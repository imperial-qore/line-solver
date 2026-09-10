/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * Delayed-hit (retrieval-system) cache analyzers, BOTH arms of mvaDispatch
 * branch 2. A miss fetches the item through the retrieval queue and returns to
 * the cache; a read arriving mid-fetch is a delayed hit.
 *
 * OPEN: the model is retrieval_simple, Source -> Cache(FIFO, 3 items, cap 1)
 * with a single-queue retrieval system (INF, Exp(2)), read distribution
 * [0.6,0.3,0.1]. The open fixed point separates hit / miss / delayed hit.
 * Reference numbers are MATLAB SolverMVA(model).getAvgTable (method fpi).
 *
 * CLOSED: Delay -> Cache -> Fetch(PS) -> Cache with two circulating jobs. There
 * is no arrival rate to run the open fixed point on, so the split comes from a
 * decomposition-aggregation sweep and the delayed-hit fraction folds into the
 * miss. Reference numbers are MATLAB's, cross-checked against SolverNC on the
 * same model.
 */

#include <cmath>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/mva/solver_mva_cacheqn_retrieval.h"
#include "line/solvers/mva/solver_mva_runner.h"

using namespace line;
using lang::Distrib;
using lang::ReplacementStrategy;
using lang::SchedStrategy;

namespace {

using D = Distrib<double>;

std::size_t station_of(const qn::NetworkStruct<double>& sn, const std::string& nm) {
    for (std::size_t i = 0; i < sn.nstations; ++i)
        if (sn.stations[i].name == nm) return i;
    FAIL("no station named ", nm);
    return 0;
}

TEST_CASE("open delayed-hit retrieval cache matches the MATLAB getAvgTable") {
    qn::Network<double> m("DelayedHits");
    const std::size_t src = m.add_source("Source");
    qn::CacheParam<double> ch;
    ch.nitems = 3;
    ch.itemcap = std::vector<int>{1};
    ch.replacestrat = ReplacementStrategy::FIFO;
    ch.pread = std::vector<std::vector<double> >{{0.6, 0.3, 0.1}, {}, {}};
    ch.hitclass = std::vector<std::size_t>{2, 0, 0};
    ch.missclass = std::vector<std::size_t>{3, 0, 0};
    const std::size_t cache = m.add_cache("Cache", ch);
    const std::size_t q = m.add_queue("Queue", SchedStrategy::INF);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t init = m.add_open_class("InitClass");
    const std::size_t hit = m.add_open_class("HitClass");
    const std::size_t miss = m.add_open_class("MissClass");
    m.set_arrival(src, init, D::exp_rate(1.0));
    m.set_service(q, init, D::exp_rate(2.0));   // read-class fetch service, inherited
    m.set_retrieval_system(cache, init, miss, std::vector<std::size_t>{q});
    qn::RoutingMatrix<double> P;
    P.set(init, init, src, cache, 1.0);
    P.set(init, init, cache, q, 1.0);
    P.set(init, init, q, cache, 1.0);
    P.set(hit, hit, cache, snk, 1.0);
    P.set(miss, miss, cache, snk, 1.0);
    m.link(P);

    mva::MvaOptions opt;
    Matrix<double> is;
    const mva::AvgResult<double> r = mva::solver_mva_run_analyzer(m.get_struct(), opt, is);
    CHECK(r.actualmethod == "fpi");
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const std::size_t Q = station_of(sn, "Queue");
    // the retrieval station carries the read class's delayed-hit occupancy
    CHECK(r.QN(Q, 0) == doctest::Approx(0.24309).epsilon(1e-4));
    CHECK(r.UN(Q, 0) == doctest::Approx(0.24309).epsilon(1e-4));
    CHECK(r.RN(Q, 0) == doctest::Approx(0.5).epsilon(1e-4));
    CHECK(r.TN(Q, 0) == doctest::Approx(0.48618).epsilon(1e-4));

    // THE CACHE ANSWER MUST RIDE WITH THE AVG RESULT. SolverMVA computed this
    // split all along and `mva_dispatch` dropped it, so `-a cache` refused
    // under `-s mva` and a host solving through the CLI kept whatever the
    // PREVIOUS solver had written onto its Cache node. Values are MATLAB's
    // `MVA(model).getAvgCacheTable` on retrieval_simple.
    REQUIRE(r.cache.caches.size() == 1);
    const solvers::CacheNodeMetrics<double>& cm = r.cache.caches[0];
    // Keyed by NAME, because the node index is this struct's own order and the
    // host's differs; see cache_metrics.h.
    CHECK(cm.name == "Cache");
    CHECK(cm.nitems == 3);
    REQUIRE(cm.hitprob.size() == sn.nclasses);
    CHECK(cm.hitprob[0] == doctest::Approx(0.413339).epsilon(1e-5));
    CHECK(cm.delayedprob[0] == doctest::Approx(0.100479).epsilon(1e-5));
    CHECK(cm.missprob[0] == doctest::Approx(0.486182).epsilon(1e-5));
    CHECK(cm.latency[0] == doctest::Approx(0.5).epsilon(1e-5));
    // The three fractions partition the reads, which is what makes
    // ArvR = arvr*(missprob + delayedprob) meaningful beside ResidT.
    CHECK(cm.hitprob[0] + cm.delayedprob[0] + cm.missprob[0] ==
          doctest::Approx(1.0).epsilon(1e-9));
    // A class that does not read the cache is NaN, not zero: absent must be
    // distinguishable from "never hits", or getAvgCacheTable prints a row.
    CHECK(std::isnan(cm.hitprob[1]));
    CHECK(std::isnan(cm.missprob[1]));
    // per-list hit fractions for the read class sum to the aggregate hit
    REQUIRE(cm.hitproblist.rows() == sn.nclasses);
    double lsum = 0.0;
    for (std::size_t l = 0; l < cm.hitproblist.cols(); ++l) lsum += cm.hitproblist(0, l);
    CHECK(lsum == doctest::Approx(0.413339).epsilon(1e-5));
    // Per-item occupancy, `[pi0(:), phit.']` of the reference: column 0 is the
    // miss and columns 1.. the per-list residency. A ROW DOES NOT SUM TO ONE
    // here, and that is the reference's shape, not a defect: the delayed-hit
    // fraction phi_i is neither a miss nor a residency, so a row sums to
    // 1 - phi_i. What ties it to the aggregates is the ACCESS-WEIGHTED sum,
    // weights being the read law [0.6, 0.3, 0.1].
    REQUIRE(cm.itemprob.rows() == 3);
    REQUIRE(cm.itemprob.cols() == 2);
    const double w[3] = {0.6, 0.3, 0.1};
    double wmiss = 0.0, whit = 0.0;
    for (std::size_t i = 0; i < 3; ++i) {
        wmiss += w[i] * cm.itemprob(i, 0);
        whit += w[i] * cm.itemprob(i, 1);
        CHECK(cm.itemprob(i, 0) + cm.itemprob(i, 1) <= 1.0);
    }
    CHECK(wmiss == doctest::Approx(cm.missprob[0]).epsilon(1e-9));
    CHECK(whit == doctest::Approx(cm.hitprob[0]).epsilon(1e-9));
    // A more popular item is cached more often; item 1 is read 6x as often as 3
    CHECK(cm.itemprob(0, 1) > cm.itemprob(2, 1));
}

/** Delay -> Cache -> Fetch -> Cache, CLOSED, with a delayed-hit retrieval
 * system. Deliberately the SAME fixture as test_nc.cpp closed_retrieval_model(),
 * so the two solvers can be compared on it. */
qn::Network<double> closed_retrieval_model() {
    qn::Network<double> m("ClosedDelayedHits");
    const std::size_t d = m.add_delay("Delay");
    qn::CacheParam<double> ch;
    ch.nitems = 3;
    ch.itemcap = std::vector<int>{1};
    ch.replacestrat = ReplacementStrategy::FIFO;
    ch.pread = std::vector<std::vector<double> >{{0.6, 0.3, 0.1}, {}, {}};
    ch.hitclass = std::vector<std::size_t>{2, 0, 0};
    ch.missclass = std::vector<std::size_t>{3, 0, 0};
    const std::size_t cn = m.add_cache("Cache", ch);
    const std::size_t q1 = m.add_queue("Fetch", SchedStrategy::PS);
    const std::size_t job = m.add_closed_class("InitClass", 2.0, d);
    const std::size_t hit = m.add_closed_class("HitClass", 0.0, d);
    const std::size_t mis = m.add_closed_class("MissClass", 0.0, d);
    m.set_service(d, job, D::exp_rate(1.0));
    m.set_service(d, hit, D::exp_rate(1.0));
    m.set_service(d, mis, D::exp_rate(1.0));
    m.set_service(q1, job, D::exp_rate(2.0));
    m.set_retrieval_system(cn, job, mis, std::vector<std::size_t>{q1});
    qn::RoutingMatrix<double> P;
    P.set(job, job, d, cn, 1.0);
    P.set(job, job, cn, q1, 1.0);
    P.set(job, job, q1, cn, 1.0);
    P.set(hit, job, cn, d, 1.0);
    P.set(mis, job, cn, d, 1.0);
    m.link(P);
    return m;
}

TEST_CASE("the closed integrated retrieval cache solves through mvaDispatch branch 2") {
    // The CLOSED variant is a different method from the OPEN one above, not the
    // same method on another topology: with no exogenous arrival rate there is
    // nothing to run `retrieval_fpi` on, so the split comes from the
    // decomposition-aggregation sweep of `da_cacheqn_retrieval` and the
    // delayed-hit fraction folds into the miss.
    qn::Network<double> m = closed_retrieval_model();
    const qn::NetworkStruct<double>& sn = m.get_struct();
    CHECK(sn.nchains == 1);
    CHECK(sn.nclasses == 6);

    const mva::MvaOptions opt;
    const mva::MvaCacheqnRetrievalSolution<double> r =
        mva::solver_mva_cacheqn_retrieval_analyzer(sn, opt);
    CHECK(r.sol.method == "fpi");
    CHECK(r.hitprob[0] == doctest::Approx(0.423425840227865).epsilon(1e-9));
    CHECK(r.missprob[0] == doctest::Approx(0.576574159772135).epsilon(1e-9));
    // The reference folds the delayed-hit fraction into the miss on this path.
    CHECK(r.delayedprob[0] == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(r.hitprob[0] + r.missprob[0] == doctest::Approx(1.0).epsilon(1e-12));
    // Quantities this path does not compute are NaN, not a fabricated number.
    CHECK(std::isnan(r.latency[0]));
    CHECK(std::isnan(r.hitproblist(0, 0)));

    // THE ORACLE IS SolverNC ON THE SAME MODEL, and it is an independent one:
    // NC solves the converged load-dependent network by normalizing constants
    // where MVA runs the load-dependent mean-value recursion. Both are exact on
    // a product-form model, so any disagreement is a defect and not a tolerance
    // -- which is how the `solver_mvald` residence-time divisor was found. The
    // literals are the MATLAB values already pinned in test_nc.cpp.
    const double fq[3] = {0.247568598687817, 0.166148415818464, 0.0717543454721912};
    const double ft[3] = {0.445312216396303, 0.298858254605033, 0.12906760706991};
    CHECK(r.sol.Q(0, 0) == doctest::Approx(1.51452864002153).epsilon(1e-9));
    CHECK(r.sol.Tp(0, 0) == doctest::Approx(1.51452864002153).epsilon(1e-9));
    double totQ = r.sol.Q(0, 0);
    for (std::size_t k = 0; k < 3; ++k) {
        CHECK(r.sol.Q(1, 3 + k) == doctest::Approx(fq[k]).epsilon(1e-9));
        CHECK(r.sol.Tp(1, 3 + k) == doctest::Approx(ft[k]).epsilon(1e-9));
        CHECK(r.sol.Q(1, k) == doctest::Approx(0.0).epsilon(1e-12));
        totQ += r.sol.Q(1, 3 + k);
    }
    // POPULATION CONSERVATION IS THE SHARPEST CHECK HERE and it is what the
    // per-visit divisor buys: before the fix the fetch row summed to 0.2799
    // instead of 0.4855 and the table lost 0.206 of the two circulating jobs,
    // while every individual number still looked plausible.
    CHECK(totQ == doctest::Approx(2.0).epsilon(1e-9));

    // The whole runner table, which is what a caller sees, and the branch it
    // took. `fpi` names the algorithm that decided the hit/miss split, as in the
    // open case and in SolverNC.
    Matrix<double> is;
    const mva::AvgResult<double> a = mva::solver_mva_run_analyzer(sn, opt, is);
    CHECK(a.actualmethod == "fpi");
    CHECK(a.QN(0, 0) == doctest::Approx(1.51452864002153).epsilon(1e-9));
    CHECK(a.TN(0, 0) == doctest::Approx(1.51452864002153).epsilon(1e-9));
    CHECK(a.WN(0, 0) == doctest::Approx(1.0).epsilon(1e-9));
    for (std::size_t k = 0; k < 3; ++k) {
        CHECK(a.QN(1, 3 + k) == doctest::Approx(fq[k]).epsilon(1e-9));
        CHECK(a.TN(1, 3 + k) == doctest::Approx(ft[k]).epsilon(1e-9));
        CHECK(a.RN(1, 3 + k) == doctest::Approx(0.555943873921247).epsilon(1e-9));
        CHECK(a.WN(1, 3 + k) == doctest::Approx(0.83391581088187).epsilon(1e-9));
    }
}

TEST_CASE("the closed retrieval decomposition refuses exact arithmetic by name") {
    // Two tolerance-stopped solves alternate and the coalescing rate is a real
    // power, so there is no exact answer to return.
    using R = num_traits<Rational>;
    const Rational one = R::from_int(1);
    qn::Network<Rational> m("ClosedDelayedHitsExact");
    const std::size_t d = m.add_delay("Delay");
    qn::CacheParam<Rational> ch;
    ch.nitems = 3;
    ch.itemcap = std::vector<int>{1};
    ch.replacestrat = ReplacementStrategy::FIFO;
    ch.pread = std::vector<std::vector<Rational> >{
        {Rational(R::from_int(3) / R::from_int(5)), Rational(R::from_int(3) / R::from_int(10)),
         Rational(one / R::from_int(10))},
        {},
        {}};
    ch.hitclass = std::vector<std::size_t>{2, 0, 0};
    ch.missclass = std::vector<std::size_t>{3, 0, 0};
    const std::size_t cn = m.add_cache("Cache", ch);
    const std::size_t q1 = m.add_queue("Fetch", SchedStrategy::PS);
    const std::size_t job = m.add_closed_class("InitClass", 2.0, d);
    const std::size_t hit = m.add_closed_class("HitClass", 0.0, d);
    const std::size_t mis = m.add_closed_class("MissClass", 0.0, d);
    m.set_service(d, job, Distrib<Rational>::exp_rate(one));
    m.set_service(d, hit, Distrib<Rational>::exp_rate(one));
    m.set_service(d, mis, Distrib<Rational>::exp_rate(one));
    m.set_service(q1, job, Distrib<Rational>::exp_rate(R::from_int(2)));
    m.set_retrieval_system(cn, job, mis, std::vector<std::size_t>{q1});
    qn::RoutingMatrix<Rational> P;
    P.set(job, job, d, cn, one);
    P.set(job, job, cn, q1, one);
    P.set(job, job, q1, cn, one);
    P.set(hit, job, cn, d, one);
    P.set(mis, job, cn, d, one);
    m.link(P);
    const mva::MvaOptions opt;
    CHECK_THROWS_AS(mva::solver_mva_cacheqn_retrieval_analyzer(m.get_struct(), opt),
                    UnsupportedError);
}

TEST_CASE("the closed retrieval branch carries its cache split into AvgResult") {
    // The analyzer's own numbers are asserted above; what this covers is the
    // WIRING, dispatch -> runner -> AvgResult::cache, which is what a host
    // solving through line-cli actually reads. Dropped here, MATLAB's Cache
    // node kept the previous solver's split and getAvgCacheTable reported it.
    qn::Network<double> m = closed_retrieval_model();
    const mva::MvaOptions opt;
    Matrix<double> is;
    const mva::AvgResult<double> r = mva::solver_mva_run_analyzer(m.get_struct(), opt, is);
    REQUIRE(r.cache.caches.size() == 1);
    const solvers::CacheNodeMetrics<double>& cm = r.cache.caches[0];
    CHECK(cm.name == "Cache");
    CHECK(cm.hitprob[0] == doctest::Approx(0.423425840227865).epsilon(1e-9));
    CHECK(cm.missprob[0] == doctest::Approx(0.576574159772135).epsilon(1e-9));
    CHECK(cm.hitprob[0] + cm.missprob[0] == doctest::Approx(1.0).epsilon(1e-12));
    // This path computes no latency, and NaN is how it says so.
    CHECK(std::isnan(cm.latency[0]));
}

}  // namespace
