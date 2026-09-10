/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The integrated caching-queueing network under the fluid solver: end-to-end
 * runs of `solver_fld_cacheqn_analyzer` and `solver_fld_cacheqn_tran` on a
 * Network built here, with nothing hand-assembled in between.
 *
 * WHY THE ORACLES ARE IDENTITIES AND NOT A REFERENCE TABLE. The fluid solver is
 * an APPROXIMATION twice over: the queueing network is replaced by its mean
 * drift, and the cache miss rates come from a refined (1/N-accurate) mean field.
 * Neither reproduces MVA or NC except in the limit, so pinning a queue length
 * against an exact solver would pin the approximation error and not the port.
 * What the model does fix exactly, whatever the drift does with it, is:
 *
 *   1. hit probability + miss probability = 1 at every cache and reading class,
 *      because the decomposition defines the hit probability as 1 - missrate /
 *      arrival rate and the mean field clips the per-item miss probability to
 *      [0,1];
 *   2. a cache whose capacity equals its item count NEVER misses -- every item
 *      is resident, the drift out of the "not cached" level is identically
 *      zero, and it stays zero for all t;
 *   3. a COLD cache misses everything at t = 0, since every item starts at the
 *      "not cached" level;
 *   4. the Source reports the arrival rate the model was given, which the
 *      decomposition cannot change: it rewrites the cache's routing, not the
 *      external arrivals.
 *
 * The queue lengths themselves are deliberately not asserted. The C++
 * `da_cacheqn` over-routes the cache switch to every connected node, as the
 * reference does, which inflates the pass-through throughputs; see the note in
 * test_mva_cacheqn.cpp.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/fluid/fluid_cacheqn.h"

using namespace line;
using lang::Distrib;
using lang::ReplacementStrategy;
using lang::SchedStrategy;

namespace {

using D = Distrib<double>;

/** A decaying read profile; a UNIFORM one makes every replacement policy agree. */
std::vector<double> fldcq_zipf(std::size_t n, double alpha = 1.0) {
    std::vector<double> p(n, 0.0);
    double tot = 0.0;
    for (std::size_t k = 0; k < n; ++k) {
        p[k] = 1.0 / std::pow(static_cast<double>(k + 1), alpha);
        tot += p[k];
    }
    for (std::size_t k = 0; k < n; ++k) p[k] /= tot;
    return p;
}

std::size_t fldcq_station_of(const qn::NetworkStruct<double>& sn, const std::string& nm) {
    for (std::size_t i = 0; i < sn.nstations; ++i)
        if (sn.stations[i].name == nm) return i;
    FAIL("no station named ", nm);
    return 0;
}

/**
 * Source -> Cache -> HitQueue / MissQueue -> Sink, three open classes: InitClass
 * reads the cache, HitClass and MissClass are the switched classes. This is the
 * gallery cache-routing shape, the same model test_mva_cacheqn.cpp uses, with
 * the replacement policy and the cache geometry left to the caller.
 */
qn::Network<double> fldcq_model(
    ReplacementStrategy strat, std::size_t nitems, int cap, double arrival = 1.0,
    const std::vector<std::vector<Matrix<double> > >& accost =
        std::vector<std::vector<Matrix<double> > >(),
    const std::vector<int>& caps = std::vector<int>(),
    const std::vector<double>& pread = std::vector<double>()) {
    qn::Network<double> m("fld-cache-routing");
    const std::size_t src = m.add_source("Source");
    qn::CacheParam<double> ch;
    ch.nitems = nitems;
    ch.itemcap = caps.empty() ? std::vector<int>{cap} : caps;
    ch.replacestrat = strat;
    ch.accost = accost;
    // InitClass reads with the given profile, uniform by default; the switched
    // classes do not read (empty rows).
    ch.pread = std::vector<std::vector<double> >{
        pread.empty() ? std::vector<double>(nitems, 1.0 / static_cast<double>(nitems)) : pread,
        {},
        {}};
    ch.hitclass = std::vector<std::size_t>{2, 0, 0};   // InitClass -> HitClass
    ch.missclass = std::vector<std::size_t>{3, 0, 0};  //           -> MissClass
    const std::size_t cache = m.add_cache("Cache", ch);
    const std::size_t hq = m.add_queue("HitQueue", SchedStrategy::PS);
    const std::size_t mq = m.add_queue("MissQueue", SchedStrategy::PS);
    const std::size_t snk = m.add_sink("Sink");
    const std::size_t init = m.add_open_class("InitClass");
    const std::size_t hit = m.add_open_class("HitClass");
    const std::size_t miss = m.add_open_class("MissClass");
    m.set_arrival(src, init, D::exp_rate(arrival));
    m.set_service(hq, hit, D::exp_rate(4.0));
    m.set_service(mq, miss, D::exp_rate(2.0));
    qn::RoutingMatrix<double> P;
    P.set(init, init, src, cache, 1.0);
    P.set(hit, hit, cache, hq, 1.0);
    P.set(hit, hit, hq, snk, 1.0);
    P.set(miss, miss, cache, mq, 1.0);
    P.set(miss, miss, mq, snk, 1.0);
    m.link(P);
    return m;
}

}  // namespace

TEST_CASE("fld cacheqn: the converged split is a probability and the arrivals survive it") {
    qn::Network<double> m = fldcq_model(ReplacementStrategy::RR, 4, 2, 1.0);
    const fluid::FluidCacheqnSolution<double> s =
        fluid::solver_fld_cacheqn_analyzer(m.get_struct(), fluid::FluidOptions());

    REQUIRE(s.hitprob.rows() == 1);
    REQUIRE(s.hitprob.cols() == 3);
    // Identity 1, on the reading class. Not a tolerance on the fluid answer:
    // the decomposition constructs the pair to sum to one.
    CHECK(s.hitprob(0, 0) + s.missprob(0, 0) == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(s.hitprob(0, 0) >= 0.0);
    CHECK(s.hitprob(0, 0) <= 1.0);
    CHECK(s.missprob(0, 0) >= 0.0);
    CHECK(s.missprob(0, 0) <= 1.0);
    // A class that never reads the cache has neither probability, and the
    // reference reports zero for both rather than the degenerate 1 and 0.
    CHECK(s.hitprob(0, 1) == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(s.missprob(0, 1) == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(s.hitprob(0, 2) == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(s.missprob(0, 2) == doctest::Approx(0.0).epsilon(1e-12));

    // Identity 4: the decomposition rewrites the cache's routing, never the
    // external arrival process, so the Source still reports 1.0.
    const qn::NetworkStruct<double>& sn = m.get_struct();
    CHECK(s.sol.TN(fldcq_station_of(sn, "Source"), 0) == doctest::Approx(1.0).epsilon(1e-6));
    CHECK(s.iter >= 1);
    CHECK(s.sol.method == "rmf");
    for (std::size_t i = 0; i < sn.nstations; ++i)
        for (std::size_t r = 0; r < sn.nclasses; ++r) {
            CHECK(std::isfinite(s.sol.QN(i, r)));
            CHECK(s.sol.QN(i, r) >= 0.0);
            CHECK(std::isfinite(s.sol.TN(i, r)));
            CHECK(s.sol.TN(i, r) >= 0.0);
        }
}

TEST_CASE("fld cacheqn: a cache that holds every item never misses") {
    // Identity 2. With capacity equal to the item count every item is resident,
    // the drift out of the not-cached level is identically zero, and the mean
    // field returns a zero miss rate however long it integrates.
    qn::Network<double> m = fldcq_model(ReplacementStrategy::RR, 3, 3, 1.0);
    const fluid::FluidCacheqnSolution<double> s =
        fluid::solver_fld_cacheqn_analyzer(m.get_struct(), fluid::FluidOptions());
    CHECK(s.missprob(0, 0) == doctest::Approx(0.0).epsilon(1e-8));
    CHECK(s.hitprob(0, 0) == doctest::Approx(1.0).epsilon(1e-8));
}

TEST_CASE("fld cacheqn: FIFO(m) on the linear chain is solved as RANDOM(m)") {
    // Gast15 Thm 1: pi_FIFO(m) = pi_RAND(m) on the linear access graph, so the
    // two are the SAME steady-state computation, not two approximations that
    // happen to agree. The oracle is that identity, checked against the RR run
    // of the same geometry rather than against a stored number.
    qn::Network<double> mr = fldcq_model(ReplacementStrategy::RR, 4, 2, 1.0);
    qn::Network<double> mf = fldcq_model(ReplacementStrategy::FIFO, 4, 2, 1.0);
    const fluid::FluidCacheqnSolution<double> sr =
        fluid::solver_fld_cacheqn_analyzer(mr.get_struct(), fluid::FluidOptions());
    const fluid::FluidCacheqnSolution<double> sf =
        fluid::solver_fld_cacheqn_analyzer(mf.get_struct(), fluid::FluidOptions());
    CHECK(sf.missprob(0, 0) == doctest::Approx(sr.missprob(0, 0)).epsilon(1e-12));
    CHECK(sf.hitprob(0, 0) + sf.missprob(0, 0) == doctest::Approx(1.0).epsilon(1e-12));
}

TEST_CASE("fld cacheqn: a policy without a drift is refused, never approximated") {
    // LRU, HLRU, CLIMB and QLRU are solved elsewhere by the characteristic-time
    // fixed point, which is not a fluid method; the reference refuses them here
    // and so does the port.
    const fluid::FluidOptions o;
    qn::Network<double> lru = fldcq_model(ReplacementStrategy::LRU, 4, 2, 1.0);
    CHECK_THROWS_AS(fluid::solver_fld_cacheqn_analyzer(lru.get_struct(), o), UnsupportedError);
    qn::Network<double> hlru = fldcq_model(ReplacementStrategy::HLRU, 4, 2, 1.0);
    CHECK_THROWS_AS(fluid::solver_fld_cacheqn_analyzer(hlru.get_struct(), o), UnsupportedError);
    qn::Network<double> climb = fldcq_model(ReplacementStrategy::CLIMB, 4, 2, 1.0);
    CHECK_THROWS_AS(fluid::solver_fld_cacheqn_analyzer(climb.get_struct(), o), UnsupportedError);
    qn::Network<double> qlru = fldcq_model(ReplacementStrategy::QLRU, 4, 2, 1.0);
    CHECK_THROWS_AS(fluid::solver_fld_cacheqn_analyzer(qlru.get_struct(), o), UnsupportedError);
}

TEST_CASE("fld cacheqn: strict FIFO(m) takes its own drift, not RANDOM(m)'s") {
    // Gast15 proves pi_FIFO(m) = pi_RAND(m) and shows strict FIFO(m) DIFFERS.
    // The oracle is therefore a disagreement: solving strict FIFO through
    // `cache_miss_sfifo_rmf` must not reproduce the RR number, or the
    // position-resolved drift is not being reached.
    //
    // THE POPULARITY MUST BE NON-UNIFORM AND THE CACHE MULTI-LIST, or the
    // oracle is vacuous. Under a uniform profile every item is exchangeable, so
    // every replacement policy holds the same expected number of items and they
    // all agree; and strict FIFO degenerates to FIFO when every list below the
    // top holds one item, because there is then no within-list order to
    // disagree about. Two lists of two, with a decaying profile, has both.
    const std::vector<int> caps{2, 2};
    const std::vector<double> pop = fldcq_zipf(6);
    qn::Network<double> mr = fldcq_model(ReplacementStrategy::RR, 6, 0, 1.0,
                                         std::vector<std::vector<Matrix<double> > >(), caps, pop);
    qn::Network<double> ms = fldcq_model(ReplacementStrategy::SFIFO, 6, 0, 1.0,
                                         std::vector<std::vector<Matrix<double> > >(), caps, pop);
    const fluid::FluidCacheqnSolution<double> sr =
        fluid::solver_fld_cacheqn_analyzer(mr.get_struct(), fluid::FluidOptions());
    const fluid::FluidCacheqnSolution<double> ss =
        fluid::solver_fld_cacheqn_analyzer(ms.get_struct(), fluid::FluidOptions());
    CHECK(ss.hitprob(0, 0) + ss.missprob(0, 0) == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(ss.missprob(0, 0) != doctest::Approx(sr.missprob(0, 0)).epsilon(1e-6));
}

TEST_CASE("fld cacheqn: a declared access graph reaches the drift") {
    // A graph that ADMITS A MISS STRAIGHT TO THE TOP LIST, instead of to list 1,
    // changes which items the cache retains under a decaying popularity profile,
    // so the miss probability must move. Solving it as the linear chain -- which
    // is what a dropped `accost` would do -- reproduces the no-graph number
    // exactly, so this separates "the graph reached the drift" from "the graph
    // was silently ignored" with no reference table.
    //
    // It cannot be checked on a ONE-list cache: with h = 1 the only edges are
    // admit and stay, and a graph that refuses admission simply freezes the
    // drift at its pre-filled initial occupancy -- which is the reference's
    // behaviour too, and not a number that distinguishes anything.
    const std::size_t nitems = 6;
    const std::size_t h = 2;
    const std::vector<int> caps{2, 2};
    const std::vector<double> pop = fldcq_zipf(nitems);
    // The graph must stay a TREE over the lists: `cache_gamma_lp` refuses a list
    // with two parents, and that is a constraint of the ISOLATION step, not of
    // the drift. Keeping admission on list 1 and splitting the promotion out of
    // list 1 between staying and moving up leaves list 1 with parent `out` and
    // list 2 with parent list 1, while differing from the chain (which promotes
    // with probability one).
    Matrix<double> g(h + 1, h + 1, 0.0);
    g(0, 1) = 1.0;  // a miss enters list 1, as in the chain
    g(1, 1) = 0.4;  // a hit in list 1 STAYS with probability 0.4
    g(1, 2) = 0.6;  // and promotes with 0.6, where the chain always promotes
    g(2, 2) = 1.0;  // a hit in the top list stays
    const std::vector<std::vector<Matrix<double> > > accost(
        3, std::vector<Matrix<double> >(nitems, g));
    qn::Network<double> mg =
        fldcq_model(ReplacementStrategy::RR, nitems, 0, 1.0, accost, caps, pop);
    qn::Network<double> ml = fldcq_model(ReplacementStrategy::RR, nitems, 0, 1.0,
                                         std::vector<std::vector<Matrix<double> > >(), caps, pop);
    const fluid::FluidCacheqnSolution<double> sg =
        fluid::solver_fld_cacheqn_analyzer(mg.get_struct(), fluid::FluidOptions());
    const fluid::FluidCacheqnSolution<double> sl =
        fluid::solver_fld_cacheqn_analyzer(ml.get_struct(), fluid::FluidOptions());
    CHECK(sg.hitprob(0, 0) + sg.missprob(0, 0) == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(sg.missprob(0, 0) != doctest::Approx(sl.missprob(0, 0)).epsilon(1e-6));
}

TEST_CASE("fld cacheqn tran: a cold cache misses everything at t = 0") {
    // Identity 3. Seeding every item at the not-cached level makes the miss
    // probability exactly one at the first point of the trajectory, whatever the
    // popularity profile or the arrival rates the decomposition converged on.
    qn::Network<double> m = fldcq_model(ReplacementStrategy::RR, 4, 2, 1.0);
    const std::size_t nitems = 4, h = 1;
    std::vector<double> cold(nitems * (h + 1), 0.0);
    for (std::size_t i = 0; i < nitems; ++i) cold[i] = 1.0;  // level 0 = not cached

    const std::vector<fluid::FluidCacheqnTranCache<double> > tr =
        fluid::solver_fld_cacheqn_tran(m.get_struct(), fluid::FluidOptions(), 0.0, 20.0,
                                       std::vector<std::vector<double> >{cold});
    REQUIRE(tr.size() == 1);
    const fluid::FluidCacheqnTranCache<double>& c = tr[0];
    REQUIRE(c.t.size() >= 2);
    CHECK(c.t.front() == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(c.t.back() == doctest::Approx(20.0).epsilon(1e-9));
    CHECK(c.missprob_t(0, 0) == doctest::Approx(1.0).epsilon(1e-9));
    CHECK(c.hitprob_t(0, 0) == doctest::Approx(0.0).epsilon(1e-9));
    CHECK(c.arate[0] > 0.0);

    // Identity 1 again, now at every point of the trajectory.
    for (std::size_t j = 0; j < c.t.size(); ++j) {
        CHECK(c.hitprob_t(0, j) + c.missprob_t(0, j) == doctest::Approx(1.0).epsilon(1e-12));
        CHECK(c.missprob_t(0, j) >= 0.0);
        CHECK(c.missprob_t(0, j) <= 1.0);
    }
    // The occupancy is the flat DDPP state, one row per (item, level).
    CHECK(c.xocc.rows() == nitems * (h + 1));
    CHECK(c.xocc.cols() == c.t.size());
}

TEST_CASE("fld cacheqn tran: a full cache misses nothing at any time") {
    // Identity 2 along the whole trajectory: the default seed puts all three
    // items in the single list of a capacity-3 cache, and the drift has no way
    // to move any of them out.
    qn::Network<double> m = fldcq_model(ReplacementStrategy::RR, 3, 3, 1.0);
    const std::vector<fluid::FluidCacheqnTranCache<double> > tr =
        fluid::solver_fld_cacheqn_tran(m.get_struct(), fluid::FluidOptions(), 0.0, 10.0);
    REQUIRE(tr.size() == 1);
    const fluid::FluidCacheqnTranCache<double>& c = tr[0];
    for (std::size_t j = 0; j < c.t.size(); ++j) {
        CHECK(c.missprob_t(0, j) == doctest::Approx(0.0).epsilon(1e-8));
        CHECK(c.hitprob_t(0, j) == doctest::Approx(1.0).epsilon(1e-8));
    }
}

TEST_CASE("fld cacheqn tran: FIFO(m) takes the position-resolved drift, not RANDOM(m)'s") {
    // RANDOM(m) and FIFO(m) agree on the stationary occupancy but NOT on the
    // path to it: FIFO evicts the deterministic tail (residence exactly m
    // insertions) where RANDOM draws a uniform victim (geometric residence).
    // The oracle is that the two trajectories differ somewhere while both start
    // from a cold cache at miss probability one, which is what identity 3 fixes.
    const std::size_t nitems = 6;
    const std::vector<int> caps{2, 2};
    const std::vector<double> pop = fldcq_zipf(nitems);
    qn::Network<double> mf = fldcq_model(ReplacementStrategy::FIFO, nitems, 0, 1.0,
                                         std::vector<std::vector<Matrix<double> > >(), caps, pop);
    qn::Network<double> mr = fldcq_model(ReplacementStrategy::RR, nitems, 0, 1.0,
                                         std::vector<std::vector<Matrix<double> > >(), caps, pop);
    // The two families do NOT share a state space: RANDOM(m) tracks
    // nitems*(h+1) per-list occupancies, FIFO(m) nitems*sum(m) per-SLOT ones.
    // Seeding each with its own cold vector is what makes t = 0 comparable.
    const std::vector<std::vector<double> > coldf(
        1, std::vector<double>(nitems * (static_cast<std::size_t>(caps[0]) +
                                         static_cast<std::size_t>(caps[1])), 0.0));
    // A cold RANDOM(m) state is NOT the zero vector: its level 0 is "not
    // cached" and carries the whole mass, whereas the positional state tracks
    // only the in-cache slots and is empty when the cache is.
    std::vector<std::vector<double> > coldr(1, std::vector<double>(nitems * (caps.size() + 1), 0.0));
    for (std::size_t k = 0; k < nitems; ++k)
        coldr[0][cache::cache_miss_rmf_index(k, 0, nitems)] = 1.0;
    const std::vector<fluid::FluidCacheqnTranCache<double> > tf =
        fluid::solver_fld_cacheqn_tran(mf.get_struct(), fluid::FluidOptions(), 0.0, 5.0, coldf);
    const std::vector<fluid::FluidCacheqnTranCache<double> > tr =
        fluid::solver_fld_cacheqn_tran(mr.get_struct(), fluid::FluidOptions(), 0.0, 5.0, coldr);
    REQUIRE(tf.size() == 1);
    REQUIRE(tr.size() == 1);
    REQUIRE(tf[0].missprob_t.cols() > 1);
    // Both start cold, so the first point is a full miss for either policy.
    CHECK(tf[0].missprob_t(0, 0) == doctest::Approx(1.0).epsilon(1e-9));
    CHECK(tr[0].missprob_t(0, 0) == doctest::Approx(1.0).epsilon(1e-9));
    // And the two policies must part company on the way to the fixed point.
    bool differs = false;
    for (std::size_t c = 1; c < tf[0].missprob_t.cols() && c < tr[0].missprob_t.cols(); ++c)
        if (std::fabs(tf[0].missprob_t(0, c) - tr[0].missprob_t(0, c)) > 1e-6) differs = true;
    CHECK(differs);
}

TEST_CASE("fld cacheqn: the moment closure answers a cache model") {
    // "minnormal" used to be refused on any model with a cache node: the closing
    // family read sn.rt station-major, which on a cache model drops the flow
    // through the cache entirely. It is now answered through this same
    // decomposition, with the closure in the network step. The cache layer is
    // the refined mean field either way, so the hit probability must be
    // IDENTICAL to the rmf route while the queueing metrics move.
    qn::Network<double> mr = fldcq_model(ReplacementStrategy::RR, 6, 3, 1.0);
    const qn::NetworkStruct<double>& sn = mr.get_struct();

    fluid::FluidOptions orm;
    orm.method = "rmf";
    const fluid::FluidCacheqnSolution<double> srm =
        fluid::solver_fld_cacheqn_analyzer(sn, orm);

    fluid::FluidOptions omn;
    omn.method = "minnormal";
    const fluid::FluidCacheqnSolution<double> smn =
        fluid::solver_fld_cacheqn_analyzer(sn, omn);

    CHECK(smn.sol.method == "minnormal");
    CHECK(smn.hitprob(0, 0) == doctest::Approx(srm.hitprob(0, 0)).epsilon(1e-9));

    const std::size_t hq = fldcq_station_of(sn, "HitQueue");
    const std::size_t mq = fldcq_station_of(sn, "MissQueue");
    // both routes answer, and neither returns the degenerate all-zero solution
    // the misread routing used to produce
    CHECK(smn.sol.QN(hq, 1) > 0.0);
    CHECK(smn.sol.QN(mq, 2) > 0.0);
    CHECK(std::isfinite(smn.sol.QN(hq, 1)));
    CHECK(std::isfinite(smn.sol.QN(mq, 2)));
}
