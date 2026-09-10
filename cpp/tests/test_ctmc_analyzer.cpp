/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * End-to-end tests of `solver_ctmc_analyzer`: a Network in, an AvgTable out,
 * with nothing hand-assembled in between.
 *
 * WHY THESE ARE NOT THE SAME AS test_ctmc_generator.cpp. Those drive the
 * generator with a state space and a synchronization list built in the test;
 * these drive the whole analyzer, so the cutoff resolution, the reducible-
 * component selection, the stationary solve and the AvgTable assembly are all
 * on the path. Each is a place a handler that passes its unit test can still
 * produce a wrong table -- a POLLING, Cache or SPN model reaching a generator
 * whose component was chosen wrongly reports plausible numbers that are the
 * stationary law of the wrong chain.
 *
 * The oracles are closed forms where one exists (M/M/1/K, a two-state marked
 * graph, a work-conserving polling system) and the ALREADY-VERIFIED exact
 * product-form solvers where one does not: a closed network that both SolverMVA
 * and SolverCTMC can represent must give the same numbers, and they arrive by
 * completely different routes.
 */
#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/num/number.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/mva/solver_mva_runner.h"

namespace qn = line::qn;
namespace ctmc = line::ctmc;
using line::lang::PollingType;
using line::lang::SchedStrategy;
using line::lang::TimingStrategy;
using Dist = line::lang::Distrib<double>;

namespace {

/** Source -> FCFS Queue (capacity K) -> Sink, one open class. */
qn::Network<double> mm1k(double lambda, double mu, int K) {
    qn::Network<double> m("mm1k");
    const std::size_t src = m.add_source("Src");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t k = m.add_sink("Sink");
    const std::size_t c = m.add_open_class("C1");
    m.set_arrival(src, c, Dist::exp_rate(lambda));
    m.set_service(q, c, Dist::exp_rate(mu));
    m.set_capacity(q, K);
    qn::RoutingMatrix<double> P;
    P.set(src, q, 1.0);
    P.set(q, k, 1.0);
    m.link(P);
    return m;
}

/** One mode of a marked graph: consume one from `from`, produce one to `to`. */
qn::TransitionParam<double> one_mode(std::size_t nnodes, std::size_t from, std::size_t to,
                                     double rate) {
    qn::TransitionParam<double> tp;
    tp.nmodes = 1;
    tp.modenames.push_back("fire");
    tp.enabling.assign(1, line::Matrix<double>(nnodes, 1, 0.0));
    tp.inhibiting.assign(1, line::Matrix<double>(nnodes, 1, std::numeric_limits<double>::infinity()));
    tp.firing.assign(1, line::Matrix<double>(nnodes, 1, 0.0));
    tp.enabling[0](from - 1, 0) = 1.0;
    tp.firing[0](to - 1, 0) = 1.0;
    tp.nmodeservers.push_back(1.0);
    tp.firingphases.push_back(1);
    tp.timing.push_back(TimingStrategy::TIMED);
    tp.fireweight.push_back(1.0);
    tp.firingproc.push_back(Dist::exp_rate(rate));
    return tp;
}

}  // namespace

TEST_CASE("the analyzer solves an M/M/1/K end to end") {
    const double lambda = 0.6, mu = 1.0;
    const int K = 4;
    qn::Network<double> m = mm1k(lambda, mu, K);
    ctmc::CtmcOptions opt;
    opt.cutoff = K;
    const line::mva::AvgResult<double> r = ctmc::solver_ctmc_run_analyzer(m.get_struct(), opt);

    const double rho = lambda / mu;
    double norm = 0, q = 0;
    for (int i = 0; i <= K; ++i) norm += std::pow(rho, i);
    for (int i = 0; i <= K; ++i) q += i * std::pow(rho, i) / norm;
    const double p0 = 1.0 / norm;

    CHECK(r.QN(1, 0) == doctest::Approx(q).epsilon(1e-9));
    CHECK(r.UN(1, 0) == doctest::Approx(1.0 - p0).epsilon(1e-9));
    CHECK(r.TN(1, 0) == doctest::Approx(mu * (1.0 - p0)).epsilon(1e-9));
    // ArvR is the OFFERED rate, not the carried one, and MATLAB reports the
    // same 0.6 here beside a Tput of 0.566. It comes from
    // `sn_get_arvr_from_tput`, which multiplies the SOURCE throughput by the
    // visit ratio, so the arrivals lost at the full buffer are still counted as
    // offered. Utilization is the quantity that must be carried instead, which
    // is what the `canDropClass` guard in `solver_ctmc_avg_from_pi` is for.
    CHECK(r.AN(1, 0) == doctest::Approx(lambda).epsilon(1e-9));
    // One visit per job, so residence time and response time coincide.
    CHECK(r.WN(1, 0) == doctest::Approx(r.RN(1, 0)).epsilon(1e-12));
    CHECK(r.actualmethod == "default");
}

TEST_CASE("an open model with no cutoff takes the reference's automatic one") {
    qn::Network<double> m = mm1k(0.5, 1.0, 20);
    ctmc::CtmcOptions opt;  // no cutoff given
    const ctmc::CtmcSolution<double> d = ctmc::solver_ctmc_analyzer(m.get_struct(), opt);

    // ceil(6000^(1/(M*K))) with M=2 stations and K=1 class is 78, capped by the
    // station's own capacity of 20; the point is that it is FINITE and reported,
    // because the answer is the stationary law of the truncated chain.
    REQUIRE(d.cutoff.size() == 1);
    CHECK(d.cutoff[0] > 0);
    CHECK(!d.chain.space.empty());
    double tot = 0;
    for (std::size_t i = 0; i < d.pi.size(); ++i) tot += d.pi[i];
    CHECK(tot == doctest::Approx(1.0).epsilon(1e-12));
}

TEST_CASE("a closed product-form network agrees with exact MVA") {
    // Delay -> FCFS Queue -> Delay, three jobs. Both solvers represent it
    // exactly and reach the answer by entirely different routes: MVA by the
    // arrival theorem, CTMC by enumerating 4 states and solving pi Q = 0.
    qn::Network<double> m("cqn");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", 3.0, d);
    m.set_service(d, c, Dist::exp_rate(1.0 / 2.0));
    m.set_service(q, c, Dist::exp_rate(3.0));
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);

    const qn::NetworkStruct<double>& sn = m.get_struct();
    const line::mva::AvgResult<double> rc =
        ctmc::solver_ctmc_run_analyzer(sn, ctmc::CtmcOptions());

    line::mva::MvaOptions mopt;
    line::Matrix<double> init;
    const line::mva::AvgResult<double> rm = line::mva::solver_mva_run_analyzer(sn, mopt, init);

    for (std::size_t i = 0; i < sn.nstations; ++i) {
        CHECK(rc.QN(i, 0) == doctest::Approx(rm.QN(i, 0)).epsilon(1e-8));
        CHECK(rc.UN(i, 0) == doctest::Approx(rm.UN(i, 0)).epsilon(1e-8));
        CHECK(rc.TN(i, 0) == doctest::Approx(rm.TN(i, 0)).epsilon(1e-8));
        CHECK(rc.RN(i, 0) == doctest::Approx(rm.RN(i, 0)).epsilon(1e-8));
    }
    // Population is conserved exactly, which no approximation guarantees.
    CHECK(rc.QN(0, 0) + rc.QN(1, 0) == doctest::Approx(3.0).epsilon(1e-9));
}

TEST_CASE("the analyzer runs in exact arithmetic and matches the double run") {
    const int K = 3;
    qn::Network<double> md = mm1k(0.5, 1.0, K);
    ctmc::CtmcOptions opt;
    opt.cutoff = K;
    const line::mva::AvgResult<double> rd = ctmc::solver_ctmc_run_analyzer(md.get_struct(), opt);

    // The same model over the rationals. Every step from the generator to the
    // means is a field operation, so this is the EXACT stationary law -- the
    // agreement below is a check on the double path's conditioning, not on the
    // exact path's accuracy.
    qn::Network<line::Rational> mr("mm1k");
    const std::size_t src = mr.add_source("Src");
    const std::size_t q = mr.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t sk = mr.add_sink("Sink");
    const std::size_t c = mr.add_open_class("C1");
    mr.set_arrival(src, c, line::lang::Distrib<line::Rational>::exp_rate(
                               line::num_traits<line::Rational>::from_double(0.5)));
    mr.set_service(q, c, line::lang::Distrib<line::Rational>::exp_rate(
                             line::num_traits<line::Rational>::from_int(1)));
    mr.set_capacity(q, K);
    qn::RoutingMatrix<line::Rational> P;
    P.set(src, q, line::num_traits<line::Rational>::from_int(1));
    P.set(q, sk, line::num_traits<line::Rational>::from_int(1));
    mr.link(P);
    const line::mva::AvgResult<line::Rational> rr =
        ctmc::solver_ctmc_run_analyzer(mr.get_struct(), opt);

    CHECK(line::num_traits<line::Rational>::to_double(rr.QN(1, 0)) ==
          doctest::Approx(rd.QN(1, 0)).epsilon(1e-12));
    CHECK(line::num_traits<line::Rational>::to_double(rr.UN(1, 0)) ==
          doctest::Approx(rd.UN(1, 0)).epsilon(1e-12));
}

TEST_CASE("SolverCTMC refuses an unlisted method, and runs gpu on the CPU") {
    qn::Network<double> m = mm1k(0.5, 1.0, 2);
    ctmc::CtmcOptions bad;
    bad.method = "amva";
    CHECK_THROWS_AS(ctmc::solver_ctmc_run_analyzer(m.get_struct(), bad), line::UnsupportedError);
    // 'gpu' FALLS BACK rather than refusing, because that is what the reference
    // does: `ctmc_solve.m` wraps the gpuArray solve in a try/catch and, with no
    // GPU, warns and runs the same `Qnnz' \ bnnz` the default takes. It must
    // therefore answer, and answer exactly what the default answers; the
    // fallback is reported through CtmcSolution::warning.
    ctmc::CtmcOptions gpu;
    gpu.method = "gpu";
    const line::mva::AvgResult<double> rgpu = ctmc::solver_ctmc_run_analyzer(m.get_struct(), gpu);
    ctmc::CtmcOptions def;
    const line::mva::AvgResult<double> rdef = ctmc::solver_ctmc_run_analyzer(m.get_struct(), def);
    CHECK(rgpu.QN(1, 0) == doctest::Approx(rdef.QN(1, 0)).epsilon(1e-12));
    CHECK(rgpu.UN(1, 0) == doctest::Approx(rdef.UN(1, 0)).epsilon(1e-12));
    const ctmc::CtmcSolution<double> sgpu =
        ctmc::solver_ctmc_analyzer(m.get_struct(), gpu);
    CHECK(sgpu.warning.find("Switching to default method") != std::string::npos);
}

TEST_CASE("an SPN marked graph reaches its two-state balance end to end") {
    // P1 -T1-> P2 -T2-> P1 with one token: a two-state chain whose stationary
    // law is pi(P1) = r2/(r1+r2). The whole analyzer runs, so the reducible
    // handling and the global-sync path are both exercised.
    const double r1 = 2.0, r2 = 3.0;
    qn::Network<double> m("marked-graph");
    const std::size_t p1 = m.add_place("P1");
    const std::size_t p2 = m.add_place("P2");
    const std::size_t c = m.add_closed_class("Tok", 1.0, p1);
    m.set_service(p1, c, Dist::exp_rate(1.0));
    m.set_service(p2, c, Dist::exp_rate(1.0));
    m.add_transition("T1", one_mode(4, p1, p2, r1));
    m.add_transition("T2", one_mode(4, p2, p1, r2));

    const ctmc::CtmcSolution<double> d =
        ctmc::solver_ctmc_analyzer(m.get_struct(), ctmc::CtmcOptions());

    // Aggregate by the marking of P1: the token is either there or at P2.
    double at_p1 = 0, at_p2 = 0;
    for (std::size_t s = 0; s < d.chain.space.size(); ++s) {
        const double n1 = d.chain.space[s].local[0][d.chain.space[s].local[0].size() - 1];
        if (n1 > 0.5) at_p1 += d.pi[s];
        else at_p2 += d.pi[s];
    }
    CHECK(at_p1 + at_p2 == doctest::Approx(1.0).epsilon(1e-9));
    CHECK(at_p1 == doctest::Approx(r2 / (r1 + r2)).epsilon(1e-6));
    // The token is never destroyed: the two places hold it between them.
    CHECK(d.avg.QN(0, 0) + d.avg.QN(1, 0) == doctest::Approx(1.0).epsilon(1e-6));
}

TEST_CASE("a polling station agrees with the MATLAB CTMC, table for table") {
    // Two open classes into one POLLING station under EXHAUSTIVE service with
    // no switchover, capacity 6, cutoff 6. The expected values below are
    // MATLAB's own `SolverCTMC(model,'cutoff',6).getAvgTable` on the identical
    // model -- the same truncation, so the two are comparable to solver
    // tolerance and not merely to the untruncated limit.
    const double l1 = 0.2, l2 = 0.3, mu = 2.0;
    qn::Network<double> m("polling");
    const std::size_t src = m.add_source("Src");
    const std::size_t q = m.add_queue("P", SchedStrategy::POLLING);
    const std::size_t sk = m.add_sink("Sink");
    const std::size_t c1 = m.add_open_class("C1");
    const std::size_t c2 = m.add_open_class("C2");
    m.set_arrival(src, c1, Dist::exp_rate(l1));
    m.set_arrival(src, c2, Dist::exp_rate(l2));
    m.set_service(q, c1, Dist::exp_rate(mu));
    m.set_service(q, c2, Dist::exp_rate(mu));
    m.set_polling(q, PollingType::EXHAUSTIVE);
    qn::RoutingMatrix<double> P;
    P.set(c1, c1, src, q, 1.0);
    P.set(c1, c1, q, sk, 1.0);
    P.set(c2, c2, src, q, 1.0);
    P.set(c2, c2, q, sk, 1.0);
    m.link(P);
    m.set_capacity(q, 6);

    ctmc::CtmcOptions opt;
    opt.cutoff = 6;
    const line::mva::AvgResult<double> r = ctmc::solver_ctmc_run_analyzer(m.get_struct(), opt);

    // MATLAB SolverCTMC on the same model and the same cutoff.
    CHECK(r.QN(1, 0) == doctest::Approx(0.1344178322).epsilon(1e-6));
    CHECK(r.QN(1, 1) == doctest::Approx(0.1984882290).epsilon(1e-6));
    CHECK(r.UN(1, 0) == doctest::Approx(0.0999816883).epsilon(1e-6));
    CHECK(r.UN(1, 1) == doctest::Approx(0.1499725325).epsilon(1e-6));
    CHECK(r.TN(1, 0) == doctest::Approx(0.1999633767).epsilon(1e-6));
    CHECK(r.TN(1, 1) == doctest::Approx(0.2999450650).epsilon(1e-6));

    // The same numbers, read as physics rather than as a reference table: no
    // work is lost except at the finite buffer, so each class's throughput is
    // its arrival rate less the loss (well under a per mille at rho = 0.25),
    // and utilization is the carried load. This holds for ANY polling order, so
    // it survives a change of discipline that would invalidate the table above.
    CHECK(r.TN(1, 0) == doctest::Approx(l1).epsilon(2e-3));
    CHECK(r.TN(1, 1) == doctest::Approx(l2).epsilon(2e-3));
    CHECK(r.UN(1, 0) + r.UN(1, 1) == doctest::Approx((l1 + l2) / mu).epsilon(2e-3));
}

TEST_CASE("a cache splits its arrival stream into hits and misses end to end") {
    // Source -> Cache -> Sink, two items and room for one. Every read leaves as
    // either a hit or a miss, so the two departure classes must carry the whole
    // arrival rate between them -- an identity independent of the replacement
    // policy, and the one a mis-wired afterEventCache breaks first.
    const double lambda = 1.0;
    qn::Network<double> m("cacheqn");
    const std::size_t src = m.add_source("Src");
    const std::size_t sk = m.add_sink("Sink");
    const std::size_t rd = m.add_open_class("Read");
    const std::size_t hit = m.add_open_class("Hit");
    const std::size_t mis = m.add_open_class("Miss");

    qn::CacheParam<double> cp;
    cp.nitems = 2;
    cp.itemcap.push_back(1);
    cp.replacestrat = line::lang::ReplacementStrategy::LRU;
    cp.pread.assign(3, std::vector<double>());
    cp.pread[rd - 1] = std::vector<double>{0.5, 0.5};
    cp.hitclass.assign(3, 0);
    cp.missclass.assign(3, 0);
    cp.hitclass[rd - 1] = hit;
    cp.missclass[rd - 1] = mis;
    const std::size_t ca = m.add_cache("C", cp);

    m.set_arrival(src, rd, Dist::exp_rate(lambda));
    qn::RoutingMatrix<double> P;
    P.set(rd, rd, src, ca, 1.0);
    P.set(hit, hit, ca, sk, 1.0);
    P.set(mis, mis, ca, sk, 1.0);
    m.link(P);

    ctmc::CtmcOptions opt;
    opt.cutoff = 1;
    const ctmc::CtmcSolution<double> d =
        ctmc::solver_ctmc_analyzer(m.get_struct(), opt);

    // The stationary probability mass is normalized, which is the invariant the
    // reducible-component selection has to preserve on a model whose encoding
    // space is wider than its dynamics -- a cache's item lists are exactly that.
    double tot = 0;
    for (std::size_t s = 0; s < d.pi.size(); ++s) tot += d.pi[s];
    CHECK(tot == doctest::Approx(1.0).epsilon(1e-9));

    // Departure rates at the cache, per class, weighted by the stationary law.
    // Every read leaves as exactly one of hit or miss, so the two together carry
    // the read arrival rate; a cache handler that dropped or duplicated a read
    // shows up here and nowhere in a per-state unit test.
    const std::size_t isf = m.get_struct().stateful_index(ca);
    REQUIRE(isf != 0);
    double dep_hit = 0, dep_miss = 0, arv_read = 0;
    for (std::size_t s = 0; s < d.chain.space.size(); ++s) {
        dep_hit += d.pi[s] * d.chain.dep_rates[s][isf - 1][hit - 1];
        dep_miss += d.pi[s] * d.chain.dep_rates[s][isf - 1][mis - 1];
        arv_read += d.pi[s] * d.chain.arv_rates[s][isf - 1][rd - 1];
    }
    CHECK(dep_hit + dep_miss == doctest::Approx(arv_read).epsilon(1e-9));
    // MATLAB's `Cache.getHitRatio` / `getMissRatio` on the identical model both
    // report 0.5: with two equally popular items and room for one, a read hits
    // exactly when it repeats the previous one. At lambda = 1 that makes each
    // departure stream carry half the arrival rate.
    CHECK(dep_hit == doctest::Approx(0.5 * lambda).epsilon(1e-6));
    CHECK(dep_miss == doctest::Approx(0.5 * lambda).epsilon(1e-6));
}
