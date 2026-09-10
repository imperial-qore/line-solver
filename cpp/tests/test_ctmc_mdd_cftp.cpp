/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The two SolverCTMC methods that never build the |S|-state generator: `mdd`
 * (Miner-Ciardo-Donatelli level aggregation over a decision diagram) and
 * `cftp` / `cftp.approx` (Kijima-Matsui perfect sampling).
 *
 * THE ORACLE IS THE SAME MODEL UNDER THE ENUMERATED CHAIN, which is the only
 * oracle that discriminates here: both methods return a complete, plausible
 * AvgTable for any closed model they accept, so a wrong descriptor or a wrong
 * estimator is invisible without a second route to the same numbers. On a
 * product-form model the paper's single approximation is an identity, so `mdd`
 * must agree with the explicit chain to solver tolerance, not merely closely.
 * `cftp` draws from the EXACT stationary law but reports sample averages, so it
 * is checked at a Monte Carlo tolerance and, more sharply, on the identities its
 * estimators are built to keep exactly -- flow balance and C = N/X -- which hold
 * per draw and so carry no sampling error at all.
 *
 * The refusals are tested too, because both methods are gated far more narrowly
 * than the solver: a model outside the gate would come back with the answer of a
 * different network rather than an error.
 */
#include <cmath>
#include <string>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/num/number.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/ctmc/solver_ctmc_cftp.h"
#include "line/solvers/ctmc/solver_ctmc_mdd_analyzer.h"

namespace qn = line::qn;
namespace ctmc = line::ctmc;
using line::lang::SchedStrategy;
using Dist = line::lang::Distrib<double>;

namespace {

/** Delay -> FCFS Queue -> Delay, one closed class of N jobs. */
qn::Network<double> cqn_delay_queue(double N, double think_rate, double serv_rate) {
    qn::Network<double> m("cqn");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", N, d);
    m.set_service(d, c, Dist::exp_rate(think_rate));
    m.set_service(q, c, Dist::exp_rate(serv_rate));
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    return m;
}

/** A cycle of three PS queues, one closed class of N jobs. */
qn::Network<double> cqn_three_ps(double N, double m1, double m2, double m3) {
    qn::Network<double> m("cqn3");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::PS);
    const std::size_t q3 = m.add_queue("Q3", SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("C1", N, q1);
    m.set_service(q1, c, Dist::exp_rate(m1));
    m.set_service(q2, c, Dist::exp_rate(m2));
    m.set_service(q3, c, Dist::exp_rate(m3));
    qn::RoutingMatrix<double> P;
    P.set(q1, q2, 1.0);
    P.set(q2, q3, 1.0);
    P.set(q3, q1, 1.0);
    m.link(P);
    return m;
}

}  // namespace

TEST_CASE("mdd reproduces the enumerated chain on a closed product-form network") {
    qn::Network<double> m = cqn_delay_queue(4.0, 1.0 / 2.0, 3.0);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    ctmc::CtmcOptions dflt;
    const line::mva::AvgResult<double> exact = ctmc::solver_ctmc_run_analyzer(sn, dflt);

    ctmc::CtmcOptions opt;
    opt.method = "mdd";
    const ctmc::CtmcMddSolution<double> s = ctmc::solver_ctmc_mdd_analyzer(sn, opt);

    for (std::size_t i = 0; i < sn.nstations; ++i) {
        CHECK(s.avg.QN(i, 0) == doctest::Approx(exact.QN(i, 0)).epsilon(1e-8));
        CHECK(s.avg.UN(i, 0) == doctest::Approx(exact.UN(i, 0)).epsilon(1e-8));
        CHECK(s.avg.TN(i, 0) == doctest::Approx(exact.TN(i, 0)).epsilon(1e-8));
        CHECK(s.avg.RN(i, 0) == doctest::Approx(exact.RN(i, 0)).epsilon(1e-8));
    }
    // The level chains are coupled only through rates, so nothing in the
    // iteration forces the marginals to describe the same population; the
    // conservation law is the test that they do.
    CHECK(s.avg.QN(0, 0) + s.avg.QN(1, 0) == doctest::Approx(4.0).epsilon(1e-6));
    CHECK(s.encoding == "np");
    CHECK(s.num_states == 5);  // n jobs at the queue, n = 0..4
    CHECK(s.iters >= 1);
}

TEST_CASE("mdd agrees with the enumerated chain on a three-station PS cycle") {
    // K = 3 is where the diagram compresses nothing, so this is the case the
    // method is NOT expected to save memory on; it must still be exact.
    qn::Network<double> m = cqn_three_ps(3.0, 2.0, 3.0, 5.0);
    const qn::NetworkStruct<double>& sn = m.get_struct();

    const line::mva::AvgResult<double> exact = ctmc::solver_ctmc_run_analyzer(sn, ctmc::CtmcOptions());
    ctmc::CtmcOptions opt;
    opt.method = "mdd";
    const ctmc::CtmcMddSolution<double> s = ctmc::solver_ctmc_mdd_analyzer(sn, opt);

    double total = 0.0;
    for (std::size_t i = 0; i < sn.nstations; ++i) {
        CHECK(s.avg.QN(i, 0) == doctest::Approx(exact.QN(i, 0)).epsilon(1e-7));
        CHECK(s.avg.TN(i, 0) == doctest::Approx(exact.TN(i, 0)).epsilon(1e-7));
        total += s.avg.QN(i, 0);
    }
    CHECK(total == doctest::Approx(3.0).epsilon(1e-6));
    CHECK(s.level_sizes.size() == 3);
}

TEST_CASE("mdd refuses the model classes its descriptor cannot encode") {
    // Open: the marking is unbounded, so the reachable set has no finite diagram.
    qn::Network<double> open("mm1");
    const std::size_t src = open.add_source("Src");
    const std::size_t q = open.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t k = open.add_sink("Sink");
    const std::size_t c = open.add_open_class("C1");
    open.set_arrival(src, c, Dist::exp_rate(0.5));
    open.set_service(q, c, Dist::exp_rate(1.0));
    qn::RoutingMatrix<double> Po;
    Po.set(src, q, 1.0);
    Po.set(q, k, 1.0);
    open.link(Po);
    ctmc::CtmcOptions opt;
    opt.method = "mdd";
    CHECK_THROWS_AS(ctmc::solver_ctmc_mdd_analyzer(open.get_struct(), opt),
                    line::UnsupportedError);

    // Multiclass: the Kronecker descriptor would need one level per
    // (station, class), which this encoding does not carry.
    qn::Network<double> two("cqn2");
    const std::size_t d = two.add_delay("Think");
    const std::size_t qq = two.add_queue("Q", SchedStrategy::PS);
    const std::size_t c1 = two.add_closed_class("A", 2.0, d);
    const std::size_t c2 = two.add_closed_class("B", 1.0, d);
    two.set_service(d, c1, Dist::exp_rate(1.0));
    two.set_service(d, c2, Dist::exp_rate(1.0));
    two.set_service(qq, c1, Dist::exp_rate(2.0));
    two.set_service(qq, c2, Dist::exp_rate(3.0));
    qn::RoutingMatrix<double> P2;
    P2.set(d, qq, 1.0);
    P2.set(qq, d, 1.0);
    two.link(P2);
    CHECK_THROWS_AS(ctmc::solver_ctmc_mdd_analyzer(two.get_struct(), opt),
                    line::UnsupportedError);
}

TEST_CASE("cftp estimates the enumerated chain and keeps its exact identities") {
    qn::Network<double> m = cqn_delay_queue(4.0, 1.0 / 2.0, 3.0);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const line::mva::AvgResult<double> exact = ctmc::solver_ctmc_run_analyzer(sn, ctmc::CtmcOptions());

    ctmc::CtmcOptions opt;
    opt.method = "cftp";
    ctmc::CtmcCftpOptions co;
    co.samples = 20000;
    co.seed = 23000;
    const ctmc::CtmcCftpSolution<double> s = ctmc::solver_ctmc_cftp(sn, opt, co);

    // A MONTE CARLO TOLERANCE, not a solver one: 2e4 draws put the standard
    // error of a mean queue length on a four-job model near 1e-2, so anything
    // tighter here would be a flaky test rather than a sharper one.
    for (std::size_t i = 0; i < sn.nstations; ++i)
        CHECK(s.avg.QN(i, 0) == doctest::Approx(exact.QN(i, 0)).epsilon(0.05));

    // These hold PER DRAW and so carry no sampling error: the population is
    // conserved by every sampled state, throughput is propagated from the
    // reference station through the visit ratios, and C = N/X is how CN is
    // formed. A drifting one of these is a wiring defect, not noise.
    CHECK(s.avg.QN(0, 0) + s.avg.QN(1, 0) == doctest::Approx(4.0).epsilon(1e-9));
    CHECK(s.avg.TN(0, 0) == doctest::Approx(s.avg.TN(1, 0)).epsilon(1e-9));
    CHECK(s.avg.CN[0] * s.avg.XN[0] == doctest::Approx(4.0).epsilon(1e-9));
    // Utilization keeps its own estimator precisely so it cannot leave [0,1].
    for (std::size_t i = 0; i < sn.nstations; ++i) {
        if (sn.stations[i].sched == SchedStrategy::INF) continue;
        CHECK(s.avg.UN(i, 0) >= 0.0);
        CHECK(s.avg.UN(i, 0) <= 1.0);
    }
    // The empirical law over the distinct draws is a distribution.
    double mass = 0.0;
    for (std::size_t j = 0; j < s.paggr.size(); ++j) mass += s.paggr[j];
    CHECK(mass == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(s.distinct_states.size() == s.paggr.size());
    CHECK(s.states.rows() == co.samples);
}

TEST_CASE("cftp.approx runs the rapidly-mixing sampler and lands on the same law") {
    qn::Network<double> m = cqn_three_ps(3.0, 2.0, 3.0, 5.0);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    const line::mva::AvgResult<double> exact = ctmc::solver_ctmc_run_analyzer(sn, ctmc::CtmcOptions());

    ctmc::CtmcOptions opt;
    opt.method = "cftp.approx";
    ctmc::CtmcCftpOptions co;
    co.samples = 20000;
    const ctmc::CtmcCftpSolution<double> s = ctmc::solver_ctmc_cftp(sn, opt, co);

    double total = 0.0;
    for (std::size_t i = 0; i < sn.nstations; ++i) {
        CHECK(s.avg.QN(i, 0) == doctest::Approx(exact.QN(i, 0)).epsilon(0.08));
        total += s.avg.QN(i, 0);
    }
    CHECK(total == doctest::Approx(3.0).epsilon(1e-9));
    CHECK(s.actualmethod == "cftp.approx");
}

TEST_CASE("cftp refuses everything outside the closed single-class product form") {
    ctmc::CtmcOptions opt;
    opt.method = "cftp";
    ctmc::CtmcCftpOptions co;
    co.samples = 100;

    // Non-exponential service: the balance function the sampler encodes is not
    // this model's, so the draws would be from a different network.
    qn::Network<double> erl("erl");
    const std::size_t d = erl.add_delay("Think");
    const std::size_t q = erl.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t c = erl.add_closed_class("C1", 2.0, d);
    erl.set_service(d, c, Dist::exp_rate(1.0));
    erl.set_service(q, c, Dist::erlang(4.0, 2));
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    erl.link(P);
    CHECK_THROWS_AS(ctmc::solver_ctmc_cftp(erl.get_struct(), opt, co), line::UnsupportedError);

    // A run length is required, not defaulted: the draw IS the answer here.
    qn::Network<double> ok = cqn_delay_queue(3.0, 1.0, 2.0);
    ctmc::CtmcCftpOptions none;
    none.samples = 0;
    CHECK_THROWS_AS(ctmc::solver_ctmc_cftp(ok.get_struct(), opt, none), line::InputError);

    // An unknown cftp variant is refused by name rather than falling back to
    // the exact sampler under a reported method that did not run.
    ctmc::CtmcOptions bad;
    bad.method = "cftp.perfect";
    CHECK_THROWS_AS(ctmc::solver_ctmc_cftp(ok.get_struct(), bad, co), line::InputError);
}

TEST_CASE("the generator-free methods are refused by the state-space analyzer") {
    // Reaching solver_ctmc_analyzer with one of these names means the caller
    // routed a generator-free method into the generator path: the answer would
    // be the enumerated one reported under a method that never ran.
    for (const char* name : {"mdd", "cftp", "cftp.approx"})
        CHECK_THROWS_AS(ctmc::check_method(name), line::UnsupportedError);
    // `exact` is an alias for the default state-space path and must pass.
    ctmc::check_method("default");
    ctmc::check_method("exact");
    // `gpu` PASSES the gate: the reference falls back to the CPU solve rather
    // than refusing, so this analyzer must run and report the fallback.
    ctmc::check_method("gpu");
    CHECK(ctmc::method_fallback_warning("gpu").find("Switching to default method") !=
          std::string::npos);
    CHECK_THROWS_AS(ctmc::check_method("nosuch"), line::UnsupportedError);

    const std::vector<std::string> valid = ctmc::list_valid_methods();
    CHECK(valid.size() == 6);
}

TEST_CASE("exact is behaviourally identical to default") {
    qn::Network<double> m = cqn_delay_queue(3.0, 1.0, 2.0);
    const qn::NetworkStruct<double>& sn = m.get_struct();
    ctmc::CtmcOptions dflt, ex;
    ex.method = "exact";
    const line::mva::AvgResult<double> a = ctmc::solver_ctmc_run_analyzer(sn, dflt);
    const line::mva::AvgResult<double> b = ctmc::solver_ctmc_run_analyzer(sn, ex);
    for (std::size_t i = 0; i < sn.nstations; ++i) {
        CHECK(a.QN(i, 0) == doctest::Approx(b.QN(i, 0)).epsilon(1e-14));
        CHECK(a.TN(i, 0) == doctest::Approx(b.TN(i, 0)).epsilon(1e-14));
    }
}
