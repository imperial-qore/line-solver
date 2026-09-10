/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * SolverCTMC and SolverFLD on a non-Markovian service law, i.e. the wiring of
 * `sn_nonmarkov_toph` into both analyzers.
 *
 * THE PH FIT IS NOT A PREFERENCE HERE, IT IS A CORRECTNESS CONDITION, and this
 * file is where that is pinned. MATLAB's CTMC carries a rational generator and
 * takes the two-moment CME by default; this port assembles an ordinary
 * generator, and an ME's D0 has off-diagonal entries that are not rates. Fed a
 * CME the chain is not a Markov chain and the answers are not probabilities:
 * measured on the model below, Util comes back 1.075 for a SINGLE server and
 * the throughput is 0.538 at the delay against 0.400 at the queue, which cannot
 * both be true around a cycle. The tests therefore assert the two properties an
 * invalid fit violates -- utilization at most one, and throughput balanced --
 * rather than only comparing a number.
 *
 * The MATLAB values below are from R2025a on the same model.
 */
#include <cmath>
#include <cstddef>
#include <vector>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/fluid/solver_fluid.h"

namespace qn = line::qn;
using line::lang::SchedStrategy;
using Dist = line::lang::Distrib<double>;

namespace {

/** Delay -> FCFS Queue -> Delay, two closed jobs, Gamma(4, 0.5) service. */
qn::Network<double> gamma_cycle(double njobs) {
    qn::Network<double> m("G");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", njobs, d);
    m.set_service(d, c, Dist::exp_rate(1.0));
    m.set_service(q, c, Dist::gamma_dist(4.0, 0.5));  // mean 2, SCV 1/4
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);
    return m;
}

}  // namespace

TEST_CASE("SolverCTMC solves a Gamma service instead of answering the raw law") {
    qn::Network<double> m = gamma_cycle(2.0);
    line::ctmc::CtmcOptions o;
    const line::ctmc::CtmcSolution<double> s =
        line::ctmc::solver_ctmc_analyzer(m.get_struct(), o);

    // The two properties an invalid surrogate breaks.
    CHECK(s.avg.UN(1, 0) <= 1.0 + 1e-9);
    CHECK(s.avg.TN(0, 0) == doctest::Approx(s.avg.TN(1, 0)).epsilon(1e-6));

    // The population is conserved.
    CHECK(s.avg.QN(0, 0) + s.avg.QN(1, 0) == doctest::Approx(2.0).epsilon(1e-6));

    // Utilization is the throughput times the mean service time.
    CHECK(s.avg.UN(1, 0) == doctest::Approx(2.0 * s.avg.TN(1, 0)).epsilon(1e-6));

    // MATLAB R2025a on this model: QLen 1.523, Util 0.95395, Tput 0.47698. The
    // gap is the fit, not a defect: MATLAB matches SCV 0.25 exactly with its
    // matrix exponential where the Bernstein phase-type lands at 0.323.
    CHECK(s.avg.QN(1, 0) == doctest::Approx(1.523).epsilon(0.01));
    CHECK(s.avg.UN(1, 0) == doctest::Approx(0.95395).epsilon(0.01));
    CHECK(s.avg.TN(1, 0) == doctest::Approx(0.47698).epsilon(0.01));
}

TEST_CASE("the CTMC phase budget is honoured and bounds the state space") {
    // A smaller budget must give a smaller chain and still a valid answer; this
    // is what makes a large model solvable at all.
    qn::Network<double> m = gamma_cycle(1.0);
    line::ctmc::CtmcOptions o;
    o.nonmkv_order = 4;
    const line::ctmc::CtmcSolution<double> s =
        line::ctmc::solver_ctmc_analyzer(m.get_struct(), o);
    CHECK(s.avg.UN(1, 0) <= 1.0 + 1e-9);
    CHECK(s.avg.TN(0, 0) == doctest::Approx(s.avg.TN(1, 0)).epsilon(1e-6));
    CHECK(s.avg.QN(0, 0) + s.avg.QN(1, 0) == doctest::Approx(1.0).epsilon(1e-6));

    line::ctmc::CtmcOptions o20;
    const line::ctmc::CtmcSolution<double> s20 =
        line::ctmc::solver_ctmc_analyzer(m.get_struct(), o20);
    CHECK(s20.chain.space.size() > s.chain.space.size());
}

TEST_CASE("SolverFLD solves a Gamma service and conserves the population") {
    qn::Network<double> m = gamma_cycle(2.0);
    line::fluid::FluidOptions o;
    const line::fluid::FluidSolution s = line::fluid::solver_fluid(m.get_struct(), o);

    CHECK(s.UN(1, 0) <= 1.0 + 1e-6);
    CHECK(s.QN(0, 0) + s.QN(1, 0) == doctest::Approx(2.0).epsilon(1e-3));
    CHECK(s.TN(0, 0) == doctest::Approx(s.TN(1, 0)).epsilon(1e-3));
    // A fluid limit at two jobs is an approximation, so this is a sanity band
    // around the exact CTMC answer rather than an agreement claim.
    CHECK(s.QN(1, 0) > 0.5);
    CHECK(s.QN(1, 0) < 2.0);
}

TEST_CASE("a Markovian model is untouched by the conversion") {
    // Exponential service everywhere: the conversion must be a no-op, so the
    // answer is the ordinary closed-form one for a two-station cycle.
    qn::Network<double> m("E");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q = m.add_queue("Q", SchedStrategy::FCFS);
    const std::size_t c = m.add_closed_class("C1", 1.0, d);
    m.set_service(d, c, Dist::exp_rate(1.0));
    m.set_service(q, c, Dist::exp_rate(0.5));
    qn::RoutingMatrix<double> P;
    P.set(d, q, 1.0);
    P.set(q, d, 1.0);
    m.link(P);

    line::ctmc::CtmcOptions o;
    const line::ctmc::CtmcSolution<double> s =
        line::ctmc::solver_ctmc_analyzer(m.get_struct(), o);
    // One job alternating between two exponential stages of means 1 and 2:
    // it sits at the queue two thirds of the time.
    CHECK(s.avg.QN(1, 0) == doctest::Approx(2.0 / 3.0).epsilon(1e-9));
    CHECK(s.avg.QN(0, 0) == doctest::Approx(1.0 / 3.0).epsilon(1e-9));
    CHECK(s.chain.space.size() == 2);  // no extra phases were introduced
}
