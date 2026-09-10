/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The fluid solver's fork-join arm (the `has_fork` branch of
 * solvers/fluid/fluid_runner.h).
 *
 * SolverFluid used to refuse every Fork and Join outright, in all four
 * codebases. It now drives the SAME solver-agnostic fixed point MVA and NC
 * drive (`mva::fj_fixed_point`), with the fluid analyzer as the inner solve:
 * the MMT transform emits only Source, Delay, Queue, Router and ClassSwitch,
 * every one of which the drift already carries, so no fork-join code was added
 * to the fluid solver itself.
 *
 * THE ORACLE IS THE OTHER THREE CODEBASES, to the digit -- and as of 2026-09-08
 * that is FOUR again. The three-port pin this file used to carry was a
 * deliberate exception recorded as BUG-100: between 2026-08-31 and 2026-09-01
 * the ports answered 1.76630545 / 4.11063759 / 1.03650606 / 3.49158828 while
 * MATLAB @@SolverFLD answered 1.76038792 / 4.09961825 / 1.02862318 / 3.48364563,
 * and the ports were pinned rather than the reference. The joint share-capacity
 * closure of 2026-09-08 moved the ports onto the reference and closed BUG-100,
 * so the pin below is MATLAB's row and every codebase meets it:
 *
 *   MATLAB 1.760390673 / 4.099634162 / 1.028637721 / 3.483656026
 *   JAR    1.760389779 / 4.099637323 / 1.028636563 / 3.483658111
 *   python 1.760389736 / 4.099637475 / 1.028636507 / 3.483658210
 *   C++    1.76039     / 4.09964     / 1.02864     / 3.48365
 *
 * All four sit within 1.2e-6 relative of each other, which is ODE-solver
 * scatter and an order below the 1e-5 asserted here.
 *
 * What is checked besides the golden is what the fixed point is responsible
 * for: flow balance across the fork (every station on the cycle carries the
 * reference station's throughput once the auxiliary classes are merged back)
 * and a synchronisation delay actually charged at the join. A regression in the
 * auxiliary-class refresh shows up as a join whose response time is zero and a
 * delay carrying twice the queues' throughput, which is exactly the state the
 * three ports were in before it was fixed.
 */
#include <cmath>
#include <string>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/fluid/fluid_runner.h"
#include "line/solvers/mva/solver_mva_runner.h"

namespace {

using DD = line::lang::Distrib<double>;
using Net = line::qn::Network<double>;
using Routing = line::qn::RoutingMatrix<double>;
using line::lang::SchedStrategy;

/** Delay -> Fork -> {Q1, Q2} -> Join -> Delay, the model the four codebases share. */
Net closed_fj(double n = 6.0) {
    Net m("fjclosed");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::PS);
    const std::size_t f = m.add_fork("Fork");
    const std::size_t j = m.add_join("Join", f);
    const std::size_t c = m.add_closed_class("C1", n, d);
    m.set_service(d, c, DD::exp_rate(1.0));
    m.set_service(q1, c, DD::exp_rate(2.0));
    m.set_service(q2, c, DD::exp_rate(3.0));
    Routing P;
    P.set(c, c, d, f, 1.0);
    P.set(c, c, f, q1, 1.0);
    P.set(c, c, f, q2, 1.0);
    P.set(c, c, q1, j, 1.0);
    P.set(c, c, q2, j, 1.0);
    P.set(c, c, j, d, 1.0);
    m.link(P);
    return m;
}

line::fluid::FluidSolution fld(Net& m, const std::string& method = "default") {
    line::fluid::FluidOptions o;
    o.method = method;
    return line::fluid::solver_fluid_run_analyzer(m.get_struct(), o);
}

}  // namespace

TEST_CASE("the fluid solver no longer refuses a fork-join model") {
    Net m = closed_fj();
    const line::fluid::FluidSolution s = fld(m);
    REQUIRE(s.QN.rows() == 4);
    for (std::size_t i = 0; i < s.QN.rows(); ++i) {
        CHECK(std::isfinite(s.QN(i, 0)));
        CHECK(s.QN(i, 0) >= 0.0);
    }
}

TEST_CASE("the fork-join fluid table matches the other three codebases") {
    Net m = closed_fj();
    const line::fluid::FluidSolution s = fld(m);
    // MATLAB @@SolverFLD on the same model, at default options, which is the
    // reference row every codebase now meets; see the file header for the four
    // measurements and for the BUG-100 pin this replaced. The tolerance is 1e-5
    // relative because the codebases integrate the drift with different ODE
    // solvers, which is the documented reason FLD is not expected to agree bit
    // for bit; everything else about the fixed point is shared, so a larger gap
    // than that is a port defect. All four sit within 1.2e-6 of these values.
    CHECK(s.QN(0, 0) == doctest::Approx(1.760390673).epsilon(1e-5));
    CHECK(s.QN(1, 0) == doctest::Approx(4.099634162).epsilon(1e-5));
    CHECK(s.QN(2, 0) == doctest::Approx(1.028637721).epsilon(1e-5));
    CHECK(s.QN(3, 0) == doctest::Approx(3.483656026).epsilon(1e-5));
}

TEST_CASE("flow is conserved across the fork") {
    Net m = closed_fj();
    const line::fluid::FluidSolution s = fld(m);
    const double ref = s.TN(0, 0);
    REQUIRE(ref > 0.0);
    for (std::size_t i = 1; i < s.TN.rows(); ++i)
        CHECK(s.TN(i, 0) == doctest::Approx(ref).epsilon(1e-5));
}

TEST_CASE("the join charges a synchronisation delay") {
    Net m = closed_fj();
    const line::fluid::FluidSolution s = fld(m);
    // Zero here is the signature of an auxiliary class that never carried
    // traffic, i.e. a rate written to sn.rates but not to the phase
    // representation the drift reads.
    CHECK(s.RN(3, 0) > 1e-3);
    CHECK(s.QN(3, 0) > 1e-3);
}

TEST_CASE("fluid and MVA solve the same transformed model to within their own gap") {
    Net m = closed_fj();
    const line::fluid::FluidSolution sf = fld(m);
    Net mm = closed_fj();
    line::mva::MvaOptions opt;
    opt.method = "amva";
    line::Matrix<double> init;
    const line::mva::AvgResult<double> sm = line::mva::solver_mva_run_analyzer(mm.get_struct(), opt, init);
    const double rel = std::abs(sf.TN(0, 0) - sm.TN(0, 0)) / sm.TN(0, 0);
    CHECK(rel < 0.15);
}

TEST_CASE("a method that cannot take a mixed model still refuses a fork") {
    Net m = closed_fj();
    // The transform rides the parallelism on auxiliary OPEN classes, so a
    // closed-only method cannot run the fixed point; the featset says so by
    // name rather than the transform failing somewhere inside.
    CHECK_THROWS(fld(m, "tbi"));
    CHECK_THROWS(fld(m, "statedep"));
}
