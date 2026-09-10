/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * `CtmcOptions::fes_stations`, the first solver consumer of `fes_aggregate`.
 *
 * Flow-equivalent aggregation collapses a station subset into one
 * load-dependent station whose rate is the isolated subnetwork's throughput.
 * The transform has existed in all four codebases with NO solver consumer at
 * all: it was exercised by examples and tests only, so nothing in the solver
 * stack depended on it. This is that consumer, and the collapsed stations'
 * own metrics come back through the Chandy-Herzog-Woo conditional sum
 * E[Q_i] = sum_n P(N_fes = n) * Q_i(n).
 *
 * THE ORACLE IS AN IDENTITY, not a golden. On a product-form model the
 * decomposition is EXACT, so the reduced solve plus the conditioning must
 * reproduce the full chain's table station by station, INCLUDING the collapsed
 * stations. That is a statement about the transform, the reduced solve and the
 * back-mapping together: a wrong FES rate moves the surviving stations, a wrong
 * conditional sum moves only the collapsed ones, and a wrong visit ratio moves
 * only their throughput.
 */
#include <cmath>

#include "doctest.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"

namespace {

using DD = line::lang::Distrib<double>;
using Net = line::qn::Network<double>;
using Routing = line::qn::RoutingMatrix<double>;
using line::lang::SchedStrategy;

/** Think -> Q1 -> Q2 -> Q3 -> Think, a closed product-form cycle. */
Net cycle(double n = 4.0) {
    Net m("fes");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::PS);
    const std::size_t q3 = m.add_queue("Q3", SchedStrategy::PS);
    const std::size_t c = m.add_closed_class("C1", n, d);
    m.set_service(d, c, DD::exp_rate(1.0));
    m.set_service(q1, c, DD::exp_rate(2.0));
    m.set_service(q2, c, DD::exp_rate(3.0));
    m.set_service(q3, c, DD::exp_rate(4.0));
    Routing P;
    P.set(c, c, d, q1, 1.0);
    P.set(c, c, q1, q2, 1.0);
    P.set(c, c, q2, q3, 1.0);
    P.set(c, c, q3, d, 1.0);
    m.link(P);
    return m;
}

line::mva::AvgResult<double> solve(Net& m, const std::vector<std::size_t>& fes) {
    line::ctmc::CtmcOptions o;
    o.fes_stations = fes;
    return line::ctmc::solver_ctmc_run_analyzer(m.get_struct(), o);
}

}  // namespace

TEST_CASE("the FES-reduced solve reproduces the exact table, collapsed stations included") {
    Net a = cycle();
    Net b = cycle();
    const line::mva::AvgResult<double> exact = solve(a, {});
    const line::mva::AvgResult<double> fes = solve(b, {3, 4});

    REQUIRE(fes.QN.rows() == exact.QN.rows());
    for (std::size_t i = 0; i < exact.QN.rows(); ++i) {
        CHECK(fes.QN(i, 0) == doctest::Approx(exact.QN(i, 0)).epsilon(1e-9));
        CHECK(fes.UN(i, 0) == doctest::Approx(exact.UN(i, 0)).epsilon(1e-9));
        CHECK(fes.TN(i, 0) == doctest::Approx(exact.TN(i, 0)).epsilon(1e-9));
    }
}

TEST_CASE("the reduction names itself") {
    Net b = cycle();
    const line::mva::AvgResult<double> fes = solve(b, {3, 4});
    // A caller reading the table has to be able to tell that two of its rows
    // were recovered by conditioning rather than enumerated.
    CHECK(fes.actualmethod == "default/fes");
}

TEST_CASE("the collapsed subnetwork holds what the FES station held") {
    Net b = cycle();
    const line::mva::AvgResult<double> fes = solve(b, {3, 4});
    // The conditional sum splits the FES population between the two collapsed
    // stations; it must not create or destroy jobs.
    double total = 0.0;
    for (std::size_t i = 0; i < fes.QN.rows(); ++i) total += fes.QN(i, 0);
    CHECK(total == doctest::Approx(4.0).epsilon(1e-9));
}

TEST_CASE("a different subset gives the same answer") {
    // The choice of subset is the caller's and changes only which stations are
    // enumerated, so an exact decomposition must be invariant to it.
    Net a = cycle();
    Net b = cycle();
    Net c = cycle();
    const line::mva::AvgResult<double> exact = solve(a, {});
    const line::mva::AvgResult<double> fes23 = solve(b, {2, 3});
    const line::mva::AvgResult<double> fes34 = solve(c, {3, 4});
    for (std::size_t i = 0; i < exact.QN.rows(); ++i) {
        CHECK(fes23.QN(i, 0) == doctest::Approx(exact.QN(i, 0)).epsilon(1e-9));
        CHECK(fes34.QN(i, 0) == doctest::Approx(exact.QN(i, 0)).epsilon(1e-9));
    }
}

TEST_CASE("a subset that saves nothing is refused by name") {
    Net a = cycle();
    // One station is not a subnetwork.
    CHECK_THROWS_AS(solve(a, {3}), line::InputError);
    Net b = cycle();
    // Every station leaves no complement to solve.
    CHECK_THROWS_AS(solve(b, {1, 2, 3, 4}), line::InputError);
    Net c = cycle();
    CHECK_THROWS_AS(solve(c, {3, 99}), line::InputError);
}
