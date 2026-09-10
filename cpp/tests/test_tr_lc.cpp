/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Load concealment as a model TRANSFORMATION (`tr::transform_solve_lc`).
 *
 * Birman-Kogan Algorithm 2 exists twice in the tree and the two are not the
 * same computation. `pfqn_bklc` is the KERNEL: it sweeps chains on a demand
 * matrix with MVA as the inner single-chain solve. `transform_solve_lc` is the
 * TRANSFORMATION: the same sweep, but each single-chain subproblem is a real
 * single-class struct solved by whichever analyzer the caller bound.
 *
 * THE ORACLE IS THE KERNEL'S OWN FIXED POINT. On a model whose per-chain
 * subproblem is single-class, PS and exponential, the chain aggregation is
 * exact, so the two must reach the SAME fixed point in the SAME number of
 * sweeps. Agreement on the fixed point alone is not enough: a Jacobi sweep
 * reaches the same point at a different sweep count, and the sweep count is
 * what the four codebases are pinned on.
 *
 * The MATLAB twin is line-test.git test_tr_lc.m, the python twin is
 * python/tests/test_tr_lc.py and the JAR twin is LcStrategyTest.java.
 */
#include "doctest.h"
#include "line/api/pfqn/pfqn_bk.h"
#include "line/lang/qn/network_builder.h"
#include "line/solvers/ctmc/solver_ctmc_analyzer.h"
#include "line/solvers/tr/transform_solve.h"

namespace {

using DD = line::lang::Distrib<double>;
using Net = line::qn::Network<double>;
using Routing = line::qn::RoutingMatrix<double>;
using line::lang::SchedStrategy;

/** Delay -> Q1 -> Q2 -> Q3 -> Delay, two closed classes, one chain each. */
Net build() {
    Net m("lc");
    const std::size_t d = m.add_delay("Think");
    const std::size_t q1 = m.add_queue("Q1", SchedStrategy::PS);
    const std::size_t q2 = m.add_queue("Q2", SchedStrategy::PS);
    const std::size_t q3 = m.add_queue("Q3", SchedStrategy::PS);
    const std::size_t c1 = m.add_closed_class("C1", 4.0, d);
    const std::size_t c2 = m.add_closed_class("C2", 3.0, d);
    m.set_service(d, c1, DD::exp_rate(1.0));
    m.set_service(d, c2, DD::exp_rate(2.0));
    m.set_service(q1, c1, DD::exp_rate(2.0));
    m.set_service(q1, c2, DD::exp_rate(3.0));
    m.set_service(q2, c1, DD::exp_rate(3.0));
    m.set_service(q2, c2, DD::exp_rate(1.5));
    m.set_service(q3, c1, DD::exp_rate(1.8));
    m.set_service(q3, c2, DD::exp_rate(2.5));
    Routing P;
    P.set(c1, c1, d, q1, 1.0);
    P.set(c1, c1, q1, q2, 1.0);
    P.set(c1, c1, q2, q3, 1.0);
    P.set(c1, c1, q3, d, 1.0);
    P.set(c2, c2, d, q1, 1.0);
    P.set(c2, c2, q1, q2, 1.0);
    P.set(c2, c2, q2, q3, 1.0);
    P.set(c2, c2, q3, d, 1.0);
    m.link(P);
    return m;
}

/** The demand matrix, populations and think times the model reduces to. */
line::Matrix<double> demands() {
    line::Matrix<double> L(4, 2, 0.0);
    L(1, 0) = 0.5;          L(1, 1) = 1.0 / 3;
    L(2, 0) = 1.0 / 3;      L(2, 1) = 2.0 / 3;
    L(3, 0) = 5.0 / 9;      L(3, 1) = 0.4;
    return L;
}

line::mva::AvgResult<double> elevated() {
    Net m = build();
    auto inner = [](const line::qn::NetworkStruct<double>& s) {
        line::ctmc::CtmcOptions o;
        return line::ctmc::solver_ctmc_run_analyzer(s, o);
    };
    return line::tr::transform_solve_lc<double>(m.get_struct(), inner, "default", 1000);
}

}  // namespace

TEST_CASE("the transformation reaches the kernel's own fixed point") {
    const std::vector<double> N = {4.0, 3.0}, Z = {1.0, 0.5};
    const line::pfqn::BkLcResult<double> k =
        line::pfqn::pfqn_bklc(demands(), N, Z, "mva", 1e-10, 1000);
    const line::mva::AvgResult<double> t = elevated();
    REQUIRE(t.XN.size() == k.X.size());
    for (std::size_t r = 0; r < k.X.size(); ++r) {
        CHECK(t.XN[r] == doctest::Approx(k.X[r]).epsilon(1e-8));
    }
}

TEST_CASE("the sweep is Gauss-Seidel, not Jacobi") {
    // Same fixed point in the same number of sweeps. Jacobi coupling reaches the
    // same point at a DIFFERENT sweep count, which is what this catches.
    const std::vector<double> N = {4.0, 3.0}, Z = {1.0, 0.5};
    const line::pfqn::BkLcResult<double> k =
        line::pfqn::pfqn_bklc(demands(), N, Z, "mva", 1e-10, 1000);
    CHECK(elevated().iter == static_cast<int>(k.it));
}

TEST_CASE("the transformed solve names itself") {
    CHECK(elevated().actualmethod == "default/lc");
}

TEST_CASE("the CtmcOptions flag reaches the same answer as the direct call") {
    // `transform_solve_lc` is reachable from the analyzer, not only from a test:
    // CtmcOptions::load_concealment is the C++ counterpart of the method name the
    // other three codebases carry in options.config.
    Net m = build();
    line::ctmc::CtmcOptions o;
    o.load_concealment = true;
    const line::mva::AvgResult<double> viaOption =
        line::ctmc::solver_ctmc_run_analyzer(m.get_struct(), o);
    const line::mva::AvgResult<double> direct = elevated();
    REQUIRE(viaOption.XN.size() == direct.XN.size());
    for (std::size_t r = 0; r < direct.XN.size(); ++r)
        CHECK(viaOption.XN[r] == doctest::Approx(direct.XN[r]).epsilon(1e-12));
    CHECK(viaOption.iter == direct.iter);
}
