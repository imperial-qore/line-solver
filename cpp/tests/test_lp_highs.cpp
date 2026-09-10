/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The HiGHS sparse LP backend, against the exact dense tableau.
 *
 * The point of this file is the AGREEMENT check. util/simplex.h is exact under
 * Rational and is what every mapqn golden was taken with; lp_highs.h exists
 * only to reach models the dense tableau cannot. If the two ever disagree on a
 * model small enough for both, the goldens are no longer trustworthy, so that
 * comparison is asserted here on the same shapes the mapqn tests use.
 */
#include <vector>

#include "doctest.h"
#include "line/api/mapqn/mapqn_qr_bounds_rsrd.h"
#include "line/num/number.h"
#include "line/util/lp_highs.h"
#include "line/util/matrix.h"
#include "line/util/simplex.h"

using line::Matrix;
using line::lp::highs_available;
using line::lp::LpModel;
using line::lp::LpSolution;
using line::lp::LpStatus;
using line::lp::simplex_solve;

TEST_CASE("the backend is compiled in") { CHECK(highs_available()); }

#ifdef LINE_MP_HAVE_HIGHS
TEST_CASE("HiGHS reproduces the dense tableau on a small LP") {
    // max x0 + x1 s.t. x0 + 2 x1 <= 1, x >= 0. Optimum 1 at (1,0).
    LpModel<double> m(2);
    m.set_bounds(0, 0.0, 1e30);
    m.set_bounds(1, 0.0, 1e30);
    m.set_free_upper(0);
    m.set_free_upper(1);
    m.row_add(0, 1.0);
    m.row_add(1, 2.0);
    m.emit_le(1.0);
    m.set_cost(0, 1.0);
    m.set_cost(1, 1.0);
    m.set_maximize(true);

    const LpSolution<double> dense = simplex_solve(m);
    const LpSolution<double> sparse = line::lp::highs_solve(m);
    REQUIRE(dense.ok());
    REQUIRE(sparse.ok());
    CHECK(sparse.objective == doctest::Approx(dense.objective).epsilon(1e-9));
    CHECK(sparse.objective == doctest::Approx(1.0).epsilon(1e-9));
}

TEST_CASE("HiGHS handles equalities, ranges and a minimize sense") {
    // min x0 + x1 s.t. x0 + x1 = 3, x0 <= 2. Optimum 3.
    LpModel<double> m(2);
    m.set_bounds(0, 0.0, 2.0);
    m.set_bounds(1, 0.0, 10.0);
    m.row_add(0, 1.0);
    m.row_add(1, 1.0);
    m.emit_eq(3.0);
    m.set_cost(0, 1.0);
    m.set_cost(1, 1.0);
    m.set_maximize(false);

    const LpSolution<double> dense = simplex_solve(m);
    const LpSolution<double> sparse = line::lp::highs_solve(m);
    REQUIRE(dense.ok());
    REQUIRE(sparse.ok());
    CHECK(sparse.objective == doctest::Approx(dense.objective).epsilon(1e-9));
    CHECK(sparse.objective == doctest::Approx(3.0).epsilon(1e-9));
}

TEST_CASE("HiGHS detects infeasibility the dense tableau also rejects") {
    LpModel<double> m(1);
    m.set_bounds(0, 0.0, 1.0);
    m.row_add(0, 1.0);
    m.emit_eq(5.0);  // x = 5 with x <= 1
    const LpSolution<double> dense = simplex_solve(m);
    const LpSolution<double> sparse = line::lp::highs_solve(m);
    CHECK_FALSE(dense.ok());
    CHECK_FALSE(sparse.ok());
    CHECK(sparse.status == LpStatus::Infeasible);
}

TEST_CASE("on the RS-RD model both backends give the product-form answer") {
    // The same collapse instance test_mapqn_qr_bounds_rsrd.cpp uses, where the
    // exact answers are 3/7 and 6/7. This exercises a real QRF polytope --
    // redundant equalities, pinned variables and all -- rather than a toy.
    line::mapqn::QrRsrdParams<double> p;
    p.M = 2;
    p.N = 2;
    p.F.push_back(2);
    p.F.push_back(2);
    p.K.assign(2, 1);
    p.mu.push_back(Matrix<double>(1, 1, 1.0));
    p.mu.push_back(Matrix<double>(1, 1, 0.5));
    p.v.push_back(Matrix<double>(1, 1, 0.0));
    p.v.push_back(Matrix<double>(1, 1, 0.0));
    p.r = Matrix<double>(2, 2, 0.0);
    p.r(0, 1) = 1.0;
    p.r(1, 0) = 1.0;

    const line::mapqn::QrRsrdResult<double> a =
        line::mapqn::mapqn_qr_bounds_rsrd(p, 0, line::mapqn::MapqnSense::Min);
    REQUIRE(a.ok);
    CHECK(a.objective == doctest::Approx(3.0 / 7.0).epsilon(1e-9));
}
#endif
