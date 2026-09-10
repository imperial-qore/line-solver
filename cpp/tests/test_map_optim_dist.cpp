/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * map_optim_dist / map_optim_dist_acf: fitting a MAP's D1 by minimizing a
 * distance to a reference MAP, with D0 held fixed.
 *
 * THE ORACLE IS A FIXED POINT. Fit the reference against ITSELF -- B0 = A0 and
 * alB = alA -- and the answer must be B1 = A1 with distance zero. That is an
 * exact identity, not a tolerance: it holds only if the constraint block, the
 * six Lyapunov solves, the Kronecker ordering and the optimizer are all right
 * together. Measured, both fitters recover A1 to machine precision.
 *
 * THE VECTOR alA IS THE EMBEDDED AT-ARRIVALS LAW, `map_pie`, NOT the phase
 * stationary law `map_prob`. The first constraint block IS the statement that
 * alB is stationary at arrivals, so passing `map_prob` makes the reference's own
 * D1 infeasible -- measured: a constraint residual of 2.9e-1 against 1.1e-16
 * with the right vector. That trap is what this file's first test pins.
 */
#include <cmath>
#include <cstddef>
#include <vector>

#include "doctest.h"
#include "line/api/mam/map_moment.h"
#include "line/api/mam/map_optim_dist.h"

namespace mam = line::mam;
using line::Matrix;

namespace {

/** A two-phase MAP with correlated arrivals. */
mam::Map<double> ref_map() {
    mam::Map<double> a;
    a.D0 = Matrix<double>(2, 2, 0.0);
    a.D1 = Matrix<double>(2, 2, 0.0);
    a.D0(0, 0) = -2.0;
    a.D0(0, 1) = 0.6;
    a.D0(1, 1) = -1.0;
    a.D1(0, 0) = 1.0;
    a.D1(0, 1) = 0.4;
    a.D1(1, 0) = 0.7;
    a.D1(1, 1) = 0.3;
    return a;
}

}  // namespace

TEST_CASE("the constraints are satisfied by the reference's own D1") {
    const mam::Map<double> a = ref_map();
    const std::vector<double> pie = mam::map_pie(a);

    Matrix<double> Aeq;
    std::vector<double> beq;
    mam::optdistdetail::build_constraints(a.D0, pie, &Aeq, &beq);
    const std::vector<double> v = mam::optdistdetail::vec(a.D1);
    REQUIRE(Aeq.cols() == v.size());

    double worst = 0.0;
    for (std::size_t i = 0; i < Aeq.rows(); ++i) {
        double s = -beq[i];
        for (std::size_t j = 0; j < v.size(); ++j) s += Aeq(i, j) * v[j];
        worst = std::max(worst, std::fabs(s));
    }
    CHECK(worst < 1e-12);

    // The WRONG law makes the reference's own D1 infeasible, which is the trap
    // this pins: map_prob is the phase stationary law, map_pie the embedded
    // at-arrivals one, and only the latter is what the constraint states.
    Matrix<double> Aeq2;
    std::vector<double> beq2;
    mam::optdistdetail::build_constraints(a.D0, mam::map_prob(a), &Aeq2, &beq2);
    double worst2 = 0.0;
    for (std::size_t i = 0; i < Aeq2.rows(); ++i) {
        double s = -beq2[i];
        for (std::size_t j = 0; j < v.size(); ++j) s += Aeq2(i, j) * v[j];
        worst2 = std::max(worst2, std::fabs(s));
    }
    CHECK(worst2 > 1e-3);
}

TEST_CASE("fitting a MAP against itself recovers its own D1, at distance zero") {
    const mam::Map<double> a = ref_map();
    const std::vector<double> pie = mam::map_pie(a);

    const mam::MapOptimDist<double> r = mam::map_optim_dist(a, pie, a.D0, pie, 1u);
    CHECK(r.global == true);  // the lag-1 objective is a convex quadratic
    CHECK(std::fabs(r.d) < 1e-9);
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j)
            CHECK(r.B1(i, j) == doctest::Approx(a.D1(i, j)).epsilon(1e-6));
}

TEST_CASE("the autocorrelation fit recovers it too") {
    const mam::Map<double> a = ref_map();
    const std::vector<double> pie = mam::map_pie(a);

    const mam::MapOptimDist<double> r = mam::map_optim_dist_acf(a, pie, a.D0, pie);
    // A general nonlinear program, so no global claim is made.
    CHECK(r.global == false);
    CHECK(std::fabs(r.d) < 1e-9);
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j)
            CHECK(r.B1(i, j) == doctest::Approx(a.D1(i, j)).epsilon(1e-5));
}

TEST_CASE("the fitted D1 is a valid MAP half, whatever the objective") {
    const mam::Map<double> a = ref_map();
    const std::vector<double> pie = mam::map_pie(a);

    // A DIFFERENT D0, so the fit is not the identity and the constraints are
    // doing real work.
    Matrix<double> B0(2, 2, 0.0);
    B0(0, 0) = -1.5;
    B0(0, 1) = 0.4;
    B0(1, 1) = -0.9;
    std::vector<double> alB(2, 0.0);
    alB[0] = 0.6;
    alB[1] = 0.4;

    const mam::MapOptimDist<double> r = mam::map_optim_dist(a, pie, B0, alB, 1u);
    for (std::size_t i = 0; i < 2; ++i) {
        // Every entry is at or above the reference's own 1e-6 floor, which is
        // what keeps the fitted MAP irreducible.
        for (std::size_t j = 0; j < 2; ++j) CHECK(r.B1(i, j) >= 1e-6 - 1e-12);
        // The row sums of B0 + B1 vanish, i.e. it is a generator.
        double s = 0.0;
        for (std::size_t j = 0; j < 2; ++j) s += B0(i, j) + r.B1(i, j);
        CHECK(std::fabs(s) < 1e-6);
    }
    CHECK(r.d >= -1e-12);
}

TEST_CASE("a lag count above one takes the nonlinear path and says so") {
    const mam::Map<double> a = ref_map();
    const std::vector<double> pie = mam::map_pie(a);
    const mam::MapOptimDist<double> r = mam::map_optim_dist(a, pie, a.D0, pie, 3u);
    // No global claim off the quadratic path: the reference calls fmincon here.
    CHECK(r.global == false);
    CHECK(r.d >= -1e-9);
    for (std::size_t i = 0; i < 2; ++i) {
        double s = 0.0;
        for (std::size_t j = 0; j < 2; ++j) s += a.D0(i, j) + r.B1(i, j);
        CHECK(std::fabs(s) < 1e-5);
    }
}

TEST_CASE("the refusals are by name") {
    const mam::Map<double> a = ref_map();
    const std::vector<double> pie = mam::map_pie(a);
    Matrix<double> B0(2, 2, 0.0);
    B0(0, 0) = -1.0;
    B0(1, 1) = -1.0;

    CHECK_THROWS_AS(mam::map_optim_dist(a, pie, Matrix<double>(2, 3, 0.0), pie, 1u),
                    line::InputError);
    CHECK_THROWS_AS(mam::map_optim_dist(a, pie, B0, std::vector<double>(3, 0.3), 1u),
                    line::InputError);
    CHECK_THROWS_AS(mam::map_optim_dist(a, std::vector<double>(5, 0.2), B0, pie, 1u),
                    line::InputError);
    CHECK_THROWS_AS(mam::map_optim_dist(a, pie, B0, pie, 0u), line::InputError);
    CHECK_THROWS_AS(mam::map_optim_dist_acf(a, pie, Matrix<double>(2, 3, 0.0), pie),
                    line::InputError);
}
