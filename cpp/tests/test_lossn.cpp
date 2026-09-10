/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Erlang fixed point and the damped iteration under it. Oracles: Erlang's loss
 * formula against its rational recursion and against hand values, the
 * single-link network where the fixed point is the loss formula itself, and
 * monotonicity of blocking in the offered load.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/lossn/lossn_erlangfp.h"

using line::Matrix;
using line::da::FpiOptions;
using line::lossn::erlang_b;
using line::lossn::lossn_erlangfp;

namespace {

/** Erlang B by the rational recursion B_k = nu B_{k-1} / (k + nu B_{k-1}). */
double erlang_b_recursive(double nu, int C) {
    double b = 1.0;
    for (int k = 1; k <= C; ++k) b = nu * b / (k + nu * b);
    return b;
}

}  // namespace

TEST_CASE("erlang_b agrees with the rational recursion") {
    for (double nu : {0.5, 1.0, 3.0, 7.5}) {
        for (int C : {1, 2, 5, 10}) {
            INFO("nu = ", nu, ", C = ", C);
            CHECK(erlang_b(nu, C) == doctest::Approx(erlang_b_recursive(nu, C)).epsilon(1e-12));
        }
    }
    // B(nu, 1) = nu/(1+nu) in closed form.
    CHECK(erlang_b(2.0, 1) == doctest::Approx(2.0 / 3.0).epsilon(1e-13));
}

TEST_CASE("erlang_b is increasing in load and decreasing in capacity") {
    CHECK(erlang_b(3.0, 5) > erlang_b(1.0, 5));
    CHECK(erlang_b(3.0, 10) < erlang_b(3.0, 5));
}

TEST_CASE("single-link loss network reproduces Erlang B") {
    // One link, one class, unit circuit requirement: the fixed point is
    // E = B(nu/(1-E), C) and the carried traffic is nu(1-E).
    Matrix<double> A(1, 1);
    A(0, 0) = 1.0;
    const std::vector<double> nu{4.0};
    const std::vector<int> C{6};
    auto r = lossn_erlangfp(nu, A, C);
    CHECK(r.converged);
    CHECK(r.E[0] > 0.0);
    CHECK(r.E[0] < 1.0);
    // The reduced load of link j already contains the factor (1-E_j)^A(j,r),
    // which the division by (1-E_j) cancels when A(j,r) = 1. So on a single
    // link with unit circuit requirement the fixed point collapses to Erlang B
    // at the raw offered load, with no reduction at all.
    CHECK(r.E[0] == doctest::Approx(erlang_b(nu[0], C[0])).epsilon(1e-7));
    // Carried traffic and loss are consistent.
    CHECK(r.QLen[0] == doctest::Approx(nu[0] * (1.0 - r.E[0])).epsilon(1e-12));
    CHECK(r.Loss[0] == doctest::Approx(r.E[0]).epsilon(1e-12));
}

TEST_CASE("two-link network: a class crossing both links loses more") {
    // Class 0 uses link 0 only, class 1 uses both links.
    Matrix<double> A(2, 2, 0.0);
    A(0, 0) = 1.0;
    A(0, 1) = 1.0;
    A(1, 1) = 1.0;
    const std::vector<double> nu{3.0, 3.0};
    const std::vector<int> C{5, 5};
    auto r = lossn_erlangfp(nu, A, C);
    CHECK(r.converged);
    CHECK(r.Loss[1] > r.Loss[0]);  // the two-link route is blocked by either link
    for (std::size_t j = 0; j < 2; ++j) {
        CHECK(r.E[j] > 0.0);
        CHECK(r.E[j] < 1.0);
    }
    for (std::size_t c = 0; c < 2; ++c) {
        CHECK(r.QLen[c] <= nu[c]);
        CHECK(r.Loss[c] >= 0.0);
    }
}

TEST_CASE("loss grows with offered load") {
    Matrix<double> A(1, 1);
    A(0, 0) = 1.0;
    const std::vector<int> C{4};
    auto light = lossn_erlangfp(std::vector<double>{1.0}, A, C);
    auto heavy = lossn_erlangfp(std::vector<double>{8.0}, A, C);
    CHECK(heavy.Loss[0] > light.Loss[0]);
}

TEST_CASE("da_fpi damping reaches the same fixed point more slowly") {
    Matrix<double> A(1, 1);
    A(0, 0) = 1.0;
    const std::vector<double> nu{4.0};
    const std::vector<int> C{6};

    FpiOptions damped;
    damped.damping = 0.5;
    auto undamped = lossn_erlangfp(nu, A, C);
    auto slow = lossn_erlangfp(nu, A, C, damped);
    CHECK(slow.converged);
    CHECK(slow.E[0] == doctest::Approx(undamped.E[0]).epsilon(1e-6));
    CHECK(slow.iterations >= undamped.iterations);
}
