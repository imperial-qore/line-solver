/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Fast adaptive uniformization, ctmc_fau.
 *
 * Oracles:
 *  - the closed-form transient of the two-state chain, which is exact;
 *  - ctmc_foxglynn at a tolerance three orders tighter, which is a genuinely
 *    different algorithm on the same generator: full-space uniformization at
 *    max_i |q_ii| against the adaptive rate sequence.
 *
 * Two structural properties are asserted alongside the numbers, because they
 * are what the method promises and what a plausible-looking wrong answer would
 * break: the result is a componentwise lower bound on the exact distribution,
 * and its missing mass IS its L1 error.
 *
 * Twin of python/tests/test_ctmc_fau.py and
 * jar/src/test/java/jline/api/mc/CtmcFauTest.java.
 */
#include <cmath>
#include <cstddef>
#include <vector>

#include "doctest.h"
#include "line/api/mc/ctmc_fau.h"
#include "line/api/mc/ctmc_foxglynn.h"
#include "line/util/matrix.h"

using line::Matrix;
using line::mc::ctmc_fau;
using line::mc::ctmc_foxglynn;

namespace {

Matrix<double> mm1k_generator(double lambda, double mu, std::size_t K) {
    const std::size_t n = K + 1;
    Matrix<double> Q(n, n, 0.0);
    for (std::size_t i = 0; i + 1 < n; ++i) {
        Q(i, i + 1) = lambda;
        Q(i + 1, i) = mu;
    }
    for (std::size_t i = 0; i < n; ++i) {
        double off = 0.0;
        if (i + 1 < n) off += lambda;
        if (i > 0) off += mu;
        Q(i, i) = -off;
    }
    return Q;
}

std::vector<double> unit_vector(std::size_t n, std::size_t i) {
    std::vector<double> v(n, 0.0);
    v[i] = 1.0;
    return v;
}

double l1(const std::vector<double>& a, const std::vector<double>& b) {
    double s = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) s += std::fabs(a[i] - b[i]);
    return s;
}

}  // namespace

TEST_CASE("ctmc_fau matches the closed-form transient of a two-state chain") {
    const double a = 3.0, b = 1.5, t = 0.7;
    Matrix<double> Q(2, 2, 0.0);
    Q(0, 0) = -a;
    Q(0, 1) = a;
    Q(1, 0) = b;
    Q(1, 1) = -b;

    const auto r = ctmc_fau(unit_vector(2, 0), Q, t);

    const double decay = std::exp(-(a + b) * t);
    CHECK(r.pit[0] == doctest::Approx(b / (a + b) + a / (a + b) * decay).epsilon(1e-6));
    CHECK(r.pit[1] == doctest::Approx(a / (a + b) * (1.0 - decay)).epsilon(1e-6));
    CHECK(r.steps > 1);
}

TEST_CASE("ctmc_fau reports its error as missing mass, and it is a lower bound") {
    const Matrix<double> Q = mm1k_generator(1.0, 2.0, 30);
    const std::vector<double> pi0 = unit_vector(31, 0);
    const auto r = ctmc_fau(pi0, Q, 5.0);
    const auto reference = ctmc_foxglynn(pi0, Q, 5.0, 1e-14, -1);

    // The three approximations only remove mass, so the defect IS the error.
    CHECK(r.errorBound == doctest::Approx(l1(r.pit, reference.pi)).epsilon(1e-6));
    CHECK(r.errorBound <= 1e-5);
    for (std::size_t i = 0; i < r.pit.size(); ++i) CHECK(r.pit[i] <= reference.pi[i] + 1e-12);
}

TEST_CASE("ctmc_fau tightens with the tolerance") {
    const Matrix<double> Q = mm1k_generator(1.0, 2.0, 20);
    const std::vector<double> pi0 = unit_vector(21, 0);
    const auto loose = ctmc_fau(pi0, Q, 4.0, 1e-4);
    const auto tight = ctmc_fau(pi0, Q, 4.0, 1e-10);
    CHECK(tight.errorBound < loose.errorBound);
    CHECK(tight.steps > loose.steps);
}

TEST_CASE("ctmc_fau does not pay for fast states that carry no mass") {
    // Six slow states, then four states of rate 1e6 that the initial
    // distribution cannot reach within the horizon. Ordinary uniformization
    // pays for the fast ones, adaptive uniformization does not.
    const std::size_t ns = 6, n = ns + 4;
    const double t = 1.0;
    Matrix<double> Q(n, n, 0.0);
    for (std::size_t i = 0; i + 1 < ns; ++i) {
        Q(i, i + 1) = 0.5;
        Q(i + 1, i) = 0.4;
    }
    for (std::size_t i = ns; i + 1 < n; ++i) {
        Q(i, i + 1) = 1e6;
        Q(i + 1, i) = 1e6;
    }
    Q(n - 1, ns) = 1e6;
    for (std::size_t i = 0; i < n; ++i) {
        double off = 0.0;
        for (std::size_t j = 0; j < n; ++j)
            if (j != i) off += Q(i, j);
        Q(i, i) = -off;
    }
    const std::vector<double> pi0 = unit_vector(n, 0);

    const auto r = ctmc_fau(pi0, Q, t);
    const auto reference = ctmc_foxglynn(pi0, Q, t, 1e-14, -1);

    CHECK(l1(r.pit, reference.pi) < 1e-6);
    CHECK(r.lambdaMax <= 1.0);
    CHECK(r.uniformRate >= 1e6);
    // The step count follows the visited rate, not the global one.
    CHECK(r.steps < 50);
}

TEST_CASE("ctmc_fau bounds the support by the occupancy threshold") {
    const std::size_t K = 400;
    const Matrix<double> Q = mm1k_generator(1.0, 3.0, K);
    const std::vector<double> pi0 = unit_vector(K + 1, 0);
    const auto r = ctmc_fau(pi0, Q, 4.0, 1e-8, 1e-10);
    const auto reference = ctmc_foxglynn(pi0, Q, 4.0, 1e-14, -1);

    CHECK(r.supportMax < K + 1);
    CHECK(l1(r.pit, reference.pi) <= r.errorBound + 1e-12);
    CHECK(l1(r.pit, reference.pi) < 1e-7);
}

TEST_CASE("ctmc_fau terminates on an absorbing chain") {
    Matrix<double> Q(3, 3, 0.0);
    Q(0, 0) = -2.0;
    Q(0, 1) = 2.0;
    Q(1, 1) = -1.0;
    Q(1, 2) = 1.0;
    const std::vector<double> pi0 = unit_vector(3, 0);
    const auto r = ctmc_fau(pi0, Q, 3.0);
    const auto reference = ctmc_foxglynn(pi0, Q, 3.0, 1e-14, -1);
    CHECK(l1(r.pit, reference.pi) < 1e-6);
    CHECK(r.absorbed);
    CHECK(r.steps == 3);
}

TEST_CASE("ctmc_fau returns the initial distribution at a zero horizon") {
    const Matrix<double> Q = mm1k_generator(1.0, 2.0, 5);
    const std::vector<double> pi0 = unit_vector(6, 2);
    const auto r = ctmc_fau(pi0, Q, 0.0);
    CHECK(r.pit == pi0);
    CHECK(r.errorBound == 0.0);
}

TEST_CASE("ctmc_fau refuses a non-square generator and a negative horizon") {
    Matrix<double> rect(2, 3, 0.0);
    CHECK_THROWS_AS(ctmc_fau(unit_vector(2, 0), rect, 1.0), line::InputError);
    const Matrix<double> Q = mm1k_generator(1.0, 2.0, 3);
    CHECK_THROWS_AS(ctmc_fau(unit_vector(4, 0), Q, -1.0), line::InputError);
}
