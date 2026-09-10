/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Nonintegral-population product-form analysis and the AMVA joint-dependence
 * hook. Oracles, in order of strength:
 *   1. Interpolation identity. Both pfqn_nintmva and pfqn_dnc are continuations
 *      through the integral points, so at every INTEGER population they must
 *      reproduce exact MVA and the exact convolution normalizing constant to
 *      solver precision. That is an equality, not a tolerance band, and it is
 *      what distinguishes a continuation from an interpolation.
 *   2. Closed forms on one- and two-station models computed by hand.
 *   3. Monotonicity and smoothness between the integral points.
 *   4. pfqn_jdfun is the reciprocal of the supplied callable by definition, and
 *      agrees term for term with pfqn_cdfun on any callable both accept.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/pfqn/pfqn_ca.h"
#include "line/api/pfqn/pfqn_cdfun.h"
#include "line/api/pfqn/pfqn_cub_evals.h"
#include "line/api/pfqn/pfqn_dnc.h"
#include "line/api/pfqn/pfqn_jdfun.h"
#include "line/api/pfqn/pfqn_mva.h"
#include "line/api/pfqn/pfqn_nintmva.h"

using line::Matrix;
using line::Rational;
using line::Real50;
using line::pfqn::CdScaling;
using line::pfqn::JdScaling;
using line::pfqn::pfqn_cdfun;
using line::pfqn::pfqn_cub_evals;
using line::pfqn::pfqn_dnc;
using line::pfqn::pfqn_jdfun;
using line::pfqn::pfqn_nintmva;

namespace {

constexpr double TOL = 1e-9;

}  // namespace

// ---------------------------------------------------------------------------
// pfqn_jdfun
// ---------------------------------------------------------------------------

TEST_CASE("pfqn_jdfun returns the reciprocal of a broadcast scaling") {
    Matrix<double> nvec(3, 2);
    nvec(0, 0) = 2; nvec(0, 1) = 1;
    nvec(1, 0) = 0; nvec(1, 1) = 4;
    nvec(2, 0) = 3; nvec(2, 1) = 3;
    // eta_i(n) = min(sum_r n_ir, 2), the flagship multiserver scaling
    std::vector<JdScaling<double>> jd(3);
    for (std::size_t i = 0; i < 3; ++i)
        jd[i] = [](const std::vector<double>& n) {
            double s = 0;
            for (double v : n) s += v;
            return std::vector<double>(1, std::min(s, 2.0));
        };
    const std::vector<double> r = pfqn_jdfun(nvec, jd);
    CHECK(r[0] == doctest::Approx(1.0 / 2.0).epsilon(TOL));
    CHECK(r[1] == doctest::Approx(1.0 / 2.0).epsilon(TOL));
    CHECK(r[2] == doctest::Approx(1.0 / 2.0).epsilon(TOL));
}

TEST_CASE("pfqn_jdfun selects the requested class from a per-class scaling") {
    Matrix<double> nvec(2, 3);
    nvec(0, 0) = 1; nvec(0, 1) = 2; nvec(0, 2) = 3;
    nvec(1, 0) = 4; nvec(1, 1) = 0; nvec(1, 2) = 1;
    // eta_{i,r}(n) = 1 + n_ir, a genuinely joint per-class rate
    std::vector<JdScaling<double>> jd(2);
    for (std::size_t i = 0; i < 2; ++i)
        jd[i] = [](const std::vector<double>& n) {
            std::vector<double> out(n.size());
            for (std::size_t r = 0; r < n.size(); ++r) out[r] = 1.0 + n[r];
            return out;
        };
    for (std::size_t cls = 0; cls < 3; ++cls) {
        const std::vector<double> r = pfqn_jdfun(nvec, jd, cls);
        CHECK(r[0] == doctest::Approx(1.0 / (1.0 + nvec(0, cls))).epsilon(TOL));
        CHECK(r[1] == doctest::Approx(1.0 / (1.0 + nvec(1, cls))).epsilon(TOL));
    }
}

TEST_CASE("pfqn_jdfun agrees with pfqn_cdfun term for term") {
    Matrix<double> nvec(3, 2);
    nvec(0, 0) = 5; nvec(0, 1) = 2;
    nvec(1, 0) = 1; nvec(1, 1) = 1;
    nvec(2, 0) = 0; nvec(2, 1) = 7;
    auto f = [](const std::vector<double>& n) {
        return std::vector<double>(1, 1.0 + 0.5 * n[0] + 0.25 * n[1]);
    };
    std::vector<JdScaling<double>> jd(3, f);
    std::vector<CdScaling<double>> cd(3, f);
    const std::vector<double> rj = pfqn_jdfun(nvec, jd);
    const std::vector<double> rc = pfqn_cdfun(nvec, cd);
    REQUIRE(rj.size() == rc.size());
    for (std::size_t i = 0; i < rj.size(); ++i) CHECK(rj[i] == doctest::Approx(rc[i]).epsilon(TOL));
}

TEST_CASE("pfqn_jdfun leaves an absent callable at one and is exact-capable") {
    Matrix<Rational> nvec(2, 1);
    nvec(0, 0) = Rational(3);
    nvec(1, 0) = Rational(4);
    std::vector<JdScaling<Rational>> jd(2);
    jd[1] = [](const std::vector<Rational>& n) {
        return std::vector<Rational>(1, Rational(1) + n[0]);
    };
    const std::vector<Rational> r = pfqn_jdfun(nvec, jd);
    CHECK(r[0] == Rational(1));
    CHECK(r[1] == Rational(1, 5));
}

TEST_CASE("pfqn_jdfun rejects a zero scaling and a mis-sized vector") {
    Matrix<double> nvec(2, 1);
    nvec(0, 0) = 1;
    nvec(1, 0) = 1;
    std::vector<JdScaling<double>> jd(2);
    jd[0] = [](const std::vector<double>&) { return std::vector<double>(1, 0.0); };
    CHECK_THROWS(pfqn_jdfun(nvec, jd));
    std::vector<JdScaling<double>> wrong(3);
    wrong[0] = [](const std::vector<double>&) { return std::vector<double>(1, 1.0); };
    CHECK_THROWS(pfqn_jdfun(nvec, wrong));
}

// ---------------------------------------------------------------------------
// pfqn_cub_evals
// ---------------------------------------------------------------------------

TEST_CASE("pfqn_cub_evals counts the Grundmann-Moeller nodes") {
    // M = 3 -> n = 2; order 0: C(2,2) = 1; order 1: + C(4,2) = 6 -> 7;
    // order 2: + C(6,2) = 15 -> 22; order 3: + C(8,2) = 28 -> 50
    CHECK(pfqn_cub_evals(3, 0) == doctest::Approx(1.0));
    CHECK(pfqn_cub_evals(3, 1) == doctest::Approx(7.0));
    CHECK(pfqn_cub_evals(3, 2) == doctest::Approx(22.0));
    CHECK(pfqn_cub_evals(3, 3) == doctest::Approx(50.0));
    // M = 1 -> n = 0, every binomial C(2d,0) = 1
    CHECK(pfqn_cub_evals(1, 4) == doctest::Approx(5.0));
}

TEST_CASE("pfqn_cub_evals prices the think-time v-quadrature") {
    CHECK(pfqn_cub_evals(3, 2, 0.0) == doctest::Approx(22.0));
    CHECK(pfqn_cub_evals(3, 2, 1.5) == doctest::Approx(22.0 * line::pfqn::CUB_V_STEPS));
    // below FineTol the think time does not count
    CHECK(pfqn_cub_evals(3, 2, 1e-12) == doctest::Approx(22.0));
}

TEST_CASE("pfqn_cub_evals is the budget pfqn_nc lowers the order against") {
    // With a think time the budget of 1e7 admits at most 1000 simplex nodes.
    // M = 4 (n = 3): order 5 costs C(3,3)+C(5,3)+C(7,3)+C(9,3)+C(11,3)+C(13,3)
    //              = 1 + 10 + 35 + 84 + 165 + 286 = 581 -> 5.81e6, fits;
    // order 6 adds C(15,3) = 455 -> 1036 -> 1.036e7, does not.
    CHECK(pfqn_cub_evals(4, 5, 1.0) < line::pfqn::CUB_MAX_EVALS);
    CHECK(pfqn_cub_evals(4, 6, 1.0) > line::pfqn::CUB_MAX_EVALS);
}

// ---------------------------------------------------------------------------
// pfqn_nintmva
// ---------------------------------------------------------------------------

TEST_CASE("pfqn_nintmva reproduces exact MVA at every integer population") {
    const std::vector<double> L = {0.8, 0.4, 0.2};
    Matrix<double> Lm(3, 1);
    for (std::size_t i = 0; i < 3; ++i) Lm(i, 0) = L[i];
    for (int n = 1; n <= 8; ++n) {
        const std::vector<int> Nv(1, n);
        Matrix<double> Zm(1, 1, 0.0);
        const auto exact = line::pfqn::pfqn_mva(Lm, Nv, Zm);
        const auto nint = pfqn_nintmva<double>(L, static_cast<double>(n), 0.0);
        CHECK(nint.X == doctest::Approx(exact.XN[0]).epsilon(TOL));
        for (std::size_t i = 0; i < 3; ++i)
            CHECK(nint.Q[i] == doctest::Approx(exact.QN(i, 0)).epsilon(TOL));
    }
}

TEST_CASE("pfqn_nintmva reproduces exact MVA with a think time") {
    const std::vector<double> L = {1.0, 0.5};
    Matrix<double> Lm(2, 1);
    Lm(0, 0) = 1.0;
    Lm(1, 0) = 0.5;
    for (int n = 1; n <= 6; ++n) {
        const std::vector<int> Nv(1, n);
        Matrix<double> Zm(1, 1, 3.0);
        const auto exact = line::pfqn::pfqn_mva(Lm, Nv, Zm);
        const auto nint = pfqn_nintmva<double>(L, static_cast<double>(n), 3.0);
        CHECK(nint.X == doctest::Approx(exact.XN[0]).epsilon(TOL));
        CHECK(nint.U[0] == doctest::Approx(exact.UN(0, 0)).epsilon(TOL));
        CHECK(nint.R[1] == doctest::Approx(exact.CN(1, 0)).epsilon(TOL));
    }
}

TEST_CASE("pfqn_nintmva one station closed form at a fractional population") {
    // Single queueing station, no delay: the base n0 is the only step, so
    // R = L, X = n0/L, Q = n0 for the fractional part, and thereafter the usual
    // unit steps. At N = 0.4: R = 0.5, X = 0.8, Q = 0.4.
    const std::vector<double> L = {0.5};
    const auto r = pfqn_nintmva<double>(L, 0.4, 0.0);
    CHECK(r.R[0] == doctest::Approx(0.5).epsilon(TOL));
    CHECK(r.X == doctest::Approx(0.8).epsilon(TOL));
    CHECK(r.Q[0] == doctest::Approx(0.4).epsilon(TOL));
    CHECK(r.U[0] == doctest::Approx(0.4).epsilon(TOL));
}

TEST_CASE("pfqn_nintmva interpolates monotonically between the integers") {
    const std::vector<double> L = {0.7, 0.3};
    double prev = 0.0;
    for (int k = 0; k <= 20; ++k) {
        const double N = 1.0 + 0.1 * k;
        const auto r = pfqn_nintmva<double>(L, N, 1.0);
        CHECK(r.X > prev - 1e-12);
        prev = r.X;
        double qsum = 0;
        for (double q : r.Q) qsum += q;
        // the queueing stations hold N minus the jobs thinking
        CHECK(qsum + r.X * 1.0 == doctest::Approx(N).epsilon(1e-9));
    }
}

TEST_CASE("pfqn_nintmva at zero population and in exact arithmetic") {
    const std::vector<double> L = {1.0, 2.0};
    const auto z = pfqn_nintmva<double>(L, 0.0, 0.0);
    CHECK(z.X == doctest::Approx(0.0));
    CHECK(z.Q[0] == doctest::Approx(0.0));
    const std::vector<Rational> Lr = {Rational(1, 2), Rational(1, 4)};
    const auto r = pfqn_nintmva<Rational>(Lr, Rational(2), Rational(0));
    // exact MVA at N = 2: R = (1/2, 1/4), X = 2/(3/4) ... stepped
    const auto r1 = pfqn_nintmva<Rational>(Lr, Rational(1), Rational(0));
    CHECK(r1.X == Rational(4, 3));
    CHECK(r.X > r1.X);
    CHECK(Rational(r.Q[0] + r.Q[1]) == Rational(2));
}

TEST_CASE("pfqn_nintmva rejects a negative population") {
    const std::vector<double> L = {1.0};
    CHECK_THROWS(pfqn_nintmva<double>(L, -1.0, 0.0));
}

// ---------------------------------------------------------------------------
// pfqn_dnc
// ---------------------------------------------------------------------------

TEST_CASE("pfqn_dnc reproduces the convolution constant at every integer") {
    const std::vector<double> L = {0.9, 0.5, 0.2};
    Matrix<double> Lm(3, 1);
    for (std::size_t i = 0; i < 3; ++i) Lm(i, 0) = L[i];
    for (int n = 1; n <= 10; ++n) {
        const std::vector<int> Nv(1, n);
        Matrix<double> Zm(1, 1, 0.0);
        const auto ca = line::pfqn::pfqn_ca(Lm, Nv, Zm);
        const auto d = pfqn_dnc<double>(L, static_cast<double>(n));
        CHECK(d.lG == doctest::Approx(ca.lG).epsilon(1e-8));
    }
}

TEST_CASE("pfqn_dnc throughput matches exact MVA at every integer") {
    const std::vector<double> L = {0.9, 0.5, 0.2};
    Matrix<double> Lm(3, 1);
    for (std::size_t i = 0; i < 3; ++i) Lm(i, 0) = L[i];
    for (int n = 1; n <= 10; ++n) {
        const std::vector<int> Nv(1, n);
        Matrix<double> Zm(1, 1, 0.0);
        const auto exact = line::pfqn::pfqn_mva(Lm, Nv, Zm);
        const auto d = pfqn_dnc<double>(L, static_cast<double>(n));
        CHECK(d.X == doctest::Approx(exact.XN[0]).epsilon(1e-8));
    }
}

TEST_CASE("pfqn_dnc handles repeated loads through the multiplicity branch") {
    // Three identical stations: G(n) = C(n+2,2) L^n exactly.
    const double L1 = 0.6;
    const std::vector<double> L = {L1, L1, L1};
    for (int n = 1; n <= 8; ++n) {
        const double want = (n + 1.0) * (n + 2.0) / 2.0 * std::pow(L1, n);
        const auto d = pfqn_dnc<double>(L, static_cast<double>(n));
        CHECK(d.G == doctest::Approx(want).epsilon(1e-7));
    }
}

TEST_CASE("pfqn_dnc handles a partially repeated load set") {
    const double a = 0.8, b = 0.3;
    const std::vector<double> L = {a, a, b};
    Matrix<double> Lm(3, 1);
    Lm(0, 0) = a; Lm(1, 0) = a; Lm(2, 0) = b;
    for (int n = 1; n <= 8; ++n) {
        const std::vector<int> Nv(1, n);
        Matrix<double> Zm(1, 1, 0.0);
        const auto ca = line::pfqn::pfqn_ca(Lm, Nv, Zm);
        const auto d = pfqn_dnc<double>(L, static_cast<double>(n));
        CHECK(d.lG == doctest::Approx(ca.lG).epsilon(1e-7));
    }
}

TEST_CASE("pfqn_dnc interpolates smoothly and monotonically in N") {
    const std::vector<double> L = {0.9, 0.4};
    // G is NOT monotone: the multiplicity factor C(N+g-1,g-1) pushes it up while
    // max(L)^N < 1 pulls it down, so it peaks (near N = 1.75 here) and then
    // decays. Only positivity and the rising throughput are invariants.
    double prevX = 0.0;
    for (int k = 1; k <= 40; ++k) {
        const double N = 0.25 * k;
        const auto d = pfqn_dnc<double>(L, N);
        CHECK(d.G > 0.0);
        CHECK(d.X > prevX - 1e-9);
        prevX = d.X;
    }
    // and stays below the bottleneck bound 1/max(L)
    const auto far = pfqn_dnc<double>(L, 60.0);
    CHECK(far.X < 1.0 / 0.9 + 1e-9);
    CHECK(far.X > 1.0 / 0.9 - 1e-2);
}

TEST_CASE("pfqn_dnc single station is the geometric constant") {
    const std::vector<double> L = {0.7};
    const auto d = pfqn_dnc<double>(L, 3.5);
    CHECK(d.G == doctest::Approx(std::pow(0.7, 3.5)).epsilon(1e-10));
    CHECK(d.X == doctest::Approx(1.0 / 0.7).epsilon(1e-10));
}

TEST_CASE("pfqn_dnc rejects an empty demand set and a negative population") {
    const std::vector<double> zeros = {0.0, 0.0};
    CHECK_THROWS(pfqn_dnc<double>(zeros, 1.0));
    const std::vector<double> L = {1.0};
    CHECK_THROWS(pfqn_dnc<double>(L, -0.5));
}

TEST_CASE("pfqn_dnc continues below one job but refuses the empty network") {
    const std::vector<double> L = {0.9, 0.4};
    // the series is analytic for n > -1, so G(N-1) exists at N = 0.5 and the
    // throughput is a genuine value, not a blank
    const auto half = pfqn_dnc<double>(L, 0.5);
    CHECK(!std::isnan(half.X));
    CHECK(half.X > 0.0);
    CHECK(half.G > 0.0);
    // at N = 0 there is no job to complete and the reference returns NaN
    const auto zero = pfqn_dnc<double>(L, 0.0);
    CHECK(std::isnan(zero.X));
    CHECK(zero.G == doctest::Approx(1.0).epsilon(1e-12));
}
