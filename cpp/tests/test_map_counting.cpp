/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * MAP distribution and counting-process descriptors: map_cdf, map_pdf,
 * map_acfc, map_count_mean, map_count_var, map_varcount, map_count_moment and
 * mmap_count_var.
 *
 * Oracles, in order of strength:
 *   - closed forms: the Erlang-2 CDF and PDF, the Poisson counting moments
 *     (Touchard), Var[N(t)] = lambda t for Poisson, rho(k) = 0 for Poisson;
 *   - defining properties: a CDF is nondecreasing, starts at 0 and tends to 1,
 *     and its density integrates to it;
 *   - internal consistency: map_count_var and map_varcount are two spellings
 *     of the same quantity, and map_count_moment must reproduce map_count_mean
 *     at order 1 and map_count_var at order 2;
 *   - MATLAB reference values, produced by R2025a with lineStart on the path.
 *
 * MATLAB's map_count_moment differentiates the counting MGF NUMERICALLY
 * (derivest) up to order 4, so its high-order values carry a large error; the
 * exact Poisson moments are used to show which of the two is right, and the
 * MATLAB comparison is asserted only at the accuracy MATLAB actually delivers.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/mam/map_acfc.h"
#include "line/api/mam/map_cdf.h"
#include "line/api/mam/map_count_mean.h"
#include "line/api/mam/map_count_moment.h"
#include "line/api/mam/map_count_var.h"
#include "line/api/mam/map_pdf.h"
#include "line/api/mam/map_varcount.h"
#include "line/api/mam/mmap_count_var.h"

using line::Matrix;
using line::mam::Map;
using line::mam::Mmap;
using namespace line::mam;

namespace {

constexpr double MATLAB_TOL = 1e-10;

/** Two-state MMPP used throughout: rates 1 and 3, switching 0.05 / 0.02. */
Map<double> mmpp() {
    Map<double> m;
    m.D0 = Matrix<double>{{-1.05, 0.05}, {0.02, -3.02}};
    m.D1 = Matrix<double>{{1.0, 0.0}, {0.0, 3.0}};
    return m;
}

/** Erlang-2 with rate 2 per phase: mean 1, SCV 1/2. */
Map<double> erlang2() {
    Map<double> m;
    m.D0 = Matrix<double>{{-2.0, 2.0}, {0.0, -2.0}};
    m.D1 = Matrix<double>{{0.0, 0.0}, {2.0, 0.0}};
    return m;
}

Map<double> poisson(double lam) {
    Map<double> m;
    m.D0 = Matrix<double>(1, 1, -lam);
    m.D1 = Matrix<double>(1, 1, lam);
    return m;
}

}  // namespace

TEST_CASE("map_cdf is a proper distribution function") {
    const Map<double> m = mmpp();
    std::vector<double> pts;
    for (int k = 0; k <= 200; ++k) pts.push_back(0.05 * k);
    const std::vector<double> F = map_cdf(m, pts);
    CHECK(F[0] == 0.0);  // F(0) = 0 exactly
    for (std::size_t k = 1; k < F.size(); ++k) {
        CHECK(F[k] >= F[k - 1] - 1e-15);  // nondecreasing
        CHECK(F[k] >= 0.0);
        CHECK(F[k] <= 1.0 + 1e-15);
    }
    const std::vector<double> far = map_cdf(m, std::vector<double>{50.0});
    CHECK(far[0] == doctest::Approx(1.0).epsilon(1e-12));  // tends to one
}

TEST_CASE("map_cdf and map_pdf match the Erlang-2 closed form") {
    const Map<double> m = erlang2();
    const std::vector<double> t{0.1, 0.5, 1.0, 2.0, 5.0};
    const std::vector<double> F = map_cdf(m, t);
    const std::vector<double> f = map_pdf(m, t);
    for (std::size_t k = 0; k < t.size(); ++k) {
        const double x = t[k];
        CHECK(F[k] == doctest::Approx(1.0 - (1.0 + 2.0 * x) * std::exp(-2.0 * x)).epsilon(1e-13));
        CHECK(f[k] == doctest::Approx(4.0 * x * std::exp(-2.0 * x)).epsilon(1e-13));
    }
    // MATLAB reference (map_cdf.m, map_pdf.m)
    CHECK(F[0] == doctest::Approx(0.017523096306421793).epsilon(MATLAB_TOL));
    CHECK(F[2] == doctest::Approx(0.59399415029016189).epsilon(MATLAB_TOL));
    CHECK(F[4] == doctest::Approx(0.99950060077261271).epsilon(MATLAB_TOL));
    CHECK(f[0] == doctest::Approx(0.32749230123119277).epsilon(MATLAB_TOL));
    CHECK(f[2] == doctest::Approx(0.54134113294645081).epsilon(MATLAB_TOL));
    CHECK(f[4] == doctest::Approx(0.000907998595249697).epsilon(MATLAB_TOL));
}

TEST_CASE("map_cdf and map_pdf of an MMPP match MATLAB") {
    const Map<double> m = mmpp();
    const std::vector<double> t{0.1, 0.5, 1.0, 2.0, 5.0};
    const std::vector<double> F = map_cdf(m, t);
    const std::vector<double> f = map_pdf(m, t);
    const double Fref[5] = {0.2397955648876845, 0.7310163373424603, 0.91211527447636631,
                            0.98202088857917091, 0.99931699781826644};
    const double fref[5] = {2.0658870429523408, 0.66110795535236977, 0.17589639588143369,
                            0.022949680653318462, 0.00071727923247247901};
    for (std::size_t k = 0; k < 5; ++k) {
        CHECK(F[k] == doctest::Approx(Fref[k]).epsilon(MATLAB_TOL));
        CHECK(f[k] == doctest::Approx(fref[k]).epsilon(MATLAB_TOL));
    }
    // The density integrates to the distribution: trapezoid on a fine grid.
    std::vector<double> grid;
    const double dt = 0.002;
    for (int k = 0; k <= 500; ++k) grid.push_back(dt * k);
    const std::vector<double> dens = map_pdf(m, grid);
    double area = 0.0;
    for (std::size_t k = 1; k < grid.size(); ++k) area += 0.5 * dt * (dens[k] + dens[k - 1]);
    const std::vector<double> Fend = map_cdf(m, std::vector<double>{grid.back()});
    CHECK(area == doctest::Approx(Fend[0]).epsilon(1e-5));  // trapezoid discretization
}

TEST_CASE("map_cdf and map_pdf of a Poisson process are the exponential law") {
    const double lam = 0.7;
    const Map<double> m = poisson(lam);
    const std::vector<double> t{0.1, 1.0, 5.0};
    const std::vector<double> F = map_cdf(m, t);
    const std::vector<double> f = map_pdf(m, t);
    for (std::size_t k = 0; k < t.size(); ++k) {
        CHECK(F[k] == doctest::Approx(1.0 - std::exp(-lam * t[k])).epsilon(1e-14));
        CHECK(f[k] == doctest::Approx(lam * std::exp(-lam * t[k])).epsilon(1e-14));
    }
    // MATLAB reference
    CHECK(F[0] == doctest::Approx(0.067606180094051727).epsilon(MATLAB_TOL));
    CHECK(F[1] == doctest::Approx(0.50341469620859047).epsilon(MATLAB_TOL));
    CHECK(F[2] == doctest::Approx(0.96980261657768152).epsilon(MATLAB_TOL));
}

TEST_CASE("map_count_mean is lambda t, exactly") {
    const Map<double> m = mmpp();
    const std::vector<double> t{1.0, 5.0, 20.0};
    const std::vector<double> mu = map_count_mean(m, t);
    const double lam = map_lambda(m);
    CHECK(lam == doctest::Approx(2.4285714285714284).epsilon(1e-14));
    for (std::size_t k = 0; k < t.size(); ++k) CHECK(mu[k] == lam * t[k]);
    // MATLAB reference
    CHECK(mu[0] == doctest::Approx(2.4285714285714284).epsilon(MATLAB_TOL));
    CHECK(mu[2] == doctest::Approx(48.571428571428569).epsilon(MATLAB_TOL));
}

TEST_CASE("map_count_var of a Poisson process is lambda t") {
    const double lam = 0.7;
    const Map<double> m = poisson(lam);
    const std::vector<double> t{1.0, 5.0, 20.0};
    const std::vector<double> v = map_count_var(m, t);
    for (std::size_t k = 0; k < t.size(); ++k)
        CHECK(v[k] == doctest::Approx(lam * t[k]).epsilon(1e-13));
}

TEST_CASE("map_count_var and map_varcount agree and match MATLAB") {
    const std::vector<double> t{1.0, 5.0, 20.0};
    SUBCASE("MMPP") {
        const Map<double> m = mmpp();
        const std::vector<double> a = map_count_var(m, t);
        const std::vector<double> b = map_varcount(m, t);
        const double ref[3] = {3.2261790607074552, 30.364627977913699, 264.01398215463757};
        for (std::size_t k = 0; k < 3; ++k) {
            CHECK(a[k] == doctest::Approx(b[k]).epsilon(1e-12));
            CHECK(a[k] == doctest::Approx(ref[k]).epsilon(MATLAB_TOL));
        }
        // Overdispersed: the MMPP is burstier than Poisson at every scale.
        const double lam = map_lambda(m);
        for (std::size_t k = 0; k < 3; ++k) CHECK(a[k] > lam * t[k]);
    }
    SUBCASE("Erlang-2 renewal process") {
        const Map<double> m = erlang2();
        const std::vector<double> a = map_count_var(m, t);
        const std::vector<double> b = map_varcount(m, t);
        const double ref[3] = {0.62271054513890811, 2.6249999997423563, 10.125000000000004};
        for (std::size_t k = 0; k < 3; ++k) {
            CHECK(a[k] == doctest::Approx(b[k]).epsilon(1e-12));
            CHECK(a[k] == doctest::Approx(ref[k]).epsilon(MATLAB_TOL));
        }
        // Underdispersed: an Erlang renewal stream is smoother than Poisson.
        // Var[N(t)] = lambda SCV t + O(1), so the SLOPE, not the ratio, is the
        // asymptotic index of dispersion lambda SCV = 1/2.
        CHECK((a[2] - a[1]) / 15.0 == doctest::Approx(0.5).epsilon(1e-9));
    }
}

TEST_CASE("map_acfc is zero for a Poisson process and matches MATLAB otherwise") {
    const std::vector<unsigned> lags{1u, 2u, 3u};
    SUBCASE("Poisson counts are uncorrelated") {
        const Map<double> m = poisson(0.7);
        const std::vector<double> r = map_acfc(m, lags, 1.0);
        for (std::size_t k = 0; k < r.size(); ++k) CHECK(std::fabs(r[k]) < 1e-14);
    }
    SUBCASE("MMPP counts are positively correlated and decay geometrically") {
        const Map<double> m = mmpp();
        const std::vector<unsigned> ks{1u, 2u, 3u, 5u};
        const std::vector<double> r = map_acfc(m, ks, 1.0);
        const double ref[4] = {0.23602182226959623, 0.22006528844711162, 0.20518751492390674,
                               0.17838145590011367};
        for (std::size_t k = 0; k < 4; ++k) {
            CHECK(r[k] == doctest::Approx(ref[k]).epsilon(MATLAB_TOL));
            CHECK(r[k] > 0.0);
        }
        // The decay factor is exp(-(q01+q10) u), the second eigenvalue of Q.
        const double decay = std::exp(-(0.05 + 0.02) * 1.0);
        CHECK(r[1] / r[0] == doctest::Approx(decay).epsilon(1e-10));
        CHECK(r[2] / r[1] == doctest::Approx(decay).epsilon(1e-10));

        const std::vector<double> r5 = map_acfc(m, lags, 5.0);
        const double ref5[3] = {0.47847779998275697, 0.33717760684266168, 0.23760504366188259};
        for (std::size_t k = 0; k < 3; ++k)
            CHECK(r5[k] == doctest::Approx(ref5[k]).epsilon(MATLAB_TOL));
    }
    SUBCASE("Erlang-2 counts are negatively correlated") {
        const Map<double> m = erlang2();
        const std::vector<double> r = map_acfc(m, lags, 1.0);
        const double ref[3] = {-0.096724733543280886, -0.0017715752911877704,
                               -3.2447533297592775e-05};
        for (std::size_t k = 0; k < 3; ++k) CHECK(r[k] == doctest::Approx(ref[k]).epsilon(1e-9));
    }
}

TEST_CASE("map_count_moment reproduces the exact Poisson moments") {
    // N(t) ~ Poisson(lambda t): E[N] = a, E[N^2] = a^2 + a,
    // E[N^3] = a^3 + 3a^2 + a, E[N^4] = a^4 + 6a^3 + 7a^2 + a, a = lambda t.
    const double lam = 0.7, t = 3.0, a = lam * t;
    const Map<double> m = poisson(lam);
    const std::vector<unsigned> ord{0u, 1u, 2u, 3u, 4u};
    const std::vector<double> M = map_count_moment(m, t, ord);
    CHECK(M[0] == doctest::Approx(1.0).epsilon(1e-14));
    CHECK(M[1] == doctest::Approx(a).epsilon(1e-12));
    CHECK(M[2] == doctest::Approx(a * a + a).epsilon(1e-12));
    CHECK(M[3] == doctest::Approx(a * a * a + 3.0 * a * a + a).epsilon(1e-12));
    CHECK(M[4] == doctest::Approx(a * a * a * a + 6.0 * a * a * a + 7.0 * a * a + a).epsilon(1e-12));

    // MATLAB (derivest) returns 107.9823181930049 for the fourth moment, which
    // is 1.8e-3 away from the exact 107.9841: numerical differentiation, not
    // the model, is the error. Orders 1 and 2 agree to 1e-10.
    CHECK(M[1] == doctest::Approx(2.0999999999999401).epsilon(1e-12));
    CHECK(M[2] == doctest::Approx(6.5100000000205114).epsilon(1e-10));
    CHECK(std::fabs(M[4] - 107.9823181930049) > 1e-4);
    CHECK(M[4] == doctest::Approx(107.9841).epsilon(1e-12));
}

TEST_CASE("map_count_moment is consistent with map_count_mean and map_count_var") {
    const std::vector<unsigned> ord{1u, 2u};
    SUBCASE("MMPP at t = 1") {
        const Map<double> m = mmpp();
        const std::vector<double> M = map_count_moment(m, 1.0, ord);
        const std::vector<double> mu = map_count_mean(m, std::vector<double>{1.0});
        const std::vector<double> v = map_count_var(m, std::vector<double>{1.0});
        CHECK(M[0] == doctest::Approx(mu[0]).epsilon(1e-12));
        CHECK(M[1] - M[0] * M[0] == doctest::Approx(v[0]).epsilon(1e-10));
        // MATLAB reference (derivest): order 1 and 2 at t = 1
        CHECK(M[0] == doctest::Approx(2.4285714285714564).epsilon(1e-12));
        CHECK(M[1] == doctest::Approx(9.124138244380072).epsilon(1e-9));
    }
    SUBCASE("MMPP at t = 5") {
        const Map<double> m = mmpp();
        const std::vector<double> M = map_count_moment(m, 5.0, ord);
        const std::vector<double> v = map_count_var(m, std::vector<double>{5.0});
        CHECK(M[0] == doctest::Approx(12.142857142857142).epsilon(1e-12));
        CHECK(M[1] - M[0] * M[0] == doctest::Approx(v[0]).epsilon(1e-9));
        CHECK(M[1] == doctest::Approx(177.81360756935092).epsilon(1e-8));
    }
    SUBCASE("Erlang-2 at t = 2") {
        const Map<double> m = erlang2();
        const std::vector<double> M = map_count_moment(m, 2.0, ord);
        const std::vector<double> v = map_count_var(m, std::vector<double>{2.0});
        CHECK(M[0] == doctest::Approx(2.0).epsilon(1e-12));  // lambda = 1
        CHECK(M[1] - M[0] * M[0] == doctest::Approx(v[0]).epsilon(1e-9));
        CHECK(M[1] == doctest::Approx(5.1249580671664834).epsilon(1e-8));
    }
}

TEST_CASE("mmap_count_var splits the counts of a marked MAP and matches MATLAB") {
    const Map<double> base = mmpp();
    Mmap<double> mm;
    mm.D0 = base.D0;
    mm.D1 = base.D1;
    Matrix<double> A = base.D1, B = base.D1;
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j) {
            A(i, j) *= 0.3;
            B(i, j) *= 0.7;
        }
    mm.Dc.push_back(A);
    mm.Dc.push_back(B);

    const std::vector<double> v1 = mmap_count_var(mm, 1.0);
    const std::vector<double> v5 = mmap_count_var(mm, 5.0);
    CHECK(v1[0] == doctest::Approx(0.80035611546367047).epsilon(MATLAB_TOL));
    CHECK(v1[1] == doctest::Approx(2.0908277397466506).epsilon(MATLAB_TOL));
    CHECK(v5[0] == doctest::Approx(5.2828165180122326).epsilon(MATLAB_TOL));
    CHECK(v5[1] == doctest::Approx(17.428667709177695).epsilon(MATLAB_TOL));

    // A class is itself a MAP (D0 + the other class matrices, Dc), and its
    // per-class count variance must equal that MAP's map_count_var.
    for (std::size_t c = 0; c < 2; ++c) {
        Map<double> mc;
        mc.D1 = mm.Dc[c];
        mc.D0 = mm.D0;
        for (std::size_t i = 0; i < 2; ++i)
            for (std::size_t j = 0; j < 2; ++j)
                mc.D0(i, j) += mm.Dc[1 - c](i, j);  // the other class is hidden
        const std::vector<double> ref = map_count_var(mc, std::vector<double>{1.0});
        CHECK(v1[c] == doctest::Approx(ref[0]).epsilon(1e-10));
    }

    // Marking a Poisson process splits it into independent Poisson streams.
    Mmap<double> pm;
    pm.D0 = Matrix<double>(1, 1, -1.0);
    pm.D1 = Matrix<double>(1, 1, 1.0);
    pm.Dc.push_back(Matrix<double>(1, 1, 0.25));
    pm.Dc.push_back(Matrix<double>(1, 1, 0.75));
    const std::vector<double> pv = mmap_count_var(pm, 4.0);
    CHECK(pv[0] == doctest::Approx(1.0).epsilon(1e-12));  // 0.25 * 4
    CHECK(pv[1] == doctest::Approx(3.0).epsilon(1e-12));  // 0.75 * 4
}

TEST_CASE("the counting descriptors instantiate at high precision") {
    // Same MMPP at 50 digits: the double answers are a prefix of these, and the
    // exponential is evaluated with a precision-aware scaling (see expm.h), so
    // the Poisson identity Var[N(t)] = lambda t holds far below double epsilon.
    using line::Real50;
    Map<Real50> m;
    m.D0 = Matrix<Real50>{{Real50("-1.05"), Real50("0.05")}, {Real50("0.02"), Real50("-3.02")}};
    m.D1 = Matrix<Real50>{{Real50(1), Real50(0)}, {Real50(0), Real50(3)}};
    const std::vector<Real50> t{Real50(1), Real50(5)};
    const std::vector<Real50> v = map_count_var(m, t);
    const std::vector<Real50> w = map_varcount(m, t);
    CHECK(static_cast<double>(v[0]) == doctest::Approx(3.2261790607074552).epsilon(MATLAB_TOL));
    CHECK(static_cast<double>(v[1]) == doctest::Approx(30.364627977913699).epsilon(MATLAB_TOL));
    // The two spellings are algebraically identical; at 50 digits they agree to
    // 1e-40, where at double they differ in the 14th digit.
    CHECK(static_cast<double>(abs(v[0] - w[0])) < 1e-40);

    Map<Real50> p;
    p.D0 = Matrix<Real50>(1, 1, Real50("-0.7"));
    p.D1 = Matrix<Real50>(1, 1, Real50("0.7"));
    const std::vector<Real50> pv = map_count_var(p, t);
    CHECK(static_cast<double>(abs(pv[1] - Real50("3.5"))) < 1e-40);
    const std::vector<Real50> F = map_cdf(p, std::vector<Real50>{Real50(1)});
    const Real50 exact = Real50(1) - exp(Real50("-0.7"));
    CHECK(static_cast<double>(abs(F[0] - exact)) < 1e-45);
}
