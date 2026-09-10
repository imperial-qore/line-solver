/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Harel-Namn-Sturm bounds, the replica multiplicity fold, and the asymptotic
 * cache miss ratio. Oracles, in order of strength:
 *   1. Cross-algorithm identity. G(n) formed by the Newton-Girard recurrence
 *      inside pfqn_harel_bounds IS the normalizing constant of the single-class
 *      closed load-independent network, so it must equal pfqn_ca term for term;
 *      and TH(n) = G(n-1)/G(n) must equal the exact pfqn_mva throughput at
 *      population n. This is what certifies the port against the reference's
 *      hardcoded G(0)..G(7) polynomials without reproducing the table.
 *   2. Bound semantics. LB <= X_exact <= UB(n) at every admissible n, checked
 *      against pfqn_mva on models whose exact throughput is available.
 *   3. Exact arithmetic. G and UB are field expressions, so the Rational
 *      instantiation must agree with the double one to full double precision.
 *   4. Definitional invariants for cache_miss_asy: the ratio lies in [0,1], it
 *      is 1 at zero capacity, and it decreases as capacity grows.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/cache/cache_miss_asy.h"
#include "line/api/pfqn/pfqn_ca.h"
#include "line/api/pfqn/pfqn_harel_bounds.h"
#include "line/api/pfqn/pfqn_mva.h"
#include "line/api/pfqn/pfqn_unique.h"

using line::Matrix;
using line::Rational;
using line::cache::cache_miss_asy;
using line::pfqn::pfqn_ca;
using line::pfqn::pfqn_combine_mi;
using line::pfqn::pfqn_harel_bounds;
using line::pfqn::pfqn_harel_lb;
using line::pfqn::pfqn_harel_ub;

namespace {

constexpr double TOL = 1e-10;

/** The single-class closed network with demands rho, as pfqn_ca takes it. */
Matrix<double> as_column(const std::vector<double>& rho) {
    Matrix<double> L(rho.size(), 1);
    for (std::size_t i = 0; i < rho.size(); ++i) L(i, 0) = rho[i];
    return L;
}

}  // namespace

// ---------------------------------------------------------------------------
// pfqn_harel_bounds: G(n) against the exact normalizing constant
// ---------------------------------------------------------------------------

TEST_CASE("harel G(n) equals the exact normalizing constant of pfqn_ca") {
    // Three shapes: distinct demands, a repeated demand, and a near-balanced set
    const std::vector<std::vector<double> > models = {
        {0.9, 0.4, 0.2}, {0.7, 0.7, 0.3, 0.1}, {0.55, 0.5, 0.45}};
    for (const std::vector<double>& rho : models) {
        const Matrix<double> L = as_column(rho);
        const auto b = pfqn_harel_bounds<double>(rho, 7);
        for (int n = 1; n <= 7; ++n) {
            // TH(n) = G(n-1)/G(n), and pfqn_ca gives G at each population
            const double Gn = pfqn_ca<double>(L, std::vector<int>(1, n)).G;
            const double Gnm1 = pfqn_ca<double>(L, std::vector<int>(1, n - 1)).G;
            CHECK(b.TH[n] == doctest::Approx(Gnm1 / Gn).epsilon(TOL));
        }
    }
}

TEST_CASE("harel TH(n) equals the exact MVA throughput at population n") {
    const std::vector<double> rho = {0.9, 0.4, 0.2};
    const auto b = pfqn_harel_bounds<double>(rho, 7);
    const Matrix<double> L = as_column(rho);
    for (int n = 1; n <= 7; ++n) {
        const auto mva = line::pfqn::pfqn_mva<double>(L, std::vector<int>(1, n));
        CHECK(b.TH[n] == doctest::Approx(mva.XN[0]).epsilon(1e-9));
    }
}

TEST_CASE("harel bounds bracket the exact throughput") {
    const std::vector<std::vector<double> > models = {
        {0.9, 0.4, 0.2}, {0.7, 0.7, 0.3, 0.1}, {0.55, 0.5, 0.45}, {1.2, 0.3}};
    for (const std::vector<double>& rho : models) {
        const Matrix<double> L = as_column(rho);
        for (int N = 2; N <= 7; ++N) {
            const auto b = pfqn_harel_bounds<double>(rho, N);
            const double Xexact = line::pfqn::pfqn_mva<double>(L, std::vector<int>(1, N)).XN[0];
            CHECK(b.LB <= doctest::Approx(Xexact).epsilon(1e-12));
            for (int n = 2; n <= b.maxUB; ++n) CHECK(b.UB[n] >= doctest::Approx(Xexact).epsilon(1e-12));
        }
    }
}

TEST_CASE("harel UB is exact at the extrapolation point n = N") {
    // At n = N the extrapolation factor (N-1)/(n-1) is one, so UB(N) collapses
    // to N / (N/TH(N)) = TH(N), the exact throughput.
    const std::vector<double> rho = {0.9, 0.4, 0.2};
    for (int N = 2; N <= 7; ++N) {
        const auto b = pfqn_harel_bounds<double>(rho, N);
        CHECK(b.UB[N] == doctest::Approx(b.TH[N]).epsilon(TOL));
    }
}

TEST_CASE("harel standalone entry points agree with the batch one") {
    const std::vector<double> rho = {0.8, 0.5, 0.35, 0.1};
    const int N = 6;
    const auto b = pfqn_harel_bounds<double>(rho, N);
    CHECK(pfqn_harel_lb<double>(rho, N) == doctest::Approx(b.LB).epsilon(TOL));
    for (int n = 2; n <= b.maxUB; ++n)
        CHECK(pfqn_harel_ub<double>(rho, N, n) == doctest::Approx(b.UB[n]).epsilon(TOL));
}

TEST_CASE("harel lower bound at N = 1 and N = 2 is the closed form") {
    const std::vector<double> rho = {0.9, 0.4, 0.2};
    const double A1 = 0.9 + 0.4 + 0.2;
    const double A2 = 0.81 + 0.16 + 0.04;
    CHECK(pfqn_harel_lb<double>(rho, 1) == doctest::Approx(1.0 / A1).epsilon(TOL));
    CHECK(pfqn_harel_lb<double>(rho, 2) == doctest::Approx(2.0 / (A1 + A2 / A1)).epsilon(TOL));
}

TEST_CASE("harel UB agrees between exact and double arithmetic") {
    // UB is a field expression and must agree to full double precision. The
    // BATCH entry cannot be used in Rational for N > 2 because it always forms
    // LB, whose (N-1)-st root has no field expression; the standalone UB does
    // not touch LB and stays exact.
    const std::vector<double> rhod = {0.5, 0.25, 0.125};
    std::vector<Rational> rhor(3);
    rhor[0] = Rational(1, 2);
    rhor[1] = Rational(1, 4);
    rhor[2] = Rational(1, 8);
    const auto bd = pfqn_harel_bounds<double>(rhod, 5);
    for (int n = 2; n <= 5; ++n) {
        const Rational ur = pfqn_harel_ub<Rational>(rhor, 5, n);
        CHECK(static_cast<double>(ur) == doctest::Approx(bd.UB[n]).epsilon(1e-14));
    }
    // and the batch entry does refuse the root in exact arithmetic
    CHECK_THROWS_AS(pfqn_harel_bounds<Rational>(rhor, 5, Rational(0), 5),
                    line::UnsupportedError);
    // at N == 2 the exponent is one, so the whole batch stays exact
    const auto b2 = pfqn_harel_bounds<Rational>(rhor, 2, Rational(0), 2);
    const auto d2 = pfqn_harel_bounds<double>(rhod, 2);
    CHECK(static_cast<double>(b2.LB) == doctest::Approx(d2.LB).epsilon(1e-14));
    CHECK(static_cast<double>(b2.UB[2]) == doctest::Approx(d2.UB[2]).epsilon(1e-14));
}

TEST_CASE("harel refuses a think time, a bad population and n above 7") {
    const std::vector<double> rho = {0.9, 0.4};
    CHECK_THROWS_AS(pfqn_harel_bounds<double>(rho, 5, 0.5, 0), line::InputError);
    CHECK_THROWS_AS(pfqn_harel_bounds<double>(rho, 0), line::InputError);
    CHECK_THROWS_AS(pfqn_harel_bounds<double>(rho, 9, 0.0, 8), line::InputError);
    CHECK_THROWS_AS(pfqn_harel_ub<double>(rho, 9, 8), line::InputError);
    CHECK_THROWS_AS(pfqn_harel_ub<double>(rho, 5, 1), line::InputError);
    CHECK_THROWS_AS(pfqn_harel_ub<double>(rho, 5, 6), line::InputError);
    const std::vector<double> bad = {0.9, 0.0};
    CHECK_THROWS_AS(pfqn_harel_bounds<double>(bad, 3), line::InputError);
}

// ---------------------------------------------------------------------------
// pfqn_combine_mi
// ---------------------------------------------------------------------------

TEST_CASE("pfqn_combine_mi sums multiplicities along the mapping") {
    const std::vector<int> mi = {1, 3, 2, 1};
    const std::vector<std::size_t> mapping = {0, 1, 0, 1};
    const std::vector<int> got = pfqn_combine_mi(mi, mapping, 2);
    REQUIRE(got.size() == 2);
    CHECK(got[0] == 3);  // stations 0 and 2
    CHECK(got[1] == 4);  // stations 1 and 3
}

TEST_CASE("pfqn_combine_mi composes with pfqn_unique on replicated demands") {
    // Two identical rows must merge, and unit multiplicities must fold to the
    // group sizes pfqn_unique itself reports.
    Matrix<double> L(3, 2);
    L(0, 0) = 0.4; L(0, 1) = 0.2;
    L(1, 0) = 0.4; L(1, 1) = 0.2;
    L(2, 0) = 0.7; L(2, 1) = 0.1;
    const auto u = line::pfqn::pfqn_unique<double>(L);
    REQUIRE(u.L.rows() == 2);
    const std::vector<int> ones(3, 1);
    const std::vector<int> got = pfqn_combine_mi(ones, u.mapping, u.L.rows());
    REQUIRE(got.size() == u.mi.size());
    for (std::size_t g = 0; g < got.size(); ++g) CHECK(got[g] == u.mi[g]);
}

TEST_CASE("pfqn_combine_mi refuses a mismatched mapping") {
    const std::vector<int> mi = {1, 1};
    CHECK_THROWS_AS(pfqn_combine_mi(mi, std::vector<std::size_t>(3, 0), 1), line::InputError);
    CHECK_THROWS_AS(pfqn_combine_mi(mi, std::vector<std::size_t>(2, 5), 1), line::InputError);
}

// ---------------------------------------------------------------------------
// cache_miss_asy
// ---------------------------------------------------------------------------

TEST_CASE("cache_miss_asy returns a ratio that falls as capacity grows") {
    // gamma is LIST-major here: 2 lists, 6 items, geometrically decaying access
    Matrix<double> gamma(2, 6);
    for (std::size_t k = 0; k < 6; ++k) {
        gamma(0, k) = 1.0 / static_cast<double>(k + 1);
        gamma(1, k) = 0.5 / static_cast<double>(k + 1);
    }
    double prev = 2.0;
    for (int c = 0; c <= 2; ++c) {
        const std::vector<int> m = {c, c};
        const double miss = cache_miss_asy<double>(gamma, m);
        CHECK(miss >= 0.0);
        CHECK(miss <= 1.0);
        CHECK(miss <= prev);
        prev = miss;
    }
}

TEST_CASE("cache_miss_asy is one at zero or negative capacity") {
    Matrix<double> gamma(2, 4);
    for (std::size_t l = 0; l < 2; ++l)
        for (std::size_t k = 0; k < 4; ++k) gamma(l, k) = 1.0 / static_cast<double>(k + 1);
    CHECK(cache_miss_asy<double>(gamma, std::vector<int>(2, 0)) == doctest::Approx(1.0));
    const std::vector<int> neg = {2, -1};
    CHECK(cache_miss_asy<double>(gamma, neg) == doctest::Approx(1.0));
}

TEST_CASE("cache_miss_asy caches everything when capacity covers every item") {
    // With as many slots per list as there are items, no item is ever evicted,
    // so every request hits and the ratio is zero.
    Matrix<double> gamma(1, 4);
    for (std::size_t k = 0; k < 4; ++k) gamma(0, k) = 1.0 / static_cast<double>(k + 1);
    const std::vector<int> m = {4};
    CHECK(cache_miss_asy<double>(gamma, m) == doctest::Approx(0.0).epsilon(1e-12));
}

TEST_CASE("cache_miss_asy refuses a capacity vector of the wrong length") {
    Matrix<double> gamma(2, 4);
    for (std::size_t l = 0; l < 2; ++l)
        for (std::size_t k = 0; k < 4; ++k) gamma(l, k) = 0.5;
    CHECK_THROWS_AS(cache_miss_asy<double>(gamma, std::vector<int>(3, 1)), line::InputError);
}
