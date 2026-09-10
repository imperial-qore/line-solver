/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The cache importance samplers. These are Monte Carlo estimators driven by a
 * pseudo-random stream, so there is no run-for-run reference value to compare
 * against MATLAB: the sampling order differs by construction and only the
 * estimate is comparable. What IS asserted is what the estimator promises --
 * unbiasedness against the exact enumeration of cache_erec / cache_prob_erec
 * on a model small enough for both, and the two conservation laws the hit
 * probabilities must satisfy exactly whatever the sample path.
 *
 * The exact values used as the target were confirmed against MATLAB:
 * cache_erec([0.9 0.4;0.7 0.3;0.5 0.2;0.3 0.1;0.2 0.05;0.1 0.02],[1 2])
 * = 1.2108, and MATLAB cache_is on the same model with 2e5 samples returns
 * 1.21076616, a relative error of 2.8e-5.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/cache/cache_erec.h"
#include "line/api/cache/cache_is.h"
#include "line/api/cache/cache_miss_is.h"
#include "line/api/cache/cache_prob_erec.h"
#include "line/api/cache/cache_prob_is.h"

using line::Matrix;
using namespace line::cache;

namespace {

Matrix<double> gamma6x2() {
    Matrix<double> g(6, 2);
    const double a[6][2] = {{0.9, 0.4}, {0.7, 0.3}, {0.5, 0.2},
                            {0.3, 0.1}, {0.2, 0.05}, {0.1, 0.02}};
    for (std::size_t i = 0; i < 6; ++i)
        for (std::size_t j = 0; j < 2; ++j) g(i, j) = a[i][j];
    return g;
}

const std::vector<int> caps{1, 2};

}  // namespace

TEST_CASE("cache_is is unbiased for the exact normalizing constant") {
    // MATLAB cache_erec on this model gives 1.2108 exactly; the sampler must
    // land on it to within its own Monte Carlo error, which at 2e5 samples is
    // a few parts in 1e4 (MATLAB measures 2.8e-5 on its own stream).
    const double exact = cache_erec(gamma6x2(), caps);
    CHECK(exact == doctest::Approx(1.2108).epsilon(1e-12));
    const auto r = cache_is(gamma6x2(), caps, static_cast<std::size_t>(200000),
                            static_cast<std::uint64_t>(20260721));
    CHECK(r.E == doctest::Approx(exact).epsilon(5e-3));
    CHECK(r.lE == doctest::Approx(std::log(exact)).epsilon(5e-3));
    CHECK(std::exp(r.lE) == doctest::Approx(r.E).epsilon(1e-12));
}

TEST_CASE("cache_is takes the exact route in the degenerate configurations") {
    // Every item cached: one configuration, so the sampler must defer to
    // cache_erec rather than sample.
    Matrix<double> g(2, 1);
    g(0, 0) = 0.7;
    g(1, 0) = 0.4;
    const std::vector<int> m2{2};
    const auto r = cache_is(g, m2, static_cast<std::size_t>(10), static_cast<std::uint64_t>(1));
    CHECK(r.E == doctest::Approx(cache_erec(g, m2)).epsilon(1e-14));

    // No capacity: the empty configuration alone, E = 1.
    const std::vector<int> m0{0};
    const auto r0 = cache_is(g, m0, static_cast<std::size_t>(10), static_cast<std::uint64_t>(1));
    CHECK(r0.E == doctest::Approx(1.0).epsilon(1e-14));
    CHECK(r0.lE == doctest::Approx(0.0).epsilon(1e-14));

    // Fewer items than slots: no valid configuration, E = 0.
    const std::vector<int> m3{3};
    const auto r3 = cache_is(g, m3, static_cast<std::size_t>(10), static_cast<std::uint64_t>(1));
    CHECK(r3.E == 0.0);
}

TEST_CASE("cache_prob_is obeys the two conservation laws exactly") {
    // Whatever the sample path, each item's row is a probability distribution
    // over {miss, list 1, ..., list h}, and each list column sums to its
    // capacity because every drawn configuration fills the list exactly.
    const auto P = cache_prob_is(gamma6x2(), caps, static_cast<std::size_t>(50000),
                                 static_cast<std::uint64_t>(7));
    REQUIRE(P.rows() == 6);
    REQUIRE(P.cols() == 3);
    for (std::size_t i = 0; i < 6; ++i) {
        double s = 0.0;
        for (std::size_t j = 0; j < 3; ++j) {
            CHECK(P(i, j) >= 0.0);
            CHECK(P(i, j) <= 1.0);
            s += P(i, j);
        }
        CHECK(s == doctest::Approx(1.0).epsilon(1e-12));
    }
    for (std::size_t j = 0; j < 2; ++j) {
        double c = 0.0;
        for (std::size_t i = 0; i < 6; ++i) c += P(i, 1 + j);
        CHECK(c == doctest::Approx(static_cast<double>(caps[j])).epsilon(1e-10));
    }
    double miss = 0.0;
    for (std::size_t i = 0; i < 6; ++i) miss += P(i, 0);
    CHECK(miss == doctest::Approx(6.0 - 3.0).epsilon(1e-10));
}

TEST_CASE("cache_prob_is converges to the exact hit probabilities") {
    // Target from cache_prob_erec, which MATLAB confirms entry by entry:
    // its first column is [0.2082920383 0.2659398745 0.3658738024 0.5675586389
    // 0.725140403 0.8671952428]. The sampler is asserted at 1e-2 absolute,
    // which is the error a 2e5-sample self-normalized estimator supports here
    // (MATLAB measures at most 2.0e-3 on its own stream).
    const Matrix<double> exact = cache_prob_erec(gamma6x2(), caps);
    CHECK(exact(0, 0) == doctest::Approx(0.2082920383).epsilon(1e-8));
    CHECK(exact(5, 0) == doctest::Approx(0.8671952428).epsilon(1e-8));
    const auto P = cache_prob_is(gamma6x2(), caps, static_cast<std::size_t>(200000),
                                 static_cast<std::uint64_t>(20260721));
    for (std::size_t i = 0; i < 6; ++i)
        for (std::size_t j = 0; j < 3; ++j) CHECK(std::fabs(P(i, j) - exact(i, j)) < 1e-2);
}

TEST_CASE("cache_miss_is reproduces the exact miss rates within the sampling error") {
    // Two users over the six items; the exact answer comes from the per-item
    // miss probabilities of cache_prob_erec, so this checks the linear map on
    // top of the sampler as well as the sampler itself.
    Matrix<double> lambda(2, 6);
    for (std::size_t k = 0; k < 6; ++k) {
        lambda(0, k) = 1.0 / static_cast<double>(k + 1);
        lambda(1, k) = 0.5;
    }
    const Matrix<double> exact = cache_prob_erec(gamma6x2(), caps);
    const auto r = cache_miss_is(gamma6x2(), caps, lambda, static_cast<std::size_t>(200000),
                                 static_cast<std::uint64_t>(20260721));
    REQUIRE(r.MU.size() == 2);
    REQUIRE(r.MI.size() == 6);
    double eM = 0.0;
    std::vector<double> eMU(2, 0.0);
    for (std::size_t v = 0; v < 2; ++v)
        for (std::size_t k = 0; k < 6; ++k) eMU[v] += lambda(v, k) * exact(k, 0);
    for (std::size_t k = 0; k < 6; ++k) eM += (lambda(0, k) + lambda(1, k)) * exact(k, 0);
    CHECK(r.M == doctest::Approx(eM).epsilon(1e-2));
    CHECK(r.MU[0] == doctest::Approx(eMU[0]).epsilon(1e-2));
    CHECK(r.MU[1] == doctest::Approx(eMU[1]).epsilon(1e-2));
    // the global rate is the sum of the per-item rates, exactly
    double sMI = 0.0;
    for (std::size_t k = 0; k < 6; ++k) sMI += r.MI[k];
    CHECK(r.M == doctest::Approx(sMI).epsilon(1e-12));
    // and the per-item rate is the aggregate rate times the miss probability
    for (std::size_t k = 0; k < 6; ++k)
        CHECK(r.MI[k] == doctest::Approx((lambda(0, k) + lambda(1, k)) * r.pi0[k]).epsilon(1e-12));
}

TEST_CASE("cache_miss_is without request rates returns the mean miss probability") {
    const auto r = cache_miss_is(gamma6x2(), caps, Matrix<double>(),
                                 static_cast<std::size_t>(20000), static_cast<std::uint64_t>(3));
    CHECK(r.MU.empty());
    CHECK(r.MI.empty());
    double s = 0.0;
    for (std::size_t k = 0; k < 6; ++k) s += r.pi0[k];
    CHECK(r.M == doctest::Approx(s / 6.0).epsilon(1e-12));
    // the miss probabilities sum to the number of uncached items
    CHECK(s == doctest::Approx(3.0).epsilon(1e-10));
}
