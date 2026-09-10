/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Decomposition-aggregation inputs (api/da).
 *
 * Oracles: the superposed SCV is a convex combination of the input SCVs, so it
 * lies between their extremes and is invariant to a common rescaling of the
 * rates; the isolated cache preserves the aggregate request rate of every
 * class; and the access factors of the default linear cache have the closed
 * form gamma(i,l) = prod_{j<=l} (lambda p_i) that the tree recursion must
 * reproduce. All of these are exact in the rational instantiation.
 *
 * Reference values were produced by running matlab/src/api/da/*.m under
 * MATLAB R2025 (lineStart; da_traffic_superpos / da_cache_isolate). Both
 * routines are finite field computations, so the port is compared to MATLAB at
 * the double-rounding level (1e-15 relative), not at a method tolerance:
 * neither has one.
 */
#include <cmath>
#include <limits>
#include <vector>

#include "doctest.h"
#include "line/api/da/da_cache_isolate.h"
#include "line/api/da/da_fpi.h"
#include "line/api/da/da_traffic_superpos.h"

using line::Rational;
using line::Real50;
using line::da::CacheParam;
using line::da::da_cache_isolate;
using line::da::da_traffic_superpos;

TEST_CASE("da_traffic_superpos is the rate-weighted SCV mixture, exactly") {
    // lambda = [1,2,3], a2 = [1/2,1,2] -> (1/2 + 2 + 6)/6 = 17/12.
    std::vector<Rational> lambda;
    lambda.push_back(Rational(1));
    lambda.push_back(Rational(2));
    lambda.push_back(Rational(3));
    std::vector<Rational> a2;
    a2.push_back(Rational(1, 2));
    a2.push_back(Rational(1));
    a2.push_back(Rational(2));
    CHECK(da_traffic_superpos(lambda, a2) == Rational(17, 12));

    // Convexity: the mixture lies between the smallest and largest input SCV.
    const Rational d2 = da_traffic_superpos(lambda, a2);
    CHECK(d2 >= Rational(1, 2));
    CHECK(d2 <= Rational(2));

    // Scale invariance: multiplying every rate by 5 cannot move the mixture.
    std::vector<Rational> scaled(lambda);
    for (std::size_t i = 0; i < scaled.size(); ++i) scaled[i] *= Rational(5);
    CHECK(da_traffic_superpos(scaled, a2) == d2);

    // Identical SCVs superpose to themselves whatever the rates.
    std::vector<Rational> same(3, Rational(7, 3));
    CHECK(da_traffic_superpos(lambda, same) == Rational(7, 3));
}

TEST_CASE("da_traffic_superpos matches MATLAB") {
    // MATLAB: da_traffic_superpos([1 2 3],[0.5 1 2]) = 1.4166666666666667
    std::vector<double> lambda;
    lambda.push_back(1.0);
    lambda.push_back(2.0);
    lambda.push_back(3.0);
    std::vector<double> a2;
    a2.push_back(0.5);
    a2.push_back(1.0);
    a2.push_back(2.0);
    CHECK(da_traffic_superpos(lambda, a2) == doctest::Approx(1.4166666666666667).epsilon(1e-15));

    // MATLAB drops non-finite rates: da_traffic_superpos([1 Inf 3],[0.5 1 2])
    // = 1.625, the mixture of flows 1 and 3 alone.
    std::vector<double> withInf(lambda);
    withInf[1] = std::numeric_limits<double>::infinity();
    CHECK(da_traffic_superpos(withInf, a2) == doctest::Approx(1.625).epsilon(1e-15));

    // The mask is a no-op in an exact instantiation: rationals are all finite.
    std::vector<Real50> lr, ar;
    for (std::size_t i = 0; i < 3; ++i) {
        lr.push_back(Real50(lambda[i]));
        ar.push_back(Real50(a2[i]));
    }
    CHECK(std::fabs(static_cast<double>(da_traffic_superpos(lr, ar)) - 1.4166666666666667) < 1e-15);
}

TEST_CASE("da_traffic_superpos rejects malformed input") {
    std::vector<double> one(1, 1.0);
    std::vector<double> two(2, 1.0);
    CHECK_THROWS_AS(da_traffic_superpos(one, two), line::InputError);
    CHECK_THROWS_AS(da_traffic_superpos(std::vector<double>(), std::vector<double>()),
                    line::InputError);
    std::vector<double> zero(2, 0.0);
    CHECK_THROWS_AS(da_traffic_superpos(zero, two), line::NumericError);
}

/** Two-list, three-item, single-class cache with the default linear routing. */
static CacheParam<Rational> exact_cache() {
    CacheParam<Rational> ch;
    ch.itemcap.push_back(1);
    ch.itemcap.push_back(1);
    ch.nitems = 3;
    std::vector<Rational> pread;
    pread.push_back(Rational(1, 2));
    pread.push_back(Rational(3, 10));
    pread.push_back(Rational(1, 5));
    ch.pread.push_back(pread);
    return ch;
}

TEST_CASE("da_cache_isolate conserves the aggregate request rate") {
    const CacheParam<Rational> ch = exact_cache();
    std::vector<Rational> lambda(1, Rational(2));
    const line::da::CacheIsolateResult<Rational> r = da_cache_isolate(ch, lambda);

    // CONSERVATION: the isolated cache is driven by exactly the rate the
    // network hands it. Summing the per-item rates at any list position must
    // return lambda, exactly, because pread is a probability distribution.
    REQUIRE(r.lambda_cache.size() == 1u);
    REQUIRE(r.lambda_cache[0].rows() == 3u);
    REQUIRE(r.lambda_cache[0].cols() == 3u);  // h + 1 positions
    for (std::size_t l = 0; l < 3; ++l) {
        Rational total(0);
        for (std::size_t k = 0; k < 3; ++k) total += r.lambda_cache[0](k, l);
        CHECK(total == Rational(2));
    }

    // A class that does not read the cache (empty pread, MATLAB's NaN) draws
    // no rate at all, so the conserved total is unchanged by adding it.
    CacheParam<Rational> ch2 = ch;
    ch2.pread.push_back(std::vector<Rational>());
    std::vector<Rational> lambda2;
    lambda2.push_back(Rational(2));
    lambda2.push_back(Rational(5));
    const line::da::CacheIsolateResult<Rational> r2 = da_cache_isolate(ch2, lambda2);
    REQUIRE(r2.lambda_cache.size() == 2u);
    for (std::size_t k = 0; k < 3; ++k)
        for (std::size_t l = 0; l < 3; ++l) CHECK(r2.lambda_cache[1](k, l) == Rational(0));
    for (std::size_t l = 0; l < 3; ++l) {
        Rational total(0);
        for (std::size_t v = 0; v < 2; ++v)
            for (std::size_t k = 0; k < 3; ++k) total += r2.lambda_cache[v](k, l);
        CHECK(total == Rational(2));
    }
}

TEST_CASE("da_cache_isolate builds the default linear cache routing") {
    const CacheParam<Rational> ch = exact_cache();
    std::vector<Rational> lambda(1, Rational(2));
    const line::da::CacheIsolateResult<Rational> r = da_cache_isolate(ch, lambda);

    REQUIRE(r.Rcost.size() == 1u);
    REQUIRE(r.Rcost[0].size() == 3u);
    // R = [0 1 0; 0 0 1; 0 0 1]: list l feeds list l+1, the last list holds.
    for (std::size_t k = 0; k < 3; ++k) {
        const line::Matrix<Rational>& R = r.Rcost[0][k];
        REQUIRE(R.rows() == 3u);
        REQUIRE(R.cols() == 3u);
        CHECK(R(0, 1) == Rational(1));
        CHECK(R(1, 2) == Rational(1));
        CHECK(R(2, 2) == Rational(1));
        Rational total(0);
        for (std::size_t a = 0; a < 3; ++a)
            for (std::size_t b = 0; b < 3; ++b) total += R(a, b);
        CHECK(total == Rational(3));  // exactly the three ones above
    }
}

TEST_CASE("da_cache_isolate access factors have the exact linear-cache form") {
    // For the default linear cache the path to list l is 0 -> 1 -> ... -> l,
    // and every edge carries the full item rate, so gamma(i,l) = (lambda p_i)^l.
    const CacheParam<Rational> ch = exact_cache();
    std::vector<Rational> lambda(1, Rational(2));
    const line::da::CacheIsolateResult<Rational> r = da_cache_isolate(ch, lambda);

    REQUIRE(r.gamma.rows() == 3u);
    REQUIRE(r.gamma.cols() == 2u);
    for (std::size_t k = 0; k < 3; ++k) {
        const Rational rate = Rational(2) * ch.pread[0][k];
        CHECK(r.gamma(k, 0) == rate);
        CHECK(r.gamma(k, 1) == rate * rate);
    }
    // MATLAB: da_cache_isolate on the same cache returns
    // gamma = [1 1; 0.6 0.36; 0.4 0.16].
    CHECK(r.gamma(0, 0) == Rational(1));
    CHECK(r.gamma(0, 1) == Rational(1));
    CHECK(r.gamma(1, 0) == Rational(3, 5));
    CHECK(r.gamma(1, 1) == Rational(9, 25));
    CHECK(r.gamma(2, 0) == Rational(2, 5));
    CHECK(r.gamma(2, 1) == Rational(4, 25));
}

TEST_CASE("da_cache_isolate matches MATLAB in double arithmetic") {
    // MATLAB reference, ch.itemcap=[1 1], ch.nitems=3, pread={[0.5 0.3 0.2]},
    // lambda=2:  gamma = [1 1; 0.59999999999999998 0.35999999999999999;
    //                     0.40000000000000002 0.16000000000000003]
    CacheParam<double> ch;
    ch.itemcap.push_back(1);
    ch.itemcap.push_back(1);
    ch.nitems = 3;
    std::vector<double> pread;
    pread.push_back(0.5);
    pread.push_back(0.3);
    pread.push_back(0.2);
    ch.pread.push_back(pread);
    std::vector<double> lambda(1, 2.0);

    const line::da::CacheIsolateResult<double> r = da_cache_isolate(ch, lambda);
    // Products of doubles only, so agreement is at the rounding level.
    CHECK(r.gamma(0, 0) == doctest::Approx(1.0).epsilon(1e-15));
    CHECK(r.gamma(0, 1) == doctest::Approx(1.0).epsilon(1e-15));
    CHECK(r.gamma(1, 0) == doctest::Approx(0.59999999999999998).epsilon(1e-15));
    CHECK(r.gamma(1, 1) == doctest::Approx(0.35999999999999999).epsilon(1e-15));
    CHECK(r.gamma(2, 0) == doctest::Approx(0.40000000000000002).epsilon(1e-15));
    CHECK(r.gamma(2, 1) == doctest::Approx(0.16000000000000003).epsilon(1e-15));
    // MATLAB: size(lambda_cache) = 1 3 3, and every list position carries the
    // same per-item rate lambda * pread.
    CHECK(r.lambda_cache[0](0, 0) == doctest::Approx(1.0).epsilon(1e-15));
    CHECK(r.lambda_cache[0](0, 2) == doctest::Approx(1.0).epsilon(1e-15));
    CHECK(r.lambda_cache[0](1, 0) == doctest::Approx(0.6).epsilon(1e-15));
    CHECK(r.lambda_cache[0](2, 2) == doctest::Approx(0.4).epsilon(1e-15));
}

TEST_CASE("da_cache_isolate honours an explicit access-cost matrix") {
    // A cache whose second list is fed directly from "not cached": the path to
    // list 2 is 0 -> 2 and gamma(i,2) becomes the single-edge rate, not its
    // square. This is the branch the default routing never exercises.
    CacheParam<Rational> ch = exact_cache();
    line::Matrix<Rational> R(3, 3, Rational(0));
    R(0, 1) = Rational(1, 2);
    R(0, 2) = Rational(1, 2);
    R(1, 1) = Rational(1);
    R(2, 2) = Rational(1);
    ch.accost.assign(1, std::vector<line::Matrix<Rational>>(3, R));
    std::vector<Rational> lambda(1, Rational(2));

    const line::da::CacheIsolateResult<Rational> r = da_cache_isolate(ch, lambda);
    for (std::size_t k = 0; k < 3; ++k) {
        const Rational rate = Rational(2) * ch.pread[0][k];
        CHECK(r.gamma(k, 0) == rate * Rational(1, 2));
        CHECK(r.gamma(k, 1) == rate * Rational(1, 2));
    }
}

TEST_CASE("da_cache_isolate rejects malformed cache parameters") {
    CacheParam<Rational> ch = exact_cache();
    std::vector<Rational> lambda(1, Rational(2));

    CacheParam<Rational> noLists = ch;
    noLists.itemcap.clear();
    CHECK_THROWS_AS(da_cache_isolate(noLists, lambda), line::InputError);

    CacheParam<Rational> noItems = ch;
    noItems.nitems = 0;
    CHECK_THROWS_AS(da_cache_isolate(noItems, lambda), line::InputError);

    CacheParam<Rational> shortPread = ch;
    shortPread.pread[0].pop_back();
    CHECK_THROWS_AS(da_cache_isolate(shortPread, lambda), line::InputError);

    CacheParam<Rational> badAccost = ch;
    badAccost.accost.assign(1, std::vector<line::Matrix<Rational>>(3, line::Matrix<Rational>(2, 2)));
    CHECK_THROWS_AS(da_cache_isolate(badAccost, lambda), line::InputError);
}

TEST_CASE("da_cache_isolate composes with the da_fpi driver") {
    // The decomposition loop the two are meant to form: the cache rate is
    // updated from the hit rate it produces, here through the closed-form
    // one-list miss ratio, and da_fpi drives it to a fixed point. This only
    // checks that the two pieces compose and that the fixed point is the one
    // the map's algebra predicts.
    CacheParam<double> ch;
    ch.itemcap.push_back(1);
    ch.nitems = 2;
    std::vector<double> pread;
    pread.push_back(0.5);
    pread.push_back(0.5);
    ch.pread.push_back(pread);

    struct Step {
        CacheParam<double> ch;
        std::pair<std::vector<double>, std::vector<double>> operator()(const std::vector<double>& x,
                                                                      std::size_t) const {
            const line::da::CacheIsolateResult<double> r = da_cache_isolate(ch, x);
            double g = 0.0;
            for (std::size_t k = 0; k < r.gamma.rows(); ++k) g += r.gamma(k, 0);
            std::vector<double> xnew(1, 1.0 + 0.5 * g / (1.0 + g));
            return std::make_pair(xnew, x);
        }
    };
    Step step;
    step.ch = ch;

    line::da::FpiOptions opts;
    opts.iter_tol = 1e-12;
    const line::da::FpiResult<double> fp =
        line::da::da_fpi<double>(step, std::vector<double>(1, 1.0), opts);
    CHECK(fp.converged);
    // At the fixed point g = sum_k lambda p_k = lambda, so lambda solves
    // lambda = 1 + lambda / (2 (1 + lambda)), i.e. 2 lambda^2 - lambda - 2 = 0.
    const double expected = (1.0 + std::sqrt(1.0 + 16.0)) / 4.0;
    CHECK(fp.x[0] == doctest::Approx(expected).epsilon(1e-9));  // iter_tol = 1e-12
}
