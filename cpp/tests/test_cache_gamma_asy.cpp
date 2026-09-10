/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The two cache routines that had no MATLAB source until the api/cache port
 * audit: cache_gamma, the access factors over a general access GRAPH, and
 * cache_miss_asy, the rank-threshold asymptotic miss ratio. Both now exist in
 * all four codebases; these cases are the ones the MATLAB and Python twins are
 * pinned against.
 *
 * Oracles, in order of strength:
 *   1. Closed forms. On a linear chain the BFS path is forced and the access
 *      factor is a product that can be written down; for a single list of size
 *      m the rank rule puts exactly the m most popular items in cache, so the
 *      miss ratio is the tail mass, independent of the iteration.
 *   2. Structural invariants: an unreachable list has a zero access factor,
 *      and a degenerate capacity misses every request.
 *   3. Cross-codebase agreement with the MATLAB and Python twins on the
 *      multi-user diamond graph, where the BFS tie-break in node order is what
 *      selects the path.
 */
#include <vector>

#include "doctest.h"
#include "line/api/cache/cache_gamma.h"
#include "line/api/cache/cache_miss_asy.h"

using line::Matrix;
using line::cache::cache_gamma;
using line::cache::cache_miss_asy;

namespace {

constexpr double TOL = 1e-12;

/** (h+1) x (h+1) adjacency, scaled, from a row-major list of 0/1 entries. */
Matrix<double> graph4(const int* a, double scale) {
    Matrix<double> R(4, 4, 0.0);
    for (std::size_t i = 0; i < 4; ++i)
        for (std::size_t j = 0; j < 4; ++j) R(i, j) = a[i * 4 + j] * scale;
    return R;
}

}  // namespace

TEST_CASE("cache_miss_asy: single list of size m misses exactly the tail mass") {
    // gamma is LIST-major here, so a single list is one row.
    Matrix<double> g(1, 5, 0.0);
    const double p[5] = {0.4, 0.3, 0.15, 0.1, 0.05};
    for (std::size_t k = 0; k < 5; ++k) g(0, k) = p[k];

    // The rank rule caches the two most popular items, so every request for
    // items 3..5 misses and no request for items 1..2 does.
    CHECK(cache_miss_asy(g, std::vector<int>{2}) == doctest::Approx(0.30).epsilon(TOL));
    CHECK(cache_miss_asy(g, std::vector<int>{1}) == doctest::Approx(0.60).epsilon(TOL));
    CHECK(cache_miss_asy(g, std::vector<int>{4}) == doctest::Approx(0.05).epsilon(TOL));
    // More slots than competitors: nothing can miss.
    CHECK(cache_miss_asy(g, std::vector<int>{5}) == doctest::Approx(0.0).epsilon(TOL));
}

TEST_CASE("cache_miss_asy: two lists, and the degenerate capacity") {
    Matrix<double> g(2, 4, 0.0);
    const double r0[4] = {0.5, 0.3, 0.15, 0.05};
    const double r1[4] = {0.4, 0.4, 0.1, 0.1};
    for (std::size_t k = 0; k < 4; ++k) {
        g(0, k) = r0[k];
        g(1, k) = r1[k];
    }
    // MATLAB cache_miss_asy(g, [1 2]) and the Python twin.
    CHECK(cache_miss_asy(g, std::vector<int>{1, 2}) == doctest::Approx(0.35).epsilon(1e-9));
    // No capacity at all: every request misses.
    CHECK(cache_miss_asy(g, std::vector<int>{0, 0}) == doctest::Approx(1.0).epsilon(TOL));
    // A negative capacity is degenerate in the same way.
    CHECK(cache_miss_asy(g, std::vector<int>{-1, 2}) == doctest::Approx(1.0).epsilon(TOL));
}

TEST_CASE("cache_gamma: linear chain has a closed-form access factor") {
    // 1 -> 2 -> 3 -> 4 with a self-loop on the top list, one user, one item.
    static const int chain[16] = {0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1};
    std::vector<Matrix<double>> lambda(1, Matrix<double>(1, 4, 0.0));
    const double lam[4] = {2, 3, 5, 7};
    for (std::size_t t = 0; t < 4; ++t) lambda[0](0, t) = lam[t];
    std::vector<std::vector<Matrix<double>>> R(1);
    R[0].push_back(graph4(chain, 1.0));

    const line::cache::CacheGammaGraphResult<double> res = cache_gamma(lambda, R);
    CHECK(res.u == 1);
    CHECK(res.n == 1);
    CHECK(res.h == 3);
    // Column 0 is the one-node path, so it is the miss-node rate alone;
    // column 1 multiplies by lambda(0) R(0,1); column 2 by lambda(1) R(1,2).
    CHECK(res.gamma(0, 0) == doctest::Approx(2.0).epsilon(TOL));
    CHECK(res.gamma(0, 1) == doctest::Approx(2.0 * 2.0).epsilon(TOL));
    CHECK(res.gamma(0, 2) == doctest::Approx(2.0 * 2.0 * 3.0).epsilon(TOL));
}

TEST_CASE("cache_gamma: an unreachable list has a zero access factor") {
    // 1 -> 2 -> 4, so node 3 is never entered.
    static const int broken[16] = {0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<Matrix<double>> lambda(1, Matrix<double>(1, 4, 0.0));
    const double lam[4] = {2, 3, 5, 7};
    for (std::size_t t = 0; t < 4; ++t) lambda[0](0, t) = lam[t];
    std::vector<std::vector<Matrix<double>>> R(1);
    R[0].push_back(graph4(broken, 1.0));

    const line::cache::CacheGammaGraphResult<double> res = cache_gamma(lambda, R);
    CHECK(res.gamma(0, 0) == doctest::Approx(2.0).epsilon(TOL));
    CHECK(res.gamma(0, 1) == doctest::Approx(4.0).epsilon(TOL));
    CHECK(res.gamma(0, 2) == doctest::Approx(0.0).epsilon(TOL));
}

TEST_CASE("cache_gamma: multi-user diamond graph agrees with MATLAB and Python") {
    // 1 -> {2,3} -> 4, so lists 2 and 3 are siblings and list 4 has two
    // parents: a structure cache_gamma_lp rejects outright.
    static const int diamond[16] = {0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1};
    std::vector<Matrix<double>> lambda(2, Matrix<double>(2, 4, 0.0));
    const double l00[4] = {1, 2, 3, 4};
    const double l01[4] = {9, 1, 2, 3};
    const double l10[4] = {5, 6, 7, 8};
    const double l11[4] = {4, 5, 6, 7};
    for (std::size_t t = 0; t < 4; ++t) {
        lambda[0](0, t) = l00[t];
        lambda[0](1, t) = l01[t];
        lambda[1](0, t) = l10[t];
        lambda[1](1, t) = l11[t];
    }
    std::vector<std::vector<Matrix<double>>> R(2);
    for (std::size_t v = 0; v < 2; ++v)
        for (std::size_t i = 0; i < 2; ++i) R[v].push_back(graph4(diamond, 0.5 * (v + 1)));

    const line::cache::CacheGammaGraphResult<double> res = cache_gamma(lambda, R);
    CHECK(res.gamma(0, 0) == doctest::Approx(6.0).epsilon(TOL));
    CHECK(res.gamma(0, 1) == doctest::Approx(33.0).epsilon(TOL));
    CHECK(res.gamma(0, 2) == doctest::Approx(33.0).epsilon(TOL));
    CHECK(res.gamma(1, 0) == doctest::Approx(13.0).epsilon(TOL));
    CHECK(res.gamma(1, 1) == doctest::Approx(110.5).epsilon(TOL));
    CHECK(res.gamma(1, 2) == doctest::Approx(110.5).epsilon(TOL));
}
