/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Mean value analysis of caches and the mean-field drift. Oracles:
 *   1. cache_prob_erec, an independent (normalizing-constant) route to the
 *      same hit probabilities; in exact arithmetic the two must agree as
 *      values, not merely to a tolerance.
 *   2. The identity x(l) = E(m - e_l)/E(m) between the MVA throughput and the
 *      cache_erec constants.
 *   3. Conservation laws asserted exactly: sum_j pij(k,j) + pi0(k) = 1 for MVA,
 *      sum_k (1 - Mk(k)) = sum(m) for cache_mva_miss, and a drift that sums to
 *      zero over the levels of every item.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/cache/cache_erec.h"
#include "line/api/cache/cache_mva.h"
#include "line/api/cache/cache_mva_miss.h"
#include "line/api/cache/cache_prob_erec.h"
#include "line/api/cache/cache_rrm_meanfield_ode.h"

using line::Matrix;
using line::Rational;
using line::Real50;
using line::cache::cache_erec;
using line::cache::cache_mva;
using line::cache::cache_mva_miss;
using line::cache::cache_prob_erec;
using line::cache::cache_rrm_meanfield_ode;

namespace {

constexpr double TOL = 1e-9;

template <class T>
Matrix<T> gamma_3x2() {
    Matrix<T> g(3, 2);
    g(0, 0) = line::num_traits<T>::from_rational(1, 2);
    g(0, 1) = line::num_traits<T>::from_rational(1, 3);
    g(1, 0) = line::num_traits<T>::from_rational(1, 4);
    g(1, 1) = line::num_traits<T>::from_rational(1, 5);
    g(2, 0) = line::num_traits<T>::from_rational(1, 6);
    g(2, 1) = line::num_traits<T>::from_rational(1, 7);
    return g;
}

}  // namespace

TEST_CASE("cache_mva single list single slot matches the closed form") {
    // Two items, one slot: item k is cached with probability g_k/(g1+g2) and
    // the list throughput is 1/(g1+g2).
    Matrix<Rational> g(2, 1);
    g(0, 0) = Rational(1, 2);
    g(1, 0) = Rational(1, 3);
    const auto r = cache_mva(g, std::vector<int>{1});
    CHECK(r.x[0] == Rational(6, 5));  // 1 / (1/2 + 1/3)
    CHECK(r.pij(0, 0) == Rational(3, 5));
    CHECK(r.pij(1, 0) == Rational(2, 5));
    CHECK(r.pi[0] == Rational(3, 5));
    CHECK(r.pi0[0] == Rational(2, 5));
    CHECK(r.u(0, 0) == Rational(1, 2) * Rational(6, 5));

    // A single item and a single slot is always cached.
    Matrix<Rational> one(1, 1);
    one(0, 0) = Rational(7, 11);
    const auto r1 = cache_mva(one, std::vector<int>{1});
    CHECK(r1.pi[0] == Rational(1));
    CHECK(r1.pi0[0] == Rational(0));
    CHECK(r1.x[0] == Rational(11, 7));
}

TEST_CASE("cache_mva reproduces cache_prob_erec exactly") {
    // sum(m) must not exceed the item count, else no placement is admissible.
    const std::vector<std::vector<int>> caps{{1, 1}, {2, 1}, {1, 2}};
    for (const std::vector<int>& m : caps) {
        const auto mva = cache_mva(gamma_3x2<Rational>(), m);
        const Matrix<Rational> p = cache_prob_erec(gamma_3x2<Rational>(), m);
        for (std::size_t k = 0; k < 3; ++k) {
            for (std::size_t l = 0; l < 2; ++l) CHECK(mva.pij(k, l) == p(k, 1 + l));
            CHECK(mva.pi0[k] == p(k, 0));
        }
    }
}

TEST_CASE("cache_mva throughput equals the ratio of cache_erec constants") {
    // Equating the MVA hit probability gamma(k,l) pi0(k;m-e_l) x(l) with the
    // constant-ratio form m(l) gamma(k,l) E(gamma\k,m-e_l)/E(m) gives
    // x(l) = m(l) E(m-e_l)/E(m). Verified against MATLAB cache_mva on this
    // very model: x = [7.8360655737704956, 5.1639344262295079].
    const std::vector<int> m{2, 1};
    const auto r = cache_mva(gamma_3x2<Rational>(), m);
    const Rational E = cache_erec(gamma_3x2<Rational>(), m);
    for (std::size_t l = 0; l < 2; ++l) {
        std::vector<int> ml = m;
        ml[l] -= 1;
        CHECK(r.x[l] == Rational(m[l]) * cache_erec(gamma_3x2<Rational>(), ml) / E);
    }
    CHECK(static_cast<double>(r.x[0]) == doctest::Approx(7.8360655737704956).epsilon(1e-12));
    CHECK(static_cast<double>(r.x[1]) == doctest::Approx(5.1639344262295079).epsilon(1e-12));
}

TEST_CASE("an over-committed cache is refused, not silently answered") {
    // Four slots for three items: no admissible placement, E = 0.
    const std::vector<int> m{2, 2};
    CHECK(cache_erec(gamma_3x2<Rational>(), m) == Rational(0));
    CHECK_THROWS_AS(cache_prob_erec(gamma_3x2<Rational>(), m), line::NumericError);
    CHECK_THROWS_AS(cache_mva(gamma_3x2<Rational>(), m), line::NumericError);
}

TEST_CASE("cache_mva conserves probability exactly and fills every list") {
    const std::vector<int> m{2, 1};
    const auto r = cache_mva(gamma_3x2<Rational>(), m);
    Rational occ0(0), occ1(0);
    for (std::size_t k = 0; k < 3; ++k) {
        CHECK(r.pij(k, 0) + r.pij(k, 1) + r.pi0[k] == Rational(1));
        occ0 += r.pij(k, 0);
        occ1 += r.pij(k, 1);
    }
    CHECK(occ0 == Rational(2));
    CHECK(occ1 == Rational(1));
    // Documented reference defect: E is never computed and is always one.
    CHECK(r.E == Rational(1));
}

TEST_CASE("cache_mva double and Real50 agree with the exact rational") {
    const std::vector<int> m{2, 1};
    const auto q = cache_mva(gamma_3x2<Rational>(), m);
    const auto d = cache_mva(gamma_3x2<double>(), m);
    const auto rr = cache_mva(gamma_3x2<Real50>(), m);
    for (std::size_t k = 0; k < 3; ++k) {
        CHECK(d.pi0[k] == doctest::Approx(static_cast<double>(q.pi0[k])).epsilon(TOL));
        CHECK(static_cast<double>(rr.pi0[k]) ==
              doctest::Approx(static_cast<double>(q.pi0[k])).epsilon(TOL));
        for (std::size_t l = 0; l < 2; ++l)
            CHECK(d.pij(k, l) == doctest::Approx(static_cast<double>(q.pij(k, l))).epsilon(TOL));
    }
    for (std::size_t l = 0; l < 2; ++l)
        CHECK(d.x[l] == doctest::Approx(static_cast<double>(q.x[l])).epsilon(TOL));
}

TEST_CASE("cache_mva_miss single list matches the closed form") {
    // h = 1, m = [1], R all ones: w(k) = p(k), x = 1/sum(p) and
    // Mk(k) = 1 - p(k)/sum(p).
    std::vector<Rational> p{Rational(1, 2), Rational(1, 3), Rational(1, 6)};
    Matrix<Rational> R(1, 3, Rational(1));
    const auto r = cache_mva_miss(p, std::vector<int>{1}, R);
    CHECK(r.Mk[0] == Rational(1, 2));
    CHECK(r.Mk[1] == Rational(2, 3));
    CHECK(r.Mk[2] == Rational(5, 6));
    // M = sum p(k) Mk(k) = 1/4 + 2/9 + 5/36 = 22/36 = 11/18.
    CHECK(r.M == Rational(11, 18));
}

TEST_CASE("cache_mva_miss fills exactly sum(m) slots") {
    std::vector<Rational> p{Rational(2, 5), Rational(3, 10), Rational(1, 5), Rational(1, 10)};
    Matrix<Rational> R(2, 4, Rational(1));
    const std::vector<std::vector<int>> caps{{1, 0}, {1, 1}, {2, 1}};
    for (const std::vector<int>& m : caps) {
        const auto r = cache_mva_miss(p, m, R);
        Rational filled(0);
        for (std::size_t k = 0; k < p.size(); ++k) filled += Rational(1) - r.Mk[k];
        long mt = 0;
        for (int v : m) mt += v;
        CHECK(filled == Rational(mt));
    }
}

TEST_CASE("cache_mva_miss empty cache misses everything") {
    std::vector<double> p{0.5, 0.5};
    Matrix<double> R(1, 2, 1.0);
    const auto r = cache_mva_miss(p, std::vector<int>{0}, R);
    CHECK(r.Mk[0] == 1.0);
    CHECK(r.Mk[1] == 1.0);
    CHECK(r.M == doctest::Approx(1.0).epsilon(TOL));
}

TEST_CASE("cache_mva_miss double agrees with the exact rational") {
    std::vector<Rational> pq{Rational(1, 2), Rational(1, 3), Rational(1, 6)};
    std::vector<double> pd{0.5, 1.0 / 3.0, 1.0 / 6.0};
    Matrix<Rational> Rq(2, 3, Rational(3, 4));
    Matrix<double> Rd(2, 3, 0.75);
    const std::vector<int> m{1, 1};
    const auto q = cache_mva_miss(pq, m, Rq);
    const auto d = cache_mva_miss(pd, m, Rd);
    CHECK(d.M == doctest::Approx(static_cast<double>(q.M)).epsilon(TOL));
    for (std::size_t k = 0; k < 3; ++k)
        CHECK(d.Mk[k] == doctest::Approx(static_cast<double>(q.Mk[k])).epsilon(TOL));
}

TEST_CASE("cache_rrm_meanfield_ode drift sums to zero over the levels of each item") {
    // Three items, two lists; an arbitrary interior point of the simplex.
    const std::size_t n = 3, h = 2;
    std::vector<Rational> lambda{Rational(1, 2), Rational(1, 3), Rational(1, 5)};
    std::vector<int> m{1, 2};
    std::vector<Rational> x(n * (h + 1));
    const Rational vals[3][3] = {{Rational(1, 2), Rational(1, 4), Rational(1, 4)},
                                 {Rational(1, 3), Rational(1, 3), Rational(1, 3)},
                                 {Rational(1, 6), Rational(1, 2), Rational(1, 3)}};
    for (std::size_t k = 0; k < n; ++k)
        for (std::size_t s = 0; s <= h; ++s) x[k + s * n] = vals[k][s];

    const std::vector<Rational> dx = cache_rrm_meanfield_ode(x, lambda, m);
    for (std::size_t k = 0; k < n; ++k) {
        Rational s(0);
        for (std::size_t l = 0; l <= h; ++l) s += dx[k + l * n];
        CHECK(s == Rational(0));  // each item stays on its simplex
    }
}

TEST_CASE("cache_rrm_meanfield_ode single item single list has the expected fixed point") {
    // One item, one slot: dx1/dt = lam x0 (1 - x1). The only fixed point on
    // x0 + x1 = 1 is x1 = 1, the item permanently cached.
    std::vector<Rational> lambda{Rational(3, 7)};
    std::vector<int> m{1};
    std::vector<Rational> cached{Rational(0), Rational(1)};
    const std::vector<Rational> d0 = cache_rrm_meanfield_ode(cached, lambda, m);
    CHECK(d0[0] == Rational(0));
    CHECK(d0[1] == Rational(0));

    // Away from it the drift pushes the item into the cache.
    std::vector<Rational> half{Rational(1, 2), Rational(1, 2)};
    const std::vector<Rational> d1 = cache_rrm_meanfield_ode(half, lambda, m);
    CHECK(d1[1] == Rational(3, 7) * Rational(1, 2) * Rational(1, 2));  // lam x0 (1 - x1)
    CHECK(d1[0] == -d1[1]);
}

TEST_CASE("cache_rrm_meanfield_ode double agrees with the exact rational") {
    std::vector<Rational> lq{Rational(1, 2), Rational(1, 4)};
    std::vector<double> ld{0.5, 0.25};
    std::vector<int> m{1};
    std::vector<Rational> xq{Rational(3, 4), Rational(1, 2), Rational(1, 4), Rational(1, 2)};
    std::vector<double> xd{0.75, 0.5, 0.25, 0.5};
    const std::vector<Rational> dq = cache_rrm_meanfield_ode(xq, lq, m);
    const std::vector<double> dd = cache_rrm_meanfield_ode(xd, ld, m);
    for (std::size_t i = 0; i < dq.size(); ++i)
        CHECK(dd[i] == doctest::Approx(static_cast<double>(dq[i])).epsilon(TOL));
}
