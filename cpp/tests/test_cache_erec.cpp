/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Exact recursive cache analysis. Oracles, in order of strength:
 *   1. Hand-computed closed forms. For a single list the normalizing constant
 *      collapses to E(gamma,[m]) = m! e_m(gamma), the m-th elementary
 *      symmetric polynomial, which gives exact values for every quantity in
 *      the family.
 *   2. Conservation laws that hold identically, asserted with == in the exact
 *      instantiation: sum_i prob(i,1+j) = m(j) (exactly m(j) items sit in list
 *      j) and sum_j prob(i,j) = 1.
 *   3. Agreement between cache_miss's pi0 and cache_prob_erec's miss column,
 *      two independent routes to the same probability.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/cache/cache_erec.h"
#include "line/api/cache/cache_gamma_lp.h"
#include "line/api/cache/cache_miss.h"
#include "line/api/cache/cache_prob_erec.h"

using line::Matrix;
using line::Rational;
using line::Real50;
using line::cache::cache_erec;
using line::cache::cache_gamma_lp;
using line::cache::cache_miss;
using line::cache::cache_prob_erec;

namespace {

constexpr double TOL = 1e-9;

/** Single-list model: gamma = [1/2, 1/3, 1/5]. */
template <class T>
Matrix<T> gamma_1list() {
    Matrix<T> g(3, 1);
    g(0, 0) = line::num_traits<T>::from_rational(1, 2);
    g(1, 0) = line::num_traits<T>::from_rational(1, 3);
    g(2, 0) = line::num_traits<T>::from_rational(1, 5);
    return g;
}

/** Two-list model, three items, no symmetry. */
template <class T>
Matrix<T> gamma_2list() {
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

TEST_CASE("cache_erec single list equals m! times the elementary symmetric polynomial") {
    // gamma = [1/2, 1/3, 1/5]:
    //   e1 = 31/30, e2 = 1/6 + 1/10 + 1/15 = 1/3, e3 = 1/30.
    const Matrix<Rational> g = gamma_1list<Rational>();
    CHECK(cache_erec(g, std::vector<int>{0}) == Rational(1));
    CHECK(cache_erec(g, std::vector<int>{1}) == Rational(31, 30));
    CHECK(cache_erec(g, std::vector<int>{2}) == Rational(2, 3));   // 2! * 1/3
    CHECK(cache_erec(g, std::vector<int>{3}) == Rational(1, 5));   // 3! * 1/30
    // More slots than items: no admissible placement.
    CHECK(cache_erec(g, std::vector<int>{4}) == Rational(0));

    // A single item and a single slot is just its own access factor.
    Matrix<Rational> one(1, 1);
    one(0, 0) = Rational(7, 11);
    CHECK(cache_erec(one, std::vector<int>{1}) == Rational(7, 11));
    CHECK(cache_erec(one, std::vector<int>{0}) == Rational(1));
}

TEST_CASE("cache_erec double agrees with the exact rational") {
    const std::vector<std::vector<int>> caps{{0}, {1}, {2}, {3}};
    for (const std::vector<int>& m : caps) {
        const double d = cache_erec(gamma_1list<double>(), m);
        const Rational q = cache_erec(gamma_1list<Rational>(), m);
        CHECK(d == doctest::Approx(static_cast<double>(q)).epsilon(TOL));
        const Real50 r = cache_erec(gamma_1list<Real50>(), m);
        CHECK(static_cast<double>(r) == doctest::Approx(static_cast<double>(q)).epsilon(TOL));
    }

    // Two lists, capacity [1,1].
    const std::vector<int> m11{1, 1};
    CHECK(cache_erec(gamma_2list<double>(), m11) ==
          doctest::Approx(static_cast<double>(cache_erec(gamma_2list<Rational>(), m11)))
              .epsilon(TOL));
}

TEST_CASE("cache_prob_erec single list matches the hand-computed probabilities") {
    // With m = [1], prob(i,2) = gamma(i) / e1 = gamma(i) * 30/31.
    const Matrix<Rational> p = cache_prob_erec(gamma_1list<Rational>(), std::vector<int>{1});
    CHECK(p(0, 1) == Rational(15, 31));
    CHECK(p(1, 1) == Rational(10, 31));
    CHECK(p(2, 1) == Rational(6, 31));
    CHECK(p(0, 0) == Rational(16, 31));
    CHECK(p(1, 0) == Rational(21, 31));
    CHECK(p(2, 0) == Rational(25, 31));
}

TEST_CASE("cache_prob_erec obeys both conservation laws exactly") {
    const Rational one(1);
    SUBCASE("one list") {
        for (int mv = 1; mv <= 3; ++mv) {
            const Matrix<Rational> p = cache_prob_erec(gamma_1list<Rational>(), std::vector<int>{mv});
            Rational occ(0);
            for (std::size_t i = 0; i < p.rows(); ++i) {
                Rational row(0);
                for (std::size_t j = 0; j < p.cols(); ++j) row += p(i, j);
                CHECK(row == one);  // miss plus every hit is exactly one
                occ += p(i, 1);
            }
            CHECK(occ == Rational(mv));  // exactly mv items occupy the list
        }
    }
    SUBCASE("two lists") {
        const std::vector<int> m{1, 1};
        const Matrix<Rational> p = cache_prob_erec(gamma_2list<Rational>(), m);
        Rational occ0(0), occ1(0);
        for (std::size_t i = 0; i < p.rows(); ++i) {
            CHECK(p(i, 0) + p(i, 1) + p(i, 2) == one);
            occ0 += p(i, 1);
            occ1 += p(i, 2);
        }
        CHECK(occ0 == Rational(1));
        CHECK(occ1 == Rational(1));
    }
}

TEST_CASE("cache_prob_erec double agrees with the exact rational") {
    const std::vector<int> m{1, 1};
    const Matrix<double> pd = cache_prob_erec(gamma_2list<double>(), m);
    const Matrix<Rational> pq = cache_prob_erec(gamma_2list<Rational>(), m);
    for (std::size_t i = 0; i < pd.rows(); ++i)
        for (std::size_t j = 0; j < pd.cols(); ++j)
            CHECK(pd(i, j) == doctest::Approx(static_cast<double>(pq(i, j))).epsilon(TOL));
}

TEST_CASE("cache_miss matches the hand-computed rate and closes with cache_prob_erec") {
    // m = [1]: M = E([2])/E([1]) = (2/3)/(31/30) = 20/31.
    Matrix<Rational> lambda(1, 3);
    lambda(0, 0) = Rational(1);
    lambda(0, 1) = Rational(1);
    lambda(0, 2) = Rational(1);
    const auto r = cache_miss(gamma_1list<Rational>(), std::vector<int>{1}, lambda);
    CHECK(r.M == Rational(20, 31));

    // pi0 from the constant ratio must equal the miss column of cache_prob_erec.
    const Matrix<Rational> p = cache_prob_erec(gamma_1list<Rational>(), std::vector<int>{1});
    Rational tot(0);
    for (std::size_t k = 0; k < 3; ++k) {
        CHECK(r.pi0[k] == p(k, 0));
        CHECK(r.MI[k] == r.pi0[k]);  // unit rates
        tot += r.pi0[k];
    }
    // Three items, one slot: exactly two are always missing.
    CHECK(tot == Rational(2));
    CHECK(r.MU[0] == tot);
}

TEST_CASE("cache_miss two lists: pi0 agrees with cache_prob_erec, exactly") {
    const std::vector<int> m{1, 1};
    Matrix<Rational> lambda(2, 3);
    lambda(0, 0) = Rational(1, 2);
    lambda(0, 1) = Rational(1, 3);
    lambda(0, 2) = Rational(1, 7);
    lambda(1, 0) = Rational(2, 5);
    lambda(1, 1) = Rational(1, 11);
    lambda(1, 2) = Rational(3, 4);
    const auto r = cache_miss(gamma_2list<Rational>(), m, lambda);
    const Matrix<Rational> p = cache_prob_erec(gamma_2list<Rational>(), m);
    for (std::size_t k = 0; k < 3; ++k) {
        CHECK(r.pi0[k] == p(k, 0));
        CHECK(r.MI[k] == (lambda(0, k) + lambda(1, k)) * r.pi0[k]);
    }
    Rational mu0(0);
    for (std::size_t k = 0; k < 3; ++k) mu0 += lambda(0, k) * p(k, 0);
    CHECK(r.MU[0] == mu0);

    // double must reproduce the same numbers.
    const auto rd = cache_miss(gamma_2list<double>(), m, line::matrix_from<double>(lambda.view()));
    CHECK(rd.M == doctest::Approx(static_cast<double>(r.M)).epsilon(TOL));
    for (std::size_t k = 0; k < 3; ++k)
        CHECK(rd.pi0[k] == doctest::Approx(static_cast<double>(r.pi0[k])).epsilon(TOL));
}

TEST_CASE("cache_miss without request rates returns only the global rate") {
    const auto r = cache_miss(gamma_1list<double>(), std::vector<int>{1});
    CHECK(r.M == doctest::Approx(20.0 / 31.0).epsilon(TOL));
    CHECK(r.pi0.empty());
    CHECK(r.MU.empty());
}

TEST_CASE("cache_gamma_lp builds the path product of a two-level chain") {
    // One user, one item, two lists in a chain 0 -> 1 -> 2 with rates
    // lambda = [a,b,c] at the three nodes and routing probability one on each
    // edge. gamma(0,0) = a and gamma(0,1) = a*b.
    const Rational a(1, 2), b(1, 3), c(1, 5);
    std::vector<Matrix<Rational>> lambda(1, Matrix<Rational>(1, 3));
    lambda[0](0, 0) = a;
    lambda[0](0, 1) = b;
    lambda[0](0, 2) = c;

    std::vector<std::vector<Matrix<Rational>>> R(1, std::vector<Matrix<Rational>>(1));
    R[0][0] = Matrix<Rational>(3, 3, Rational(0));
    R[0][0](0, 1) = Rational(1);
    R[0][0](1, 2) = Rational(1);

    const auto r = cache_gamma_lp(lambda, R);
    CHECK(r.u == 1u);
    CHECK(r.n == 1u);
    CHECK(r.h == 2u);
    CHECK(r.gamma(0, 0) == a);
    CHECK(r.gamma(0, 1) == a * b);
}

TEST_CASE("cache_gamma_lp aggregates the users and rejects a non-tree structure") {
    // Two users on the same one-list chain: gamma is the summed flow.
    std::vector<Matrix<Rational>> lambda(2, Matrix<Rational>(1, 2));
    lambda[0](0, 0) = Rational(1, 2);
    lambda[0](0, 1) = Rational(0);
    lambda[1](0, 0) = Rational(1, 3);
    lambda[1](0, 1) = Rational(0);

    std::vector<std::vector<Matrix<Rational>>> R(2, std::vector<Matrix<Rational>>(1));
    for (int v = 0; v < 2; ++v) {
        R[v][0] = Matrix<Rational>(2, 2, Rational(0));
        R[v][0](0, 1) = Rational(1);
    }
    CHECK(cache_gamma_lp(lambda, R).gamma(0, 0) == Rational(1, 2) + Rational(1, 3));

    // Node 2 reachable from both 0 and 1: not a tree.
    std::vector<Matrix<Rational>> l3(1, Matrix<Rational>(1, 3, Rational(1)));
    std::vector<std::vector<Matrix<Rational>>> R3(1, std::vector<Matrix<Rational>>(1));
    R3[0][0] = Matrix<Rational>(3, 3, Rational(0));
    R3[0][0](0, 1) = Rational(1);
    R3[0][0](0, 2) = Rational(1);
    R3[0][0](1, 2) = Rational(1);
    CHECK_THROWS_AS(cache_gamma_lp(l3, R3), line::InputError);
}
