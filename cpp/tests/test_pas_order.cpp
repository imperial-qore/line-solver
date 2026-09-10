/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Placement order of pass-and-swap networks. Oracles, in order of strength:
 *   1. The defining property: an ordering is feasible iff it is a linear
 *      extension of the precedence relation, so the placeable-next set of a
 *      full population must be exactly the minimal classes, and repeatedly
 *      placing a placeable class must consume every job. Checked exactly.
 *   2. Transitive closure is idempotent and contains the graph.
 *   3. MATLAB, on the swap graph that actually forces a precedence.
 */
#include <algorithm>
#include <functional>
#include <vector>

#include "doctest.h"
#include "line/api/pfqn/pas_placement.h"
#include "line/api/pfqn/pas_swap2order.h"

using line::Matrix;
using line::Rational;
using line::pfqn::PasRateFun;
using line::pfqn::pas_placement;
using line::pfqn::pas_swap2order;

namespace {

/** Head-of-line service: only the first job in the order is served. */
template <class T>
PasRateFun<T> headRate(const std::vector<T>& mu) {
    return [mu](const std::vector<int>& c) -> T {
        if (c.empty()) return line::num_traits<T>::from_int(0);
        return mu[static_cast<std::size_t>(c[0]) - 1];
    };
}

/** Processor-sharing style: every job in the order contributes its rate. */
template <class T>
PasRateFun<T> sumRate(const std::vector<T>& mu) {
    return [mu](const std::vector<int>& c) -> T {
        T s = line::num_traits<T>::from_int(0);
        for (int r : c) s += mu[static_cast<std::size_t>(r) - 1];
        return s;
    };
}

/** The swap graph of the MATLAB reference run: classes 1 and 2 chase each other. */
template <class T>
Matrix<T> swapPair() {
    Matrix<T> G(3, 3, line::num_traits<T>::from_int(0));
    G(0, 1) = line::num_traits<T>::from_int(1);
    G(1, 0) = line::num_traits<T>::from_int(1);
    return G;
}

}  // namespace

TEST_CASE("pas_swap2order recovers the forced precedence MATLAB reports") {
    const std::vector<Rational> mu{Rational(1), Rational(2), Rational(3)};
    const Matrix<Rational> G = swapPair<Rational>();
    const std::vector<PasRateFun<Rational>> rates{headRate(mu), headRate(mu)};
    const Matrix<Rational> H = pas_swap2order<Rational>({G, G}, rates, {1, 1, 1});
    REQUIRE(H.rows() == 3);
    // MATLAB pas_swap2order({G,G}, {svc,svc}, ones(1,3)) -> H(1,2) = 1 only
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) {
            INFO("H(", i + 1, ",", j + 1, ")");
            CHECK(H(i, j) == Rational(i == 0 && j == 1 ? 1 : 0));
        }

    // the placement order is a class-level property, so a rate function that
    // serves every position instead of only the head must give the same DAG
    const std::vector<PasRateFun<Rational>> ps{sumRate(mu), sumRate(mu)};
    const Matrix<Rational> Hps = pas_swap2order<Rational>({G, G}, ps, {1, 1, 1});
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) CHECK(Hps(i, j) == H(i, j));
}

TEST_CASE("pas_swap2order returns no constraint for a pure OI network") {
    const std::vector<Rational> mu{Rational(1), Rational(2), Rational(3)};
    const Matrix<Rational> Z(3, 3, Rational(0));
    const std::vector<PasRateFun<Rational>> rates{headRate(mu), headRate(mu)};
    const Matrix<Rational> H = pas_swap2order<Rational>({Z, Z}, rates, {1, 1, 1});
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) CHECK(H(i, j) == Rational(0));
    // a single graph is replicated to both queues
    const Matrix<Rational> H1 = pas_swap2order<Rational>({Z}, rates, {1, 1, 1});
    CHECK(H1.rows() == 3);
    CHECK_THROWS_AS(pas_swap2order<Rational>({}, rates, {1, 1, 1}), line::InputError);
    CHECK_THROWS_AS(pas_swap2order<Rational>({Z, Z}, {headRate(mu)}, {1, 1, 1}), line::InputError);
}

TEST_CASE("pas_swap2order with a one-sided swap graph leaves the classes incomparable") {
    // With the swap graph on queue 1 only, both orders of classes 1 and 2 stay
    // reachable, so nothing is forced. MATLAB agrees (H all zero).
    const std::vector<Rational> mu{Rational(1), Rational(2), Rational(3)};
    Matrix<Rational> G(3, 3, Rational(0));
    G(0, 1) = Rational(1);
    const std::vector<PasRateFun<Rational>> rates{headRate(mu), headRate(mu)};
    const Matrix<Rational> H =
        pas_swap2order<Rational>({G, Matrix<Rational>(3, 3, Rational(0))}, rates, {1, 1, 1});
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) CHECK(H(i, j) == Rational(0));
}

TEST_CASE("pas_placement closes the relation and is idempotent") {
    Matrix<Rational> H(3, 3, Rational(0));
    H(0, 1) = Rational(1);  // 1 before 2
    H(1, 2) = Rational(1);  // 2 before 3
    const line::pfqn::PasPlacement<Rational> pl = pas_placement(H);
    CHECK(pl.P(0, 1) == Rational(1));
    CHECK(pl.P(1, 2) == Rational(1));
    CHECK(pl.P(0, 2) == Rational(1));  // the closure adds 1 before 3
    CHECK(pl.P(2, 0) == Rational(0));
    // closing an already closed relation changes nothing
    const line::pfqn::PasPlacement<Rational> again = pas_placement(pl.P);
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) CHECK(again.P(i, j) == pl.P(i, j));
}

TEST_CASE("pas_placement placeable is the set of minimal remaining classes") {
    const Matrix<Rational> H = swapPair<Rational>();  // as a precedence: 1 before 2, 2 before 1
    // Use the DAG the reference actually produces instead, 1 before 2 only.
    Matrix<Rational> D(3, 3, Rational(0));
    D(0, 1) = Rational(1);
    const line::pfqn::PasPlacement<Rational> pl = pas_placement(D);
    const std::vector<Rational> full{Rational(1), Rational(1), Rational(1)};
    CHECK(pl.placeable(full) == std::vector<std::size_t>{0, 2});  // MATLAB: [1 3]
    const std::vector<Rational> noFirst{Rational(0), Rational(1), Rational(1)};
    CHECK(pl.placeable(noFirst) == std::vector<std::size_t>{1, 2});  // MATLAB: [2 3]
    CHECK(H(0, 1) == Rational(1));

    // Placing repeatedly must consume every job: a feasible ordering exists.
    std::vector<Rational> x{Rational(2), Rational(1), Rational(3)};
    std::size_t placed = 0;
    while (true) {
        const std::vector<std::size_t> next = pl.placeable(x);
        if (next.empty()) break;
        x[next.front()] -= Rational(1);
        ++placed;
    }
    CHECK(placed == 6);
    for (const Rational& v : x) CHECK(v == Rational(0));
}

TEST_CASE("pas_placement with an empty graph admits every present class") {
    const line::pfqn::PasPlacement<Rational> pl = pas_placement(Matrix<Rational>());
    CHECK(pl.P.rows() == 0);
    const std::vector<Rational> x{Rational(1), Rational(0), Rational(2)};
    CHECK(pl.placeable(x) == std::vector<std::size_t>{0, 2});
    CHECK_THROWS_AS(pas_placement(Matrix<Rational>(2, 3, Rational(0))), line::InputError);
}

TEST_CASE("pas_placement agrees between exact and double") {
    Matrix<double> H(4, 4, 0.0);
    H(0, 1) = 1.0;
    H(1, 2) = 1.0;
    H(2, 3) = 1.0;
    Matrix<Rational> He(4, 4, Rational(0));
    He(0, 1) = Rational(1);
    He(1, 2) = Rational(1);
    He(2, 3) = Rational(1);
    const line::pfqn::PasPlacement<double> pd = pas_placement(H);
    const line::pfqn::PasPlacement<Rational> pe = pas_placement(He);
    for (std::size_t i = 0; i < 4; ++i)
        for (std::size_t j = 0; j < 4; ++j)
            CHECK(pd.P(i, j) == doctest::Approx(static_cast<double>(pe.P(i, j))));
    // a total order: only the first class is ever placeable
    const std::vector<double> x{1.0, 1.0, 1.0, 1.0};
    CHECK(pd.placeable(x) == std::vector<std::size_t>{0});
}
