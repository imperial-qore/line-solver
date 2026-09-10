/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Station consolidation. The oracle is pfqn_ca: consolidating identical
 * stations and solving the consolidated model with multiplicities must return
 * the IDENTICAL rational normalizing constant as solving the expanded one, and
 * expansion must be the exact inverse of the merge.
 */
#include <vector>

#include "doctest.h"
#include "line/api/pfqn/pfqn_mva.h"
#include "line/api/pfqn/pfqn_expand.h"
#include "line/api/pfqn/pfqn_mushift.h"
#include "line/api/pfqn/pfqn_unique.h"

using line::Matrix;
using line::Rational;
using namespace line::pfqn;

TEST_CASE("pfqn_unique merges exactly equal rows and counts multiplicities") {
    Matrix<Rational> L(4, 2);
    L(0, 0) = Rational(1, 2); L(0, 1) = Rational(1, 3);
    L(1, 0) = Rational(1, 2); L(1, 1) = Rational(1, 3);  // duplicate of row 0
    L(2, 0) = Rational(1, 5); L(2, 1) = Rational(1, 7);
    L(3, 0) = Rational(1, 2); L(3, 1) = Rational(1, 3);  // duplicate of row 0

    auto u = pfqn_unique(L);
    CHECK(u.L.rows() == 2);
    CHECK(u.mi[0] == 3);
    CHECK(u.mi[1] == 1);
    CHECK(u.mapping[0] == 0);
    CHECK(u.mapping[1] == 0);
    CHECK(u.mapping[2] == 1);
    CHECK(u.mapping[3] == 0);
}

TEST_CASE("pfqn_unique does NOT merge rows that only nearly agree") {
    // The MATLAB reference merges within 1e-14; the port requires exact
    // equality, which is the only meaningful test in the rational field.
    Matrix<Rational> L(2, 1);
    L(0, 0) = Rational(1, 3);
    L(1, 0) = Rational(1, 3) + Rational(1, 1000000000000000000LL);
    auto u = pfqn_unique(L);
    CHECK(u.L.rows() == 2);
}

TEST_CASE("consolidation preserves the exact normalizing constant") {
    // Three stations, two of them identical: the consolidated model with
    // multiplicity 2 must give the same G as the expanded one under MVA.
    Matrix<Rational> L(3, 1);
    L(0, 0) = Rational(1, 2);
    L(1, 0) = Rational(1, 2);
    L(2, 0) = Rational(1, 5);
    const std::vector<int> N{4};

    auto full = pfqn_mva(L, N);
    auto u = pfqn_unique(L);
    auto cons = pfqn_mva(u.L, N, Matrix<Rational>(), u.mi);

    // The consolidated MVA reproduces the throughput of the expanded model.
    CHECK(cons.XN[0] == full.XN[0]);
    CHECK(cons.G == full.G);
}

TEST_CASE("pfqn_expand is the exact inverse of the merge") {
    Matrix<Rational> L(3, 1);
    L(0, 0) = Rational(1, 2);
    L(1, 0) = Rational(1, 2);
    L(2, 0) = Rational(1, 5);
    auto u = pfqn_unique(L);
    auto cons = pfqn_mva(u.L, std::vector<int>{4}, Matrix<Rational>(), u.mi);

    auto e = pfqn_expand(cons.QN, cons.UN, cons.CN, u.mapping);
    CHECK(e.QN.rows() == 3);
    // The two merged stations receive identical rows, by construction.
    CHECK(e.QN(0, 0) == e.QN(1, 0));
    CHECK(e.UN(0, 0) == e.UN(1, 0));
    CHECK(e.CN(0, 0) == e.CN(1, 0));
    CHECK(e.QN(2, 0) == cons.QN(1, 0));
}

TEST_CASE("pfqn_mushift shifts only the selected rows") {
    Matrix<Rational> mu(2, 4);
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t k = 0; k < 4; ++k)
            mu(i, k) = Rational(static_cast<long>(10 * (i + 1) + k), 1);

    auto s = pfqn_mushift(mu, std::vector<std::size_t>{0});
    CHECK(s.cols() == 3);
    CHECK(s(0, 0) == mu(0, 1));  // shifted
    CHECK(s(0, 2) == mu(0, 3));
    CHECK(s(1, 0) == mu(1, 0));  // untouched
    CHECK(s(1, 2) == mu(1, 2));
}
