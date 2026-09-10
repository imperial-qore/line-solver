/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Regression for lu_solve under a pivot sequence that moves an
 * already-eliminated row.
 *
 * lu_factor swaps whole rows, multiplier columns included, so LU(i,k) belongs
 * to the row that ends up at position i. A forward substitution that
 * interleaves the swaps with the updates therefore pairs a step-k row identity
 * with a final-order multiplier and returns a wrong solution, while leaving
 * LU = PA intact and every pivot healthy. The matrix below is the transposed
 * balance system of a two-class DPS chain truncated at level two, which is the
 * instance that exposed it; diagonally dominant systems do not, which is why
 * the defect survived the rest of the suite.
 */
#include <vector>

#include "doctest.h"
#include "line/num/number.h"
#include "line/util/linalg.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

using line::Matrix;
using line::Rational;

namespace {
/** Generator of a symmetric two-class DPS chain truncated at total level 2. */
Matrix<Rational> dps_generator() {
    const Rational q[6][6] = {
        {Rational(-1, 2), Rational(1, 4), Rational(0), Rational(1, 4), Rational(0), Rational(0)},
        {Rational(1), Rational(-3, 2), Rational(1, 4), Rational(0), Rational(1, 4), Rational(0)},
        {Rational(0), Rational(1), Rational(-1), Rational(0), Rational(0), Rational(0)},
        {Rational(1), Rational(0), Rational(0), Rational(-3, 2), Rational(1, 4), Rational(1, 4)},
        {Rational(0), Rational(1, 2), Rational(0), Rational(1, 2), Rational(-1), Rational(0)},
        {Rational(0), Rational(0), Rational(0), Rational(1), Rational(0), Rational(-1)}};
    Matrix<Rational> Q(6, 6);
    for (std::size_t i = 0; i < 6; ++i)
        for (std::size_t j = 0; j < 6; ++j) Q(i, j) = q[i][j];
    return Q;
}
}  // namespace

TEST_CASE("lu_solve is correct when the pivot sequence moves eliminated rows") {
    const Matrix<Rational> Q = dps_generator();
    // A x = b with the last balance equation replaced by the normalization.
    Matrix<Rational> A(6, 6);
    for (std::size_t i = 0; i < 6; ++i)
        for (std::size_t j = 0; j < 6; ++j) A(i, j) = (i == 5) ? Rational(1) : Q(j, i);
    std::vector<Rational> b(6, Rational(0));
    b[5] = Rational(1);

    const std::vector<Rational> x = line::solve(A, b);

    // The residual is exactly zero in rational arithmetic; nothing here is
    // approximate, so an inexact solve cannot hide behind a tolerance.
    for (std::size_t i = 0; i < 6; ++i) {
        Rational s(0);
        for (std::size_t j = 0; j < 6; ++j) s += A(i, j) * x[j];
        CHECK(s == b[i]);
    }
    // The pivot sequence really does move rows after they were eliminated.
    Matrix<Rational> LU = A;
    const std::vector<std::size_t> piv = line::lu_factor(LU);
    bool late_move = false;
    for (std::size_t k = 1; k < piv.size(); ++k)
        if (piv[k] != k) late_move = true;
    CHECK(late_move);

    // The stationary distribution of the chain: a symmetric two-class DPS at
    // equal weights is M/M/1 in the total, with a binomial split at each level.
    const Rational expect[6] = {Rational(4, 7),  Rational(1, 7),  Rational(1, 28),
                                Rational(1, 7),  Rational(1, 14), Rational(1, 28)};
    for (std::size_t i = 0; i < 6; ++i) CHECK(x[i] == expect[i]);
}

TEST_CASE("lu_solve and inverse agree on the same pivoted system") {
    const Matrix<Rational> Q = dps_generator();
    Matrix<Rational> A(6, 6);
    for (std::size_t i = 0; i < 6; ++i)
        for (std::size_t j = 0; j < 6; ++j) A(i, j) = (i == 5) ? Rational(1) : Q(j, i);
    std::vector<Rational> b(6, Rational(0));
    b[5] = Rational(1);
    const std::vector<Rational> x = line::solve(A, b);
    const std::vector<Rational> y = line::mulvec(line::inverse(A), b);
    for (std::size_t i = 0; i < 6; ++i) CHECK(x[i] == y[i]);
    // A A^-1 = I, exactly.
    const Matrix<Rational> P = line::matmul(A, line::inverse(A));
    for (std::size_t i = 0; i < 6; ++i)
        for (std::size_t j = 0; j < 6; ++j) CHECK(P(i, j) == Rational(i == j ? 1 : 0));
}
