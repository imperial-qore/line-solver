/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * dtmc_stochcomp (WS-E): the stochastic complement S = P11 + P12 (I-P22)^-1 P21
 * of a DTMC over a retained set. Exact under Rational, so it is checked both in
 * double against a hand computation and in exact rational arithmetic.
 */

#include "doctest.h"
#include "line/api/mc/dtmc_stochcomp.h"
#include "line/num/number.h"

using namespace line;

namespace {

TEST_CASE("retaining every state returns the matrix unchanged") {
    Matrix<double> P(2, 2);
    P(0, 0) = 0.3; P(0, 1) = 0.7;
    P(1, 0) = 0.6; P(1, 1) = 0.4;
    const Matrix<double> S = mc::dtmc_stochcomp(P, {0, 1});
    CHECK(S(0, 0) == doctest::Approx(0.3));
    CHECK(S(1, 1) == doctest::Approx(0.4));
}

TEST_CASE("complement over one state censors the transient middle state") {
    // 3-state chain, retain {0, 2}, censor state 1.
    // P = [[0, 1, 0], [0.5, 0, 0.5], [0, 1, 0]].
    // From 0 you always pass through 1, which splits 0.5/0.5 to {0,2}; likewise
    // from 2. So S = [[0.5,0.5],[0.5,0.5]].
    Matrix<double> P(3, 3, 0.0);
    P(0, 1) = 1.0;
    P(1, 0) = 0.5; P(1, 2) = 0.5;
    P(2, 1) = 1.0;
    const Matrix<double> S = mc::dtmc_stochcomp(P, {0, 2});
    CHECK(S.rows() == 2);
    CHECK(S(0, 0) == doctest::Approx(0.5));
    CHECK(S(0, 1) == doctest::Approx(0.5));
    CHECK(S(1, 0) == doctest::Approx(0.5));
    CHECK(S(1, 1) == doctest::Approx(0.5));
    // rows of a stochastic complement stay stochastic
    CHECK(S(0, 0) + S(0, 1) == doctest::Approx(1.0));
}

TEST_CASE("the complement is exact under Rational") {
    Matrix<Rational> P(3, 3, Rational(0));
    P(0, 1) = Rational(1);
    P(1, 0) = Rational(1, 2); P(1, 2) = Rational(1, 2);
    P(2, 1) = Rational(1);
    const Matrix<Rational> S = mc::dtmc_stochcomp(P, {0, 2});
    CHECK(S(0, 0) == Rational(1, 2));
    CHECK(S(0, 1) == Rational(1, 2));
}

}  // namespace
