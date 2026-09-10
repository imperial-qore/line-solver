/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * dtmc_solve and ctmc_stochcomp. The defining properties are checked exactly:
 * pi P = pi for the DTMC, and for the stochastic complement both that S is a
 * generator (zero row sums) and that its stationary distribution is the
 * conditional distribution of the full chain on the selected states, which is
 * the theorem stochastic complementation rests on.
 */
#include <vector>

#include "doctest.h"
#include "line/api/mc/ctmc_solve.h"
#include "line/api/mc/dtmc_solve.h"

using line::Matrix;
using line::Rational;
using line::mc::ctmc_makeinfgen;
using line::mc::ctmc_solve;
using line::mc::ctmc_stochcomp;
using line::mc::dtmc_solve;

TEST_CASE("dtmc_solve two-state chain matches the closed form, exactly") {
    // P = [[1-a, a], [b, 1-b]] has pi = (b, a)/(a+b).
    Matrix<Rational> P(2, 2);
    P(0, 0) = Rational(3, 4);
    P(0, 1) = Rational(1, 4);
    P(1, 0) = Rational(1, 3);
    P(1, 1) = Rational(2, 3);
    std::vector<Rational> pi = dtmc_solve(P);
    // a = 1/4, b = 1/3 -> pi = (1/3, 1/4)/(7/12) = (4/7, 3/7)
    CHECK(pi[0] == Rational(4, 7));
    CHECK(pi[1] == Rational(3, 7));
}

TEST_CASE("dtmc_solve satisfies pi P = pi identically in exact arithmetic") {
    Matrix<Rational> P(3, 3);
    P(0, 0) = Rational(1, 2); P(0, 1) = Rational(1, 4); P(0, 2) = Rational(1, 4);
    P(1, 0) = Rational(1, 3); P(1, 1) = Rational(1, 3); P(1, 2) = Rational(1, 3);
    P(2, 0) = Rational(1, 6); P(2, 1) = Rational(1, 2); P(2, 2) = Rational(1, 3);
    std::vector<Rational> pi = dtmc_solve(P);

    Rational total(0);
    for (const Rational& v : pi) total += v;
    CHECK(total == Rational(1));

    for (std::size_t j = 0; j < 3; ++j) {
        Rational s(0);
        for (std::size_t i = 0; i < 3; ++i) s += pi[i] * P(i, j);
        CHECK(s == pi[j]);
    }
}

TEST_CASE("ctmc_stochcomp produces a generator whose stationary vector is the conditional one") {
    Matrix<Rational> Q(4, 4, Rational(0));
    Q(0, 1) = Rational(1, 3); Q(0, 2) = Rational(2, 5);
    Q(1, 0) = Rational(1, 7); Q(1, 3) = Rational(3, 4);
    Q(2, 3) = Rational(5, 6); Q(2, 0) = Rational(1, 2);
    Q(3, 1) = Rational(2, 3); Q(3, 2) = Rational(1, 9);
    Matrix<Rational> G = ctmc_makeinfgen(Q);

    const std::vector<std::size_t> I{0, 1};
    auto sc = ctmc_stochcomp(G, I);

    // S is a generator: every row sums to zero, exactly.
    for (std::size_t i = 0; i < 2; ++i) {
        Rational s(0);
        for (std::size_t j = 0; j < 2; ++j) s += sc.S(i, j);
        CHECK(s == Rational(0));
    }

    // Its stationary vector equals the full chain's, conditioned on I.
    std::vector<Rational> piFull = ctmc_solve(G);
    std::vector<Rational> piS = ctmc_solve(sc.S);
    const Rational mass = piFull[0] + piFull[1];
    CHECK(piS[0] == piFull[0] / mass);
    CHECK(piS[1] == piFull[1] / mass);
}

TEST_CASE("ctmc_stochcomp default subset is the first half of the states") {
    Matrix<double> Q(4, 4, 0.0);
    Q(0, 1) = 1.0; Q(1, 2) = 1.0; Q(2, 3) = 1.0; Q(3, 0) = 1.0;
    auto sc = ctmc_stochcomp(ctmc_makeinfgen(Q));
    CHECK(sc.S.rows() == 2);
    CHECK(sc.Q22.rows() == 2);
}

TEST_CASE("ctmc_stochcomp on the whole state set returns the generator unchanged") {
    Matrix<Rational> Q(3, 3, Rational(0));
    Q(0, 1) = Rational(1); Q(1, 2) = Rational(2); Q(2, 0) = Rational(3);
    Matrix<Rational> G = ctmc_makeinfgen(Q);
    auto sc = ctmc_stochcomp(G, std::vector<std::size_t>{0, 1, 2});
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) CHECK(sc.S(i, j) == G(i, j));
}
