/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Oracles: closed-form stationary distributions (two-state chain, M/M/1/K
 * birth-death), the defining equation pi Q = 0 checked exactly, and agreement
 * between the double and exact instantiations.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/mc/ctmc_solve.h"

using line::Matrix;
using line::Rational;
using line::Real50;
using line::mc::ctmc_makeinfgen;
using line::mc::ctmc_solve;

namespace {

constexpr double TOL = 1e-12;

/** Birth-death chain with birth rate lambda and death rate mu, K+1 states. */
template <class T>
Matrix<T> birth_death(int K, long lam_num, long lam_den, long mu_num, long mu_den) {
    Matrix<T> Q(K + 1, K + 1, line::num_traits<T>::from_int(0));
    for (int i = 0; i < K; ++i) {
        Q(i, i + 1) = line::num_traits<T>::from_rational(lam_num, lam_den);
        Q(i + 1, i) = line::num_traits<T>::from_rational(mu_num, mu_den);
    }
    return Q;
}

}  // namespace

TEST_CASE("ctmc_makeinfgen zeroes the row sums and ignores the old diagonal") {
    Matrix<double> Q{{99.0, 2.0}, {3.0, -7.0}};
    Matrix<double> G = ctmc_makeinfgen(Q);
    CHECK(G(0, 0) == doctest::Approx(-2.0));
    CHECK(G(1, 1) == doctest::Approx(-3.0));
    for (std::size_t i = 0; i < 2; ++i) {
        double s = 0.0;
        for (std::size_t j = 0; j < 2; ++j) s += G(i, j);
        CHECK(s == doctest::Approx(0.0).epsilon(TOL));
    }
}

TEST_CASE("ctmc_solve two-state chain matches the closed form") {
    // States 0,1 with rates a (0->1) and b (1->0): pi = (b, a)/(a+b).
    const double a = 2.0, b = 3.0;
    Matrix<double> Q{{0.0, a}, {b, 0.0}};
    std::vector<double> pi = ctmc_solve(Q);
    CHECK(pi[0] == doctest::Approx(b / (a + b)).epsilon(TOL));
    CHECK(pi[1] == doctest::Approx(a / (a + b)).epsilon(TOL));

    // Exact arithmetic returns the rational stationary distribution itself.
    Matrix<Rational> Qq{{Rational(0), Rational(2)}, {Rational(3), Rational(0)}};
    std::vector<Rational> piq = ctmc_solve(Qq);
    CHECK(piq[0] == Rational(3, 5));
    CHECK(piq[1] == Rational(2, 5));
}

TEST_CASE("ctmc_solve M/M/1/K matches the geometric closed form, exactly") {
    // lambda = 1/2, mu = 1, rho = 1/2: pi_k = rho^k (1-rho)/(1-rho^(K+1)).
    const int K = 6;
    Matrix<Rational> Q = birth_death<Rational>(K, 1, 2, 1, 1);
    std::vector<Rational> pi = ctmc_solve(Q);

    Rational rho(1, 2), num(1), denom(0);
    std::vector<Rational> pw(K + 1);
    for (int k = 0; k <= K; ++k) {
        pw[k] = num;
        denom += num;
        num *= rho;
    }
    for (int k = 0; k <= K; ++k) CHECK(pi[k] == pw[k] / denom);
}

TEST_CASE("ctmc_solve residual pi Q = 0 is exactly zero in exact arithmetic") {
    Matrix<Rational> Q(4, 4, Rational(0));
    Q(0, 1) = Rational(1, 3);
    Q(0, 2) = Rational(2, 5);
    Q(1, 0) = Rational(1, 7);
    Q(1, 3) = Rational(3, 4);
    Q(2, 3) = Rational(5, 6);
    Q(2, 0) = Rational(1, 2);
    Q(3, 1) = Rational(2, 3);
    Q(3, 2) = Rational(1, 9);
    Matrix<Rational> G = ctmc_makeinfgen(Q);
    std::vector<Rational> pi = ctmc_solve(G);

    Rational total(0);
    for (const Rational& v : pi) total += v;
    CHECK(total == Rational(1));

    for (std::size_t j = 0; j < 4; ++j) {
        Rational r(0);
        for (std::size_t i = 0; i < 4; ++i) r += pi[i] * G(i, j);
        CHECK(r == Rational(0));  // not "close to zero": identically zero
    }
}

TEST_CASE("ctmc_solve agrees across arithmetics") {
    Matrix<double> Qd(4, 4, 0.0);
    Matrix<Rational> Qq(4, 4, Rational(0));
    Matrix<Real50> Qr(4, 4, Real50(0));
    const long num[4][4] = {{0, 1, 2, 0}, {3, 0, 0, 1}, {0, 4, 0, 5}, {1, 0, 2, 0}};
    for (std::size_t i = 0; i < 4; ++i)
        for (std::size_t j = 0; j < 4; ++j) {
            Qd(i, j) = static_cast<double>(num[i][j]) / 7.0;
            Qq(i, j) = Rational(num[i][j], 7);
            Qr(i, j) = Real50(num[i][j]) / Real50(7);
        }
    std::vector<double> pd = ctmc_solve(Qd);
    std::vector<Rational> pq = ctmc_solve(Qq);
    std::vector<Real50> pr = ctmc_solve(Qr);
    for (std::size_t i = 0; i < 4; ++i) {
        CHECK(static_cast<double>(pq[i]) == doctest::Approx(pd[i]).epsilon(1e-9));
        CHECK(static_cast<double>(pr[i]) == doctest::Approx(pd[i]).epsilon(1e-9));
    }
}

TEST_CASE("ctmc_solve reducible generator solves each component and renormalizes") {
    // Two disconnected two-state chains: each contributes half the mass.
    Matrix<Rational> Q(4, 4, Rational(0));
    Q(0, 1) = Rational(2);
    Q(1, 0) = Rational(3);
    Q(2, 3) = Rational(1);
    Q(3, 2) = Rational(1);
    std::vector<Rational> pi = ctmc_solve(Q);
    Rational total(0);
    for (const Rational& v : pi) total += v;
    CHECK(total == Rational(1));
    // Within a component the ratio is preserved: pi0/pi1 = 3/2, pi2 = pi3.
    CHECK(pi[0] / pi[1] == Rational(3, 2));
    CHECK(pi[2] == pi[3]);
}

TEST_CASE("ctmc_solve answers an absorbing chain with the point mass") {
    // State 1 is ABSORBING and state 0 feeds it, so all the mass ends up in 1 and
    // pi = [0 1] is the unique solution of pi Q = 0. This used to THROW: the trim
    // dropped a state with an all-zero row, which is precisely the absorbing one,
    // and then cascaded through its feeders until nothing was left. Only an
    // ISOLATED state -- nothing in, nothing out -- is trimmed now, and since
    // ctmc_makeinfgen gives any state with an outgoing rate a nonzero diagonal, an
    // isolated state is one with no rate at all: the all-zero generator, which is
    // answered uniformly above. The trim can therefore no longer empty the generator.
    Matrix<double> Q(2, 2, 0.0);
    Q(0, 1) = 1.0;
    const std::vector<double> pi = ctmc_solve(Q);
    REQUIRE(pi.size() == 2u);
    CHECK(std::abs(pi[0]) < 1e-12);
    CHECK(pi[1] == doctest::Approx(1.0));
}
