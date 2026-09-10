/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Port of jar/src/test/java/jline/api/pfqn/nc/PfqnNcApiTest.java, extended to
 * the exact and high-precision instantiations. The oracles are the same: an
 * inline convolution for the single-class case and a direct enumeration over
 * the state space for the multiclass case, both independent of the algorithm
 * under test.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/pfqn/pfqn_ca.h"

using line::Matrix;
using line::Rational;
using line::Real50;
using line::pfqn::pfqn_ca;

namespace {

constexpr double EXACT_TOL = 1e-8;

// Test network: two queueing stations, one delay; single class.
constexpr double L1 = 0.5, L2 = 0.4, Z0 = 0.3;

/** G(N) = sum_{k0+k1+k2=N} Z^k0/k0! L1^k1 L2^k2, exactly as in the Java test. */
double exact_log_g(double l1, double l2, double z, int n) {
    double g = 0.0;
    for (int k0 = 0; k0 <= n; ++k0) {
        double zTerm = std::pow(z, k0);
        for (int f = 2; f <= k0; ++f) zTerm /= f;
        for (int k1 = 0; k1 + k0 <= n; ++k1) {
            int k2 = n - k0 - k1;
            g += zTerm * std::pow(l1, k1) * std::pow(l2, k2);
        }
    }
    return std::log(g);
}

const double LM[2][2] = {{0.5, 0.3}, {0.4, 0.6}};
const double ZM[2] = {0.3, 0.2};
const int NM[2] = {2, 1};

double binom(int n, int k) {
    double b = 1.0;
    for (int i = 1; i <= k; ++i) b = b * (n - k + i) / i;
    return b;
}

double exact_multiclass_log_g() {
    double g = 0.0;
    for (int d1 = 0; d1 <= NM[0]; ++d1) {
        for (int d2 = 0; d2 <= NM[1]; ++d2) {
            double zTerm = std::pow(ZM[0], d1) * std::pow(ZM[1], d2);
            for (int f = 2; f <= d1; ++f) zTerm /= f;
            for (int f = 2; f <= d2; ++f) zTerm /= f;
            int r1 = NM[0] - d1, r2 = NM[1] - d2;
            for (int a1 = 0; a1 <= r1; ++a1) {
                for (int a2 = 0; a2 <= r2; ++a2) {
                    int b1 = r1 - a1, b2 = r2 - a2;
                    double st1 = binom(a1 + a2, a1) * std::pow(LM[0][0], a1) * std::pow(LM[0][1], a2);
                    double st2 = binom(b1 + b2, b1) * std::pow(LM[1][0], b1) * std::pow(LM[1][1], b2);
                    g += zTerm * st1 * st2;
                }
            }
        }
    }
    return std::log(g);
}

/** Build the single-class model in the arithmetic under test. */
template <class T>
void single_class_model(Matrix<T>& L, Matrix<T>& Z) {
    L = Matrix<T>(2, 1);
    L(0, 0) = line::num_traits<T>::from_rational(1, 2);  // 0.5
    L(1, 0) = line::num_traits<T>::from_rational(2, 5);  // 0.4
    Z = Matrix<T>(1, 1);
    Z(0, 0) = line::num_traits<T>::from_rational(3, 10);  // 0.3
}

template <class T>
void multiclass_model(Matrix<T>& L, Matrix<T>& Z) {
    L = Matrix<T>(2, 2);
    L(0, 0) = line::num_traits<T>::from_rational(1, 2);
    L(0, 1) = line::num_traits<T>::from_rational(3, 10);
    L(1, 0) = line::num_traits<T>::from_rational(2, 5);
    L(1, 1) = line::num_traits<T>::from_rational(3, 5);
    Z = Matrix<T>(1, 2);
    Z(0, 0) = line::num_traits<T>::from_rational(3, 10);
    Z(0, 1) = line::num_traits<T>::from_rational(1, 5);
}

}  // namespace

TEST_CASE("pfqn_ca single class matches inline convolution, all arithmetics") {
    const std::vector<int> N{4};
    const double expected = exact_log_g(L1, L2, Z0, 4);

    Matrix<double> Ld, Zd;
    single_class_model(Ld, Zd);
    CHECK(pfqn_ca(Ld, N, Zd).lG == doctest::Approx(expected).epsilon(EXACT_TOL));

    Matrix<Rational> Lq, Zq;
    single_class_model(Lq, Zq);
    CHECK(pfqn_ca(Lq, N, Zq).lG == doctest::Approx(expected).epsilon(EXACT_TOL));

    Matrix<Real50> Lr, Zr;
    single_class_model(Lr, Zr);
    CHECK(pfqn_ca(Lr, N, Zr).lG == doctest::Approx(expected).epsilon(EXACT_TOL));
}

TEST_CASE("pfqn_ca multiclass matches state-space enumeration, all arithmetics") {
    const std::vector<int> N{NM[0], NM[1]};
    const double expected = exact_multiclass_log_g();

    Matrix<double> Ld, Zd;
    multiclass_model(Ld, Zd);
    CHECK(pfqn_ca(Ld, N, Zd).lG == doctest::Approx(expected).epsilon(EXACT_TOL));

    Matrix<Rational> Lq, Zq;
    multiclass_model(Lq, Zq);
    CHECK(pfqn_ca(Lq, N, Zq).lG == doctest::Approx(expected).epsilon(EXACT_TOL));

    Matrix<Real50> Lr, Zr;
    multiclass_model(Lr, Zr);
    CHECK(pfqn_ca(Lr, N, Zr).lG == doctest::Approx(expected).epsilon(EXACT_TOL));
}

TEST_CASE("pfqn_ca exact result is a rational, reported without rounding") {
    // Integer demands make G a rational with a small denominator, so the exact
    // path can be checked against a hand-computed value: with L = [[1]] and
    // N = 2, Z = 0, G = 1 and the delay term never contributes.
    Matrix<Rational> L(1, 1);
    L(0, 0) = Rational(1);
    auto r = pfqn_ca(L, std::vector<int>{2}, Matrix<Rational>());
    CHECK(r.G == Rational(1));
    CHECK(r.lG == doctest::Approx(0.0));

    // Single station, single class, L = 1/2, N = 3, no delay: G = (1/2)^3.
    Matrix<Rational> L2m(1, 1);
    L2m(0, 0) = Rational(1, 2);
    auto r2 = pfqn_ca(L2m, std::vector<int>{3}, Matrix<Rational>());
    CHECK(r2.G == Rational(1, 8));
    CHECK(r2.lG == doctest::Approx(std::log(0.125)));
}

TEST_CASE("pfqn_ca edge cases follow the MATLAB contract") {
    // Zero population: G = 1, lG = 0.
    Matrix<double> L(2, 1, 1.0);
    CHECK(pfqn_ca(L, std::vector<int>{0}).G == 1.0);
    CHECK(pfqn_ca(L, std::vector<int>{0}).lG == 0.0);

    // Negative population: G = 0, lG = -inf.
    auto neg = pfqn_ca(L, std::vector<int>{-1});
    CHECK(neg.G == 0.0);
    CHECK(std::isinf(neg.lG));
    CHECK(neg.lG < 0);

    // Delay-only network (M = 0): G = Z^N / N!.
    Matrix<double> noL;
    Matrix<double> Z(1, 1);
    Z(0, 0) = 2.0;
    auto dly = pfqn_ca(noL, std::vector<int>{3}, Z);
    CHECK(dly.G == doctest::Approx(8.0 / 6.0));
}

TEST_CASE("pfqn_ca scaling keeps lG finite where the unscaled recursion overflows") {
    // Demands large enough that G(N) leaves the double range: the Lam scaling
    // must still return a finite lG (the JAR and MATLAB contract).
    Matrix<double> L(2, 1);
    L(0, 0) = 1e30;
    L(1, 0) = 1e30;
    auto r = pfqn_ca(L, std::vector<int>{40});
    CHECK(std::isfinite(r.lG));
    CHECK(r.lG > 2700.0);  // ~ 40*log(1e30) + log(41)

    // The exact path needs no scaling at all and must agree.
    const line::BigInt e30("1000000000000000000000000000000");
    Matrix<Rational> Lq(2, 1);
    Lq(0, 0) = Rational(e30);
    Lq(1, 0) = Lq(0, 0);
    auto rq = pfqn_ca(Lq, std::vector<int>{40});
    CHECK(rq.lG == doctest::Approx(r.lG).epsilon(1e-12));
}

TEST_CASE("pfqn_ca scaling survives a class with no think time beside one with") {
    // An LQN layer reaches SolverNC as a think-time class beside a zero-Z call
    // class. The per-configuration estimate this scaling replaced asked ONE
    // station -- or the delay -- to hold every class at once, so it discarded
    // the delay entirely here and read kscale off the all-at-the-queue term:
    // -2051.6 against a true log G of -359.134. The scaling then went the WRONG
    // WAY, Z/2^-30 = 1.07e9, and the delay column Z^n/n! overflowed at n=[40,0].
    Matrix<double> L(1, 2);
    L(0, 0) = 1e-9;
    L(0, 1) = 1.0;
    Matrix<double> Z(1, 2);
    Z(0, 0) = 1.0;
    Z(0, 1) = 0.0;
    auto r = pfqn_ca(L, std::vector<int>{99, 1}, Z);
    // 120-digit mpmath convolution of the same recursion
    CHECK(std::isfinite(r.lG));
    CHECK(r.lG == doctest::Approx(-359.13420517157538927).epsilon(1e-12));
    CHECK(r.G > 0.0);
    CHECK(std::isfinite(r.G));
}
