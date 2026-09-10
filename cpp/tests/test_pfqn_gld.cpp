/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Tests for the load-dependent normalizing constant. Three independent
 * oracles are used:
 *
 *  1. pfqn_ca. Setting every rate to one turns pfqn_gld into the classical
 *     convolution, so the two must agree as rationals with no rounding at all,
 *     not merely to a tolerance. The same holds when a think time is expressed
 *     the way pfqn_gld requires it, as an infinite-server station with rates
 *     1, 2, ..., Nt: its balance function collapses to Z^n/n!, which is
 *     exactly pfqn_ca's delay term.
 *  2. Closed-form values for one-station models, where G reduces to a single
 *     balance function and can be written down by hand.
 *  3. mp_pfqn's GMP implementation (bin/gld -e), whose exact numerator and
 *     denominator are embedded here for three of its .qn models.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/pfqn/pfqn_ca.h"
#include "line/api/pfqn/pfqn_gld.h"

using line::BigInt;
using line::InputError;
using line::Matrix;
using line::Rational;
using line::Real50;
using line::pfqn::pfqn_ca;
using line::pfqn::pfqn_gld;

namespace {

constexpr double EXACT_TOL = 1e-8;

/** Two queueing stations, two classes, the pfqn_ca test model. */
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

/** (M x Nt) rate matrix of all ones: every station a single server. */
template <class T>
Matrix<T> unit_rates(std::size_t M, std::size_t Nt) {
    return Matrix<T>(M, Nt, line::num_traits<T>::from_int(1));
}

/** Appends the think-time row and its infinite-server rates 1, 2, ..., Nt. */
template <class T>
void append_delay(const Matrix<T>& L, const Matrix<T>& Z, std::size_t Nt, Matrix<T>& Lout,
                  Matrix<T>& muout) {
    const std::size_t M = L.rows(), R = L.cols();
    Lout = Matrix<T>(M + 1, R);
    muout = Matrix<T>(M + 1, Nt);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t r = 0; r < R; ++r) Lout(i, r) = L(i, r);
    for (std::size_t r = 0; r < R; ++r) Lout(M, r) = Z(0, r);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t k = 0; k < Nt; ++k) muout(i, k) = line::num_traits<T>::from_int(1);
    for (std::size_t k = 0; k < Nt; ++k)
        muout(M, k) = line::num_traits<T>::from_int(static_cast<long>(k) + 1);
}

/** Rate row of an S-server station: mu(k) = min(k, S). */
template <class T>
void multiserver_row(Matrix<T>& mu, std::size_t i, long S, std::size_t Nt) {
    for (std::size_t k = 1; k <= Nt; ++k)
        mu(i, k - 1) = line::num_traits<T>::from_int(static_cast<long>(k) < S ? static_cast<long>(k) : S);
}

}  // namespace

TEST_CASE("pfqn_gld with unit rates is pfqn_ca, exactly in rational arithmetic") {
    const std::vector<int> N{2, 1};

    Matrix<Rational> Lq, Zq;
    multiclass_model(Lq, Zq);
    const Rational gca = pfqn_ca(Lq, N, Matrix<Rational>()).G;
    const Rational ggld = pfqn_gld(Lq, N, unit_rates<Rational>(2, 3)).G;
    CHECK(ggld == gca);
    // The convenience overload builds the same all-ones rate matrix.
    CHECK(pfqn_gld(Lq, N).G == gca);

    // Inexact arithmetics run the identical recursion, so they agree to rounding.
    Matrix<double> Ld, Zd;
    multiclass_model(Ld, Zd);
    CHECK(pfqn_gld(Ld, N, unit_rates<double>(2, 3)).lG ==
          doctest::Approx(pfqn_ca(Ld, N, Matrix<double>()).lG).epsilon(EXACT_TOL));

    Matrix<Real50> Lr, Zr;
    multiclass_model(Lr, Zr);
    CHECK(pfqn_gld(Lr, N, unit_rates<Real50>(2, 3)).lG ==
          doctest::Approx(pfqn_ca(Lr, N, Matrix<Real50>()).lG).epsilon(EXACT_TOL));
}

TEST_CASE("pfqn_gld reproduces the pfqn_ca delay term as an infinite server") {
    const std::vector<int> N{2, 1};
    const std::size_t Nt = 3;

    Matrix<Rational> Lq, Zq;
    multiclass_model(Lq, Zq);
    Matrix<Rational> Lz, muz;
    append_delay(Lq, Zq, Nt, Lz, muz);
    CHECK(pfqn_gld(Lz, N, muz).G == pfqn_ca(Lq, N, Zq).G);

    Matrix<Real50> Lr, Zr;
    multiclass_model(Lr, Zr);
    Matrix<Real50> Lzr, muzr;
    append_delay(Lr, Zr, Nt, Lzr, muzr);
    CHECK(pfqn_gld(Lzr, N, muzr).lG == doctest::Approx(pfqn_ca(Lr, N, Zr).lG).epsilon(EXACT_TOL));

    // A pure delay: G must be Z^N / N!, the same value pfqn_ca returns for an
    // empty demand matrix.
    Matrix<Rational> Ld(1, 1);
    Ld(0, 0) = Rational(1, 2);
    Matrix<Rational> mud(1, 3);
    for (std::size_t k = 0; k < 3; ++k) mud(0, k) = Rational(static_cast<long>(k) + 1);
    CHECK(pfqn_gld(Ld, std::vector<int>{3}, mud).G == Rational(1, 48));
    Matrix<Rational> Zd(1, 1);
    Zd(0, 0) = Rational(1, 2);
    CHECK(pfqn_ca(Matrix<Rational>(), std::vector<int>{3}, Zd).G == Rational(1, 48));
}

TEST_CASE("pfqn_gld matches the closed form of a single load-dependent station") {
    // One station, one class, L = 1/2, N = 3. With a single balance function
    // G = L^N / prod_{k=1}^{N} mu(k).
    Matrix<Rational> L(1, 1);
    L(0, 0) = Rational(1, 2);

    // Single server: G = (1/2)^3.
    CHECK(pfqn_gld(L, std::vector<int>{3}, unit_rates<Rational>(1, 3)).G == Rational(1, 8));

    // Two servers, mu = 1, 2, 2: G = (1/8) / 4.
    Matrix<Rational> mu2(1, 3);
    multiserver_row(mu2, 0, 2, 3);
    CHECK(pfqn_gld(L, std::vector<int>{3}, mu2).G == Rational(1, 32));

    // Three servers, mu = 1, 2, 3: G = (1/8) / 6.
    Matrix<Rational> mu3(1, 3);
    multiserver_row(mu3, 0, 3, 3);
    CHECK(pfqn_gld(L, std::vector<int>{3}, mu3).G == Rational(1, 48));
    CHECK(pfqn_gld(L, std::vector<int>{3}, mu3).lG == doctest::Approx(std::log(1.0 / 48.0)));

    // Multiclass single station: G = (n1+n2)!/(n1! n2!) L1^n1 L2^n2 / prod mu.
    Matrix<Rational> Lm(1, 2);
    Lm(0, 0) = Rational(1, 2);
    Lm(0, 1) = Rational(1, 3);
    Matrix<Rational> mum(1, 3);
    multiserver_row(mum, 0, 2, 3);
    // n = (2,1): 3!/(2!1!) = 3, (1/2)^2 (1/3) = 1/12, prod mu = 1*2*2 = 4.
    CHECK(pfqn_gld(Lm, std::vector<int>{2, 1}, mum).G == Rational(3, 48));
}

TEST_CASE("pfqn_gld agrees with mp_pfqn's GMP reference on its .qn models") {
    // models/13_gld_small.qn: no MU section, so mp_pfqn's gld_auto_mu builds
    // mu(k) = min(k, mi) = 1 for both stations, a load-independent model.
    // bin/gld -e reports 14560000 / 1.
    Matrix<Rational> L13(2, 2);
    L13(0, 0) = Rational(10);
    L13(0, 1) = Rational(20);
    L13(1, 0) = Rational(30);
    L13(1, 1) = Rational(40);
    const std::vector<int> N13{2, 2};
    const Rational g13 = pfqn_gld(L13, N13, unit_rates<Rational>(2, 4)).G;
    CHECK(BigInt(numerator(g13)).str() == "14560000");
    CHECK(BigInt(denominator(g13)).str() == "1");

    // models/14_ld_multi.qn: same demands, MU = [1 1 1 1; 2 3 3 3].
    // bin/gld -e reports 3440000 / 3.
    Matrix<Rational> mu14(2, 4);
    for (std::size_t k = 0; k < 4; ++k) mu14(0, k) = Rational(1);
    mu14(1, 0) = Rational(2);
    mu14(1, 1) = Rational(3);
    mu14(1, 2) = Rational(3);
    mu14(1, 3) = Rational(3);
    const Rational g14 = pfqn_gld(L13, N13, mu14).G;
    CHECK(BigInt(numerator(g14)).str() == "3440000");
    CHECK(BigInt(denominator(g14)).str() == "3");

    // models/15_repairman.qn: N = (3,2), Z = (10,20), two queueing stations
    // with MU = [1 2 3 3 3; 1 2 3 4 5], and the think time carried by the
    // appended infinite server. bin/gld -e reports 13541728 / 27.
    Matrix<Rational> L15(3, 2);
    L15(0, 0) = Rational(5);
    L15(0, 1) = Rational(8);
    L15(1, 0) = Rational(3);
    L15(1, 1) = Rational(4);
    L15(2, 0) = Rational(10);  // think time as an infinite-server station
    L15(2, 1) = Rational(20);
    Matrix<Rational> mu15(3, 5);
    const long r0[5] = {1, 2, 3, 3, 3};
    const long r1[5] = {1, 2, 3, 4, 5};
    for (std::size_t k = 0; k < 5; ++k) {
        mu15(0, k) = Rational(r0[k]);
        mu15(1, k) = Rational(r1[k]);
        mu15(2, k) = Rational(static_cast<long>(k) + 1);
    }
    const Rational g15 = pfqn_gld(L15, std::vector<int>{3, 2}, mu15).G;
    CHECK(BigInt(numerator(g15)).str() == "13541728");
    CHECK(BigInt(denominator(g15)).str() == "27");
}

TEST_CASE("pfqn_gld edge cases follow the MATLAB contract") {
    Matrix<double> L(2, 1, 1.0);

    // Zero population: G = 1, lG = 0, whatever the rates.
    auto empty = pfqn_gld(L, std::vector<int>{0});
    CHECK(empty.G == 1.0);
    CHECK(empty.lG == 0.0);

    // Negative population: G = 0, lG = -inf, as in pfqn_ca.
    auto neg = pfqn_gld(L, std::vector<int>{-1});
    CHECK(neg.G == 0.0);
    CHECK(std::isinf(neg.lG));
    CHECK(neg.lG < 0);

    // No station and a positive population: MATLAB returns G = 0.
    auto nostation = pfqn_gld(Matrix<double>(), std::vector<int>{3});
    CHECK(nostation.G == 0.0);
    CHECK(std::isinf(nostation.lG));

    // A station that cannot serve one of the classes still admits the states
    // in which that class is elsewhere.
    Matrix<Rational> Lz(2, 2);
    Lz(0, 0) = Rational(1, 2);
    Lz(0, 1) = Rational(0);
    Lz(1, 0) = Rational(1, 3);
    Lz(1, 1) = Rational(1, 4);
    Matrix<Rational> muz(2, 2);
    multiserver_row(muz, 0, 2, 2);
    multiserver_row(muz, 1, 1, 2);
    // Class 1 must sit at station 2: G = sum over the states, checked against
    // the convolution written out by hand.
    //   n = (1,1): Y1(1,0) Y2(0,1) + Y1(0,0) Y2(1,1)
    //            = (1/2)(1/4) + 1 * (2!/(1!1!)) (1/3)(1/4)
    //            = 1/8 + 1/6 = 7/24
    CHECK(pfqn_gld(Lz, std::vector<int>{1, 1}, muz).G == Rational(7, 24));
}

TEST_CASE("pfqn_gld rejects malformed rate matrices") {
    Matrix<double> L(2, 1, 1.0);

    // Fewer rate columns than jobs.
    Matrix<double> few(2, 2, 1.0);
    CHECK_THROWS_AS(pfqn_gld(L, std::vector<int>{3}, few), InputError);

    // Wrong number of stations.
    Matrix<double> wrongM(3, 3, 1.0);
    CHECK_THROWS_AS(pfqn_gld(L, std::vector<int>{3}, wrongM), InputError);

    // A zero rate leaves the balance function undefined.
    Matrix<double> zeroRate(2, 3, 1.0);
    zeroRate(1, 1) = 0.0;
    CHECK_THROWS_AS(pfqn_gld(L, std::vector<int>{3}, zeroRate), InputError);

    // Demands and population disagreeing on the class count.
    Matrix<double> L2(2, 2, 1.0);
    CHECK_THROWS_AS(pfqn_gld(L2, std::vector<int>{3}, Matrix<double>(2, 3, 1.0)), InputError);
}

TEST_CASE("pfqn_gld scaling keeps lG finite where the unscaled recursion overflows") {
    // Demands far outside the range in which G(N) is representable: the
    // power-of-two rescaling must still return a finite lG, and the exact path,
    // which needs no scaling, must agree with it.
    Matrix<double> L(2, 1);
    L(0, 0) = 1e30;
    L(1, 0) = 1e30;
    Matrix<double> mu(2, 40);
    multiserver_row(mu, 0, 4, 40);
    multiserver_row(mu, 1, 1, 40);
    auto r = pfqn_gld(L, std::vector<int>{40}, mu);
    CHECK(std::isfinite(r.lG));

    const BigInt e30("1000000000000000000000000000000");
    Matrix<Rational> Lq(2, 1);
    Lq(0, 0) = Rational(e30);
    Lq(1, 0) = Lq(0, 0);
    Matrix<Rational> muq(2, 40);
    multiserver_row(muq, 0, 4, 40);
    multiserver_row(muq, 1, 1, 40);
    auto rq = pfqn_gld(Lq, std::vector<int>{40}, muq);
    CHECK(rq.lG == doctest::Approx(r.lG).epsilon(1e-12));
}
