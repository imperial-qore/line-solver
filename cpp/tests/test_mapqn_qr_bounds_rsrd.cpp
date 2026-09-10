/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Quadratic-reduction bounds under RS-RD blocking (qrf_rsrd).
 *
 * Oracle 1, and it is exact: with every capacity at or above N the blocking
 * never binds, so the polytope collapses onto the product-form solution. For
 * M = 2, N = 2, mu = [1, 0.5], F = [2, 2] the Gordon-Newell constant is
 * G = x1^2 + x1 x2 + x2^2 = 1 + 2 + 4 = 7 with x_i = 1/mu_i, so
 * U1 = 3/7 and U2 = 6/7 EXACTLY, in both senses. qrf_rsrd.m returns
 * 0.428571428571 and 0.857142857143 with min == max, confirming the collapse.
 * Oracle 2: pb must be zero there, since nothing is ever blocked.
 */
#include <vector>

#include "doctest.h"
#include "line/api/mapqn/mapqn_params.h"
#include "line/api/mapqn/mapqn_qr_bounds_rsrd.h"
#include "line/num/number.h"
#include "line/util/matrix.h"

using line::Matrix;
using line::Rational;
using line::mapqn::mapqn_qr_bounds_rsrd;
using line::mapqn::MapqnSense;
using line::mapqn::QrRsrdParams;
using line::mapqn::QrRsrdResult;

namespace {
template <class T>
QrRsrdParams<T> tandem(int N, int F0, int F1, double mu0, double mu1) {
    QrRsrdParams<T> p;
    p.M = 2;
    p.N = N;
    p.F.push_back(F0);
    p.F.push_back(F1);
    p.K.assign(2, 1);
    p.mu.push_back(Matrix<T>(1, 1, line::num_traits<T>::from_double(mu0)));
    p.mu.push_back(Matrix<T>(1, 1, line::num_traits<T>::from_double(mu1)));
    p.v.push_back(Matrix<T>(1, 1, T()));
    p.v.push_back(Matrix<T>(1, 1, T()));
    p.r = Matrix<T>(2, 2, T());
    p.r(0, 1) = line::num_traits<T>::from_int(1);
    p.r(1, 0) = line::num_traits<T>::from_int(1);
    return p;
}
}  // namespace

TEST_CASE("RS-RD with non-binding capacity collapses onto the product form, exactly") {
    const QrRsrdParams<Rational> p = tandem<Rational>(2, 2, 2, 1.0, 0.5);
    const Rational three_sevenths = line::num_traits<Rational>::from_rational(3, 7);
    const Rational six_sevenths = line::num_traits<Rational>::from_rational(6, 7);
    for (int s = 0; s < 2; ++s) {
        const MapqnSense sense = s == 0 ? MapqnSense::Max : MapqnSense::Min;
        const QrRsrdResult<Rational> a = mapqn_qr_bounds_rsrd(p, 0, sense);
        const QrRsrdResult<Rational> b = mapqn_qr_bounds_rsrd(p, 1, sense);
        REQUIRE(a.ok);
        REQUIRE(b.ok);
        CHECK(a.objective == three_sevenths);
        CHECK(b.objective == six_sevenths);
    }
}

TEST_CASE("RS-RD: nothing is blocked when capacity does not bind") {
    const QrRsrdParams<Rational> p = tandem<Rational>(2, 2, 2, 1.0, 0.5);
    const QrRsrdResult<Rational> r = mapqn_qr_bounds_rsrd(p, 0, MapqnSense::Min);
    REQUIRE(r.ok);
    CHECK(r.pb[0] == Rational());
    CHECK(r.pb[1] == Rational());
    CHECK(r.U[0] == r.Ueff[0]);
    CHECK(r.U[1] == r.Ueff[1]);
}

TEST_CASE("RS-RD: MATLAB agreement in double") {
    const QrRsrdParams<double> p = tandem<double>(2, 2, 2, 1.0, 0.5);
    const QrRsrdResult<double> a = mapqn_qr_bounds_rsrd(p, 0, MapqnSense::Min);
    const QrRsrdResult<double> b = mapqn_qr_bounds_rsrd(p, 1, MapqnSense::Min);
    REQUIRE(a.ok);
    REQUIRE(b.ok);
    CHECK(a.objective == doctest::Approx(0.428571428571).epsilon(1e-9));
    CHECK(b.objective == doctest::Approx(0.857142857143).epsilon(1e-9));
}

TEST_CASE("RS-RD: argument validation") {
    const QrRsrdParams<double> p = tandem<double>(2, 2, 2, 1.0, 0.5);
    CHECK_THROWS(mapqn_qr_bounds_rsrd(p, -1, MapqnSense::Min));
    CHECK_THROWS(mapqn_qr_bounds_rsrd(p, 2, MapqnSense::Min));
    QrRsrdParams<double> bad = p;
    bad.F[0] = 99;
    CHECK_THROWS(mapqn_qr_bounds_rsrd(bad, 0, MapqnSense::Min));
}
