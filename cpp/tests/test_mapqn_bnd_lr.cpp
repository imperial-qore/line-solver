/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * General linear-reduction bounds (mapqn_bnd_lr), the p1-only model.
 *
 * Oracles, in decreasing order of strength:
 *  1. THE K == 1 COLLAPSE. With one phase per station the phase structure
 *     disappears and this model must reproduce mapqn_bnd_lr_pf on the same
 *     instance in the same sense. Both are ported here, so at line::Rational
 *     the two must agree EXACTLY, not to a tolerance. This is the oracle the
 *     MATLAB reference itself uses, where it holds to 1e-7.
 *  2. The symmetric tandem is product form, so U1 = 3/4 at N = 3, exactly.
 *  3. MATLAB reference values from mapqn_bnd_lr.m, R2025a, lpAlgorithm
 *     'interior-point', measured 2026-08-01.
 *  4. A cross-model invariant: LR is a strictly weaker relaxation than the
 *     general QR model, since it keeps no joint variables. Its interval must
 *     therefore CONTAIN the QR interval on the same instance. This catches a
 *     family accidentally shared between the two entry points in the wrong
 *     direction, which the shared mapqn_p1_common.h now makes possible.
 */
#include <string>
#include <vector>

#include "doctest.h"
#include "line/api/mapqn/mapqn_bnd_lr.h"
#include "line/api/mapqn/mapqn_bnd_lr_pf.h"
#include "line/api/mapqn/mapqn_bnd_qr.h"
#include "line/api/mapqn/mapqn_params.h"
#include "line/num/number.h"
#include "line/util/matrix.h"

using line::Matrix;
using line::Rational;
using line::mapqn::LrIndex;
using line::mapqn::LrPfParams;
using line::mapqn::LrPfResult;
using line::mapqn::mapqn_bnd_lr;
using line::mapqn::mapqn_bnd_lr_pf;
using line::mapqn::mapqn_bnd_qr;
using line::mapqn::MapqnBndLrResult;
using line::mapqn::MapqnBndQrResult;
using line::mapqn::MapqnParams;
using line::mapqn::MapqnSense;

namespace {

template <class T>
T num(double v) {
    return line::num_traits<T>::from_double(v);
}

/** Two single-phase queues in a cycle with the given rates. */
template <class T>
MapqnParams<T> single_phase(int N, double mu1, double mu2) {
    MapqnParams<T> p;
    p.M = 2;
    p.N = N;
    p.K.assign(2, 1);
    p.mu.push_back(Matrix<T>(1, 1, num<T>(mu1)));
    p.mu.push_back(Matrix<T>(1, 1, num<T>(mu2)));
    p.v.push_back(Matrix<T>(1, 1, T()));
    p.v.push_back(Matrix<T>(1, 1, T()));
    p.r = Matrix<T>(2, 2, T());
    p.r(0, 1) = line::num_traits<T>::from_int(1);
    p.r(1, 0) = line::num_traits<T>::from_int(1);
    return p;
}

/** The same network as an lr_pf product-form instance. */
template <class T>
LrPfParams<T> single_phase_pf(int N, double mu1, double mu2) {
    LrPfParams<T> p;
    p.M = 2;
    p.N = N;
    p.mu.push_back(num<T>(mu1));
    p.mu.push_back(num<T>(mu2));
    p.r = Matrix<T>(2, 2, T());
    p.r(0, 1) = line::num_traits<T>::from_int(1);
    p.r(1, 0) = line::num_traits<T>::from_int(1);
    return p;
}

/** Test 5 of test_mapqn_bnd.m: two queues, two phases each, N = 2. */
template <class T>
MapqnParams<T> twophase() {
    MapqnParams<T> p;
    p.M = 2;
    p.N = 2;
    p.K.assign(2, 2);
    Matrix<T> mu1(2, 2), mu2(2, 2), v1(2, 2, T()), v2(2, 2, T());
    mu1(0, 0) = num<T>(0.8);
    mu1(0, 1) = num<T>(0.2);
    mu1(1, 0) = num<T>(0.1);
    mu1(1, 1) = num<T>(0.6);
    mu2(0, 0) = num<T>(0.5);
    mu2(0, 1) = num<T>(0.1);
    mu2(1, 0) = num<T>(0.2);
    mu2(1, 1) = num<T>(0.7);
    v1(0, 1) = num<T>(0.1);
    v1(1, 0) = num<T>(0.05);
    v2(0, 1) = num<T>(0.05);
    v2(1, 0) = num<T>(0.1);
    p.mu.push_back(mu1);
    p.mu.push_back(mu2);
    p.v.push_back(v1);
    p.v.push_back(v2);
    p.r = Matrix<T>(2, 2, T());
    p.r(0, 1) = line::num_traits<T>::from_int(1);
    p.r(1, 0) = line::num_traits<T>::from_int(1);
    return p;
}

}  // namespace

TEST_CASE("K == 1: the LR bound collapses onto mapqn_bnd_lr_pf, exactly") {
    const int senses[2] = {0, 1};
    for (int s = 0; s < 2; ++s) {
        const MapqnSense sense = senses[s] == 0 ? MapqnSense::Max : MapqnSense::Min;
        for (int N = 2; N <= 3; ++N) {
            const MapqnBndLrResult<Rational> lr =
                mapqn_bnd_lr(single_phase<Rational>(N, 1.0, 1.0), 0, 0, sense);
            const LrPfResult<Rational> pf =
                mapqn_bnd_lr_pf(single_phase_pf<Rational>(N, 1.0, 1.0), 1, sense);
            REQUIRE(lr.ok);
            REQUIRE(pf.ok);
            CHECK(lr.objective == pf.objective);
        }
    }
}

TEST_CASE("K == 1 asymmetric: the LR bound still collapses onto mapqn_bnd_lr_pf") {
    for (int s = 0; s < 2; ++s) {
        const MapqnSense sense = s == 0 ? MapqnSense::Max : MapqnSense::Min;
        const MapqnBndLrResult<Rational> lr =
            mapqn_bnd_lr(single_phase<Rational>(3, 2.0, 1.0), 0, 0, sense);
        const LrPfResult<Rational> pf =
            mapqn_bnd_lr_pf(single_phase_pf<Rational>(3, 2.0, 1.0), 1, sense);
        REQUIRE(lr.ok);
        REQUIRE(pf.ok);
        CHECK(lr.objective == pf.objective);
    }
}

TEST_CASE("symmetric tandem: the LR polytope collapses onto U1 = 3/4 exactly") {
    const MapqnParams<Rational> p = single_phase<Rational>(3, 1.0, 1.0);
    const MapqnBndLrResult<Rational> hi = mapqn_bnd_lr(p, 0, 0, MapqnSense::Max);
    const MapqnBndLrResult<Rational> lo = mapqn_bnd_lr(p, 0, 0, MapqnSense::Min);
    REQUIRE(hi.ok);
    REQUIRE(lo.ok);
    const Rational exact = line::num_traits<Rational>::from_rational(3, 4);
    CHECK(hi.objective == exact);
    CHECK(lo.objective == exact);
}

TEST_CASE("asymmetric single-phase tandem: MATLAB reference values") {
    // mapqn_bnd_lr.m, mu = [2, 1], N = 3, R2025a 'interior-point', 2026-08-01:
    //   max = 0.499999797767, min = 0.454545454567.
    // The reference's interior-point residual is ~2e-7 on the max; the exact
    // optimum this port returns is the vertex it is approaching.
    const MapqnParams<double> p = single_phase<double>(3, 2.0, 1.0);
    const MapqnBndLrResult<double> hi = mapqn_bnd_lr(p, 0, 0, MapqnSense::Max);
    const MapqnBndLrResult<double> lo = mapqn_bnd_lr(p, 0, 0, MapqnSense::Min);
    REQUIRE(hi.ok);
    REQUIRE(lo.ok);
    CHECK(hi.objective == doctest::Approx(0.499999797767).epsilon(1e-5));
    CHECK(lo.objective == doctest::Approx(0.454545454567).epsilon(1e-6));
}

TEST_CASE("two-phase instance: MATLAB reference values, both senses") {
    // mapqn_bnd_lr.m on the Test 5 instance, R2025a 'interior-point',
    // measured 2026-08-01: max = 0.218225418927, min = 0.196969697432.
    const MapqnParams<double> p = twophase<double>();
    const MapqnBndLrResult<double> hi = mapqn_bnd_lr(p, 0, 0, MapqnSense::Max);
    const MapqnBndLrResult<double> lo = mapqn_bnd_lr(p, 0, 0, MapqnSense::Min);
    REQUIRE(hi.ok);
    REQUIRE(lo.ok);
    CHECK(hi.objective == doctest::Approx(0.218225418927).epsilon(1e-6));
    CHECK(lo.objective == doctest::Approx(0.196969697432).epsilon(1e-6));
    // A vacuous bound is the [0,1] box; an inventory gap shows up here first.
    CHECK(hi.objective < 1.0 - 1e-9);
    CHECK(lo.objective > 1e-9);
}

TEST_CASE("LR is the weaker relaxation: its interval contains the QR interval") {
    const MapqnParams<double> p = twophase<double>();
    const MapqnBndLrResult<double> lr_hi = mapqn_bnd_lr(p, 0, 0, MapqnSense::Max);
    const MapqnBndLrResult<double> lr_lo = mapqn_bnd_lr(p, 0, 0, MapqnSense::Min);
    const MapqnBndQrResult<double> qr_hi = mapqn_bnd_qr(p, 0, 0, MapqnSense::Max);
    const MapqnBndQrResult<double> qr_lo = mapqn_bnd_qr(p, 0, 0, MapqnSense::Min);
    REQUIRE(lr_hi.ok);
    REQUIRE(lr_lo.ok);
    REQUIRE(qr_hi.ok);
    REQUIRE(qr_lo.ok);
    CHECK(lr_lo.objective <= qr_lo.objective + 1e-9);
    CHECK(lr_hi.objective >= qr_hi.objective - 1e-9);
}

TEST_CASE("a bound solution satisfies the aggregate identities it is built on") {
    const MapqnParams<Rational> p = twophase<Rational>();
    const MapqnBndLrResult<Rational> r = mapqn_bnd_lr(p, 0, 0, MapqnSense::Max);
    REQUIRE(r.ok);

    const Rational one = line::num_traits<Rational>::from_int(1);
    for (int i = 0; i < p.M; ++i) {
        Rational busy_plus_idle = Rational();
        for (int k = 0; k < p.K[i]; ++k)
            busy_plus_idle = busy_plus_idle +
                             r.U(static_cast<std::size_t>(i), static_cast<std::size_t>(k)) +
                             r.IT(static_cast<std::size_t>(i), static_cast<std::size_t>(k));
        CHECK(busy_plus_idle == one);  // ONE
    }

    Rational total = Rational();
    for (int i = 0; i < p.M; ++i)
        for (int k = 0; k < p.K[i]; ++k)
            total = total + r.Q(static_cast<std::size_t>(i), static_cast<std::size_t>(k));
    CHECK(total == line::num_traits<Rational>::from_int(p.N));  // POPC
}

TEST_CASE("bounds are monotone in the direction they are taken") {
    const MapqnParams<double> p = twophase<double>();
    for (int k = 0; k < 2; ++k) {
        const MapqnBndLrResult<double> hi = mapqn_bnd_lr(p, 1, k, MapqnSense::Max);
        const MapqnBndLrResult<double> lo = mapqn_bnd_lr(p, 1, k, MapqnSense::Min);
        REQUIRE(hi.ok);
        REQUIRE(lo.ok);
        CHECK(lo.objective <= hi.objective + 1e-12);
    }
}

TEST_CASE("exact and double instantiations agree") {
    const MapqnBndLrResult<Rational> exact =
        mapqn_bnd_lr(twophase<Rational>(), 0, 0, MapqnSense::Max);
    const MapqnBndLrResult<double> dbl = mapqn_bnd_lr(twophase<double>(), 0, 0, MapqnSense::Max);
    REQUIRE(exact.ok);
    REQUIRE(dbl.ok);
    CHECK(line::num_traits<Rational>::to_double(exact.objective) ==
          doctest::Approx(dbl.objective).epsilon(1e-9));
}

TEST_CASE("the LR layout allocates no joint variables") {
    std::vector<int> K;
    K.push_back(2);
    K.push_back(1);
    const LrIndex lr(2, 3, K, false);
    const LrIndex qr(2, 3, K, true);
    CHECK(lr.num_vars() < qr.num_vars());
    CHECK(qr.num_vars() - lr.num_vars() == lr.block * lr.block);
}

TEST_CASE("argument validation") {
    const MapqnParams<double> p = twophase<double>();
    CHECK_THROWS(mapqn_bnd_lr(p, -1, 0, MapqnSense::Max));
    CHECK_THROWS(mapqn_bnd_lr(p, 2, 0, MapqnSense::Max));
    CHECK_THROWS(mapqn_bnd_lr(p, 0, 2, MapqnSense::Max));
    MapqnParams<double> zero = p;
    zero.N = 0;
    CHECK_THROWS(mapqn_bnd_lr(zero, 0, 0, MapqnSense::Max));
}
