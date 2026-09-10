/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * General quadratic-reduction bounds (mapqn_bnd_qr), the p1/p2 hybrid model.
 *
 * Oracles, in decreasing order of strength:
 *  1. The symmetric single-phase tandem of Test 1 in
 *     matlab/lib/qrf/test_mapqn_bnd.m. Two exponential queues of equal rate
 *     with N = 3 is product form with four equally likely states, so
 *     U1 = P(n1 >= 1) = 3/4 EXACTLY, and the reference records that the QR
 *     polytope collapses onto it (min == max == 0.75). At line::Rational the
 *     port must return the fraction 3/4 in both senses, not a rounded decimal.
 *     This is also the sharpest available check that no family is missing:
 *     omit THM3 or THM30 and the polytope opens back up to the [0,1] box.
 *  2. Structural properties a valid relaxation must have: min <= max, both in
 *     [0,1], U + IT = 1 per station, sum Q = N, and Q(i,k) <= N U(i,k).
 *  3. Agreement between the Double and Rational instantiations.
 *
 * The single-point property in oracle 1 is what makes this file worth having:
 * the 2026-07-20 defect in mapqn_bnd_qr.m returned [0, 1/3] on an instance
 * whose true range is [0.1970, 0.2182] while every subscript in the file was
 * internally consistent, so only a bound compared against a known value
 * detects an inventory gap.
 */
#include <string>
#include <vector>

#include "doctest.h"
#include "line/api/mapqn/mapqn_bnd_qr.h"
#include "line/api/mapqn/mapqn_params.h"
#include "line/num/number.h"
#include "line/util/matrix.h"

using line::Matrix;
using line::Rational;
using line::mapqn::mapqn_bnd_qr;
using line::mapqn::MapqnBndQrResult;
using line::mapqn::MapqnParams;
using line::mapqn::MapqnSense;
using line::mapqn::QrIndex;

namespace {

template <class T>
T num(double v) {
    return line::num_traits<T>::from_double(v);
}

/**
 * Test 1 of test_mapqn_bnd.m: two single-phase exponential queues of unit rate
 * in a cycle, N = 3. Product form, four equally likely states, U1 = 3/4.
 */
template <class T>
MapqnParams<T> tandem_k1(int N) {
    MapqnParams<T> p;
    p.M = 2;
    p.N = N;
    p.K.assign(2, 1);
    p.mu.push_back(Matrix<T>(1, 1, line::num_traits<T>::from_int(1)));
    p.mu.push_back(Matrix<T>(1, 1, line::num_traits<T>::from_int(1)));
    p.v.push_back(Matrix<T>(1, 1, T()));
    p.v.push_back(Matrix<T>(1, 1, T()));
    p.r = Matrix<T>(2, 2, T());
    p.r(0, 1) = line::num_traits<T>::from_int(1);
    p.r(1, 0) = line::num_traits<T>::from_int(1);
    return p;
}

/**
 * Test 5 of test_mapqn_bnd.m: two queues, two phases each, N = 2. The same
 * network instance_A() of test_mapqn.cpp carries, minus the alpha table, since
 * this model has no load dependence.
 */
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

TEST_CASE("symmetric single-phase tandem: the QR polytope collapses onto U1 = 3/4 exactly") {
    const MapqnParams<Rational> p = tandem_k1<Rational>(3);
    const MapqnBndQrResult<Rational> hi = mapqn_bnd_qr(p, 0, 0, MapqnSense::Max);
    const MapqnBndQrResult<Rational> lo = mapqn_bnd_qr(p, 0, 0, MapqnSense::Min);

    REQUIRE(hi.ok);
    REQUIRE(lo.ok);
    const Rational exact = line::num_traits<Rational>::from_rational(3, 4);
    CHECK(hi.objective == exact);
    CHECK(lo.objective == exact);
}

TEST_CASE("symmetric tandem: the collapsed solution is the product-form distribution") {
    const MapqnParams<Rational> p = tandem_k1<Rational>(3);
    const MapqnBndQrResult<Rational> r = mapqn_bnd_qr(p, 0, 0, MapqnSense::Max);
    REQUIRE(r.ok);

    const Rational three_quarters = line::num_traits<Rational>::from_rational(3, 4);
    const Rational one_quarter = line::num_traits<Rational>::from_rational(1, 4);
    // Both stations are symmetric, so both read the same utilization.
    CHECK(r.U(0, 0) == three_quarters);
    CHECK(r.U(1, 0) == three_quarters);
    CHECK(r.IT(0, 0) == one_quarter);
    CHECK(r.IT(1, 0) == one_quarter);
    // Mean length is 3/2 at each station by symmetry, and the two carry N = 3.
    const Rational three_halves = line::num_traits<Rational>::from_rational(3, 2);
    CHECK(r.Q(0, 0) == three_halves);
    CHECK(r.Q(1, 0) == three_halves);
}

TEST_CASE("symmetric tandem at N = 2: U1 = 2/3 exactly") {
    const MapqnParams<Rational> p = tandem_k1<Rational>(2);
    const MapqnBndQrResult<Rational> hi = mapqn_bnd_qr(p, 0, 0, MapqnSense::Max);
    const MapqnBndQrResult<Rational> lo = mapqn_bnd_qr(p, 0, 0, MapqnSense::Min);
    REQUIRE(hi.ok);
    REQUIRE(lo.ok);
    const Rational exact = line::num_traits<Rational>::from_rational(2, 3);
    CHECK(hi.objective == exact);
    CHECK(lo.objective == exact);
}

TEST_CASE("two-phase instance: MATLAB reference values, both senses") {
    // matlab/lib/qrf/mapqn_bnd_qr.m on the Test 5 instance, R2025a,
    // lpAlgorithm 'interior-point', measured 2026-08-01:
    //   max = 0.218225419307 (exitflag 1), min = 0.196970294459 (exitflag 1).
    // This is the instance whose true range _kb/03-api-layer.md quotes as
    // [0.1970, 0.2182] and on which the 2026-07-20 MATLAB defect returned
    // [0, 1/3].
    const MapqnParams<double> p = twophase<double>();
    const MapqnBndQrResult<double> hi = mapqn_bnd_qr(p, 0, 0, MapqnSense::Max);
    const MapqnBndQrResult<double> lo = mapqn_bnd_qr(p, 0, 0, MapqnSense::Min);
    REQUIRE(hi.ok);
    REQUIRE(lo.ok);

    CHECK(hi.objective == doctest::Approx(0.218225419307).epsilon(1e-6));
    CHECK(lo.objective == doctest::Approx(0.196970294459).epsilon(1e-6));
    // A vacuous bound is the [0,1] box; an inventory gap shows up here first.
    CHECK(hi.objective < 1.0 - 1e-9);
}

TEST_CASE("a bound solution satisfies the aggregate identities it is built on") {
    // Exact arithmetic on the tandem, where the whole polytope is one point.
    // The two-phase instance is checked the same way below in double: one
    // exact solve of its 260-variable model costs ~70 s, which does not belong
    // in a unit test.
    const MapqnParams<Rational> p = tandem_k1<Rational>(3);
    const MapqnBndQrResult<Rational> r = mapqn_bnd_qr(p, 0, 0, MapqnSense::Max);
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

TEST_CASE("the two-phase bound solution satisfies ONE, POPC and QUB1") {
    const MapqnParams<double> p = twophase<double>();
    const MapqnBndQrResult<double> r = mapqn_bnd_qr(p, 0, 0, MapqnSense::Max);
    REQUIRE(r.ok);

    for (int i = 0; i < p.M; ++i) {
        double busy_plus_idle = 0.0;
        for (int k = 0; k < p.K[i]; ++k)
            busy_plus_idle += r.U(static_cast<std::size_t>(i), static_cast<std::size_t>(k)) +
                              r.IT(static_cast<std::size_t>(i), static_cast<std::size_t>(k));
        CHECK(busy_plus_idle == doctest::Approx(1.0).epsilon(1e-9));  // ONE
    }

    double total = 0.0;
    for (int i = 0; i < p.M; ++i)
        for (int k = 0; k < p.K[i]; ++k)
            total += r.Q(static_cast<std::size_t>(i), static_cast<std::size_t>(k));
    CHECK(total == doctest::Approx(static_cast<double>(p.N)).epsilon(1e-9));  // POPC

    for (int i = 0; i < p.M; ++i) {  // QUB1
        for (int k = 0; k < p.K[i]; ++k) {
            const std::size_t ii = static_cast<std::size_t>(i), kk = static_cast<std::size_t>(k);
            CHECK(r.Q(ii, kk) <= static_cast<double>(p.N) * r.U(ii, kk) + 1e-9);
        }
    }
}

TEST_CASE("bounds are monotone in the direction they are taken") {
    const MapqnParams<double> p = twophase<double>();
    for (int k = 0; k < 2; ++k) {
        const MapqnBndQrResult<double> hi = mapqn_bnd_qr(p, 1, k, MapqnSense::Max);
        const MapqnBndQrResult<double> lo = mapqn_bnd_qr(p, 1, k, MapqnSense::Min);
        REQUIRE(hi.ok);
        REQUIRE(lo.ok);
        CHECK(lo.objective <= hi.objective + 1e-12);
    }
}

TEST_CASE("exact and double instantiations agree on the tandem") {
    // On the tandem, not the two-phase instance: agreement between arithmetics
    // is a property of the solver, not of the instance, and the exact solve of
    // the larger model costs ~70 s.
    for (int N = 2; N <= 3; ++N) {
        const MapqnBndQrResult<Rational> exact =
            mapqn_bnd_qr(tandem_k1<Rational>(N), 0, 0, MapqnSense::Max);
        const MapqnBndQrResult<double> dbl =
            mapqn_bnd_qr(tandem_k1<double>(N), 0, 0, MapqnSense::Max);
        REQUIRE(exact.ok);
        REQUIRE(dbl.ok);
        CHECK(line::num_traits<Rational>::to_double(exact.objective) ==
              doctest::Approx(dbl.objective).epsilon(1e-9));
    }
}

TEST_CASE("the variable layout is a bijection onto its own range") {
    std::vector<int> K;
    K.push_back(2);
    K.push_back(1);
    K.push_back(3);
    const QrIndex x(3, 2, K, true);

    std::vector<char> seen(x.num_vars(), 0);
    for (int i = 0; i < 3; ++i) {
        for (int k = 0; k < K[i]; ++k) {
            for (int t = 0; t < 3; ++t) {
                REQUIRE(x.C(i, k, t) < x.num_vars());
                CHECK(seen[x.C(i, k, t)] == 0);
                seen[x.C(i, k, t)] = 1;
            }
            REQUIRE(x.U(i, k) < x.num_vars());
            CHECK(seen[x.U(i, k)] == 0);
            seen[x.U(i, k)] = 1;
            REQUIRE(x.IT(i, k) < x.num_vars());
            CHECK(seen[x.IT(i, k)] == 0);
            seen[x.IT(i, k)] = 1;
            REQUIRE(x.Q(i, k) < x.num_vars());
            CHECK(seen[x.Q(i, k)] == 0);
            seen[x.Q(i, k)] = 1;
            for (int t = 0; t < 3; ++t) {
                for (int nt = 0; nt <= 2; ++nt) {
                    for (int h = 0; h < K[t]; ++h) {
                        REQUIRE(x.p1(i, k, t, nt, h) < x.num_vars());
                        CHECK(seen[x.p1(i, k, t, nt, h)] == 0);
                        seen[x.p1(i, k, t, nt, h)] = 1;
                        REQUIRE(x.p1c(i, k, t, nt, h) < x.num_vars());
                        CHECK(seen[x.p1c(i, k, t, nt, h)] == 0);
                        seen[x.p1c(i, k, t, nt, h)] = 1;
                    }
                }
            }
        }
    }
    for (int j = 0; j < 3; ++j)
        for (int nj = 0; nj <= 2; ++nj)
            for (int k = 0; k < K[j]; ++k)
                for (int i = 0; i < 3; ++i)
                    for (int ni = 0; ni <= 2; ++ni)
                        for (int h = 0; h < K[i]; ++h) {
                            const std::size_t idx = x.p2(j, nj, k, i, ni, h);
                            REQUIRE(idx < x.num_vars());
                            CHECK(seen[idx] == 0);
                            seen[idx] = 1;
                        }

    std::size_t used = 0;
    for (std::size_t i = 0; i < seen.size(); ++i) used += static_cast<std::size_t>(seen[i]);
    CHECK(used == x.num_vars());
}

TEST_CASE("argument validation") {
    const MapqnParams<double> p = twophase<double>();
    CHECK_THROWS(mapqn_bnd_qr(p, -1, 0, MapqnSense::Max));
    CHECK_THROWS(mapqn_bnd_qr(p, 2, 0, MapqnSense::Max));
    CHECK_THROWS(mapqn_bnd_qr(p, 0, 2, MapqnSense::Max));
    MapqnParams<double> zero = p;
    zero.N = 0;
    CHECK_THROWS(mapqn_bnd_qr(zero, 0, 0, MapqnSense::Max));
}
