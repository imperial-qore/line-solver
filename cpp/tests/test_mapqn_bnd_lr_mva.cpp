/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * MVA-shaped linear-reduction bounds (mapqn_bnd_lr_mva).
 *
 * Oracles, in decreasing order of strength:
 *  1. MATLAB reference values from mapqn_bnd_lr_mva.m on Test 4 of
 *     matlab/lib/qrf/test_mapqn_bnd.m, R2025a, measured 2026-08-01, for every
 *     one of the six (queue, level) objectives in both senses. Twelve numbers
 *     is enough to pin the whole polytope: a dropped family moves at least one
 *     of them, and the two that sit at zero pin the lower face.
 *  2. Structural properties: min <= max, the bounds lie in [0,1], the mean
 *     lengths of a bound solution carry exactly N, and QN <= N UN.
 *  3. Agreement between the Double and Rational instantiations.
 *
 * There is no collapse oracle here, unlike mapqn_bnd_lr: this model has no
 * product-form special case to fall back onto, since the MAP queue's levels do
 * not disappear at K == 1 (the network is then plain exponential, but the model
 * still carries the B variables and the MCC cuts). The MATLAB values are
 * therefore the primary check, which is why every objective is asserted rather
 * than a representative one.
 */
#include <string>
#include <vector>

#include "doctest.h"
#include "line/api/mapqn/mapqn_bnd_lr_mva.h"
#include "line/api/mapqn/mapqn_params.h"
#include "line/num/number.h"
#include "line/util/matrix.h"

using line::Matrix;
using line::Rational;
using line::mapqn::LrMvaIndex;
using line::mapqn::LrMvaParams;
using line::mapqn::mapqn_bnd_lr_mva;
using line::mapqn::MapqnBndLrMvaResult;
using line::mapqn::MapqnSense;

namespace {

template <class T>
T num(double v) {
    return line::num_traits<T>::from_double(v);
}

/**
 * Test 4 of test_mapqn_bnd.m: three queues, N = 5, the MAP at queue 3 with two
 * levels, queues 1 and 2 exponential at rates 2.0 and 1.5.
 */
template <class T>
LrMvaParams<T> instance_test4() {
    LrMvaParams<T> p;
    p.M = 3;
    p.N = 5;
    p.K = 2;
    p.muM.push_back(num<T>(2.0));
    p.muM.push_back(num<T>(1.5));
    p.muMAP = Matrix<T>(2, 2);
    p.muMAP(0, 0) = num<T>(0.5);
    p.muMAP(0, 1) = num<T>(0.1);
    p.muMAP(1, 0) = num<T>(0.2);
    p.muMAP(1, 1) = num<T>(0.8);
    p.v = Matrix<T>(2, 2, T());
    p.v(0, 1) = num<T>(0.05);
    p.v(1, 0) = num<T>(0.03);
    p.r = Matrix<T>(3, 3, T());
    p.r(0, 1) = num<T>(0.5);
    p.r(0, 2) = num<T>(0.5);
    p.r(1, 0) = line::num_traits<T>::from_int(1);
    p.r(2, 0) = line::num_traits<T>::from_int(1);
    return p;
}

/** A two-queue instance, N = 2, the MAP at queue 2. */
template <class T>
LrMvaParams<T> instance_small() {
    LrMvaParams<T> p;
    p.M = 2;
    p.N = 2;
    p.K = 2;
    p.muM.push_back(line::num_traits<T>::from_int(1));
    p.muMAP = Matrix<T>(2, 2);
    p.muMAP(0, 0) = num<T>(0.6);
    p.muMAP(0, 1) = num<T>(0.2);
    p.muMAP(1, 0) = num<T>(0.3);
    p.muMAP(1, 1) = num<T>(0.5);
    p.v = Matrix<T>(2, 2, T());
    p.v(0, 1) = num<T>(0.1);
    p.v(1, 0) = num<T>(0.05);
    p.r = Matrix<T>(2, 2, T());
    p.r(0, 1) = line::num_traits<T>::from_int(1);
    p.r(1, 0) = line::num_traits<T>::from_int(1);
    return p;
}

}  // namespace

TEST_CASE("Test 4 instance: MATLAB reference values for every objective and sense") {
    // mapqn_bnd_lr_mva.m, R2025a, measured 2026-08-01. Queue and level are
    // 1-based in MATLAB and 0-based here.
    const double ref_max[3][2] = {{0.646953785916, 0.477479257083},
                                  {0.505263157895, 0.454976303318},
                                  {0.605263157895, 0.394736842105}};
    const double ref_min[3][2] = {{0.176807585998, 0.011181497387},
                                  {0.000000000210, 0.000000000197},
                                  {0.483193277487, 0.315126050542}};

    const LrMvaParams<double> p = instance_test4<double>();
    for (int i = 0; i < 3; ++i) {
        for (int k = 0; k < 2; ++k) {
            const MapqnBndLrMvaResult<double> hi = mapqn_bnd_lr_mva(p, i, k, MapqnSense::Max);
            const MapqnBndLrMvaResult<double> lo = mapqn_bnd_lr_mva(p, i, k, MapqnSense::Min);
            REQUIRE(hi.ok);
            REQUIRE(lo.ok);
            CHECK(hi.objective == doctest::Approx(ref_max[i][k]).epsilon(1e-6));
            // The two queue-2 lower bounds are zero; the reference's ~2e-10 is
            // its interior-point residual, so compare absolutely there.
            if (ref_min[i][k] < 1e-8) {
                CHECK(lo.objective < 1e-9);
            } else {
                CHECK(lo.objective == doctest::Approx(ref_min[i][k]).epsilon(1e-6));
            }
        }
    }
}

TEST_CASE("Test 4 instance: the MAP queue's level bounds sum to one") {
    // UN(3,1) + UN(3,2) = 23/38 + 15/38 = 1 in the reference. Both are upper
    // bounds taken separately, so their sum being exactly 1 is a property of
    // this instance's polytope, not an identity, and it is worth pinning.
    const LrMvaParams<Rational> p = instance_test4<Rational>();
    const MapqnBndLrMvaResult<Rational> a = mapqn_bnd_lr_mva(p, 2, 0, MapqnSense::Max);
    const MapqnBndLrMvaResult<Rational> b = mapqn_bnd_lr_mva(p, 2, 1, MapqnSense::Max);
    REQUIRE(a.ok);
    REQUIRE(b.ok);
    CHECK(a.objective + b.objective == line::num_traits<Rational>::from_int(1));
}

TEST_CASE("small instance: MATLAB reference values") {
    // mapqn_bnd_lr_mva.m, M = 2, N = 2, K = 2, measured 2026-08-01.
    const LrMvaParams<double> p = instance_small<double>();
    const MapqnBndLrMvaResult<double> hi = mapqn_bnd_lr_mva(p, 0, 0, MapqnSense::Max);
    const MapqnBndLrMvaResult<double> lo = mapqn_bnd_lr_mva(p, 0, 0, MapqnSense::Min);
    REQUIRE(hi.ok);
    REQUIRE(lo.ok);
    CHECK(hi.objective == doctest::Approx(0.441783100368).epsilon(1e-6));
    CHECK(lo.objective == doctest::Approx(0.214899713470).epsilon(1e-6));
}

TEST_CASE("a bound solution carries the population and respects QMAX") {
    const LrMvaParams<Rational> p = instance_test4<Rational>();
    const MapqnBndLrMvaResult<Rational> r = mapqn_bnd_lr_mva(p, 0, 0, MapqnSense::Max);
    REQUIRE(r.ok);

    Rational total = Rational();
    for (int i = 0; i < p.M; ++i)
        for (int k = 0; k < p.K; ++k)
            total = total + r.QN(static_cast<std::size_t>(i), static_cast<std::size_t>(k));
    CHECK(total == line::num_traits<Rational>::from_int(p.N));  // POPCONSTR

    for (int i = 0; i < p.M; ++i) {  // QMAX
        for (int k = 0; k < p.K; ++k) {
            const std::size_t ii = static_cast<std::size_t>(i), kk = static_cast<std::size_t>(k);
            CHECK(r.QN(ii, kk) <= line::num_traits<Rational>::from_int(p.N) * r.UN(ii, kk));
        }
    }

    for (int i = 0; i < p.M; ++i) {  // UMAX
        Rational busy = Rational();
        for (int k = 0; k < p.K; ++k)
            busy = busy + r.UN(static_cast<std::size_t>(i), static_cast<std::size_t>(k));
        CHECK(busy <= line::num_traits<Rational>::from_int(1));
    }
}

TEST_CASE("bounds are monotone and lie in [0,1]") {
    const LrMvaParams<double> p = instance_test4<double>();
    for (int i = 0; i < 3; ++i) {
        for (int k = 0; k < 2; ++k) {
            const MapqnBndLrMvaResult<double> hi = mapqn_bnd_lr_mva(p, i, k, MapqnSense::Max);
            const MapqnBndLrMvaResult<double> lo = mapqn_bnd_lr_mva(p, i, k, MapqnSense::Min);
            REQUIRE(hi.ok);
            REQUIRE(lo.ok);
            CHECK(lo.objective <= hi.objective + 1e-12);
            CHECK(lo.objective >= -1e-12);
            CHECK(hi.objective <= 1.0 + 1e-12);
        }
    }
}

TEST_CASE("exact and double instantiations agree") {
    const MapqnBndLrMvaResult<Rational> exact =
        mapqn_bnd_lr_mva(instance_test4<Rational>(), 0, 0, MapqnSense::Max);
    const MapqnBndLrMvaResult<double> dbl =
        mapqn_bnd_lr_mva(instance_test4<double>(), 0, 0, MapqnSense::Max);
    REQUIRE(exact.ok);
    REQUIRE(dbl.ok);
    CHECK(line::num_traits<Rational>::to_double(exact.objective) ==
          doctest::Approx(dbl.objective).epsilon(1e-9));
}

TEST_CASE("the variable layout is a bijection onto its own range") {
    const LrMvaIndex x(3, 2);
    std::vector<char> seen(x.num_vars(), 0);
    for (int i = 0; i < 3; ++i) {
        for (int k = 0; k < 2; ++k) {
            REQUIRE(x.UN(i, k) < x.num_vars());
            CHECK(seen[x.UN(i, k)] == 0);
            seen[x.UN(i, k)] = 1;
            REQUIRE(x.QN(i, k) < x.num_vars());
            CHECK(seen[x.QN(i, k)] == 0);
            seen[x.QN(i, k)] = 1;
            for (int j = 0; j < 3; ++j) {
                REQUIRE(x.B(i, k, j) < x.num_vars());
                CHECK(seen[x.B(i, k, j)] == 0);
                seen[x.B(i, k, j)] = 1;
            }
        }
    }
    std::size_t used = 0;
    for (std::size_t i = 0; i < seen.size(); ++i) used += static_cast<std::size_t>(seen[i]);
    CHECK(used == x.num_vars());
    CHECK(x.num_vars() == 2u * 3u * 2u + 3u * 2u * 3u);
}

TEST_CASE("argument validation") {
    const LrMvaParams<double> p = instance_test4<double>();
    CHECK_THROWS(mapqn_bnd_lr_mva(p, -1, 0, MapqnSense::Max));
    CHECK_THROWS(mapqn_bnd_lr_mva(p, 3, 0, MapqnSense::Max));
    CHECK_THROWS(mapqn_bnd_lr_mva(p, 0, 2, MapqnSense::Max));
    LrMvaParams<double> bad = p;
    bad.muM.pop_back();
    CHECK_THROWS(mapqn_bnd_lr_mva(bad, 0, 0, MapqnSense::Max));
    LrMvaParams<double> one = p;
    one.M = 1;
    CHECK_THROWS(mapqn_bnd_lr_mva(one, 0, 0, MapqnSense::Max));
}
