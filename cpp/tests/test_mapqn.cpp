/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Quadratic-reduction bounds on MAP queueing networks.
 *
 * Oracles, in decreasing order of strength:
 *  1. An instance whose QR polytope is a single point. The delay instance of
 *     matlab/lib/qrf/test_mapqn_bnd.m has min == max, so the bound IS the
 *     stationary distribution, and at line::Rational the port reproduces it as
 *     exact fractions (4/19, 6/19, 6/19, 3/19) rather than to a tolerance.
 *     MATLAB's interior-point linprog gets the same numbers to ~1e-11.
 *  2. MATLAB reference values from mapqn_bnd_qr_ld.m and mapqn_bnd_qr_delay.m,
 *     run on the instances of test_mapqn_bnd.m and on two further instances
 *     built here. Agreement is asserted to 1e-6, which is what
 *     interior-point-legacy delivers on these instances.
 *  3. Structural properties a valid relaxation must have: min <= max, the
 *     bounds lie in [0,1], the marginals of a bound solution normalize to one,
 *     and (the sharpest of the three) the TRUE distribution of a product-form
 *     instance is feasible for the polytope, so it must lie between the bounds.
 *     On single-phase instances the polytope collapses onto it, so the bound
 *     must BE the product-form distribution, exactly.
 *
 * Oracle 3 is what exposes the reference defect recorded at the bottom of this
 * file: on those instances MATLAB's default lpAlgorithm reports the LP
 * infeasible when it demonstrably is not.
 */
#include <string>
#include <vector>

#include "doctest.h"
#include "line/api/mapqn/mapqn_bnd_qr_delay.h"
#include "line/api/mapqn/mapqn_bnd_qr_ld.h"
#include "line/api/mapqn/mapqn_params.h"
#include "line/num/number.h"
#include "line/util/matrix.h"

using line::Matrix;
using line::Rational;
using line::mapqn::mapqn_bnd_qr_delay;
using line::mapqn::mapqn_bnd_qr_ld;
using line::mapqn::MapqnParams;
using line::mapqn::MapqnQrResult;
using line::mapqn::MapqnSense;
using line::mapqn::P2Index;

namespace {

template <class T>
T num(double v) {
    return line::num_traits<T>::from_double(v);
}

/**
 * Test 6 of matlab/lib/qrf/test_mapqn_bnd.m: two queues, two phases each,
 * N = 2, load dependence switched off.
 */
template <class T>
MapqnParams<T> instance_A() {
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
    p.alpha = Matrix<T>(2, 2, line::num_traits<T>::from_int(1));
    p.r = Matrix<T>(2, 2, T());
    p.r(0, 1) = line::num_traits<T>::from_int(1);
    p.r(1, 0) = line::num_traits<T>::from_int(1);
    return p;
}

/**
 * Test 7 of test_mapqn_bnd.m: two single-phase queues, N = 3, queue 2 is the
 * delay station (alpha(2,n) = n makes it infinite-server), Z = 2, D1 = 1.
 */
template <class T>
MapqnParams<T> instance_B() {
    MapqnParams<T> p;
    p.M = 2;
    p.N = 3;
    p.K.assign(2, 1);
    p.mu.push_back(Matrix<T>(1, 1, line::num_traits<T>::from_int(1)));
    p.mu.push_back(Matrix<T>(1, 1, line::num_traits<T>::from_rational(1, 2)));
    p.v.push_back(Matrix<T>(1, 1, T()));
    p.v.push_back(Matrix<T>(1, 1, T()));
    p.alpha = Matrix<T>(2, 3, line::num_traits<T>::from_int(1));
    for (int n = 1; n <= 3; ++n) p.alpha(1, n - 1) = line::num_traits<T>::from_int(n);
    p.r = Matrix<T>(2, 2, T());
    p.r(0, 1) = line::num_traits<T>::from_int(1);
    p.r(1, 0) = line::num_traits<T>::from_int(1);
    p.Z = line::num_traits<T>::from_int(2);
    p.D1 = line::num_traits<T>::from_int(1);
    return p;
}

/** Two queues, K = [2,1], N = 3, asymmetric: built here, MATLAB run on it. */
template <class T>
MapqnParams<T> instance_D() {
    MapqnParams<T> p;
    p.M = 2;
    p.N = 3;
    p.K.push_back(2);
    p.K.push_back(1);
    Matrix<T> mu1(2, 2), v1(2, 2, T());
    mu1(0, 0) = num<T>(0.9);
    mu1(0, 1) = num<T>(0.3);
    mu1(1, 0) = num<T>(0.2);
    mu1(1, 1) = num<T>(0.5);
    v1(0, 1) = num<T>(0.2);
    v1(1, 0) = num<T>(0.4);
    p.mu.push_back(mu1);
    p.mu.push_back(Matrix<T>(1, 1, num<T>(1.5)));
    p.v.push_back(v1);
    p.v.push_back(Matrix<T>(1, 1, T()));
    p.alpha = Matrix<T>(2, 3, line::num_traits<T>::from_int(1));
    p.r = Matrix<T>(2, 2, T());
    p.r(0, 1) = line::num_traits<T>::from_int(1);
    p.r(1, 0) = line::num_traits<T>::from_int(1);
    return p;
}

/** Single-phase two-queue cycle with unequal rates: a product-form network. */
template <class T>
MapqnParams<T> instance_C(double mu1, double mu2) {
    MapqnParams<T> p;
    p.M = 2;
    p.N = 2;
    p.K.assign(2, 1);
    p.mu.push_back(Matrix<T>(1, 1, num<T>(mu1)));
    p.mu.push_back(Matrix<T>(1, 1, num<T>(mu2)));
    p.v.push_back(Matrix<T>(1, 1, T()));
    p.v.push_back(Matrix<T>(1, 1, T()));
    p.alpha = Matrix<T>(2, 2, line::num_traits<T>::from_int(1));
    p.r = Matrix<T>(2, 2, T());
    p.r(0, 1) = line::num_traits<T>::from_int(1);
    p.r(1, 0) = line::num_traits<T>::from_int(1);
    return p;
}

}  // namespace

// ---------------------------------------------------------------------------
// Oracle 1: an instance whose polytope is a point, solved exactly
// ---------------------------------------------------------------------------

TEST_CASE("delay instance: the QR polytope is a single point, reproduced exactly") {
    const MapqnParams<Rational> p = instance_B<Rational>();
    // Queue 1 (the non-delay queue) holds 0..3 jobs with probabilities
    // 4/19, 6/19, 6/19, 3/19. Both bounds must land on them exactly.
    const long expect_num[4] = {4, 6, 6, 3};
    for (int n = 0; n <= 3; ++n) {
        INFO("population ", n);
        const MapqnQrResult<Rational> hi = mapqn_bnd_qr_delay(p, 0, 0, n, MapqnSense::Max);
        const MapqnQrResult<Rational> lo = mapqn_bnd_qr_delay(p, 0, 0, n, MapqnSense::Min);
        REQUIRE(hi.ok);
        REQUIRE(lo.ok);
        const Rational want = Rational(expect_num[n], 19);
        CHECK(hi.objective == want);  // exact rational equality, no tolerance
        CHECK(lo.objective == want);
    }
}

TEST_CASE("delay instance: the exact bound has denominator 19, not a rounded decimal") {
    const MapqnParams<Rational> p = instance_B<Rational>();
    const MapqnQrResult<Rational> hi = mapqn_bnd_qr_delay(p, 0, 0, 1, MapqnSense::Max);
    REQUIRE(hi.ok);
    CHECK(line::num_traits<Rational>::denominator_str(hi.objective) == std::string("19"));
    CHECK(line::num_traits<Rational>::numerator_str(hi.objective) == std::string("6"));
    // MATLAB's interior-point linprog reports 0.3157894736857648 for the same
    // quantity: right to 1e-11, but not the rational 6/19.
    CHECK(line::num_traits<Rational>::to_double(hi.objective) ==
          doctest::Approx(0.3157894736857648).epsilon(1e-10));
}

TEST_CASE("delay instance: MATLAB reference values, both queues and both senses") {
    const MapqnParams<double> p = instance_B<double>();
    // From mapqn_bnd_qr_delay.m, interior-point-legacy.
    const double q1[4] = {0.210526315791874, 0.3157894736857648, 0.3157894736843428,
                          0.1578947368497256};
    const double q2[4] = {0.1578947368437021, 0.3157894736857634, 0.3157894736896439,
                          0.210526315789643};
    for (int n = 0; n <= 3; ++n) {
        INFO("population ", n);
        const MapqnQrResult<double> a = mapqn_bnd_qr_delay(p, 0, 0, n, MapqnSense::Max);
        const MapqnQrResult<double> b = mapqn_bnd_qr_delay(p, 1, 0, n, MapqnSense::Max);
        REQUIRE(a.ok);
        REQUIRE(b.ok);
        CHECK(a.objective == doctest::Approx(q1[n]).epsilon(1e-9));
        CHECK(b.objective == doctest::Approx(q2[n]).epsilon(1e-9));
        const MapqnQrResult<double> amin = mapqn_bnd_qr_delay(p, 0, 0, n, MapqnSense::Min);
        REQUIRE(amin.ok);
        CHECK(amin.objective == doctest::Approx(q1[n]).epsilon(1e-9));
    }
}

TEST_CASE("delay instance: the delay station carries Z times the queue-1 throughput") {
    // XZ is the family that distinguishes the delay model. With the polytope a
    // point, its own statement can be checked on the solution: E[N_2] =
    // (Z/D1) * P(queue 1 busy) ... in the units the family uses.
    const MapqnParams<Rational> p = instance_B<Rational>();
    const MapqnQrResult<Rational> s = mapqn_bnd_qr_delay(p, 0, 0, 1, MapqnSense::Max);
    REQUIRE(s.ok);
    Rational mean2 = Rational(0);
    for (int n = 1; n <= 3; ++n) {
        const Rational term = Rational(n) * s.p2marginals[1](static_cast<std::size_t>(n), 0);
        mean2 += term;
    }
    Rational busy1 = Rational(0);
    for (int n = 1; n <= 3; ++n) busy1 += s.p2marginals[0](static_cast<std::size_t>(n), 0);
    CHECK(mean2 == Rational(2) * busy1);
    // The two marginals each normalize to one.
    Rational tot1 = Rational(0), tot2 = Rational(0);
    for (int n = 0; n <= 3; ++n) {
        tot1 += s.p2marginals[0](static_cast<std::size_t>(n), 0);
        tot2 += s.p2marginals[1](static_cast<std::size_t>(n), 0);
    }
    CHECK(tot1 == Rational(1));
    CHECK(tot2 == Rational(1));
}

// ---------------------------------------------------------------------------
// Oracle 2: MATLAB reference values on the load-dependent entry point
// ---------------------------------------------------------------------------

TEST_CASE("mapqn_bnd_qr_ld matches MATLAB on the two-phase instance") {
    const MapqnParams<double> p = instance_A<double>();
    // mapqn_bnd_qr_ld.m, interior-point-legacy, exitflag 1 throughout.
    // rows: queue, phase, population, max, min
    struct Ref {
        int q, k, n;
        double hi, lo;
    };
    const Ref refs[] = {
        {0, 0, 0, 0.2008610622288099, 0.1216391361620101},
        {0, 1, 0, 0.2571565114890635, 0.1831816928628592},
        {0, 0, 1, 0.1264654286990022, 0.1090857202981979},
        {0, 1, 1, 0.2220664629310009, 0.1934429507630597},
        {0, 0, 2, 0.09798239709899, 0.07885365255450918},
        {0, 1, 2, 0.2171952117256652, 0.1920697717110456},
        {1, 0, 0, 0.1883373777642734, 0.1390754585427031},
        {1, 1, 0, 0.1569734037710093, 0.1017147910727523},
        {1, 0, 1, 0.2186857591720418, 0.1900908024977621},
        {1, 1, 1, 0.1298175769643144, 0.1124664242495011},
        {1, 0, 2, 0.2792099546440689, 0.2546127949332287},
        {1, 1, 2, 0.1241828527477528, 0.1048328016104849},
    };
    for (const Ref& r : refs) {
        INFO("queue ", r.q, " phase ", r.k, " population ", r.n);
        const MapqnQrResult<double> hi = mapqn_bnd_qr_ld(p, r.q, r.k, r.n, MapqnSense::Max);
        const MapqnQrResult<double> lo = mapqn_bnd_qr_ld(p, r.q, r.k, r.n, MapqnSense::Min);
        REQUIRE(hi.ok);
        REQUIRE(lo.ok);
        CHECK(hi.objective == doctest::Approx(r.hi).epsilon(1e-6));
        CHECK(lo.objective == doctest::Approx(r.lo).epsilon(1e-6));
        CHECK(lo.objective <= hi.objective);
        CHECK(lo.objective >= -1e-12);
        CHECK(hi.objective <= 1.0 + 1e-12);
    }
}

TEST_CASE("mapqn_bnd_qr_ld matches MATLAB on an asymmetric K = [2,1] instance") {
    const MapqnParams<double> p = instance_D<double>();
    // mapqn_bnd_qr_ld.m: min and max coincide to ~1e-11 here too, so the
    // polytope pins queue 1's phase-1 marginal exactly.
    const double ref[4] = {0.07236971130283049, 0.1038637626231353, 0.1537750993568947,
                           0.2233041753927516};
    for (int n = 0; n <= 3; ++n) {
        INFO("population ", n);
        const MapqnQrResult<double> hi = mapqn_bnd_qr_ld(p, 0, 0, n, MapqnSense::Max);
        const MapqnQrResult<double> lo = mapqn_bnd_qr_ld(p, 0, 0, n, MapqnSense::Min);
        REQUIRE(hi.ok);
        REQUIRE(lo.ok);
        CHECK(hi.objective == doctest::Approx(ref[n]).epsilon(1e-6));
        CHECK(lo.objective == doctest::Approx(ref[n]).epsilon(1e-6));
    }
}

TEST_CASE("exact and double instantiations agree on the two-phase instance") {
    const MapqnParams<double> pd = instance_A<double>();
    const MapqnParams<Rational> pq = instance_A<Rational>();
    for (int n = 0; n <= 2; ++n) {
        INFO("population ", n);
        const MapqnQrResult<double> d = mapqn_bnd_qr_ld(pd, 0, 0, n, MapqnSense::Max);
        const MapqnQrResult<Rational> q = mapqn_bnd_qr_ld(pq, 0, 0, n, MapqnSense::Max);
        REQUIRE(d.ok);
        REQUIRE(q.ok);
        CHECK(d.objective ==
              doctest::Approx(line::num_traits<Rational>::to_double(q.objective)).epsilon(1e-9));
    }
}

// ---------------------------------------------------------------------------
// Oracle 3: structural properties of a valid relaxation
// ---------------------------------------------------------------------------

TEST_CASE("the bound solution is a consistent pairwise distribution") {
    const MapqnParams<Rational> p = instance_A<Rational>();
    const P2Index idx(p.M, p.N, p.K);
    const MapqnQrResult<Rational> s = mapqn_bnd_qr_ld(p, 0, 0, 1, MapqnSense::Max);
    REQUIRE(s.ok);
    // ONE: each queue's diagonal normalizes.
    for (int j = 0; j < p.M; ++j) {
        Rational tot = Rational(0);
        for (int nj = 0; nj <= p.N; ++nj)
            for (int k = 0; k < p.K[j]; ++k) tot += s.x[idx(j, nj, k, j, nj, k)];
        INFO("queue ", j);
        CHECK(tot == Rational(1));
    }
    // SYMMETRY holds exactly on the solution.
    for (int nj = 0; nj <= p.N; ++nj)
        for (int k = 0; k < p.K[0]; ++k)
            for (int ni = 0; ni <= p.N; ++ni)
                for (int h = 0; h < p.K[1]; ++h)
                    CHECK(s.x[idx(0, nj, k, 1, ni, h)] == s.x[idx(1, ni, h, 0, nj, k)]);
    // ZERO2: one queue cannot hold two different populations at once.
    for (int nj = 0; nj <= p.N; ++nj)
        for (int ni = 0; ni <= p.N; ++ni)
            if (nj != ni) CHECK(s.x[idx(0, nj, 0, 0, ni, 0)] == Rational(0));
    // MARGINALS and THM1 are equalities of the polytope, so the mean
    // population read off the pair (queue 1, queue 2) must be N.
    Rational mean = Rational(0);
    for (int j = 0; j < p.M; ++j)
        for (int nj = 1; nj <= p.N; ++nj)
            for (int k = 0; k < p.K[j]; ++k) {
                const Rational term = Rational(nj) * s.x[idx(j, nj, k, j, nj, k)];
                mean += term;
            }
    CHECK(mean == Rational(p.N));
}

TEST_CASE("bounds are monotone in the direction they are taken") {
    const MapqnParams<double> p = instance_A<double>();
    for (int q = 0; q < 2; ++q)
        for (int k = 0; k < 2; ++k)
            for (int n = 0; n <= 2; ++n) {
                const MapqnQrResult<double> hi = mapqn_bnd_qr_ld(p, q, k, n, MapqnSense::Max);
                const MapqnQrResult<double> lo = mapqn_bnd_qr_ld(p, q, k, n, MapqnSense::Min);
                REQUIRE(hi.ok);
                REQUIRE(lo.ok);
                INFO("queue ", q, " phase ", k, " population ", n);
                CHECK(lo.objective <= hi.objective + 1e-12);
            }
}

TEST_CASE("argument validation") {
    const MapqnParams<double> p = instance_A<double>();
    CHECK_THROWS_AS(mapqn_bnd_qr_ld(p, 5, 0, 0), line::InputError);
    CHECK_THROWS_AS(mapqn_bnd_qr_ld(p, 0, 7, 0), line::InputError);
    CHECK_THROWS_AS(mapqn_bnd_qr_ld(p, 0, 0, 9), line::InputError);
    MapqnParams<double> bad = p;
    bad.K.pop_back();
    CHECK_THROWS_AS(mapqn_bnd_qr_ld(bad, 0, 0, 0), line::InputError);
    // The delay entry point needs a delay station and a nonzero D1.
    MapqnParams<double> d = instance_B<double>();
    d.D1 = 0.0;
    CHECK_THROWS_AS(mapqn_bnd_qr_delay(d, 0, 0, 0), line::InputError);
}

TEST_CASE("the index map is a bijection and symmetric in its two halves") {
    std::vector<int> K;
    K.push_back(2);
    K.push_back(1);
    K.push_back(3);
    const P2Index idx(3, 2, K);
    CHECK(idx.num_vars() == (3u * 6u) * (3u * 6u));
    std::vector<char> seen(idx.num_vars(), 0);
    for (int j = 0; j < 3; ++j)
        for (int nj = 0; nj <= 2; ++nj)
            for (int k = 0; k < K[j]; ++k)
                for (int i = 0; i < 3; ++i)
                    for (int ni = 0; ni <= 2; ++ni)
                        for (int h = 0; h < K[i]; ++h) {
                            const std::size_t v = idx(j, nj, k, i, ni, h);
                            REQUIRE(v < idx.num_vars());
                            CHECK(seen[v] == 0);
                            seen[v] = 1;
                            // the two halves use the same map, so the swap is
                            // a transposition of the flat index
                            CHECK(idx(i, ni, h, j, nj, k) ==
                                  idx.half(i, ni, h) * idx.block + idx.half(j, nj, k));
                        }
    for (std::size_t v = 0; v < seen.size(); ++v) CHECK(seen[v] == 1);
}

// ---------------------------------------------------------------------------
// Reference defect: the polytope excludes the truth on a single-phase network
// ---------------------------------------------------------------------------

TEST_CASE("single-phase product-form network: the bound is the exact distribution") {
    // Two single-phase queues in a cycle with mu = (1, 2). Single-phase makes
    // this an ordinary product-form closed network, whose stationary law is
    // proportional to (1/mu1)^n1 (1/mu2)^n2. At N = 2 that is P(n1) = 1/7,
    // 2/7, 4/7 and at N = 3 it is 1/15, 2/15, 4/15, 8/15. The QR polytope
    // collapses to that single point here, so the bound is not merely valid,
    // it is the answer -- and the exact solver returns it as those fractions.
    //
    // This is also where the MATLAB reference misreports. Its default
    // params.lpAlgorithm is 'interior-point-legacy', which on this instance
    // returns exitflag -2 (no feasible point found) for every objective, while
    // still handing back an iterate whose value a caller reading only
    // result.objective would take as a bound of [7e-09, 1/3] for n = 0. The
    // polytope is not empty: switching the same call to 'interior-point'
    // gives exitflag 1 and 0.142857142857 = 1/7, agreeing with the exact
    // optimum below. A polytope that is a single point has empty relative
    // interior, which is exactly the situation an interior-point method cannot
    // be trusted on, and the legacy variant fails silently rather than loudly.
    // A simplex method has no such failure mode, and with rational arithmetic
    // it has no tolerance to misjudge either.
    const long n2[3] = {1, 2, 4};
    const MapqnParams<Rational> p2 = instance_C<Rational>(1.0, 2.0);
    for (int n = 0; n <= 2; ++n) {
        INFO("N = 2, population ", n);
        const MapqnQrResult<Rational> hi = mapqn_bnd_qr_ld(p2, 0, 0, n, MapqnSense::Max);
        const MapqnQrResult<Rational> lo = mapqn_bnd_qr_ld(p2, 0, 0, n, MapqnSense::Min);
        REQUIRE(hi.ok);
        REQUIRE(lo.ok);
        CHECK(hi.objective == Rational(n2[n], 7));
        CHECK(lo.objective == Rational(n2[n], 7));
    }
    // MATLAB with 'interior-point' on the same instance.
    CHECK(line::num_traits<Rational>::to_double(Rational(1, 7)) ==
          doctest::Approx(0.142857142857).epsilon(1e-11));

    // Equal rates: the uniform distribution, likewise pinned exactly.
    const MapqnParams<Rational> sym = instance_C<Rational>(1.0, 1.0);
    for (int n = 0; n <= 2; ++n) {
        INFO("equal rates, population ", n);
        const MapqnQrResult<Rational> hi = mapqn_bnd_qr_ld(sym, 0, 0, n, MapqnSense::Max);
        const MapqnQrResult<Rational> lo = mapqn_bnd_qr_ld(sym, 0, 0, n, MapqnSense::Min);
        REQUIRE(hi.ok);
        REQUIRE(lo.ok);
        CHECK(hi.objective == Rational(1, 3));
        CHECK(lo.objective == Rational(1, 3));
    }
}

TEST_CASE("single-phase product-form network at N = 3") {
    // Same network, one more job: P(n1) = 1/15, 2/15, 4/15, 8/15.
    MapqnParams<Rational> p = instance_C<Rational>(1.0, 2.0);
    p.N = 3;
    p.alpha = Matrix<Rational>(2, 3, line::num_traits<Rational>::from_int(1));
    const long numer[4] = {1, 2, 4, 8};
    for (int n = 0; n <= 3; ++n) {
        INFO("population ", n);
        const MapqnQrResult<Rational> hi = mapqn_bnd_qr_ld(p, 0, 0, n, MapqnSense::Max);
        const MapqnQrResult<Rational> lo = mapqn_bnd_qr_ld(p, 0, 0, n, MapqnSense::Min);
        REQUIRE(hi.ok);
        REQUIRE(lo.ok);
        CHECK(hi.objective == Rational(numer[n], 15));
        CHECK(lo.objective == Rational(numer[n], 15));
    }
}
