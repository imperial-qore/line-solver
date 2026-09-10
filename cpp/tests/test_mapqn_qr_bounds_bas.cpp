/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Quadratic-reduction bounds under blocking-after-service (qrf_bas).
 *
 * Oracle: qrf_bas.m in R2025a, 2026-08-01, on a tiny instance (three
 * single-phase queues, N = 2, finite queue 0 with capacity 1, three blocking
 * configurations), for every queue in both senses. The instance is small on
 * purpose: this model is MR * B^2 columns and lp::simplex_solve uses a dense
 * tableau, so example_bas_small.m at 3894 columns does not finish.
 *
 * The MATLAB oracle is WEAKER here than for the three unblocked bounds: on the
 * paper instance its residual is conditioning-bound at 8.5e-07 and linprog can
 * return finite-but-infeasible points, which is why qrf_lp_residual exists.
 * The tolerance below is therefore 1e-5, not the 1e-6 used elsewhere in this
 * family, and the structural checks carry proportionally more weight.
 */
#include <string>
#include <vector>

#include "doctest.h"
#include "line/api/mapqn/mapqn_params.h"
#include "line/api/mapqn/mapqn_qr_bounds_bas.h"
#include "line/num/number.h"
#include "line/util/matrix.h"

using line::Matrix;
using line::mapqn::mapqn_qr_bounds_bas;
using line::mapqn::MapqnSense;
using line::mapqn::QrBasIndex;
using line::mapqn::QrBasParams;
using line::mapqn::QrBasResult;

namespace {

template <class T>
T num(double v) {
    return line::num_traits<T>::from_double(v);
}

/**
 * A tiny BAS network: three single-phase queues, N = 2, the finite queue is 0
 * with capacity 1, three blocking configurations. Converted to 0-based, so the
 * blocking-order entries MM and MM1 carry -1 where MATLAB writes 0 for
 * "absent".
 *
 * Deliberately smaller than example_bas_small.m. This model is
 * MR * B^2 + sum_i K(i) columns with B = (N+1) sum_i K(i), so that example is
 * 3894 columns, which lp::simplex_solve's DENSE tableau does not finish in
 * reasonable time. Here it is 246. See the cost note in the header.
 */
template <class T>
QrBasParams<T> instance_small() {
    QrBasParams<T> p;
    p.M = 3;
    p.N = 2;
    p.f = 0;
    p.F.push_back(1);
    p.F.push_back(2);
    p.F.push_back(2);
    p.K.assign(3, 1);
    const double rate[3] = {1.0, 0.8, 0.6};
    for (int i = 0; i < 3; ++i) {
        p.mu.push_back(Matrix<T>(1, 1, num<T>(rate[i])));
        p.v.push_back(Matrix<T>(1, 1, T()));
    }
    p.r = Matrix<T>(3, 3, T());
    p.r(0, 1) = num<T>(0.5);
    p.r(0, 2) = num<T>(0.5);
    p.r(1, 0) = num<T>(0.5);
    p.r(1, 2) = num<T>(0.5);
    p.r(2, 0) = num<T>(0.5);
    p.r(2, 1) = num<T>(0.5);

    p.MR = 3;
    const int bb[3][3] = {{0, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    const int mm[3][2] = {{-1, -1}, {1, -1}, {2, -1}};
    for (int m = 0; m < 3; ++m) {
        std::vector<int> row;
        for (int i = 0; i < 3; ++i) row.push_back(bb[m][i]);
        p.BB.push_back(row);
        std::vector<int> ord;
        ord.push_back(mm[m][0]);
        ord.push_back(mm[m][1]);
        p.MM.push_back(ord);
        std::vector<int> ext(3, -1);
        p.MM1.push_back(ext);
    }
    p.ZZ.push_back(0);
    p.ZZ.push_back(1);
    p.ZZ.push_back(1);
    p.ZM = 1;
    return p;
}

}  // namespace

TEST_CASE("BAS tiny instance: MATLAB reference bounds for every queue") {
    // qrf_bas.m, R2026a, re-measured 2026-08-22 on this instance after the
    // objective moved from occupancy to UTILIZATION (it now optimises the e
    // variables, which UEFF restricts to the configurations where the station
    // is not blocked). Superseding values measured 2026-08-01 under R2025a:
    //   ref_min = {0.374247090377, 0.512782760339, 0.670647149738}
    //   ref_max = {0.383667177081, 0.531587056203, 0.687993574882}
    //
    // Queue 0 is f, the finite-capacity queue. f is the blocking DESTINATION and
    // is never itself blocked, so occupancy and utilization coincide there and
    // its two constants are UNCHANGED to every digit -- which is the check that
    // this regeneration tracked a real change of quantity rather than drift.
    //
    // The MINIMA are where the distinction bites, since occupancy >= utilization
    // pointwise: queue 1 falls 0.5128 -> 0.4196 and queue 2 0.6706 -> 0.5701.
    // The maxima move by 1.3e-8 and 8.0e-7, i.e. LP tolerance, because the
    // maximising vertex carries no blocking mass at those stations.
    const double ref_min[3] = {0.374247090377, 0.419623624125, 0.570107980685};
    const double ref_max[3] = {0.383667177081, 0.531587043308, 0.687992773604};
    const QrBasParams<double> p = instance_small<double>();
    for (int i = 0; i < 3; ++i) {
        const QrBasResult<double> lo = mapqn_qr_bounds_bas(p, i, MapqnSense::Min);
        const QrBasResult<double> hi = mapqn_qr_bounds_bas(p, i, MapqnSense::Max);
        REQUIRE(lo.ok);
        REQUIRE(hi.ok);
        CHECK(lo.objective == doctest::Approx(ref_min[i]).epsilon(1e-5));
        CHECK(hi.objective == doctest::Approx(ref_max[i]).epsilon(1e-5));
        CHECK(lo.objective <= hi.objective + 1e-12);
        // A vacuous bound is the [0,1] box; an inventory gap shows here first.
        CHECK(hi.objective < 1.0 - 1e-6);
        CHECK(lo.objective > 1e-6);
    }
}

TEST_CASE("BAS: example_bas_small.m, the instance the dense tableau cannot carry") {
    // matlab/lib/qrf/example_bas_small.m: M = 3, N = 5, f = 0 with capacity 3,
    // two phases each, MR = 3. 3894 columns, which util/simplex.h's dense
    // tableau does not finish; via lp_solve it goes to HiGHS and returns in
    // well under a second. qrf_bas.m, R2025a, 2026-08-01: min 0.571555060034,
    // max 0.797487083957. This case is what proves the sparse seam is wired.
    QrBasParams<double> p;
    p.M = 3;
    p.N = 5;
    p.f = 0;
    p.F.push_back(3);
    p.F.push_back(5);
    p.F.push_back(5);
    p.K.assign(3, 2);
    for (int i = 0; i < 3; ++i) {
        Matrix<double> mu(2, 2);
        mu(0, 0) = 1.0;
        mu(0, 1) = 0.1;
        mu(1, 0) = 0.1;
        mu(1, 1) = 0.5;
        p.mu.push_back(mu);
        p.v.push_back(Matrix<double>(2, 2, 0.0));
    }
    p.r = Matrix<double>(3, 3, 0.0);
    p.r(0, 1) = 0.5;
    p.r(0, 2) = 0.5;
    p.r(1, 0) = 0.5;
    p.r(1, 2) = 0.5;
    p.r(2, 0) = 0.5;
    p.r(2, 1) = 0.5;
    p.MR = 3;
    const int bb[3][3] = {{0, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    const int mm[3][2] = {{-1, -1}, {1, -1}, {2, -1}};
    for (int m = 0; m < 3; ++m) {
        std::vector<int> row;
        for (int i = 0; i < 3; ++i) row.push_back(bb[m][i]);
        p.BB.push_back(row);
        std::vector<int> ord;
        ord.push_back(mm[m][0]);
        ord.push_back(mm[m][1]);
        p.MM.push_back(ord);
        p.MM1.push_back(std::vector<int>(3, -1));
    }
    p.ZZ.push_back(0);
    p.ZZ.push_back(1);
    p.ZZ.push_back(1);
    p.ZM = 1;

    const QrBasResult<double> lo = mapqn_qr_bounds_bas(p, 0, MapqnSense::Min);
    const QrBasResult<double> hi = mapqn_qr_bounds_bas(p, 0, MapqnSense::Max);
    REQUIRE(lo.ok);
    REQUIRE(hi.ok);
    CHECK(lo.num_vars == 3894u);
    CHECK(lo.objective == doctest::Approx(0.571555060034).epsilon(1e-6));
    CHECK(hi.objective == doctest::Approx(0.797487083957).epsilon(1e-6));
}

TEST_CASE("BAS: the variable layout is a bijection onto its own range") {
    std::vector<int> K;
    K.push_back(2);
    K.push_back(1);
    const QrBasIndex x(2, 2, K, 2);
    std::vector<char> seen(x.num_vars(), 0);
    for (int j = 0; j < 2; ++j)
        for (int nj = 0; nj <= 2; ++nj)
            for (int kj = 0; kj < K[j]; ++kj)
                for (int i = 0; i < 2; ++i)
                    for (int ni = 0; ni <= 2; ++ni)
                        for (int hi = 0; hi < K[i]; ++hi)
                            for (int m = 0; m < 2; ++m) {
                                const std::size_t idx = x.p2(j, nj, kj, i, ni, hi, m);
                                REQUIRE(idx < x.num_vars());
                                CHECK(seen[idx] == 0);
                                seen[idx] = 1;
                            }
    for (int i = 0; i < 2; ++i)
        for (int ki = 0; ki < K[i]; ++ki) {
            REQUIRE(x.e(i, ki) < x.num_vars());
            CHECK(seen[x.e(i, ki)] == 0);
            seen[x.e(i, ki)] = 1;
        }
    std::size_t used = 0;
    for (std::size_t i = 0; i < seen.size(); ++i) used += static_cast<std::size_t>(seen[i]);
    CHECK(used == x.num_vars());
}

TEST_CASE("BAS: argument validation") {
    const QrBasParams<double> p = instance_small<double>();
    CHECK_THROWS(mapqn_qr_bounds_bas(p, -1, MapqnSense::Min));
    CHECK_THROWS(mapqn_qr_bounds_bas(p, 3, MapqnSense::Min));
    QrBasParams<double> badf = p;
    badf.f = 3;
    CHECK_THROWS(mapqn_qr_bounds_bas(badf, 0, MapqnSense::Min));
    QrBasParams<double> badmr = p;
    badmr.MR = 2;
    CHECK_THROWS(mapqn_qr_bounds_bas(badmr, 0, MapqnSense::Min));
    QrBasParams<double> badcap = p;
    badcap.F[0] = 99;
    CHECK_THROWS(mapqn_qr_bounds_bas(badcap, 0, MapqnSense::Min));
}
