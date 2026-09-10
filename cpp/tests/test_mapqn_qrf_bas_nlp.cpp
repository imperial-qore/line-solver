/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * mapqn_qrf_bas_mmi / _mem / _bethe: the nonlinear bounds on the BAS-blocking
 * polytope.
 *
 * THE POINT OF THESE ENTRY POINTS IS THAT THEY SHARE THE LP'S POLYTOPE, so
 * that is what the tests check: the answer must lie inside the linear bounds
 * `mapqn_qr_bounds_bas` reports for the same instance, since it is a feasible
 * point of the same set. A port that re-derived the constraint families and got
 * one wrong would show up as an answer outside that interval, which no amount
 * of comparing the objective against itself would catch.
 *
 * MEM and MMI are again different KINDS of problem -- MEM convex, MMI not (see
 * mapqn_qrf_common.h) -- so MEM is expected to be reproducible and MMI is
 * checked for feasibility and for the properties a utilization has, not against
 * a number.
 *
 * THE ITERATION CAP IS SMALL ON PURPOSE, and it costs the tests nothing they
 * assert. Every Frank-Wolfe iterate is a convex combination of two points of
 * the polytope, hence feasible, so a truncated run returns a WORSE bound but
 * never an invalid one -- and feasibility is exactly what is being checked. The
 * BAS polytope here is 246 columns and each iteration is one LP plus a line
 * search, so the uncapped default would make this file minutes long for no
 * additional assurance.
 */
#include <cstddef>
#include <vector>

#include "doctest.h"
#include "line/api/mapqn/mapqn_qr_bounds_bas.h"
#include "line/api/mapqn/mapqn_qrf_bas_nlp.h"

using namespace line;
using namespace line::mapqn;

namespace {

/**
 * The same tiny BAS network `test_mapqn_qr_bounds_bas.cpp` uses: three
 * single-phase queues, N = 2, the finite queue is 0 with capacity 1, three
 * blocking configurations. Small on purpose -- 246 columns.
 */
QrBasParams<double> instance_small() {
    QrBasParams<double> p;
    p.M = 3;
    p.N = 2;
    p.f = 0;
    p.F.push_back(1);
    p.F.push_back(2);
    p.F.push_back(2);
    p.K.assign(3, 1);
    const double rate[3] = {1.0, 0.8, 0.6};
    for (int i = 0; i < 3; ++i) {
        p.mu.push_back(Matrix<double>(1, 1, rate[i]));
        p.v.push_back(Matrix<double>(1, 1, 0.0));
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
    return p;
}

/** See the header note: feasibility holds at every iterate. */
const unsigned kIter = 6;

}  // namespace

TEST_CASE("the MEM bound is a feasible point of the LP's own polytope") {
    const QrBasParams<double> p = instance_small();
    const QrfMetrics<double> r = mapqn_qrf_bas_mem(p, kIter);

    REQUIRE(r.UN.size() == 3u);
    REQUIRE(r.QN.size() == 3u);
    for (int i = 0; i < 3; ++i) {
        // The polytope is shared, so the nonlinear answer must sit inside the
        // linear bounds on the same station. That is the check that a
        // re-derived constraint family would fail.
        const QrBasResult<double> lo = mapqn_qr_bounds_bas(p, i, MapqnSense::Min);
        const QrBasResult<double> hi = mapqn_qr_bounds_bas(p, i, MapqnSense::Max);
        REQUIRE(lo.ok);
        REQUIRE(hi.ok);
        CHECK(r.UN[i] >= lo.objective - 1e-6);
        CHECK(r.UN[i] <= hi.objective + 1e-6);
        // A utilization is a probability, and a queue length is bounded by the
        // station's own capacity.
        CHECK(r.UN[i] >= -1e-9);
        CHECK(r.UN[i] <= 1.0 + 1e-9);
        CHECK(r.QN[i] >= -1e-9);
        CHECK(r.QN[i] <= static_cast<double>(p.F[i]) + 1e-6);
    }
}

TEST_CASE("the MMI bound is feasible too, on the same polytope") {
    const QrBasParams<double> p = instance_small();
    const QrfMetrics<double> r = mapqn_qrf_bas_mmi(p, kIter);

    REQUIRE(r.UN.size() == 3u);
    for (int i = 0; i < 3; ++i) {
        const QrBasResult<double> lo = mapqn_qr_bounds_bas(p, i, MapqnSense::Min);
        const QrBasResult<double> hi = mapqn_qr_bounds_bas(p, i, MapqnSense::Max);
        REQUIRE(lo.ok);
        REQUIRE(hi.ok);
        CHECK(r.UN[i] >= lo.objective - 1e-6);
        CHECK(r.UN[i] <= hi.objective + 1e-6);
        CHECK(r.QN[i] >= -1e-9);
        CHECK(r.QN[i] <= static_cast<double>(p.F[i]) + 1e-6);
    }
}

TEST_CASE("MEM is reproducible where MMI need not be") {
    // MEM is convex on the polytope, so two runs from the same input give the
    // same point. MMI is not convex, so this is asserted only of MEM -- see
    // mapqn_qrf_common.h for the measurement behind that distinction.
    const QrBasParams<double> p = instance_small();
    const QrfMetrics<double> a = mapqn_qrf_bas_mem(p, kIter);
    const QrfMetrics<double> b = mapqn_qrf_bas_mem(p, kIter);
    for (int i = 0; i < 3; ++i) {
        CHECK(a.UN[i] == doctest::Approx(b.UN[i]).epsilon(1e-12));
        CHECK(a.QN[i] == doctest::Approx(b.QN[i]).epsilon(1e-12));
    }
}

TEST_CASE("the BETHE bound is feasible too, on the same polytope") {
    const QrBasParams<double> p = instance_small();
    const QrfMetrics<double> r = mapqn_qrf_bas_bethe(p, kIter);

    REQUIRE(r.UN.size() == 3u);
    for (int i = 0; i < 3; ++i) {
        const QrBasResult<double> lo = mapqn_qr_bounds_bas(p, i, MapqnSense::Min);
        const QrBasResult<double> hi = mapqn_qr_bounds_bas(p, i, MapqnSense::Max);
        REQUIRE(lo.ok);
        REQUIRE(hi.ok);
        CHECK(r.UN[i] >= lo.objective - 1e-6);
        CHECK(r.UN[i] <= hi.objective + 1e-6);
        CHECK(r.QN[i] >= -1e-9);
        CHECK(r.QN[i] <= static_cast<double>(p.F[i]) + 1e-6);
    }
}

TEST_CASE("BETHE is reproducible, as its convexity says it should be") {
    // lambda = 1/M is chosen so the tree-reweighted entropy is concave and the
    // program convex, which is the property MMI does not have. The BAS
    // polytope adds the blocking families on top of the marginal ones, so this
    // is asserted as a MEASUREMENT and not inherited from the no-blocking twin.
    const QrBasParams<double> p = instance_small();
    const QrfMetrics<double> a = mapqn_qrf_bas_bethe(p, kIter);
    const QrfMetrics<double> b = mapqn_qrf_bas_bethe(p, kIter);
    for (int i = 0; i < 3; ++i) {
        CHECK(a.UN[i] == doctest::Approx(b.UN[i]).epsilon(1e-12));
        CHECK(a.QN[i] == doctest::Approx(b.QN[i]).epsilon(1e-12));
    }
}

TEST_CASE("the BETHE term sets span the idle cells that MMI and MEM skip") {
    // Both of its blocks run from n = 0, which is what separates it from
    // lam*MMI + MEM. A regression restoring only one block would still produce
    // a plausible number, so the ranges are asserted directly.
    const QrBasParams<double> p = instance_small();
    const QrBasIndex x(p.M, p.N, p.K, p.MR);

    std::vector<std::size_t> ij0, ii0, jj0, ij1, ii1, jj1;
    qrfbas::mmi_terms(x, p.F, &ij0, &ii0, &jj0, 0);
    qrfbas::mmi_terms(x, p.F, &ij1, &ii1, &jj1, 1);
    CHECK(ij0.size() > ij1.size());
    CHECK(qrfbas::mem_terms(x, p.F, 0).size() > qrfbas::mem_terms(x, p.F, 1).size());
}

TEST_CASE("the parameter validation is the LP bound's own") {
    QrBasParams<double> bad = instance_small();
    bad.f = 7;  // out of range
    CHECK_THROWS_AS(mapqn_qrf_bas_mem(bad), line::InputError);

    QrBasParams<double> badmr = instance_small();
    badmr.MR = 0;
    CHECK_THROWS_AS(mapqn_qrf_bas_mmi(badmr), line::InputError);

    QrBasParams<double> badcap = instance_small();
    badcap.F[0] = 9;  // above N
    CHECK_THROWS_AS(mapqn_qrf_bas_mem(badcap), line::InputError);

    QrBasParams<double> badbethe = instance_small();
    badbethe.f = 7;
    CHECK_THROWS_AS(mapqn_qrf_bas_bethe(badbethe), line::InputError);
}
