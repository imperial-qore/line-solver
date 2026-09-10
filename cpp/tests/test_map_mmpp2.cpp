/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * map_mmpp2 and map2_fit_idc.
 *
 * THE ORACLE IS INVERSION. map_mmpp2 is the closed-form inverse of the moment
 * and autocorrelation map, so the test applies the FORWARD map to its output
 * and checks it lands back on the request: map_mean, map_scv, map_skew and
 * map_acf are computed by code that shares nothing with the fit. A transcription
 * error anywhere in the 18 KB of Maple algebra moves at least one of the four.
 *
 * The structural checks matter as much: the result must be an MMPP(2), i.e. D1
 * diagonal and non-negative, D0 with non-negative off-diagonals, rows of
 * D0 + D1 summing to zero.
 */
#include <cmath>
#include <cstddef>

#include "doctest.h"
#include "line/api/mam/map2_fit_idc.h"
#include "line/api/mam/map_mmpp2.h"
#include "line/api/mam/map_moment.h"
#include "line/api/mam/map_transform.h"

namespace mam = line::mam;
using line::Matrix;

namespace {

void check_is_mmpp2(const mam::Map<double>& m) {
    REQUIRE(m.order() == 2);
    CHECK(m.D1(0, 1) == doctest::Approx(0.0));  // D1 diagonal: no phase change on arrival
    CHECK(m.D1(1, 0) == doctest::Approx(0.0));
    CHECK(m.D1(0, 0) >= -1e-12);
    CHECK(m.D1(1, 1) >= -1e-12);
    CHECK(m.D0(0, 1) >= -1e-12);  // the modulating rates
    CHECK(m.D0(1, 0) >= -1e-12);
    for (std::size_t i = 0; i < 2; ++i) {
        double s = 0.0;
        for (std::size_t j = 0; j < 2; ++j) s += m.D0(i, j) + m.D1(i, j);
        CHECK(std::fabs(s) < 1e-9);
    }
}

}  // namespace

TEST_CASE("map_mmpp2 inverts the moments it was given") {
    // Feasible interior points: SCV > 1, a skewness above the minimum, and a
    // lag-1 autocorrelation inside [0, (1-1/SCV)/2].
    struct Case { double mean, scv, skew, acf1; };
    const Case cases[] = {
        {1.0, 2.0, 4.0, 0.10},
        {0.5, 4.0, 12.0, 0.20},
        {2.0, 1.5, 3.0, 0.05},
    };
    for (std::size_t c = 0; c < 3; ++c) {
        const Case& k = cases[c];
        const mam::Map<double> m = mam::map_mmpp2(k.mean, k.scv, k.skew, k.acf1);
        check_is_mmpp2(m);
        CHECK(mam::map_mean(m) == doctest::Approx(k.mean).epsilon(1e-8));
        CHECK(mam::map_scv(m) == doctest::Approx(k.scv).epsilon(1e-7));
        CHECK(mam::map_skew(m) == doctest::Approx(k.skew).epsilon(1e-6));
        CHECK(mam::map_acf(m, std::vector<unsigned>{1})[0] ==
              doctest::Approx(k.acf1).epsilon(1e-6));
    }
}

TEST_CASE("the ACF1 = -1 sentinel asks for the maximum feasible autocorrelation") {
    const double mean = 1.0, scv = 3.0, skew = 6.0;
    const mam::Map<double> m = mam::map_mmpp2(mean, scv, skew, -1.0);
    check_is_mmpp2(m);
    CHECK(mam::map_mean(m) == doctest::Approx(mean).epsilon(1e-8));
    CHECK(mam::map_scv(m) == doctest::Approx(scv).epsilon(1e-7));
    // The bound is (1 - 1/SCV)/2; the sentinel must approach it from below.
    const double bound = 0.5 * (1.0 - 1.0 / scv);
    const double got = mam::map_acf(m, std::vector<unsigned>{1})[0];
    CHECK(got > 0.9 * bound);
    CHECK(got <= bound + 1e-9);
}

TEST_CASE("the SKEW = -1 sentinel asks for the minimum third moment") {
    const mam::Map<double> m = mam::map_mmpp2(1.0, 3.0, -1.0, 0.1);
    check_is_mmpp2(m);
    CHECK(mam::map_mean(m) == doctest::Approx(1.0).epsilon(1e-8));
    CHECK(mam::map_scv(m) == doctest::Approx(3.0).epsilon(1e-7));
    // e3 sits just above its lower limit (3/2) e2^2 / e1.
    const double e1 = 1.0, e2 = (1.0 + 3.0) * e1 * e1;
    const double e3 = mam::map_moment(m, 3);
    CHECK(e3 >= 1.5 * e2 * e2 / e1);
    CHECK(e3 < 1.51 * e2 * e2 / e1);
}

TEST_CASE("map_mmpp2 refuses what an MMPP(2) cannot represent, by name") {
    CHECK_THROWS_AS(mam::map_mmpp2(1.0, 0.5, 2.0, 0.1), line::InputError);   // under-dispersed
    CHECK_THROWS_AS(mam::map_mmpp2(1.0, 1.0, 2.0, 0.1), line::InputError);   // Poisson boundary
    CHECK_THROWS_AS(mam::map_mmpp2(1.0, 2.0, 4.0, -0.2), line::InputError);  // negative acf
    CHECK_THROWS_AS(mam::map_mmpp2(1.0, 2.0, 4.0, 0.9), line::InputError);   // acf above bound
    CHECK_THROWS_AS(mam::map_mmpp2(1.0, 2.0, 0.1, 0.1), line::InputError);   // e3 below minimum
}

TEST_CASE("map2_fit_idc matches the index of dispersion it was asked for") {
    // SCV 3 and I 6 are inside the MAP(2) region; g2 = (I - SCV)/(I - 1).
    const double e1 = 1.0, scv = 3.0;
    const double e2 = (1.0 + scv) * e1 * e1;
    const double e3 = 3.0 * e2 * e2 / e1;
    const double I = 6.0;
    const mam::Map2FitIdcResult<double> r = mam::map2_fit_idc(e1, e2, e3, I);
    CHECK(r.status <= 3);
    CHECK(mam::map_mean(r.map) == doctest::Approx(e1).epsilon(1e-8));
    CHECK(mam::map_scv(r.map) == doctest::Approx(scv).epsilon(1e-6));
    if (r.status == 0)
        CHECK(mam::map_idc(r.map) == doctest::Approx(I).epsilon(1e-6));
}

TEST_CASE("map2_fit_idc returns the exponential where burstiness is not representable") {
    // The reference's Norton's-theorem rule: SCV <= 1 or I < SCV must give an
    // exponential, not a MAP(2) fitted to a marginal that is not independent of
    // the rest of the model.
    const double e1 = 1.0;
    const mam::Map2FitIdcResult<double> a = mam::map2_fit_idc(e1, 2.0 * e1 * e1, 6.0, 5.0);
    CHECK(a.status == 1);
    CHECK(a.map.order() == 1);
    CHECK(mam::map_mean(a.map) == doctest::Approx(e1).epsilon(1e-12));

    const double scv = 3.0, e2 = (1.0 + scv) * e1 * e1;
    const mam::Map2FitIdcResult<double> b = mam::map2_fit_idc(e1, e2, 3.0 * e2 * e2 / e1, 2.0);
    CHECK(b.status == 1);  // I = 2 < SCV = 3
    CHECK(b.map.order() == 1);
}
