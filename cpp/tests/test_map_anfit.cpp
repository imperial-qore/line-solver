/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * map_anfit: the Andersen-Nielsen superposition of interrupted Poisson
 * processes fitted to a Hurst parameter.
 *
 * ORACLES. The construction has no closed-form answer to compare against, so
 * the checks are the properties it is built to have:
 *  - the result is a valid MAP of order 2^d (or 2^(d+1) with a Poisson
 *    remainder), since superposing d two-state IPPs doubles the state space
 *    each time;
 *  - the arrival rate is the target ls, which is what the eta normalization
 *    exists to achieve;
 *  - the switching ladder is geometric with the documented ratio, and the
 *    autocorrelation is positive and DECAYING, which is the long-range
 *    dependence the fit is for;
 *  - a larger Hurst parameter gives slower decay, the ordering being the whole
 *    point of the parameter.
 */
#include <cmath>
#include <cstddef>
#include <vector>

#include "doctest.h"
#include "line/api/mam/map_anfit.h"
#include "line/api/mam/map_moment.h"
#include "line/api/mam/map_transform.h"

namespace mam = line::mam;

namespace {

void check_is_map(const mam::Map<double>& m) {
    for (std::size_t i = 0; i < m.order(); ++i) {
        double s = 0.0;
        for (std::size_t j = 0; j < m.order(); ++j) {
            if (i != j) CHECK(m.D0(i, j) >= -1e-9);
            CHECK(m.D1(i, j) >= -1e-9);
            s += m.D0(i, j) + m.D1(i, j);
        }
        CHECK(std::fabs(s) < 1e-8);
    }
}

}  // namespace

TEST_CASE("map_anfit builds a valid MAP at the requested arrival rate") {
    const double ls = 1.0, rho = 0.3, H = 0.7;
    const mam::MapAnfitResult<double> r = mam::map_anfit(ls, rho, H, 3.0, 3);
    check_is_map(r.map);
    CHECK(r.d >= 3);
    // d IPPs plus the Poisson remainder: 2^d * 1.
    CHECK(r.map.order() == (std::size_t(1) << r.d));
    CHECK(mam::map_lambda(r.map) == doctest::Approx(ls).epsilon(1e-6));
}

TEST_CASE("the switching ladder is geometric") {
    const mam::MapAnfitResult<double> r = mam::map_anfit(1.0, 0.3, 0.7, 3.0, 4);
    REQUIRE(r.switching.size() == r.d);
    CHECK(r.switching[0] == doctest::Approx(0.8).epsilon(1e-12));
    // k(2,i) = a^(1-i) k(2,1), so consecutive ratios are constant.
    const double ratio = r.switching[1] / r.switching[0];
    for (std::size_t i = 1; i + 1 < r.switching.size(); ++i)
        CHECK(r.switching[i + 1] / r.switching[i] == doctest::Approx(ratio).epsilon(1e-9));
    CHECK(ratio < 1.0);  // the ladder descends
}

TEST_CASE("the autocorrelation is positive and decaying") {
    const mam::MapAnfitResult<double> r = mam::map_anfit(1.0, 0.3, 0.75, 3.0, 3);
    std::vector<unsigned> lags;
    for (unsigned k = 1; k <= 6; ++k) lags.push_back(k);
    const std::vector<double> acf = mam::map_acf(r.map, lags);
    for (std::size_t i = 0; i < acf.size(); ++i) CHECK(acf[i] > 0.0);
    for (std::size_t i = 0; i + 1 < acf.size(); ++i) CHECK(acf[i + 1] <= acf[i] + 1e-12);
}

TEST_CASE("a larger Hurst parameter decays more slowly") {
    std::vector<unsigned> lags;
    for (unsigned k = 1; k <= 8; ++k) lags.push_back(k);
    const mam::MapAnfitResult<double> lo = mam::map_anfit(1.0, 0.3, 0.6, 3.0, 3);
    const mam::MapAnfitResult<double> hi = mam::map_anfit(1.0, 0.3, 0.85, 3.0, 3);
    const std::vector<double> alo = mam::map_acf(lo.map, lags);
    const std::vector<double> ahi = mam::map_acf(hi.map, lags);
    // Compare the ratio of the last lag to the first: slower decay is a larger
    // ratio, and that ordering is what the Hurst parameter controls.
    const double dlo = alo.back() / alo.front();
    const double dhi = ahi.back() / ahi.front();
    CHECK(dhi > dlo);
}

TEST_CASE("map_anfit refuses inputs the construction cannot serve") {
    CHECK_THROWS_AS(mam::map_anfit(1.0, 0.3, 0.7, 3.0, 1), line::InputError);
    CHECK_THROWS_AS(mam::map_anfit(0.0, 0.3, 0.7, 3.0, 3), line::InputError);
    CHECK_THROWS_AS(mam::map_anfit(-1.0, 0.3, 0.7, 3.0, 3), line::InputError);
}

TEST_CASE("the least-squares variant returns a valid MAP and refuses bad targets") {
    std::vector<unsigned> lags;
    for (unsigned k = 1; k <= 4; ++k) lags.push_back(k);
    // Target the autocorrelation of the deterministic fit itself, so a feasible
    // point exists and the search has somewhere to land.
    const mam::MapAnfitResult<double> base = mam::map_anfit(1.0, 0.3, 0.7, 3.0, 3);
    const std::vector<double> target = mam::map_acf(base.map, lags);

    const mam::MapAnfitResult<double> r =
        mam::map_anfit_lsq(1.0, 0.3, 0.7, 3.0, 3, target, lags, 20);
    check_is_map(r.map);
    CHECK(r.map.order() == base.map.order());

    std::vector<double> short_target(2, 0.1);
    CHECK_THROWS_AS(mam::map_anfit_lsq(1.0, 0.3, 0.7, 3.0, 3, short_target, lags),
                    line::InputError);
    CHECK_THROWS_AS(mam::map_anfit_lsq(1.0, 0.3, 0.7, 3.0, 3, std::vector<double>(),
                                       std::vector<unsigned>()),
                    line::InputError);
}
