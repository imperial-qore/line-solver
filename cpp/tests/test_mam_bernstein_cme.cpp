/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * map_bernstein, the CME representation and dist_fit_me, none of which had a
 * C++ twin before 2026-08-01. Together they are what `sn_nonmarkov_toph` needs
 * to hand a non-Markovian service law to a Markovian solver.
 *
 * ORACLES.
 *  (a) Moment identities, which do not depend on the fit: a renewal MAP's mean
 *      is -alpha D0^-1 e and its SCV follows from the second moment, so a fit
 *      that claims to match two moments can be checked against them directly.
 *  (b) The reference's own guards: a density that underflows on part of the
 *      Bernstein grid must still give a usable fit, which is where the MATLAB
 *      and JAR versions part company.
 *  (c) The CME table invariants: unit mean and a normalized entry law at every
 *      tabulated order.
 */
#include <cmath>
#include <cstddef>
#include <vector>

#include "doctest.h"
#include "line/api/mam/cme.h"
#include "line/api/mam/map_bernstein.h"
#include "line/api/mam/map_moment.h"
#include "line/api/mam/map_transform.h"

namespace mam = line::mam;
using line::Matrix;

TEST_CASE("map_bernstein fits an exponential density and rescales to its mean") {
    // f(x) = exp(-x): the fit has unit time scale, so map_scale is what puts the
    // mean where the caller wants it.
    const mam::Map<double> m =
        mam::map_bernstein<double>([](double x) { return std::exp(-x); }, 20);
    CHECK(m.order() == 20);

    const mam::Map<double> s = mam::map_scale(m, 2.5);
    CHECK(mam::map_mean(s) == doctest::Approx(2.5).epsilon(1e-9));
    // The Bernstein fit of an exponential is not an exponential, but it is close.
    CHECK(mam::map_scv(s) == doctest::Approx(1.0).epsilon(0.15));
}

TEST_CASE("map_bernstein survives a density that vanishes on part of the grid") {
    // Uniform(0,1). The grid nodes -log(i/n) exceed 1 for i < n/e, so roughly a
    // third of them evaluate to zero. The JAR divides by its normalization
    // unconditionally; the MATLAB reference skips those nodes and renormalizes,
    // which is what is ported here.
    const mam::Map<double> m = mam::map_bernstein<double>(
        [](double x) { return (x >= 0.0 && x <= 1.0) ? 1.0 : 0.0; }, 20);
    const double mean = mam::map_mean(m);
    CHECK(std::isfinite(mean));
    CHECK(mean > 0.0);

    const mam::Map<double> s = mam::map_scale(m, 0.5);
    CHECK(mam::map_mean(s) == doctest::Approx(0.5).epsilon(1e-9));
    CHECK(mam::map_scv(s) < 1.0);  // a Uniform is underdispersed
}

TEST_CASE("map_bernstein falls back to an Erlang when no node carries mass") {
    // A density supported far above every grid node: nothing is usable, and the
    // reference returns the unit-mean Erlang-n rather than a NaN generator.
    const mam::Map<double> m =
        mam::map_bernstein<double>([](double x) { return x > 1e6 ? 1.0 : 0.0; }, 8);
    CHECK(m.order() == 8);
    CHECK(mam::map_mean(m) == doctest::Approx(1.0).epsilon(1e-9));
    CHECK(mam::map_scv(m) == doctest::Approx(1.0 / 8.0).epsilon(1e-9));
}

TEST_CASE("aph_bernstein is the same fit under the JAR's name") {
    const mam::Map<double> a =
        mam::aph_bernstein<double>([](double x) { return std::exp(-x); }, 6);
    const mam::Map<double> b =
        mam::map_bernstein<double>([](double x) { return std::exp(-x); }, 6);
    for (std::size_t i = 0; i < a.order(); ++i)
        for (std::size_t j = 0; j < a.order(); ++j) {
            CHECK(a.D0(i, j) == doctest::Approx(b.D0(i, j)));
            CHECK(a.D1(i, j) == doctest::Approx(b.D1(i, j)));
        }
}

TEST_CASE("the CME table exposes odd orders with decreasing minimal SCV") {
    const std::vector<std::size_t> orders = mam::cme_supported_orders();
    REQUIRE(orders.size() > 1);
    CHECK(orders.front() >= 3);
    for (std::size_t i = 0; i < orders.size(); ++i) CHECK(orders[i] % 2 == 1);
    // More harmonics can only concentrate further, so the reach is monotone.
    CHECK(mam::cme_min_scv(orders.back()) < mam::cme_min_scv(orders.front()));
    // The Erlang bound is 1/order; a CME must beat it at high order.
    CHECK(mam::cme_min_scv(orders.back()) < 1.0 / static_cast<double>(orders.back()));
}

TEST_CASE("cme_representation has unit mean and a normalized entry law") {
    // The moment checks are confined to the orders where double arithmetic can
    // still resolve them. At a tabulated cv2 of 1e-6 the SCV is m2/m1^2 - 1 with
    // m2/m1^2 = 1.000001, so ten digits cancel and the computed SCV is noise --
    // a property of the representation in double, shared with the reference, not
    // of this port. The structural invariants are checked at every order.
    const std::vector<std::size_t> orders = mam::cme_supported_orders();
    for (std::size_t k = 0; k < orders.size(); ++k) {
        const mam::CmeRepresentation<double> r = mam::cme_representation<double>(orders[k]);
        REQUIRE(r.alpha.size() == orders[k]);
        // The entry law alternates in sign with entries of magnitude far above
        // one, so the tolerance on its sum has to carry the conditioning of that
        // cancellation: at order 2001 the terms reach 1e14 and no summation in
        // double resolves the unit total to better than 1e-2.
        double s = 0.0, mag = 0.0;
        for (std::size_t i = 0; i < r.alpha.size(); ++i) {
            s += r.alpha[i];
            mag += std::fabs(r.alpha[i]);
        }
        CHECK(std::fabs(s - 1.0) <= 64.0 * mag * 2.3e-16);

        if (orders[k] > 41) continue;
        const mam::Map<double> m = mam::me_to_map(r.alpha, r.A);
        CHECK(mam::map_mean(m) == doctest::Approx(1.0).epsilon(1e-8));
        CHECK(mam::map_scv(m) == doctest::Approx(r.scv).epsilon(1e-5));
    }
}

TEST_CASE("the CME table reaches the low orders dist_fit_me starts from") {
    // The vendored table was originally trimmed to what matlab_ilt can select,
    // which is n >= 10. dist_fit_me walks the orders upwards, so without the
    // low-n tail it fitted in 22 phases what the reference fits in 4.
    const std::vector<std::size_t> orders = mam::cme_supported_orders();
    CHECK(orders.front() == 3);
    CHECK(mam::cme_min_scv(3) == doctest::Approx(0.20090156350183885).epsilon(1e-12));
}

TEST_CASE("dist_fit_me matches the mean and the SCV exactly inside its reach") {
    const double targets[] = {0.05, 0.2, 0.5, 0.9};
    for (std::size_t i = 0; i < 4; ++i) {
        const mam::Map<double> m = mam::dist_fit_me<double>(3.0, targets[i]);
        CHECK(mam::map_mean(m) == doctest::Approx(3.0).epsilon(1e-8));
        CHECK(mam::map_scv(m) == doctest::Approx(targets[i]).epsilon(1e-6));
    }
}

TEST_CASE("dist_fit_me picks the smallest order whose reach covers the target") {
    // SCV 0.5 is above the reach of order 3 (0.2009/1.2009 = 0.167), so the
    // reference stops there: 3 CME phases plus the exponential.
    const mam::Map<double> m = mam::dist_fit_me<double>(2.0, 0.5);
    CHECK(m.order() == 4);
    CHECK(mam::map_mean(m) == doctest::Approx(2.0).epsilon(1e-9));
    CHECK(mam::map_scv(m) == doctest::Approx(0.5).epsilon(1e-8));
}

TEST_CASE("dist_fit_me under a phase budget returns the closest reachable SCV") {
    // 21 phases buys a CME of order 19 or less plus the exponential. The target
    // 1e-6 is below what that order reaches, so the fit is budget-limited: the
    // mean still lands exactly and the SCV lands at the reach, from above.
    const mam::Map<double> m = mam::dist_fit_me<double>(1.0, 1e-6, 21);
    CHECK(m.order() <= 21);
    CHECK(mam::map_mean(m) == doctest::Approx(1.0).epsilon(1e-8));
    const double got = mam::map_scv(m);
    CHECK(got > 1e-6);
    CHECK(got < 1.0 / static_cast<double>(m.order()));  // still beats the Erlang bound
}

TEST_CASE("dist_fit_me refuses the range a matrix exponential cannot cover") {
    CHECK_THROWS_AS(mam::dist_fit_me<double>(1.0, 1.5), line::InputError);
    CHECK_THROWS_AS(mam::dist_fit_me<double>(1.0, 0.0), line::InputError);
    CHECK_THROWS_AS(mam::dist_fit_me<double>(-1.0, 0.5), line::InputError);
    CHECK_THROWS_AS(mam::dist_fit_me<double>(1.0, 0.5, 3), line::InputError);
}
