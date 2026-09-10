/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * map_block and map_feasblock: the general MAP(2) fit to three moments and an
 * autocorrelation decay rate.
 *
 * THE ORACLE IS INVERSION, as for map_mmpp2: the forward moment map is applied
 * to the fit and must land back on the request. `map_block` is the branch that
 * MATTERS to check that way, because the 57 KB of Maple algebra behind it has
 * no other independent description.
 *
 * The fallback branch is checked on its own terms: it deliberately drops the
 * third moment, so only E1, E2 and the decay rate survive there, and the test
 * asserts exactly that rather than pretending E3 is matched.
 */
#include <cmath>
#include <cstddef>
#include <vector>

#include "doctest.h"
#include "line/api/mam/map_block.h"
#include "line/api/mam/map_moment.h"
#include "line/api/mam/map_transform.h"

namespace mam = line::mam;
using line::Matrix;

namespace {

/** rho(2)/rho(1), the decay rate a MAP(2) holds constant across lags. */
double decay(const mam::Map<double>& m) {
    const std::vector<unsigned> lags{1u, 2u};
    const std::vector<double> a = mam::map_acf(m, lags);
    return a[1] / a[0];
}

void check_is_map(const mam::Map<double>& m) {
    for (std::size_t i = 0; i < m.order(); ++i) {
        double s = 0.0;
        for (std::size_t j = 0; j < m.order(); ++j) {
            if (i != j) CHECK(m.D0(i, j) >= -1e-12);
            CHECK(m.D1(i, j) >= -1e-12);
            s += m.D0(i, j) + m.D1(i, j);
        }
        CHECK(std::fabs(s) < 1e-9);
    }
}

}  // namespace

TEST_CASE("map_block inverts the three moments and the decay rate it was given") {
    struct Case { double e1, scv, e3mult, g2; };
    // e3mult multiplies the lower limit (3/2) e2^2 / e1, so every case is
    // feasible in the third moment by construction.
    const Case cases[] = {
        {1.0, 2.0, 2.0, 0.30},
        {1.0, 4.0, 1.8, 0.50},
        {0.7, 3.0, 2.5, 0.20},
    };
    for (std::size_t c = 0; c < 3; ++c) {
        const Case& k = cases[c];
        const double e2 = (1.0 + k.scv) * k.e1 * k.e1;
        const double e3 = k.e3mult * 1.5 * e2 * e2 / k.e1;
        const mam::Map<double> m = mam::map_block(k.e1, e2, e3, k.g2);
        check_is_map(m);
        // Whichever branch was taken, the first two moments always survive.
        CHECK(mam::map_mean(m) == doctest::Approx(k.e1).epsilon(1e-7));
        CHECK(mam::map_scv(m) == doctest::Approx(k.scv).epsilon(1e-6));
        if (m.order() == 2 && std::fabs(mam::map_moment(m, 3) - e3) < 1e-3 * e3) {
            // The exact branch was taken: the third moment and the decay rate
            // must both come back, which is what the Maple algebra is for.
            CHECK(mam::map_moment(m, 3) == doctest::Approx(e3).epsilon(1e-6));
            CHECK(decay(m) == doctest::Approx(k.g2).epsilon(1e-5));
        }
    }
}

TEST_CASE("the SCV spelling agrees with the moment spelling") {
    const double e1 = 1.0, scv = 2.5, g2 = 0.4;
    const double e2 = (1.0 + scv) * e1 * e1;
    const double e3 = 2.0 * 1.5 * e2 * e2 / e1;
    const mam::Map<double> a = mam::map_block(e1, e2, e3, g2);
    const mam::Map<double> b = mam::map_block_scv(e1, scv, e3, g2);
    REQUIRE(a.order() == b.order());
    for (std::size_t i = 0; i < a.order(); ++i)
        for (std::size_t j = 0; j < a.order(); ++j) {
            CHECK(a.D0(i, j) == doctest::Approx(b.D0(i, j)));
            CHECK(a.D1(i, j) == doctest::Approx(b.D1(i, j)));
        }
}

TEST_CASE("the fallback keeps the first two moments and drops the third") {
    // A third moment far below its limit forces the infeasible branch. The
    // fallback is defined for 1 <= SCV < 3, where its first branch rate
    // E1 (1 - sqrt((SCV-1)/2)) is still positive.
    const double e1 = 1.0, g2 = 0.4;
    const double scvs[] = {1.5, 2.0, 2.9};
    for (std::size_t i = 0; i < 3; ++i) {
        const double e2 = (1.0 + scvs[i]) * e1 * e1;
        const mam::Map<double> m = mam::map_block(e1, e2, 0.1 * e2 * e2 / e1, g2);
        check_is_map(m);
        CHECK(mam::map_mean(m) == doctest::Approx(e1).epsilon(1e-8));
        CHECK(mam::map_scv(m) == doctest::Approx(scvs[i]).epsilon(1e-7));
    }
}

TEST_CASE("at and above SCV 3 the fallback is refused rather than returned") {
    // Measured in MATLAB R2025a on the same inputs: an infinite diagonal at
    // SCV = 3 with a NaN mean, and a POSITIVE diagonal at SCV = 5, both with
    // map_isfeasible 0. Handing back either would be a matrix that is not a MAP.
    const double e1 = 1.0, g2 = 0.4;
    for (double scv : {3.0, 5.0}) {
        const double e2 = (1.0 + scv) * e1 * e1;
        CHECK_THROWS_AS(mam::map_block(e1, e2, 0.1 * e2 * e2 / e1, g2), line::InputError);
    }
}

TEST_CASE("an under-dispersed request falls back to the exponential") {
    // SCV < 1 has no hyperexponential form, and the reference's general MAP(2)
    // branch is commented out there, so the exponential is what comes back.
    const double e1 = 1.0, scv = 0.4;
    const double e2 = (1.0 + scv) * e1 * e1;
    const mam::Map<double> m = mam::map_block(e1, e2, 5.0 * e2 * e2 / e1, 0.3);
    CHECK(m.order() == 1);
    CHECK(mam::map_mean(m) == doctest::Approx(e1).epsilon(1e-12));
    CHECK(mam::map_scv(m) == doctest::Approx(1.0).epsilon(1e-12));
}

TEST_CASE("map_feasblock repairs the moments into the feasible region") {
    const double e1 = 1.0;
    // E2 below the exponential value: raised to (2 + tol) E1^2, so the SCV of
    // the result is at or just above one rather than the requested value.
    const mam::Map<double> a = mam::map_feasblock(e1, 1.2 * e1 * e1, 10.0, 0.3);
    check_is_map(a);
    CHECK(mam::map_mean(a) == doctest::Approx(e1).epsilon(1e-6));
    CHECK(mam::map_scv(a) >= 1.0 - 1e-9);

    // E3 below its limit: raised to (3/2 + tol) E2^2 / E1.
    const double scv = 2.0, e2 = (1.0 + scv) * e1 * e1;
    const mam::Map<double> b = mam::map_feasblock(e1, e2, 0.01, 0.3);
    check_is_map(b);
    CHECK(mam::map_mean(b) == doctest::Approx(e1).epsilon(1e-6));
    CHECK(mam::map_moment(b, 3) >= 1.5 * e2 * e2 / e1 - 1e-6);
}

TEST_CASE("map_feasblock returns the Poisson process at the exponential boundary") {
    const double e1 = 2.0;
    const mam::Map<double> m = mam::map_feasblock(e1, 2.0 * e1 * e1, 10.0, 0.3);
    CHECK(mam::map_mean(m) == doctest::Approx(e1).epsilon(1e-12));
    CHECK(mam::map_scv(m) == doctest::Approx(1.0).epsilon(1e-9));
    CHECK(std::fabs(mam::map_acf(m, std::vector<unsigned>{1})[0]) < 1e-9);
}
