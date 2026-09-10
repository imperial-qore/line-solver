/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * mamap2m_can1_coefficients and mamap2m_can2_coefficients.
 *
 * THE ORACLE IS MATLAB ITSELF: every value below was printed from R2025a at
 * twelve significant digits on the same inputs. That is possible here, and was
 * not for the fitters built on these tables, because the coefficients are a
 * pure function of (h1, h2, r1, r2) with no solver in the way.
 *
 * G(9) IS THE ONE DELIBERATE DIVERGENCE. MATLAB assigns G(10) twice in
 * consecutive lines, so the value intended for G(9) is discarded and G(9) stays
 * zero; the JAR writes it to index 9 and so does this port. The test asserts the
 * CORRECTED value and names the reference's zero, so the divergence cannot be
 * mistaken for a transcription slip later.
 */
#include <cmath>
#include <cstddef>
#include <vector>

#include "doctest.h"
#include "line/api/mam/mamap2m_coefficients.h"

namespace mam = line::mam;

TEST_CASE("the first canonical coefficients match MATLAB") {
    const mam::Mamap2mCoefficients<double> c =
        mam::mamap2m_can1_coefficients(1.0, 2.0, 0.5, 0.4);
    REQUIRE(c.G.size() == 15);
    const double G[15] = {0.375, 0.375, 0.25,  0.1875, 0.3,  0.075, 0.1125, 0.225,
                          0.1,   0.375, 0.975, 0.65,   0.75, 0.75,  0.5};
    for (std::size_t i = 0; i < 15; ++i) {
        if (i == 8) continue;  // checked below, against the corrected value
        CHECK(c.G[i] == doctest::Approx(G[i]).epsilon(1e-10));
    }
    // MATLAB prints 0 here because of the duplicated G(10) assignment.
    CHECK(c.G[8] == doctest::Approx(0.1).epsilon(1e-10));

    const double U[12] = {0.64, -2.56, 0.0,   2.56, 0.0,   0.0,
                          0.64, -2.176, -0.768, 2.56, 0.384, -0.768};
    REQUIRE(c.U.size() == 12);
    for (std::size_t i = 0; i < 12; ++i) CHECK(c.U[i] == doctest::Approx(U[i]).epsilon(1e-10));

    REQUIRE(c.Y.size() == 3);
    CHECK(std::fabs(c.Y[0]) < 1e-15);
    CHECK(std::fabs(c.Y[1]) < 1e-15);
    CHECK(c.Y[2] == doctest::Approx(-0.15).epsilon(1e-10));
}

TEST_CASE("the second canonical coefficients match MATLAB") {
    const mam::Mamap2mCoefficients<double> c =
        mam::mamap2m_can2_coefficients(1.0, 2.0, 0.5, 0.4);
    REQUIRE(c.G.size() == 14);
    const double E[14] = {0.230769230769, 0.461538461538, 0.307692307692, 0.369230769231,
                          0.0923076923077, 0.138461538462, 0.276923076923, 0.123076923077,
                          0.230769230769, 1.06153846154,  0.707692307692, 0.461538461538,
                          0.923076923077, 0.615384615385};
    for (std::size_t i = 0; i < 14; ++i) CHECK(c.G[i] == doctest::Approx(E[i]).epsilon(1e-9));

    const double V[12] = {-1.69,  6.76,   -6.76,  0.0,     0.0,    0.0,
                          -1.69,  6.084,  -4.394, -1.014,  -0.676, 1.352};
    for (std::size_t i = 0; i < 12; ++i) CHECK(c.U[i] == doctest::Approx(V[i]).epsilon(1e-9));

    CHECK(std::fabs(c.Y[0]) < 1e-15);
    CHECK(std::fabs(c.Y[1]) < 1e-15);
    CHECK(c.Y[2] == doctest::Approx(0.138461538462).epsilon(1e-9));
}

TEST_CASE("a second parameter set matches MATLAB in both forms") {
    const mam::Mamap2mCoefficients<double> a =
        mam::mamap2m_can1_coefficients(0.7, 1.3, 0.3, 0.6);
    CHECK(a.G[0] == doctest::Approx(0.48275862069).epsilon(1e-9));
    CHECK(a.G[8] == doctest::Approx(0.186206896552).epsilon(1e-9));  // corrected G(9)
    CHECK(a.G[14] == doctest::Approx(0.403448275862).epsilon(1e-9));
    CHECK(a.U[1] == doctest::Approx(-0.791816).epsilon(1e-9));
    CHECK(a.Y[0] == doctest::Approx(0.005728352946).epsilon(1e-8));
    CHECK(a.Y[2] == doctest::Approx(-0.131843043995).epsilon(1e-9));

    const mam::Mamap2mCoefficients<double> b =
        mam::mamap2m_can2_coefficients(0.7, 1.3, 0.3, 0.6);
    CHECK(b.G[0] == doctest::Approx(0.21875).epsilon(1e-10));
    CHECK(b.G[13] == doctest::Approx(0.609375).epsilon(1e-10));
    CHECK(b.U[8] == doctest::Approx(-1.6608256).epsilon(1e-9));
    CHECK(b.Y[2] == doctest::Approx(0.0467578125).epsilon(1e-9));
}

TEST_CASE("a third parameter set matches MATLAB in both forms") {
    const mam::Mamap2mCoefficients<double> a =
        mam::mamap2m_can1_coefficients(2.0, 0.5, 0.8, 0.2);
    CHECK(a.G[6] == doctest::Approx(0.426666666667).epsilon(1e-9));
    CHECK(a.G[8] == doctest::Approx(0.0333333333333).epsilon(1e-9));  // corrected G(9)
    CHECK(a.U[4] == doctest::Approx(1.16736).epsilon(1e-9));
    CHECK(a.Y[1] == doctest::Approx(0.0527777777778).epsilon(1e-9));

    const mam::Mamap2mCoefficients<double> b =
        mam::mamap2m_can2_coefficients(2.0, 0.5, 0.8, 0.2);
    CHECK(b.G[9] == doctest::Approx(1.2275862069).epsilon(1e-9));
    CHECK(b.U[2] == doctest::Approx(-2.32).epsilon(1e-9));
    CHECK(b.Y[0] == doctest::Approx(-0.00685554963303).epsilon(1e-8));
}

TEST_CASE("the degenerate denominators are refused by name") {
    // Form 1's denominators are r2(r1-1)+1 and r1 r2 - r2 + 1; both vanish at
    // (r1, r2) = (0, 1).
    CHECK_THROWS_AS(mam::mamap2m_can1_coefficients(1.0, 2.0, 0.0, 1.0), line::NumericError);
    // Form 2 adds r1 + r2 - r1 r2 - 2, which vanishes when (1-r1)(1-r2) = -1;
    // (r1, r2) = (2, 0) kills that and r2(r1-1) - r1 + 2 at once.
    CHECK_THROWS_AS(mam::mamap2m_can2_coefficients(1.0, 2.0, 2.0, 0.0), line::NumericError);
}
