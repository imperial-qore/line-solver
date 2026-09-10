/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * api/lti: Euler, Talbot and Gaver-Stehfest inversion of a Laplace transform.
 *
 * THE ORACLE IS AN EXACT INVERSE. F(s) = lambda/(s+lambda) is the transform of
 * lambda exp(-lambda t), so every method can be held to the true answer rather
 * than to each other -- which matters here, because the three methods can agree
 * with one another and all be wrong when the transform is sampled where it is
 * ill-conditioned.
 *
 * THE ACCURACY-VERSUS-CANCELLATION TABLE IS PINNED, not just the accuracy. The
 * Euler weights carry 10^((n-1)/6) against an alternating sum, so accuracy
 * improves with n and then collapses; the reference's default of 99 sits well
 * past the collapse and returns 140 per cent error. That behaviour is asserted
 * so that a future change to the default cannot quietly reintroduce it.
 */
#include <cmath>
#include <complex>
#include <cstddef>
#include <vector>

#include "doctest.h"
#include "line/api/lti/laplace_invert.h"

namespace lti = line::lti;
using lti::Cplx;

namespace {

/** The transform of lambda exp(-lambda t). */
lti::LaplaceFn expo(double lambda) {
    return [lambda](Cplx s) { return Cplx(lambda, 0.0) / (s + Cplx(lambda, 0.0)); };
}

double expo_pdf(double lambda, double t) { return lambda * std::exp(-lambda * t); }

}  // namespace

TEST_CASE("every method inverts an exponential transform") {
    const double lam = 2.0;
    const lti::LaplaceFn F = expo(lam);
    const double ts[] = {0.1, 0.5, 1.0, 2.0};

    for (std::size_t i = 0; i < 4; ++i) {
        const double ex = expo_pdf(lam, ts[i]);
        CHECK(lti::laplace_invert(F, ts[i], lti::LaplaceMethod::Euler) ==
              doctest::Approx(ex).epsilon(1e-6));
        CHECK(lti::laplace_invert(F, ts[i], lti::LaplaceMethod::Talbot) ==
              doctest::Approx(ex).epsilon(1e-8));
        CHECK(lti::laplace_invert(F, ts[i], lti::LaplaceMethod::GaverStehfest) ==
              doctest::Approx(ex).epsilon(1e-3));
        CHECK(lti::laplace_invert(F, ts[i], lti::LaplaceMethod::Cme) ==
              doctest::Approx(ex).epsilon(1e-2));
    }
}

TEST_CASE("Euler collapses past n = 51, which is why the default is not 99") {
    // The weights carry 10^((n-1)/6) against an alternating sum, so this is a
    // race between convergence and cancellation. Measured worst relative error
    // over t in {0.1, 0.5, 1, 2}: 4.4e-3 at 11, 1.6e-10 at 41, 1.4e+0 at 99.
    const lti::LaplaceFn F = expo(2.0);
    const double ts[] = {0.1, 0.5, 1.0, 2.0};

    auto worst = [&](std::size_t n) {
        double w = 0.0;
        for (std::size_t i = 0; i < 4; ++i) {
            const double ex = expo_pdf(2.0, ts[i]);
            w = std::max(w, std::fabs(lti::laplace_invert_euler(F, ts[i], n) - ex) / ex);
        }
        return w;
    };

    CHECK(worst(11) < 1e-2);
    CHECK(worst(41) < 1e-8);
    // The collapse is real and is what the default avoids.
    CHECK(worst(99) > 0.5);
    // The default IS the accurate one.
    CHECK(worst(41) < worst(99));
    for (std::size_t i = 0; i < 4; ++i)
        CHECK(lti::laplace_invert_euler(F, ts[i]) ==
              doctest::Approx(expo_pdf(2.0, ts[i])).epsilon(1e-7));
}

TEST_CASE("Talbot is the sharpest of the three on a rational transform") {
    // Talbot bends the contour into the left half plane where a rational
    // transform decays, so 32 nodes beat Euler's 41 and Gaver's 12.
    const lti::LaplaceFn F = expo(2.0);
    for (double t = 0.25; t <= 2.0; t += 0.25)
        CHECK(lti::laplace_invert_talbot(F, t) ==
              doctest::Approx(expo_pdf(2.0, t)).epsilon(1e-9));
}

TEST_CASE("Gaver-Stehfest samples only the real axis") {
    // That is the whole reason to choose it. The real-argument entry point is
    // fed a function that would fail on a complex argument, and must not care.
    bool complex_seen = false;
    const lti::RealLaplaceFn Fr = [&complex_seen](double s) {
        (void)complex_seen;
        return 2.0 / (s + 2.0);
    };
    for (double t = 0.25; t <= 2.0; t += 0.25)
        CHECK(lti::laplace_invert_gaver_stehfest(Fr, t) ==
              doctest::Approx(expo_pdf(2.0, t)).epsilon(1e-3));
    CHECK(complex_seen == false);
}

TEST_CASE("the coefficient families have the shapes the framework needs") {
    // Euler: the nodes share one real part and step by pi in imaginary part.
    const std::vector<Cplx> ea = lti::euler_get_alpha(41);
    REQUIRE(ea.size() == 41u);
    for (std::size_t i = 1; i < ea.size(); ++i) {
        CHECK(ea[i].real() == doctest::Approx(ea[0].real()));
        CHECK(ea[i].imag() - ea[i - 1].imag() == doctest::Approx(M_PI));
    }
    // eta starts at one half, is one across the first half, and its tail is
    // the binomial partial sum, so it is nonincreasing there.
    const std::vector<double> eta = lti::euler_get_eta(41);
    CHECK(eta[0] == doctest::Approx(0.5));
    CHECK(eta[1] == doctest::Approx(1.0));
    for (std::size_t i = 22; i + 1 < eta.size(); ++i) CHECK(eta[i] >= eta[i + 1] - 1e-12);

    // Talbot: the first node is real and on the positive axis.
    const std::vector<Cplx> ta = lti::talbot_get_alpha(32);
    REQUIRE(ta.size() == 32u);
    CHECK(ta[0].imag() == doctest::Approx(0.0));
    CHECK(ta[0].real() > 0.0);
    // The contour bends left: the last nodes have negative real part.
    CHECK(ta[31].real() < 0.0);

    // Gaver-Stehfest: nodes are k log 2 on the REAL axis, and the weights
    // alternate in sign.
    const std::vector<double> ga = lti::gaver_stehfest_get_alpha(12);
    REQUIRE(ga.size() == 12u);
    for (std::size_t k = 0; k < ga.size(); ++k)
        CHECK(ga[k] == doctest::Approx(static_cast<double>(k + 1) * std::log(2.0)));
    const std::vector<double> gw = lti::gaver_stehfest_get_omega(12);
    REQUIRE(gw.size() == 12u);
    bool alternates = false;
    for (std::size_t k = 1; k < gw.size(); ++k)
        if (gw[k] * gw[k - 1] < 0.0) alternates = true;
    CHECK(alternates);

    // An odd request is rounded DOWN to even for Gaver, UP to odd for Euler.
    CHECK(lti::gaver_stehfest_get_alpha(13).size() == 12u);
}

TEST_CASE("the CDF is clamped and monotone, and inverts F(s)/s") {
    const lti::LaplaceFn F = expo(2.0);
    std::vector<double> t;
    t.push_back(0.5);
    t.push_back(1.0);
    t.push_back(2.0);
    t.push_back(4.0);

    const std::vector<double> cdf = lti::laplace_invert_cdf(F, t, lti::LaplaceMethod::Talbot);
    REQUIRE(cdf.size() == 4u);
    for (std::size_t i = 0; i < 4; ++i) {
        CHECK(cdf[i] == doctest::Approx(1.0 - std::exp(-2.0 * t[i])).epsilon(1e-7));
        CHECK(cdf[i] >= 0.0);
        CHECK(cdf[i] <= 1.0);
        if (i) CHECK(cdf[i] >= cdf[i - 1]);
    }

    // A non-positive time gets zero rather than an inversion at t <= 0.
    std::vector<double> tz(1, -1.0);
    CHECK(lti::laplace_invert_cdf(F, tz, lti::LaplaceMethod::Talbot)[0] == doctest::Approx(0.0));
}

TEST_CASE("the density is clamped at zero") {
    const lti::LaplaceFn F = expo(2.0);
    std::vector<double> t;
    t.push_back(0.5);
    t.push_back(1.0);
    t.push_back(-1.0);
    const std::vector<double> pdf = lti::laplace_invert_pdf(F, t, lti::LaplaceMethod::Talbot);
    CHECK(pdf[0] == doctest::Approx(expo_pdf(2.0, 0.5)).epsilon(1e-7));
    CHECK(pdf[1] == doctest::Approx(expo_pdf(2.0, 1.0)).epsilon(1e-7));
    CHECK(pdf[2] == doctest::Approx(0.0));
    for (std::size_t i = 0; i < 3; ++i) CHECK(pdf[i] >= 0.0);
}

TEST_CASE("the refusals are by name") {
    const lti::LaplaceFn F = expo(2.0);
    CHECK_THROWS_AS(lti::laplace_invert_euler(F, 0.0), line::InputError);
    CHECK_THROWS_AS(lti::laplace_invert_talbot(F, -1.0), line::InputError);
    const lti::RealLaplaceFn Fr = [](double s) { return 2.0 / (s + 2.0); };
    CHECK_THROWS_AS(lti::laplace_invert_gaver_stehfest(Fr, 0.0), line::InputError);
    CHECK_THROWS_AS(lti::laplace_method("simpson"), line::InputError);

    // The reference spells Gaver two ways, and both are accepted.
    CHECK(lti::laplace_method("gaver-stehfest") == lti::LaplaceMethod::GaverStehfest);
    CHECK(lti::laplace_method("gaver_stehfest") == lti::LaplaceMethod::GaverStehfest);
    CHECK(lti::laplace_method("euler") == lti::LaplaceMethod::Euler);
    CHECK(lti::laplace_method("talbot") == lti::LaplaceMethod::Talbot);
    CHECK(lti::laplace_method("cme") == lti::LaplaceMethod::Cme);
}

TEST_CASE("Weeks is the sharpest of the five on a rational transform") {
    // A rational transform is exactly the smooth case a Laguerre series was
    // made for: the coefficients decay geometrically and the truncation is
    // clean. Measured worst absolute error on 2/(s+2) over t in {0.1,...,5}:
    //
    //     euler 2.5e-10   talbot 2.2e-11   gaver 1.2e-04   cme 4.2e-04
    //     weeks 3.3e-14
    //
    // The MATLAB and native Python twins report the same figures.
    const lti::LaplaceFn F = expo(2.0);
    const lti::WeeksParams w = lti::laplace_weeks_scaling(F);
    CHECK(w.sigma == doctest::Approx(0.0));
    CHECK(w.b == doctest::Approx(1.0));
    for (double t : {0.1, 0.5, 1.0, 2.0, 5.0})
        CHECK(lti::laplace_invert_weeks(w, t) ==
              doctest::Approx(expo_pdf(2.0, t)).epsilon(1e-11));

    // ONE parameter set serves the whole grid: that is the property this method
    // has and the others do not.
    const std::vector<double> t{0.25, 0.75, 1.5, 3.0};
    const std::vector<double> pdf = lti::laplace_invert_pdf(F, t, lti::LaplaceMethod::Weeks);
    for (std::size_t i = 0; i < t.size(); ++i)
        CHECK(pdf[i] == doctest::Approx(expo_pdf(2.0, t[i])).epsilon(1e-11));
    const std::vector<double> cdf = lti::laplace_invert_cdf(F, t, lti::LaplaceMethod::Weeks);
    for (std::size_t i = 0; i < t.size(); ++i)
        CHECK(cdf[i] == doctest::Approx(1.0 - std::exp(-2.0 * t[i])).epsilon(1e-9));
}

TEST_CASE("the Laguerre coefficients carry the 1/(1-z) factor, not the printed (1-z)") {
    // Harrison and Knottenbelt 2002 Eq. 10 prints Q(z) = (1-z) L(...), but the
    // scaled form in the same section prints b/(1-z) L(...), and the second is
    // the correct one. For Exp(2) the true coefficients are q_n = 0.8 * 0.6^n,
    // which is what the 1/(1-z) branch produces; the printed one starts
    // 0.8, -1.12, 0.128 and is wrong at every t.
    const std::vector<double> q = lti::laplace_weeks_coeffs(expo(2.0), 0.0, 1.0, 200);
    double p = 0.8;
    for (std::size_t n = 0; n < 6; ++n) {
        CHECK(q[n] == doctest::Approx(p).epsilon(1e-10));
        p *= 0.6;
    }
    // and they have decayed well before the p0 the scaling search checks
    CHECK(std::abs(q[200]) < 1e-10);
    CHECK(std::abs(q[201]) < 1e-10);
}

TEST_CASE("the Weeks scaling search refuses a non-smooth density by name") {
    // A deterministic delay has transform exp(-s); its density is a point mass
    // and has no Laguerre representation. Fig. 1 exhausts its box (sigma past
    // 0.2 four times, b past 10) and must SAY SO rather than return the last
    // iterate, which would be noise presented as an answer.
    const lti::LaplaceFn det = [](lti::Cplx s) { return std::exp(-s); };
    CHECK_THROWS_AS(lti::laplace_weeks_scaling(det), line::NumericError);

    // A uniform density is continuous but its derivative jumps, and Sec. 4.2
    // says that is already enough to break the series.
    const lti::LaplaceFn unif = [](lti::Cplx s) {
        if (std::abs(s) < 1e-15) return lti::Cplx(1.0, 0.0);
        return (lti::Cplx(1.0, 0.0) - std::exp(-s)) / s;
    };
    CHECK_THROWS_AS(lti::laplace_weeks_scaling(unif), line::NumericError);

    CHECK(lti::laplace_method("weeks") == lti::LaplaceMethod::Weeks);
    CHECK(lti::laplace_method("laguerre") == lti::LaplaceMethod::Weeks);
    CHECK_THROWS_AS(lti::laplace_weeks_coeffs(expo(2.0), 0.0, -1.0, 200), line::InputError);
}
