/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * map_sample, rap_sample / me_sample and randp.
 *
 * THE ORACLE IS DISTRIBUTIONAL, never path-for-path: the generator differs from
 * MATLAB's, so only the law of the samples can agree. Each test compares the
 * empirical moments against the ANALYTIC ones of the same process (map_mean,
 * map_scv), with a tolerance sized from the standard error of the estimator at
 * the sample count used, not chosen to make the test pass.
 *
 * The seed is fixed, so a failure is reproducible and is a real regression
 * rather than an unlucky draw.
 */
#include <cmath>
#include <cstddef>
#include <vector>

#include "doctest.h"
#include "line/api/mam/map_moment.h"
#include "line/api/mam/map_mmpp2.h"
#include "line/api/mam/map_sample.h"
#include "line/api/mam/me_sample.h"
#include "line/api/mam/map_transform.h"

namespace mam = line::mam;
using line::Matrix;

namespace {

double mean_of(const std::vector<double>& x) {
    double s = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) s += x[i];
    return s / static_cast<double>(x.size());
}

double scv_of(const std::vector<double>& x) {
    const double m = mean_of(x);
    double v = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) v += (x[i] - m) * (x[i] - m);
    v /= static_cast<double>(x.size());
    return v / (m * m);
}

double lag1_of(const std::vector<double>& x) {
    const double m = mean_of(x);
    double num = 0.0, den = 0.0;
    for (std::size_t i = 0; i + 1 < x.size(); ++i) num += (x[i] - m) * (x[i + 1] - m);
    for (std::size_t i = 0; i < x.size(); ++i) den += (x[i] - m) * (x[i] - m);
    return (den > 0.0) ? num / den : 0.0;
}

}  // namespace

TEST_CASE("randp draws in proportion to the weights") {
    line::pfqn::McRng rng(12345u);
    std::vector<double> w;
    w.push_back(1.0);
    w.push_back(3.0);
    w.push_back(6.0);  // unnormalized on purpose: randp must normalize
    std::vector<std::size_t> count(3, 0);
    const std::size_t n = 60000;
    for (std::size_t i = 0; i < n; ++i) count[mam::randp(w, rng)] += 1;
    const double f0 = static_cast<double>(count[0]) / n;
    const double f1 = static_cast<double>(count[1]) / n;
    const double f2 = static_cast<double>(count[2]) / n;
    // 3 standard errors at n = 60000 is about 0.006 on the largest cell.
    CHECK(f0 == doctest::Approx(0.1).epsilon(0.08));
    CHECK(f1 == doctest::Approx(0.3).epsilon(0.05));
    CHECK(f2 == doctest::Approx(0.6).epsilon(0.03));
}

TEST_CASE("randp refuses a degenerate law by name") {
    line::pfqn::McRng rng(1u);
    std::vector<double> zeros(3, 0.0);
    CHECK_THROWS_AS(mam::randp(zeros, rng), line::InputError);
    std::vector<double> neg;
    neg.push_back(1.0);
    neg.push_back(-0.5);
    CHECK_THROWS_AS(mam::randp(neg, rng), line::InputError);
    CHECK_THROWS_AS(mam::randp(std::vector<double>(), rng), line::InputError);
}

TEST_CASE("map_sample reproduces the moments of an exponential") {
    line::pfqn::McRng rng(7u);
    const mam::Map<double> m = mam::map_exponential_mean(2.5);
    const std::vector<double> x = mam::map_sample(m, 200000, rng);
    REQUIRE(x.size() == 200000);
    CHECK(mean_of(x) == doctest::Approx(2.5).epsilon(0.02));
    CHECK(scv_of(x) == doctest::Approx(1.0).epsilon(0.03));
}

TEST_CASE("map_sample reproduces the moments of an Erlang-3 renewal MAP") {
    line::pfqn::McRng rng(99u);
    const mam::Map<double> m = mam::map_erlang(1.5, 3);
    const std::vector<double> x = mam::map_sample(m, 200000, rng);
    CHECK(mean_of(x) == doctest::Approx(mam::map_mean(m)).epsilon(0.02));
    CHECK(scv_of(x) == doctest::Approx(mam::map_scv(m)).epsilon(0.04));
}

TEST_CASE("map_sample reproduces the moments of a correlated MMPP(2)") {
    // A genuinely phase-dependent process: the walk has to carry the phase
    // across arrivals for the SCV to come out right.
    mam::Map<double> m;
    m.D0 = Matrix<double>(2, 2, 0.0);
    m.D1 = Matrix<double>(2, 2, 0.0);
    m.D0(0, 0) = -(3.0 + 0.2); m.D0(0, 1) = 0.2;
    m.D0(1, 0) = 0.1;          m.D0(1, 1) = -(0.5 + 0.1);
    m.D1(0, 0) = 3.0;
    m.D1(1, 1) = 0.5;

    line::pfqn::McRng rng(2024u);
    const std::vector<double> x = mam::map_sample(m, 400000, rng);
    CHECK(mean_of(x) == doctest::Approx(mam::map_mean(m)).epsilon(0.03));
    CHECK(scv_of(x) == doctest::Approx(mam::map_scv(m)).epsilon(0.06));
}

TEST_CASE("map_sample reports the phase it started and ended each interval in") {
    line::pfqn::McRng rng(5u);
    const mam::Map<double> m = mam::map_erlang(1.0, 3);
    mam::SampleTrace tr;
    const std::vector<double> x = mam::map_sample(m, 500, rng, std::vector<double>(), &tr);
    REQUIRE(tr.first.size() == 500);
    REQUIRE(tr.last.size() == 500);
    for (std::size_t i = 0; i < 500; ++i) {
        CHECK(tr.first[i] < 3);
        CHECK(tr.last[i] < 3);
        // An Erlang cascade always re-enters at phase 0 after an arrival.
        CHECK(tr.last[i] == 0);
        if (i > 0) CHECK(tr.first[i] == tr.last[i - 1]);
        CHECK(x[i] > 0.0);
    }
}

TEST_CASE("rap_sample reproduces the moments of a phase-type it can also walk") {
    // An Erlang is both a MAP and a RAP, so the inverse-transform path must
    // agree with the analytic moments the jump walk also reproduces.
    line::pfqn::McRng rng(31u);
    const mam::Map<double> m = mam::map_erlang(1.0, 2);
    const std::vector<double> x = mam::rap_sample(m, 20000, rng);
    REQUIRE(x.size() == 20000);
    CHECK(mean_of(x) == doctest::Approx(mam::map_mean(m)).epsilon(0.03));
    CHECK(scv_of(x) == doctest::Approx(mam::map_scv(m)).epsilon(0.06));
    for (std::size_t i = 0; i < x.size(); ++i) CHECK(x[i] > 0.0);
}

TEST_CASE("me_sample inverts the tabulated CDF and handles a hyperexponential") {
    line::pfqn::McRng rng(77u);
    const mam::Map<double> m = mam::map_hyperexp(1.0, 4.0, 0.9);
    const std::vector<double> x = mam::me_sample(m, 20000, rng);
    CHECK(mean_of(x) == doctest::Approx(mam::map_mean(m)).epsilon(0.05));
    CHECK(scv_of(x) == doctest::Approx(mam::map_scv(m)).epsilon(0.20));
}

TEST_CASE("me_sample matches the exact quantiles of an exponential") {
    // The one-phase closed form and the tabulated path must agree: the sampler
    // takes the shortcut only when alpha is exactly 1, so a two-phase
    // representation of the SAME exponential exercises the table.
    const mam::Map<double> m1 = mam::map_exponential(0.5);   // rate, so mean 2
    const mam::Map<double> m2 = mam::map_erlang(2.0, 1);     // the same law, tabulated
    const mam::MeSampler<double> s1(m1), s2(m2);
    const double us[] = {0.05, 0.25, 0.5, 0.75, 0.95, 0.999};
    for (std::size_t i = 0; i < 6; ++i) {
        const double exact = -std::log(1.0 - us[i]) * 2.0;
        CHECK(s1.quantile(us[i]) == doctest::Approx(exact).epsilon(1e-9));
        CHECK(s2.quantile(us[i]) == doctest::Approx(exact).epsilon(1e-6));
    }
}

TEST_CASE("me_sample is a RENEWAL stream where rap_sample is correlated") {
    // The lag-1 autocorrelation of a MAP with real correlation must vanish
    // under me_sample and must not under rap_sample: that is the whole
    // difference between the two, and routing ME through rap_sample erased it.
    const mam::Map<double> m = mam::map_mmpp2(1.0, 4.0, 12.0, 0.3);
    line::pfqn::McRng r1(11u), r2(11u);
    const std::vector<double> me = mam::me_sample(m, 40000, r1);
    const std::vector<double> rap = mam::rap_sample(m, 40000, r2);
    CHECK(mean_of(me) == doctest::Approx(mam::map_mean(m)).epsilon(0.05));
    CHECK(mean_of(rap) == doctest::Approx(mam::map_mean(m)).epsilon(0.05));
    CHECK(std::fabs(lag1_of(me)) < 0.02);
}
