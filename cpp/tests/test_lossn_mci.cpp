/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Monte Carlo summation for loss networks. A random estimator cannot be pinned
 * to a recorded number, so nothing here asserts one. The oracles are:
 *   1. A CLOSED FORM. One link of capacity C with unit circuit requirements is
 *      an Erlang loss system: the blocking probability is Erlang B and the
 *      normalizing constant is sum_n nu^n/n!, both known exactly. The estimate
 *      must approach them.
 *   2. CONVERGENCE, not a single draw: the root-mean-square error over several
 *      independent seeds must fall when the sample count rises. A single lucky
 *      sample would satisfy an absolute tolerance and tell us nothing.
 *   3. COVERAGE: the true value must lie inside the reported confidence
 *      interval.
 *   4. REPRODUCIBILITY: one seed, one answer, bit for bit, and a different
 *      seed a different one. There is no global generator to reseed.
 *   5. A second, genuinely multi-link instance whose exact g is enumerated in
 *      the test itself, so the estimator is checked away from the degenerate
 *      case where every sampled state is feasible.
 */
#include <cmath>
#include <cstdint>
#include <random>
#include <vector>

#include "doctest.h"
#include "line/api/lossn/lossn_mci.h"

using line::Matrix;
using line::Real50;
using line::lossn::LossnMciOptions;
using line::lossn::LossnMciResult;
using line::lossn::lossn_mci;

namespace {

/** Erlang B by its rational recursion, exact for integer capacity. */
double erlangB(double nu, int C) {
    double b = 1.0;
    for (int k = 1; k <= C; ++k) b = nu * b / (k + nu * b);
    return b;
}

/** Single link, capacity C, one route needing one circuit. */
LossnMciResult<double> singleLink(double nu, double C, std::size_t S, std::uint64_t seed) {
    LossnMciOptions<double> opt;
    opt.samples = S;
    return lossn_mci<double>({nu}, Matrix<double>{{1.0}}, {C}, opt, seed);
}

}  // namespace

TEST_CASE("lossn_mci converges to the Erlang B blocking probability") {
    const double nu = 4.0;
    const int C = 5;
    const double exact = erlangB(nu, C);
    CHECK(exact == doctest::Approx(0.199066874027994).epsilon(1e-12));  // MATLAB

    double rmsSmall = 0.0, rmsLarge = 0.0;
    const int seeds = 6;
    for (int s = 0; s < seeds; ++s) {
        const double eSmall = singleLink(nu, C, 2000, 100 + s).Loss[0] - exact;
        const double eLarge = singleLink(nu, C, 50000, 100 + s).Loss[0] - exact;
        rmsSmall += eSmall * eSmall;
        rmsLarge += eLarge * eLarge;
    }
    rmsSmall = std::sqrt(rmsSmall / seeds);
    rmsLarge = std::sqrt(rmsLarge / seeds);
    INFO("rms at S=2000 ", rmsSmall, ", at S=50000 ", rmsLarge);
    CHECK(rmsLarge < rmsSmall);   // twenty-five times the samples, less error
    CHECK(rmsLarge < 5e-3);
    CHECK(rmsSmall < 5e-2);
}

TEST_CASE("lossn_mci recovers the exact normalizing constant of the Erlang system") {
    const double nu = 4.0;
    const int C = 5;
    double g = 0.0, fact = 1.0;
    for (int n = 0; n <= C; ++n) {
        if (n > 0) fact *= n;
        g += std::pow(nu, n) / fact;
    }
    CHECK(std::log(g) == doctest::Approx(3.75809452313541).epsilon(1e-12));  // MATLAB
    const LossnMciResult<double> r = singleLink(nu, C, 50000, 2026);
    CHECK(r.lG == doctest::Approx(std::log(g)).epsilon(1e-2));
    CHECK(r.nsamples == 50000);
    CHECK(r.level == doctest::Approx(0.95));
    // carried load and blocking are complementary by construction
    CHECK(r.QLen[0] == doctest::Approx(nu * (1.0 - r.Loss[0])).epsilon(1e-12));
    CHECK(r.Loss[0] == doctest::Approx(r.lossPoint[0]));
    CHECK(r.acceptPoint[0] == doctest::Approx(1.0 - r.Loss[0]));
}

TEST_CASE("lossn_mci reports a confidence interval that covers the truth") {
    const double exact = erlangB(4.0, 5);
    int covered = 0;
    const int seeds = 8;
    for (int s = 0; s < seeds; ++s) {
        const LossnMciResult<double> r = singleLink(4.0, 5.0, 20000, 7000 + s);
        CHECK(r.lossCI(0, 0) <= r.Loss[0]);
        CHECK(r.Loss[0] <= r.lossCI(0, 1));
        if (r.lossCI(0, 0) <= exact && exact <= r.lossCI(0, 1)) ++covered;
    }
    INFO("covered ", covered, " of ", seeds);
    CHECK(covered >= 6);  // nominal 95 per cent, allowing for the small sample
}

TEST_CASE("lossn_mci is reproducible in its seed and carries no global state") {
    const LossnMciResult<double> a = singleLink(4.0, 5.0, 5000, 42);
    const LossnMciResult<double> b = singleLink(4.0, 5.0, 5000, 42);
    CHECK(a.Loss[0] == b.Loss[0]);  // bit for bit, not to a tolerance
    CHECK(a.lG == b.lG);
    const LossnMciResult<double> c = singleLink(4.0, 5.0, 5000, 43);
    CHECK(a.Loss[0] != c.Loss[0]);

    // the caller's own engine is advanced in place, so consecutive calls on
    // one engine differ while a freshly seeded engine repeats the first
    std::mt19937_64 gen(42);
    LossnMciOptions<double> opt;
    opt.samples = 5000;
    const LossnMciResult<double> d =
        lossn_mci<double>({4.0}, Matrix<double>{{1.0}}, {5.0}, opt, gen);
    const LossnMciResult<double> e =
        lossn_mci<double>({4.0}, Matrix<double>{{1.0}}, {5.0}, opt, gen);
    CHECK(d.Loss[0] == a.Loss[0]);
    CHECK(d.Loss[0] != e.Loss[0]);
}

TEST_CASE("lossn_mci on a two-link two-route network approaches the enumerated truth") {
    // Route 1 uses link 1 only, route 2 uses both. Feasible states are
    // n1 + n2 <= 6 and n2 <= 4, small enough to enumerate exactly.
    const double nu[2] = {3.0, 2.0};
    const Matrix<double> A{{1.0, 1.0}, {0.0, 1.0}};
    const std::vector<double> C{6.0, 4.0};
    double g = 0.0, g1 = 0.0, g2 = 0.0;
    std::vector<double> f(8, 1.0);
    for (std::size_t n = 1; n < f.size(); ++n) f[n] = f[n - 1] * static_cast<double>(n);
    for (int n1 = 0; n1 <= 6; ++n1)
        for (int n2 = 0; n2 <= 4; ++n2) {
            if (n1 + n2 > 6) continue;
            const double t = std::pow(nu[0], n1) / f[n1] * std::pow(nu[1], n2) / f[n2];
            g += t;
            if (n1 + n2 <= 5) g1 += t;
            if (n1 + n2 <= 5 && n2 <= 3) g2 += t;
        }
    // MATLAB, by the same enumeration
    CHECK(std::log(g) == doctest::Approx(4.71816399380131).epsilon(1e-12));
    CHECK(1.0 - g1 / g == doctest::Approx(0.185888132187116).epsilon(1e-12));
    CHECK(1.0 - g2 / g == doctest::Approx(0.209705630605486).epsilon(1e-12));

    LossnMciOptions<double> opt;
    opt.samples = 100000;
    const LossnMciResult<double> r =
        lossn_mci<double>({nu[0], nu[1]}, A, C, opt, static_cast<std::uint64_t>(11));
    CHECK(r.lG == doctest::Approx(std::log(g)).epsilon(2e-2));
    CHECK(r.Loss[0] == doctest::Approx(1.0 - g1 / g).epsilon(5e-2));
    CHECK(r.Loss[1] == doctest::Approx(1.0 - g2 / g).epsilon(5e-2));
    // route 2 crosses both links, so it must block at least as often
    CHECK(r.Loss[1] > r.Loss[0]);
}

TEST_CASE("lossn_mci accepts explicit importance parameters and rejects bad input") {
    LossnMciOptions<double> opt;
    opt.samples = 20000;
    opt.gamma = std::vector<double>{4.0};  // the offered load itself
    const LossnMciResult<double> r =
        lossn_mci<double>({4.0}, Matrix<double>{{1.0}}, {5.0}, opt, static_cast<std::uint64_t>(5));
    CHECK(r.Loss[0] == doctest::Approx(erlangB(4.0, 5)).epsilon(5e-2));

    LossnMciOptions<double> bad;
    bad.samples = 1;
    CHECK_THROWS_AS(lossn_mci<double>({4.0}, Matrix<double>{{1.0}}, {5.0}, bad,
                                      static_cast<std::uint64_t>(1)),
                    line::InputError);
    LossnMciOptions<double> ok;
    ok.samples = 100;
    CHECK_THROWS_AS(lossn_mci<double>({4.0, 1.0}, Matrix<double>{{1.0}}, {5.0}, ok,
                                      static_cast<std::uint64_t>(1)),
                    line::InputError);
    ok.gamma = std::vector<double>{1.0, 2.0};
    CHECK_THROWS_AS(
        lossn_mci<double>({4.0}, Matrix<double>{{1.0}}, {5.0}, ok, static_cast<std::uint64_t>(1)),
        line::InputError);
}

TEST_CASE("lossn_mci runs at high precision and gives the same estimate") {
    // The estimate is a function of the same uniforms, so raising the working
    // precision must not move it beyond the double rounding of the weights.
    LossnMciOptions<Real50> opt;
    opt.samples = 5000;
    const LossnMciResult<Real50> r =
        lossn_mci<Real50>({Real50(4)}, Matrix<Real50>{{Real50(1)}}, {Real50(5)}, opt,
                          static_cast<std::uint64_t>(42));
    const LossnMciResult<double> d = singleLink(4.0, 5.0, 5000, 42);
    CHECK(static_cast<double>(r.Loss[0]) == doctest::Approx(d.Loss[0]).epsilon(1e-12));
}
