/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Transient CTMC analysis. Oracles: the two-state chain has a closed-form
 * transient solution, pi(t) must converge to the stationary distribution as
 * t grows, and probability must be conserved at every horizon.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/mc/ctmc_solve.h"
#include "line/api/mc/ctmc_uniformization.h"

using line::Matrix;
using line::Real50;
using line::mc::ctmc_makeinfgen;
using line::mc::ctmc_solve;
using line::mc::ctmc_timeaverage;
using line::mc::ctmc_uniformization;

namespace {

/** Two-state chain, rates a (0->1) and b (1->0). */
Matrix<double> two_state(double a, double b) {
    Matrix<double> Q(2, 2, 0.0);
    Q(0, 1) = a;
    Q(1, 0) = b;
    return ctmc_makeinfgen(Q);
}

}  // namespace

TEST_CASE("ctmc_uniformization matches the two-state closed form") {
    const double a = 2.0, b = 3.0, t = 0.35;
    Matrix<double> Q = two_state(a, b);
    std::vector<double> pi0{1.0, 0.0};
    auto r = ctmc_uniformization(pi0, Q, t);

    // pi_0(t) = b/(a+b) + a/(a+b) exp(-(a+b) t) when starting in state 0.
    const double s = a + b;
    const double p0 = b / s + (a / s) * std::exp(-s * t);
    CHECK(r.pi[0] == doctest::Approx(p0).epsilon(1e-11));
    CHECK(r.pi[1] == doctest::Approx(1.0 - p0).epsilon(1e-11));
    CHECK(r.kmax > 0);
}

TEST_CASE("ctmc_uniformization conserves probability and converges to stationarity") {
    Matrix<double> Q(3, 3, 0.0);
    Q(0, 1) = 1.0; Q(1, 2) = 2.0; Q(2, 0) = 3.0; Q(1, 0) = 0.5;
    Q = ctmc_makeinfgen(Q);
    std::vector<double> pi0{1.0, 0.0, 0.0};

    for (double t : {0.1, 1.0, 10.0}) {
        auto r = ctmc_uniformization(pi0, Q, t);
        double s = 0.0;
        for (double v : r.pi) s += v;
        CHECK(s == doctest::Approx(1.0).epsilon(1e-11));
        for (double v : r.pi) CHECK(v >= -1e-12);
    }

    // Long horizon: the transient solution must approach pi Q = 0.
    auto rlong = ctmc_uniformization(pi0, Q, 200.0);
    std::vector<double> pist = ctmc_solve(Q);
    for (std::size_t i = 0; i < 3; ++i)
        CHECK(rlong.pi[i] == doctest::Approx(pist[i]).epsilon(1e-8));
}

TEST_CASE("ctmc_uniformization segments long horizons without losing accuracy") {
    // q t well past the MAXQT = 500 split, so the segmented path is exercised.
    Matrix<double> Q = two_state(2.0, 3.0);
    std::vector<double> pi0{1.0, 0.0};
    auto r = ctmc_uniformization(pi0, Q, 400.0);
    std::vector<double> pist = ctmc_solve(Q);
    CHECK(r.pi[0] == doctest::Approx(pist[0]).epsilon(1e-10));
    CHECK(r.pi[1] == doctest::Approx(pist[1]).epsilon(1e-10));
}

TEST_CASE("ctmc_timeaverage integrates the transient solution") {
    const double a = 2.0, b = 3.0, t = 0.7;
    Matrix<double> Q = two_state(a, b);
    std::vector<double> pi0{1.0, 0.0};
    auto ta = ctmc_timeaverage(pi0, Q, t);

    // (1/t) int_0^t pi_0(u) du = b/(a+b) + a/((a+b)^2 t) (1 - exp(-(a+b) t)).
    const double s = a + b;
    const double avg0 = b / s + (a / (s * s * t)) * (1.0 - std::exp(-s * t));
    CHECK(ta.piTimeAvg[0] == doctest::Approx(avg0).epsilon(1e-10));
    CHECK(ta.piTimeAvg[0] + ta.piTimeAvg[1] == doctest::Approx(1.0).epsilon(1e-10));

    // piExit must agree with a direct transient solve at the same horizon.
    auto tr = ctmc_uniformization(pi0, Q, t);
    CHECK(ta.piExit[0] == doctest::Approx(tr.pi[0]).epsilon(1e-11));
}

TEST_CASE("ctmc_uniformization runs in high precision as well as double") {
    Matrix<Real50> Q(2, 2, Real50(0));
    Q(0, 1) = Real50(2);
    Q(1, 0) = Real50(3);
    Q = ctmc_makeinfgen(Q);
    std::vector<Real50> pi0{Real50(1), Real50(0)};
    auto r = ctmc_uniformization(pi0, Q, Real50(1) / Real50(2));

    Matrix<double> Qd = two_state(2.0, 3.0);
    std::vector<double> pi0d{1.0, 0.0};
    auto rd = ctmc_uniformization(pi0d, Qd, 0.5);
    CHECK(static_cast<double>(r.pi[0]) == doctest::Approx(rd.pi[0]).epsilon(1e-12));
}
