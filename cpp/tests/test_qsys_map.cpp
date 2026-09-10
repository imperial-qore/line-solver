/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Queues driven by a non-Poisson arrival process: PH/M/1, PH/M/c and D/M/c.
 *
 * The governing check is the Markovian collapse. A phase-type arrival process
 * whose Laplace transform is lambda/(s+lambda) is a Poisson process, so
 * PH/M/1 must return the M/M/1 values and PH/M/c the M/M/c ones, to the
 * accuracy of the root finder and of the R iteration respectively. Beyond
 * that the oracles are Little's law, W = Wq + 1/mu, and MATLAB reference
 * values obtained by running matlab/src/api/qsys through
 *   matlab -singleCompThread -batch "addpath(genpath('matlab/src')); ..."
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/qsys/qsys_dmc.h"
#include "line/api/qsys/qsys_gm1.h"
#include "line/api/qsys/qsys_mm1.h"
#include "line/api/qsys/qsys_mmk.h"
#include "line/api/qsys/qsys_phm1.h"
#include "line/api/qsys/qsys_phmc.h"

using line::Matrix;

namespace {
/** One-phase PH with rate lambda, i.e. a Poisson arrival process. */
Matrix<double> exp_ph(double lambda) {
    Matrix<double> T(1, 1);
    T(0, 0) = -lambda;
    return T;
}
/** Erlang-2 with the given phase rate: mean 2/rate, scv 1/2. */
Matrix<double> erlang2(double rate) {
    Matrix<double> T(2, 2, 0.0);
    T(0, 0) = -rate;
    T(0, 1) = rate;
    T(1, 1) = -rate;
    return T;
}
}  // namespace

TEST_CASE("PH/M/1 with a Poisson arrival process is M/M/1") {
    // lambda = 2, mu = 3: sigma must be rho = 2/3 and L = rho/(1-rho) = 2.
    auto r = line::qsys::qsys_phm1(std::vector<double>{1.0}, exp_ph(2.0), 3.0);
    auto mm1 = line::qsys::qsys_mm1(2.0, 3.0);
    CHECK(r.sigma == doctest::Approx(2.0 / 3.0).epsilon(1e-13));
    CHECK(r.utilization == doctest::Approx(2.0 / 3.0).epsilon(1e-15));
    CHECK(r.meanQueueLength == doctest::Approx(2.0).epsilon(1e-12));
    CHECK(r.meanSojournTime == doctest::Approx(mm1.W).epsilon(1e-12));
    // MATLAB qsys_phm1(1,-2,3) -> L 1.9999999999999991, W 0.99999999999999956,
    // sigma 0.66666666666666652; the bisection is asserted at 1e-12 relative.
    CHECK(r.meanQueueLength == doctest::Approx(1.9999999999999991).epsilon(1e-12));
    CHECK(r.meanSojournTime == doctest::Approx(0.99999999999999956).epsilon(1e-12));
}

TEST_CASE("PH/M/1 with Erlang-2 arrivals matches MATLAB and the exact root") {
    // Erlang-2 rate 4 (mean 1/2, lambda = 2), mu = 4, rho = 1/2. The root of
    // sigma = (4/(4+4(1-sigma)))^2 is sigma = (3-sqrt(5))/2.
    std::vector<double> alpha = {1.0, 0.0};
    auto r = line::qsys::qsys_phm1(alpha, erlang2(4.0), 4.0);
    const double exact_sigma = (3.0 - std::sqrt(5.0)) / 2.0;
    CHECK(r.sigma == doctest::Approx(exact_sigma).epsilon(1e-13));
    CHECK(r.utilization == doctest::Approx(0.5).epsilon(1e-15));
    // MATLAB references, bisection vs fzero: 1e-12 relative.
    CHECK(r.meanQueueLength == doctest::Approx(0.80901699437494745).epsilon(1e-12));
    CHECK(r.meanWaitingQueue == doctest::Approx(0.30901699437494751).epsilon(1e-12));
    CHECK(r.meanWaitingTime == doctest::Approx(0.15450849718747375).epsilon(1e-12));
    CHECK(r.meanSojournTime == doctest::Approx(0.40450849718747373).epsilon(1e-12));
    // W = Wq + 1/mu, and Lq = lambda Wq with lambda = 2.
    CHECK(r.meanSojournTime == doctest::Approx(r.meanWaitingTime + 0.25).epsilon(1e-14));
    CHECK(r.meanWaitingQueue == doctest::Approx(2.0 * r.meanWaitingTime).epsilon(1e-14));
    // Less arrival variability than Poisson, so less delay than M/M/1.
    CHECK(r.meanSojournTime < line::qsys::qsys_mm1(2.0, 4.0).W);
}

TEST_CASE("PH/M/1 agrees with the G/M/1 closed form at its own root") {
    std::vector<double> alpha = {1.0, 0.0};
    auto r = line::qsys::qsys_phm1(alpha, erlang2(4.0), 4.0);
    CHECK(r.meanSojournTime == doctest::Approx(line::qsys::qsys_gm1(r.sigma, 4.0)).epsilon(1e-12));
}

TEST_CASE("PH/M/1 rejects an unstable or malformed instance") {
    CHECK_THROWS_AS(line::qsys::qsys_phm1(std::vector<double>{1.0}, exp_ph(2.0), 1.0),
                    line::InputError);  // rho = 2
    CHECK_THROWS_AS(line::qsys::qsys_phm1(std::vector<double>{1.0, 0.0}, exp_ph(2.0), 3.0),
                    line::InputError);  // alpha length
    CHECK_THROWS_AS(line::qsys::qsys_phm1(std::vector<double>{1.0}, exp_ph(2.0), -1.0),
                    line::InputError);  // mu <= 0
}

TEST_CASE("PH/M/c with a Poisson arrival process is M/M/c") {
    // c = 1: M/M/1 with lambda = 2, mu = 3.
    auto r1 = line::qsys::qsys_phmc(std::vector<double>{1.0}, exp_ph(2.0), 3.0, 1u);
    CHECK(r1.meanQueueLength == doctest::Approx(2.0).epsilon(1e-9));
    CHECK(r1.meanSojournTime == doctest::Approx(1.0).epsilon(1e-9));
    // c = 2: M/M/2 with lambda = 1, mu = 1, so L = 4/3 and W = 4/3.
    auto r2 = line::qsys::qsys_phmc(std::vector<double>{1.0}, exp_ph(1.0), 1.0, 2u);
    auto mmk = line::qsys::qsys_mmk(1.0, 1.0, 2u);
    CHECK(r2.meanQueueLength == doctest::Approx(4.0 / 3.0).epsilon(1e-9));
    CHECK(r2.meanSojournTime == doctest::Approx(mmk.W).epsilon(1e-9));
    CHECK(r2.utilization == doctest::Approx(0.5).epsilon(1e-15));
    // MATLAB qsys_phmc(1,-1,1,2) -> L 1.3333333333333162, Lq 0.33333333333332071.
    CHECK(r2.meanQueueLength == doctest::Approx(1.3333333333333162).epsilon(1e-9));
    CHECK(r2.meanWaitingQueue == doctest::Approx(0.33333333333332071).epsilon(1e-8));
}

TEST_CASE("PH/M/c with Erlang-2 arrivals matches MATLAB") {
    // Erlang-2 rate 4 (lambda = 2), mu = 1.5, c = 2, rho = 2/3.
    std::vector<double> alpha = {1.0, 0.0};
    auto r = line::qsys::qsys_phmc(alpha, erlang2(4.0), 1.5, 2u);
    CHECK(r.utilization == doctest::Approx(2.0 / 3.0).epsilon(1e-15));
    // The R fixed point is driven to 1e-14, so 1e-9 relative is claimed.
    CHECK(r.meanQueueLength == doctest::Approx(2.0220560277297017).epsilon(1e-9));
    CHECK(r.meanWaitingQueue == doctest::Approx(0.68872269439637501).epsilon(1e-9));
    CHECK(r.meanWaitingTime == doctest::Approx(0.34436134719818751).epsilon(1e-9));
    CHECK(r.meanSojournTime == doctest::Approx(1.0110280138648542).epsilon(1e-9));
    CHECK(r.meanSojournTime == doctest::Approx(r.meanWaitingTime + 1.0 / 1.5).epsilon(1e-12));
    CHECK(r.meanWaitingQueue == doctest::Approx(2.0 * r.meanWaitingTime).epsilon(1e-12));
    // Erlang-2 arrivals are less bursty than Poisson at the same rate.
    CHECK(r.meanSojournTime < line::qsys::qsys_mmk(2.0, 1.5, 2u).W);
}

TEST_CASE("PH/M/c rejects an unstable instance") {
    CHECK_THROWS_AS(line::qsys::qsys_phmc(std::vector<double>{1.0}, exp_ph(4.0), 1.0, 2u),
                    line::InputError);  // rho = 2
    CHECK_THROWS_AS(line::qsys::qsys_phmc(std::vector<double>{1.0}, exp_ph(1.0), 1.0, 0u),
                    line::InputError);  // c = 0
}

TEST_CASE("D/M/1 agrees with the G/M/1 closed form at the deterministic root") {
    // For deterministic interarrivals the G/M/1 root solves
    // sigma = exp(-mu(1-sigma)/lambda); with lambda = 0.5, mu = 1 it is
    // found by bisection here and fed to qsys_gm1.
    const double lambda = 0.5, mu = 1.0;
    double lo = 1e-12, hi = 1.0 - 1e-12;
    for (int i = 0; i < 200; ++i) {
        const double m = 0.5 * (lo + hi);
        (m - std::exp(-mu * (1.0 - m) / lambda) < 0.0 ? lo : hi) = m;
    }
    const double sigma = 0.5 * (lo + hi);
    const double W_exact = line::qsys::qsys_gm1(sigma, mu);

    auto r = line::qsys::qsys_dmc(lambda, mu, 1u, 60u, 200u);
    CHECK(r.utilization == doctest::Approx(0.5).epsilon(1e-15));
    // The reference integrates the cycle with 200 trapezoid steps, which is
    // the dominant error; agreement with the exact G/M/1 value is only 1e-5.
    CHECK(r.meanSojournTime == doctest::Approx(W_exact).epsilon(1e-4));
    // Against MATLAB, which runs the same truncation and the same 200 steps,
    // the agreement is that of expm against expm: 1e-10 relative.
    CHECK(r.meanQueueLength == doctest::Approx(0.62750380750500567).epsilon(1e-10));
    CHECK(r.meanWaitingQueue == doctest::Approx(0.12750116205125941).epsilon(1e-10));
    CHECK(r.meanWaitingTime == doctest::Approx(0.25500232410251883).epsilon(1e-10));
    CHECK(r.meanSojournTime == doctest::Approx(1.2550023241025188).epsilon(1e-10));
}

TEST_CASE("D/M/c matches MATLAB and beats M/M/c") {
    // lambda = 1.2, mu = 1, c = 2, rho = 0.6.
    auto r = line::qsys::qsys_dmc(1.2, 1.0, 2u, 60u, 200u);
    CHECK(r.utilization == doctest::Approx(0.6).epsilon(1e-14));
    CHECK(r.meanQueueLength == doctest::Approx(1.3823737803185796).epsilon(1e-10));
    CHECK(r.meanWaitingQueue == doctest::Approx(0.18237288274621721).epsilon(1e-10));
    CHECK(r.meanWaitingTime == doctest::Approx(0.15197740228851436).epsilon(1e-10));
    CHECK(r.meanSojournTime == doctest::Approx(1.1519774022885143).epsilon(1e-10));
    CHECK(r.meanWaitingTime == doctest::Approx(r.meanWaitingQueue / 1.2).epsilon(1e-13));
    // Deterministic arrivals are the least bursty renewal stream at a given
    // rate, so the delay is below the M/M/c one.
    CHECK(r.meanSojournTime < line::qsys::qsys_mmk(1.2, 1.0, 2u).W);
}

TEST_CASE("D/M/c rejects an unstable or degenerate instance") {
    CHECK_THROWS_AS(line::qsys::qsys_dmc(3.0, 1.0, 2u, 40u, 50u), line::InputError);  // rho = 1.5
    CHECK_THROWS_AS(line::qsys::qsys_dmc(0.5, -1.0, 1u, 40u, 50u), line::InputError);
    CHECK_THROWS_AS(line::qsys::qsys_dmc(0.5, 1.0, 1u, 40u, 0u), line::InputError);
}
