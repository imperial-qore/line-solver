/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
/**
 * Piecewise-linear Lyapunov bound of a multitype Markovian queueing network.
 *
 * Reference values come from the MATLAB reference
 * (matlab/src/api/npfqn/npfqn_bnd_bgt.m) and agree with the JAR and python
 * ports to the digits printed here. On M/M/1 the gamma of the LP is checked
 * against its closed form (mu-lambda)/(lambda+mu), which is what the drift of
 * L = 1 gives once the rates are uniformized.
 *
 * The tolerances below are ABSOLUTE, so they are written as an explicit
 * fabs difference rather than doctest::Approx, whose epsilon is relative and
 * would loosen the Lu-Kumar cases (Qub ~ 8e6 at 1e-4) by orders of magnitude.
 */

#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/npfqn/npfqn_bnd_bgt.h"
#include "line/num/number.h"

using namespace line;

TEST_CASE("npfqn_bnd_bgt on M/M/1 matches the closed form") {
    // J = 1, so the exception parameter vanishes
    const double rho[4] = {0.3, 0.5, 0.7, 0.9};
    const double gref[4] = {0.538461538461538, 1.0 / 3.0, 0.176470588235294, 0.0526315789473684};
    const double qref[4] = {23.9340659340659, 32.6666666666667, 53.6862745098039,
                            160.105263157895};
    for (int i = 0; i < 4; ++i) {
        CAPTURE(rho[i]);
        const std::vector<double> lam(1, rho[i]);
        const std::vector<std::vector<double> > mu(1, std::vector<double>(1, 1.0));
        const std::vector<std::vector<std::size_t> > sg(1, std::vector<std::size_t>(1, 0));
        const npfqn::BndBgt<double> r = npfqn::npfqn_bnd_bgt(lam, mu, sg, 1u);
        CHECK(std::fabs(r.gamma - gref[i]) <= 1e-9);
        // the closed form of the uniformized drift of L = 1
        CHECK(std::fabs(r.gamma - (1.0 - rho[i]) / (1.0 + rho[i])) <= 1e-9);
        CHECK(std::fabs(r.Qub[0][0] - qref[i]) <= 1e-8);
        // it IS an upper bound on the exact mean queue length
        CHECK(r.Qub[0][0] >= rho[i] / (1.0 - rho[i]));
    }
}

TEST_CASE("npfqn_bnd_bgt on a two-station tandem") {
    const std::vector<double> lam(1, 0.5);
    std::vector<std::vector<double> > mu(1);
    mu[0].push_back(1.0);
    mu[0].push_back(1.0);
    std::vector<std::vector<std::size_t> > sg(1);
    sg[0].push_back(0);
    sg[0].push_back(1);
    const npfqn::BndBgt<double> r = npfqn::npfqn_bnd_bgt(lam, mu, sg, 2u);
    CHECK(std::fabs(r.gamma - 0.2) <= 1e-9);
    CHECK(std::fabs(r.Lmax - 1.0) <= 1e-9);
    CHECK(std::fabs(r.B - 5529.6) <= 1e-6);
    CHECK(std::fabs(r.U - 5578.0) <= 1e-6);
    CHECK(std::fabs(r.Qub[0][0] - 5578.0) <= 1e-6);
    CHECK(std::fabs(r.Qub[0][1] - 5578.0) <= 1e-6);
    CHECK(std::fabs(r.tailRatio - 1.1 / 1.15) <= 1e-9);
    CHECK(std::fabs(r.tailStep - 2.2) <= 1e-9);
}

TEST_CASE("npfqn_bnd_bgt on the globally stable Lu-Kumar re-entrant line") {
    const std::vector<double> lam(1, 1.0);
    std::vector<std::vector<double> > mu(1);
    mu[0].push_back(1 / 0.3);
    mu[0].push_back(1 / 0.6);
    mu[0].push_back(1 / 0.3);
    mu[0].push_back(1 / 0.1);
    std::vector<std::vector<std::size_t> > sg(1);
    sg[0].push_back(0);
    sg[0].push_back(1);
    sg[0].push_back(1);
    sg[0].push_back(0);
    const npfqn::BndBgt<double> r = npfqn::npfqn_bnd_bgt(lam, mu, sg, 2u);
    CHECK(std::fabs(r.gamma - 0.00574712643678161) <= 1e-12);
    CHECK(std::fabs(r.rhoStation[0] - 0.4) <= 1e-12);
    CHECK(std::fabs(r.rhoStation[1] - 0.9) <= 1e-12);
    CHECK(std::fabs(r.Qub[0][0] - 7886457.48275862) <= 1e-4);
    CHECK(std::fabs(r.Qub[0][2] - 23659372.4482759) <= 1e-3);
}

TEST_CASE("npfqn_bnd_bgt refuses a globally unstable network at loads below one") {
    // rho_2 + rho_4 = 1.2 > 1 is the Rybko-Stolyar instability; every station
    // is at 0.7, so a per-station load test would accept it and the LP must not
    const std::vector<double> lam(1, 1.0);
    std::vector<std::vector<double> > mu(1);
    mu[0].push_back(1 / 0.1);
    mu[0].push_back(1 / 0.6);
    mu[0].push_back(1 / 0.1);
    mu[0].push_back(1 / 0.6);
    std::vector<std::vector<std::size_t> > sg(1);
    sg[0].push_back(0);
    sg[0].push_back(1);
    sg[0].push_back(1);
    sg[0].push_back(0);
    CHECK_THROWS_AS(npfqn::npfqn_bnd_bgt(lam, mu, sg, 2u), std::exception);
}

TEST_CASE("npfqn_bnd_bgt with two types crossing two stations") {
    const std::vector<double> lam(2, 0.3);
    std::vector<std::vector<double> > mu(2, std::vector<double>(2, 1.0));
    std::vector<std::vector<std::size_t> > sg(2, std::vector<std::size_t>(2, 0));
    sg[0][0] = 0;
    sg[0][1] = 1;
    sg[1][0] = 1;
    sg[1][1] = 0;
    const npfqn::BndBgt<double> r = npfqn::npfqn_bnd_bgt(lam, mu, sg, 2u);
    CHECK(std::fabs(r.gamma - 0.152173913043478) <= 1e-12);
    CHECK(std::fabs(r.U - 16969.7098491564) <= 1e-6);
}

TEST_CASE("npfqn_bnd_bgt on M/M/1 in exact arithmetic") {
    typedef Rational R;
    const std::vector<R> lam(1, R(7) / R(10));
    const std::vector<std::vector<R> > mu(1, std::vector<R>(1, R(1)));
    const std::vector<std::vector<std::size_t> > sg(1, std::vector<std::size_t>(1, 0));
    const npfqn::BndBgt<R> r = npfqn::npfqn_bnd_bgt(lam, mu, sg, 1u);
    // gamma = (mu-lambda)/(lambda+mu) = (3/10)/(17/10) = 3/17
    CHECK(r.gamma == R(3) / R(17));
}
