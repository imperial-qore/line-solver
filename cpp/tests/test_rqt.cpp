/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
/**
 * Robust Queueing Theory: the worst-case system time of a G/G/k queue, the
 * service adaptation of Section 7.1, and the network characterization of
 * Theorem 10.
 *
 * Reference values come from the MATLAB reference (matlab/src/api/qsys/
 * qsys_gigk_rqt.m and matlab/src/api/npfqn/npfqn_traffic_rqt.m) and, for the
 * calculus operators, from closed forms the paper states directly: passage
 * through a queue leaves the uncertainty set unchanged, a split by f scales
 * Gamma by f^(-1/alpha), and a merge of p identical streams by p^((1-alpha)/alpha).
 *
 * The tolerances below are ABSOLUTE, so they are written as an explicit
 * fabs difference rather than doctest::Approx, whose epsilon is relative and
 * would loosen the large-magnitude cases by orders of magnitude.
 */

#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/npfqn/npfqn_traffic_rqt.h"
#include "line/api/qsys/qsys_gig1_rqt.h"
#include "line/api/qsys/qsys_gigk_rqt.h"
#include "line/api/qsys/qsys_gigk_rqt_gamma.h"
#include "line/util/matrix.h"

using namespace line;

TEST_CASE("qsys_gig1_rqt matches the MATLAB M/M/1 reference") {
    const double mu = 1.0;
    const double rho[4] = {0.5, 0.7, 0.9, 0.95};
    const double Wref[4] = {3.04, 3.85523809523809, 10.4711111111111, 20.8126315789474};
    const double Sref[4] = {2.056013895625, 3.54825301974245, 10.410484762801, 20.7854598611252};
    for (int i = 0; i < 4; ++i) {
        CAPTURE(rho[i]);
        const double lambda = rho[i] * mu;
        const double sa = 1.0 / lambda, ss = 1.0 / mu;
        const double Gs = qsys::qsys_gigk_rqt_gamma(rho[i], mu, sa, ss, 1u, 2.0, "independent");
        const qsys::GigkRqtResult<double> r = qsys::qsys_gig1_rqt(lambda, mu, sa, Gs);
        CHECK(std::fabs(r.W - Wref[i]) <= 1e-9);
        CHECK(std::fabs(r.Sworst - Sref[i]) <= 1e-8);
        // the closed form bounds the exact worst case
        CHECK(r.W >= r.Sworst - 1e-9);
    }
}

TEST_CASE("qsys_gigk_rqt matches the MATLAB M/M/k reference") {
    const std::size_t ks[4] = {1, 3, 6, 10};
    const double Wref[4] = {10.4711111111111, 4.05111111111111, 2.44611111111, 1.80411111111111};
    for (int i = 0; i < 4; ++i) {
        const std::size_t k = ks[i];
        CAPTURE(k);
        const double mu = 1.0, rho = 0.9;
        const double lambda = rho * k * mu;
        const double sa = 1.0 / lambda, ss = 1.0 / mu;
        const double Gs = qsys::qsys_gigk_rqt_gamma(rho, mu, sa, ss, k, 2.0, "independent");
        const qsys::GigkRqtResult<double> r = qsys::qsys_gigk_rqt(lambda, mu, sa, Gs, k, 2.0, 2.0);
        CHECK(std::fabs(r.W - Wref[i]) <= 1e-9);
    }
}

TEST_CASE("a heavier tail raises the robust system time") {
    double prev = 0.0;
    const double al[3] = {2.0, 1.8, 1.5};
    for (int i = 0; i < 3; ++i) {
        CAPTURE(al[i]);
        const qsys::GigkRqtResult<double> r = qsys::qsys_gig1_rqt(0.9, 1.0, 1.0, 1.0, al[i], al[i]);
        if (i > 0) {
            CHECK(r.W > prev);
        }
        prev = r.W;
    }
}

TEST_CASE("npfqn_traffic_rqt: passage through a queue leaves the set unchanged") {
    // tandem: the robust Burke theorem leaves the set unchanged
    std::vector<double> l0(2, 0.0), G0(2, 0.0), a0(2, 2.0);
    l0[0] = 1.0;
    G0[0] = 0.5;
    Matrix<double> F(2, 2, 0.0);
    F(0, 1) = 1.0;
    npfqn::TrafficRqt<double> t = npfqn::npfqn_traffic_rqt(l0, G0, a0, F);
    CHECK(std::fabs(t.lambda[1] - 1.0) <= 1e-12);
    CHECK(std::fabs(t.Gamma[1] - 0.5) <= 1e-12);
}

TEST_CASE("npfqn_traffic_rqt: a split scales Gamma by f^(-1/alpha)") {
    std::vector<double> l0(3, 0.0), G0(3, 0.0), a0(3, 2.0);
    l0[0] = 1.0;
    G0[0] = 0.5;
    Matrix<double> F(3, 3, 0.0);
    F(0, 1) = 0.25;
    F(0, 2) = 0.75;
    npfqn::TrafficRqt<double> t = npfqn::npfqn_traffic_rqt(l0, G0, a0, F);
    CHECK(std::fabs(t.lambda[1] - 0.25) <= 1e-12);
    CHECK(std::fabs(t.Gamma[1] - 0.5 / std::sqrt(0.25)) <= 1e-12);
    CHECK(std::fabs(t.Gamma[2] - 0.5 / std::sqrt(0.75)) <= 1e-12);
}

TEST_CASE("npfqn_traffic_rqt: a merge of two identical streams") {
    std::vector<double> l0(3, 0.0), G0(3, 0.0), a0(3, 2.0);
    l0[0] = l0[1] = 1.0;
    G0[0] = G0[1] = 0.5;
    Matrix<double> F(3, 3, 0.0);
    F(0, 2) = 1.0;
    F(1, 2) = 1.0;
    npfqn::TrafficRqt<double> t = npfqn::npfqn_traffic_rqt(l0, G0, a0, F);
    CHECK(std::fabs(t.lambda[2] - 2.0) <= 1e-12);
    CHECK(std::fabs(t.Gamma[2] - std::sqrt(2.0 * 0.25) / 2.0) <= 1e-12);
}

TEST_CASE("npfqn_traffic_rqt: the heaviest upstream tail dominates downstream") {
    std::vector<double> l0(3, 0.0), G0(3, 0.0), a0(3, 2.0);
    l0[0] = l0[1] = 1.0;
    G0[0] = G0[1] = 0.5;
    a0[0] = 1.5;
    Matrix<double> F(3, 3, 0.0);
    F(0, 2) = 1.0;
    F(1, 2) = 1.0;
    npfqn::TrafficRqt<double> t = npfqn::npfqn_traffic_rqt(l0, G0, a0, F);
    CHECK(std::fabs(t.alpha[2] - 1.5) <= 1e-12);
    // only the alpha=1.5 stream contributes, Theorem 9
    const double z = std::pow(1.0 * 0.5, 1.5 / 0.5);
    CHECK(std::fabs(t.Gamma[2] - std::pow(z, 0.5 / 1.5) / 2.0) <= 1e-12);
}
