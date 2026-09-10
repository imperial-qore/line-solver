/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
/**
 * Morrison's heavy-usage expansion for a closed think+DPS network
 * (line/api/npfqn/npfqn_dps_morrison.h).
 *
 * The pinned values are the MATLAB reference (matlab/src/api/npfqn/npfqn_dps_morrison.m) and agree
 * with the JAR and python ports to the digits written here; the same four models appear in the
 * test of each codebase, so a divergence in any of them fails somewhere.
 *
 * Reference: J.A. Morrison, Queueing Systems 9 (1991) 191-214.
 */

#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/npfqn/npfqn_dps_morrison.h"
#include "line/num/number.h"

using namespace line;

namespace {

// M2: two classes, K = 100 each, rho = 0.9 -- the regime the expansion is derived for.
const std::vector<double> M2_N = {100.0, 100.0};
const std::vector<double> M2_Z = {1.0, 0.5};
const std::vector<double> M2_S = {0.006, 0.0015};
const std::vector<double> M2_W = {1.0, 4.0};

// M3: three classes, unequal populations and weights, rho = 1.06 (appendix A's saturated side).
const std::vector<double> M3_N = {20.0, 12.0, 28.0};
const std::vector<double> M3_Z = {1.0, 0.5, 2.0};
const std::vector<double> M3_S = {0.02, 0.01, 0.03};
const std::vector<double> M3_W = {1.0, 4.0, 2.0};

void check_close(const std::vector<double>& got, const std::vector<double>& want, double tol) {
    REQUIRE(got.size() == want.size());
    for (std::size_t i = 0; i < want.size(); ++i) {
        CAPTURE(i);
        CAPTURE(got[i]);
        CAPTURE(want[i]);
        CHECK(std::fabs(got[i] - want[i]) <= tol);
    }
}

}  // namespace

TEST_CASE("npfqn_dps_morrison matches the MATLAB reference, two classes") {
    const auto res = npfqn::npfqn_dps_morrison(M2_N, M2_Z, M2_S, M2_W);
    CHECK(std::fabs(res.rho - 0.9) < 1e-12);
    CHECK(std::fabs(res.a - 0.1) < 1e-12);
    check_close(res.Q, {4.053964005705, 0.930378969124}, 1e-9);
    check_close(res.R, {0.042452089190, 0.004666835854}, 1e-9);
    check_close(res.X, {95.946035994295, 198.139242061752}, 1e-8);
    check_close(res.Qlead, {4.373155763276, 0.546644470409}, 1e-9);
    check_close(res.sigma, {0.084297520661, -0.337190082645}, 1e-9);
}

TEST_CASE("npfqn_dps_morrison matches the MATLAB reference, three classes") {
    const auto res = npfqn::npfqn_dps_morrison(M3_N, M3_Z, M3_S, M3_W);
    CHECK(std::fabs(res.rho - 1.06) < 1e-12);
    check_close(res.Q, {3.213021210443, 0.918773941158, 2.645208279899}, 1e-9);
    check_close(res.R, {0.208463525546, 0.039776387080, 0.202390704352}, 1e-9);
    check_close(res.sigma, {0.497718870654, -0.218459289315, -0.258992817331}, 1e-9);
}

TEST_CASE("Little's law is exact, and the two response-time forms converge") {
    // Throughput closes the think station EXACTLY, which is why the analyzers report R = Q/T.
    // Morrison's RESULT 2 (4.17) is the EXPANDED ratio (4.11)/(4.15), so it differs from Q/X by
    // higher-order terms -- 22% at K=6, 0.5% at K=100, 0.09% at K=400. The two agree only
    // asymptotically, and that convergence is what is asserted here.
    const auto res = npfqn::npfqn_dps_morrison(M3_N, M3_Z, M3_S, M3_W);
    for (std::size_t j = 0; j < M3_N.size(); ++j) {
        CAPTURE(j);
        CHECK(std::fabs(res.X[j] - (M3_N[j] - res.Q[j]) / M3_Z[j]) < 1e-10);
    }

    auto gap = [](const npfqn::DpsMorrisonResult<double>& k) {
        double g = 0;
        for (std::size_t j = 0; j < k.R.size(); ++j) {
            g = std::max(g, std::fabs(k.R[j] - k.Q[j] / k.X[j]) / k.R[j]);
        }
        return g;
    };
    const double g6 = gap(npfqn::npfqn_dps_morrison(std::vector<double>{6.0, 6.0}, M2_Z,
                                                    std::vector<double>{0.075, 0.0375}, M2_W));
    const double g100 = gap(npfqn::npfqn_dps_morrison(M2_N, M2_Z, M2_S, M2_W));
    const double g400 = gap(npfqn::npfqn_dps_morrison(std::vector<double>{400.0, 400.0}, M2_Z,
                                                      std::vector<double>{0.0015, 0.000375}, M2_W));
    CHECK(g400 < g100);
    CHECK(g100 < g6);
    CHECK(g100 < 0.01);
}

TEST_CASE("equal weights collapse the discriminatory correction, Morrison (4.18)-(4.21)") {
    // w = 1 makes the network product-form: D = B, Q = C, M = H = I, L = J, hence R = S = V = 0,
    // A = -K and sigma = 0. Every codebase must reproduce these identities exactly.
    const std::vector<double> w1 = {1.0, 1.0};
    const auto res = npfqn::npfqn_dps_morrison(M2_N, M2_Z, M2_S, w1);
    CHECK(std::fabs(res.cD - res.cB) < 1e-9 * res.cB);
    CHECK(std::fabs(res.cQ - res.cC) < 1e-9 * res.cC);
    CHECK(std::fabs(res.cH - res.cI) < 1e-9 * res.cH);
    CHECK(std::fabs(res.cM - res.cH) < 1e-9 * res.cH);
    CHECK(std::fabs(res.cL - res.cJ) < 1e-9 * res.cL);
    CHECK(std::fabs(res.cR) < 1e-9);
    CHECK(std::fabs(res.cS) < 1e-9);
    CHECK(std::fabs(res.cV) < 1e-9);
    CHECK(std::fabs(res.cU) < 1e-9);
    CHECK(std::fabs(res.cA + res.cK) < 1e-9 * res.cK);
    CHECK(std::fabs(res.delta) < 1e-9);
    for (std::size_t j = 0; j < w1.size(); ++j) {
        CAPTURE(j);
        CHECK(std::fabs(res.sigma[j]) < 1e-9);
    }
    check_close(res.Q, {3.737598812660, 1.937595955666}, 1e-9);
}

TEST_CASE("the kernel refuses inputs it has no derivation for") {
    // A refusal is the right answer where the expansion has none; returning a number would be
    // worse. The solver-level gates (a DPS model off the shape, and "morrison" named on a model
    // with no DPS station at all) are covered by the MATLAB, JAR and python suites, which have a
    // model builder to hand; here the kernel's own guards stand for them.
    CHECK_THROWS(npfqn::npfqn_dps_morrison(M2_N, M2_Z, M2_S, std::vector<double>{1.0, -4.0}));
    CHECK_THROWS(npfqn::npfqn_dps_morrison(std::vector<double>{0.0, 100.0}, M2_Z, M2_S, M2_W));
    CHECK_THROWS(npfqn::npfqn_dps_morrison(M2_N, M2_Z, M2_S, std::vector<double>{1.0}));
    CHECK_THROWS(npfqn::npfqn_dps_morrison(M2_N, std::vector<double>{1.0, 0.0}, M2_S, M2_W));
}

TEST_CASE("the approximation improves with the population at fixed usage") {
    // Both models run at rho = 0.9 with the same weights and think times; only the populations
    // differ. The exact values are from the CTMC of eq. (2.1). An asymptotic expansion has to get
    // relatively better as N grows -- that is the whole claim being made.
    const std::vector<double> smallN = {6.0, 6.0};
    const std::vector<double> smallS = {0.075, 0.0375};
    const std::vector<double> exact_small = {1.178868777433, 0.746116560046};
    const std::vector<double> exact_big = {4.082236085180, 0.862829749692};

    const auto small = npfqn::npfqn_dps_morrison(smallN, M2_Z, smallS, M2_W);
    const auto big = npfqn::npfqn_dps_morrison(M2_N, M2_Z, M2_S, M2_W);
    for (std::size_t j = 0; j < 2; ++j) {
        CAPTURE(j);
        const double e_small = std::fabs(small.Q[j] - exact_small[j]) / exact_small[j];
        const double e_big = std::fabs(big.Q[j] - exact_big[j]) / exact_big[j];
        CHECK(e_big < e_small);
        CHECK(e_big < 0.1);
    }
}
