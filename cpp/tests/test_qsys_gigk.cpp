/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Multiserver G/I/G/k approximations. The structural checks are that the
 * Kingman (Lee-Longton) scaling collapses onto its G/I/G/1 counterpart at
 * k = 1 in exact arithmetic, and that Kingman, Cosmetatos and Whitt all
 * return the exact M/M/k answer at ca = cs = 1.
 */
#include <cmath>

#include "doctest.h"
#include "line/api/qsys/qsys_gig1_approx_allencunneen.h"
#include "line/api/qsys/qsys_gigk_approx.h"
#include "line/api/qsys/qsys_gigk_approx_cosmetatos.h"
#include "line/api/qsys/qsys_gigk_approx_kingman.h"
#include "line/api/qsys/qsys_gigk_approx_whitt.h"
#include "line/api/qsys/qsys_mmk.h"

using line::Rational;
using line::Real50;

namespace {
constexpr double REL_TOL = 1e-9;

// Working point: lambda = 2, mu = 1, k = 3 (rho = 2/3), ca = 3/2, cs = 1/2.
constexpr double LAM = 2.0, MU = 1.0, CA = 1.5, CS = 0.5;
constexpr unsigned K = 3;
}  // namespace

TEST_CASE("qsys_gigk_approx_kingman with k=1 equals its G/I/G/1 counterpart exactly") {
    const Rational lam(1, 2), mu(1), ca(3, 2), cs(1, 2);
    auto k1 = line::qsys::qsys_gigk_approx_kingman(lam, mu, ca, cs, 1u);
    auto g1 = line::qsys::qsys_gig1_approx_allencunneen(lam, mu, ca, cs);
    CHECK(k1.W == g1.W);
    CHECK(k1.rhohat == g1.rhohat);

    // Second operating point, and the degenerate ca = cs = 1 case.
    const Rational lam2(1, 3), one(1);
    CHECK(line::qsys::qsys_gigk_approx_kingman(lam2, mu, one, one, 1u).W ==
          line::qsys::qsys_gig1_approx_allencunneen(lam2, mu, one, one).W);
}

TEST_CASE("qsys_gigk_approx_kingman at ca=cs=1 is the exact M/M/k answer") {
    const Rational lam(2), mu(1), one(1);
    CHECK(line::qsys::qsys_gigk_approx_kingman(lam, mu, one, one, K).W ==
          line::qsys::qsys_mmk(lam, mu, K).W);
}

TEST_CASE("qsys_gigk_approx_kingman exact and double agree") {
    const double d = line::qsys::qsys_gigk_approx_kingman(LAM, MU, CA, CS, K).W;
    const Rational q =
        line::qsys::qsys_gigk_approx_kingman(Rational(2), Rational(1), Rational(3, 2),
                                             Rational(1, 2), K)
            .W;
    CHECK(d == doctest::Approx(static_cast<double>(q)).epsilon(REL_TOL));

    const double r = static_cast<double>(
        line::qsys::qsys_gigk_approx_kingman(Real50(2), Real50(1), Real50(3) / Real50(2),
                                             Real50(1) / Real50(2), K)
            .W);
    CHECK(d == doctest::Approx(r).epsilon(REL_TOL));

    // Independent evaluation of the formula.
    const double Wmmk = line::qsys::qsys_mmk(LAM, MU, K).W;
    CHECK(d == doctest::Approx((CA * CA + CS * CS) / 2 * (Wmmk - 1.0 / MU) + 1.0 / MU).epsilon(REL_TOL));
}

TEST_CASE("qsys_gigk_approx matches the MATLAB branch structure") {
    // rho = 2/3 <= 0.7, so alpha = rho^((k+1)/2) = rho^2.
    const double rho = LAM / (MU * K);
    const double alpha_lo = std::pow(rho, (K + 1) / 2.0);
    const double Wlo =
        (alpha_lo / MU) * (1 / (1 - rho)) * (CA * CA + CS * CS) / (2.0 * K) + 1.0 / MU;
    CHECK(line::qsys::qsys_gigk_approx(LAM, MU, CA, CS, K).W == doctest::Approx(Wlo).epsilon(REL_TOL));

    // rho > 0.7 branch: lambda = 2.4, k = 3 gives rho = 0.8.
    const double lam2 = 2.4, rho2 = lam2 / (MU * K);
    const double alpha_hi = (std::pow(rho2, static_cast<double>(K)) + rho2) / 2;
    const double Whi =
        (alpha_hi / MU) * (1 / (1 - rho2)) * (CA * CA + CS * CS) / (2.0 * K) + 1.0 / MU;
    CHECK(line::qsys::qsys_gigk_approx(lam2, MU, CA, CS, K).W == doctest::Approx(Whi).epsilon(REL_TOL));

    // Real50 agrees with double.
    CHECK(static_cast<double>(line::qsys::qsys_gigk_approx(Real50(2), Real50(1), Real50(3) / Real50(2),
                                                           Real50(1) / Real50(2), K)
                                  .W) == doctest::Approx(Wlo).epsilon(REL_TOL));
}

TEST_CASE("qsys_gigk_approx_cosmetatos at ca=cs=1 is the exact M/M/k answer") {
    const double Wmmk = line::qsys::qsys_mmk(LAM, MU, K).W;
    CHECK(line::qsys::qsys_gigk_approx_cosmetatos(LAM, MU, 1.0, 1.0, K).W ==
          doctest::Approx(Wmmk).epsilon(REL_TOL));
}

TEST_CASE("qsys_gigk_approx_cosmetatos matches the MATLAB formula, both branches") {
    const double rho = LAM / (K * MU);
    const double Wq_mmk = line::qsys::qsys_mmk(LAM, MU, K).W - 1.0 / MU;
    const double gamma =
        std::min(0.24, (1 - rho) * (K - 1) * (std::sqrt(4.0 + 5.0 * K) - 2.0) / (16.0 * K * rho));
    const double phi1 = 1 + gamma;
    const double phi3 = (1 - 4 * gamma) * std::exp(-2 * (1 - rho) / (3 * rho));

    // Inside the unit box: ca = 0.8, cs = 0.5.
    const double ca_in = 0.8, ca2 = 0.64, cs2 = 0.25;
    const double Wq_in = (ca2 * cs2 + ca2 * (1 - cs2) * phi1 / 2 + (1 - ca2) * cs2 * phi3 / 2) * Wq_mmk;
    CHECK(line::qsys::qsys_gigk_approx_cosmetatos(LAM, MU, ca_in, CS, K).W ==
          doctest::Approx(Wq_in + 1.0 / MU).epsilon(REL_TOL));

    // Outside the unit box (ca = 1.5): Lee-Longton fallback.
    const double Wq_out = ((CA * CA + CS * CS) / 2) * Wq_mmk;
    CHECK(line::qsys::qsys_gigk_approx_cosmetatos(LAM, MU, CA, CS, K).W ==
          doctest::Approx(Wq_out + 1.0 / MU).epsilon(REL_TOL));

    // Real50 agrees with double.
    CHECK(static_cast<double>(
              line::qsys::qsys_gigk_approx_cosmetatos(Real50(2), Real50(1), Real50(3) / Real50(2),
                                                      Real50(1) / Real50(2), K)
                  .W) == doctest::Approx(Wq_out + 1.0 / MU).epsilon(REL_TOL));
}

TEST_CASE("qsys_gigk_approx_whitt at ca=cs=1 is the exact M/M/k answer") {
    const double Wmmk = line::qsys::qsys_mmk(LAM, MU, K).W;
    CHECK(line::qsys::qsys_gigk_approx_whitt(LAM, MU, 1.0, 1.0, K).W ==
          doctest::Approx(Wmmk).epsilon(REL_TOL));
}

TEST_CASE("qsys_gigk_approx_whitt matches the MATLAB formula on all three phi branches") {
    const double rho = LAM / (K * MU);
    const double Wq_mmk = line::qsys::qsys_mmk(LAM, MU, K).W - 1.0 / MU;
    const double gamma =
        std::min(0.24, (1 - rho) * (K - 1) * (std::sqrt(4.0 + 5.0 * K) - 2.0) / (16.0 * K * rho));
    const double phi1 = 1 + gamma;
    const double phi3 = (1 - 4 * gamma) * std::exp(-2 * (1 - rho) / (3 * rho));
    const double phi4 = std::min(1.0, (phi1 + phi3) / 2);

    auto whitt_ref = [&](double ca, double cs) {
        const double ca2 = ca * ca, cs2 = cs * cs, c2 = (ca2 + cs2) / 2;
        const double psi = c2 >= 1 ? 1.0 : std::pow(phi4, 2 * (1 - c2));
        double phi;
        if (std::fabs(ca2 - cs2) < 1e-12) {
            phi = psi;
        } else if (ca2 > cs2) {
            phi = (4 * (ca2 - cs2) / (4 * ca2 - 3 * cs2)) * phi1 + (cs2 / (4 * ca2 - 3 * cs2)) * psi;
        } else {
            phi = ((cs2 - ca2) / (2 * (ca2 + cs2))) * phi3 + ((cs2 + 3 * ca2) / (2 * (ca2 + cs2))) * psi;
        }
        return phi * c2 * Wq_mmk + 1.0 / MU;
    };

    // ca2 > cs2
    CHECK(line::qsys::qsys_gigk_approx_whitt(LAM, MU, CA, CS, K).W ==
          doctest::Approx(whitt_ref(CA, CS)).epsilon(REL_TOL));
    // ca2 < cs2
    CHECK(line::qsys::qsys_gigk_approx_whitt(LAM, MU, 0.5, 1.5, K).W ==
          doctest::Approx(whitt_ref(0.5, 1.5)).epsilon(REL_TOL));
    // ca2 == cs2, and below the c2 = 1 threshold so psi is a real power
    CHECK(line::qsys::qsys_gigk_approx_whitt(LAM, MU, 0.5, 0.5, K).W ==
          doctest::Approx(whitt_ref(0.5, 0.5)).epsilon(REL_TOL));

    // Real50 agrees with double.
    CHECK(static_cast<double>(
              line::qsys::qsys_gigk_approx_whitt(Real50(2), Real50(1), Real50(3) / Real50(2),
                                                 Real50(1) / Real50(2), K)
                  .W) == doctest::Approx(whitt_ref(CA, CS)).epsilon(REL_TOL));
}

TEST_CASE("gigk W increases with rho for every approximation") {
    double prevK = 0, prevA = 0, prevC = 0, prevW = 0;
    for (double lam = 0.6; lam < 2.9; lam += 0.4) {
        const double kg = line::qsys::qsys_gigk_approx_kingman(lam, MU, CA, CS, K).W;
        const double ap = line::qsys::qsys_gigk_approx(lam, MU, CA, CS, K).W;
        const double co = line::qsys::qsys_gigk_approx_cosmetatos(lam, MU, CA, CS, K).W;
        const double wh = line::qsys::qsys_gigk_approx_whitt(lam, MU, CA, CS, K).W;
        CHECK(kg > prevK);
        CHECK(ap > prevA);
        CHECK(co > prevC);
        CHECK(wh > prevW);
        prevK = kg;
        prevA = ap;
        prevC = co;
        prevW = wh;
    }
}
