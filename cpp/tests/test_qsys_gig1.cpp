/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * G/I/G/1 approximations and bounds. The load-bearing checks are the exact
 * reductions: at ca = cs = 1 the Allen-Cunneen, Heyman, Kimura, Marchal and
 * Kraemer/Langenbach-Belz formulas must return the M/M/1 response time, and
 * the first four do so in exact rational arithmetic. Kobayashi, Gelenbe and
 * the two bounds do not reduce, so they are checked against the formulas
 * evaluated inline.
 */
#include <cmath>

#include "doctest.h"
#include "line/api/qsys/qsys_gig1_approx_allencunneen.h"
#include "line/api/qsys/qsys_gig1_approx_gelenbe.h"
#include "line/api/qsys/qsys_gig1_approx_heyman.h"
#include "line/api/qsys/qsys_gig1_approx_kimura.h"
#include "line/api/qsys/qsys_gig1_approx_klb.h"
#include "line/api/qsys/qsys_gig1_approx_kobayashi.h"
#include "line/api/qsys/qsys_gig1_approx_marchal.h"
#include "line/api/qsys/qsys_gig1_lbnd.h"
#include "line/api/qsys/qsys_gig1_ubnd_kingman.h"
#include "line/api/qsys/qsys_mm1.h"

using line::Rational;
using line::Real50;

namespace {
constexpr double REL_TOL = 1e-9;

// Working point used throughout: lambda = 1/2, mu = 1, ca = 3/2, cs = 1/2.
constexpr double LAM = 0.5, MU = 1.0, CA = 1.5, CS = 0.5;

Rational lam_q() { return Rational(1, 2); }
Rational mu_q() { return Rational(1); }
Rational ca_q() { return Rational(3, 2); }
Rational cs_q() { return Rational(1, 2); }
}  // namespace

TEST_CASE("gig1 approximations reduce to M/M/1 exactly at ca=cs=1") {
    const Rational lam(1, 3), mu(1), one(1);
    const Rational Wmm1 = line::qsys::qsys_mm1(lam, mu).W;

    CHECK(line::qsys::qsys_gig1_approx_allencunneen(lam, mu, one, one).W == Wmm1);
    CHECK(line::qsys::qsys_gig1_approx_heyman(lam, mu, one, one).W == Wmm1);
    CHECK(line::qsys::qsys_gig1_approx_kimura(lam, mu, one, one).W == Wmm1);
    CHECK(line::qsys::qsys_gig1_approx_marchal(lam, mu, one, one).W == Wmm1);

    // KLB needs exp for its correction factor g, which is exactly 1 at ca = 1,
    // so it reduces too, but only in a transcendental arithmetic.
    const double Wd = line::qsys::qsys_mm1(1.0 / 3.0, 1.0).W;
    CHECK(line::qsys::qsys_gig1_approx_klb(1.0 / 3.0, 1.0, 1.0, 1.0).W ==
          doctest::Approx(Wd).epsilon(REL_TOL));
}

TEST_CASE("allencunneen and heyman are the same formula, exact and double") {
    auto a = line::qsys::qsys_gig1_approx_allencunneen(lam_q(), mu_q(), ca_q(), cs_q());
    auto h = line::qsys::qsys_gig1_approx_heyman(lam_q(), mu_q(), ca_q(), cs_q());
    CHECK(a.W == h.W);

    // rho/(1-rho)/mu*(ca^2+cs^2)/2 + 1/mu = 1*(2.25+0.25)/2 + 1 = 2.25.
    CHECK(a.W == Rational(9, 4));
    CHECK(line::qsys::qsys_gig1_approx_allencunneen(LAM, MU, CA, CS).W ==
          doctest::Approx(2.25).epsilon(REL_TOL));
    CHECK(static_cast<double>(
              line::qsys::qsys_gig1_approx_allencunneen(Real50(1) / Real50(2), Real50(1),
                                                        Real50(3) / Real50(2), Real50(1) / Real50(2))
                  .W) == doctest::Approx(2.25).epsilon(REL_TOL));
}

TEST_CASE("gig1 exact and double agree to 1e-9 relative on every field-only function") {
    struct Case {
        const char* name;
        double d;
        Rational q;
    };
    const Case cases[] = {
        {"allencunneen", line::qsys::qsys_gig1_approx_allencunneen(LAM, MU, CA, CS).W,
         line::qsys::qsys_gig1_approx_allencunneen(lam_q(), mu_q(), ca_q(), cs_q()).W},
        {"heyman", line::qsys::qsys_gig1_approx_heyman(LAM, MU, CA, CS).W,
         line::qsys::qsys_gig1_approx_heyman(lam_q(), mu_q(), ca_q(), cs_q()).W},
        {"kimura", line::qsys::qsys_gig1_approx_kimura(LAM, MU, CA, CS).W,
         line::qsys::qsys_gig1_approx_kimura(lam_q(), mu_q(), ca_q(), cs_q()).W},
        {"marchal", line::qsys::qsys_gig1_approx_marchal(LAM, MU, CA, CS).W,
         line::qsys::qsys_gig1_approx_marchal(lam_q(), mu_q(), ca_q(), cs_q()).W},
        {"ubnd_kingman", line::qsys::qsys_gig1_ubnd_kingman(LAM, MU, CA, CS).W,
         line::qsys::qsys_gig1_ubnd_kingman(lam_q(), mu_q(), ca_q(), cs_q()).W},
        {"lbnd", line::qsys::qsys_gig1_lbnd(LAM, MU, CA, CS).W,
         line::qsys::qsys_gig1_lbnd(lam_q(), mu_q(), ca_q(), cs_q()).W},
    };
    for (const Case& c : cases) {
        INFO("function: ", c.name);
        CHECK(c.d == doctest::Approx(static_cast<double>(c.q)).epsilon(REL_TOL));
    }
}

TEST_CASE("gig1 transcendental approximations agree between double and Real50") {
    CHECK(static_cast<double>(line::qsys::qsys_gig1_approx_kobayashi(
                                  Real50(1) / Real50(2), Real50(1), Real50(3) / Real50(2),
                                  Real50(1) / Real50(2))
                                  .W) ==
          doctest::Approx(line::qsys::qsys_gig1_approx_kobayashi(LAM, MU, CA, CS).W).epsilon(REL_TOL));
    CHECK(static_cast<double>(line::qsys::qsys_gig1_approx_gelenbe(
                                  Real50(1) / Real50(2), Real50(1), Real50(3) / Real50(2),
                                  Real50(1) / Real50(2))
                                  .W) ==
          doctest::Approx(line::qsys::qsys_gig1_approx_gelenbe(LAM, MU, CA, CS).W).epsilon(REL_TOL));
    CHECK(static_cast<double>(
              line::qsys::qsys_gig1_approx_klb(Real50(1) / Real50(2), Real50(1),
                                               Real50(3) / Real50(2), Real50(1) / Real50(2))
                  .W) ==
          doctest::Approx(line::qsys::qsys_gig1_approx_klb(LAM, MU, CA, CS).W).epsilon(REL_TOL));
}

TEST_CASE("gig1 approximations match the MATLAB formulas evaluated inline") {
    const double rho = LAM / MU, ca2 = CA * CA, cs2 = CS * CS;

    // Kobayashi
    const double rhohat = std::exp(-2 * (1 - rho) / (rho * (ca2 + cs2 / rho)));
    CHECK(line::qsys::qsys_gig1_approx_kobayashi(LAM, MU, CA, CS).W ==
          doctest::Approx(rhohat / (1 - rhohat) / LAM).epsilon(REL_TOL));

    // Gelenbe
    const double rhat = std::exp(-2 * (1 - rho) / (rho * ca2 + cs2));
    CHECK(line::qsys::qsys_gig1_approx_gelenbe(LAM, MU, CA, CS).W ==
          doctest::Approx(1.0 / (MU * (1 - rhat))).epsilon(REL_TOL));

    // KLB, ca > 1 branch
    const double g = std::exp(-(1 - rho) * (ca2 - 1) / (ca2 + 4 * cs2));
    CHECK(line::qsys::qsys_gig1_approx_klb(LAM, MU, CA, CS).W ==
          doctest::Approx(1.0 / MU * ((rho / (1 - rho)) * ((cs2 + ca2) / 2) * g + 1)).epsilon(REL_TOL));

    // KLB, ca <= 1 branch
    const double ca_lo = 0.5, ca_lo2 = 0.25;
    const double g2 = std::exp(-2 * (1 - rho) * (1 - ca_lo2) * (1 - ca_lo2) / (3 * rho * (ca_lo2 + cs2)));
    CHECK(line::qsys::qsys_gig1_approx_klb(LAM, MU, ca_lo, CS).W ==
          doctest::Approx(1.0 / MU * ((rho / (1 - rho)) * ((cs2 + ca_lo2) / 2) * g2 + 1)).epsilon(REL_TOL));

    // Kimura
    CHECK(line::qsys::qsys_gig1_approx_kimura(LAM, MU, CA, CS).W ==
          doctest::Approx(rho * (ca2 + cs2) / MU / (1 - rho) / (1 + ca2) + 1.0 / MU).epsilon(REL_TOL));

    // Marchal (ca appears unsquared in the correction ratio, as in MATLAB)
    const double Wmm1 = rho / (1 - rho);
    CHECK(line::qsys::qsys_gig1_approx_marchal(LAM, MU, CA, CS).W ==
          doctest::Approx(Wmm1 * (1 + cs2) / 2 / MU * (CA + rho * rho * cs2) / (1 + rho * rho * cs2) +
                          1.0 / MU)
              .epsilon(REL_TOL));

    // Kingman upper bound
    CHECK(line::qsys::qsys_gig1_ubnd_kingman(LAM, MU, CA, CS).W ==
          doctest::Approx(LAM * (ca2 / (LAM * LAM) + cs2 / (MU * MU)) / (2 * (1 - rho)) + 1.0 / MU)
              .epsilon(REL_TOL));

    // Lower bound is the mean service time
    CHECK(line::qsys::qsys_gig1_lbnd(LAM, MU, CA, CS).W == doctest::Approx(1.0 / MU).epsilon(REL_TOL));
}

TEST_CASE("gig1 W increases with rho for every approximation") {
    // The Kingman upper bound is excluded here and checked separately: it
    // carries a term ca^2/lambda that diverges as lambda -> 0, so the bound is
    // U-shaped in rho at fixed mu rather than monotone.
    double prevA = 0, prevH = 0, prevK = 0, prevM = 0, prevKo = 0, prevG = 0, prevL = 0;
    for (double lam = 0.1; lam < 0.95; lam += 0.1) {
        const double a = line::qsys::qsys_gig1_approx_allencunneen(lam, MU, CA, CS).W;
        const double h = line::qsys::qsys_gig1_approx_heyman(lam, MU, CA, CS).W;
        const double ki = line::qsys::qsys_gig1_approx_kimura(lam, MU, CA, CS).W;
        const double m = line::qsys::qsys_gig1_approx_marchal(lam, MU, CA, CS).W;
        const double ko = line::qsys::qsys_gig1_approx_kobayashi(lam, MU, CA, CS).W;
        const double g = line::qsys::qsys_gig1_approx_gelenbe(lam, MU, CA, CS).W;
        const double l = line::qsys::qsys_gig1_approx_klb(lam, MU, CA, CS).W;
        CHECK(a > prevA);
        CHECK(h > prevH);
        CHECK(ki > prevK);
        CHECK(m > prevM);
        CHECK(ko > prevKo);
        CHECK(g > prevG);
        CHECK(l > prevL);
        prevA = a;
        prevH = h;
        prevK = ki;
        prevM = m;
        prevKo = ko;
        prevG = g;
        prevL = l;
    }
}

TEST_CASE("gig1 Kingman upper bound increases with rho above its minimum") {
    // Wq = (ca^2/lambda + lambda cs^2/mu^2)/(2(1-rho)) has a 1/lambda term, so
    // it is only increasing once lambda is past the interior minimum. With
    // ca = 3/2, cs = 1/2, mu = 1 that minimum sits below lambda = 0.5.
    double prev = 0.0;
    for (double lam = 0.5; lam < 0.95; lam += 0.1) {
        const double u = line::qsys::qsys_gig1_ubnd_kingman(lam, MU, CA, CS).W;
        CHECK(u > prev);
        prev = u;
    }
}

TEST_CASE("gig1 bounds bracket the approximations at the working point") {
    const double lo = line::qsys::qsys_gig1_lbnd(LAM, MU, CA, CS).W;
    const double up = line::qsys::qsys_gig1_ubnd_kingman(LAM, MU, CA, CS).W;
    const double ac = line::qsys::qsys_gig1_approx_allencunneen(LAM, MU, CA, CS).W;
    CHECK(lo <= ac);
    CHECK(ac <= up);
}
