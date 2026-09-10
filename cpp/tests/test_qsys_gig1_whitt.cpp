/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * qsys_gig1_approx_whitt: Whitt's correction to the Kingman diffusion form.
 *
 * Goldens are native Python's on the same inputs. Beyond matching them, the
 * tests pin what the correction IS: a discount on Kingman that operates only
 * where Kingman is loose. With bursty service (cs^2 > 1) the correction factor
 * is exactly one and the two must coincide -- so `phi` is checked by asserting
 * that coincidence rather than by reading it out.
 *
 * The M/M/1 case is the anchor: at ca = cs = 1 the exact response time is
 * 1/(mu - lambda), and Whitt returns it while Kingman, being an upper bound,
 * does not.
 */
#include <cmath>

#include "doctest.h"
#include "line/api/qsys/qsys_gig1_approx_whitt.h"
#include "line/api/qsys/qsys_gig1_ubnd_kingman.h"

namespace qs = line::qsys;

TEST_CASE("M/M/1 is recovered exactly, where Kingman only bounds") {
    // lambda = 0.5, mu = 1: the exact mean response time is 1/(1-0.5) = 2.
    const qs::QsysResult<double> w = qs::qsys_gig1_approx_whitt(0.5, 1.0, 1.0, 1.0);
    CHECK(w.W == doctest::Approx(2.0).epsilon(1e-10));
    CHECK(w.rhohat == doctest::Approx(0.5).epsilon(1e-10));

    const qs::QsysResult<double> k = qs::qsys_gig1_ubnd_kingman(0.5, 1.0, 1.0, 1.0);
    CHECK(k.W == doctest::Approx(3.5).epsilon(1e-10));
    CHECK(w.W < k.W);  // the correction is a discount, and here a large one
}

TEST_CASE("the goldens match native Python") {
    const qs::QsysResult<double> a = qs::qsys_gig1_approx_whitt(0.5, 1.0, 0.5, 0.5);
    CHECK(a.W == doctest::Approx(1.118091638).epsilon(1e-8));
    CHECK(a.rhohat == doctest::Approx(0.3585820328).epsilon(1e-8));

    const qs::QsysResult<double> b = qs::qsys_gig1_approx_whitt(0.8, 1.0, 2.0, 0.5);
    CHECK(b.W == doctest::Approx(8.538823712).epsilon(1e-8));

    const qs::QsysResult<double> c = qs::qsys_gig1_approx_whitt(0.5, 1.0, 0.5, 2.0);
    CHECK(c.W == doctest::Approx(3.125).epsilon(1e-8));

    const qs::QsysResult<double> d = qs::qsys_gig1_approx_whitt(0.9, 1.0, 1.5, 1.5);
    CHECK(d.W == doctest::Approx(21.25).epsilon(1e-8));
}

TEST_CASE("with BURSTY SERVICE the correction is inert and this IS Kingman") {
    // cs^2 > 1 leaves phi at one in both remaining arms, so the queue-length
    // term is Kingman's. Checking the coincidence is how phi gets tested
    // without reading it out of the implementation.
    const double cs = 2.0;
    for (double ca = 0.5; ca <= 2.0; ca += 0.5) {
        const qs::QsysResult<double> w = qs::qsys_gig1_approx_whitt(0.5, 1.0, ca, cs);
        // Kingman's W is the waiting term plus one service time; Whitt's is
        // (Lq + rho)/lambda, which is the same quantity when phi = 1.
        const double rho = 0.5;
        const double Lq = rho * rho * (ca * ca + cs * cs) / (2.0 * (1.0 - rho));
        CHECK(w.W == doctest::Approx((Lq + rho) / 0.5).epsilon(1e-9));
    }
}

TEST_CASE("the correction only ever discounts, never inflates") {
    for (double ca = 0.25; ca <= 2.0; ca += 0.25)
        for (double cs = 0.25; cs <= 2.0; cs += 0.25)
            for (double lam = 0.2; lam <= 0.8; lam += 0.3) {
                const qs::QsysResult<double> w = qs::qsys_gig1_approx_whitt(lam, 1.0, ca, cs);
                const double rho = lam;
                const double Lq = rho * rho * (ca * ca + cs * cs) / (2.0 * (1.0 - rho));
                // phi <= 1 always, so the answer never exceeds the uncorrected
                // diffusion form.
                CHECK(w.W <= (Lq + rho) / lam + 1e-9);
                CHECK(w.W > 0.0);
                CHECK(w.rhohat > 0.0);
                CHECK(w.rhohat < 1.0);
            }
}

TEST_CASE("an unstable queue reports infinity, not a division") {
    const qs::QsysResult<double> u = qs::qsys_gig1_approx_whitt(2.0, 1.0, 1.0, 1.0);
    CHECK(std::isinf(u.W));
    CHECK(u.rhohat == doctest::Approx(1.0));
    const qs::QsysResult<double> e = qs::qsys_gig1_approx_whitt(1.0, 1.0, 1.0, 1.0);
    CHECK(std::isinf(e.W));  // rho exactly one is unstable too
}
