/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * G/G/1 dispatch and the remaining G/I/G/1 approximations: qsys_gg1, the two
 * Myskja third-moment forms, and the Whitt-You robust-queueing bound.
 *
 * Oracles: the exact branches of qsys_gg1 must reproduce qsys_mm1, qsys_mg1
 * and qsys_gm1 bit for bit; Myskja is exact for M/G/1 by construction, so at
 * ca = 1 it must equal the Pollaczek-Khinchine value; the robust-queueing
 * workload must satisfy its own variational inequality at every trial point.
 * MATLAB reference values come from running matlab/src/api/qsys through
 *   matlab -singleCompThread -batch "addpath(genpath('matlab/src')); ..."
 */
#include <cmath>

#include "doctest.h"
#include "line/api/qsys/qsys_gg1.h"
#include "line/api/qsys/qsys_gig1_approx_allencunneen.h"
#include "line/api/qsys/qsys_gig1_approx_myskja.h"
#include "line/api/qsys/qsys_gig1_approx_myskja2.h"
#include "line/api/qsys/qsys_gig1_rq.h"
#include "line/api/qsys/qsys_mg1.h"
#include "line/api/qsys/qsys_mm1.h"

using line::Real50;

TEST_CASE("qsys_gg1 reproduces its exact branches") {
    // ca2 = cs2 = 1 must be M/M/1 exactly.
    auto mm = line::qsys::qsys_gg1(0.6, 1.0, 1.0, 1.0);
    auto mm1 = line::qsys::qsys_mm1(0.6, 1.0);
    CHECK(mm.W == mm1.W);
    CHECK(mm.rhohat == mm1.rhohat);
    // ca2 = 1 must be M/G/1 with cs = sqrt(cs2), exactly.
    auto mg = line::qsys::qsys_gg1(0.6, 1.0, 1.0, 4.0);
    auto mg1 = line::qsys::qsys_mg1(0.6, 1.0, 2.0);
    CHECK(mg.W == mg1.W);
    CHECK(mg.rhohat == mg1.rhohat);
    // MATLAB qsys_gg1(0.6,1,1,4) -> W 4.75, rhohat 0.74025974025974028.
    CHECK(mg.W == doctest::Approx(4.75).epsilon(1e-15));
    CHECK(mg.rhohat == doctest::Approx(0.74025974025974028).epsilon(1e-15));
}

TEST_CASE("qsys_gg1 general branch is Allen-Cunneen") {
    auto gen = line::qsys::qsys_gg1(0.6, 1.0, 2.0, 3.0);
    auto ac = line::qsys::qsys_gig1_approx_allencunneen(0.6, 1.0, std::sqrt(2.0), std::sqrt(3.0));
    CHECK(gen.W == ac.W);
    // MATLAB qsys_gg1(0.6,1,2,3) -> W 4.75, rhohat 0.74025974025974028.
    CHECK(gen.W == doctest::Approx(4.75).epsilon(1e-14));
    CHECK(gen.rhohat == doctest::Approx(0.74025974025974028).epsilon(1e-14));
}

TEST_CASE("qsys_gg1 G/M/1 branch matches MATLAB on both sides of ca2 = 1") {
    // The fixed point stops at an absolute step of 1e-13, so 1e-11 relative is
    // the tolerance claimed against MATLAB's identical iteration.
    auto h2 = line::qsys::qsys_gg1(0.6, 1.0, 2.0, 1.0);  // H2 interarrival fit
    CHECK(h2.W == doctest::Approx(3.2019410160095809).epsilon(1e-11));
    CHECK(h2.rhohat == doctest::Approx(0.6576707807866522).epsilon(1e-11));
    auto er = line::qsys::qsys_gg1(0.6, 1.0, 0.5, 1.0);  // Erlang mixture fit
    CHECK(er.W == doctest::Approx(1.9834994352913979).epsilon(1e-11));
    CHECK(er.rhohat == doctest::Approx(0.54339977411641249).epsilon(1e-11));
    // Less arrival variability, less delay.
    CHECK(er.W < h2.W);
    // And both bracket the M/M/1 value at ca2 = 1.
    CHECK(er.W < line::qsys::qsys_mm1(0.6, 1.0).W);
    CHECK(h2.W > line::qsys::qsys_mm1(0.6, 1.0).W);
}

TEST_CASE("Myskja is exact for M/G/1 at ca = 1") {
    // The correction term carries the factor (ca^2-1), so at ca = 1 the
    // formula is the Pollaczek-Khinchine mean plus the service time.
    const double lambda = 0.6, mu = 1.0, cs = 0.8;
    auto m = line::qsys::qsys_gig1_approx_myskja(lambda, mu, 1.0, cs, 0.5, 0.9);
    auto mg1 = line::qsys::qsys_mg1(lambda, mu, cs);
    CHECK(m.W == doctest::Approx(mg1.W).epsilon(1e-14));
    // Myskja2 branches to qsys_mg1 outright at ca = 1.
    auto m2 = line::qsys::qsys_gig1_approx_myskja2(lambda, mu, 1.0, cs, 0.5, 0.9);
    CHECK(m2.W == mg1.W);
    CHECK(m2.rhohat == mg1.rhohat);
}

TEST_CASE("Myskja matches MATLAB") {
    // MATLAB qsys_gig1_approx_myskja(0.6,1,1.5,0.8,0.5,0.9)
    //   -> W 3.0646979634944174, rhohat 0.64774081117271165.
    auto m = line::qsys::qsys_gig1_approx_myskja(0.6, 1.0, 1.5, 0.8, 0.5, 0.9);
    CHECK(m.W == doctest::Approx(3.0646979634944174).epsilon(1e-14));
    CHECK(m.rhohat == doctest::Approx(0.64774081117271165).epsilon(1e-14));
    // MATLAB qsys_gig1_approx_myskja2(0.6,1,1.5,0.8,0.5,0.9)
    //   -> W 10.756604689105092, rhohat 0.86584317295038826.
    auto m2 = line::qsys::qsys_gig1_approx_myskja2(0.6, 1.0, 1.5, 0.8, 0.5, 0.9);
    CHECK(m2.W == doctest::Approx(10.756604689105092).epsilon(1e-14));
    CHECK(m2.rhohat == doctest::Approx(0.86584317295038826).epsilon(1e-14));
}

TEST_CASE("Myskja reads ca and cs as coefficients of variation, not SCVs") {
    // The JAR writes (1+cs) and (ca-1) where MATLAB writes (1+cs^2) and
    // (ca^2-1). This pins the MATLAB convention: feeding the squares gives a
    // different answer, so a port that silently adopted the JAR reading would
    // fail here.
    auto cv = line::qsys::qsys_gig1_approx_myskja(0.6, 1.0, 1.5, 0.8, 0.5, 0.9);
    auto scv = line::qsys::qsys_gig1_approx_myskja(0.6, 1.0, 1.5 * 1.5, 0.8 * 0.8, 0.5, 0.9);
    CHECK(cv.W != scv.W);
    CHECK(cv.W == doctest::Approx(3.0646979634944174).epsilon(1e-14));
}

TEST_CASE("Myskja agrees between double and Real50") {
    auto d = line::qsys::qsys_gig1_approx_myskja(0.6, 1.0, 1.5, 0.8, 0.5, 0.9);
    auto h = line::qsys::qsys_gig1_approx_myskja(Real50(6) / Real50(10), Real50(1),
                                                 Real50(15) / Real50(10), Real50(8) / Real50(10),
                                                 Real50(5) / Real50(10), Real50(9) / Real50(10));
    CHECK(static_cast<double>(h.W) == doctest::Approx(d.W).epsilon(1e-13));
}

TEST_CASE("robust queueing reproduces the exact M/M/1 workload") {
    // With I_a = c_s^2 = 1 the supremum is attained in closed form and gives
    // Z = rho/(1-rho)/mu, hence W = Z/rho - 1/mu = the M/M/1 waiting time.
    // MATLAB qsys_gig1_rq(0.5,2,1,@(x)1) -> Z 0.5, W 0.5, Q 0.5, X 1.
    auto r = line::qsys::qsys_gig1_rq<double>(0.5, 2.0, 1.0, [](double) { return 1.0; });
    CHECK(r.Z == doctest::Approx(0.5).epsilon(1e-9));
    CHECK(r.W == doctest::Approx(0.5).epsilon(1e-8));
    CHECK(r.Q == doctest::Approx(0.5).epsilon(1e-8));
    CHECK(r.X == doctest::Approx(1.0).epsilon(1e-8));
    // The M/M/1 waiting time at lambda = 1, mu = 2 is rho/(mu(1-rho)) = 0.5.
    CHECK(r.W == doctest::Approx(line::qsys::qsys_mm1(1.0, 2.0).W - 0.5).epsilon(1e-8));
}

TEST_CASE("robust queueing matches MATLAB on a variable instance") {
    // MATLAB qsys_gig1_rq(0.7,1,2,@(x)2) -> Z 4.6666666666666652,
    // W 5.1666666666666652, Q 3.6166666666666654, X 4.3166666666666655.
    // The optimizer is driven to TolX 1e-10 on the same bracket, so 1e-9
    // relative is the tolerance claimed.
    auto r = line::qsys::qsys_gig1_rq<double>(0.7, 1.0, 2.0, [](double) { return 2.0; });
    CHECK(r.Z == doctest::Approx(4.6666666666666652).epsilon(1e-9));
    CHECK(r.W == doctest::Approx(5.1666666666666652).epsilon(1e-9));
    CHECK(r.Q == doctest::Approx(3.6166666666666654).epsilon(1e-9));
    CHECK(r.X == doctest::Approx(4.3166666666666655).epsilon(1e-9));
    CHECK(r.X == doctest::Approx(r.Q + 0.7).epsilon(1e-14));
}

TEST_CASE("robust queueing returns a genuine supremum") {
    // Z must dominate the objective at every trial point; that is the defining
    // property and does not depend on the optimizer at all.
    const double rho = 0.8, mu = 1.5, cs2 = 3.0;
    auto ia = [](double) { return 4.0; };
    auto r = line::qsys::qsys_gig1_rq<double>(rho, mu, cs2, ia);
    for (double x = 1e-3; x < 1e5; x *= 2.0) {
        const double f = -(1.0 - rho) * x + std::sqrt(2.0 * rho * x * (ia(x) + cs2) / mu);
        CHECK(r.Z >= f - 1e-9);
    }
    CHECK(r.Z > 0.0);
    // Zero load gives zero everywhere, and an unstable load is rejected.
    auto z = line::qsys::qsys_gig1_rq<double>(0.0, 1.0, 1.0, ia);
    CHECK(z.Z == 0.0);
    CHECK(z.X == 0.0);
    CHECK_THROWS_AS(line::qsys::qsys_gig1_rq<double>(1.0, 1.0, 1.0, ia), line::InputError);
}
