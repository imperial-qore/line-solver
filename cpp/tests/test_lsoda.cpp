/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

/**
 * LSODA, the integrator the fluid solver runs on.
 *
 * THE EXPECTED VALUES ARE NOT THIS PORT'S OWN OUTPUT. They are the C reference
 * figures already carried by `python/tests/test_lsoda.py`, produced by
 * compiling and running the C LSODA at full precision, and the native-Python
 * port is held to the same numbers. Checking the C++ against them is what
 * makes "all three codebases run one algorithm" a measured claim rather than a
 * provenance argument -- see the lineage note at the top of
 * `third_party/lsoda.hpp`.
 *
 * The tolerances below are the ones the Python suite uses. They are looser
 * than the integration tolerance on purpose: two ports of an ADAPTIVE solver
 * take the same steps only until the first difference in floating-point
 * evaluation order, after which the step sequences drift apart while both stay
 * within the requested error. Demanding bit equality would be testing the
 * compiler, not the algorithm.
 */

#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/solvers/fluid/fluid_stiff.h"
#include "line/util/lsoda.h"

using namespace line;

namespace {

/** Relative check with an absolute floor, as the Python suite applies it. */
void check_close(double got, double want, double rtol, double atol = 0.0) {
    const double tol = atol + rtol * std::fabs(want);
    CHECK(std::fabs(got - want) <= tol);
}

TEST_CASE("LSODA integrates Robertson to the C reference values") {
    // The classic stiff benchmark, over twelve decades of time. atol is per
    // component: the middle species is ~1e-8 at its peak and would be lost in
    // the noise of a scalar 1e-6.
    const LsodaRhs f = [](double, const double* y, double* d) {
        d[0] = 1.0e4 * y[1] * y[2] - 0.04 * y[0];
        d[2] = 3.0e7 * y[1] * y[1];
        d[1] = -(d[0] + d[2]);
    };
    const std::vector<double> y0{1.0, 0.0, 0.0};

    std::vector<double> t_eval{0.0};
    for (int i = 0; i < 12; ++i) t_eval.push_back(0.4 * std::pow(10.0, i));

    LsodaOptions opt;
    opt.rtol = 1e-4;
    opt.atol_vec = std::vector<double>{1e-6, 1e-10, 1e-6};

    const LsodaSolution s = lsoda_integrate(f, y0, t_eval, opt);
    REQUIRE(s.success);
    REQUIRE(s.y.size() == t_eval.size());

    static const double expected[12][3] = {
        {9.851712049052344e-01, 3.386380058450371e-05, 1.479493129418126e-02},
        {9.055332934421416e-01, 2.240655242773938e-05, 9.444430000543071e-02},
        {7.158403474138746e-01, 9.186333848221152e-06, 2.841504662522774e-01},
        {4.505250299518244e-01, 3.222963773719375e-06, 5.494717470844018e-01},
        {1.831975696546021e-01, 8.941773261377527e-07, 8.168015361680721e-01},
        {3.898729244522436e-02, 1.621939581471608e-07, 9.610125453607797e-01},
        {4.936362159709025e-03, 1.984220853834811e-08, 9.950636179980765e-01},
        {5.161832921225141e-04, 2.065786542448266e-09, 9.994838146420849e-01},
        {5.179803644367760e-05, 2.072026913293681e-10, 9.999482017563474e-01},
        {5.283686292063946e-06, 2.113485306887202e-11, 9.999947162925692e-01},
        {4.658646268319763e-07, 1.863459326277758e-12, 9.999995341335010e-01},
        {1.431681921554496e-08, 5.726732909838437e-14, 9.999999856831133e-01}};

    for (std::size_t i = 0; i < 12; ++i)
        for (std::size_t k = 0; k < 3; ++k)
            check_close(s.y[i + 1][k], expected[i][k], 5e-3, 1e-8);

    // Mass is conserved exactly in the model, so it is a check the reference
    // values cannot influence: the three species must always sum to one.
    for (const std::vector<double>& yi : s.y)
        check_close(yi[0] + yi[1] + yi[2], 1.0, 1e-6);

    // Robertson is stiff, so LSODA must have switched away from Adams.
    CHECK(s.method == "bdf");
}

TEST_CASE("LSODA integrates HIRES to the C reference values") {
    // Eight-equation photochemistry from Hairer and Wanner, moderately stiff.
    const LsodaRhs f = [](double, const double* y, double* d) {
        d[0] = -1.71 * y[0] + 0.43 * y[1] + 8.32 * y[2] + 0.0007;
        d[1] = 1.71 * y[0] - 8.75 * y[1];
        d[2] = -10.03 * y[2] + 0.43 * y[3] + 0.035 * y[4];
        d[3] = 8.32 * y[1] + 1.71 * y[2] - 1.12 * y[3];
        d[4] = -1.745 * y[4] + 0.43 * y[5] + 0.43 * y[6];
        d[5] = -280.0 * y[5] * y[7] + 0.69 * y[3] + 1.71 * y[4] - 0.43 * y[5] + 0.69 * y[6];
        d[6] = 280.0 * y[5] * y[7] - 1.81 * y[6];
        d[7] = -280.0 * y[5] * y[7] + 1.81 * y[6];
    };
    const std::vector<double> y0{1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0057};
    const std::vector<double> t_eval{0.0, 100.0, 200.0, 321.8122};

    LsodaOptions opt;
    opt.rtol = 1e-6;
    opt.atol = 1e-8;

    const LsodaSolution s = lsoda_integrate(f, y0, t_eval, opt);
    REQUIRE(s.success);

    static const double expected[3][8] = {
        {4.520862869231245e-03, 8.839063180882597e-04, 7.971949677087163e-04,
         7.811332285106224e-03, 1.323853727051267e-01, 5.301681739823277e-01,
         5.631339591581244e-03, 6.866040841874939e-05},
        {2.736525896095168e-03, 5.351908569598239e-04, 4.485119404720507e-04,
         4.688161408919465e-03, 7.083443534644768e-02, 2.804641383048971e-01,
         5.571599220744725e-03, 1.284007792552683e-04},
        {7.371423243758239e-04, 1.442507544489246e-04, 5.888935680491721e-05,
         1.175672019019430e-03, 2.386687478330434e-03, 6.239970967822740e-03,
         2.850266805796065e-03, 2.849733194203930e-03}};

    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t k = 0; k < 8; ++k) check_close(s.y[i + 1][k], expected[i][k], 1e-4, 1e-12);
}

TEST_CASE("LSODA integrates the stiff Van der Pol oscillator inside its limit cycle") {
    // mu = 1000 is the standard very-stiff setting. There is no published
    // reference vector, so the invariant is the one the physics guarantees:
    // the orbit is attracted to a limit cycle with |y1| just above 2.
    const double mu = 1000.0;
    const LsodaRhs f = [mu](double, const double* y, double* d) {
        d[0] = y[1];
        d[1] = mu * (1.0 - y[0] * y[0]) * y[1] - y[0];
    };
    const std::vector<double> y0{2.0, 0.0};

    std::vector<double> t_eval;
    for (int i = 0; i <= 100; ++i) t_eval.push_back(3000.0 * i / 100.0);

    LsodaOptions opt;
    opt.rtol = 1e-6;
    opt.atol = 1e-8;

    const LsodaSolution s = lsoda_integrate(f, y0, t_eval, opt);
    REQUIRE(s.success);

    double ymax = 0.0;
    for (const std::vector<double>& yi : s.y) ymax = std::max(ymax, std::fabs(yi[0]));
    CHECK(ymax <= 2.5);
    CHECK(ymax > 1.5);
    CHECK(s.method == "bdf");  // a non-stiff method cannot hold this problem
}

TEST_CASE("LSODA reproduces an exactly solvable decay") {
    // y' = -2y has the closed form y = y0 exp(-2t); no reference file needed
    // and it pins the sign and scaling conventions of the wrapper itself.
    const LsodaRhs f = [](double, const double* y, double* d) { d[0] = -2.0 * y[0]; };
    const std::vector<double> y0{3.0};
    const std::vector<double> t_eval{0.0, 0.5, 1.0, 2.0, 5.0};

    LsodaOptions opt;
    opt.rtol = 1e-10;
    opt.atol = 1e-12;

    const LsodaSolution s = lsoda_integrate(f, y0, t_eval, opt);
    REQUIRE(s.success);
    for (std::size_t i = 0; i < t_eval.size(); ++i)
        check_close(s.y[i][0], 3.0 * std::exp(-2.0 * t_eval[i]), 1e-8, 1e-12);

    // A linear non-stiff problem: LSODA should stay on Adams throughout.
    CHECK(s.method == "adams");
    CHECK(s.final_time() == doctest::Approx(5.0));
}

TEST_CASE("force_stiff pins the BDF half and keeps HIRES accurate") {
    // The pin exists because the Adams half loses its stability bound at a fixed
    // point (see LsodaOptions::force_stiff); the JAR carries the same fork as
    // LSODA.setForceStiff(true), MATLAB and native Python as forceStiff. There
    // is no C reference for the pinned solver -- the C cannot pin -- so it is
    // checked against the reference row the auto-switcher reproduces exactly.
    const LsodaRhs f = [](double, const double* y, double* d) {
        d[0] = -1.71 * y[0] + 0.43 * y[1] + 8.32 * y[2] + 0.0007;
        d[1] = 1.71 * y[0] - 8.75 * y[1];
        d[2] = -10.03 * y[2] + 0.43 * y[3] + 0.035 * y[4];
        d[3] = 8.32 * y[1] + 1.71 * y[2] - 1.12 * y[3];
        d[4] = -1.745 * y[4] + 0.43 * y[5] + 0.43 * y[6];
        d[5] = -280.0 * y[5] * y[7] + 0.69 * y[3] + 1.71 * y[4] - 0.43 * y[5] + 0.69 * y[6];
        d[6] = 280.0 * y[5] * y[7] - 1.81 * y[6];
        d[7] = -280.0 * y[5] * y[7] + 1.81 * y[6];
    };
    const std::vector<double> y0{1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0057};
    const std::vector<double> t_eval{0.0, 100.0, 200.0, 321.8122};
    static const double expected[8] = {7.371423243758239e-04, 1.442507544489246e-04,
                                       5.888935680491721e-05, 1.175672019019430e-03,
                                       2.386687478330434e-03, 6.239970967822740e-03,
                                       2.850266805796065e-03, 2.849733194203930e-03};

    LsodaOptions opt;
    opt.rtol = 1e-6;
    opt.atol = 1e-8;
    opt.force_stiff = true;

    const LsodaSolution s = lsoda_integrate(f, y0, t_eval, opt);
    REQUIRE(s.success);
    // Never on Adams, not even for the first step, which is what "pinned" means.
    CHECK(s.method == "bdf");
    for (std::size_t k = 0; k < 8; ++k) check_close(s.y[3][k], expected[k], 1e-4, 1e-12);
}

TEST_CASE("force_stiff never leaves BDF on a problem LSODA would call nonstiff") {
    // dy/dt = -y is where the auto-switcher stays on Adams for the whole solve
    // (the decay case above asserts exactly that), so this separates the pin
    // from the detector rather than agreeing with it by accident.
    const LsodaRhs f = [](double, const double* y, double* d) { d[0] = -y[0]; };
    const std::vector<double> y0{1.0};
    const std::vector<double> t_eval{0.0, 1.0, 2.0, 5.0};

    LsodaOptions opt;
    opt.rtol = 1e-8;
    opt.atol = 1e-10;
    opt.force_stiff = true;

    const LsodaSolution s = lsoda_integrate(f, y0, t_eval, opt);
    REQUIRE(s.success);
    CHECK(s.method == "bdf");
    for (std::size_t i = 1; i < t_eval.size(); ++i)
        check_close(s.y[i][0], std::exp(-t_eval[i]), 1e-6);
}

TEST_CASE("force_stiff reaches the one-step stepper too") {
    // LsodaStepper is the arm the NonNegative clip-and-reset rule drives, so the
    // pin has to be in force there as well and not only in lsoda_integrate.
    const LsodaRhs f = [](double, const double* y, double* d) { d[0] = -1000.0 * (y[0] - 1.0); };
    const std::vector<double> y0{0.0};

    LsodaOptions opt;
    opt.rtol = 1e-8;
    opt.atol = 1e-10;
    opt.force_stiff = true;

    LsodaStepper stepper(f, y0, 0.0, 1.0, opt);
    std::size_t taken = 0;
    while (stepper.step() && taken < 10000) ++taken;
    REQUIRE_FALSE(stepper.failed());
    stepper.settle_at_end();
    CHECK(taken > 0);
    check_close(stepper.y()[0], 1.0, 1e-6);
}

TEST_CASE("force_stiff survives a cold start with a component pinned at zero") {
    // THE FAILURE THIS ASSERTS AGAINST IS REAL AND WAS MEASURED: pinning BDF
    // takes a finite-difference Jacobian at t0, where the auto-switcher would
    // still be on Adams and take none, and ODEPACK's increment
    // max(sqrt(eps)*|y_j|, r0/ewt_j) collapses to ~1e-19 for y2 = 0 under atol
    // 1e-10. The column is then noise, the corrector converges to nonsense, and
    // Robertson ran away to y1 = -1.5e7 while reporting success -- flipping on a
    // last-bit change to the right-hand side. prja's numjac floor now applies on
    // every step under the pin, which is what ode15s does.
    const LsodaRhs f = [](double, const double* y, double* d) {
        d[0] = 1.0e4 * y[1] * y[2] - 0.04 * y[0];
        d[2] = 3.0e7 * y[1] * y[1];
        d[1] = -(d[0] + d[2]);
    };
    const std::vector<double> y0{1.0, 0.0, 0.0};
    std::vector<double> t_eval{0.0};
    for (int i = 0; i < 12; ++i) t_eval.push_back(0.4 * std::pow(10.0, i));

    LsodaOptions opt;
    opt.rtol = 1e-4;
    opt.rtol_vec = std::vector<double>(3, 1e-4);
    opt.atol_vec = std::vector<double>{1e-6, 1e-10, 1e-6};
    opt.force_stiff = true;

    const LsodaSolution s = lsoda_integrate(f, y0, t_eval, opt);
    REQUIRE(s.success);
    CHECK(s.method == "bdf");
    const std::vector<double>& yend = s.final_state();
    // y1 -> 1.43e-8 and y3 -> 1, both well inside their own atol of 1e-6, and
    // the mass y1+y2+y3 = 1 the system conserves exactly.
    CHECK(std::fabs(yend[0]) < 1e-6);
    CHECK(std::fabs(yend[2] - 1.0) < 1e-6);
    CHECK(std::fabs(yend[0] + yend[1] + yend[2] - 1.0) < 1e-8);
}

TEST_CASE("force_stiff reaches the fluid integration arm") {
    // In-solver coverage: fluid_integrate_leg is what solver_fluid drives, and
    // it takes the same LsodaOptions, so the pin has to arrive through it and
    // not only through the two wrapper entry points. The drift is a two-stage
    // relaxation with a fast mode, i.e. the shape of a fluid leg, and its fixed
    // point (1, 1) is known in closed form.
    const std::function<void(double, const double*, double*)> drift =
        [](double, const double* x, double* dx) {
            dx[0] = -1000.0 * (x[0] - 1.0);
            dx[1] = x[0] - x[1];
        };
    const std::vector<double> x0{0.0, 0.0};

    LsodaOptions lopt;
    lopt.rtol = 1e-8;
    lopt.atol = 1e-10;
    LsodaOptions pinned = lopt;
    pinned.force_stiff = true;

    const std::vector<double> autoswitch = line::fluid::fluid_integrate_leg(drift, 0.0, 20.0, x0, lopt);
    const std::vector<double> forced = line::fluid::fluid_integrate_leg(drift, 0.0, 20.0, x0, pinned);
    for (std::size_t k = 0; k < 2; ++k) {
        check_close(forced[k], 1.0, 1e-6);
        check_close(forced[k], autoswitch[k], 1e-6);
    }
}

TEST_CASE("lsoda_final returns the end state and the wrapper refuses bad input") {
    const LsodaRhs f = [](double, const double* y, double* d) { d[0] = -y[0]; };
    const std::vector<double> y0{1.0};

    const std::vector<double> yT = lsoda_final(f, y0, 0.0, 1.0);
    check_close(yT[0], std::exp(-1.0), 1e-5);

    // Refused BY NAME rather than silently reinterpreted.
    CHECK_THROWS_AS(lsoda_integrate(f, y0, std::vector<double>{}), InputError);
    CHECK_THROWS_AS(lsoda_integrate(f, std::vector<double>{}, std::vector<double>{0.0, 1.0}),
                    InputError);
    CHECK_THROWS_AS(lsoda_integrate(f, y0, std::vector<double>{0.0, 1.0, 0.5}), InputError);

    LsodaOptions bad;
    bad.atol_vec = std::vector<double>{1e-8, 1e-8};  // two entries for one equation
    CHECK_THROWS_AS(lsoda_integrate(f, y0, std::vector<double>{0.0, 1.0}, bad), InputError);
}

}  // namespace
