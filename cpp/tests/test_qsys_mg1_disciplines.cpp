/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * M/G/1 scheduling disciplines: HOL priorities, FB/LAS, SETF, LRPT, PSJF,
 * SRPT, and multiclass DPS.
 *
 * On the exponential-service reduction. FB is work conserving and blind to the
 * job size, so for exponential service its mean response time is the M/M/1
 * value; that identity is asserted. The other disciplines here are size or
 * age based and deliberately are NOT M/M/1: PSJF and SRPT beat it, LRPT and
 * SETF lose to it. Asserting M/M/1 for those would assert that the discipline
 * does nothing. What is asserted instead is the ordering that must hold for
 * any work-conserving M/G/1,
 *
 *     SRPT <= PSJF <= FB = M/M/1 <= LRPT,      SETF >= FB,
 *
 * together with MATLAB reference values for each.
 *
 * 2026-07-29: mva_feature_set was narrowed to MATLAB's 59 declared names, so a
 * size-based model no longer reaches these closed forms through SolverMVA in
 * either codebase. Every case here calls line::qsys::qsys_mg1_* directly, which
 * is an api-level entry point and is unaffected; the closed forms remain
 * supported, only the SolverMVA route to them is gone.
 *
 * References were obtained by running matlab/src/api/qsys through
 *   matlab -singleCompThread -batch "addpath(genpath('matlab/src')); ..."
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/qsys/qsys_mg1_fb.h"
#include "line/api/qsys/qsys_mg1_lrpt.h"
#include "line/api/qsys/qsys_mg1_prio.h"
#include "line/api/qsys/qsys_mg1_psjf.h"
#include "line/api/qsys/qsys_mg1_setf.h"
#include "line/api/qsys/qsys_mg1_srpt.h"
#include "line/api/qsys/qsys_mm1.h"
#include "line/api/qsys/qsys_mm1_dps.h"

using line::Rational;

namespace {
std::vector<double> v3(double a, double b) { return std::vector<double>{a, b}; }
}  // namespace

TEST_CASE("M/G/1 HOL priority is exact and matches MATLAB") {
    // MATLAB qsys_mg1_prio([0.2;0.3],[1;1],[0.5;1]) -> W = [1.53125; 2.0625],
    // rhohat = 0.48051948051948051. Rational arithmetic, so 1e-14 relative.
    std::vector<Rational> lam = {Rational(1, 5), Rational(3, 10)};
    std::vector<Rational> mu = {Rational(1), Rational(1)};
    std::vector<Rational> cs = {Rational(1, 2), Rational(1)};
    auto r = line::qsys::qsys_mg1_prio(lam, mu, cs);
    CHECK(r.W[0] == Rational(49, 32));  // 1.53125 exactly
    CHECK(static_cast<double>(r.W[0]) == doctest::Approx(1.53125).epsilon(1e-15));
    CHECK(static_cast<double>(r.W[1]) == doctest::Approx(2.0625).epsilon(1e-15));
    CHECK(static_cast<double>(r.rhohat) == doctest::Approx(0.48051948051948051).epsilon(1e-14));
    // Higher priority is served first, so its response time is smaller.
    CHECK(r.W[0] < r.W[1]);
}

TEST_CASE("M/G/1 priority with one exponential class is M/M/1") {
    // A single class with cs = 1 makes the Cobham formula the P-K formula.
    std::vector<Rational> lam = {Rational(1, 2)};
    std::vector<Rational> mu = {Rational(1)};
    std::vector<Rational> cs = {Rational(1)};
    auto r = line::qsys::qsys_mg1_prio(lam, mu, cs);
    auto mm1 = line::qsys::qsys_mm1(Rational(1, 2), Rational(1));
    CHECK(r.W[0] == mm1.W);  // exactly 2
    CHECK(r.W[0] == Rational(2));
    CHECK_THROWS_AS(line::qsys::qsys_mg1_prio(std::vector<Rational>{Rational(2)},
                                              std::vector<Rational>{Rational(1)},
                                              std::vector<Rational>{Rational(1)}),
                    line::InputError);  // unstable
}

TEST_CASE("FB with exponential service is exactly the M/M/1 mean response time") {
    // This is the discipline-independence of the mean under a size-blind,
    // work-conserving policy. MATLAB itself lands on 1.9999999093092409
    // because it truncates the outer integral at 20 mean service times; the
    // port reproduces that same truncation, so both sit 4.5e-8 below 2.
    auto r = line::qsys::qsys_mg1_fb(std::vector<double>{0.5}, std::vector<double>{1.0},
                                     std::vector<double>{1.0});
    auto mm1 = line::qsys::qsys_mm1(0.5, 1.0);
    CHECK(r.W[0] == doctest::Approx(mm1.W).epsilon(1e-7));
    CHECK(r.W[0] == doctest::Approx(1.9999999093092409).epsilon(1e-9));
}

TEST_CASE("the size-based disciplines order as SRPT <= PSJF <= FB <= LRPT") {
    const std::vector<double> lam{0.5}, mu{1.0}, cs{1.0};
    const double fb = line::qsys::qsys_mg1_fb(lam, mu, cs).W[0];
    const double psjf = line::qsys::qsys_mg1_psjf(lam, mu, cs).W[0];
    const double srpt = line::qsys::qsys_mg1_srpt(lam, mu, cs).W[0];
    const double lrpt = line::qsys::qsys_mg1_lrpt(lam, mu, cs).W[0];
    const double setf = line::qsys::qsys_mg1_setf(lam, mu, cs).W[0];
    CHECK(srpt <= psjf);
    CHECK(psjf <= fb);
    CHECK(fb <= lrpt);
    CHECK(setf >= fb);
    // MATLAB reference values for the same instance.
    CHECK(psjf == doctest::Approx(1.5313831926055728).epsilon(1e-8));
    CHECK(srpt == doctest::Approx(1.4253720785260464).epsilon(1e-5));
    CHECK(lrpt == doctest::Approx(3.9999999093092402).epsilon(1e-9));
    CHECK(setf == doctest::Approx(3.2066798723311383).epsilon(1e-13));
}

TEST_CASE("two-class exponential disciplines match MATLAB") {
    // lambda = [0.2 0.3], mu = [1 2], cs = [1 1]; rho = 0.35.
    const std::vector<double> lam = v3(0.2, 0.3), mu = v3(1.0, 2.0), cs = v3(1.0, 1.0);

    auto setf = line::qsys::qsys_mg1_setf(lam, mu, cs);
    CHECK(setf.W[0] == doctest::Approx(2.461340796845878).epsilon(1e-13));
    CHECK(setf.W[1] == doctest::Approx(1.5073502149086924).epsilon(1e-13));
    CHECK(setf.rhohat == doctest::Approx(0.48572189745855437).epsilon(1e-13));

    // FB and PSJF go through the quadrature; MATLAB asks integral() for
    // RelTol 1e-8, so that is the tolerance claimed.
    auto fb = line::qsys::qsys_mg1_fb(lam, mu, cs);
    CHECK(fb.W[0] == doctest::Approx(1.5732010984766482).epsilon(1e-8));
    CHECK(fb.W[1] == doctest::Approx(0.72291123066388818).epsilon(1e-8));
    CHECK(fb.rhohat == doctest::Approx(0.34705117391623186).epsilon(1e-8));

    auto lrpt = line::qsys::qsys_mg1_lrpt(lam, mu, cs);
    CHECK(lrpt.W[0] == doctest::Approx(2.1893490444933392).epsilon(1e-9));
    CHECK(lrpt.W[1] == doctest::Approx(1.4201183085581286).epsilon(1e-9));
    CHECK(lrpt.rhohat == doctest::Approx(0.46349205659030945).epsilon(1e-9));

    auto psjf = line::qsys::qsys_mg1_psjf(lam, mu, cs);
    CHECK(psjf.W[0] == doctest::Approx(1.3707251721087317).epsilon(1e-8));
    CHECK(psjf.W[1] == doctest::Approx(0.60363165980002575).epsilon(1e-8));
    CHECK(psjf.rhohat == doctest::Approx(0.31282554271368002).epsilon(1e-8));

    // SRPT is a fixed trapezoid grid in the reference, accurate to about 1e-5
    // relative; the port uses the same grid, so the two agree far more closely
    // than either agrees with the true value.
    auto srpt = line::qsys::qsys_mg1_srpt(lam, mu, cs);
    CHECK(srpt.W[0] == doctest::Approx(1.2867477903534605).epsilon(1e-5));
    CHECK(srpt.W[1] == doctest::Approx(0.59194260398252196).epsilon(1e-5));
}

TEST_CASE("the non-exponential surrogate paths match MATLAB") {
    // cs = [2 0.5] takes every discipline off its exponential branch.
    const std::vector<double> lam = v3(0.2, 0.3), mu = v3(1.0, 2.0), cs = v3(2.0, 0.5);

    auto setf = line::qsys::qsys_mg1_setf(lam, mu, cs);
    CHECK(setf.W[0] == doctest::Approx(3.7795857988165675).epsilon(1e-13));
    CHECK(setf.W[1] == doctest::Approx(2.3611111111111107).epsilon(1e-13));

    auto fb = line::qsys::qsys_mg1_fb(lam, mu, cs);
    CHECK(fb.W[0] == doctest::Approx(1.9526627218934909).epsilon(1e-13));
    CHECK(fb.W[1] == doctest::Approx(0.77777777777777768).epsilon(1e-13));

    auto lrpt = line::qsys::qsys_mg1_lrpt(lam, mu, cs);
    CHECK(lrpt.W[0] == doctest::Approx(1.25).epsilon(1e-13));
    CHECK(lrpt.W[1] == doctest::Approx(1.0288461538461537).epsilon(1e-13));

    auto psjf = line::qsys::qsys_mg1_psjf(lam, mu, cs);
    CHECK(psjf.W[0] == doctest::Approx(2.8328402366863905).epsilon(1e-13));
    CHECK(psjf.W[1] == doctest::Approx(0.65311418685121114).epsilon(1e-13));
}

TEST_CASE("the discipline family rejects unstable and malformed input") {
    const std::vector<double> bad = v3(2.0, 2.0), mu = v3(1.0, 1.0), cs = v3(1.0, 1.0);
    CHECK_THROWS_AS(line::qsys::qsys_mg1_fb(bad, mu, cs), line::InputError);
    CHECK_THROWS_AS(line::qsys::qsys_mg1_setf(bad, mu, cs), line::InputError);
    CHECK_THROWS_AS(line::qsys::qsys_mg1_lrpt(bad, mu, cs), line::InputError);
    CHECK_THROWS_AS(line::qsys::qsys_mg1_psjf(bad, mu, cs), line::InputError);
    CHECK_THROWS_AS(line::qsys::qsys_mg1_srpt(bad, mu, cs), line::InputError);
    const std::vector<double> one{1.0};
    CHECK_THROWS_AS(line::qsys::qsys_mg1_fb(one, mu, cs), line::InputError);  // length mismatch
}

TEST_CASE("M/M/1 DPS with equal weights and rates is M/M/1 per class") {
    // Two symmetric classes, lambda = 0.25 each, mu = 1: the total is M/M/1 at
    // rho = 0.5, and by symmetry each class sees W = 2. The truncation at
    // level 40 leaves a residual of order 1e-10, which is the assertion.
    auto r = line::qsys::qsys_mm1_dps(v3(0.25, 0.25), v3(1.0, 1.0), v3(1.0, 1.0), 1e-10, 40u);
    CHECK(r.rho == doctest::Approx(0.5).epsilon(1e-15));
    CHECK(r.T_[0] == doctest::Approx(2.0).epsilon(1e-9));
    CHECK(r.T_[1] == doctest::Approx(2.0).epsilon(1e-9));
    CHECK(r.T_[0] == doctest::Approx(r.T_[1]).epsilon(1e-12));
    // MATLAB qsys_mm1_dps([0.25 0.25],[1 1],[1 1],1e-10,40) -> 1.9999999999627098.
    CHECK(r.T_[0] == doctest::Approx(1.9999999999627098).epsilon(1e-10));
}

TEST_CASE("M/M/1 DPS matches MATLAB on an asymmetric instance") {
    // lambda = [0.2 0.3], mu = [1 2], w = [1 3], maxCutoff = 40.
    auto r = line::qsys::qsys_mm1_dps(v3(0.2, 0.3), v3(1.0, 2.0), v3(1.0, 3.0), 1e-10, 40u);
    CHECK(r.rho == doctest::Approx(0.35).epsilon(1e-15));
    CHECK(r.T_[0] == doctest::Approx(1.6166883963494127).epsilon(1e-9));
    CHECK(r.T_[1] == doctest::Approx(0.66492829204693593).epsilon(1e-9));
    // The heavier weight gets more capacity, so its response time is smaller.
    CHECK(r.T_[1] < r.T_[0]);
}

TEST_CASE("M/M/1 DPS rejects unstable and non-positive input") {
    CHECK_THROWS_AS(line::qsys::qsys_mm1_dps(v3(2.0, 2.0), v3(1.0, 1.0), v3(1.0, 1.0), 1e-10, 20u),
                    line::InputError);
    CHECK_THROWS_AS(line::qsys::qsys_mm1_dps(v3(0.1, 0.1), v3(1.0, 1.0), v3(1.0, 0.0), 1e-10, 20u),
                    line::InputError);
}
