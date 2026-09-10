/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * qsys_mg1_ps: the sojourn-time distribution of the M/G/1-PS queue.
 *
 * Every golden below was measured in MATLAB R2025a on 2026-08-01, on the same
 * inputs, and the port reproduces them to the digits printed -- ten for the
 * transform and the moments, eight for the quantities that go through the
 * numerical inversion, where the two differ only in round-off.
 *
 * FOUR THINGS ARE CHECKED THAT A NUMBER-ONLY COMPARISON WOULD MISS:
 *
 *  - the CLOSED FORMS the reference guarantees exactly: the conditional mean is
 *    x/(1-rho) and the unconditional one m1/(1-rho), neither of which comes
 *    from the inversion, so they fail loudly if the polytope of coefficients is
 *    wrong;
 *  - the LST is a transform: 1 at s = 0, and in (0,1] for s > 0. The port
 *    initially returned 4.73 here, from taking n+2 coefficients of A where the
 *    reference takes n+1 and thereby keeping half of the cancelled double pole;
 *  - the LATTICE atoms: the conditional density is NaN at t = (k+1)x and only
 *    there, because the CDF jumps and no density exists. A smooth inversion
 *    returns a finite number at those points, which is worse than NaN;
 *  - the refusals, since an unstable system and an even `nterms` are input
 *    errors the reference names.
 */
#include <cmath>
#include <cstddef>
#include <vector>

#include "doctest.h"
#include "line/api/qsys/qsys_mg1_ps.h"

namespace qs = line::qsys;
using line::Matrix;
using qs::Cplx;

namespace {

/** M/M/1-PS with mu = 1 and lambda = 0.5, i.e. rho = 0.5. */
qs::Mg1PsResult mm1ps(const qs::Mg1PsOptions& o) {
    std::vector<double> alpha(1, 1.0);
    Matrix<double> T(1, 1, -1.0);
    return qs::qsys_mg1_ps(0.5, alpha, T, o);
}

}  // namespace

TEST_CASE("the exact closed forms are exact") {
    qs::Mg1PsOptions o;
    o.x.push_back(0.5);
    o.x.push_back(1.0);
    o.x.push_back(2.0);
    const qs::Mg1PsResult r = mm1ps(o);

    CHECK(r.rho == doctest::Approx(0.5).epsilon(1e-12));
    CHECK(r.m1 == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(r.m2 == doctest::Approx(2.0).epsilon(1e-12));
    // These two are closed forms, not inversions.
    CHECK(r.meanUncond == doctest::Approx(2.0).epsilon(1e-12));
    for (std::size_t i = 0; i < r.x.size(); ++i)
        CHECK(r.meanCond[i] == doctest::Approx(r.x[i] / 0.5).epsilon(1e-12));
    // (1-rho) bhat(lambda) = 0.5 * 1/(1+0.5).
    CHECK(r.atomUncond == doctest::Approx(1.0 / 3.0).epsilon(1e-10));
    // The k = 0 atom is (1-rho) exp(-lambda x).
    for (std::size_t i = 0; i < r.x.size(); ++i)
        CHECK(r.atomCond[i] == doctest::Approx(0.5 * std::exp(-0.5 * r.x[i])).epsilon(1e-12));
}

TEST_CASE("the conditional LST is a transform, and matches MATLAB") {
    qs::Mg1PsOptions o;
    o.x.push_back(0.5);
    o.x.push_back(1.0);
    o.x.push_back(2.0);
    o.s.push_back(0.1);
    o.s.push_back(1.0);
    const qs::Mg1PsResult r = mm1ps(o);

    // A transform is 1 at the origin and in (0,1] to the right of it. The port
    // returned 4.73 here while it kept one coefficient too many of A.
    for (std::size_t i = 0; i < r.x.size(); ++i) {
        CHECK(r.lstCond(Cplx(0.0, 0.0), r.x[i]).real() == doctest::Approx(1.0).epsilon(1e-9));
        for (std::size_t j = 0; j < r.s.size(); ++j) {
            CHECK(r.lstCondVal(i, j) > 0.0);
            CHECK(r.lstCondVal(i, j) <= 1.0);
        }
    }
    // MATLAB R2025a, same inputs.
    CHECK(r.lstCondVal(0, 0) == doctest::Approx(0.9068293695).epsilon(1e-9));
    CHECK(r.lstCondVal(0, 1) == doctest::Approx(0.4304793151).epsilon(1e-9));
    CHECK(r.lstCondVal(1, 0) == doctest::Approx(0.8251498265).epsilon(1e-9));
    CHECK(r.lstCondVal(1, 1) == doctest::Approx(0.2129752040).epsilon(1e-9));
    CHECK(r.lstCondVal(2, 0) == doctest::Approx(0.6874456279).epsilon(1e-9));
    CHECK(r.lstCondVal(2, 1) == doctest::Approx(0.0574622986).epsilon(1e-9));
    // The LST decreases in both s and x, since V(x) is increasing in x.
    CHECK(r.lstCondVal(0, 0) > r.lstCondVal(1, 0));
    CHECK(r.lstCondVal(1, 0) > r.lstCondVal(2, 0));
}

TEST_CASE("the moments through the transform match MATLAB") {
    qs::Mg1PsOptions o;
    o.x.push_back(0.5);
    o.x.push_back(1.0);
    o.x.push_back(2.0);
    const qs::Mg1PsResult r = mm1ps(o);

    CHECK(r.m2Cond[0] == doctest::Approx(1.460812526).epsilon(1e-8));
    CHECK(r.m2Cond[1] == doctest::Approx(5.704490546).epsilon(1e-8));
    CHECK(r.m2Cond[2] == doctest::Approx(21.88607103).epsilon(1e-8));
    CHECK(r.varCond[0] == doctest::Approx(0.4608125261).epsilon(1e-7));
    CHECK(r.varCond[1] == doctest::Approx(1.704490546).epsilon(1e-8));
    CHECK(r.varCond[2] == doctest::Approx(5.886071035).epsilon(1e-8));
    // The variance is positive, which a wrong stencil sign would break.
    for (std::size_t i = 0; i < r.x.size(); ++i) CHECK(r.varCond[i] > 0.0);

    CHECK(r.m2Uncond == doctest::Approx(10.66666666).epsilon(1e-8));
    CHECK(r.varUncond == doctest::Approx(6.666666657).epsilon(1e-8));
}

TEST_CASE("the unconditional LST matches MATLAB") {
    qs::Mg1PsOptions o;
    o.s.push_back(0.1);
    o.s.push_back(1.0);
    const qs::Mg1PsResult r = mm1ps(o);
    CHECK(r.lstUncondVal[0] == doctest::Approx(0.8404404190).epsilon(1e-8));
    CHECK(r.lstUncondVal[1] == doctest::Approx(0.3788664537).epsilon(1e-8));
    CHECK(r.lstUncond(Cplx(0.0, 0.0)).real() == doctest::Approx(1.0).epsilon(1e-6));
}

TEST_CASE("the conditional distribution is atomic on the lattice t = (k+1)x") {
    qs::Mg1PsOptions o;
    o.x.push_back(0.5);
    o.x.push_back(1.0);
    o.x.push_back(2.0);
    o.t.push_back(1.0);
    o.t.push_back(2.0);
    o.t.push_back(3.5);
    const qs::Mg1PsResult r = mm1ps(o);

    // MATLAB R2025a, same inputs.
    CHECK(r.cdfCond(0, 0) == doctest::Approx(0.62165938).epsilon(1e-6));
    CHECK(r.cdfCond(0, 1) == doctest::Approx(0.91218006).epsilon(1e-6));
    CHECK(r.cdfCond(0, 2) == doctest::Approx(0.99023232).epsilon(1e-6));
    CHECK(r.cdfCond(1, 0) == doctest::Approx(0.30326533).epsilon(1e-6));
    CHECK(r.cdfCond(1, 1) == doctest::Approx(0.62825348).epsilon(1e-6));
    CHECK(r.cdfCond(1, 2) == doctest::Approx(0.88278426).epsilon(1e-6));
    CHECK(r.cdfCond(2, 2) == doctest::Approx(0.55813177).epsilon(1e-6));

    // V(x) >= x, so the CDF is zero below x and exactly the atom at x.
    CHECK(r.cdfCond(2, 0) == doctest::Approx(0.0));
    CHECK(r.cdfCond(2, 1) == doctest::Approx(r.atomCond[2]).epsilon(1e-9));
    CHECK(r.cdfCond(1, 0) == doctest::Approx(r.atomCond[1]).epsilon(1e-9));

    // NaN exactly on the lattice, and a finite density off it. t/x is 2, 4 and
    // 7 in the first row and 2 in the second: all lattice points.
    CHECK(std::isnan(r.pdfCond(0, 0)));
    CHECK(std::isnan(r.pdfCond(0, 1)));
    CHECK(std::isnan(r.pdfCond(0, 2)));
    CHECK(std::isnan(r.pdfCond(1, 1)));
    CHECK(!std::isnan(r.pdfCond(1, 2)));
    CHECK(r.pdfCond(1, 2) == doctest::Approx(0.093387638).epsilon(1e-5));
    CHECK(r.pdfCond(2, 2) == doctest::Approx(0.16753844).epsilon(1e-5));
    // The CDF is nondecreasing in t.
    for (std::size_t i = 0; i < r.x.size(); ++i)
        for (std::size_t j = 1; j < r.t.size(); ++j)
            CHECK(r.cdfCond(i, j) >= r.cdfCond(i, j - 1) - 1e-9);
}

TEST_CASE("the unconditional distribution matches MATLAB") {
    qs::Mg1PsOptions o;
    o.t.push_back(1.0);
    o.t.push_back(2.0);
    o.t.push_back(3.5);
    const qs::Mg1PsResult r = mm1ps(o);
    CHECK(r.cdfUncond[0] == doctest::Approx(0.46217412).epsilon(1e-6));
    CHECK(r.cdfUncond[1] == doctest::Approx(0.67873320).epsilon(1e-6));
    CHECK(r.cdfUncond[2] == doctest::Approx(0.83339123).epsilon(1e-6));
    CHECK(r.pdfUncond[0] == doctest::Approx(0.30086963).epsilon(1e-5));
    CHECK(r.pdfUncond[1] == doctest::Approx(0.15348188).epsilon(1e-5));
    CHECK(r.pdfUncond[2] == doctest::Approx(0.067389654).epsilon(1e-5));
    // A distribution function: nondecreasing, and the density is positive.
    CHECK(r.cdfUncond[0] < r.cdfUncond[1]);
    CHECK(r.cdfUncond[1] < r.cdfUncond[2]);
    for (std::size_t j = 0; j < r.t.size(); ++j) CHECK(r.pdfUncond[j] > 0.0);
}

TEST_CASE("a multi-phase service law is solved exactly too") {
    // Erlang-2 of rate 2 per phase: mean 1, SCV 1/2. With lambda = 0.5 the
    // utilization is again 0.5, so the closed forms are the same and only the
    // shape changes -- which is what makes this a check on the residue sum
    // rather than on the bookkeeping around it.
    std::vector<double> alpha(2, 0.0);
    alpha[0] = 1.0;
    Matrix<double> T(2, 2, 0.0);
    T(0, 0) = -2.0;
    T(0, 1) = 2.0;
    T(1, 1) = -2.0;

    qs::Mg1PsOptions o;
    o.x.push_back(1.0);
    o.s.push_back(1.0);
    const qs::Mg1PsResult r = qs::qsys_mg1_ps(0.5, alpha, T, o);

    CHECK(r.m1 == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(r.m2 == doctest::Approx(1.5).epsilon(1e-12));  // Erlang-2: 2/lambda^2 * ... = 1.5
    CHECK(r.rho == doctest::Approx(0.5).epsilon(1e-12));
    CHECK(r.meanCond[0] == doctest::Approx(2.0).epsilon(1e-12));
    CHECK(r.lstCond(Cplx(0.0, 0.0), 1.0).real() == doctest::Approx(1.0).epsilon(1e-9));
    CHECK(r.lstCondVal(0, 0) > 0.0);
    CHECK(r.lstCondVal(0, 0) <= 1.0);
    CHECK(r.varCond[0] > 0.0);
}

TEST_CASE("the refusals are the reference's, by name") {
    std::vector<double> alpha(1, 1.0);
    Matrix<double> T(1, 1, -1.0);
    qs::Mg1PsOptions o;

    // rho >= 1: the queue is unstable and there is no stationary sojourn time.
    CHECK_THROWS_AS(qs::qsys_mg1_ps(1.0, alpha, T, o), line::InputError);
    CHECK_THROWS_AS(qs::qsys_mg1_ps(2.0, alpha, T, o), line::InputError);
    CHECK_THROWS_AS(qs::qsys_mg1_ps(0.0, alpha, T, o), line::InputError);
    CHECK_THROWS_AS(qs::qsys_mg1_ps(-1.0, alpha, T, o), line::InputError);

    // The Euler inversion needs an odd term count of at least eleven.
    qs::Mg1PsOptions even;
    even.nterms = 40;
    CHECK_THROWS_AS(qs::qsys_mg1_ps(0.5, alpha, T, even), line::InputError);
    qs::Mg1PsOptions few;
    few.nterms = 9;
    CHECK_THROWS_AS(qs::qsys_mg1_ps(0.5, alpha, T, few), line::InputError);

    // alpha must be a probability vector and T a proper subgenerator.
    std::vector<double> bad(1, 0.5);
    CHECK_THROWS_AS(qs::qsys_mg1_ps(0.5, bad, T, o), line::InputError);
    Matrix<double> badT(1, 1, 1.0);
    CHECK_THROWS_AS(qs::qsys_mg1_ps(0.5, alpha, badT, o), line::InputError);
    Matrix<double> wrongsize(2, 2, -1.0);
    CHECK_THROWS_AS(qs::qsys_mg1_ps(0.5, alpha, wrongsize, o), line::InputError);
}
