/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Conformance of the vendored RODAS (third_party/rodas.hpp) against the
 * unmodified Hairer-Wanner Fortran.
 *
 * THE ORACLE IS THE ORIGINAL, NOT A CLOSED FORM. rodas.hpp is an f2c
 * translation, so the property that matters is not that it solves these
 * problems well but that it solves them IDENTICALLY to the Fortran it came
 * from: a translation defect shows up as a last-digit drift long before it
 * shows up as a wrong answer. The expected values below were produced by
 *
 *     gfortran -O2 -std=legacy rodas.f dc_decsol.f decsol.f
 *
 * on the sources named in rodas.hpp's provenance block, printed at 17
 * significant digits so the decimal round-trips through a double exactly.
 *
 * THE STEP COUNTERS ARE CHECKED TOO, and are the more sensitive test: nfcn,
 * njac and nstep diverge on any change to the step-size controller or the
 * error estimate, while the final value can absorb a small perturbation and
 * still look right.
 *
 * The cases reach three different IJOB paths in DECOMR/SLVROD -- identity,
 * banded and full mass matrix -- crossed with analytic and numerical
 * Jacobians, because those are separate code paths and a merge error can hit
 * one and not the others. The singular mass matrix (index-1 DAE) is the case
 * LINE needs; the rest guard it against regressions elsewhere in the file.
 */
#include <cstring>

#include "doctest.h"
#include "rodas.hpp"

using namespace line::rodas_impl;

namespace {

/* Robertson's problem as an index-1 DAE:  M y' = f(y),  M = diag(1,1,0).
 * The third equation is the algebraic constraint y1+y2+y3 = 1. */
int frob_(integer *, doublereal *, doublereal *y, doublereal *f,
          doublereal *, integer *) {
    f[0] = -0.04 * y[0] + 1.0e4 * y[1] * y[2];
    f[1] = 0.04 * y[0] - 1.0e4 * y[1] * y[2] - 3.0e7 * y[1] * y[1];
    f[2] = y[0] + y[1] + y[2] - 1.0;
    return 0;
}
int jrob_(integer *, doublereal *, doublereal *y, doublereal *dfy,
          integer *ldfy, doublereal *, integer *) {
    const int ld = *ldfy;
    dfy[0 + 0 * ld] = -0.04;
    dfy[0 + 1 * ld] = 1.0e4 * y[2];
    dfy[0 + 2 * ld] = 1.0e4 * y[1];
    dfy[1 + 0 * ld] = 0.04;
    dfy[1 + 1 * ld] = -1.0e4 * y[2] - 6.0e7 * y[1];
    dfy[1 + 2 * ld] = -1.0e4 * y[1];
    dfy[2 + 0 * ld] = 1.0;
    dfy[2 + 1 * ld] = 1.0;
    dfy[2 + 2 * ld] = 1.0;
    return 0;
}
/* banded storage, MLMAS=MUMAS=0, i.e. the diagonal -- reaches IJOB=3 */
int mdiag3_(integer *, doublereal *am, integer *lmas, doublereal *, integer *) {
    const int ld = *lmas;
    am[0 + 0 * ld] = 1.0;
    am[0 + 1 * ld] = 1.0;
    am[0 + 2 * ld] = 0.0;
    return 0;
}
/* the same mass matrix stored full -- reaches IJOB=5, and must agree */
int mfull3_(integer *n, doublereal *am, integer *lmas, doublereal *, integer *) {
    const int ld = *lmas;
    for (int j = 0; j < *n; ++j)
        for (int i = 0; i < *n; ++i) am[i + j * ld] = 0.0;
    am[0 + 0 * ld] = 1.0;
    am[1 + 1 * ld] = 1.0;
    return 0;
}
/* Van der Pol at eps=1e-6: stiff, M = identity, reaches IJOB=1 */
int fvdp_(integer *, doublereal *, doublereal *y, doublereal *f,
          doublereal *, integer *) {
    const double eps = 1.0e-6;
    f[0] = y[1];
    f[1] = ((1.0 - y[0] * y[0]) * y[1] - y[0]) / eps;
    return 0;
}
int jvdp_(integer *, doublereal *, doublereal *y, doublereal *dfy,
          integer *ldfy, doublereal *, integer *) {
    const double eps = 1.0e-6;
    const int ld = *ldfy;
    dfy[0 + 0 * ld] = 0.0;
    dfy[0 + 1 * ld] = 1.0;
    dfy[1 + 0 * ld] = (-2.0 * y[0] * y[1] - 1.0) / eps;
    dfy[1 + 1 * ld] = (1.0 - y[0] * y[0]) / eps;
    return 0;
}
int sout0_(integer *, doublereal *, doublereal *, doublereal *, doublereal *,
           integer *, integer *, doublereal *, integer *, integer *) { return 0; }
int dfxd_(integer *, doublereal *, doublereal *, doublereal *, doublereal *,
          integer *) { return 0; }

struct Result {
    integer idid, nfcn, njac, nstep, nacc;
    doublereal y[3];
};

/* LWORK per the documented formula N*(LJAC+LMAS+LE1+14)+20; 400 is comfortably
 * above it for these sizes. Both work arrays must be zeroed for the defaults to
 * apply -- RODAS reads a zero as "use the default", not as a value. */
Result run_rob(int ijac_in, bool mass_full, int mlmas, int mumas,
               double rt, double at) {
    const int nd = 3, lwork = 400, liwork = 200;
    doublereal y[nd] = {1.0, 0.0, 0.0}, work[lwork], rpar[1], rtol[1], atol[1];
    integer iwork[liwork], ipar[1];
    integer n = nd, itol = 0, mljac = nd, mujac = nd, idid = 0;
    doublereal x = 0.0, xend = 0.4, h = 1.0e-6;
    rtol[0] = rt;
    atol[0] = at;
    std::memset(work, 0, sizeof work);
    std::memset(iwork, 0, sizeof iwork);
    integer ifcn = 0, ijac = ijac_in, idfx = 0, imas = 1, iout = 0;
    integer mlm = mlmas, mum = mumas, lw = lwork, liw = liwork;
    rodas_(&n, (U_fp)frob_, &ifcn, &x, y, &xend, &h, rtol, atol, &itol,
           (U_fp)jrob_, &ijac, &mljac, &mujac, (U_fp)dfxd_, &idfx,
           (U_fp)(mass_full ? mfull3_ : mdiag3_), &imas, &mlm, &mum,
           (U_fp)sout0_, &iout, work, &lw, iwork, &liw, rpar, ipar, &idid);
    Result r{idid, iwork[13], iwork[14], iwork[15], iwork[16], {y[0], y[1], y[2]}};
    return r;
}

Result run_vdp(int ijac_in, double rt, double at) {
    const int nd = 2, lwork = 400, liwork = 200;
    doublereal y[nd] = {2.0, -0.6}, work[lwork], rpar[1], rtol[1], atol[1];
    integer iwork[liwork], ipar[1];
    integer n = nd, itol = 0, mljac = nd, mujac = nd, idid = 0;
    doublereal x = 0.0, xend = 2.0, h = 1.0e-6;
    rtol[0] = rt;
    atol[0] = at;
    std::memset(work, 0, sizeof work);
    std::memset(iwork, 0, sizeof iwork);
    integer ifcn = 0, ijac = ijac_in, idfx = 0, imas = 0, iout = 0;
    integer mlm = 0, mum = 0, lw = lwork, liw = liwork;
    rodas_(&n, (U_fp)fvdp_, &ifcn, &x, y, &xend, &h, rtol, atol, &itol,
           (U_fp)jvdp_, &ijac, &mljac, &mujac, (U_fp)dfxd_, &idfx,
           (U_fp)mdiag3_, &imas, &mlm, &mum,
           (U_fp)sout0_, &iout, work, &lw, iwork, &liw, rpar, ipar, &idid);
    Result r{idid, iwork[13], iwork[14], iwork[15], iwork[16], {y[0], y[1], 0.0}};
    return r;
}

/* Bit equality, not a tolerance. The claim under test is that the translation
 * reproduces the Fortran exactly, and a tolerance would let a real drift pass. */
void check_exact(doublereal got, doublereal want) { CHECK(got == want); }

/* The two van der Pol cases at the bottom cannot be asserted bit exactly, and
 * the reason is libm rather than this translation. RODAS's step size
 * controller calls pow(err, 1/4) -- `pow_dd` in rodas.hpp, at the FAC1 and
 * FACGUS lines -- and glibc's pow is not correctly rounded, so it is not
 * portable across glibc releases. Logging every pow call of the first vdp
 * solve on the cluster's two images, call 274 takes the same input
 * 0x3fe4ff512a0207cf and returns 0x3fecccfcb80b2bd7 under glibc 2.35 (22.04)
 * against 0x3fecccfcb80b2bd6 under glibc 2.31 (20.04): one ulp. That perturbs
 * the accepted step, and from there the step sequence, the work counters and
 * the last digits of y all depend on the host. Measured across both images:
 *
 *   analytic  J: counters identical, y agrees to 1.4e-12 relative
 *   numerical J: nfcn 32755 vs 32663, nstep 5467 vs 5449, y to 1.1e-9
 *
 * Both solves are run at rtol 1e-8, so agreement is asserted at 1e-8. Asserting
 * tighter than the accuracy actually asked of the integrator would be pinning
 * the host's libm rather than this port -- which is what the bit exact form of
 * these two cases was doing, and it is why they failed on every 20.04 worker.
 * Do NOT tighten these back without removing the pow call first: replacing
 * pow(x,1/4) by sqrt(sqrt(x)) is bit portable (sqrt is correctly rounded by
 * IEEE 754) and was measured to give identical results on both images, but it
 * moves the answer off the vendored Fortran's, so it is a separate decision.
 *
 * The four Robertson cases above stay bit exact: they are short enough that
 * the controller never diverges, and that was verified on both images. */
void check_close(doublereal got, doublereal want) {
    CHECK(got == doctest::Approx(want).epsilon(1e-8));
}

/* Work counters move with the step sequence, so they are bounded rather than
 * fixed. The measured spread is under 0.4%; 5% still catches a real regression
 * in the work a solve costs. */
void check_count(integer got, integer want) {
    CHECK(static_cast<double>(got) >= 0.95 * static_cast<double>(want));
    CHECK(static_cast<double>(got) <= 1.05 * static_cast<double>(want));
}

}  // namespace

TEST_CASE("rodas: index-1 DAE, banded (diagonal) mass matrix, numerical Jacobian") {
    Result r = run_rob(0, false, 0, 0, 1.0e-8, 1.0e-10);
    CHECK(r.idid == 1);
    CHECK(r.nfcn == 281);
    CHECK(r.njac == 46);
    CHECK(r.nstep == 47);
    CHECK(r.nacc == 46);
    check_exact(r.y[0], 0.98517211385739623);
    check_exact(r.y[1], 0.33863953610620744e-4);
    check_exact(r.y[2], 0.14794022188993190e-1);
}

TEST_CASE("rodas: index-1 DAE, banded mass matrix, analytic Jacobian") {
    Result r = run_rob(1, false, 0, 0, 1.0e-8, 1.0e-10);
    CHECK(r.idid == 1);
    CHECK(r.nfcn == 281);
    CHECK(r.njac == 46);
    CHECK(r.nstep == 47);
    check_exact(r.y[0], 0.98517211385689285);
    check_exact(r.y[1], 0.33863953595594975e-4);
    check_exact(r.y[2], 0.14794022189511580e-1);
}

TEST_CASE("rodas: index-1 DAE at a tight tolerance") {
    Result r = run_rob(1, false, 0, 0, 1.0e-11, 1.0e-13);
    CHECK(r.idid == 1);
    CHECK(r.nfcn == 1510);
    CHECK(r.njac == 250);
    CHECK(r.nstep == 252);
    check_exact(r.y[0], 0.98517211386097781);
    check_exact(r.y[1], 0.33863953787836347e-4);
    check_exact(r.y[2], 0.14794022185234504e-1);
}

TEST_CASE("rodas: full mass matrix reaches the same answer as banded") {
    /* Same problem, different storage, therefore a different IJOB path through
     * DECOMR and SLVROD. Agreement here is what says the two paths were merged
     * correctly, without needing an external reference for either. */
    Result full = run_rob(1, true, 3, 3, 1.0e-8, 1.0e-10);
    Result band = run_rob(1, false, 0, 0, 1.0e-8, 1.0e-10);
    CHECK(full.idid == 1);
    CHECK(full.nfcn == band.nfcn);
    CHECK(full.nstep == band.nstep);
    for (int i = 0; i < 3; ++i) check_exact(full.y[i], band.y[i]);
    check_exact(full.y[0], 0.98517211385689285);
}

TEST_CASE("rodas: stiff ODE, identity mass matrix, analytic Jacobian") {
    Result r = run_vdp(1, 1.0e-8, 1.0e-10);
    CHECK(r.idid == 1);
    check_count(r.nfcn, 32446);
    check_count(r.njac, 5401);
    check_count(r.nstep, 5409);
    check_close(r.y[0], 0.17061674639889246e1);
    check_close(r.y[1], -0.89280998821562574);
}

TEST_CASE("rodas: stiff ODE, identity mass matrix, numerical Jacobian") {
    Result r = run_vdp(0, 1.0e-8, 1.0e-10);
    CHECK(r.idid == 1);
    check_count(r.nfcn, 32755);
    check_count(r.njac, 5420);
    check_count(r.nstep, 5467);
    check_close(r.y[0], 0.17061674637556288e1);
    check_close(r.y[1], -0.89280998834260528);
}
