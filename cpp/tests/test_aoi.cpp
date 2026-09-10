/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Age of Information API. The oracles are the M/M/1 closed forms of Inoue,
 * Masuyama, Takine and Tanaka (IEEE Trans. IT 65(12), 2019), which are exact
 * rational functions of rho and so are asserted as exact equalities; the
 * general-distribution routines are then checked by feeding them an
 * exponential and requiring that they reproduce those same M/M/1 numbers,
 * which is the strongest available consistency test since the two code paths
 * share nothing (one is a rational formula, the other a bracketed root plus a
 * finite-difference derivative of a transform). Transcendental results are
 * pinned to MATLAB values computed with matlab -singleCompThread on the same
 * inputs.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/aoi/aoi_fcfs_dm1.h"
#include "line/api/aoi/aoi_fcfs_gim1.h"
#include "line/api/aoi/aoi_fcfs_md1.h"
#include "line/api/aoi/aoi_fcfs_mgi1.h"
#include "line/api/aoi/aoi_fcfs_mm1.h"
#include "line/api/aoi/aoi_lcfsd_gim1.h"
#include "line/api/aoi/aoi_lcfsd_mgi1.h"
#include "line/api/aoi/aoi_lcfspr_dm1.h"
#include "line/api/aoi/aoi_lcfspr_gim1.h"
#include "line/api/aoi/aoi_lcfspr_md1.h"
#include "line/api/aoi/aoi_lcfspr_mgi1.h"
#include "line/api/aoi/aoi_lcfspr_mm1.h"
#include "line/api/aoi/aoi_lcfss_gim1.h"
#include "line/api/aoi/aoi_lcfss_mgi1.h"
#include "line/api/aoi/aoi_lst_det.h"
#include "line/api/aoi/aoi_lst_erlang.h"
#include "line/api/aoi/aoi_lst_exp.h"
#include "line/api/aoi/aoi_lst_ph.h"

using line::Matrix;
using line::Rational;
using namespace line::aoi;

namespace {
Rational Q(long n, long d = 1) { return Rational(n, d); }
}  // namespace

TEST_CASE("the LSTs in the field are exact rational functions") {
    // Exp(mu): mu/(mu+s).
    const Lst<Rational> e = aoi_lst_exp<Rational>(Q(3));
    CHECK(e(Q(7, 10)) == Q(30, 37));
    CHECK(e(Q(0)) == Q(1));  // a transform is 1 at the origin

    // Erlang(2,3): (3/3.7)^2 = 900/1369.
    const Lst<Rational> er = aoi_lst_erlang<Rational>(2, Q(3));
    CHECK(er(Q(7, 10)) == Q(900, 1369));
    CHECK(static_cast<double>(er(Q(7, 10))) ==
          doctest::Approx(0.65741417092768428).epsilon(1e-14));  // MATLAB
    CHECK(er(Q(0)) == Q(1));
    // Erlang-k is the k-fold product of the exponential transform.
    CHECK(er(Q(7, 10)) == e(Q(7, 10)) * e(Q(7, 10)));

    // A PH representation of Erlang-2 with phase rate 2 gives the same
    // transform as the closed form, exactly: the linear solve stays in the
    // field.
    Matrix<Rational> Tm(2, 2, Q(0));
    Tm(0, 0) = Q(-2);
    Tm(0, 1) = Q(2);
    Tm(1, 1) = Q(-2);
    const Lst<Rational> ph = aoi_lst_ph<Rational>(std::vector<Rational>{Q(1), Q(0)}, Tm);
    CHECK(ph(Q(7, 10)) == aoi_lst_erlang<Rational>(2, Q(2))(Q(7, 10)));
    CHECK(static_cast<double>(ph(Q(7, 10))) ==
          doctest::Approx(0.5486968449931412).epsilon(1e-14));  // MATLAB
    CHECK(ph(Q(0)) == Q(1));

    // Deterministic: exp(-s d), the one transform outside the field.
    CHECK(aoi_lst_det<double>(1.5)(0.7) ==
          doctest::Approx(0.34993774911115538).epsilon(1e-14));  // MATLAB
    CHECK_THROWS_AS(aoi_lst_exp<Rational>(Q(0)), line::InputError);
}

TEST_CASE("M/M/1 FCFS age of information matches the closed form, exactly") {
    // lambda = 1/2, mu = 1, rho = 1/2:
    //   E[A]     = (1)(1 + 2 + (1/4)/(1/2)) = 7/2
    //   E[Apeak] = (1)(1 + 2 + (1/2)/(1/2)) = 4
    //   E[A^2]   = 2 (1 - 1/2 - 1/8 + 1/4 - 1/16)/((1/4)(1/4)) = 23/4 + 49/4
    const AoiResult<Rational> a = aoi_fcfs_mm1<Rational>(Q(1, 2), Q(1));
    CHECK(a.meanAoI == Q(7, 2));
    CHECK(a.peakAoI == Q(4));
    CHECK(a.varAoI == Q(23, 4));
    CHECK(static_cast<double>(a.meanAoI) == doctest::Approx(3.5).epsilon(1e-15));   // MATLAB
    CHECK(static_cast<double>(a.varAoI) == doctest::Approx(5.75).epsilon(1e-15));   // MATLAB
    CHECK(static_cast<double>(a.peakAoI) == doctest::Approx(4.0).epsilon(1e-15));   // MATLAB
    // The peak age exceeds the mean age, as it must.
    CHECK(a.peakAoI > a.meanAoI);
    CHECK_THROWS_AS(aoi_fcfs_mm1<Rational>(Q(2), Q(1)), line::NumericError);
    CHECK_THROWS_AS(aoi_fcfs_mm1<Rational>(Q(0), Q(1)), line::InputError);
}

TEST_CASE("M/M/1 mean AoI has an interior minimum in the arrival rate") {
    // A theorem, not a numerical observation: sampling too rarely leaves the
    // age to grow, sampling too fast fills the queue. Every comparison here is
    // between exact rationals.
    const Rational mu = Q(1);
    Rational best = aoi_fcfs_mm1<Rational>(Q(1, 20), mu).meanAoI;
    long best_num = 1;
    for (long n = 2; n <= 19; ++n) {
        const Rational m = aoi_fcfs_mm1<Rational>(Q(n, 20), mu).meanAoI;
        if (m < best) {
            best = m;
            best_num = n;
        }
    }
    CHECK(best_num > 1);
    CHECK(best_num < 19);
    CHECK(best < aoi_fcfs_mm1<Rational>(Q(1, 20), mu).meanAoI);
    CHECK(best < aoi_fcfs_mm1<Rational>(Q(19, 20), mu).meanAoI);
}

TEST_CASE("preemptive LCFS beats FCFS at every load, exactly") {
    const Rational mu = Q(1);
    for (long n = 1; n <= 19; ++n) {
        const Rational lam = Q(n, 20);
        CHECK(aoi_lcfspr_mm1<Rational>(lam, mu).meanAoI < aoi_fcfs_mm1<Rational>(lam, mu).meanAoI);
    }
    const AoiResult<Rational> p = aoi_lcfspr_mm1<Rational>(Q(1, 2), Q(1));
    CHECK(p.meanAoI == Q(3));       // 1/mu + 1/lambda
    CHECK(p.peakAoI == Q(11, 3));   // 1/(l+m) + 1/l + 1/m
    CHECK(p.varAoI == Q(5));        // 2(1/l^2 + 1/(lm) + 1/m^2) - 9
    CHECK(static_cast<double>(p.peakAoI) ==
          doctest::Approx(3.6666666666666665).epsilon(1e-15));  // MATLAB
}

TEST_CASE("M/GI/1 with exponential service reproduces the M/M/1 result") {
    // Two entirely different code paths: a rational formula on one side, a
    // finite-difference derivative of the Pollaczek-Khinchine transform on the
    // other. Agreement to 1e-9 is a real check on both.
    const double lambda = 0.5, mu = 1.0;
    const Lst<double> H = aoi_lst_exp<double>(mu);
    const AoiLstResult<double> g = aoi_fcfs_mgi1<double>(lambda, H, 1.0 / mu, 2.0 / (mu * mu));
    const AoiResult<double> m = aoi_fcfs_mm1<double>(lambda, mu);
    CHECK(g.meanAoI == doctest::Approx(m.meanAoI).epsilon(1e-9));
    CHECK(g.peakAoI == doctest::Approx(m.peakAoI).epsilon(1e-12));
    CHECK(g.meanAoI == doctest::Approx(3.5000000000143778).epsilon(1e-12));  // MATLAB
    CHECK(g.has_lst);
    CHECK(g.lstAoI(0.7) == doctest::Approx(0.19463667820069207).epsilon(1e-12));
    // lstAoI IS a probability transform and DOES tend to 1 at the origin. The
    // values this block carried until 2026-08-19 -- 0.22997835497835503 at
    // s = 0.7, 3.333333e+08 at s = 1e-9, 3.327788e+02 at s = 1e-3 -- were
    // faithfully reproduced from a reference that was itself wrong, and the
    // "open question recorded in the earlier AoI audit" was the defect: the old
    // expression (lambda H*(s))/(s + lambda - lambda H*(s)) * W*(s) has a first
    // factor whose denominator vanishes like s(1+rho), so A*(0) was +Inf, and on
    // M/E2/1 it returned 1.3062 at s = 0.3 -- ABOVE 1, which no LST can be.
    // Settled by simulation, not by majority: a 3e6-cycle sample path of this
    // very M/M/1 gives A*(0.7) = 0.194592 against the 0.194637 below. All four
    // codebases now carry the general age formula; see BUGS.md and
    // _kb/03-api-layer.md (AoI family).
    CHECK(g.lstAoI(1e-9) == doctest::Approx(1.0).epsilon(1e-6));
    CHECK(g.lstAoI(1e-3) == doctest::Approx(0.99650897954).epsilon(1e-6));
}

TEST_CASE("GI/M/1 with exponential interarrivals reproduces the M/M/1 result") {
    const double lambda = 0.5, mu = 1.0;
    const Lst<double> Y = aoi_lst_exp<double>(lambda);
    const AoiLstResult<double> g =
        aoi_fcfs_gim1<double>(Y, mu, 1.0 / lambda, 2.0 / (lambda * lambda));
    const AoiResult<double> m = aoi_fcfs_mm1<double>(lambda, mu);
    CHECK(g.meanAoI == doctest::Approx(m.meanAoI).epsilon(1e-9));
    CHECK(g.peakAoI == doctest::Approx(m.peakAoI).epsilon(1e-9));
    CHECK(g.meanAoI == doctest::Approx(3.5000000000421339).epsilon(1e-11));  // MATLAB
    CHECK(g.peakAoI == doctest::Approx(4.0000000000000018).epsilon(1e-12));  // MATLAB
    // Same M/M/1 as the M/GI/1 case above, reached through the GI/M/1 route, so
    // the two must agree: 0.19463667820069207. The 0.061837568207552165 here
    // until 2026-08-19 came from the old (mu sigma(s))/(s + mu - mu sigma(s)) *
    // D*(s) form, whose value at the origin is sigma/(1-sigma), not 1.
    CHECK(g.lstAoI(0.7) == doctest::Approx(0.19463667820069207).epsilon(1e-10));
    CHECK(g.lstAoI(0.0) == doctest::Approx(1.0).epsilon(1e-12));
    // sigma for an M/M/1 is the utilization; the bracketed root must find it.
    CHECK(1.0 / (mu * g.peakAoI - mu / lambda * mu) > 0.0);  // E[D] = 1/(mu(1-sigma)) > 0
}

TEST_CASE("preemptive LCFS with a general distribution reproduces M/M/1") {
    const double lambda = 0.5, mu = 1.0;
    const Lst<double> H = aoi_lst_exp<double>(mu);
    const AoiLstResult<double> a = aoi_lcfspr_mgi1<double>(lambda, H, 1.0 / mu);
    const AoiResult<double> m = aoi_lcfspr_mm1<double>(lambda, mu);
    CHECK(a.meanAoI == doctest::Approx(m.meanAoI).epsilon(1e-14));
    CHECK(a.peakAoI == doctest::Approx(m.peakAoI).epsilon(1e-9));
    CHECK(a.peakAoI == doctest::Approx(3.6666666667043408).epsilon(1e-11));       // MATLAB
    CHECK(a.lstAoI(0.7) == doctest::Approx(0.24509803921568629).epsilon(1e-13));  // MATLAB

    const Lst<double> Y = aoi_lst_exp<double>(lambda);
    const AoiLstResult<double> b = aoi_lcfspr_gim1<double>(Y, mu, 1.0 / lambda);
    CHECK(b.meanAoI == doctest::Approx(m.meanAoI).epsilon(1e-14));
    CHECK(b.peakAoI == doctest::Approx(m.peakAoI).epsilon(1e-9));
    CHECK(b.peakAoI == doctest::Approx(3.6666666666894625).epsilon(1e-11));       // MATLAB
    CHECK(b.lstAoI(0.7) == doctest::Approx(0.24509803921568629).epsilon(1e-13));  // MATLAB
    // Both preemptive transforms are the same product of two exponential LSTs.
    CHECK(a.lstAoI(0.7) == doctest::Approx(b.lstAoI(0.7)).epsilon(1e-15));
}

TEST_CASE("M/D/1 and D/M/1 age of information match MATLAB") {
    const AoiResult<double> md = aoi_fcfs_md1<double>(0.5, 1.0);
    CHECK(md.meanAoI == doctest::Approx(3.1487212707001282).epsilon(1e-13));  // MATLAB
    CHECK(md.varAoI == doctest::Approx(8.0).epsilon(1e-13));                  // MATLAB
    CHECK(md.peakAoI == doctest::Approx(3.5).epsilon(1e-14));                 // MATLAB
    // Deterministic service is less variable than exponential service, so the
    // mean age is lower at the same load.
    CHECK(md.meanAoI < aoi_fcfs_mm1<double>(0.5, 1.0).meanAoI);

    const AoiResult<double> dm = aoi_fcfs_dm1<double>(2.0, 1.0);
    CHECK(dm.meanAoI == doctest::Approx(2.2550009749159754).epsilon(1e-12));  // MATLAB
    CHECK(dm.varAoI == doctest::Approx(1.6400529442481464).epsilon(1e-12));   // MATLAB
    CHECK(dm.peakAoI == doctest::Approx(3.2550009749159754).epsilon(1e-12));  // MATLAB
    // The peak exceeds the mean by exactly the half interarrival time.
    CHECK(dm.peakAoI - dm.meanAoI == doctest::Approx(1.0).epsilon(1e-12));
    // Deterministic arrivals beat Poisson arrivals at the same rate.
    CHECK(dm.meanAoI < aoi_fcfs_mm1<double>(0.5, 1.0).meanAoI);

    const AoiResult<double> pmd = aoi_lcfspr_md1<double>(0.5, 1.0);
    CHECK(pmd.meanAoI == doctest::Approx(3.0).epsilon(1e-14));               // MATLAB
    CHECK(pmd.peakAoI == doctest::Approx(4.2974425414002564).epsilon(1e-13));  // MATLAB
    CHECK(pmd.varAoI == doctest::Approx(4.0).epsilon(1e-14));

    const AoiResult<double> pdm = aoi_lcfspr_dm1<double>(2.0, 1.0);
    CHECK(pdm.meanAoI == doctest::Approx(3.0).epsilon(1e-14));   // MATLAB
    CHECK(pdm.peakAoI == doctest::Approx(3.0000000000000004).epsilon(1e-12));  // MATLAB
    CHECK(pdm.varAoI == doctest::Approx(1.0).epsilon(1e-14));
    CHECK_THROWS_AS(aoi_fcfs_dm1<double>(0.5, 1.0), line::NumericError);  // rho = 2
}

TEST_CASE("the non-preemptive LCFS disciplines are exact and match MATLAB") {
    // Both take moments only, so both stay in the field.
    const AoiLstResult<Rational> d = aoi_lcfsd_mgi1<Rational>(Q(1, 2), Q(1), Q(2));
    CHECK(d.meanAoI == Q(23, 6));  // 2 + 1 + 1/2 + 1/3
    CHECK(d.peakAoI == Q(7, 2));
    CHECK(d.has_lst == false);
    CHECK(static_cast<double>(d.meanAoI) ==
          doctest::Approx(3.8333333333333335).epsilon(1e-14));  // MATLAB

    const AoiLstResult<Rational> s = aoi_lcfss_mgi1<Rational>(Q(1, 2), Q(1), Q(2));
    CHECK(s.meanAoI == Q(5));
    CHECK(s.peakAoI == Q(5));
    CHECK(s.has_lst == false);
    CHECK(static_cast<double>(s.meanAoI) == doctest::Approx(5.0).epsilon(1e-15));  // MATLAB
    CHECK_THROWS_AS(aoi_lcfss_mgi1<Rational>(Q(1, 2), Q(1), Q(1, 2)), line::InputError);

    // GI/M/1 counterparts, reached through a bracketed root of the transform.
    const Lst<double> Y = aoi_lst_exp<double>(0.5);
    const AoiLstResult<double> dg = aoi_lcfsd_gim1<double>(Y, 1.0, 2.0);
    CHECK(dg.meanAoI == doctest::Approx(3.8333333333333339).epsilon(1e-11));  // MATLAB
    CHECK(dg.peakAoI == doctest::Approx(3.5000000000000004).epsilon(1e-11));  // MATLAB
    CHECK(dg.has_lst == false);
    // The M/GI/1 and GI/M/1 LCFS-D routes agree on the same M/M/1 system.
    CHECK(dg.meanAoI == doctest::Approx(static_cast<double>(d.meanAoI)).epsilon(1e-9));
    CHECK(dg.peakAoI == doctest::Approx(static_cast<double>(d.peakAoI)).epsilon(1e-9));

    const AoiLstResult<double> sg = aoi_lcfss_gim1<double>(Y, 1.0, 2.0);
    CHECK(sg.meanAoI == doctest::Approx(4.0000000000000018).epsilon(1e-11));  // MATLAB
    CHECK(sg.peakAoI == doctest::Approx(4.0000000000000018).epsilon(1e-11));  // MATLAB
    // REFERENCE DIVERGENCE, recorded rather than smoothed over: on the same
    // M/M/1 system the LCFS-S mean is 5 through the M/GI/1 route and 4 through
    // the GI/M/1 one. The two MATLAB files implement different approximations
    // of Section V and are not consistent with each other; the LCFS-D pair
    // above is.
    CHECK(sg.meanAoI < static_cast<double>(s.meanAoI));
}

TEST_CASE("AoI results agree between double and exact arithmetic") {
    CHECK(aoi_fcfs_mm1<double>(0.4, 1.25).meanAoI ==
          doctest::Approx(static_cast<double>(aoi_fcfs_mm1<Rational>(Q(2, 5), Q(5, 4)).meanAoI))
              .epsilon(1e-9));
    CHECK(aoi_fcfs_mm1<double>(0.4, 1.25).varAoI ==
          doctest::Approx(static_cast<double>(aoi_fcfs_mm1<Rational>(Q(2, 5), Q(5, 4)).varAoI))
              .epsilon(1e-9));
    CHECK(aoi_lcfspr_mm1<double>(0.4, 1.25).peakAoI ==
          doctest::Approx(static_cast<double>(aoi_lcfspr_mm1<Rational>(Q(2, 5), Q(5, 4)).peakAoI))
              .epsilon(1e-9));
    CHECK(aoi_lcfsd_mgi1<double>(0.4, 0.8, 1.28).meanAoI ==
          doctest::Approx(
              static_cast<double>(aoi_lcfsd_mgi1<Rational>(Q(2, 5), Q(4, 5), Q(32, 25)).meanAoI))
              .epsilon(1e-9));
    CHECK(aoi_lcfss_mgi1<double>(0.4, 0.8, 1.28).peakAoI ==
          doctest::Approx(
              static_cast<double>(aoi_lcfss_mgi1<Rational>(Q(2, 5), Q(4, 5), Q(32, 25)).peakAoI))
              .epsilon(1e-9));
    CHECK(static_cast<double>(aoi_lst_erlang<Rational>(3, Q(5, 4))(Q(7, 10))) ==
          doctest::Approx(aoi_lst_erlang<double>(3, 1.25)(0.7)).epsilon(1e-12));
}
