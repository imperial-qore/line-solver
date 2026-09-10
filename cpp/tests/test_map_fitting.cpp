/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * MAP / PH closed-form moment matching.
 *
 * Three kinds of oracle are used.
 *
 * 1. Round-tripping. A fit is correct when the process it produces reproduces
 *    the characteristics it was fitted to: map_moment 1..3 must return the
 *    target moments, and the autocorrelation decay rate, recovered as
 *    (acf[k+1] - 1)/(acf[k] - 1) from the unnormalized kpctoolbox acf, must
 *    return the target gamma. This is checked to the precision the closed form
 *    allows (1e-9 relative), not to an eyeballed tolerance.
 *
 * 2. Structural feasibility. Every fitted (D0, D1) is put through
 *    map_isfeasible, which checks the off-diagonal signs, the row sums of
 *    D0 + D1 and the stochasticity of the embedded chain. A fit handed an
 *    infeasible moment set (an SCV below the Erlang bound 1/n) must relax the
 *    moments and still return a feasible generator, never a generator with
 *    negative rates.
 *
 * 3. MATLAB reference values, obtained by running the reference implementation
 *    directly (matlab -singleCompThread -batch over matlab/lib/kpctoolbox and
 *    matlab/lib/m3a). Ten fits are pinned entry by entry.
 *
 * The predicate-only and assembly-only routines are additionally instantiated
 * at Rational and asserted with exact equality, since they contain no root.
 */
#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/mam/amap2_assemble.h"
#include "line/api/mam/amap2_fit_gamma.h"
#include "line/api/mam/amap2_fitall_gamma.h"
#include "line/api/mam/aph2_adjust.h"
#include "line/api/mam/aph2_assemble.h"
#include "line/api/mam/aph2_fit.h"
#include "line/api/mam/aph2_fitall.h"
#include "line/api/mam/aph_fit.h"
#include "line/api/mam/m3pp2m_fitc.h"
#include "line/api/mam/map2_fit.h"
#include "line/api/mam/map2mmpp.h"
#include "line/api/mam/map_mark.h"
#include "line/api/mam/mmpp2_fit.h"
#include "line/api/mam/mmpp2_fit1.h"
#include "line/api/mam/mmpp2_fit2.h"
#include "line/api/mam/mmpp2_fit3.h"
#include "line/api/mam/mmpp2_fit4.h"
#include "line/api/mam/mmpp2_fitc.h"

using line::Matrix;
using line::Rational;
using line::mam::Map;
using namespace line::mam;

namespace {

const double kRel = 1e-9;

/** Autocorrelation decay rate recovered from the fitted process. */
double decay(const Map<double>& m, unsigned k) {
    // map_acf now carries map_acf.m's own (x - 1)/scv normalization, so the
    // ratio is taken directly; scv cancels in it either way, but the -1 does
    // not, which is why this used to be written out by hand here.
    const std::vector<double> a = map_acf(m, std::vector<unsigned>{k, k + 1});
    return a[1] / a[0];
}

/** Lag-1 autocorrelation coefficient of the fitted process. */
double acf1(const Map<double>& m) {
    return map_acf(m, std::vector<unsigned>{1})[0];
}

void check_moments(const Map<double>& m, double e1, double e2, double e3) {
    CHECK(map_moment(m, 1) == doctest::Approx(e1).epsilon(kRel));
    CHECK(map_moment(m, 2) == doctest::Approx(e2).epsilon(kRel));
    CHECK(map_moment(m, 3) == doctest::Approx(e3).epsilon(kRel));
}

void check_feasible(const Map<double>& m) {
    CHECK(map_isfeasible(m, 1e-10));
}

void check_entries(const Map<double>& m, const std::vector<double>& d0,
                   const std::vector<double>& d1, double tol) {
    const std::size_t n = m.order();
    REQUIRE(d0.size() == n * n);
    REQUIRE(d1.size() == n * n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) {
            CHECK(m.D0(i, j) == doctest::Approx(d0[i * n + j]).epsilon(tol).scale(1.0));
            CHECK(m.D1(i, j) == doctest::Approx(d1[i * n + j]).epsilon(tol).scale(1.0));
        }
}

}  // namespace

// ---------------------------------------------------------------------------
// Exact (rational) assembly and predicates
// ---------------------------------------------------------------------------

TEST_CASE("aph2_assemble is exact in rational arithmetic") {
    const Rational l1 = line::num_traits<Rational>::from_rational(1, 2);
    const Rational l2 = line::num_traits<Rational>::from_rational(1, 3);
    const Rational p1 = line::num_traits<Rational>::from_rational(1, 4);
    const Map<Rational> m = aph2_assemble(l1, l2, p1);

    CHECK(m.D0(0, 0) == Rational(-2));
    CHECK(m.D0(0, 1) == Rational(1, 2));
    CHECK(m.D0(1, 0) == Rational(0));
    CHECK(m.D0(1, 1) == Rational(-3));
    CHECK(m.D1(0, 0) == Rational(3, 2));
    CHECK(m.D1(0, 1) == Rational(0));
    CHECK(m.D1(1, 0) == Rational(3));
    CHECK(m.D1(1, 1) == Rational(0));

    // Every row of D0 + D1 vanishes exactly, and the MAP is feasible with a
    // zero tolerance -- no rounding residual to absorb.
    for (std::size_t i = 0; i < 2; ++i) {
        Rational s(0);
        for (std::size_t j = 0; j < 2; ++j) s += m.D0(i, j) + m.D1(i, j);
        CHECK(s == Rational(0));
    }
    CHECK(map_isfeasible(m));
}

TEST_CASE("amap2_assemble is exact in both canonical forms") {
    const Rational l1 = line::num_traits<Rational>::from_rational(1, 2);
    const Rational l2 = line::num_traits<Rational>::from_rational(1, 5);
    const Rational p1 = line::num_traits<Rational>::from_rational(1, 4);
    const Rational p2 = line::num_traits<Rational>::from_rational(2, 5);

    const Map<Rational> f1 = amap2_assemble(l1, l2, p1, p2, 1);
    CHECK(f1.D0(0, 0) == Rational(-2));
    CHECK(f1.D0(0, 1) == Rational(1, 2));
    CHECK(f1.D1(0, 0) == Rational(3, 2));
    CHECK(f1.D1(0, 1) == Rational(0));
    CHECK(f1.D1(1, 0) == Rational(3));
    CHECK(f1.D1(1, 1) == Rational(2));
    CHECK(map_isfeasible(f1));

    const Map<Rational> f2 = amap2_assemble(l1, l2, p1, p2, 2);
    CHECK(f2.D1(0, 0) == Rational(0));
    CHECK(f2.D1(0, 1) == Rational(3, 2));
    CHECK(map_isfeasible(f2));

    CHECK_THROWS_AS(amap2_assemble(l1, l2, p1, p2, 3), line::InputError);
}

TEST_CASE("map_mark splits D1 exactly and preserves the inter-arrival process") {
    Map<Rational> base;
    base.D0 = Matrix<Rational>(2, 2, Rational(0));
    base.D1 = Matrix<Rational>(2, 2, Rational(0));
    base.D0(0, 0) = Rational(-3);
    base.D0(0, 1) = Rational(1);
    base.D0(1, 1) = Rational(-4);
    base.D1(0, 0) = Rational(2);
    base.D1(1, 0) = Rational(4);

    std::vector<Rational> prob;
    prob.push_back(Rational(1, 3));
    prob.push_back(Rational(2, 3));
    const Mmap<Rational> mm = map_mark(base, prob);

    REQUIRE(mm.classes() == 2u);
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j) {
            CHECK(mm.Dc[0](i, j) + mm.Dc[1](i, j) == base.D1(i, j));
            CHECK(mm.D1(i, j) == base.D1(i, j));
            CHECK(mm.D0(i, j) == base.D0(i, j));
        }
    CHECK(mm.Dc[0](0, 0) == Rational(2, 3));
    CHECK(mm.Dc[1](0, 0) == Rational(4, 3));

    // Unnormalized weights are renormalized, degenerate ones rejected.
    std::vector<Rational> weights;
    weights.push_back(Rational(1));
    weights.push_back(Rational(2));
    const Mmap<Rational> mw = map_mark(base, weights);
    CHECK(mw.Dc[0](0, 0) == Rational(2, 3));

    std::vector<Rational> zero_w(2, Rational(0));
    CHECK_THROWS_AS(map_mark(base, zero_w), line::InputError);
}

TEST_CASE("map2mmpp separates the generator from the arrival rates") {
    Map<Rational> mmpp;
    mmpp.D0 = Matrix<Rational>(2, 2, Rational(0));
    mmpp.D1 = Matrix<Rational>(2, 2, Rational(0));
    mmpp.D0(0, 0) = Rational(-3);
    mmpp.D0(0, 1) = Rational(1);
    mmpp.D0(1, 0) = Rational(2);
    mmpp.D0(1, 1) = Rational(-5);
    mmpp.D1(0, 0) = Rational(2);
    mmpp.D1(1, 1) = Rational(3);

    const Map2mmppResult<Rational> r = map2mmpp(mmpp);
    CHECK(r.is_mmpp);
    CHECK(r.offdiag_norm == Rational(0));
    CHECK(r.Q(0, 0) == Rational(-1));
    CHECK(r.Q(0, 1) == Rational(1));
    CHECK(r.Q(1, 0) == Rational(2));
    CHECK(r.Q(1, 1) == Rational(-2));
    CHECK(r.LAMBDA(0, 0) == Rational(2));

    Map<Rational> notmmpp = mmpp;
    notmmpp.D1(0, 1) = Rational(1);
    notmmpp.D0(0, 1) = Rational(0);
    const Map2mmppResult<Rational> r2 = map2mmpp(notmmpp);
    CHECK_FALSE(r2.is_mmpp);
    CHECK(r2.offdiag_norm == Rational(1));
}

// ---------------------------------------------------------------------------
// APH(2)
// ---------------------------------------------------------------------------

TEST_CASE("aph2_fitall reproduces the moments it was given") {
    const std::vector<Map<double>> a = aph2_fitall<double>(1.0, 2.5, 12.0);
    REQUIRE(a.size() == 1u);
    check_moments(a[0], 1.0, 2.5, 12.0);
    check_feasible(a[0]);

    // MATLAB: aph2_fitall(1,2.5,12) -> a single APH(2).
    check_entries(a[0], {-1.2612038749637415, 0.11834673210659863, 0.0, -0.45308183932197288},
                  {1.1428571428571428, 0.0, 0.45308183932197288, 0.0}, 1e-12);
}

TEST_CASE("aph2_fitall enumerates both branches when they are feasible") {
    // A hyperexponential moment set: both roots of the discriminant give
    // positive phase means and a probability in [0, 1].
    const std::vector<Map<double>> a = aph2_fitall<double>(1.0, 8.0, 200.0);
    REQUIRE(a.size() == 1u);
    check_moments(a[0], 1.0, 8.0, 200.0);
    check_feasible(a[0]);
    // MATLAB: aph2_fitall(1,8,200) -> a single APH(2); the second root of the
    // discriminant gives a branching probability outside [0, 1].
    check_entries(a[0], {-1.5829709254103086, 0.063740156179539112, 0.0, -0.10933676689738396},
                  {1.5192307692307694, 0.0, 0.10933676689738396, 0.0}, 1e-12);
}

TEST_CASE("aph2_fit refuses an infeasible moment set instead of fabricating rates") {
    // SCV = 0.2 is below the APH(2) bound of 1/2: no two-phase acyclic PH can
    // match it. The reference falls back to aph_fit with nmax = 2, which
    // relaxes the moments to the Erlang-2 corner rather than returning a
    // generator with negative entries.
    const std::vector<Map<double>> a = aph2_fitall<double>(1.0, 1.2, 3.0);
    REQUIRE(a.size() == 1u);
    check_feasible(a[0]);
    CHECK(a[0].order() == 2u);
    // Relaxed to the n = 2 corner n2 = 3/2, n3 = 2 n2 - 1 = 2.
    CHECK(map_moment(a[0], 2) == doctest::Approx(1.5).epsilon(kRel));
    CHECK(map_moment(a[0], 3) == doctest::Approx(3.0).epsilon(kRel));
    // Off-diagonal signs and row sums are intact.
    CHECK(a[0].D0(0, 1) >= 0.0);
    CHECK(a[0].D0(0, 0) <= 0.0);
    CHECK(a[0].D0(1, 1) <= 0.0);
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j) CHECK(a[0].D1(i, j) >= 0.0);

    const Aph2FitResult<double> f = aph2_fit<double>(1.0, 1.2, 3.0);
    check_feasible(f.aph);
}

TEST_CASE("aph2_adjust lifts infeasible characteristics onto the boundary") {
    // SCV below 1/2 is lifted to exactly 1/2, i.e. M2 = 3/2 M1^2.
    const Aph2AdjustResult<double> a = aph2_adjust<double>(1.0, 1.2, 3.0);
    CHECK(a.M2a == doctest::Approx(1.5).epsilon(kRel));
    // With scva = 1/2 the third moment is clamped into [lb, ub].
    const double lb = 3.0 * (3.0 * 0.5 - 1.0 + std::sqrt(2.0) * std::pow(0.5, 1.5));
    const double ub = 6.0 * 0.5;
    CHECK(a.M3a >= lb - 1e-12);
    CHECK(a.M3a <= ub + 1e-12);
    // A feasible set is returned untouched.
    const Aph2AdjustResult<double> b = aph2_adjust<double>(1.0, 2.5, 12.0);
    CHECK(b.M2a == doctest::Approx(2.5).epsilon(kRel));
    CHECK(b.M3a == doctest::Approx(12.0).epsilon(kRel));
    // SCV > 1 pushes the third moment strictly above 3/2 M1^3 (1 + scv)^2.
    // scv = 2, lb = 3/2 M1^3 (1 + scv)^2 = 13.5, raised by the 1e-4 slack.
    const Aph2AdjustResult<double> c = aph2_adjust<double>(1.0, 3.0, 1.0);
    CHECK(c.M3a == doctest::Approx(13.5 * (1.0 + 1e-4)).epsilon(kRel));
}

// ---------------------------------------------------------------------------
// APH(n) -- Bobbio, Horvath and Telek
// ---------------------------------------------------------------------------

TEST_CASE("aph_fit matches three moments with a minimal-order APH") {
    // MATLAB: aph_fit(1,3,20,10) -> exact APH(2), second fitting case.
    const AphFitResult<double> a = aph_fit<double>(1.0, 3.0, 20.0, 10u);
    CHECK(a.isexact);
    CHECK(a.order == 2u);
    check_moments(a.aph, 1.0, 3.0, 20.0);
    check_feasible(a.aph);
    check_entries(a.aph, {-0.34173549813061532, 0.34173549813061532, 0.0, -1.3505721941770767},
                  {0.0, 0.0, 0.11980296340784612, 1.2307692307692306}, 1e-12);

    // MATLAB: aph_fit(1,1.5,3.5,10) -> exact APH(4), first fitting case.
    const AphFitResult<double> b = aph_fit<double>(1.0, 1.5, 3.5, 10u);
    CHECK(b.isexact);
    CHECK(b.order == 4u);
    check_moments(b.aph, 1.0, 1.5, 3.5);
    check_feasible(b.aph);
    CHECK(b.aph.D0(0, 0) == doctest::Approx(-0.78000353558888513).epsilon(1e-12));
    CHECK(b.aph.D0(0, 1) == doctest::Approx(0.78000353558888513).epsilon(1e-12));
    CHECK(b.aph.D0(1, 1) == doctest::Approx(-3.2993527191071674).epsilon(1e-12));
    CHECK(b.aph.D1(3, 0) == doctest::Approx(0.23349617929173699).epsilon(1e-12));
    CHECK(b.aph.D1(3, 1) == doctest::Approx(3.0658565398154303).epsilon(1e-12));
}

TEST_CASE("aph_fit reports inexactness rather than returning an invalid PH") {
    // SCV = 0.05 needs at least 20 phases; with nmax = 10 the reference
    // relaxes the moment set to the n = 10 corner and flags it.
    const AphFitResult<double> a = aph_fit<double>(1.0, 1.05, 1.2, 10u);
    CHECK_FALSE(a.isexact);
    CHECK(a.order == 10u);
    check_feasible(a.aph);
    // Relaxed to n2 = (n+1)/n, n3 = 2 n2 - 1.
    CHECK(map_moment(a.aph, 2) == doctest::Approx(1.1).epsilon(kRel));
    CHECK(map_moment(a.aph, 3) == doctest::Approx(1.32).epsilon(kRel));
    // Nothing negative leaked into the representation.
    for (std::size_t i = 0; i < a.aph.order(); ++i) {
        CHECK(a.aph.D0(i, i) <= 0.0);
        for (std::size_t j = 0; j < a.aph.order(); ++j) {
            if (i != j) CHECK(a.aph.D0(i, j) >= -1e-14);
            CHECK(a.aph.D1(i, j) >= -1e-14);
        }
    }
    CHECK_THROWS_AS(aph_fit<double>(1.0, 3.0, 20.0, 1u), line::InputError);
    CHECK_THROWS_AS(aph_fit<double>(0.0, 3.0, 20.0, 10u), line::InputError);
}

TEST_CASE("aph_fit returns the exponential on the degenerate moment set") {
    const AphFitResult<double> a = aph_fit<double>(2.0, 8.0, 48.0);
    CHECK(a.order == 1u);
    CHECK(a.isexact);
    CHECK(map_moment(a.aph, 1) == doctest::Approx(2.0).epsilon(kRel));
    CHECK(map_scv(a.aph) == doctest::Approx(1.0).epsilon(kRel));
}

// ---------------------------------------------------------------------------
// MAP(2) / AMAP(2)
// ---------------------------------------------------------------------------

TEST_CASE("map2_fit matches three moments and the lag-1 decay rate") {
    // MATLAB: map2_fit(1,4,30,0.3) -> ERR = 0.
    const Map2FitResult<double> r = map2_fit<double>(1.0, 4.0, 30.0, 0.3);
    REQUIRE(r.has_map);
    CHECK(r.err == 0);
    check_moments(r.map, 1.0, 4.0, 30.0);
    CHECK(decay(r.map, 1) == doctest::Approx(0.3).epsilon(kRel));
    check_feasible(r.map);
    check_entries(r.map, {-0.3819660112501051, 0.0, 0.0, -2.6180339887498949},
                  {0.18849076967509043, 0.19347524157501467, 0.50652475842498534,
                   2.1115092303249097},
                  1e-12);
}

TEST_CASE("map2_fit rejects out-of-range correlation with the reference ERR code") {
    // g2 >= 1 is outside the representable range of an AMAP(2).
    // MATLAB returns ERR = 51 here. Note that the JAR (Map2_fit.java, the
    // b >= 0 hyperexponential branch) drops MATLAB's `&& g2 < 1` guard and
    // returns a MAP instead; MATLAB is ground truth, so this port refuses.
    const Map2FitResult<double> r = map2_fit<double>(1.0, 4.0, 30.0, 1.5);
    CHECK_FALSE(r.has_map);
    CHECK(r.err == 51);
    // A hypoexponential marginal (h2 = -0.2, h3 = -0.05) with a correlation
    // past its bound -(h2 + sqrt(-h3))^2/h2 = 0.00279.
    const Map2FitResult<double> s = map2_fit<double>(1.0, 1.6, 3.54, 0.9);
    CHECK_FALSE(s.has_map);
    CHECK(s.err == 53);
    // h3 = 0.06 with h2 = -0.2 sits in neither region: h3 out of bounds.
    const Map2FitResult<double> u = map2_fit<double>(1.0, 1.6, 4.2, 0.9);
    CHECK_FALSE(u.has_map);
    CHECK(u.err == 30);
}

TEST_CASE("map2_fit selects a third moment when given a sentinel") {
    // e3 = -1 maximizes the range of feasible correlations.
    const Map2FitResult<double> r = map2_fit<double>(1.0, 4.0, 0.4);
    REQUIRE(r.has_map);
    CHECK(r.e3_used == doctest::Approx(1.501 * 16.0).epsilon(kRel));
    check_moments(r.map, 1.0, 4.0, r.e3_used);
    CHECK(decay(r.map, 1) == doctest::Approx(0.4).epsilon(kRel));
    check_feasible(r.map);
}

TEST_CASE("amap2_fitall_gamma enumerates the AMAP(2)s matching gamma") {
    // MATLAB: amap2_fitall_gamma(1,4,30,0.5) -> two solutions.
    const std::vector<Map<double>> a = amap2_fitall_gamma<double>(1.0, 4.0, 30.0, 0.5);
    REQUIRE(a.size() == 2u);
    for (std::size_t k = 0; k < a.size(); ++k) {
        check_moments(a[k], 1.0, 4.0, 30.0);
        CHECK(decay(a[k], 1) == doctest::Approx(0.5).epsilon(kRel));
        check_feasible(a[k]);
    }
    check_entries(a[0], {-2.6180339887498891, 0.33725758234547876, 0.0, -0.38196601125010526},
                  {2.2807764064044105, 0.0, 0.16274241765452038, 0.21922359359558488}, 1e-11);
    check_entries(a[1], {-0.38196601125010526, 0.16274241765452038, 0.0, -2.6180339887498891},
                  {0.21922359359558491, 0.0, 0.33725758234547903, 2.28077640640441}, 1e-11);
}

TEST_CASE("amap2_fitall_gamma uses the second canonical form for negative gamma") {
    const std::vector<Map<double>> a = amap2_fitall_gamma<double>(1.0, 4.0, 30.0, -0.3);
    REQUIRE(a.size() >= 1u);
    for (std::size_t k = 0; k < a.size(); ++k) {
        check_moments(a[k], 1.0, 4.0, 30.0);
        CHECK(decay(a[k], 1) == doctest::Approx(-0.3).epsilon(kRel));
        check_feasible(a[k]);
    }
    // Form 2 puts the escape mass of phase 1 in the (0,1) entry of D1.
    REQUIRE(a.size() == 1u);
    // MATLAB: amap2_fitall_gamma(1,4,30,-0.3).
    check_entries(a[0], {-0.38196601125010526, 0.055180724783188841, 0.0, -2.6180339887498891},
                  {0.0, 0.3267852864669164, 0.91803398874989228, 1.6999999999999968}, 1e-11);
}

TEST_CASE("amap2_fitall_gamma returns nothing for infeasible characteristics") {
    // SCV below 1/2 makes the moment discriminant negative; unlike
    // aph2_fitall there is no approximate fallback.
    const std::vector<Map<double>> a = amap2_fitall_gamma<double>(1.0, 1.2, 3.0, 0.4);
    CHECK(a.empty());
}

TEST_CASE("amap2_fit_gamma falls back to Poisson when no exact fit exists") {
    const Amap2FitGammaResult<double> r = amap2_fit_gamma<double>(1.0, 4.0, 30.0, 0.5);
    CHECK_FALSE(r.poisson_fallback);
    CHECK(r.amaps.size() == 2u);
    check_moments(r.amap, 1.0, 4.0, 30.0);
    check_feasible(r.amap);

    // Unit SCV: the canonical forms degenerate, a Poisson process is returned.
    const Amap2FitGammaResult<double> p = amap2_fit_gamma<double>(2.0, 8.0, 48.0, 0.3);
    CHECK(p.poisson_fallback);
    CHECK(p.amap.order() == 1u);
    CHECK(map_moment(p.amap, 1) == doctest::Approx(2.0).epsilon(kRel));

    // Infeasible: no exact solution, so the documented Poisson fallback fires.
    const Amap2FitGammaResult<double> q = amap2_fit_gamma<double>(1.0, 1.2, 3.0, 0.4);
    CHECK(q.poisson_fallback);
    check_feasible(q.amap);
}

// ---------------------------------------------------------------------------
// MMPP(2)
// ---------------------------------------------------------------------------

TEST_CASE("mmpp2_fit3 matches three moments and the decay rate") {
    // MATLAB: mmpp2_fit3(1,4,30,0.5).
    const Map<double> m = mmpp2_fit3<double>(1.0, 4.0, 30.0, 0.5);
    check_moments(m, 1.0, 4.0, 30.0);
    CHECK(decay(m, 1) == doctest::Approx(0.5).epsilon(kRel));
    check_feasible(m);
    check_entries(m, {-0.40858968733650131, 0.18936609374091656, 0.31063390625908333,
                      -2.5914103126634984},
                  {0.21922359359558474, 0.0, 0.0, 2.2807764064044149}, 1e-11);
    // The result really is an MMPP: D1 is diagonal.
    const Map2mmppResult<double> mm = map2mmpp(m);
    CHECK(mm.offdiag_norm == doctest::Approx(0.0).epsilon(1e-14).scale(1.0));
}

TEST_CASE("mmpp2_fit3 takes the uncorrelated branch at vanishing G2") {
    const Map<double> m = mmpp2_fit3<double>(1.0, 4.0, 30.0, 0.0);
    check_moments(m, 1.0, 4.0, 30.0);
    // MATLAB: mmpp2_fit3(1,4,30,0) -> D0 = [-5/2 1/2; 1/2 -1/2], D1 = diag(2,0).
    check_entries(m, {-2.5, 0.5, 0.5, -0.5}, {2.0, 0.0, 0.0, 0.0}, 1e-12);
}

TEST_CASE("mmpp2_fit matches the lag-1 autocorrelation") {
    // MATLAB: mmpp2_fit(1,4,30,0.2).
    const Map<double> m = mmpp2_fit<double>(1.0, 4.0, 30.0, 0.2);
    check_moments(m, 1.0, 4.0, 30.0);
    CHECK(acf1(m) == doctest::Approx(0.2).epsilon(kRel));
    check_feasible(m);
    check_entries(m, {-0.39849977199567577, 0.1425304228867306, 0.2574695771132694,
                      -2.6015002280043245},
                  {0.2559693491089452, 0.0, 0.0, 2.3440306508910549}, 1e-11);
}

TEST_CASE("mmpp2_fit1 converts SCV, skewness and IDC into a MAP(2)") {
    // MATLAB: mmpp2_fit1(1,3,-1,4). skew = -1 defers the third moment.
    const Map2FitResult<double> r = mmpp2_fit1<double>(1.0, 3.0, -1.0, 4.0);
    REQUIRE(r.has_map);
    CHECK(r.err == 0);
    CHECK(map_moment(r.map, 1) == doctest::Approx(1.0).epsilon(kRel));
    CHECK(map_scv(r.map) == doctest::Approx(3.0).epsilon(kRel));
    // g2 = -(scv - idc)/(idc - 1) = 1/3.
    CHECK(decay(r.map, 1) == doctest::Approx(1.0 / 3.0).epsilon(kRel));
    check_entries(r.map, {-0.49966666681477934, 0.0, 0.0, -750.50033333326792},
                  {0.33288903733311193, 0.16677762948166747, 249.83322237054591,
                   500.66711096272201},
                  1e-11);
    CHECK_THROWS_AS(mmpp2_fit1<double>(1.0, 3.0, -1.0, 1.0), line::InputError);
}

TEST_CASE("mmpp2_fit2 and mmpp2_fit4 agree with mmpp2_fit3 on the same data") {
    // scv = 3, skew chosen so that E3 = 30 for E1 = 1, E2 = 4:
    // E3 = -(2 - 3*4 - skew * 3^(3/2)) => skew = (30 - 10)/3^1.5.
    const double skew = (30.0 - 10.0) / std::pow(3.0, 1.5);
    const Mmpp2FitResult<double> a = mmpp2_fit2<double>(1.0, 3.0, skew, 0.5);
    CHECK(a.feasible);
    check_moments(a.map, 1.0, 4.0, 30.0);
    CHECK(decay(a.map, 1) == doctest::Approx(0.5).epsilon(kRel));

    // mmpp2_fit4 takes the lag-1 coefficient; rho0 = (1 - 1/scv)/2 = 1/3.
    const Mmpp2FitResult<double> b = mmpp2_fit4<double>(1.0, 3.0, skew, 0.5 / 3.0);
    CHECK(b.feasible);
    check_moments(b.map, 1.0, 4.0, 30.0);
    CHECK(acf1(b.map) == doctest::Approx(0.5 / 3.0).epsilon(kRel));
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j)
            CHECK(a.map.D0(i, j) == doctest::Approx(b.map.D0(i, j)).epsilon(1e-11).scale(1.0));

    // Unit SCV degenerates to a Poisson process.
    const Mmpp2FitResult<double> c = mmpp2_fit2<double>(2.0, 1.0, 2.0, 0.5);
    CHECK(c.map.order() == 1u);
    CHECK(map_moment(c.map, 1) == doctest::Approx(2.0).epsilon(kRel));

    CHECK_THROWS_AS(mmpp2_fit4<double>(1.0, 3.0, -1.0, 0.1), line::InputError);
}

// ---------------------------------------------------------------------------
// Counting-process fits
// ---------------------------------------------------------------------------

TEST_CASE("mmpp2_fitc recovers the MMPP its characteristics came from") {
    // Characteristics of the MMPP(2) with r1 = 0.3, r2 = 0.7, l1 = 2, l2 = 0.5,
    // computed in MATLAB with map_count_mean / map_count_var / map_count_moment
    // at t1 = 1, t2 = 10, tinf = 1e8.
    const double a = 1.5499999999999998;
    const double bt1 = 1.2242877883271068;
    const double bt2 = 1.5487124453505583;
    const double binf = 1.6096774132580665;
    const double m3t2 = 27.4061949935558;
    const Mmpp2FitcResult<double> r = mmpp2_fitc<double>(a, bt1, bt2, binf, m3t2, 1.0, 10.0);
    CHECK_FALSE(r.degenerate);
    CHECK(r.third_moment_ok);
    check_feasible(r.map);
    // The fit recovers the generating process with the two phases relabelled
    // so that r1 >= r2.
    check_entries(r.map, {-1.2, 0.7, 0.3, -2.3}, {0.5, 0.0, 0.0, 2.0}, 1e-6);
    // The arrival rate is matched to full precision.
    CHECK(map_lambda(r.map) == doctest::Approx(a).epsilon(1e-7));
    // MATLAB reference (its fsolve leaves about 6e-8 on each entry).
    CHECK(r.map.D0(0, 0) == doctest::Approx(-1.1999999350962054).epsilon(1e-6));
    CHECK(r.map.D1(1, 1) == doctest::Approx(1.9999999490323699).epsilon(1e-6));
}

TEST_CASE("mmpp2_fitc degenerates to Poisson on an infeasible IDC profile") {
    // Constant unit IDC has no MMPP(2) representation.
    const Mmpp2FitcResult<double> r = mmpp2_fitc<double>(1.5, 1.0, 1.0, 1.0, 0.0, 1.0, 10.0);
    CHECK(r.degenerate);
    CHECK(r.map.order() == 1u);
    CHECK(map_lambda(r.map) == doctest::Approx(1.5).epsilon(kRel));
    // binf <= bt1 also has no representation.
    const Mmpp2FitcResult<double> s = mmpp2_fitc<double>(1.5, 2.0, 2.0, 1.5, 0.0, 1.0, 10.0);
    CHECK(s.degenerate);
    CHECK_THROWS_AS(mmpp2_fitc<double>(0.0, 1.2, 1.5, 1.6, 1.0, 1.0, 10.0), line::InputError);
}

TEST_CASE("lambertw0 inverts w exp(w) on the principal branch") {
    using line::mam::fitdetail::lambertw0;
    const double pts[] = {-0.3, -0.2, -0.05, 0.0, 0.5, 2.0, 10.0};
    for (double z : pts) {
        const double w = lambertw0<double>(z, 200u);
        CHECK(w * std::exp(w) == doctest::Approx(z).epsilon(1e-13).scale(1.0));
        CHECK(w >= -1.0);
    }
    // Right at the branch point z = -1/e the value is exactly -1.
    const double zb = -1.0 / std::exp(1.0);
    CHECK(lambertw0<double>(zb, 200u) == doctest::Approx(-1.0).epsilon(1e-7));
    CHECK_THROWS_AS(lambertw0<double>(-0.5, 200u), line::NumericError);
}

TEST_CASE("m3pp2m_fitc splits the MMPP(2) into per-class matrices") {
    const double a = 1.5499999999999998;
    const double bt1 = 1.2242877883271068;
    const double bt2 = 1.5487124453505583;
    const double binf = 1.6096774132580665;
    const double m3t2 = 27.4061949935558;
    std::vector<double> ai;
    ai.push_back(a * 0.4);
    ai.push_back(a * 0.6);
    std::vector<double> dv(2, 0.05);
    const M3pp2mFitcResult<double> r =
        m3pp2m_fitc<double>(a, bt1, bt2, binf, m3t2, 1.0, 10.0, ai, dv, 5.0);
    CHECK_FALSE(r.degenerate);
    REQUIRE(r.mmap.classes() == 2u);

    // The class matrices reconstitute D1 exactly.
    for (std::size_t i = 0; i < 2; ++i)
        for (std::size_t j = 0; j < 2; ++j)
            CHECK(r.mmap.Dc[0](i, j) + r.mmap.Dc[1](i, j) ==
                  doctest::Approx(r.mmap.D1(i, j)).epsilon(1e-12).scale(1.0));

    // MATLAB: m3pp2m_fitc(a,bt1,bt2,binf,m3t2,1,10,[0.4a;0.6a],[0.05;0.05],5).
    // dvt3 = 0.05 is outside the feasible region for this MMPP, so the first
    // class picks up a negative rate; the reference does the same, and the
    // approximate variant (skipped here: it is a quadratic program) exists
    // precisely to project such inputs back.
    CHECK(r.mmap.Dc[0](0, 0) == doctest::Approx(-0.12684861032473554).epsilon(1e-6));
    CHECK(r.mmap.Dc[0](1, 1) == doctest::Approx(0.94007790209311837).epsilon(1e-6));
    CHECK(r.mmap.Dc[1](0, 0) == doctest::Approx(0.6268484872820238).epsilon(1e-6));
    CHECK(r.mmap.Dc[1](1, 1) == doctest::Approx(1.0599220469392514).epsilon(1e-6));
}

TEST_CASE("m3pp2m_fitc marks a degenerate Poisson process by class rate") {
    std::vector<double> ai;
    ai.push_back(0.6);
    ai.push_back(0.9);
    std::vector<double> dv(2, 0.01);
    const M3pp2mFitcResult<double> r =
        m3pp2m_fitc<double>(1.5, 1.0, 1.0, 1.0, 0.0, 1.0, 10.0, ai, dv, 5.0);
    CHECK(r.degenerate);
    REQUIRE(r.mmap.classes() == 2u);
    CHECK(r.mmap.Dc[0](0, 0) == doctest::Approx(0.6).epsilon(kRel));
    CHECK(r.mmap.Dc[1](0, 0) == doctest::Approx(0.9).epsilon(kRel));

    // Class rates inconsistent with the total are rejected.
    std::vector<double> bad;
    bad.push_back(0.6);
    bad.push_back(0.1);
    CHECK_THROWS_AS(m3pp2m_fitc<double>(1.5, 1.0, 1.0, 1.0, 0.0, 1.0, 10.0, bad, dv, 5.0),
                    line::InputError);
}
