/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * The optimization-driven members of the mam fitting family:
 *   aph2_adjust_opt.h        aph2_adjust methods opt_param and opt_char
 *   amap2_adjust_gamma.h     methods 1..4
 *   mmpp2_fitc_approx.h      the optimized MMPP(2) counting fit
 *   m3pp2m_fitc_approx.h     m3pp2m_fitc_approx, _ag, _ag_multiclass
 *
 * ACCEPTANCE. None of these check an iterate against MATLAB, which would be
 * meaningless: the optimizer is a different one. They check the SPECIFICATION:
 *   - the fitted process is a valid MAP/MMAP (map_isfeasible / mmap_isfeasible,
 *     plus an explicit check that sum_i Dc_i = D1);
 *   - the characteristics the fit targets are reproduced, recomputed
 *     INDEPENDENTLY with map_count_var / map_count_moment / mmap_count_var
 *     rather than with the closed forms the objective uses;
 *   - an input that is already feasible comes back unchanged;
 *   - the achieved objective is no worse than a closed-form feasible point
 *     that is available for free (the 'simple' adjustment).
 */

#include <cmath>
#include <vector>

#include "doctest.h"
#include "line/api/mam/amap2_adjust_gamma.h"
#include "line/api/mam/amap2_fitall_gamma.h"
#include "line/api/mam/aph2_adjust.h"
#include "line/api/mam/aph2_adjust_opt.h"
#include "line/api/mam/aph2_fitall.h"
#include "line/api/mam/m3pp2m_fitc_approx.h"
#include "line/api/mam/map_count_mean.h"
#include "line/api/mam/map_count_moment.h"
#include "line/api/mam/map_count_var.h"
#include "line/api/mam/map_moment.h"
#include "line/api/mam/map_transform.h"
#include "line/api/mam/mmap_count_var.h"
#include "line/api/mam/mmpp2_fitc_approx.h"
#include "line/num/number.h"
#include "line/util/matrix.h"

using namespace line;
using namespace line::mam;

namespace {

using T = double;

/** An MMPP(2) with clearly distinct phases, used as the ground truth. */
Map<T> reference_mmpp2() {
    Map<T> m;
    m.D0 = Matrix<T>(2, 2, 0.0);
    m.D1 = Matrix<T>(2, 2, 0.0);
    const T l1 = 2.0, l2 = 0.5, r1 = 0.3, r2 = 0.7;
    m.D0(0, 0) = -(l1 + r1);
    m.D0(0, 1) = r1;
    m.D0(1, 0) = r2;
    m.D0(1, 1) = -(l2 + r2);
    m.D1(0, 0) = l1;
    m.D1(1, 1) = l2;
    return m;
}

/** IDC at window t, computed from the ported counting descriptors. */
T idc_at(const Map<T>& m, const T& t) {
    std::vector<T> tv(1, t);
    const std::vector<T> v = map_count_var(m, tv);
    const std::vector<T> mu = map_count_mean(m, tv);
    return v[0] / mu[0];
}

/** Third CENTRAL moment of the counts at window t. */
T m3_counts_at(const Map<T>& m, const T& t) {
    std::vector<unsigned> ord(3);
    ord[0] = 1;
    ord[1] = 2;
    ord[2] = 3;
    const std::vector<T> mt = map_count_moment(m, t, ord);
    return mt[2] - 3.0 * mt[1] * mt[0] + 2.0 * mt[0] * mt[0] * mt[0];
}

/** sum_i Dc_i - D1, the defining identity of a marked MAP. */
T mark_defect(const Mmap<T>& mm) {
    T worst = 0.0;
    for (std::size_t r = 0; r < mm.D1.rows(); ++r)
        for (std::size_t c = 0; c < mm.D1.cols(); ++c) {
            T s = 0.0;
            for (std::size_t k = 0; k < mm.Dc.size(); ++k) s += mm.Dc[k](r, c);
            const T d = std::fabs(s - mm.D1(r, c));
            if (d > worst) worst = d;
        }
    return worst;
}

}  // namespace

// ---------------------------------------------------------------------------
// aph2_adjust, methods opt_param and opt_char
// ---------------------------------------------------------------------------

TEST_CASE("aph2_adjust opt_param leaves an already feasible moment set alone") {
    // an Erlang-2-like triple, comfortably inside the APH(2) region
    const T M1 = 1.0;
    const T M2 = 1.6;   // scv = 0.6, above the 1/2 floor
    const T M3 = 3.5;   // inside [3.4735, 3.6], the APH(2) window for that scv
    REQUIRE(fitdetail::aph2_moments_feasible(M1, M2, M3, T(0.0)));

    const Aph2AdjustOptResult<T> r = aph2_adjust_opt_param(M1, M2, M3);
    CHECK(r.converged);
    CHECK(r.objective < 1e-6);
    CHECK(r.M2a == doctest::Approx(M2).epsilon(1e-6));
    CHECK(r.M3a == doctest::Approx(M3).epsilon(1e-6));
    CHECK(fitdetail::aph2_moments_feasible(M1, r.M2a, r.M3a, T(1e-6)));
}

TEST_CASE("aph2_adjust opt_param repairs an SCV below the APH(2) floor") {
    // scv = 0.3 < 1/2 is not representable by any APH(2)
    const T M1 = 1.0;
    const T M2 = 1.3;
    const T M3 = 2.0;
    REQUIRE(!fitdetail::aph2_moments_feasible(M1, M2, M3, T(0.0)));

    const Aph2AdjustOptResult<T> r = aph2_adjust_opt_param(M1, M2, M3);
    CHECK(r.converged);
    CHECK(fitdetail::aph2_moments_feasible(M1, r.M2a, r.M3a, T(1e-6)));

    // an APH(2) with the adjusted moments actually exists
    const std::vector<Map<T>> fits = aph2_fitall(M1, r.M2a, r.M3a);
    CHECK(!fits.empty());
    CHECK(map_moment(fits[0], 1) == doctest::Approx(M1).epsilon(1e-6));
    CHECK(map_moment(fits[0], 2) == doctest::Approx(r.M2a).epsilon(1e-6));

    // and the adjustment is no larger than the closed-form 'simple' one,
    // which is a feasible point available for free
    const Aph2AdjustResult<T> s = aph2_adjust(M1, M2, M3);
    const T dsimple =
        std::sqrt((s.M2a - M2) * (s.M2a - M2) + (s.M3a - M3) * (s.M3a - M3));
    CHECK(r.objective <= dsimple * (1.0 + 1e-6));
}

TEST_CASE("aph2_adjust opt_char leaves an already feasible moment set alone") {
    const T M1 = 1.0;
    const T M2 = 1.6;
    const T M3 = 3.5;
    const Aph2AdjustOptResult<T> r = aph2_adjust_opt_char(M1, M2, M3);
    CHECK(r.converged);
    CHECK(r.objective < 1e-6);
    CHECK(r.M2a == doctest::Approx(M2).epsilon(1e-6));
    CHECK(r.M3a == doctest::Approx(M3).epsilon(1e-6));
}

TEST_CASE("aph2_adjust opt_char repairs an infeasible third moment") {
    // scv = 0.6 is fine but M3 is far below the lower bound for that scv
    const T M1 = 1.0;
    const T M2 = 1.6;
    const T M3 = 2.0;
    REQUIRE(!fitdetail::aph2_moments_feasible(M1, M2, M3, T(0.0)));

    const Aph2AdjustOptResult<T> r = aph2_adjust_opt_char(M1, M2, M3);
    CHECK(r.converged);
    CHECK(fitdetail::aph2_moments_feasible(M1, r.M2a, r.M3a, T(1e-5)));
    const std::vector<Map<T>> fits = aph2_fitall(M1, r.M2a, r.M3a);
    CHECK(!fits.empty());

    const Aph2AdjustResult<T> s = aph2_adjust(M1, M2, M3);
    const T dsimple =
        std::sqrt((s.M2a - M2) * (s.M2a - M2) + (s.M3a - M3) * (s.M3a - M3));
    CHECK(r.objective <= dsimple * (1.0 + 1e-6));
}

TEST_CASE("the two aph2_adjust optimization methods agree on the feasible region") {
    // both must return a feasible pair; opt_char, which searches the moment
    // space directly, cannot do worse than opt_param on the same input
    const T M1 = 2.0;
    const T M2 = 9.0;   // scv = 1.25 > 1
    const T M3 = 30.0;  // below the lower bound 3/2 M1^3 (1 + scv)^2 = 60.75
    const Aph2AdjustOptResult<T> a = aph2_adjust_opt_param(M1, M2, M3);
    const Aph2AdjustOptResult<T> b = aph2_adjust_opt_char(M1, M2, M3);
    CHECK(fitdetail::aph2_moments_feasible(M1, a.M2a, a.M3a, T(1e-5)));
    CHECK(fitdetail::aph2_moments_feasible(M1, b.M2a, b.M3a, T(1e-5)));
    CHECK(b.objective <= a.objective * (1.0 + 1e-4));
}

TEST_CASE("aph2_adjust optimization methods reject a non-positive mean") {
    CHECK_THROWS_AS(aph2_adjust_opt_param<T>(0.0, 1.0, 1.0), InputError);
    CHECK_THROWS_AS(aph2_adjust_opt_char<T>(-1.0, 1.0, 1.0), InputError);
}

// ---------------------------------------------------------------------------
// amap2_adjust_gamma
// ---------------------------------------------------------------------------

TEST_CASE("amap2_adjust_gamma method 3 is the closed-form priority rule") {
    const T M1 = 1.0;
    const T M2 = 3.0;   // scv = 2
    const T M3 = 30.0;
    const T GAMMA = 0.5;
    const Amap2AdjustGammaResult<T> r = amap2_adjust_gamma(M1, M2, M3, GAMMA);
    // M2 and M3 come straight from the 'simple' APH(2) adjustment
    const Aph2AdjustResult<T> s = aph2_adjust(M1, M2, M3);
    CHECK(r.M2a == doctest::Approx(s.M2a));
    CHECK(r.M3a == doctest::Approx(s.M3a));
    CHECK(r.feasible);
    // and an AMAP(2) with those characteristics exists
    const std::vector<Map<T>> all = amap2_fitall_gamma(M1, r.M2a, r.M3a, r.GAMMAa);
    CHECK(!all.empty());
}

TEST_CASE("amap2_adjust_gamma clamps an out-of-range decay rate") {
    const T M1 = 1.0;
    const T M2 = 3.0;
    const T M3 = 30.0;
    const T GAMMA = 0.999;  // above 1 - tol with tol = 1e-2
    const Amap2AdjustGammaResult<T> r = amap2_adjust_gamma(M1, M2, M3, GAMMA);
    CHECK(r.GAMMAa <= 0.99 + 1e-12);
    CHECK(r.feasible);
}

TEST_CASE("amap2_adjust_gamma methods 1, 2 and 4 all return feasible characteristics") {
    const T M1 = 1.0;
    const T M2 = 3.0;
    const T M3 = 30.0;
    const T GAMMA = 0.5;
    for (int method = 1; method <= 4; ++method) {
        const Amap2AdjustGammaResult<T> r = amap2_adjust_gamma(M1, M2, M3, GAMMA, method);
        INFO("method ", method);
        CHECK(r.feasible);
        const std::vector<Map<T>> all = amap2_fitall_gamma(M1, r.M2a, r.M3a, r.GAMMAa);
        CHECK(!all.empty());
    }
}

TEST_CASE("amap2_adjust_gamma leaves a feasible triple essentially unchanged") {
    // characteristics read off an actual AMAP(2), so nothing needs adjusting
    const T M1 = 1.0;
    const Map<T> amap = map_scale(amap2_assemble<T>(1.0, 1.0 / 3.0, 0.5, 2.0 / 3.0, 1), M1);
    const T M2 = map_moment(amap, 2);
    const T M3 = map_moment(amap, 3);
    std::vector<unsigned> lags(2);
    lags[0] = 3;
    lags[1] = 4;
    const std::vector<T> acf = map_acf(amap, lags);
    const T GAMMA = acf[1] / acf[0];

    // The reference keeps every characteristic strictly inside its region by a
    // slack tol = 1e-2, so a triple sitting ON the boundary -- this GAMMA is
    // 0.9955, above the upper bound 1 - tol -- is legitimately pulled in by up
    // to that slack. The invariance claim is therefore "within the reference's
    // own strict-inequality slack", not "bitwise unchanged".
    for (int method = 1; method <= 4; ++method) {
        const Amap2AdjustGammaResult<T> r = amap2_adjust_gamma(M1, M2, M3, GAMMA, method);
        INFO("method ", method);
        CHECK(r.feasible);
        CHECK(std::fabs(r.M2a / M2 - 1.0) < 1e-2);
        CHECK(std::fabs(r.M3a / M3 - 1.0) < 1e-2);
        CHECK(std::fabs(r.GAMMAa / GAMMA - 1.0) < 1e-2);
        const std::vector<Map<T>> all = amap2_fitall_gamma(M1, r.M2a, r.M3a, r.GAMMAa);
        CHECK(!all.empty());
    }
}

TEST_CASE("amap2_adjust_gamma is deterministic across repeated calls") {
    const T M1 = 1.0, M2 = 3.0, M3 = 30.0, GAMMA = 0.5;
    for (int method = 1; method <= 4; ++method) {
        const Amap2AdjustGammaResult<T> a = amap2_adjust_gamma(M1, M2, M3, GAMMA, method);
        const Amap2AdjustGammaResult<T> b = amap2_adjust_gamma(M1, M2, M3, GAMMA, method);
        INFO("method ", method);
        CHECK(a.M2a == b.M2a);
        CHECK(a.M3a == b.M3a);
        CHECK(a.GAMMAa == b.GAMMAa);
    }
}

TEST_CASE("amap2_adjust_gamma refuses the unported constraint set and bad methods") {
    std::vector<T> w(3, 1.0);
    CHECK_THROWS_AS(amap2_adjust_gamma<T>(1.0, 3.0, 30.0, 0.5, w, 1, 1, 1e-2), UnsupportedError);
    CHECK_THROWS_AS(amap2_adjust_gamma<T>(1.0, 3.0, 30.0, 0.5, w, 5, 2, 1e-2), InputError);
    CHECK_THROWS_AS(amap2_adjust_gamma<T>(1.0, 3.0, 30.0, 0.0, w, 1, 2, 1e-2), InputError);
}

// ---------------------------------------------------------------------------
// mmpp2_fitc_approx
// ---------------------------------------------------------------------------

TEST_CASE("mmpp2_fitc_approx recovers the characteristics of an actual MMPP(2)") {
    const Map<T> src = reference_mmpp2();
    const T t1 = 1.0, t2 = 5.0, t3 = 1.0;
    const T a = map_lambda(src);
    const T bt1 = idc_at(src, t1);
    const T bt2 = idc_at(src, t2);
    const T binf = idc_at(src, T(1e6));
    const T m3t2 = m3_counts_at(src, t2);
    (void)t3;

    const Mmpp2FitcApproxResult<T> r = mmpp2_fitc_approx(a, bt1, bt2, binf, m3t2, t1, t2);

    // 1. a valid MAP
    CHECK(map_isfeasible(r.map, T(1e-10)));
    // 2. the rate is exact by construction
    CHECK(map_lambda(r.map) == doctest::Approx(a).epsilon(1e-12));
    // 3. the objective reached zero: the target IS attainable
    CHECK(r.objective < 1e-12);
    // 4. the characteristics, recomputed independently, reproduce the targets
    CHECK(idc_at(r.map, t1) == doctest::Approx(bt1).epsilon(1e-6));
    CHECK(idc_at(r.map, t2) == doctest::Approx(bt2).epsilon(1e-6));
    CHECK(idc_at(r.map, T(1e6)) == doctest::Approx(binf).epsilon(1e-6));
    CHECK(m3_counts_at(r.map, t2) == doctest::Approx(m3t2).epsilon(1e-5));
}

TEST_CASE("mmpp2_fitc_approx is deterministic and rejects degenerate targets") {
    const Map<T> src = reference_mmpp2();
    const T t1 = 1.0, t2 = 5.0;
    const T a = map_lambda(src);
    const T bt1 = idc_at(src, t1);
    const T bt2 = idc_at(src, t2);
    const T binf = idc_at(src, T(1e6));
    const T m3t2 = m3_counts_at(src, t2);

    const Mmpp2FitcApproxResult<T> x = mmpp2_fitc_approx(a, bt1, bt2, binf, m3t2, t1, t2);
    const Mmpp2FitcApproxResult<T> y = mmpp2_fitc_approx(a, bt1, bt2, binf, m3t2, t1, t2);
    CHECK(x.objective == y.objective);
    CHECK(x.map.D1(0, 0) == y.map.D1(0, 0));

    CHECK_THROWS_AS(mmpp2_fitc_approx<T>(0.0, bt1, bt2, binf, m3t2, t1, t2), InputError);
    CHECK_THROWS_AS(mmpp2_fitc_approx<T>(a, bt1, bt2, 0.0, m3t2, t1, t2), InputError);
    CHECK_THROWS_AS(mmpp2_fitc_approx<T>(a, bt1, bt2, binf, m3t2, T(-1.0), t2), InputError);
}

TEST_CASE("mmpp2_fitc_approx accepts a single time scale, which the reference cannot") {
    // matlab/lib/kpctoolbox/mmpp/mmpp2_fitc_approx.m leaves xbt2 unassigned
    // when t1 == t2 and then uses it, so the reference errors out here.
    const Map<T> src = reference_mmpp2();
    const T t = 2.0;
    const T a = map_lambda(src);
    const T bt = idc_at(src, t);
    const T binf = idc_at(src, T(1e6));
    const T m3t = m3_counts_at(src, t);
    const Mmpp2FitcApproxResult<T> r = mmpp2_fitc_approx(a, bt, bt, binf, m3t, t, t);
    CHECK(map_isfeasible(r.map, T(1e-10)));
    CHECK(r.objective < 1e-12);
    CHECK(idc_at(r.map, t) == doctest::Approx(bt).epsilon(1e-6));
}

// ---------------------------------------------------------------------------
// m3pp2m_fitc_approx and the 'ag' variants
// ---------------------------------------------------------------------------

TEST_CASE("m3pp2m_fitc_approx recovers the per-class variance differences") {
    const Map<T> src = reference_mmpp2();
    const T t1 = 1.0, t2 = 5.0, t3 = 1.0;
    const T a = map_lambda(src);
    const T bt1 = idc_at(src, t1);
    const T bt2 = idc_at(src, t2);
    const T binf = idc_at(src, T(1e6));
    const T m3t2 = m3_counts_at(src, t2);

    // ground truth: a three-class marking of the source MMPP, with per-phase
    // probabilities that sum to one in each phase
    const T q1[3] = {0.5, 0.3, 0.2};
    const T q2[3] = {0.2, 0.5, 0.3};
    std::vector<T> qq1(q1, q1 + 3), qq2(q2, q2 + 3);
    const Mmap<T> truth = m3pp2m_assemble(src, qq1, qq2);
    REQUIRE(mark_defect(truth) < 1e-14);

    const std::vector<T> ai = mmap_count_lambda(truth);
    std::vector<T> dvt3(3);
    for (std::size_t i = 0; i < 3; ++i) {
        Mmap<T> two;
        two.D0 = truth.D0;
        two.D1 = truth.D1;
        Matrix<T> rest(2, 2, 0.0);
        for (std::size_t r = 0; r < 2; ++r)
            for (std::size_t c = 0; c < 2; ++c) rest(r, c) = truth.D1(r, c) - truth.Dc[i](r, c);
        two.Dc.push_back(truth.Dc[i]);
        two.Dc.push_back(rest);
        const std::vector<T> v = mmap_count_var(two, t3);
        dvt3[i] = v[0] - v[1];
    }

    const M3pp2mFitcApproxResult<T> r =
        m3pp2m_fitc_approx(a, bt1, bt2, binf, m3t2, t1, t2, ai, dvt3, t3);

    CHECK(r.mmpp_objective < 1e-12);
    CHECK(mark_defect(r.mmap) < 1e-10);
    CHECK(r.feasible);
    CHECK(r.class_violation < 1e-8);
    // the target IS attainable, so the per-class objective reaches zero
    CHECK(r.class_objective < 1e-12);

    // per-class rates matched exactly
    const std::vector<T> ai_fit = mmap_count_lambda(r.mmap);
    for (std::size_t i = 0; i < 3; ++i) CHECK(ai_fit[i] == doctest::Approx(ai[i]).epsilon(1e-8));

    // and the variance differences, recomputed independently
    for (std::size_t i = 0; i < 3; ++i) {
        Mmap<T> two;
        two.D0 = r.mmap.D0;
        two.D1 = r.mmap.D1;
        Matrix<T> rest(2, 2, 0.0);
        for (std::size_t p = 0; p < 2; ++p)
            for (std::size_t c = 0; c < 2; ++c) rest(p, c) = r.mmap.D1(p, c) - r.mmap.Dc[i](p, c);
        two.Dc.push_back(r.mmap.Dc[i]);
        two.Dc.push_back(rest);
        const std::vector<T> v = mmap_count_var(two, t3);
        CHECK(v[0] - v[1] == doctest::Approx(dvt3[i]).epsilon(1e-5));
    }
}

TEST_CASE("m3pp2m_fitc_approx_ag_multiclass recovers the per-class variance-plus-covariance") {
    const Map<T> src = reference_mmpp2();
    const T t3 = 1.0;
    const T q1[3] = {0.5, 0.3, 0.2};
    const T q2[3] = {0.2, 0.5, 0.3};
    std::vector<T> qq1(q1, q1 + 3), qq2(q2, q2 + 3);
    const Mmap<T> truth = m3pp2m_assemble(src, qq1, qq2);
    const std::vector<T> ai = mmap_count_lambda(truth);

    // g_i = Var[N_i] + Cov[N_i, N_rest] = (V_i - V_rest + V_total)/2
    std::vector<T> tv(1, t3);
    const T vtot = map_count_var(src, tv)[0];
    std::vector<T> gt3(3);
    for (std::size_t i = 0; i < 3; ++i) {
        Mmap<T> two;
        two.D0 = truth.D0;
        two.D1 = truth.D1;
        Matrix<T> rest(2, 2, 0.0);
        for (std::size_t p = 0; p < 2; ++p)
            for (std::size_t c = 0; c < 2; ++c) rest(p, c) = truth.D1(p, c) - truth.Dc[i](p, c);
        two.Dc.push_back(truth.Dc[i]);
        two.Dc.push_back(rest);
        const std::vector<T> v = mmap_count_var(two, t3);
        gt3[i] = (v[0] - v[1] + vtot) / 2.0;
    }

    const M3pp2mFitcApproxResult<T> r = m3pp2m_fitc_approx_ag_multiclass(src, ai, gt3, t3);
    CHECK(mark_defect(r.mmap) < 1e-10);
    CHECK(r.feasible);
    CHECK(r.class_violation < 1e-8);
    CHECK(r.class_objective < 1e-12);

    const std::vector<T> ai_fit = mmap_count_lambda(r.mmap);
    for (std::size_t i = 0; i < 3; ++i) CHECK(ai_fit[i] == doctest::Approx(ai[i]).epsilon(1e-8));

    for (std::size_t i = 0; i < 3; ++i) {
        Mmap<T> two;
        two.D0 = r.mmap.D0;
        two.D1 = r.mmap.D1;
        Matrix<T> rest(2, 2, 0.0);
        for (std::size_t p = 0; p < 2; ++p)
            for (std::size_t c = 0; c < 2; ++c) rest(p, c) = r.mmap.D1(p, c) - r.mmap.Dc[i](p, c);
        two.Dc.push_back(r.mmap.Dc[i]);
        two.Dc.push_back(rest);
        const std::vector<T> v = mmap_count_var(two, t3);
        const T g = (v[0] - v[1] + vtot) / 2.0;
        CHECK(g == doctest::Approx(gt3[i]).epsilon(1e-5));
    }
}

TEST_CASE("m3pp2m_fitc_approx_ag chains the MMPP fit and the ag split") {
    const Map<T> src = reference_mmpp2();
    const T t1 = 1.0, t2 = 5.0, t3 = 1.0;
    const T a = map_lambda(src);
    const T bt1 = idc_at(src, t1);
    const T bt2 = idc_at(src, t2);
    const T binf = idc_at(src, T(1e6));
    const T m3t2 = m3_counts_at(src, t2);

    const T q1[2] = {0.6, 0.4};
    const T q2[2] = {0.25, 0.75};
    std::vector<T> qq1(q1, q1 + 2), qq2(q2, q2 + 2);
    const Mmap<T> truth = m3pp2m_assemble(src, qq1, qq2);
    const std::vector<T> ai = mmap_count_lambda(truth);

    std::vector<T> tv(1, t3);
    const T vtot = map_count_var(src, tv)[0];
    std::vector<T> gt3(2);
    for (std::size_t i = 0; i < 2; ++i) {
        Mmap<T> two;
        two.D0 = truth.D0;
        two.D1 = truth.D1;
        Matrix<T> rest(2, 2, 0.0);
        for (std::size_t p = 0; p < 2; ++p)
            for (std::size_t c = 0; c < 2; ++c) rest(p, c) = truth.D1(p, c) - truth.Dc[i](p, c);
        two.Dc.push_back(truth.Dc[i]);
        two.Dc.push_back(rest);
        const std::vector<T> v = mmap_count_var(two, t3);
        gt3[i] = (v[0] - v[1] + vtot) / 2.0;
    }

    const M3pp2mFitcApproxResult<T> r =
        m3pp2m_fitc_approx_ag(a, bt1, bt2, binf, m3t2, t1, t2, ai, gt3, t3);
    CHECK(r.mmpp_objective < 1e-12);
    CHECK(r.class_objective < 1e-10);
    CHECK(mark_defect(r.mmap) < 1e-10);
    CHECK(r.feasible);
}

TEST_CASE("the M3PP splits short-circuit on one class and on a Poisson process") {
    const Map<T> src = reference_mmpp2();
    const T a = map_lambda(src);
    std::vector<T> ai(1, a);
    std::vector<T> gt3(1, 1.0);
    const M3pp2mFitcApproxResult<T> one = m3pp2m_fitc_approx_ag_multiclass(src, ai, gt3, T(1.0));
    CHECK(one.mmap.classes() == 1);
    CHECK(mark_defect(one.mmap) < 1e-14);
    CHECK(one.class_objective == 0.0);

    const Map<T> poi = map_exponential<T>(2.0);
    std::vector<T> ai2(2);
    ai2[0] = 0.5;
    ai2[1] = 1.5;
    std::vector<T> gt32(2, 1.0);
    const M3pp2mFitcApproxResult<T> deg = m3pp2m_fitc_approx_ag_multiclass(poi, ai2, gt32, T(1.0));
    CHECK(deg.degenerate);
    CHECK(deg.mmap.classes() == 2);
    CHECK(mark_defect(deg.mmap) < 1e-14);
    const std::vector<T> lam = mmap_count_lambda(deg.mmap);
    CHECK(lam[0] == doctest::Approx(0.5));
    CHECK(lam[1] == doctest::Approx(1.5));
}

TEST_CASE("the M3PP splits reject inconsistent per-class rates and zero targets") {
    const Map<T> src = reference_mmpp2();
    std::vector<T> ai(2);
    ai[0] = 0.1;
    ai[1] = 0.1;  // does not sum to the rate of src
    std::vector<T> gt3(2, 1.0);
    CHECK_THROWS_AS(m3pp2m_fitc_approx_ag_multiclass(src, ai, gt3, T(1.0)), InputError);

    const std::vector<T> good = mmap_count_lambda(m3pp2m_assemble(
        src, std::vector<T>(2, 0.5), std::vector<T>(2, 0.5)));
    std::vector<T> zero_target(2, 0.0);
    CHECK_THROWS_AS(m3pp2m_fitc_approx_ag_multiclass(src, good, zero_target, T(1.0)), InputError);
}
