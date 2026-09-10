/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_M3PP2M_INTERLEAVE_H
#define LINE_API_MAM_M3PP2M_INTERLEAVE_H

/**
 * LUMPED interleaving of several M3PP(2, m), and the two fitters built on it
 * (matlab/lib/m3a/m3a/m3pp/m3pp2m_interleave.m,
 *  matlab/lib/m3a/m3a/m3pp/m3pp2m_fitc_theoretical.m,
 *  matlab/lib/m3a/m3a/m3pp/m3pp22_interleave_fitc.m).
 *
 * Interleaving is the cheap alternative to superposition. Superposing L
 * two-phase processes gives the PRODUCT chain, of order 2^L; interleaving
 * instead lays the L phase processes on a single BIRTH-DEATH chain of order
 * L + 1, where phase h means "the first h components are in their fast phase".
 * The off-diagonal rates are recovered by DIFFERENCING: the upper rate out of
 * level j is component j's r1 minus the rates already spent on levels above it,
 * and symmetrically downwards. That differencing is only meaningful when the
 * components' rates are ordered, which is what m3pp22_interleave_fitc's linear
 * program arranges before it ever fits a component -- it does not fit L
 * processes and then hope they interleave, it SOLVES for off-diagonal rates
 * that admit the interleaving and fits the components to those.
 *
 * The class matrices carry no cross terms: class j of component i fires at its
 * phase-1 rate on levels h <= i and at its phase-2 rate above, so sum_c Dc = D1
 * holds level by level.
 *
 * THE LINEAR PROGRAM IS A FEASIBILITY PROBLEM. The reference passes a ZERO
 * objective to linprog, so every feasible point is optimal and which vertex
 * comes back is the solver's choice; different LP backends hand different
 * MMPP(2)s to the per-pair covariance split, and a covariance one accepts
 * another can report infeasible. That is a property of the reference, not of
 * this port, and it is why m3pp22_interleave_fitc reports the realised
 * covariance per pair rather than asserting the requested one.
 *
 * Gated on transcendental arithmetic, and on double alone for the LP.
 */

#include <cstddef>
#include <string>
#include <vector>

#include "line/api/mam/m3pp22_fitc_cov.h"
#include "line/api/mam/m3pp2m_fitc.h"
#include "line/api/mam/m3pp2m_fitc_approx.h"
#include "line/api/mam/m3pp_superpos_fitc.h"
#include "line/api/mam/map_count_mean.h"
#include "line/api/mam/map_count_moment.h"
#include "line/api/mam/map_count_var.h"
#include "line/api/mam/mmap_count_var.h"
#include "line/api/mam/mmap_lambda.h"
#include "line/api/mam/mmap_stats.h"
#include "line/api/mam/mmpp2_fitc.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"
#include "line/util/simplex.h"

namespace line {
namespace mam {

/**
 * Interleave L M3PP(2, m_i) into one M3PP of order L + 1 whose class list is
 * the concatenation of theirs.
 */
template <class T>
Mmap<T> m3pp2m_interleave(const std::vector<Mmap<T>>& parts) {
    const T zero = num_traits<T>::from_int(0);
    const std::size_t L = parts.size();
    if (L == 0) throw InputError("m3pp2m_interleave: no components");
    for (std::size_t i = 0; i < L; ++i)
        if (parts[i].order() != 2)
            throw InputError("m3pp2m_interleave: every component must have order 2");

    // r(0, i): upper off-diagonal rate contributed by component i, differenced
    // from the top down; r(1, i): lower rate, differenced from the bottom up.
    std::vector<T> r0(L, zero), r1(L, zero);
    r0[L - 1] = parts[L - 1].D0(0, 1);
    for (std::size_t k = L - 1; k-- > 0;) {
        T acc = zero;
        for (std::size_t j = k + 1; j < L; ++j) acc += r0[j];
        r0[k] = parts[k].D0(0, 1) - acc;
    }
    r1[0] = parts[0].D0(1, 0);
    for (std::size_t i = 1; i < L; ++i) {
        T acc = zero;
        for (std::size_t j = 0; j < i; ++j) acc += r1[j];
        r1[i] = parts[i].D0(1, 0) - acc;
    }

    std::size_t M = 0;
    for (std::size_t i = 0; i < L; ++i) M += parts[i].classes();
    const std::size_t n = 2 + (L - 1);

    Mmap<T> s;
    s.D0 = Matrix<T>(n, n, zero);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) {
            if (j > i)
                s.D0(i, j) = r0[j - 1];
            else if (j < i)
                s.D0(i, j) = r1[j];
        }

    for (std::size_t i = 0; i < L; ++i) {
        for (std::size_t c = 0; c < parts[i].classes(); ++c) {
            Matrix<T> Dc(n, n, zero);
            for (std::size_t h = 0; h < n; ++h)
                Dc(h, h) = h <= i ? parts[i].Dc[c](0, 0) : parts[i].Dc[c](1, 1);
            s.Dc.push_back(Dc);
        }
    }

    s.D1 = Matrix<T>(n, n, zero);
    for (std::size_t c = 0; c < M; ++c)
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) s.D1(i, j) += s.Dc[c](i, j);

    for (std::size_t h = 0; h < n; ++h) {
        T acc = zero;
        for (std::size_t j = 0; j < n; ++j) acc += s.D0(h, j) + s.D1(h, j);
        s.D0(h, h) = -acc;
    }
    return s;
}

/**
 * Fit the counting characteristics of a GIVEN MMAP with an M3PP(2, m).
 *
 * @param mm     the MMAP(n, m) to fit
 * @param method 'exact_delta', 'approx_delta', 'approx_cov' or 'approx_ag'
 * @param t      the single finite time scale (the reference sets t1 = t2 = t3 = t)
 * @param tinf   the near-infinite time scale
 */
template <class T>
Mmap<T> m3pp2m_fitc_theoretical(const Mmap<T>& mm, const std::string& method, const T& t,
                                const T& tinf) {
    static_assert(num_traits<T>::has_transcendental,
                  "m3pp2m_fitc_theoretical requires transcendental arithmetic");
    const std::size_t m = mm.classes();
    if (m == 0) throw InputError("m3pp2m_fitc_theoretical: the MMAP has no classes");
    if (method == "approx_cov" && m > 2)
        throw InputError("m3pp2m_fitc_theoretical: approximate covariance fitting only supports "
                         "two classes");

    const T t1 = t, t2 = t, t3 = t;
    const Map<T> joint = mm.map();

    std::vector<T> ts;
    ts.push_back(t1);
    ts.push_back(t2);
    ts.push_back(tinf);
    ts.push_back(t3);
    const std::vector<T> mu = map_count_mean(joint, ts);
    const std::vector<T> vr = map_count_var(joint, ts);

    const T a = mu[0] / t1;
    const T bt1 = vr[0] / (a * t1);
    const T bt2 = vr[1] / (a * t2);
    const T binf = vr[2] / (a * tinf);

    std::vector<unsigned> orders;
    orders.push_back(1);
    orders.push_back(2);
    orders.push_back(3);
    const std::vector<T> mt2 = map_count_moment(joint, t2, orders);
    const T m3t2 = fitdetail::m3_from_raw(mt2[0], mt2[1], mt2[2]);

    const std::vector<T> ai = mmap_count_mean(mm, num_traits<T>::from_int(1));

    if (method == "exact_delta" || method == "approx_delta") {
        std::vector<T> dvt3(m);
        for (std::size_t i = 0; i < m; ++i) {
            Mmap<T> two;
            two.D0 = mm.D0;
            two.D1 = mm.D1;
            two.Dc.push_back(mm.Dc[i]);
            Matrix<T> rest = mm.D1;
            for (std::size_t r = 0; r < rest.rows(); ++r)
                for (std::size_t c = 0; c < rest.cols(); ++c) rest(r, c) -= mm.Dc[i](r, c);
            two.Dc.push_back(rest);
            const std::vector<T> V = mmap_count_var(two, t3);
            dvt3[i] = V[0] - V[1];
        }
        if (method == "exact_delta")
            return m3pp2m_fitc(a, bt1, bt2, binf, m3t2, t1, t2, ai, dvt3, t3).mmap;
        return m3pp2m_fitc_approx(a, bt1, bt2, binf, m3t2, t1, t2, ai, dvt3, t3).mmap;
    }
    if (method == "approx_cov") {
        const std::vector<T> vi = mmap_count_var(mm, t3);
        T sum = num_traits<T>::from_int(0);
        for (std::size_t i = 0; i < vi.size(); ++i) sum += vi[i];
        const T s = (vr[3] - sum) / num_traits<T>::from_int(2);
        return m3pp22_fitc_approx_cov(a, bt1, bt2, binf, m3t2, t1, t2, ai, s, t3).mmap;
    }
    if (method == "approx_ag") {
        std::vector<T> gt3(m);
        for (std::size_t i = 0; i < m; ++i) {
            Mmap<T> two;
            two.D0 = mm.D0;
            two.D1 = mm.D1;
            two.Dc.push_back(mm.Dc[i]);
            Matrix<T> rest = mm.D1;
            for (std::size_t r = 0; r < rest.rows(); ++r)
                for (std::size_t c = 0; c < rest.cols(); ++c) rest(r, c) -= mm.Dc[i](r, c);
            two.Dc.push_back(rest);
            const std::vector<T> V = mmap_count_var(two, t3);
            const Matrix<T> S = mmap_count_mcov(two, t3);
            gt3[i] = V[0] + S(0, 1);
        }
        return m3pp2m_fitc_approx_ag(a, bt1, bt2, binf, m3t2, t1, t2, ai, gt3, t3).mmap;
    }
    throw InputError("m3pp2m_fitc_theoretical: unknown method '" + method + "'");
}

/** m3pp2m_fitc_theoretical with the reference's default scales t = 10, tinf = 1e4. */
template <class T>
Mmap<T> m3pp2m_fitc_theoretical(const Mmap<T>& mm, const std::string& method) {
    return m3pp2m_fitc_theoretical(mm, method, num_traits<T>::from_int(10),
                                   num_traits<T>::from_double(1e4));
}

/** Result of m3pp22_interleave_fitc. */
template <class T>
struct M3pp22InterleaveResult {
    Mmap<T> mmap;                          ///< the lumped process, of order L + 1
    std::vector<Mmap<T>> parts;            ///< the L M3PP(2, 2) components
    std::vector<T> sigma;                  ///< the covariance realised for each pair
};

/**
 * Fit L PAIRS of classes into one MMAP by lumped interleaving of L M3PP(2, 2).
 *
 * @param av    (L x 2) per-class rates
 * @param btv   per-pair IDC at t
 * @param binfv per-pair asymptotic IDC
 * @param stv   per-pair count covariance at t
 * @param t     the time scale
 */
template <class T>
M3pp22InterleaveResult<T> m3pp22_interleave_fitc(const Matrix<T>& av, const std::vector<T>& btv,
                                                 const std::vector<T>& binfv,
                                                 const std::vector<T>& stv, const T& t) {
    static_assert(num_traits<T>::has_transcendental,
                  "m3pp22_interleave_fitc requires transcendental arithmetic");
    using fitdetail::lambertw0;
    using fitdetail::num_exp;
    using fitdetail::num_sqrt;
    using fitdetail::pw;

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const std::size_t L = av.rows();
    if (L == 0) throw InputError("m3pp22_interleave_fitc: no pairs");
    if (av.cols() != 2) throw InputError("m3pp22_interleave_fitc: av must have two columns");
    if (btv.size() != L || binfv.size() != L || stv.size() != L)
        throw InputError("m3pp22_interleave_fitc: btv, binfv and stv must have one entry per pair");

    // bounds on the upper off-diagonal element of each MMPP(2)
    std::vector<T> uv(L, zero), dv(L, zero);
    for (std::size_t i = 0; i < L; ++i) {
        const T a = av(i, 0) + av(i, 1);
        if (!(binfv[i] > btv[i] && btv[i] > one))
            throw InputError("m3pp22_interleave_fitc: infeasible IDC pair, IDC(inf) must exceed "
                             "IDC(t) and IDC(t) must exceed one");
        const T c = (binfv[i] - one) / (binfv[i] - btv[i]);
        const T w = lambertw0(T(-c * num_exp(T(-c))), 200u);
        const T d = (w + c) / t;
        const T z = (binfv[i] - one) * pw(d, 3) * a;
        uv[i] = d * z / (two * a * a * d * d + z);
        dv[i] = d;
    }

    // feasibility LP over the 2L per-level off-diagonal rates
    lp::LpModel<T> model(2 * L);
    model.set_maximize(false);
    for (std::size_t j = 0; j < 2 * L; ++j) {
        model.set_lower(j, zero);
        model.set_free_upper(j);
        model.set_cost(j, zero);
    }
    const T eps = num_traits<T>::from_double(1e-6);
    for (std::size_t i = 0; i < L; ++i) {
        model.row_clear();
        for (std::size_t j = i; j < L; ++j) model.row_add(j, one);
        model.emit_le(T(dv[i] - eps));
        model.row_clear();
        for (std::size_t j = i; j < L; ++j) model.row_add(j, T(-one));
        model.emit_le(T(-uv[i] - eps));
    }
    for (std::size_t i = 0; i < L; ++i) {
        model.row_clear();
        for (std::size_t j = i; j < L; ++j) model.row_add(j, one);
        for (std::size_t j = 0; j <= i; ++j) model.row_add(L + j, one);
        model.emit_eq(dv[i]);
    }
    const lp::LpSolution<T> sol = lp::simplex_solve(model);
    if (!sol.ok())
        throw NumericError("m3pp22_interleave_fitc: no feasible set of off-diagonal rates (" +
                           std::string(lp::lp_status_name(sol.status)) + ")");

    M3pp22InterleaveResult<T> out;
    for (std::size_t i = 0; i < L; ++i) {
        T r1 = zero, r2 = zero;
        for (std::size_t j = i; j < L; ++j) r1 += sol.x[j];
        for (std::size_t j = 0; j <= i; ++j) r2 += sol.x[L + j];
        const T a = av(i, 0) + av(i, 1);
        const T d = r1 + r2;
        const T z = (binfv[i] - one) * pw(d, 3) * a;
        const T delta = num_sqrt(T(z / (two * r1 * r2)));
        const T l2 = a - r2 / d * delta;
        const T l1 = l2 + delta;

        Map<T> base;
        base.D0 = Matrix<T>(2, 2, zero);
        base.D1 = Matrix<T>(2, 2, zero);
        base.D0(0, 1) = r1;
        base.D0(1, 0) = r2;
        base.D1(0, 0) = l1;
        base.D1(1, 1) = l2;
        for (std::size_t h = 0; h < 2; ++h)
            base.D0(h, h) = -(base.D0(h, 0) + base.D0(h, 1) + base.D1(h, 0) + base.D1(h, 1));

        std::vector<T> ai;
        ai.push_back(av(i, 0));
        ai.push_back(av(i, 1));
        const M3pp22FitcCovResult<T> f = m3pp22_fitc_approx_cov_multiclass(base, ai, stv[i], t);
        out.parts.push_back(f.mmap);
        out.sigma.push_back(f.sigma);
    }

    out.mmap = m3pp2m_interleave(out.parts);
    return out;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_M3PP2M_INTERLEAVE_H
