/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MAMAP2M_FIT_H
#define LINE_API_MAM_MAMAP2M_FIT_H

/**
 * Fit a MAMAP(2,m): a second-order acyclic MAP marked with m classes, matching
 * the forward and backward per-class moments.
 *
 * Templated port of matlab/lib/m3a/m3a/mamap2m/mamap2m_fit_fb_multiclass.m,
 * mamap2m_fit_gamma_fb.m and mamap2m_fit_trace.m.
 *
 * As in `maph2m_fit.h`, the TIMING and the MARKING separate: an AMAP(2) fixed by
 * (M1, M2, M3, gamma) carries the inter-arrival law and its autocorrelation
 * decay, and the class marking splits the THREE arrival flows of the canonical
 * acyclic form among the m classes. The per-class forward and backward moments
 * are affine in that split,
 *
 *   q(j,c) = fF(c) q_f(j,c) + fB(c) q_b(j,c) + q_0(j,c),
 *
 * so the fit is again a convex quadratic program, here in 2k variables (a
 * forward and a backward moment per class) under three equality constraints, one
 * per flow, and 6k inequalities keeping every q in [0,1].
 *
 * SIX BRANCHES, and which one fires is decided by the AMAP's own degeneracies,
 * not by the data. With h1, h2 the phase means, r1 the branch probability out of
 * phase one and r2 the restart probability into phase two:
 *
 *  1. POISSON. The AMAP has collapsed to one state, or a denominator of the
 *     coefficients has vanished; only the class probabilities survive.
 *  2. DEGENERATE PHASE-TYPE (form 2, r2 = 0, r1 = 1). All three flows see the
 *     same class law p.
 *  3. CANONICAL PHASE-TYPE (form 1, r2 = 0). The MAP is really an APH(2), so the
 *     problem IS `maph2m_fit_multiclass` and is delegated to it.
 *  4. NON-CANONICAL PHASE-TYPE (r1 = 1). Only the forward moments are
 *     identifiable; the first flow is split uniformly.
 *  5. DEGENERATE MMAP (form 2, r2 = 0). Either the forward or the backward
 *     moments are fitted, whichever the weights prefer, and the third flow is
 *     split uniformly.
 *  6. GENERAL. The joint (F, B) program above.
 *
 * WHICH CANONICAL FORM. `form 1` (D1(1,2) = 0) carries a positive
 * autocorrelation decay and `form 2` (D1(1,1) = 0) a negative one; the
 * coefficient sets differ and are not interchangeable. Anything else is refused,
 * because the coefficients were derived for these two forms only.
 *
 * THE SOLVER IS NOT quadprog; see the same note in `maph2m_fit.h`. The
 * acceptance is the specification: the class probabilities are reproduced, the
 * inter-arrival law is the AMAP's, and the forward and backward moments approach
 * their targets as far as the feasibility of the split allows.
 *
 * ARITHMETIC: transcendental, through the fitters and the solver.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/mam/amap2_fit_gamma.h"
#include "line/api/mam/mamap22_fit_fs.h"
#include "line/api/mam/mamap_marked_poisson.h"
#include "line/api/mam/map_moment.h"
#include "line/api/mam/map_transform.h"
#include "line/api/mam/maph2m_fit.h"
#include "line/api/mam/mmap_compress.h"
#include "line/api/mam/mmap_lambda.h"
#include "line/api/mam/mmap_stats.h"
#include "line/api/trace/mtrace_backward_moment.h"
#include "line/api/trace/mtrace_forward_moment.h"
#include "line/api/trace/mtrace_sigma.h"
#include "line/api/trace/trace_gamma.h"
#include "line/api/trace/mtrace_pc.h"
#include "line/num/number.h"
#include "line/util/auglag.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/** The fitted MAMAP and the moments it achieved. */
template <class T>
struct Mamap2mFitResult {
    Mmap<T> mmap;
    std::vector<T> fF;  ///< achieved per-class forward moments
    std::vector<T> fB;  ///< achieved per-class backward moments
};

namespace mamapdetail {

/**
 * The shared quadratic program: minimize sum_v w(v) (x(v)/target(v) - 1)^2
 * subject to the equality rows summing each flow's split to one and the
 * inequality rows keeping every q in [0,1].
 */
template <class T, class QFun>
std::vector<T> solve_split(const std::vector<T>& target, const std::vector<T>& w,
                           std::size_t nflows, std::size_t k, QFun qof) {
    const T one = num_traits<T>::from_int(1), zero = num_traits<T>::from_int(0);
    const std::size_t nv = target.size();

    auto fobj = [&](const std::vector<T>& x) {
        T s = zero;
        for (std::size_t v = 0; v < nv; ++v) {
            const T r = T(x[v] / target[v] - one);
            s += w[v] * r * r;
        }
        return s;
    };
    auto heq = [&](const std::vector<T>& x) {
        std::vector<T> v(nflows, zero);
        for (std::size_t j = 0; j < nflows; ++j) {
            T s = zero;
            for (std::size_t c = 0; c < k; ++c) s += qof(x, j, c);
            v[j] = T(s - one);
        }
        return v;
    };
    auto gineq = [&](const std::vector<T>& x) {
        std::vector<T> v;
        v.reserve(2 * nflows * k);
        for (std::size_t c = 0; c < k; ++c)
            for (std::size_t j = 0; j < nflows; ++j) {
                const T qq = qof(x, j, c);
                v.push_back(T(qq - one));
                v.push_back(-qq);
            }
        return v;
    };

    std::vector<T> x0 = target;
    std::vector<Bound<T>> bounds(nv);
    for (std::size_t v = 0; v < nv; ++v) {
        bounds[v].lo = num_traits<T>::from_double(1e-6);
        bounds[v].hi = num_traits<T>::from_double(1e6);
        if (x0[v] < bounds[v].lo) x0[v] = bounds[v].lo;
        if (x0[v] > bounds[v].hi) x0[v] = bounds[v].hi;
    }
    return auglag(fobj, heq, gineq, x0, bounds).x;
}

/** The reference's clamp-and-renormalize on each flow's split. */
template <class T>
void fix_split(std::vector<std::vector<T>>& q, std::size_t k) {
    const T zero = num_traits<T>::from_int(0);
    for (std::size_t j = 0; j < q.size(); ++j) {
        T s = zero;
        for (std::size_t c = 0; c < k; ++c) {
            if (q[j][c] < zero) q[j][c] = zero;
            s += q[j][c];
        }
        if (!(num_traits<T>::to_double(s) > 0.0))
            throw NumericError("mamap2m_fit: a flow split lost all its mass");
        for (std::size_t c = 0; c < k; ++c) q[j][c] = T(q[j][c] / s);
    }
}

}  // namespace mamapdetail

/**
 * Mark a canonical acyclic AMAP(2) with m classes, matching the forward and
 * backward moments.
 *
 * @param map          the AMAP(2), in one of the two canonical acyclic forms
 * @param p            per-class probabilities
 * @param F            per-class target forward moments
 * @param B            per-class target backward moments
 * @param classWeights per-class weights; empty means uniform
 * @param fbWeights    the forward and backward weights; empty means (1, 1)
 */
template <class T>
Mamap2mFitResult<T> mamap2m_fit_fb_multiclass(const Map<T>& map, const std::vector<T>& p,
                                              const std::vector<T>& F, const std::vector<T>& B,
                                              const std::vector<T>& classWeights = std::vector<T>(),
                                              const std::vector<T>& fbWeights = std::vector<T>()) {
    static_assert(num_traits<T>::has_transcendental,
                  "mamap2m_fit_fb_multiclass solves a quadratic program");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    if (map.D0.rows() != 2)
        throw InputError("mamap2m_fit_fb_multiclass: the underlying MAP must be second order");
    if (!(num_abs(T(map.D0(1, 0))) <= zero))
        throw InputError("mamap2m_fit_fb_multiclass: the underlying MAP must be acyclic");

    int form;
    if (map.D1(0, 1) == zero)
        form = 1;
    else if (map.D1(0, 0) == zero)
        form = 2;
    else
        throw InputError(
            "mamap2m_fit_fb_multiclass: the underlying MAP must be in canonical acyclic form "
            "(D1(1,2) = 0 for a positive decay, D1(1,1) = 0 for a negative one)");

    const std::size_t k = p.size();
    if (k == 0) throw InputError("mamap2m_fit_fb_multiclass: no classes given");
    if (F.size() != k || B.size() != k)
        throw InputError(
            "mamap2m_fit_fb_multiclass: one forward and one backward moment per class is required");
    std::vector<T> cw = classWeights;
    if (cw.empty()) cw.assign(k, one);
    std::vector<T> fbw = fbWeights;
    if (fbw.empty()) fbw.assign(2, one);

    const T h1 = T(-one / map.D0(0, 0));
    const T h2 = T(-one / map.D0(1, 1));
    const T r1 = T(map.D0(0, 1) * h1);
    const T r2 = T(map.D1(1, 1) * h2);
    const double dt = 1e-8;
    const double dr1 = num_traits<T>::to_double(r1), dr2 = num_traits<T>::to_double(r2);
    const double dh1 = num_traits<T>::to_double(h1), dh2 = num_traits<T>::to_double(h2);

    Mamap2mFitResult<T> out;
    out.mmap.D0 = map.D0;
    out.mmap.D1 = map.D1;
    out.mmap.Dc.assign(k, Matrix<T>(2, 2, zero));

    auto finish = [&](std::vector<std::vector<T>>& q) {
        mamapdetail::fix_split(q, k);
        for (std::size_t c = 0; c < k; ++c) {
            if (form == 1) {
                out.mmap.Dc[c](0, 0) = T(out.mmap.D1(0, 0) * q[0][c]);
                out.mmap.Dc[c](1, 0) = T(out.mmap.D1(1, 0) * q[1][c]);
                out.mmap.Dc[c](1, 1) = T(out.mmap.D1(1, 1) * q[2][c]);
            } else {
                out.mmap.Dc[c](0, 1) = T(out.mmap.D1(0, 1) * q[0][c]);
                out.mmap.Dc[c](1, 0) = T(out.mmap.D1(1, 0) * q[1][c]);
                out.mmap.Dc[c](1, 1) = T(out.mmap.D1(1, 1) * q[2][c]);
            }
        }
        const std::vector<unsigned> one_order(1, 1u);
        const Matrix<T> fm = mmap_forward_moment(out.mmap, one_order, true);
        const std::vector<std::vector<T>> bm = mmap_backward_moment(out.mmap, one_order, true);
        out.fF.assign(k, zero);
        out.fB.assign(k, zero);
        for (std::size_t c = 0; c < k; ++c) {
            out.fF[c] = fm(c, 0);
            out.fB[c] = bm[c][0];
        }
    };

    // 1. POISSON: the coefficients have a vanishing denominator, so nothing but
    // the class probabilities is identifiable.
    const bool poisson1 =
        form == 1 && (dr1 < dt || dr2 > 1.0 - dt || std::fabs(dh2 - dh1 * dr2) < dt ||
                      std::fabs(dh1 - dh2 + dh2 * dr1) < dt);
    const bool poisson2 =
        form == 2 && (dr2 > 1.0 - dt || std::fabs(dh1 - dh2 + dh2 * dr1) < dt ||
                      std::fabs(dh1 - dh2 - dh1 * dr1 + dh1 * dr1 * dr2) < dt);
    if (poisson1 || poisson2) {
        out.mmap = mamapdetail::marked_poisson(map_mean(map), p);
        const std::vector<unsigned> one_order(1, 1u);
        const Matrix<T> fm = mmap_forward_moment(out.mmap, one_order, true);
        const std::vector<std::vector<T>> bm = mmap_backward_moment(out.mmap, one_order, true);
        out.fF.assign(k, zero);
        out.fB.assign(k, zero);
        for (std::size_t c = 0; c < k; ++c) {
            out.fF[c] = fm(c, 0);
            out.fB[c] = bm[c][0];
        }
        return out;
    }

    // 2. DEGENERATE PHASE-TYPE: all three flows carry the same class law.
    if (form == 2 && dr2 < dt && std::fabs(1.0 - dr1) < dt) {
        std::vector<std::vector<T>> q(3, std::vector<T>(k, zero));
        for (std::size_t c = 0; c < k; ++c) q[0][c] = q[1][c] = q[2][c] = p[c];
        finish(q);
        return out;
    }

    // 3. CANONICAL PHASE-TYPE: the MAP is an APH(2); delegate.
    if (form == 1 && dr2 < dt) {
        Map<T> aph = map;
        aph.D1(1, 1) = zero;
        aph = map_normalize(aph);
        const Maph2mFitResult<T> r = maph2m_fit_multiclass(aph, p, B, cw);
        out.mmap = r.maph;
        const std::vector<unsigned> one_order(1, 1u);
        const Matrix<T> fm = mmap_forward_moment(out.mmap, one_order, true);
        const std::vector<std::vector<T>> bm = mmap_backward_moment(out.mmap, one_order, true);
        out.fF.assign(k, zero);
        out.fB.assign(k, zero);
        for (std::size_t c = 0; c < k; ++c) {
            out.fF[c] = fm(c, 0);
            out.fB[c] = bm[c][0];
        }
        return out;
    }

    // 4. NON-CANONICAL PHASE-TYPE: only the forward moments are identifiable.
    if (std::fabs(1.0 - dr1) < dt) {
        std::vector<std::vector<T>> qf(2, std::vector<T>(k, zero)),
            q0(2, std::vector<T>(k, zero));
        for (std::size_t c = 0; c < k; ++c) {
            qf[0][c] = T(p[c] * (-one / ((h1 + h2 * (r1 - one)) * (r2 - one) *
                                         (r1 + r2 - r1 * r2))));
            q0[0][c] = T(p[c] * (h2 / ((r2 - one) * (r1 + r2 - r1 * r2) *
                                       (h1 - h2 + h2 * r1))));
            qf[1][c] = T(p[c] * (-one / (r2 * (h1 + h2 * (r1 - one)) * (r1 + r2 - r1 * r2))));
            q0[1][c] = T(p[c] * ((h1 + h2 * r1) /
                                 (r2 * (r1 + r2 - r1 * r2) * (h1 - h2 + h2 * r1))));
        }
        std::vector<T> w(k, zero);
        for (std::size_t c = 0; c < k; ++c) w[c] = T(cw[c] * fbw[0]);
        auto qof = [&](const std::vector<T>& x, std::size_t j, std::size_t c) {
            return T(x[c] * qf[j][c] + q0[j][c]);
        };
        const std::vector<T> x = mamapdetail::solve_split(F, w, 2, k, qof);
        std::vector<std::vector<T>> q(3, std::vector<T>(k, zero));
        const T uni = T(one / num_traits<T>::from_int(static_cast<long>(k)));
        for (std::size_t c = 0; c < k; ++c) {
            q[0][c] = uni;
            q[1][c] = T(x[c] * qf[0][c] + q0[0][c]);
            q[2][c] = T(x[c] * qf[1][c] + q0[1][c]);
        }
        finish(q);
        return out;
    }

    // 5. DEGENERATE MMAP (gamma < 0): forward or backward, per the weights.
    if (form == 2 && dr2 < dt) {
        const bool doForward = !(fbw[0] < fbw[1]);
        std::vector<std::vector<T>> qc(2, std::vector<T>(k, zero)),
            q0(2, std::vector<T>(k, zero));
        for (std::size_t c = 0; c < k; ++c) {
            if (doForward) {
                qc[0][c] = T(p[c] * (-(r1 - two) / ((h1 + h2 * (r1 - one)) * (r1 - one))));
                q0[0][c] = T(p[c] * (one - (h1 + h2) / ((r1 - one) * (h1 - h2 + h2 * r1))));
                qc[1][c] = T(p[c] * (-(r1 - two) / (h1 + h2 * (r1 - one))));
                q0[1][c] = T(p[c] * ((h2 * (r1 - two)) / (h1 - h2 + h2 * r1)));
            } else {
                qc[0][c] = T(p[c] * (-(r1 - two) / ((h2 + h1 * (r1 - one)) * (r1 - one))));
                q0[0][c] = T(p[c] * (one - (h1 + h2) / ((r1 - one) * (h2 - h1 + h1 * r1))));
                qc[1][c] = T(p[c] * (-(r1 - two) / (h2 + h1 * (r1 - one))));
                q0[1][c] = T(p[c] * ((h1 * (r1 - two)) / (h2 - h1 + h1 * r1)));
            }
        }
        std::vector<T> w(k, zero);
        for (std::size_t c = 0; c < k; ++c) w[c] = T(cw[c] * (doForward ? fbw[0] : fbw[1]));
        auto qof = [&](const std::vector<T>& x, std::size_t j, std::size_t c) {
            return T(x[c] * qc[j][c] + q0[j][c]);
        };
        const std::vector<T> x = mamapdetail::solve_split(doForward ? F : B, w, 2, k, qof);
        std::vector<std::vector<T>> q(3, std::vector<T>(k, zero));
        const T uni = T(one / num_traits<T>::from_int(static_cast<long>(k)));
        for (std::size_t c = 0; c < k; ++c) {
            q[0][c] = T(x[c] * qc[0][c] + q0[0][c]);
            q[1][c] = T(x[c] * qc[1][c] + q0[1][c]);
            q[2][c] = uni;
        }
        finish(q);
        return out;
    }

    // 6. GENERAL: the joint (F, B) program in 2k variables.
    std::vector<std::vector<T>> qf(3, std::vector<T>(k, zero)), qb(3, std::vector<T>(k, zero)),
        q0(3, std::vector<T>(k, zero));
    for (std::size_t c = 0; c < k; ++c) {
        if (form == 1) {
            const T z = T(r1 * r2 - r2 + one);
            qf[0][c] = zero;
            qb[0][c] = T(-(p[c] * z) / ((h2 - h1 * r2) * (r1 - one) * (r2 - one)));
            q0[0][c] = T((p[c] * (h1 + h2 - h1 * r2) * z) /
                         ((h2 - h1 * r2) * (r1 - one) * (r2 - one)));
            qf[1][c] = T(-(p[c] * z) / (r1 * (h1 + h2 * (r1 - one)) * (r2 - one)));
            qb[1][c] = T(-(p[c] * z) / (r1 * (h2 - h1 * r2) * (r2 - one)));
            q0[1][c] = T((p[c] * z) / ((r1 - one) * (r2 - one)) +
                         (h1 * p[c] * z) / (r1 * (h2 - h1 * r2) * (r2 - one)) -
                         (h1 * p[c] * z) / (r1 * (h1 + h2 * (r1 - one)) * (r1 - one) * (r2 - one)));
            qf[2][c] = T(-(p[c] * z) / (r1 * r2 * (h1 - h2 + h2 * r1)));
            qb[2][c] = zero;
            q0[2][c] = T((p[c] * (h1 + h2 * r1) * z) / (r1 * r2 * (h1 - h2 + h2 * r1)));
        } else {
            const T z = T(r1 + r2 - r1 * r2 - two);
            const T d2 = T(h1 - h2 - h1 * r1 + h1 * r1 * r2);
            qf[0][c] = zero;
            qb[0][c] = T(-(p[c] * z) / ((r1 - one) * (r2 - one) * d2));
            q0[0][c] = T((p[c] * (h2 + h1 * r1 - h1 * r1 * r2) * z) /
                         ((r1 - one) * (r2 - one) * d2));
            qf[1][c] = T((p[c] * z) / ((r2 - one) * (h1 - h2 + h2 * r1)));
            qb[1][c] = zero;
            q0[1][c] = T(-(h2 * p[c] * z) / ((r2 - one) * (h1 - h2 + h2 * r1)));
            qf[2][c] = T((p[c] * z) / (r2 * (h1 + h2 * (r1 - one))));
            qb[2][c] = T((p[c] * z) / (r2 * d2));
            q0[2][c] = T((h1 * p[c] * z) / (r2 * (h1 + h2 * (r1 - one)) * (r1 - one)) -
                         (h1 * p[c] * z) / (r2 * d2) - (p[c] * z) / (r2 * (r1 - one)));
        }
    }

    // The variables interleave: x[2c] is F(c), x[2c+1] is B(c).
    std::vector<T> target(2 * k, zero), w(2 * k, zero);
    for (std::size_t c = 0; c < k; ++c) {
        target[2 * c] = F[c];
        target[2 * c + 1] = B[c];
        w[2 * c] = T(cw[c] * fbw[0]);
        w[2 * c + 1] = T(cw[c] * fbw[1]);
    }
    auto qof = [&](const std::vector<T>& x, std::size_t j, std::size_t c) {
        return T(x[2 * c] * qf[j][c] + x[2 * c + 1] * qb[j][c] + q0[j][c]);
    };
    const std::vector<T> x = mamapdetail::solve_split(target, w, 3, k, qof);
    std::vector<std::vector<T>> q(3, std::vector<T>(k, zero));
    for (std::size_t c = 0; c < k; ++c)
        for (std::size_t j = 0; j < 3; ++j) q[j][c] = qof(x, j, c);
    finish(q);
    return out;
}

/**
 * Fit a MAMAP(2,m) to three moments, the decay rate, the class probabilities and
 * the forward and backward moments, over every AMAP(2) form.
 */
template <class T>
Mmap<T> mamap2m_fit_gamma_fb(const T& M1, const T& M2, const T& M3, const T& GAMMA,
                             const std::vector<T>& p, const std::vector<T>& F,
                             const std::vector<T>& B) {
    const Amap2FitGammaResult<T> a = amap2_fit_gamma(M1, M2, M3, GAMMA);
    if (a.amaps.empty() || (a.amaps.size() == 1 && a.amaps[0].order() == 1))
        return mamapdetail::marked_poisson(M1, p);

    Mmap<T> best;
    double bestErr = 0.0;
    bool have = false;
    for (std::size_t j = 0; j < a.amaps.size(); ++j) {
        Mamap2mFitResult<T> r;
        try {
            r = mamap2m_fit_fb_multiclass(a.amaps[j], p, F, B);
        } catch (const Error&) {
            continue;
        }
        double err = 0.0;
        for (std::size_t c = 0; c < p.size(); ++c) {
            const double df = num_traits<T>::to_double(T(r.fF[c] / F[c])) - 1.0;
            const double db = num_traits<T>::to_double(T(r.fB[c] / B[c])) - 1.0;
            err += df * df + db * db;
        }
        if (!have || err < bestErr) {
            bestErr = err;
            best = r.mmap;
            have = true;
        }
    }
    if (!have)
        throw NumericError(
            "mamap2m_fit_gamma_fb: no AMAP(2) form admits a valid class split for the requested "
            "class probabilities and forward/backward moments");
    return best;
}

/**
 * The full `mamap2m_fit` dispatcher.
 *
 * Port of matlab/lib/m3a/m3a/mamap2m/mamap2m_fit.m. It chooses WHICH pair of
 * descriptors to match from the AMAP's degeneracies and from the caller's
 * weights over (forward, backward, sigma):
 *
 *  - more than two classes: F+B, since the sigma fitters are two-class only;
 *  - a negligible gamma: the process is renewal, so a MAPH is fitted instead;
 *  - otherwise one AMAP(2) form at a time, taking F+B, F+S or B+S as the
 *    degeneracy allows and the weights prefer, and keeping the form whose
 *    achieved descriptors land closest.
 *
 * The sigma branches call `mamap22_fit_fs_multiclass` and
 * `mamap22_fit_bs_multiclass` (mamap22_fit_fs.h, mamap22_fit_bs.h); they fire
 * when the weights prefer sigma over forward or backward, and on the
 * degenerate AMAP shapes. With the reference's DEFAULT weights (1,1,1) F+B is
 * preferred.
 *
 * @param fbsWeights the (forward, backward, sigma) weights; empty means (1,1,1)
 */
template <class T>
Mmap<T> mamap2m_fit(const T& M1, const T& M2, const T& M3, const T& GAMMA,
                    const std::vector<T>& p, const std::vector<T>& F, const std::vector<T>& B,
                    const Matrix<T>& S, const std::vector<T>& fbsWeights = std::vector<T>()) {
    const T one = num_traits<T>::from_int(1);
    std::vector<T> w = fbsWeights;
    if (w.empty()) w.assign(3, one);
    if (w.size() != 3) throw InputError("mamap2m_fit: three (F, B, S) weights are required");
    const Matrix<T>& S_or_empty = S;

    const double gammatol = 1e-4, degentol = 1e-8;
    if (p.size() > 2) return mamap2m_fit_gamma_fb(M1, M2, M3, GAMMA, p, F, B);
    if (std::fabs(num_traits<T>::to_double(GAMMA)) < gammatol) return maph2m_fit(M1, M2, M3, p, B);

    const Amap2FitGammaResult<T> a = amap2_fit_gamma(M1, M2, M3, GAMMA);
    if (a.amaps.empty() || (a.amaps.size() == 1 && a.amaps[0].order() == 1))
        return mamapdetail::marked_poisson(M1, p);

    const bool preferFB = !(w[0] < w[2]) && !(w[1] < w[2]);
    const bool preferFS = !preferFB && !(w[0] < w[1]);

    Mmap<T> best;
    double bestErr = 0.0;
    bool have = false;
    for (std::size_t j = 0; j < a.amaps.size(); ++j) {
        const Map<T>& mp = a.amaps[j];
        if (mp.order() != 2) continue;
        const T h1 = T(-one / mp.D0(0, 0)), h2 = T(-one / mp.D0(1, 1));
        const T r1 = T(h1 * mp.D0(0, 1)), r2 = T(h2 * mp.D1(1, 1));
        const double dh1 = num_traits<T>::to_double(h1), dh2 = num_traits<T>::to_double(h2);
        const double dr1 = num_traits<T>::to_double(r1), dr2 = num_traits<T>::to_double(r2);
        const bool posGamma = num_traits<T>::to_double(GAMMA) > 0.0;

        // The reference's degeneracy ladder: several shapes identify only one
        // of the two moments, and there sigma is the second descriptor.
        bool needsFs = false, needsBs = false;
        if (posGamma) {
            if (std::fabs(dh2 - dh1 * dr2) < degentol) needsFs = true;
            else if (std::fabs(dh1 - dh2 + dh2 * dr1) < degentol) needsBs = true;
            else if (1.0 - dr1 < degentol) needsFs = true;
        } else {
            if (std::fabs(dh1 - dh2 - dh1 * dr1 + dh1 * dr1 * dr2) < degentol) needsFs = true;
            else if (std::fabs(dh1 - dh2 + dh2 * dr1) < degentol) needsBs = true;
        }
        // Outside the degeneracies the weights decide, as in the reference.
        if (!needsFs && !needsBs && !preferFB) {
            if (preferFS) needsFs = true;
            else needsBs = true;
        }

        Mmap<T> cand;
        std::vector<T> cF, cB;
        try {
            if (needsFs) {
                const Mamap22FsFitResult<T> rf = mamap22_fit_fs_multiclass(mp, p, F, S_or_empty);
                cand = rf.mmap;
            } else if (needsBs) {
                const Mamap22FitResult<T> rb = mamap22_fit_bs_multiclass(mp, p, B, S_or_empty);
                cand = rb.mmap;
            } else {
                const Mamap2mFitResult<T> r = mamap2m_fit_fb_multiclass(mp, p, F, B);
                cand = r.mmap;
            }
        } catch (const Error&) {
            continue;
        }
        // Score every candidate the same way, on the descriptors it achieved.
        const std::vector<unsigned> ord1(1, 1u);
        const Matrix<T> fmc = mmap_forward_moment(cand, ord1, true);
        const std::vector<std::vector<T>> bmc = mmap_backward_moment(cand, ord1, true);
        const Matrix<T> fsc = mmap_sigma(cand);
        double err = 0.0;
        for (std::size_t c = 0; c < p.size(); ++c) {
            const double df = num_traits<T>::to_double(T(F[c] / fmc(c, 0))) - 1.0;
            const double db = num_traits<T>::to_double(T(B[c] / bmc[c][0])) - 1.0;
            err += num_traits<T>::to_double(w[0]) * df * df +
                   num_traits<T>::to_double(w[1]) * db * db;
        }
        if (S_or_empty.rows() > 0 && fsc.rows() > 0) {
            const double ds = num_traits<T>::to_double(T(S_or_empty(0, 0) / fsc(0, 0))) - 1.0;
            err += num_traits<T>::to_double(w[2]) * ds * ds;
        }
        if (!have || err < bestErr) {
            bestErr = err;
            best = cand;
            have = true;
        }
    }
    if (!have)
        throw NumericError(
            "mamap2m_fit: no AMAP(2) form admits a valid class split for the requested "
            "descriptors");
    return best;
}

/**
 * Fit a MAMAP(2,m) from a marked trace through the (F, B) pair alone.
 *
 * Port of matlab/lib/m3a/m3a/mamap2m/mamap2m_fit_gamma_fb_trace.m: the
 * descriptors are the class probabilities and the per-class FORWARD and
 * BACKWARD first moments, with no sigma and no descriptor-pair selection.
 * `mamap2m_fit_trace` is the full dispatcher of the reference.
 *
 * @param Tv the inter-arrival times
 * @param A  the class of each arrival, 1-based
 */
template <class T>
Mmap<T> mamap2m_fit_gamma_fb_trace(const std::vector<T>& Tv, const std::vector<int>& A) {
    if (Tv.empty() || Tv.size() != A.size())
        throw InputError(
            "mamap2m_fit_gamma_fb_trace: the trace and its labels must agree in length");
    const T zero = num_traits<T>::from_int(0);
    T m1 = zero, m2 = zero, m3 = zero;
    for (std::size_t i = 0; i < Tv.size(); ++i) {
        const T x = Tv[i];
        m1 += x;
        m2 += x * x;
        m3 += x * x * x;
    }
    const T n = num_traits<T>::from_int(static_cast<long>(Tv.size()));
    const std::vector<unsigned> one_order(1, 1u);
    const std::vector<T> p = trace::mtrace_pc<T>(A);
    const Matrix<T> fm = trace::mtrace_forward_moment(Tv, A, one_order);
    const Matrix<T> bm = trace::mtrace_backward_moment(Tv, A, one_order);
    std::vector<T> F(p.size(), zero), B(p.size(), zero);
    for (std::size_t c = 0; c < p.size(); ++c) {
        F[c] = fm(c, 0);
        B[c] = bm(c, 0);
    }
    return mamap2m_fit_gamma_fb(T(m1 / n), T(m2 / n), T(m3 / n),
                                line::trace::trace_gamma(Tv).gamma, p, F, B);
}

/**
 * Fit a MAPH(2,m) or MAMAP(2,m) matching the characteristics of a marked
 * trace.
 *
 * Port of matlab/lib/m3a/m3a/mamap2m/mamap2m_fit_trace.m: the class
 * probabilities are always matched exactly; the remaining two characteristics
 * default to the forward and backward moments unless the underlying AMAP(2)
 * is degenerate, and `fbsWeights` moves the preference among (forward,
 * backward, sigma). Unlike `mamap2m_fit_gamma_fb_trace` this computes the
 * class transition probabilities (`mtrace_sigma`) and dispatches through the
 * full `mamap2m_fit` descriptor selection.
 *
 * @param Tv the inter-arrival times
 * @param A  the class of each arrival, 1-based
 * @param fbsWeights the (forward, backward, sigma) weights; empty means (1,1,1)
 */
template <class T>
Mmap<T> mamap2m_fit_trace(const std::vector<T>& Tv, const std::vector<int>& A,
                          const std::vector<T>& fbsWeights = std::vector<T>()) {
    if (Tv.empty() || Tv.size() != A.size())
        throw InputError("mamap2m_fit_trace: the trace and its labels must agree in length");
    const T zero = num_traits<T>::from_int(0);
    T m1 = zero, m2 = zero, m3 = zero;
    for (std::size_t i = 0; i < Tv.size(); ++i) {
        const T x = Tv[i];
        m1 += x;
        m2 += x * x;
        m3 += x * x * x;
    }
    const T n = num_traits<T>::from_int(static_cast<long>(Tv.size()));
    const std::vector<unsigned> one_order(1, 1u);
    const std::vector<T> p = trace::mtrace_pc<T>(A);
    const Matrix<T> fm = trace::mtrace_forward_moment(Tv, A, one_order);
    const Matrix<T> bm = trace::mtrace_backward_moment(Tv, A, one_order);
    std::vector<T> F(p.size(), zero), B(p.size(), zero);
    for (std::size_t c = 0; c < p.size(); ++c) {
        F[c] = fm(c, 0);
        B[c] = bm(c, 0);
    }
    const Matrix<T> S = trace::mtrace_sigma<T>(A);
    return mamap2m_fit(T(m1 / n), T(m2 / n), T(m3 / n), line::trace::trace_gamma(Tv).gamma, p,
                       F, B, S, fbsWeights);
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MAMAP2M_FIT_H
