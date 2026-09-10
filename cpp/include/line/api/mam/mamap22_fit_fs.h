/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MAMAP22_FIT_FS_H
#define LINE_API_MAM_MAMAP22_FIT_FS_H

/**
 * Fit a MAMAP(2,2) matching the FORWARD moment and the class TRANSITION
 * probability sigma.
 *
 * Templated port of matlab/lib/m3a/m3a/mamap22/mamap22_fit_fs_multiclass.m, the
 * twin of `mamap22_fit_bs.h`. The structure is the same -- a closed-form inverse
 * over the coefficient tables of `mamap2m_coefficients.h`, with five degenerate
 * branches ahead of it -- and only the inverse itself differs, reading the
 * FORWARD block of the tables (U(1..6), G(13..15), Y(2)) where the backward
 * fitter reads U(7..12), G(10..12) and Y(3).
 *
 * TWO CLASSES ONLY, as in the backward twin, and for the same reason.
 *
 * ONE BRANCH DIFFERS SUBSTANTIVELY from the backward fitter and is worth
 * naming: the CANONICAL PHASE-TYPE case (form 1, r2 = 0). The reference falls
 * back to `maph2m_fit_multiclass`, which fits BACKWARD moments, and it has no
 * forward targets to give it -- so it sets both backward targets to the ordinary
 * mean and WARNS that the caller should have used B+S instead. That warning is
 * reproduced as a named diagnostic on the result rather than dropped, because a
 * caller who asked for a forward fit and silently got a mean-matched backward
 * one has no other way to find out. It also re-fits the timing as an APH(2) when
 * the SCV has fallen to one, which is the hypoexponential boundary where the
 * canonical form stops being informative.
 *
 * WHAT IS NOT PORTED, refused by name: the same nonconvex repair as the backward
 * twin -- YALMIP `bmibnb`, a spatial branch-and-bound returning a GLOBAL optimum
 * of a bilinear program. `line/util/auglag.h` finds a local KKT point, which on a
 * nonconvex set is a different answer; substituting it would report a fit the
 * reference would not have chosen. The gamma < 0 degeneracy's sigma arm IS
 * ported: the reference states it as a YALMIP program, but its feasible set is
 * an interval and its objective is (x - S11)^2, so the projection is its exact
 * global optimum.
 *
 * ARITHMETIC: transcendental.
 */

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "line/api/mam/m3a_fit_from.h"
#include "line/api/mam/mamap22_fit_bs.h"
#include "line/api/mam/amap2_fit_gamma.h"
#include "line/api/mam/amap2_fitall_gamma.h"
#include "line/api/mam/mamap2m_coefficients.h"
#include "line/api/mam/mamap_marked_poisson.h"
#include "line/api/trace/mtrace_forward_moment.h"
#include "line/api/trace/mtrace_pc.h"
#include "line/api/trace/mtrace_sigma.h"
#include "line/api/trace/trace_gamma.h"
#include "line/api/mam/map_moment.h"
#include "line/api/mam/map_transform.h"
#include "line/api/mam/maph2m_fit.h"
#include "line/api/mam/mmap_compress.h"
#include "line/api/mam/mmap_lambda.h"
#include "line/api/mam/mmap_stats.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/** The forward-plus-sigma result; `warning` carries the reference's diagnostic. */
template <class T>
struct Mamap22FsFitResult {
    Mmap<T> mmap;
    std::vector<T> fF;  ///< achieved per-class forward moments
    Matrix<T> fS;       ///< achieved class transition probabilities
    bool exact = false;
    std::string warning;  ///< empty unless a branch substituted its targets
};

namespace fsdetail {

/** The closed-form F+S inverse of the first canonical form. */
template <class T>
void fit_can1(const Mamap2mCoefficients<T>& c, const T& p1, const T& vF1, const T& vS11,
              double denumtol, T* q1, T* q2, T* q3) {
    const T den = T(p1 * (c.U[4] * vF1 + c.U[5]));
    if (std::fabs(num_traits<T>::to_double(den)) < denumtol) {
        *q1 = *q2 = *q3 = p1;
        return;
    }
    *q2 = T((c.U[0] * vF1 * vF1 * p1 * p1 + c.U[1] * vF1 * p1 * p1 + c.U[2] * vS11 +
             c.U[3] * p1 * p1) /
            den);
    *q1 = T(-(c.G[14] * p1 - vF1 * c.G[2] * p1 +
              (c.G[2] * c.G[13] - c.G[1] * c.G[14]) * (*q2)) /
            c.Y[1]);
    *q3 = T((c.G[12] * p1 - vF1 * c.G[0] * p1 +
             (c.G[0] * c.G[13] - c.G[1] * c.G[12]) * (*q2)) /
            c.Y[1]);
}

/** The closed-form F+S inverse of the second canonical form. */
template <class T>
void fit_can2(const Mamap2mCoefficients<T>& c, const T& p1, const T& vF1, const T& vS11,
              double denumtol, T* q1, T* q2, T* q3) {
    const T den = T(c.U[4] * vF1 * p1 + c.U[5] * p1);
    if (std::fabs(num_traits<T>::to_double(den)) < denumtol) {
        *q1 = *q2 = *q3 = p1;
        return;
    }
    *q3 = T((c.U[0] * vF1 * vF1 * p1 * p1 + c.U[1] * vF1 * p1 * p1 + c.U[2] * p1 * p1 +
             c.U[3] * vS11) /
            den);
    *q1 = T(-(c.G[12] * p1 - vF1 * c.G[1] * p1 +
              (c.G[1] * c.G[13] - c.G[2] * c.G[12]) * (*q3)) /
            c.Y[1]);
    *q2 = T((c.G[11] * p1 - vF1 * c.G[0] * p1 +
             (c.G[0] * c.G[13] - c.G[2] * c.G[11]) * (*q3)) /
            c.Y[1]);
}

}  // namespace fsdetail

/**
 * @param map          the AMAP(2), in one of the two canonical acyclic forms
 * @param p            the two class probabilities
 * @param F            the two target forward moments
 * @param S            the target class transition matrix; only S(0,0) is used
 * @param classWeights per-class weights; empty means uniform
 * @param fsWeights    the (forward, sigma) weights; empty means (1, 1)
 * @param adjust       repair an infeasible closed form; the repair is unported
 */
template <class T>
Mamap22FsFitResult<T> mamap22_fit_fs_multiclass(
    const Map<T>& map, const std::vector<T>& p, const std::vector<T>& F, const Matrix<T>& S,
    const std::vector<T>& classWeights = std::vector<T>(),
    const std::vector<T>& fsWeights = std::vector<T>(), bool adjust = true) {
    static_assert(num_traits<T>::has_transcendental,
                  "mamap22_fit_fs_multiclass inverts a moment system");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    if (map.D0.rows() != 2)
        throw InputError("mamap22_fit_fs_multiclass: the underlying MAP must be second order");
    if (!(num_abs(T(map.D0(1, 0))) <= zero))
        throw InputError("mamap22_fit_fs_multiclass: the underlying MAP must be acyclic");
    int form;
    if (map.D1(0, 1) == zero) form = 1;
    else if (map.D1(0, 0) == zero) form = 2;
    else
        throw InputError(
            "mamap22_fit_fs_multiclass: the underlying MAP must be in canonical acyclic form");
    if (p.size() != 2)
        throw InputError(
            "mamap22_fit_fs_multiclass: fitting the forward moment and the transition "
            "probabilities supports two classes only");
    if (F.size() != 2)
        throw InputError("mamap22_fit_fs_multiclass: one forward moment per class is required");
    if (S.rows() < 1 || S.cols() < 1)
        throw InputError("mamap22_fit_fs_multiclass: the transition matrix is empty");

    std::vector<T> cw = classWeights;
    if (cw.empty()) cw.assign(2, one);
    std::vector<T> fw = fsWeights;
    if (fw.empty()) fw.assign(2, one);

    // feastol is 1e-3 here and 1e-4 in the backward twin; that is the
    // reference's own asymmetry, not a transcription slip.
    const double degentol = 1e-6, feastol = 1e-3, denumtol = 1e-12;
    Map<T> mp = map;
    T h1 = T(-one / mp.D0(0, 0)), h2 = T(-one / mp.D0(1, 1));
    T r1 = T(mp.D0(0, 1) * h1), r2 = T(mp.D1(1, 1) * h2);
    auto dv = [](const T& v) { return num_traits<T>::to_double(v); };

    Mamap22FsFitResult<T> out;
    out.mmap.D0 = mp.D0;
    out.mmap.D1 = mp.D1;
    out.mmap.Dc.assign(2, Matrix<T>(2, 2, zero));

    auto finish = [&](T q1, T q2, T q3) {
        auto fix = [&](const T& q) {
            T v = q;
            if (v < zero) v = zero;
            if (v > one) v = one;
            return v;
        };
        q1 = fix(q1);
        q2 = fix(q2);
        q3 = fix(q3);
        if (form == 1) {
            out.mmap.Dc[0](0, 0) = T(out.mmap.D1(0, 0) * q1);
            out.mmap.Dc[1](0, 0) = T(out.mmap.D1(0, 0) * (one - q1));
        } else {
            out.mmap.Dc[0](0, 1) = T(out.mmap.D1(0, 1) * q1);
            out.mmap.Dc[1](0, 1) = T(out.mmap.D1(0, 1) * (one - q1));
        }
        out.mmap.Dc[0](1, 0) = T(out.mmap.D1(1, 0) * q2);
        out.mmap.Dc[1](1, 0) = T(out.mmap.D1(1, 0) * (one - q2));
        out.mmap.Dc[0](1, 1) = T(out.mmap.D1(1, 1) * q3);
        out.mmap.Dc[1](1, 1) = T(out.mmap.D1(1, 1) * (one - q3));
        const std::vector<unsigned> ord(1, 1u);
        const Matrix<T> fm = mmap_forward_moment(out.mmap, ord, true);
        out.fF.assign(2, zero);
        for (std::size_t c = 0; c < 2; ++c) out.fF[c] = fm(c, 0);
        out.fS = mmap_sigma(out.mmap);
    };
    auto feasible = [&](const T& q) { return dv(q) >= -feastol && dv(q) <= 1.0 + feastol; };
    auto adopt = [&](const Mmap<T>& m) {
        out.mmap = m;
        const std::vector<unsigned> ord(1, 1u);
        const Matrix<T> fm = mmap_forward_moment(out.mmap, ord, true);
        out.fF.assign(2, zero);
        for (std::size_t c = 0; c < 2; ++c) out.fF[c] = fm(c, 0);
        out.fS = mmap_sigma(out.mmap);
    };

    // ---- the Poisson perturbation ---------------------------------------
    // The forward fitter's degeneracy is h1 - h2 + h2*r1 in BOTH forms, and it
    // perturbs h1, not h2. The backward twin's is a different expression and
    // perturbs a different parameter; do not share this block with it.
    const bool degen1 = form == 1 && (dv(r1) < degentol || dv(r2) > 1.0 - degentol ||
                                      std::fabs(dv(h1) - dv(h2) + dv(h2) * dv(r1)) < degentol);
    const bool degen2 =
        form == 2 && (dv(r2) > 1.0 - degentol ||
                      std::fabs(dv(h1) - dv(h2) + dv(h2) * dv(r1)) < degentol);
    if (degen1 || degen2) {
        if (dv(r1) < degentol) r1 = num_traits<T>::from_double(degentol);
        if (dv(r2) > 1.0 - degentol) r2 = num_traits<T>::from_double(1.0 - degentol);
        if (std::fabs(dv(h1) - dv(h2) + dv(h2) * dv(r1)) < degentol)
            h1 = T(h2 * (one - r1) + num_traits<T>::from_double(degentol));
        mp.D0 = Matrix<T>(2, 2, zero);
        mp.D1 = Matrix<T>(2, 2, zero);
        mp.D0(0, 0) = T(-one / h1);
        mp.D0(0, 1) = T(r1 / h1);
        mp.D0(1, 1) = T(-one / h2);
        if (form == 1) {
            mp.D1(0, 0) = T((one - r1) / h1);
            mp.D1(1, 0) = T(r2 / h2);
            mp.D1(1, 1) = T((one - r2) / h2);
        } else {
            mp.D1(0, 1) = T((one - r1) / h1);
            mp.D1(1, 0) = T(r2 / h2);
            mp.D1(1, 1) = T((one - r2) / h2);
        }
        mp = map_normalize(mp);
        out.mmap.D0 = mp.D0;
        out.mmap.D1 = mp.D1;
    }

    // ---- the degenerate ladder ------------------------------------------
    if (form == 2 && dv(r2) < degentol && std::fabs(1.0 - dv(r1)) < degentol) {
        finish(p[0], p[0], p[0]);
        return out;
    }
    if (form == 1 && dv(r2) < degentol) {
        // The reference has no forward targets for the MAPH fitter, so it uses
        // the ordinary mean as both backward targets and WARNS; the warning is
        // carried out rather than dropped.
        Map<T> aph = mp;
        aph.D1(1, 1) = zero;
        aph = map_normalize(aph);
        if (dv(map_scv(aph)) < 1.0 + degentol) aph = aph2_fit_map(mp).aph;
        const T mean = map_mean(aph);
        const std::vector<T> Bsub(2, mean);
        const Maph2mFitResult<T> r = maph2m_fit_multiclass(aph, p, Bsub, cw);
        adopt(r.maph);
        out.warning =
            "mamap22_fit_fs_multiclass: the canonical phase-type branch fits BACKWARD moments and "
            "had no forward targets to use, so both were set to the ordinary mean; fit with "
            "mamap22_fit_bs_multiclass instead if the forward moments matter";
        return out;
    }
    if (std::fabs(1.0 - dv(r1)) < degentol) {
        // Non-canonical: only the forward moment is identifiable, through a
        // one-variable inverse whose feasible set is an interval.
        auto degen_forward = [&](const T& vF1, T* q1, T* q2, T* q3) {
            *q1 = p[0];  // meaningless when the form is truly degenerate
            *q2 = T(p[0] * (h2 - vF1) / (h1 * (r2 - one)));
            *q3 = T(p[0] * (h1 + h2 - vF1) / (h1 * r2));
        };
        T q1 = zero, q2 = zero, q3 = zero;
        degen_forward(F[0], &q1, &q2, &q3);
        if (!(feasible(q1) && feasible(q2) && feasible(q3))) {
            const T q2F = T(-p[0] / (h1 * (r2 - one)));
            const T q2_0 = T(p[0] * h2 / (h1 * (r2 - one)));
            const T q3F = T(-p[0] / (h1 * r2));
            const T q3_0 = T(p[0] * (h1 + h2) / (h1 * r2));
            double lo = 1e-6, hi = 1e6;
            const double coefs[2] = {dv(q2F), dv(q3F)};
            const double offs[2] = {dv(q2_0), dv(q3_0)};
            for (int i = 0; i < 2; ++i) {
                if (std::fabs(coefs[i]) < denumtol) continue;
                const double a = -offs[i] / coefs[i], b = (1.0 - offs[i]) / coefs[i];
                lo = std::max(lo, std::min(a, b));
                hi = std::min(hi, std::max(a, b));
            }
            if (!(lo <= hi))
                throw NumericError(
                    "mamap22_fit_fs_multiclass: the non-canonical forward fit has an empty "
                    "feasible interval for this (p, F)");
            double x = dv(F[0]);
            if (x < lo) x = lo;
            if (x > hi) x = hi;
            degen_forward(num_traits<T>::from_double(x), &q1, &q2, &q3);
        }
        finish(q1, q2, q3);
        return out;
    }
    if (form == 2 && dv(r2) < degentol) {
        auto degen_forward2 = [&](const T& vF1, T* q1, T* q2, T* q3) {
            *q1 = T(p[0] * (r1 - two) * (h1 + h2 * r1 - vF1) / ((r1 - one) * (h1 - h2 + h2 * r1)));
            *q2 = T(-p[0] * (vF1 - h2) * (r1 - two) / (h1 - h2 + h2 * r1));
            *q3 = p[0];
        };
        T q1 = zero, q2 = zero, q3 = zero;
        if (dv(fw[0]) > dv(fw[1])) {
            degen_forward2(F[0], &q1, &q2, &q3);
            if (!(feasible(q1) && feasible(q2) && feasible(q3))) {
                const T q1F = T(-p[0] * (r1 - two) / ((r1 - one) * (h1 - h2 + h2 * r1)));
                const T q1_0 =
                    T(p[0] * (r1 - two) * (h1 + h2 * r1) / ((r1 - one) * (h1 - h2 + h2 * r1)));
                const T q2F = T(-p[0] * (r1 - two) / (h1 - h2 + h2 * r1));
                const T q2_0 = T(p[0] * (r1 - two) * h2 / (h1 - h2 + h2 * r1));
                double lo = 1e-6, hi = 1e6;
                const double coefs[2] = {dv(q1F), dv(q2F)};
                const double offs[2] = {dv(q1_0), dv(q2_0)};
                for (int i = 0; i < 2; ++i) {
                    if (std::fabs(coefs[i]) < denumtol) continue;
                    const double a = -offs[i] / coefs[i], b = (1.0 - offs[i]) / coefs[i];
                    lo = std::max(lo, std::min(a, b));
                    hi = std::min(hi, std::max(a, b));
                }
                if (!(lo <= hi))
                    throw NumericError(
                        "mamap22_fit_fs_multiclass: the degenerate forward fit has an empty "
                        "feasible interval for this (p, F)");
                double x = dv(F[0]);
                if (x < lo) x = lo;
                if (x > hi) x = hi;
                degen_forward2(num_traits<T>::from_double(x), &q1, &q2, &q3);
            }
            finish(q1, q2, q3);
            return out;
        }
        // The sigma-weighted arm. The reference states this repair as a YALMIP
        // program, but it is a ONE-VARIABLE CONVEX QP -- minimize (x - S11)^2
        // over {0 <= x <= p1^2, x >= p1^2(1-(1-r1)^2), x >= p1^2-(1-p1)^2} --
        // so the projection below IS its global optimum, not a local stand-in.
        // sqrt(p1^2 - S11) is COMPLEX above p1^2, where the reference's own
        // feasibility test then fails on the real part and sends it to the
        // repair. Report that as "not feasible" instead of raising, so the
        // repair below is reachable from the same inputs.
        auto degen_transition = [&](const T& vS11, T* a, T* b, T* c2) {
            const double root = dv(T(p[0] * p[0] - vS11));
            const T s = num_traits<T>::from_double(std::sqrt(std::max(root, 0.0)));
            *a = T(p[0] + s / (r1 - one));
            *b = T(p[0] + s);
            *c2 = p[0];
            return root >= 0.0;
        };
        const bool real = degen_transition(S(0, 0), &q1, &q2, &q3);
        if (!real || !(feasible(q1) && feasible(q2) && feasible(q3))) {
            const double p1 = dv(p[0]), rr1 = dv(r1);
            double lo = std::max(0.0, std::max(p1 * p1 * (1.0 - (1.0 - rr1) * (1.0 - rr1)),
                                               p1 * p1 - (1.0 - p1) * (1.0 - p1)));
            const double hi = p1 * p1;
            if (lo > hi)
                throw NumericError(
                    "mamap22_fit_fs_multiclass: the sigma repair of the gamma < 0 degeneracy is "
                    "infeasible for this (p, r1)");
            double s11 = dv(S(0, 0));
            if (s11 < lo) s11 = lo;
            if (s11 > hi) s11 = hi;
            degen_transition(num_traits<T>::from_double(s11), &q1, &q2, &q3);
        }
        finish(q1, q2, q3);
        return out;
    }

    // ---- the full form ---------------------------------------------------
    const Mamap2mCoefficients<T> c = form == 1
                                         ? mamap2m_can1_coefficients(h1, h2, r1, r2)
                                         : mamap2m_can2_coefficients(h1, h2, r1, r2);
    T q1 = zero, q2 = zero, q3 = zero;
    if (form == 1)
        fsdetail::fit_can1(c, p[0], F[0], S(0, 0), denumtol, &q1, &q2, &q3);
    else
        fsdetail::fit_can2(c, p[0], F[0], S(0, 0), denumtol, &q1, &q2, &q3);

    if (feasible(q1) && feasible(q2) && feasible(q3)) {
        out.exact = true;
        finish(q1, q2, q3);
        return out;
    }
    if (!adjust) {
        finish(q1, q2, q3);
        return out;
    }
    throw UnsupportedError(
        "mamap22_fit_fs_multiclass: the closed-form forward-plus-sigma inverse is infeasible for "
        "these targets, and the reference's repair solves a NONCONVEX bilinear program with "
        "YALMIP's bmibnb, a spatial branch-and-bound returning a GLOBAL optimum. That solver is "
        "not ported; a local method would report a different fit under the same name. Pass "
        "adjust = false to take the clamped closed form, or use mamap2m_fit_fb_multiclass");
}


/**
 * `mamap22_fit_gamma_fs`: fit over every AMAP(2) form and keep the closest.
 *
 * Port of matlab/lib/m3a/m3a/mamap22/mamap22_fit_gamma_fs.m. When the moment
 * set admits only a one-state process, the reference perturbs the second and
 * third moments slightly above the exponential to recover a two-state form, and
 * falls back to a marked Poisson only if that also fails; both steps are here.
 */
template <class T>
Mmap<T> mamap22_fit_gamma_fs(const T& M1, const T& M2, const T& M3, const T& GAMMA,
                             const std::vector<T>& p, const std::vector<T>& F,
                             const Matrix<T>& S) {
    const T one = num_traits<T>::from_int(1);
    Amap2FitGammaResult<T> a = amap2_fit_gamma(M1, M2, M3, GAMMA);
    if (a.amaps.size() == 1 && a.amaps[0].order() == 1) {
        // Perturb just above the exponential to recover a second-order form.
        const T M2a = T(M2 * (one + num_traits<T>::from_double(1e-4)));
        const T ratio = T(M2a / M2);
        const T M3a = T(M3 * num_traits<T>::from_double(
                                 std::pow(num_traits<T>::to_double(ratio), 1.5)));
        const std::vector<Map<T>> alt = amap2_fitall_gamma(M1, M2a, M3a, GAMMA);
        if (!alt.empty()) {
            a.amaps.clear();
            for (std::size_t j = 0; j < alt.size(); ++j) a.amaps.push_back(map_normalize(alt[j]));
        } else {
            return mamapdetail::marked_poisson(M1, p);
        }
    }

    Mmap<T> best;
    double bestErr = 0.0;
    bool have = false;
    for (std::size_t j = 0; j < a.amaps.size(); ++j) {
        Mamap22FsFitResult<T> r;
        try {
            r = mamap22_fit_fs_multiclass(a.amaps[j], p, F, S);
        } catch (const Error&) {
            continue;
        }
        // the reference scores on the FIRST class alone: fF(1) and fS(1,1)
        const double df = num_traits<T>::to_double(T(F[0] / r.fF[0])) - 1.0;
        const double ds = num_traits<T>::to_double(T(S(0, 0) / r.fS(0, 0))) - 1.0;
        const double err = df * df + ds * ds;
        if (!have || err < bestErr) {
            bestErr = err;
            best = r.mmap;
            have = true;
        }
    }
    if (!have)
        throw NumericError(
            "mamap22_fit_gamma_fs: no AMAP(2) form admits a feasible forward-plus-sigma marking "
            "for these targets");
    return best;
}

/** `mamap22_fit_gamma_fs` driven from a marked trace. */
template <class T>
Mmap<T> mamap22_fit_gamma_fs_trace(const std::vector<T>& Tv, const std::vector<int>& A) {
    if (Tv.empty() || Tv.size() != A.size())
        throw InputError(
            "mamap22_fit_gamma_fs_trace: the trace and its labels must agree in length");
    const T zero = num_traits<T>::from_int(0);
    T m1 = zero, m2 = zero, m3 = zero;
    for (std::size_t i = 0; i < Tv.size(); ++i) {
        const T x = Tv[i];
        m1 += x;
        m2 += x * x;
        m3 += x * x * x;
    }
    const T n = num_traits<T>::from_int(static_cast<long>(Tv.size()));
    const std::vector<unsigned> ord(1, 1u);
    const std::vector<T> p = trace::mtrace_pc<T>(A);
    const Matrix<T> fm = trace::mtrace_forward_moment(Tv, A, ord);
    std::vector<T> F(p.size(), zero);
    for (std::size_t c = 0; c < p.size(); ++c) F[c] = fm(c, 0);
    const Matrix<T> S = trace::mtrace_sigma<T>(A);
    return mamap22_fit_gamma_fs(T(m1 / n), T(m2 / n), T(m3 / n),
                                line::trace::trace_gamma(Tv).gamma, p, F, S);
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MAMAP22_FIT_FS_H
