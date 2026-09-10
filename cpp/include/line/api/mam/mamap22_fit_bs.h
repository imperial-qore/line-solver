/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MAMAP22_FIT_BS_H
#define LINE_API_MAM_MAMAP22_FIT_BS_H

/**
 * Fit a MAMAP(2,2) matching the BACKWARD moment and the class TRANSITION
 * probability sigma.
 *
 * Templated port of matlab/lib/m3a/m3a/mamap22/mamap22_fit_bs_multiclass.m.
 *
 * TWO CLASSES ONLY, by construction: sigma is a class-to-class transition
 * probability, and the closed forms below invert it for the (1,1) entry alone.
 * The reference says so and refuses otherwise, and so does this.
 *
 * Unlike the F+B fitter, which is affine in the marking and therefore a
 * quadratic program, the B+S relation is a RATIO: q2 (form 1) or q3 (form 2) is
 * a quadratic in the backward moment over a linear denominator, and the other
 * two follow linearly from it. That inverse is closed form -- `fit_can1` and
 * `fit_can2` here -- and it is EXACT whenever its answer lands in the unit box.
 * The coefficient tables it reads are `mamap2m_coefficients.h`.
 *
 * WHAT IS NOT PORTED, and is refused by name rather than approximated. When the
 * closed form lands outside the box and the caller asked for `adjust`, the
 * reference repairs it by solving a NONCONVEX program -- the marking equality
 * q2 * den(B1) = num(B1, S11) is bilinear in the unknowns -- with YALMIP's
 * `bmibnb`, a spatial branch-and-bound that returns a GLOBAL optimum, over two
 * sides (B1 below and above the mean) and keeps the better. `line/util/auglag.h`
 * finds a local KKT point, which on a nonconvex feasible set is a different
 * answer, not a slower one; substituting it would report a fit the reference
 * would not have chosen. Every other branch IS ported, including all four
 * degeneracies and the one-variable quadratic repair of the degenerate MMAP
 * form, so the exact path and the degenerate paths are complete.
 *
 * ARITHMETIC: transcendental, through the fitters it delegates to.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/mam/amap2_fit_gamma.h"
#include "line/api/mam/m3a_fit_from.h"
#include "line/api/mam/mamap2m_coefficients.h"
#include "line/api/mam/mamap_marked_poisson.h"
#include "line/api/mam/map_gamma.h"
#include "line/api/mam/map_moment.h"
#include "line/api/trace/mtrace_backward_moment.h"
#include "line/api/trace/mtrace_pc.h"
#include "line/api/trace/mtrace_sigma.h"
#include "line/api/trace/trace_gamma.h"
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

/** The fitted MAMAP(2,2), what it achieved, and whether the fit was exact. */
template <class T>
struct Mamap22FitResult {
    Mmap<T> mmap;
    std::vector<T> fB;  ///< achieved per-class backward moments
    Matrix<T> fS;       ///< achieved class transition probabilities
    bool exact = false;
};

namespace bsdetail {

/** The closed-form B+S inverse of the first canonical form. */
template <class T>
void fit_can1(const Mamap2mCoefficients<T>& c, const T& p1, const T& vB1, const T& vS11,
              double denumtol, T* q1, T* q2, T* q3) {
    const T den = T(c.U[10] * vB1 * p1 + c.U[11] * p1);
    if (std::fabs(num_traits<T>::to_double(den)) < denumtol) {
        *q1 = *q2 = *q3 = p1;  // the inverse degenerates; only p is identifiable
        return;
    }
    *q2 = T((c.U[6] * vB1 * vB1 * p1 * p1 + c.U[7] * vB1 * p1 * p1 + c.U[8] * vS11 +
             c.U[9] * p1 * p1) /
            den);
    *q1 = T(-(c.G[11] * p1 - vB1 * c.G[2] * p1 +
              (c.G[2] * c.G[10] - c.G[1] * c.G[11]) * (*q2)) /
            c.Y[2]);
    *q3 = T((c.G[9] * p1 - vB1 * c.G[0] * p1 +
             (c.G[0] * c.G[10] - c.G[1] * c.G[9]) * (*q2)) /
            c.Y[2]);
}

/** The closed-form B+S inverse of the second canonical form. */
template <class T>
void fit_can2(const Mamap2mCoefficients<T>& c, const T& p1, const T& vB1, const T& vS11,
              double denumtol, T* q1, T* q2, T* q3) {
    const T den = T(c.U[10] * vB1 * p1 + c.U[11] * p1);
    if (std::fabs(num_traits<T>::to_double(den)) < denumtol) {
        *q1 = *q2 = *q3 = p1;
        return;
    }
    *q3 = T((c.U[6] * vB1 * vB1 * p1 * p1 + c.U[7] * vB1 * p1 * p1 + c.U[8] * p1 * p1 +
             c.U[9] * vS11) /
            den);
    *q1 = T((c.G[9] * p1 - vB1 * c.G[1] * p1 +
             (c.G[1] * c.G[10] - c.G[2] * c.G[9]) * (*q3)) /
            c.Y[2]);
    *q2 = T(-(c.G[8] * p1 - vB1 * c.G[0] * p1 +
              (c.G[0] * c.G[10] - c.G[2] * c.G[8]) * (*q3)) /
            c.Y[2]);
}

}  // namespace bsdetail

/**
 * @param map          the AMAP(2), in one of the two canonical acyclic forms
 * @param p            the two class probabilities
 * @param B            the two target backward moments
 * @param S            the target class transition matrix; only S(0,0) is used
 * @param classWeights per-class weights; empty means uniform
 * @param bsWeights    the (backward, sigma) weights; empty means (1, 1)
 * @param adjust       repair an infeasible closed form; the repair is unported
 */
template <class T>
Mamap22FitResult<T> mamap22_fit_bs_multiclass(const Map<T>& map, const std::vector<T>& p,
                                              const std::vector<T>& B, const Matrix<T>& S,
                                              const std::vector<T>& classWeights = std::vector<T>(),
                                              const std::vector<T>& bsWeights = std::vector<T>(),
                                              bool adjust = true) {
    static_assert(num_traits<T>::has_transcendental,
                  "mamap22_fit_bs_multiclass inverts a moment system");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    if (map.D0.rows() != 2)
        throw InputError("mamap22_fit_bs_multiclass: the underlying MAP must be second order");
    if (!(num_abs(T(map.D0(1, 0))) <= zero))
        throw InputError("mamap22_fit_bs_multiclass: the underlying MAP must be acyclic");
    int form;
    if (map.D1(0, 1) == zero) form = 1;
    else if (map.D1(0, 0) == zero) form = 2;
    else
        throw InputError(
            "mamap22_fit_bs_multiclass: the underlying MAP must be in canonical acyclic form");
    if (p.size() != 2)
        throw InputError(
            "mamap22_fit_bs_multiclass: fitting the backward moment and the transition "
            "probabilities supports two classes only");
    if (B.size() != 2)
        throw InputError("mamap22_fit_bs_multiclass: one backward moment per class is required");
    if (S.rows() < 1 || S.cols() < 1)
        throw InputError("mamap22_fit_bs_multiclass: the transition matrix is empty");

    std::vector<T> cw = classWeights;
    if (cw.empty()) cw.assign(2, one);
    std::vector<T> bw = bsWeights;
    if (bw.empty()) bw.assign(2, one);

    const double degentol = 1e-6, feastol = 1e-4, denumtol = 1e-12;
    Map<T> mp = map;
    T h1 = T(-one / mp.D0(0, 0)), h2 = T(-one / mp.D0(1, 1));
    T r1 = T(mp.D0(0, 1) * h1), r2 = T(mp.D1(1, 1) * h2);
    auto dv = [](const T& v) { return num_traits<T>::to_double(v); };

    Mamap22FitResult<T> out;
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
            out.mmap.Dc[0](1, 0) = T(out.mmap.D1(1, 0) * q2);
            out.mmap.Dc[0](1, 1) = T(out.mmap.D1(1, 1) * q3);
            out.mmap.Dc[1](0, 0) = T(out.mmap.D1(0, 0) * (one - q1));
            out.mmap.Dc[1](1, 0) = T(out.mmap.D1(1, 0) * (one - q2));
            out.mmap.Dc[1](1, 1) = T(out.mmap.D1(1, 1) * (one - q3));
        } else {
            out.mmap.Dc[0](0, 1) = T(out.mmap.D1(0, 1) * q1);
            out.mmap.Dc[0](1, 0) = T(out.mmap.D1(1, 0) * q2);
            out.mmap.Dc[0](1, 1) = T(out.mmap.D1(1, 1) * q3);
            out.mmap.Dc[1](0, 1) = T(out.mmap.D1(0, 1) * (one - q1));
            out.mmap.Dc[1](1, 0) = T(out.mmap.D1(1, 0) * (one - q2));
            out.mmap.Dc[1](1, 1) = T(out.mmap.D1(1, 1) * (one - q3));
        }
        const std::vector<unsigned> ord(1, 1u);
        const std::vector<std::vector<T>> bm = mmap_backward_moment(out.mmap, ord, true);
        out.fB.assign(2, zero);
        for (std::size_t c = 0; c < 2; ++c) out.fB[c] = bm[c][0];
        out.fS = mmap_sigma(out.mmap);
    };
    auto feasible = [&](const T& q) {
        return dv(q) >= -feastol && dv(q) <= 1.0 + feastol;
    };

    // ---- the Poisson perturbation, which keeps the two-state structure ----
    const bool degen1 = form == 1 && (dv(r1) < degentol || dv(r2) > 1.0 - degentol ||
                                      std::fabs(dv(h2) - dv(h1) * dv(r2)) < degentol);
    const bool degen2 =
        form == 2 && (dv(r2) > 1.0 - degentol ||
                      std::fabs(dv(h1) - dv(h2) - dv(h1) * dv(r1) + dv(h1) * dv(r1) * dv(r2)) <
                          degentol);
    if (degen1 || degen2) {
        if (dv(r1) < degentol) r1 = num_traits<T>::from_double(degentol);
        if (dv(r2) > 1.0 - degentol) r2 = num_traits<T>::from_double(1.0 - degentol);
        if (form == 1 && std::fabs(dv(h2) - dv(h1) * dv(r2)) < degentol)
            h2 = T(h1 * r2 + num_traits<T>::from_double(degentol));
        if (form == 2 &&
            std::fabs(dv(h1) - dv(h2) - dv(h1) * dv(r1) + dv(h1) * dv(r1) * dv(r2)) < degentol)
            h1 = T((h2 + num_traits<T>::from_double(degentol)) / (one - r1 + r1 * r2));
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

    // ---- the degenerate ladder, in the reference's order ------------------
    if (form == 2 && dv(r2) < degentol && std::fabs(1.0 - dv(r1)) < degentol) {
        finish(p[0], p[0], p[0]);  // only the class probabilities are identifiable
        return out;
    }
    if (form == 1 && dv(r2) < degentol) {
        // A canonical APH(2): the problem IS the MAPH fit.
        Map<T> aph = mp;
        aph.D1(1, 1) = zero;
        aph = map_normalize(aph);
        const Maph2mFitResult<T> r = maph2m_fit_multiclass(aph, p, B, cw);
        out.mmap = r.maph;
        const std::vector<unsigned> ord(1, 1u);
        const std::vector<std::vector<T>> bm = mmap_backward_moment(out.mmap, ord, true);
        out.fB.assign(2, zero);
        for (std::size_t c = 0; c < 2; ++c) out.fB[c] = bm[c][0];
        out.fS = mmap_sigma(out.mmap);
        return out;
    }
    if (std::fabs(1.0 - dv(r1)) < degentol) {
        // Non-canonical: refit the timing as an APH(2), then mark it.
        const Map<T> aph = aph2_fit_map(mp).aph;
        const Maph2mFitResult<T> r = maph2m_fit_multiclass(aph, p, B, cw);
        out.mmap = r.maph;
        const std::vector<unsigned> ord(1, 1u);
        const std::vector<std::vector<T>> bm = mmap_backward_moment(out.mmap, ord, true);
        out.fB.assign(2, zero);
        for (std::size_t c = 0; c < 2; ++c) out.fB[c] = bm[c][0];
        out.fS = mmap_sigma(out.mmap);
        return out;
    }
    if (form == 2 && dv(r2) < degentol) {
        // The gamma < 0 degeneracy: fit B or sigma, whichever is weighted higher.
        auto degen_backward = [&](const T& vB1, T* q1, T* q2, T* q3) {
            *q1 = T((p[0] * (r1 - two) * (h2 - vB1 + h1 * r1)) /
                    ((r1 - one) * (h2 - h1 + h1 * r1)));
            *q2 = T(-(p[0] * (vB1 - h1) * (r1 - two)) / (h2 - h1 + h1 * r1));
            *q3 = p[0];  // meaningless when the form is truly degenerate
        };
        T q1 = zero, q2 = zero, q3 = zero;
        if (dv(bw[0]) > dv(bw[1])) {
            degen_backward(B[0], &q1, &q2, &q3);
            if (!(feasible(q1) && feasible(q2) && feasible(q3))) {
                // One variable, four inequalities: the box the reference solves
                // a quadratic program over reduces to an interval, and the
                // minimizer of (x/B - 1)^2 on it is the projection of B onto it.
                const T q1B = T(-p[0] * (r1 - two) / ((r1 - one) * (h2 - h1 + h1 * r1)));
                const T q1_0 = T(p[0] * (r1 - two) * (h2 + h1 * r1) /
                                 ((r1 - one) * (h2 - h1 + h1 * r1)));
                const T q2B = T(-p[0] * (r1 - two) / (h2 - h1 + h1 * r1));
                const T q2_0 = T(p[0] * (r1 - two) * h1 / (h2 - h1 + h1 * r1));
                double lo = 1e-6, hi = 1e6;
                const double coefs[2] = {dv(q1B), dv(q2B)};
                const double offs[2] = {dv(q1_0), dv(q2_0)};
                for (int i = 0; i < 2; ++i) {
                    if (std::fabs(coefs[i]) < denumtol) continue;
                    // 0 <= coef x + off <= 1
                    const double a = -offs[i] / coefs[i], b = (1.0 - offs[i]) / coefs[i];
                    lo = std::max(lo, std::min(a, b));
                    hi = std::min(hi, std::max(a, b));
                }
                if (!(lo <= hi))
                    throw NumericError(
                        "mamap22_fit_bs_multiclass: the degenerate backward fit has an empty "
                        "feasible interval for this (p, B)");
                double x = dv(B[0]);
                if (x < lo) x = lo;
                if (x > hi) x = hi;
                degen_backward(num_traits<T>::from_double(x), &q1, &q2, &q3);
            }
        } else {
            // sqrt(p1^2 - S11) is COMPLEX above p1^2, where the reference's own
            // feasibility test then fails on the real part and sends it to the
            // clamp. Report that as "not feasible" instead of raising, so the
            // clamp below is reachable from the same inputs.
            auto degen_transition = [&](const T& vS11, T* a, T* b, T* c) {
                const double root = dv(T(p[0] * p[0] - vS11));
                const T s = num_traits<T>::from_double(std::sqrt(std::max(root, 0.0)));
                *a = T(p[0] + s / (r1 - one));
                *b = T(p[0] + s);
                *c = p[0];
                return root >= 0.0;
            };
            const bool real = degen_transition(S(0, 0), &q1, &q2, &q3);
            if (!real || !(feasible(q1) && feasible(q2) && feasible(q3))) {
                // The reference clamps S11 into the interval the form admits.
                const double safety = 1e-10;
                const double p1 = dv(p[0]), rr1 = dv(r1);
                const double lb = p1 * p1 * (1.0 - (1.0 - rr1) * (1.0 - rr1));
                const double ub = p1 * p1 - (1.0 - p1) * (1.0 - p1);
                double s11 = dv(S(0, 0));
                if (s11 <= lb) s11 = lb + safety;
                else if (s11 >= ub) s11 = ub - safety;
                degen_transition(num_traits<T>::from_double(s11), &q1, &q2, &q3);
            }
        }
        finish(q1, q2, q3);
        return out;
    }

    // ---- the full form: the closed-form inverse -------------------------
    const Mamap2mCoefficients<T> c = form == 1
                                         ? mamap2m_can1_coefficients(h1, h2, r1, r2)
                                         : mamap2m_can2_coefficients(h1, h2, r1, r2);
    T q1 = zero, q2 = zero, q3 = zero;
    if (form == 1)
        bsdetail::fit_can1(c, p[0], B[0], S(0, 0), denumtol, &q1, &q2, &q3);
    else
        bsdetail::fit_can2(c, p[0], B[0], S(0, 0), denumtol, &q1, &q2, &q3);

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
        "mamap22_fit_bs_multiclass: the closed-form backward-plus-sigma inverse is infeasible for "
        "these targets, and the reference's repair solves a NONCONVEX bilinear program with "
        "YALMIP's bmibnb, a spatial branch-and-bound returning a GLOBAL optimum. That solver is "
        "not ported; a local method would report a different fit under the same name. Pass "
        "adjust = false to take the clamped closed form, weight the forward moment above sigma so "
        "mamap2m_fit_fb_multiclass applies, or relax the targets");
}

/**
 * `mamap22_fit_gamma_bs`: fit over every AMAP(2) form and keep the closest.
 *
 * Port of matlab/lib/m3a/m3a/mamap22/mamap22_fit_gamma_bs.m: the class
 * probabilities are fitted exactly, the BACKWARD moments and the one-step
 * class transition probabilities approximately. When the moment set admits
 * only a one-state process, the reference perturbs the second and third
 * moments slightly above the exponential to recover a two-state form, and
 * falls back to a marked Poisson only if that also fails; both steps are
 * here, as in the `mamap22_fit_gamma_fs` twin.
 */
template <class T>
Mmap<T> mamap22_fit_gamma_bs(const T& M1, const T& M2, const T& M3, const T& GAMMA,
                             const std::vector<T>& p, const std::vector<T>& B,
                             const Matrix<T>& S) {
    const T one = num_traits<T>::from_int(1);
    Amap2FitGammaResult<T> a = amap2_fit_gamma(M1, M2, M3, GAMMA);
    if (a.amaps.size() == 1 && a.amaps[0].order() == 1) {
        // perturb just above the exponential to recover a second-order form
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
        Mamap22FitResult<T> r;
        try {
            r = mamap22_fit_bs_multiclass(a.amaps[j], p, B, S);
        } catch (const Error&) {
            continue;
        }
        // the reference scores on the FIRST class alone: fB(1) and fS(1,1)
        const double db = num_traits<T>::to_double(T(r.fB[0] / B[0])) - 1.0;
        const double ds = num_traits<T>::to_double(T(r.fS(0, 0) / S(0, 0))) - 1.0;
        const double err = db * db + ds * ds;
        if (!have || err < bestErr) {
            bestErr = err;
            best = r.mmap;
            have = true;
        }
    }
    if (!have)
        throw NumericError(
            "mamap22_fit_gamma_bs: no AMAP(2) form admits a feasible backward-plus-sigma "
            "marking for these targets");
    return best;
}

/** `mamap22_fit_gamma_bs` driven from an MMAP[2] of arbitrary order. */
template <class T>
Mmap<T> mamap22_fit_gamma_bs_mmap(const Mmap<T>& mm) {
    const Map<T> mp = mm.map();
    const std::vector<T> p = mmap_pc(mm);
    const std::vector<std::vector<T>> bm =
        mmap_backward_moment(mm, std::vector<unsigned>(1, 1u), true);
    std::vector<T> B(p.size());
    for (std::size_t c = 0; c < p.size(); ++c) B[c] = bm[c][0];
    return mamap22_fit_gamma_bs(map_moment(mp, 1u), map_moment(mp, 2u), map_moment(mp, 3u),
                                map_gamma(mp), p, B, mmap_sigma(mm));
}

/** `mamap22_fit_gamma_bs` driven from a marked trace. */
template <class T>
Mmap<T> mamap22_fit_gamma_bs_trace(const std::vector<T>& Tv, const std::vector<int>& A) {
    if (Tv.empty() || Tv.size() != A.size())
        throw InputError(
            "mamap22_fit_gamma_bs_trace: the trace and its labels must agree in length");
    const T zero = num_traits<T>::from_int(0);
    T m1 = zero, m2 = zero, m3 = zero;
    for (std::size_t i = 0; i < Tv.size(); ++i) {
        const T x = Tv[i];
        m1 += x;
        m2 += x * x;
        m3 += x * x * x;
    }
    const T n = num_traits<T>::from_int(static_cast<long>(Tv.size()));
    const std::vector<T> p = trace::mtrace_pc<T>(A);
    const Matrix<T> bm =
        trace::mtrace_backward_moment(Tv, A, std::vector<unsigned>(1, 1u));
    std::vector<T> B(p.size(), zero);
    for (std::size_t c = 0; c < p.size(); ++c) B[c] = bm(c, 0);
    const Matrix<T> S = trace::mtrace_sigma<T>(A);
    return mamap22_fit_gamma_bs(T(m1 / n), T(m2 / n), T(m3 / n),
                                line::trace::trace_gamma(Tv).gamma, p, B, S);
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MAMAP22_FIT_BS_H
