/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MAPH2M_FIT_H
#define LINE_API_MAM_MAPH2M_FIT_H

/**
 * Fit a MAPH(2,m): a second-order acyclic phase-type marked with m classes.
 *
 * Templated port of matlab/lib/m3a/m3a/maph2/maph2m_fit.m,
 * maph2m_fit_multiclass.m, maph2m_fit_mmap.m and maph2m_fit_trace.m.
 *
 * The construction separates the TIMING from the MARKING. An APH(2) in
 * canonical acyclic form fixes the inter-arrival law from (M1, M2, M3); the
 * class marking then splits each of its two exit flows among the m classes with
 * probabilities q(j,c), and those are the only free parameters left. Writing
 * h1, h2 for the two phase means and r1 for the branch probability out of phase
 * one, the per-class BACKWARD moment is affine in the split, so
 *
 *   q(j,c) = fB(c) q_b(j,c) + q_0(j,c),
 *
 * with q_b and q_0 the coefficients the reference derives. The fit is then a
 * QUADRATIC PROGRAM in the achieved backward moments fB: minimize
 * sum_c w(c) (fB(c)/B(c) - 1)^2 subject to each q(j,.) being a probability
 * vector. Its Hessian is diagonal and positive, so the program is convex and
 * its unconstrained minimizer is fB = B; the constraints are what make the
 * answer differ from the target.
 *
 * THE DEGENERATE FORM HAS NO FREEDOM AT ALL. When r1 = 1 the second phase is
 * unreachable except through the first, both exit flows see the same class law,
 * and the only thing that can be matched is the class probability vector p.
 * The reference detects that at |1 - r1| < 1e-6 and sets q(1,c) = q(2,c) = p(c)
 * without solving anything; reproduced here, because solving the program on a
 * singular coefficient set returns whatever the solver's regularization
 * happens to give.
 *
 * THE SOLVER IS NOT quadprog. MATLAB runs an interior-point QP; this port uses
 * `line/util/auglag.h`, whose header states the acceptance contract for exactly
 * this substitution. The acceptance here is the specification: the fitted MAPH
 * reproduces p exactly (it is an equality constraint), its inter-arrival moments
 * are the APH's, and its backward moments approach B as closely as the
 * feasibility of the split allows. Iterates and multipliers are NOT comparable
 * with MATLAB's.
 *
 * `maph2m_fit` runs the whole thing once per APH(2) form that `aph2_fit`
 * returns and keeps the one whose backward moments land closest, which is why
 * the fitter needs `aph2_fit`'s full list and not just its selected form.
 *
 * ARITHMETIC: transcendental, through aph2_fit and the solver.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/mam/aph2_fit.h"
#include "line/api/mam/map_moment.h"
#include "line/api/mam/mmap_compress.h"
#include "line/api/mam/mmap_lambda.h"
#include "line/api/trace/mtrace_backward_moment.h"
#include "line/api/trace/mtrace_pc.h"
#include "line/num/number.h"
#include "line/util/auglag.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/** The fitted MAPH and the backward moments it actually achieved. */
template <class T>
struct Maph2mFitResult {
    Mmap<T> maph;
    std::vector<T> fB;  ///< achieved per-class backward moments
};

/**
 * Mark a canonical acyclic APH(2) with m classes.
 *
 * @param aph           the APH(2), in canonical acyclic form
 * @param p             per-class probabilities, summing to one
 * @param B             per-class target backward moments
 * @param classWeights  per-class weights in the objective; empty means uniform
 */
template <class T>
Maph2mFitResult<T> maph2m_fit_multiclass(const Map<T>& aph, const std::vector<T>& p,
                                         const std::vector<T>& B,
                                         const std::vector<T>& classWeights = std::vector<T>()) {
    static_assert(num_traits<T>::has_transcendental,
                  "maph2m_fit_multiclass solves a quadratic program");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    if (aph.D0.rows() != 2) throw InputError("maph2m_fit_multiclass: the APH must be second order");
    if (!(num_abs(T(aph.D0(1, 0))) <= zero))
        throw InputError("maph2m_fit_multiclass: the APH must be acyclic");
    if (!(num_abs(T(aph.D1(0, 1))) <= zero) || !(num_abs(T(aph.D1(1, 1))) <= zero))
        throw InputError("maph2m_fit_multiclass: the APH must be in canonical acyclic form");

    const std::size_t k = p.size();
    if (k == 0) throw InputError("maph2m_fit_multiclass: no classes given");
    if (B.size() != k)
        throw InputError("maph2m_fit_multiclass: one backward moment per class is required");
    std::vector<T> w = classWeights;
    if (w.empty()) w.assign(k, one);
    if (w.size() != k) throw InputError("maph2m_fit_multiclass: one weight per class is required");

    const T h1 = T(-one / aph.D0(0, 0));
    const T h2 = T(-one / aph.D0(1, 1));
    const T r1 = T(aph.D0(0, 1) * h1);

    Maph2mFitResult<T> out;
    out.maph.D0 = aph.D0;
    out.maph.D1 = aph.D1;
    out.maph.Dc.assign(k, Matrix<T>(2, 2, zero));

    std::vector<std::vector<T>> q(2, std::vector<T>(k, zero));
    bool solved = false;

    const double degentol = 1e-6, feastol = 1e-8;
    if (std::fabs(1.0 - num_traits<T>::to_double(r1)) < degentol) {
        // Degenerate: one degree of freedom, so only the class probabilities.
        for (std::size_t c = 0; c < k; ++c) {
            q[0][c] = p[c];
            q[1][c] = p[c];
        }
    } else {
        std::vector<std::vector<T>> qb(2, std::vector<T>(k, zero));
        std::vector<std::vector<T>> q0(2, std::vector<T>(k, zero));
        for (std::size_t c = 0; c < k; ++c) {
            qb[0][c] = T(p[c] * (one / (h2 * (r1 - one))));
            q0[0][c] = T(p[c] * (-(h1 + h2) / (h2 * (r1 - one))));
            qb[1][c] = T(p[c] * (one / (h2 * r1)));
            q0[1][c] = T(p[c] * (-h1 / (h2 * r1)));
        }

        // min sum_c w(c) (x(c)/B(c) - 1)^2, i.e. the reference's (1/2)x'Hx+h'x
        // shifted by a constant, over the two equality and 4k inequality rows.
        auto fobj = [&](const std::vector<T>& x) {
            T s = zero;
            for (std::size_t c = 0; c < k; ++c) {
                const T r = T(x[c] / B[c] - one);
                s += w[c] * r * r;
            }
            return s;
        };
        auto heq = [&](const std::vector<T>& x) {
            std::vector<T> v(2, zero);
            for (std::size_t j = 0; j < 2; ++j) {
                T s = zero;
                for (std::size_t c = 0; c < k; ++c) s += qb[j][c] * x[c] + q0[j][c];
                v[j] = T(s - one);  // each row of q must sum to one
            }
            return v;
        };
        auto gineq = [&](const std::vector<T>& x) {
            std::vector<T> v;
            v.reserve(4 * k);
            for (std::size_t c = 0; c < k; ++c)
                for (std::size_t j = 0; j < 2; ++j) {
                    const T qq = T(qb[j][c] * x[c] + q0[j][c]);
                    v.push_back(T(qq - one));  // q <= 1
                    v.push_back(-qq);          // q >= 0
                }
            return v;
        };

        std::vector<T> x0 = B;
        std::vector<Bound<T>> bounds(k);
        for (std::size_t c = 0; c < k; ++c) {
            bounds[c].lo = num_traits<T>::from_double(1e-6);
            bounds[c].hi = num_traits<T>::from_double(1e6);
            if (x0[c] < bounds[c].lo) x0[c] = bounds[c].lo;
            if (x0[c] > bounds[c].hi) x0[c] = bounds[c].hi;
        }
        const AugLagResult<T> r = auglag(fobj, heq, gineq, x0, bounds);
        out.fB = r.x;
        solved = true;
        for (std::size_t c = 0; c < k; ++c)
            for (std::size_t j = 0; j < 2; ++j) q[j][c] = T(r.x[c] * qb[j][c] + q0[j][c]);
    }

    // The reference's feasibility gate, then its clamp-and-renormalize.
    for (std::size_t j = 0; j < 2; ++j) {
        T lo = q[j][0], s = zero;
        for (std::size_t c = 0; c < k; ++c) {
            if (q[j][c] < lo) lo = q[j][c];
            s += q[j][c];
        }
        if (num_traits<T>::to_double(lo) < -feastol ||
            num_traits<T>::to_double(s) > 1.0 + feastol)
            throw NumericError(
                "maph2m_fit_multiclass: feasibility could not be restored; the requested class "
                "probabilities and backward moments admit no valid split of the APH(2) exit flows");
    }
    for (std::size_t j = 0; j < 2; ++j) {
        T s = zero;
        for (std::size_t c = 0; c < k; ++c) {
            if (q[j][c] < zero) q[j][c] = zero;
            s += q[j][c];
        }
        if (!(num_traits<T>::to_double(s) > 0.0))
            throw NumericError("maph2m_fit_multiclass: a split lost all its mass");
        for (std::size_t c = 0; c < k; ++c) q[j][c] = T(q[j][c] / s);
    }

    // Dc = D1 .* [q(1,c) 0; q(2,c) 0]: only the first column of D1 is non-zero
    // in canonical acyclic form, so the split acts on the restart flow alone.
    for (std::size_t c = 0; c < k; ++c) {
        out.maph.Dc[c](0, 0) = T(out.maph.D1(0, 0) * q[0][c]);
        out.maph.Dc[c](1, 0) = T(out.maph.D1(1, 0) * q[1][c]);
    }
    if (!solved) {
        const std::vector<std::vector<T>> bm =
            mmap_backward_moment(out.maph, std::vector<unsigned>(1, 1u), true);
        out.fB.assign(k, zero);
        for (std::size_t c = 0; c < k; ++c) out.fB[c] = bm[c][0];
    }
    return out;
}

/**
 * Fit a MAPH(2,m) to three moments, the class probabilities and the per-class
 * backward moments, trying every APH(2) form and keeping the closest.
 */
template <class T>
Mmap<T> maph2m_fit(const T& M1, const T& M2, const T& M3, const std::vector<T>& p,
                   const std::vector<T>& B) {
    const Aph2FitResult<T> a = aph2_fit(M1, M2, M3);
    if (a.aphs.empty()) throw NumericError("maph2m_fit: no APH(2) fits the given moments");
    Mmap<T> best;
    double bestErr = 0.0;
    bool have = false;
    for (std::size_t j = 0; j < a.aphs.size(); ++j) {
        Maph2mFitResult<T> r;
        try {
            r = maph2m_fit_multiclass(a.aphs[j], p, B);
        } catch (const Error&) {
            continue;  // this form admits no valid split; the next one may
        }
        double err = 0.0;
        for (std::size_t c = 0; c < p.size(); ++c) {
            const double d = num_traits<T>::to_double(T(r.fB[c] / B[c])) - 1.0;
            err += d * d;
        }
        if (!have || err < bestErr) {
            bestErr = err;
            best = r.maph;
            have = true;
        }
    }
    if (!have)
        throw NumericError(
            "maph2m_fit: no APH(2) form admits a valid class split for the requested class "
            "probabilities and backward moments");
    return best;
}

/** Fit a MAPH(2,m) to the descriptors measured on a marked MAP. */
template <class T>
Mmap<T> maph2m_fit_mmap(const Mmap<T>& m) {
    const std::vector<T> p = mmap_pc(m);
    const std::vector<std::vector<T>> bm =
        mmap_backward_moment(m, std::vector<unsigned>(1, 1u), true);
    std::vector<T> B(p.size(), num_traits<T>::from_int(0));
    for (std::size_t c = 0; c < p.size(); ++c) B[c] = bm[c][0];
    return maph2m_fit(map_moment(m.map(), 1), map_moment(m.map(), 2), map_moment(m.map(), 3), p, B);
}

/**
 * Fit a MAPH(2,m) to the descriptors measured on a marked trace.
 *
 * @param Tv the inter-arrival times
 * @param A  the class of each arrival, 1-based as the reference indexes them
 */
template <class T>
Mmap<T> maph2m_fit_trace(const std::vector<T>& Tv, const std::vector<int>& A) {
    if (Tv.empty() || Tv.size() != A.size())
        throw InputError("maph2m_fit_trace: the trace and its labels must agree in length");
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
    const Matrix<T> bm = trace::mtrace_backward_moment(Tv, A, std::vector<unsigned>(1, 1u));
    std::vector<T> B(p.size(), zero);
    for (std::size_t c = 0; c < p.size(); ++c) B[c] = bm(c, 0);
    return maph2m_fit(T(m1 / n), T(m2 / n), T(m3 / n), p, B);
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MAPH2M_FIT_H
