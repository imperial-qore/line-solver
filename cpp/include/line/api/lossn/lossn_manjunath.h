/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_LOSSN_LOSSN_MANJUNATH_H
#define LINE_API_LOSSN_LOSSN_MANJUNATH_H

/**
 * Exact analysis of a loss network by the Manjunath-Sikdar transform.
 *
 * Templated port of matlab/src/api/lossn/lossn_manjunath.m. Calls on route r arrive
 * Poisson at rate nu_r with unit mean holding time and are admitted only while
 * every constraint holds, sum_r A(j,r) n_r <= C(j). The admissible set is
 * coordinate convex, so Kelly's truncation theorem gives the truncated product
 * form p(n) = nu^n / n! / g(C) and every metric is a ratio of normalizing
 * constants:
 *
 *     g(C)   = sum_{A n <= C} prod_r nu_r^{n_r} / n_r!
 *     E[n_r] = nu_r g(C - A e_r) / g(C)
 *     Loss_r = 1 - g(C - A e_r) / g(C)
 *
 * because a class r call is blocked exactly when the state cannot absorb one
 * more unit of its own requirement vector.
 *
 * WHY IT IS A COEFFICIENT COMPUTATION AND NOT A QUADRATURE. Writing each
 * indicator as a contour integral turns g(C) into a J-fold integral over the
 * unit circle whose integrand factorizes into the per-route z-transforms. Inside
 * the circle the only pole in z_j sits at the origin with order C_j+1, so each
 * integration is a residue, i.e. a Taylor coefficient. The routine therefore
 * never evaluates an integral: it builds the generating function as a
 * multivariate power series truncated at degree C_j in z_j, one
 * shift-and-accumulate convolution per route, and discharges each '<='
 * constraint by summing the coefficients of degrees 0..C_j along that
 * dimension. Truncation is exact because A is nonnegative -- a monomial above
 * degree C_j can never contribute to an extracted coefficient.
 *
 * THE ELIMINATION ORDER IS THE MEMORY BOUND. Contour integrations are
 * interleaved with the product rather than deferred: variable z_j is created
 * when the first route with A(j,r) != 0 is multiplied in and integrated out
 * immediately after the last one. Peak memory is therefore the product of
 * (C_j+1) over the SIMULTANEOUSLY LIVE links, an induced width of the
 * route-link incidence, not over all J links. That product is bounded by
 * `LossnManjunathOptions::max_live_states` and a region above it is refused by name
 * rather than allowed to exhaust the machine: the algorithm is exact but not
 * unconditionally cheap, and `lossn_mci` answers the same question at any size.
 *
 * WHY THIS ONE RUNS AT EXACT ARITHMETIC AND ITS TWO SIBLINGS DO NOT.
 * `lossn_erlangfp` stops on a tolerance and `lossn_mci` returns a random
 * variable, so both are transcendental-gated. Here every operation on the series
 * is an addition or a multiplication of terms nu_r^n / n!, which is rational
 * whenever nu is, and the reported quantities are RATIOS of series values, so
 * QLen and Loss come out exact under Arith::Exact. Only `lG` is transcendental,
 * and it is a double in the result for all three arithmetics anyway, obtained
 * through `num_traits<T>::log_as_double` (which is defined for Rational via
 * log_bigint, so no logarithm is ever taken of a value that has to be
 * representable as a double).
 *
 * OVERFLOW, AND WHY THE SCALING IS ARITHMETIC-DEPENDENT. The term nu^n/n! peaks
 * near n = nu at roughly e^nu / sqrt(2 pi nu), which overflows a double for
 * loads above ~700. The reference builds the sequence in log space and divides
 * it by its largest entry, accumulating the discarded logarithm and adding it
 * back in log g; the port does the same under an inexact T. Under an exact T
 * there is no overflow to avoid and the scaling would only inflate the
 * denominators of the rationals, so it is skipped. The scale cancels in every
 * ratio, so QLen and Loss are unaffected and lG is identical either way.
 *
 * A and C must be integer valued, since the residue argument counts whole units
 * of capacity; a fractional entry is refused rather than rounded, naming
 * `lossn_mci` which compares in real arithmetic. Each row is divided by the
 * greatest common divisor of its entries together with its right-hand side,
 * which is exact and shrinks the truncation degree. Routes appearing in no
 * constraint never block and contribute a factor exp(nu_r) to g(C).
 *
 * Reference: D. Manjunath and B. Sikdar, Integral Expressions for the Numerical
 * Evaluation of Product Form Expressions Over Irregular Multidimensional
 * Integer Spaces.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace lossn {

/** Controls of lossn_manjunath. The reference has none; both fields are port-local. */
struct LossnManjunathOptions {
    /**
     * Cap on the product of (C_j+1) over the simultaneously live links, i.e. on
     * the number of series coefficients held at once. The default is 2^26
     * coefficients, half a gigabyte at double, which no region a finite capacity
     * region can express reaches by accident. Raise it deliberately.
     */
    std::size_t max_live_states = static_cast<std::size_t>(1) << 26;
};

/** Result of lossn_manjunath. */
template <class T>
struct LossnManjunathResult {
    std::vector<T> QLen;  ///< mean carried load E[n_r] per route
    std::vector<T> Loss;  ///< blocking probability per route
    double lG = 0.0;      ///< log of the EXACT normalizing constant g(C)
    /** Always 1: the transform is direct, and the field exists for the shared
     *  analyzer contract that the iterative siblings fill in. */
    std::size_t iterations = 1;
    /** Peak number of live series coefficients, the realised cost. */
    std::size_t peak_states = 0;
};

namespace detail {

/** gcd on nonnegative long, with gcd(0, a) = a as MATLAB's gcd has it. */
inline long lossn_gcd(long a, long b) {
    while (b != 0) {
        const long t = a % b;
        a = b;
        b = t;
    }
    return a < 0 ? -a : a;
}

/**
 * The reduced admission rule: rows that constrain nothing dropped, every
 * remaining row divided by the gcd of its entries and its right-hand side.
 *
 * Both steps are exact. The second is what makes the cost tractable on a memory
 * budget whose class sizes share a factor -- a row (4, 8) <= 20 becomes
 * (1, 2) <= 5 and its dimension shrinks from 21 coefficients to 6.
 */
struct LossnManjunathRule {
    std::vector<std::vector<long>> A;  ///< (J x R) after dropping and reduction
    std::vector<long> C;               ///< (J)
};

/**
 * Integer view of (A, C), refusing a fractional or negative entry by name.
 *
 * The tolerance is the reference's 1e-9 on the distance to the nearest integer.
 * It is applied to the value as a double even under exact arithmetic, where a
 * fractional entry is exactly representable: the question asked is whether the
 * caller MEANT an integer, and a rational 3/2 answers no just as 1.5 does.
 */
template <class T>
LossnManjunathRule lossn_manjunath_integralize(const Matrix<T>& A, const std::vector<T>& C) {
    const std::size_t J = C.size(), R = A.cols();
    LossnManjunathRule rule;
    for (std::size_t j = 0; j < J; ++j) {
        const double c = num_traits<T>::to_double(C[j]);
        if (c < 0.0 || std::fabs(c - std::round(c)) > 1e-9)
            throw InputError(
                "lossn_manjunath: C must contain nonnegative integers -- the residue argument counts "
                "whole units of capacity. Use lossn_mci, which compares in real arithmetic, or "
                "lossn_erlangfp");
        bool any = false;
        std::vector<long> row(R, 0);
        for (std::size_t r = 0; r < R; ++r) {
            const double a = num_traits<T>::to_double(A(j, r));
            if (a < 0.0 || std::fabs(a - std::round(a)) > 1e-9)
                throw InputError(
                    "lossn_manjunath: A must contain nonnegative integers -- the residue argument counts "
                    "whole units of capacity. Use lossn_mci, which compares in real arithmetic, or "
                    "lossn_erlangfp");
            row[r] = std::lround(a);
            if (row[r] != 0) any = true;
        }
        // A row of zeros bounds nothing and is dropped, not carried as a
        // one-coefficient dimension: keeping it would make `first`/`last`
        // undefined for that link.
        if (!any) continue;
        long g = std::lround(c);
        for (std::size_t r = 0; r < R; ++r)
            if (row[r] != 0) g = lossn_gcd(g, row[r]);
        if (g > 1) {
            for (std::size_t r = 0; r < R; ++r) row[r] /= g;
            rule.C.push_back(std::lround(c) / g);  // floor, both nonnegative
        } else {
            rule.C.push_back(std::lround(c));
        }
        rule.A.push_back(row);
    }
    return rule;
}

/**
 * Coefficient-domain evaluation of the J-fold contour integral for the
 * right-hand side `C`, which is the full rule for g(C) and the rule shifted by
 * one class requirement for g(C - A e_r).
 *
 * The series lives in a flat vector indexed column-major over the live links,
 * `stride[k] = prod_{i<k} curdim[i]`, with `curdim[j] == 1` while link j is not
 * live and `C[j] + 1` while it is. `f[r]` holds the (possibly scaled) terms
 * nu_r^n / n! for n = 0..nmaxFull[r].
 */
template <class T>
T lossn_manjunath_series(const std::vector<std::vector<T>>& f, const std::vector<std::vector<long>>& A,
                  const std::vector<long>& C, const std::vector<long>& nmaxFull,
                  const LossnManjunathOptions& options, std::size_t& peak) {
    const std::size_t R = f.size(), J = C.size();
    const T zero = num_traits<T>::from_int(0);

    // The elimination order: link j is created at its first route and summed
    // out after its last, so only an induced width of links is ever live.
    std::vector<std::size_t> first(J, 0), last(J, 0);
    for (std::size_t j = 0; j < J; ++j) {
        bool seen = false;
        for (std::size_t r = 0; r < R; ++r) {
            if (A[j][r] == 0) continue;
            if (!seen) {
                first[j] = r;
                seen = true;
            }
            last[j] = r;
        }
        if (!seen)
            throw NumericError("lossn_manjunath: a constraint row with no nonzero entry reached the "
                               "series; the rule was not reduced");
    }

    std::vector<std::size_t> curdim(J, 1);
    std::vector<T> ser(1, num_traits<T>::from_int(1));

    for (std::size_t r = 0; r < R; ++r) {
        // 1. Create the links whose first route is this one.
        for (std::size_t j = 0; j < J; ++j) {
            if (first[j] != r) continue;
            const std::size_t newdim = static_cast<std::size_t>(C[j]) + 1;
            std::size_t pre = 1, post = 1;
            for (std::size_t k = 0; k < j; ++k) pre *= curdim[k];
            for (std::size_t k = j + 1; k < J; ++k) post *= curdim[k];
            if (newdim != 0 && pre * post > options.max_live_states / newdim)
                throw UnsupportedError(
                    "lossn_manjunath: the exact transform would hold more than " +
                    std::to_string(options.max_live_states) +
                    " series coefficients at once. Peak memory is the product of (C_j+1) over the "
                    "links live at the same time, so a wide constraint row with a large capacity "
                    "is what costs; raise LossnManjunathOptions::max_live_states deliberately, or use "
                    "lossn_mci, which is unbiased at any size");
            std::vector<T> grown(pre * newdim * post, zero);
            // The existing content keeps its coefficients and enters at degree
            // zero in the new variable.
            for (std::size_t q = 0; q < post; ++q)
                for (std::size_t p = 0; p < pre; ++p)
                    grown[p + q * pre * newdim] = ser[p + q * pre];
            ser.swap(grown);
            curdim[j] = newdim;
            if (ser.size() > peak) peak = ser.size();
        }

        // 2. Multiply in route r.
        bool constrained = false;
        for (std::size_t j = 0; j < J; ++j)
            if (A[j][r] != 0) constrained = true;

        if (!constrained) {
            // The route is bounded by no live link, so its z-transform is a
            // constant: sum the whole truncated sequence into the series. For a
            // route absent from every row this is the factor exp(nu_r) below,
            // already folded into `lGfree` by the caller, hence a multiply by 1.
            T s = zero;
            for (long n = 0; n <= nmaxFull[r] && static_cast<std::size_t>(n) < f[r].size(); ++n)
                s += f[r][static_cast<std::size_t>(n)];
            for (T& v : ser) v *= s;
            continue;
        }

        // The degree of route r is capped by every row it appears in, evaluated
        // at THIS right-hand side: the shifted series for g(C - A e_r) admits
        // strictly fewer calls than g(C).
        long nmax = std::min<long>(nmaxFull[r], static_cast<long>(f[r].size()) - 1);
        for (std::size_t j = 0; j < J; ++j)
            if (A[j][r] > 0) nmax = std::min<long>(nmax, C[j] / A[j][r]);

        std::vector<std::size_t> stride(J, 1);
        for (std::size_t k = 1; k < J; ++k) stride[k] = stride[k - 1] * curdim[k - 1];
        const std::size_t P = ser.size();

        std::vector<T> next(P, zero);
        std::vector<std::size_t> sub(J, 0);
        for (long n = 0; n <= nmax; ++n) {
            const T c = f[r][static_cast<std::size_t>(n)];
            if (c == zero) continue;
            if (n == 0) {
                for (std::size_t i = 0; i < P; ++i) next[i] += T(c * ser[i]);
                continue;
            }
            // Shift by n requirement vectors, dropping the coefficients the
            // shift would push past the capacity. Those monomials can never
            // contribute to an extracted coefficient, which is exactly why the
            // truncation is exact rather than an approximation.
            std::fill(sub.begin(), sub.end(), static_cast<std::size_t>(0));
            bool anyok = false;
            for (std::size_t i = 0; i < P; ++i) {
                bool ok = true;
                std::size_t tgt = 0;
                for (std::size_t j = 0; j < J && ok; ++j) {
                    const std::size_t d = sub[j] + static_cast<std::size_t>(A[j][r] * n);
                    if (static_cast<long>(d) > C[j])
                        ok = false;
                    else
                        tgt += d * stride[j];
                }
                if (ok) {
                    next[tgt] += T(c * ser[i]);
                    anyok = true;
                }
                // Odometer over the live grid, in the same column-major order
                // the strides encode.
                for (std::size_t j = 0; j < J; ++j) {
                    if (++sub[j] < curdim[j]) break;
                    sub[j] = 0;
                }
            }
            // The shift only grows with n, so once nothing fits nothing will.
            if (!anyok) break;
        }
        ser.swap(next);

        // 3. Integrate out the links whose last route was this one. The
        // multiplier (z^{C+1}-1)/(z-1) of a '<=' constraint turns the residue
        // into the partial sum of the coefficients of degrees 0..C_j, which is
        // the sum along that dimension.
        for (std::size_t j = 0; j < J; ++j) {
            if (last[j] != r) continue;
            std::size_t pre = 1, post = 1;
            for (std::size_t k = 0; k < j; ++k) pre *= curdim[k];
            for (std::size_t k = j + 1; k < J; ++k) post *= curdim[k];
            const std::size_t dj = curdim[j];
            std::vector<T> summed(pre * post, zero);
            for (std::size_t q = 0; q < post; ++q)
                for (std::size_t d = 0; d < dj; ++d)
                    for (std::size_t p = 0; p < pre; ++p)
                        summed[p + q * pre] += ser[p + d * pre + q * pre * dj];
            ser.swap(summed);
            curdim[j] = 1;
        }
    }

    if (ser.size() != 1)
        throw NumericError("lossn_manjunath: a link was never integrated out; the elimination order is "
                           "inconsistent with the constraint rows");
    return ser[0];
}

}  // namespace detail

/**
 * Exact normalizing constant, carried load and blocking of a loss network.
 *
 * @param nu      offered load of route r, nonnegative (R)
 * @param A       (J x R) nonnegative integer circuit requirements
 * @param C       (J) nonnegative integer capacities
 * @param options the live-coefficient cap
 */
template <class T>
LossnManjunathResult<T> lossn_manjunath(const std::vector<T>& nu, const Matrix<T>& A, const std::vector<T>& C,
                          const LossnManjunathOptions& options = LossnManjunathOptions()) {
    const std::size_t R = nu.size();
    if (A.cols() != R || A.rows() != C.size())
        throw InputError("lossn_manjunath: A must be J x R, matching C and nu");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    for (std::size_t r = 0; r < R; ++r)
        if (nu[r] < zero) throw InputError("lossn_manjunath: nu must be nonnegative");

    LossnManjunathResult<T> out;
    out.QLen.assign(R, zero);
    out.Loss.assign(R, zero);

    const detail::LossnManjunathRule rule = detail::lossn_manjunath_integralize(A, C);
    const std::size_t J = rule.C.size();

    // A route absent from every remaining row never blocks: its marginal is an
    // untruncated Poisson, so it carries its full offered load and factors
    // exp(nu_r) out of g(C).
    std::vector<bool> freeRoute(R, true);
    for (std::size_t j = 0; j < J; ++j)
        for (std::size_t r = 0; r < R; ++r)
            if (rule.A[j][r] != 0) freeRoute[r] = false;

    double lGfree = 0.0;
    bool allFree = true;
    for (std::size_t r = 0; r < R; ++r) {
        if (!freeRoute[r]) {
            allFree = false;
            continue;
        }
        out.QLen[r] = nu[r];
        lGfree += num_traits<T>::to_double(nu[r]);
    }
    if (J == 0 || allFree) {
        out.lG = lGfree;
        return out;
    }

    // Per-route truncation, and the terms f_r(n) = nu_r^n / n!.
    std::vector<long> nmax(R, 0);
    std::vector<std::vector<T>> f(R, std::vector<T>(1, one));
    double logscale = 0.0;
    for (std::size_t r = 0; r < R; ++r) {
        if (freeRoute[r]) continue;
        long v = std::numeric_limits<long>::max();
        for (std::size_t j = 0; j < J; ++j)
            if (rule.A[j][r] > 0) v = std::min<long>(v, rule.C[j] / rule.A[j][r]);
        nmax[r] = v;
        const std::size_t len = static_cast<std::size_t>(v) + 1;

        if (nu[r] == zero) {
            // A route with no offered load contributes only its empty term. The
            // log-space branch below would take log(0), which is why this case
            // is separated rather than clamped to a tiny load.
            f[r].assign(len, zero);
            f[r][0] = one;
            continue;
        }
        if constexpr (num_traits<T>::is_exact) {
            // No overflow to guard against, and scaling would only inflate the
            // denominators; the ratios below are scale free either way.
            f[r].assign(len, zero);
            f[r][0] = one;
            for (std::size_t n = 1; n < len; ++n)
                f[r][n] = T(f[r][n - 1] * nu[r] / num_traits<T>::from_int(static_cast<long>(n)));
        } else {
            using std::exp;
            using std::log;
            // Built in log space, then shifted by its own maximum so the peak
            // term is exactly 1: nu^n/n! peaks at e^nu/sqrt(2 pi nu) and would
            // overflow a double for a load above ~700.
            const T lnu = log(nu[r]);
            std::vector<T> lf(len, zero);
            T lfact = zero, m = zero;
            for (std::size_t n = 0; n < len; ++n) {
                if (n > 0) lfact += log(num_traits<T>::from_int(static_cast<long>(n)));
                lf[n] = T(num_traits<T>::from_int(static_cast<long>(n)) * lnu - lfact);
                if (n == 0 || lf[n] > m) m = lf[n];
            }
            f[r].assign(len, zero);
            for (std::size_t n = 0; n < len; ++n) f[r][n] = exp(T(lf[n] - m));
            logscale += num_traits<T>::to_double(m);
        }
    }

    const T G = detail::lossn_manjunath_series(f, rule.A, rule.C, nmax, options, out.peak_states);
    if (G <= zero)
        throw NumericError(
            "lossn_manjunath: the admissible set is empty -- no state satisfies A n <= C, so the loss "
            "network has no stationary distribution");
    out.lG = num_traits<T>::log_as_double(G) + logscale + lGfree;

    for (std::size_t r = 0; r < R; ++r) {
        if (freeRoute[r]) continue;
        std::vector<long> Cr = rule.C;
        bool overflows = false;
        for (std::size_t j = 0; j < J; ++j) {
            Cr[j] -= rule.A[j][r];
            if (Cr[j] < 0) overflows = true;
        }
        if (overflows) {
            // A single class r call already exceeds a capacity, so the route is
            // blocked in every state, including the empty one.
            out.Loss[r] = one;
            out.QLen[r] = zero;
            continue;
        }
        const T Gr = detail::lossn_manjunath_series(f, rule.A, Cr, nmax, options, out.peak_states);
        // The scale cancels here, which is what lets the terms be normalised.
        const T ratio = T(Gr / G);
        out.QLen[r] = T(nu[r] * ratio);
        out.Loss[r] = T(one - ratio);
    }
    return out;
}

}  // namespace lossn
}  // namespace line

#endif  // LINE_API_LOSSN_LOSSN_MANJUNATH_H
