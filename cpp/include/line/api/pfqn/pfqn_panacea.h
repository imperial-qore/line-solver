/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_PFQN_PANACEA_H
#define LINE_API_PFQN_PFQN_PANACEA_H

/**
 * PANACEA normal-usage asymptotic expansion of the normalizing constant
 * (Ramakrishnan and Mitra, BSTJ 61(10):2849-2872, 1982).
 *
 * Templated port of matlab/src/api/pfqn/pfqn_panacea.m. In the normal usage
 * regime, where alpha_r = 1 - sum_i N_i r_ir > 0 for every station, the
 * constant admits the expansion
 *
 *   log G = -sum_r factln(N_r) + sum_r N_r log Z_r + log(sum_k I_k) - sum_i log alpha_i
 *
 * whose coefficients I_2, I_3 are assembled from convolution-algorithm values
 * of the scaled demand matrix gammatilde at small auxiliary populations. The
 * expansion is available at 1, 2 or 3 terms, as in the original package.
 *
 * NOT-NORMAL-USAGE. MATLAB returns NaN when min(alpha) < 0. The port reports
 * it through a flag on the result rather than a NaN, so a caller that ignores
 * the flag gets a value it can recognize as unusable instead of a quiet NaN
 * travelling into a solver.
 *
 * ARITHMETIC. The result is the logarithm of a truncated asymptotic series, so
 * the routine is gated on num_traits<T>::has_transcendental. The pfqn_ca calls
 * inside it are exact and would remain so at exact arithmetic; it is the
 * expansion and the logarithms that are not.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_asympt_common.h"
#include "line/api/pfqn/pfqn_ca.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** Return value of pfqn_panacea, mirroring [Gn, lGn]. */
template <class T>
struct PanaceaResult {
    T G;
    T lG;
    bool normalUsage;  ///< false where MATLAB returns NaN (min alpha < 0)
};

/**
 * @param L     (M x R) demands
 * @param N     (R) population
 * @param Z     (R) think times; empty means MATLAB's 1e-8 placeholder
 * @param terms 1, 2 or 3
 */
template <class T>
PanaceaResult<T> pfqn_panacea(const Matrix<T>& L, const std::vector<int>& N,
                              const std::vector<T>& Z, int terms) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_panacea requires transcendental arithmetic (asymptotic expansion of log G)");
    using std::exp;
    using std::log;
    const std::size_t q = L.rows(), p = L.cols();
    if (N.size() != p) throw InputError("pfqn_panacea: L and N disagree on the class count");
    if (terms < 1 || terms > 3)
        throw InputError("pfqn_panacea: the terms parameter must be 1, 2 or 3");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    std::vector<T> Zv = Z;
    if (Zv.empty()) Zv.assign(p, num_traits<T>::from_double(1e-8));
    if (Zv.size() != p) throw InputError("pfqn_panacea: Z has the wrong length");

    PanaceaResult<T> res;
    res.normalUsage = true;

    // Degenerate: no station carries any demand, the delay carries everything.
    bool anyDemand = false;
    for (std::size_t i = 0; i < q; ++i)
        for (std::size_t r = 0; r < p; ++r)
            if (L(i, r) != zero) anyDemand = true;
    if (q == 0 || !anyDemand) {
        T lG = zero;
        for (std::size_t r = 0; r < p; ++r) {
            lG -= detail::num_factln<T>(num_traits<T>::from_int(N[r]));
            lG += num_traits<T>::from_int(N[r]) * log(Zv[r]);
        }
        res.lG = lG;
        res.G = exp(lG);
        return res;
    }

    // scaled-load rationale: see _kb/03-api-layer.md (cpp port notes: pfqn)
    Matrix<T> rr(q, p, zero);
    T invmax = zero;
    bool first = true;
    for (std::size_t i = 0; i < q; ++i)
        for (std::size_t s = 0; s < p; ++s) {
            if (Zv[s] == zero) throw InputError("pfqn_panacea: zero think time");
            rr(i, s) = T(L(i, s) / Zv[s]);
            if (rr(i, s) > zero) {
                const T inv = T(one / rr(i, s));
                if (first || inv > invmax) {
                    invmax = inv;
                    first = false;
                }
            }
        }
    if (first) throw InputError("pfqn_panacea: no positive demand ratio");
    const T Nt = invmax;

    std::vector<T> beta(p);
    for (std::size_t s = 0; s < p; ++s) beta[s] = T(num_traits<T>::from_int(N[s]) / Nt);
    Matrix<T> gamma(q, p);
    for (std::size_t i = 0; i < q; ++i)
        for (std::size_t s = 0; s < p; ++s) gamma(i, s) = T(rr(i, s) * Nt);
    std::vector<T> alpha(q);
    for (std::size_t i = 0; i < q; ++i) {
        T a = one;
        for (std::size_t s = 0; s < p; ++s) a -= num_traits<T>::from_int(N[s]) * rr(i, s);
        alpha[i] = a;
        if (a < zero) res.normalUsage = false;
    }
    if (!res.normalUsage) {
        res.G = zero;
        res.lG = zero;
        return res;
    }
    Matrix<T> gt(q, p);
    for (std::size_t i = 0; i < q; ++i)
        for (std::size_t s = 0; s < p; ++s) gt(i, s) = T(gamma(i, s) / alpha[i]);

    const Matrix<T> noZ;
    std::vector<T> I;
    I.push_back(one);  // A0
    if (terms >= 2) {
        T A1 = zero;
        for (std::size_t j = 0; j < p; ++j) {
            std::vector<int> m(p, 0);
            m[j] = 2;
            A1 -= beta[j] * pfqn_ca(gt, m, noZ).G;
        }
        I.push_back(T(A1 / Nt));
    }
    if (terms >= 3) {
        T A2 = zero;
        for (std::size_t j = 0; j < p; ++j) {
            std::vector<int> m(p, 0);
            m[j] = 3;
            A2 += num_traits<T>::from_int(2) * beta[j] * pfqn_ca(gt, m, noZ).G;
            m.assign(p, 0);
            m[j] = 4;
            A2 += num_traits<T>::from_int(3) * T(beta[j] * beta[j]) * pfqn_ca(gt, m, noZ).G;
            for (std::size_t k = 0; k < p; ++k) {
                if (k == j) continue;
                m.assign(p, 0);
                m[j] = 2;
                m[k] = 2;
                A2 += num_traits<T>::from_rational(1, 2) * beta[j] * beta[k] * pfqn_ca(gt, m, noZ).G;
            }
        }
        I.push_back(T(A2 / T(Nt * Nt)));
    }

    T Isum = zero;
    for (const T& v : I) Isum += v;
    if (Isum <= zero) throw NumericError("pfqn_panacea: non-positive asymptotic series");
    T lG = zero;
    for (std::size_t s = 0; s < p; ++s) {
        lG -= detail::num_factln<T>(num_traits<T>::from_int(N[s]));
        lG += num_traits<T>::from_int(N[s]) * log(Zv[s]);
    }
    lG += log(Isum);
    for (std::size_t i = 0; i < q; ++i) lG -= log(alpha[i]);
    res.lG = lG;
    res.G = exp(lG);
    return res;
}

template <class T>
PanaceaResult<T> pfqn_panacea(const Matrix<T>& L, const std::vector<int>& N,
                              const std::vector<T>& Z) {
    return pfqn_panacea(L, N, Z, 3);
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_PFQN_PANACEA_H
