/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_PFQN_MCUB_H
#define LINE_API_PFQN_PFQN_MCUB_H

/**
 * Kerola's multiclass composite bound (Perf. Eval. 6:1-9, eqs. 10-16) on the
 * per-class throughput of a closed product-form network.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_mcub.m. The multiclass Balanced
 * Job Bound gives the per-class lower bound
 *
 *   X_r^- = N_r / (sum_k L_kr + Z_r + (Ntot - 1) max_k L_kr)
 *
 * and the residual-utilization argument then gives the composite upper bound
 *
 *   X_r^+ = min_k [1 - sum_{s != r} X_s^- L_ks] / L_kr
 *
 * at O(MR). Despite the name it has nothing to do with pfqn_cub, which is the
 * Grundmann-Moeller cubature normalizing constant.
 *
 * ARITHMETIC. Sums, products, maxima and divisions only, so the bound is EXACT
 * in rational arithmetic and is deliberately left ungated.
 */

#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** Return value of pfqn_mcub, mirroring [Xub, Xlb]. */
template <class T>
struct McubBounds {
    std::vector<T> Xub;
    std::vector<T> Xlb;
};

/**
 * @param L (M x R) demands, @param N (R) population, @param Z (R) think times
 *          (empty for zero think time)
 */
template <class T>
McubBounds<T> pfqn_mcub(const Matrix<T>& L, const std::vector<T>& N, const std::vector<T>& Z) {
    const std::size_t M = L.rows(), R = L.cols();
    if (N.size() != R) throw InputError("pfqn_mcub: L and N disagree on the class count");
    if (!Z.empty() && Z.size() != R) throw InputError("pfqn_mcub: Z has the wrong length");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    T Ntot = zero;
    for (const T& v : N) Ntot += v;

    McubBounds<T> r;
    r.Xlb.assign(R, zero);
    r.Xub.assign(R, zero);
    for (std::size_t s = 0; s < R; ++s) {
        T R0 = zero, Lb = zero;
        for (std::size_t k = 0; k < M; ++k) {
            R0 += L(k, s);
            if (L(k, s) > Lb) Lb = L(k, s);
        }
        const T zr = Z.empty() ? zero : Z[s];
        const T den = T(R0 + zr + T(Ntot - one) * Lb);
        if (den == zero) throw NumericError("pfqn_mcub: zero cycle time in the balanced job bound");
        r.Xlb[s] = T(N[s] / den);
    }

    for (std::size_t s = 0; s < R; ++s) {
        bool any = false;
        T best = zero;
        for (std::size_t k = 0; k < M; ++k) {
            if (!(L(k, s) > zero)) continue;  // MATLAB leaves the device at Inf
            T Uoth = zero;
            for (std::size_t t = 0; t < R; ++t)
                if (t != s) Uoth += r.Xlb[t] * L(k, t);
            const T dev = T(T(one - Uoth) / L(k, s));
            if (!any || dev < best) {
                best = dev;
                any = true;
            }
        }
        // MATLAB's min over an all-Inf device vector is Inf; a class with no
        // demand anywhere is unconstrained, which the caller must handle.
        if (!any) throw InputError("pfqn_mcub: a class has zero demand at every station");
        r.Xub[s] = best;
    }
    return r;
}

template <class T>
McubBounds<T> pfqn_mcub(const Matrix<T>& L, const std::vector<T>& N) {
    return pfqn_mcub(L, N, std::vector<T>());
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_PFQN_MCUB_H
