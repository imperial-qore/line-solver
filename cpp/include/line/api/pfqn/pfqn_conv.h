/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_CONV_H
#define LINE_API_PFQN_CONV_H

/**
 * Multichain convolution algorithm with class-dependent service rates
 * (Sauer 1983, "Computational Algorithms for State-Dependent Queueing
 * Networks", ACM TOCS 1(1):67-92, Section 5.2).
 *
 * Templated port of matlab/src/api/pfqn/pfqn_conv.m.
 *
 * G(N) is the multivariate discrete convolution of the M station factors and
 * the delay factor,
 *
 *   G_0(n) = F_Z(n) = prod_r Z_r^{n_r} / n_r!,
 *   G_m(n) = sum_{0 <= i <= n} X_m(i) G_{m-1}(n - i),
 *
 * where a class-dependent station builds its factor from Sauer eq. (40),
 *
 *   X_m(n) = (|n| / n_r) (L(m,r) / beta_{m,r}(n)) X_m(n - e_r),  X_m(0) = 1,
 *
 * for the first class r with n_r > 0. beta is the DIMENSIONLESS scaling of the
 * service demand supplied by the caller, so beta = 1 means "no correction" and
 * the recurrence collapses to the load-independent multinomial form. A station
 * with no scaling callable is folded in by the classical in-place Buzen update
 *
 *   G_m(n) = G_{m-1}(n) + sum_r L(m,r) G_m(n - e_r)
 *
 * in O(P R) rather than O(P^2), P = prod_r (N_r + 1), so a model with no
 * class-dependent station reduces operation for operation to pfqn_ca on the
 * same demands and returns the identical value.
 *
 * Arithmetic: EXACT-CAPABLE. Every operation is an addition, a multiplication
 * or a division in the field of the inputs; the reference's use of log/exp
 * inside its local Fz is a range-management device for the delay factor and is
 * replaced here by the same detail::pff_delay that pfqn_ca uses. Whether the
 * result is exact for a class-dependent station is a property of the supplied
 * beta callables, which are evaluated but never inspected. Note that unlike
 * pfqn_ca this routine applies NO power-of-two rescaling in double: the
 * reference does not, and a class-dependent station factor is not homogeneous
 * in the demands once beta is state dependent, so no exact exponent shift
 * exists in general. Use T = Real<D> or T = Rational when the constant leaves
 * the double range.
 */

#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_ca.h"
#include "line/api/pfqn/pfqn_cdfun.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"
#include "line/util/population.h"

namespace line {
namespace pfqn {

/**
 * @param L         (M x R) service demands
 * @param N         (R) population per class, finite
 * @param Z         (K x R) think times, summed over rows; may be empty
 * @param cdscaling (M) class-dependence callables; an empty entry marks a
 *                  load-independent station. Pass an empty vector for none.
 */
template <class T>
NcResult<T> pfqn_conv(const Matrix<T>& L, const std::vector<int>& N, const Matrix<T>& Z,
                      const std::vector<CdScaling<T>>& cdscaling) {
    const std::size_t M = L.empty() ? 0 : L.rows();
    const std::size_t R = N.size();
    if (!L.empty() && L.cols() != R)
        throw InputError("pfqn_conv: L and N disagree on the class count");
    if (!cdscaling.empty() && cdscaling.size() != M)
        throw InputError("pfqn_conv: the scaling vector has the wrong station count");
    for (int v : N)
        if (v < 0) throw InputError("pfqn_conv: the convolution algorithm requires finite, "
                                    "nonnegative (closed) populations");

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    std::vector<T> Zsum(R, zero);
    if (!Z.empty()) {
        if (Z.cols() != R) throw InputError("pfqn_conv: Z and N disagree on the class count");
        for (std::size_t k = 0; k < Z.rows(); ++k)
            for (std::size_t r = 0; r < R; ++r) Zsum[r] += Z(k, r);
    }

    const std::vector<std::size_t> prods = plane_sizes(N);
    const std::size_t total = population_count(N);

    // ---- G_0(n) = F_Z(n) ----------------------------------------------------
    std::vector<T> G(total, zero);
    {
        std::vector<int> n(R, 0);
        bool more = true;
        while (more) {
            G[pop_index(n, prods)] = detail::pff_delay(Zsum, n);
            more = next_pop(n, N);
        }
    }

    // ---- fold in one station at a time --------------------------------------
    for (std::size_t ist = 0; ist < M; ++ist) {
        const bool isCd = !cdscaling.empty() && static_cast<bool>(cdscaling[ist]);
        if (!isCd) {
            // Load-independent: in-place Buzen update, lexicographic order
            // guarantees n - e_r is already at its new value.
            std::vector<int> n(R, 0);
            bool more = true;
            while (more) {
                const std::size_t idx = pop_index(n, prods);
                T acc = G[idx];
                for (std::size_t r = 0; r < R; ++r)
                    if (n[r] >= 1) acc += L(ist, r) * G[idx - prods[r]];
                G[idx] = acc;
                more = next_pop(n, N);
            }
            continue;
        }

        // Class-dependent: build the station factor X_m, then convolve.
        std::vector<T> Xm(total, zero);
        Xm[0] = one;
        {
            std::vector<int> n(R, 0);
            bool more = next_pop(n, N);  // X_m(0) is already set
            while (more) {
                const std::size_t idx = pop_index(n, prods);
                std::size_t r = 0;
                while (r < R && n[r] == 0) ++r;
                // next_pop only yields nonzero vectors past the origin, so r < R.
                std::vector<T> row(R);
                for (std::size_t s = 0; s < R; ++s)
                    row[s] = num_traits<T>::from_int(n[s]);
                const std::vector<T> bval = cdscaling[ist](row);
                if (bval.empty())
                    throw InputError("pfqn_conv: a class-dependence callable returned nothing");
                const T beta = bval.size() > 1 ? bval.at(r) : bval[0];
                if (beta > zero) {
                    int tot = 0;
                    for (int v : n) tot += v;
                    const T fac = num_traits<T>::from_int(tot) / num_traits<T>::from_int(n[r]);
                    Xm[idx] = fac * (L(ist, r) / beta) * Xm[idx - prods[r]];
                } else {
                    Xm[idx] = zero;  // a nonpositive scaling kills the state
                }
                more = next_pop(n, N);
            }
        }

        std::vector<T> Gold(G);
        std::vector<int> n(R, 0);
        bool more = true;
        while (more) {
            const std::size_t idxn = pop_index(n, prods);
            T acc = zero;
            // Inner sweep over 0 <= i <= n.
            std::vector<int> i(R, 0);
            bool more_i = true;
            while (more_i) {
                std::size_t idx_i = 0, idx_nmi = 0;
                for (std::size_t r = 0; r < R; ++r) {
                    idx_i += prods[r] * static_cast<std::size_t>(i[r]);
                    idx_nmi += prods[r] * static_cast<std::size_t>(n[r] - i[r]);
                }
                acc += Xm[idx_i] * Gold[idx_nmi];
                more_i = next_pop(i, n);
            }
            G[idxn] = acc;
            more = next_pop(n, N);
        }
    }

    const T Gn = G[total - 1];
    return {Gn, num_traits<T>::log_as_double(Gn)};
}

/** Overload with no class dependence, i.e. plain multichain convolution. */
template <class T>
NcResult<T> pfqn_conv(const Matrix<T>& L, const std::vector<int>& N, const Matrix<T>& Z) {
    return pfqn_conv(L, N, Z, std::vector<CdScaling<T>>());
}

template <class T>
NcResult<T> pfqn_conv(const Matrix<T>& L, const std::vector<int>& N) {
    return pfqn_conv(L, N, Matrix<T>(), std::vector<CdScaling<T>>());
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_CONV_H
