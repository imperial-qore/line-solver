/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_GLDSINGLE_H
#define LINE_API_PFQN_GLDSINGLE_H

/**
 * Exact normalizing constant of a SINGLE-CLASS closed network whose stations
 * are load dependent.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_gldsingle.m.
 *
 * The recursion carries a rate offset t alongside the station index m and the
 * population n, so that the jobs already placed at station m shift its rate
 * lattice without materializing a separate shifted matrix:
 *
 *   g(0, n, t)   = 0            for n >= 1
 *   g(m, 0, t)   = 1
 *   g(m, n, t)   = g(m-1, n, 1) + L(m) g(m, n-1, t+1) / mu(m, t)
 *   G            = g(M, N, 1)
 *
 * This is the single-class specialization of pfqn_gld; the two agree to the
 * last bit on every model both accept, and the specialization is kept because
 * pfqn_ncld dispatches to it directly for R = 1 and because it costs O(M N^2)
 * rather than going through the general convolution.
 *
 * Arithmetic: EXACT-CAPABLE. MATLAB carries TWO implementations of the same
 * recursion, a linear one and a log-space one selected by
 *
 *   useLog = isreal(L) && isreal(mu) && all(L>=0) && all(mu>0)
 *
 * with a pairwise log-sum-exp replacing the addition. That branch is purely a
 * range-management device for IEEE double -- the comment in the reference says
 * so explicitly, citing underflow of the delay term to realmin at N >= 190 --
 * and the two branches compute the same mathematical quantity. This port keeps
 * only the linear recursion, which is exact in any field: at T = Rational
 * there is no underflow to manage, and at T = Real<D> the exponent range is
 * wide enough that the models which drove the reference into log space stay in
 * range. Callers that genuinely need the double path on such a model should
 * raise the arithmetic rather than reintroduce the logs, since the log-sum-exp
 * form cannot represent the negative intermediate rates that the reference's
 * own `useLog` guard exists to fall back from.
 *
 * An infinite rate is accepted the same way the reference accepts it: the term
 * L/mu vanishes. In an exact field there is no infinity, so a caller expressing
 * "this station cannot hold this many jobs" must pass a zero DEMAND instead.
 */

#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_ca.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/**
 * @param L  (M x 1) service demands, one class
 * @param N  population
 * @param mu (M x >=N) load-dependent rates, mu(i,k) with k jobs at station i
 */
template <class T>
NcResult<T> pfqn_gldsingle(const Matrix<T>& L, int N, const Matrix<T>& mu) {
    if (!L.empty() && L.cols() != 1)
        throw InputError("pfqn_gldsingle: multiclass model detected, this routine is single class");
    if (N < 0) throw InputError("pfqn_gldsingle: negative population");

    const std::size_t M = L.empty() ? 0 : L.rows();
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    if (N == 0) return {one, 0.0};
    if (M == 0) return {zero, num_traits<T>::log_as_double(zero)};
    if (mu.rows() != M) throw InputError("pfqn_gldsingle: mu has the wrong station count");
    if (static_cast<int>(mu.cols()) < N)
        throw InputError("pfqn_gldsingle: mu has fewer rate columns than the population");

    const std::size_t Nu = static_cast<std::size_t>(N);
    // g[m][n][t], t = 1 .. N+2 stored at index t; index 0 unused. The station
    // index m runs 0 .. M with m = 0 the empty network.
    const std::size_t Tdim = Nu + 3;
    std::vector<T> g(static_cast<std::size_t>(M + 1) * (Nu + 1) * Tdim, zero);
    const auto at = [&](std::size_t m, std::size_t n, std::size_t t) -> T& {
        return g[(m * (Nu + 1) + n) * Tdim + t];
    };

    // g(0, n, t) = 0 for n >= 1 (already zero); g(0, 0, t) is never read.
    for (std::size_t m = 1; m <= M; ++m) {
        for (std::size_t t = 1; t <= Nu + 1; ++t) at(m, 0, t) = one;
        for (std::size_t n = 1; n <= Nu; ++n) {
            for (std::size_t t = 1; t + n <= Nu + 1; ++t) {
                const T& rate = mu(m - 1, t - 1);
                if (rate == zero)
                    throw NumericError(
                        "pfqn_gldsingle: a load-dependent rate is zero, the station cannot serve "
                        "and the normalizing constant diverges");
                at(m, n, t) = at(m - 1, n, 1) + L(m - 1, 0) * at(m, n - 1, t + 1) / rate;
            }
        }
    }

    const T G = at(M, Nu, 1);
    return {G, num_traits<T>::log_as_double(G)};
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_GLDSINGLE_H
