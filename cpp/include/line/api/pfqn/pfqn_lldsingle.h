/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_LLDSINGLE_H
#define LINE_API_PFQN_LLDSINGLE_H

/**
 * Exact normalizing constant of a SINGLE-CLASS closed network whose stations
 * are LIMITED load dependent, i.e. whose rate functions stay constant past a
 * per-station threshold.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_lldsingle.m.
 *
 * Same recursion, same arithmetic and bit-identical results to
 * pfqn_gldsingle, but with the rate-offset axis truncated at that threshold
 * instead of at the population. Unrolling the recursion of pfqn_gldsingle,
 *
 *   g(m, n, t) = g(m-1, n, 1) + L(m) g(m, n-1, t+1) / mu(m, t)
 *
 * shows that its third index is an offset into station m's rate function,
 *
 *   g(m,n,t) = sum_{j=0..n} prod_{i=0..j-1} L(m)/mu(m,t+i) * g(m-1,n-j,1)
 *
 * so once t >= s_m, where s_m is the population past which mu(m, .) stays
 * constant, every factor is mu(m, s_m), the product collapses to
 * (L(m)/mu(m,s_m))^j and
 *
 *   g(m, n, t) = g(m, n, s_m)   for all t >= s_m
 *
 * The N - s_m upper slices that pfqn_gldsingle materializes are duplicates of
 * one another. Capping the offset at s_m and reading g(m, n-1, min(t+1, s_m))
 * keeps every value the answer reads.
 *
 * COST. O(N sum_k s_k) time against O(M N^2) for pfqn_gldsingle, and
 * O(N max_k s_k) space against O(M N^2), the station levels being rolled. On a
 * multiserver model, where s_k is the server count, this is LINEAR in the
 * population rather than quadratic. Unlike pfqn_explicit_ld, which reaches the
 * same asymptotics through Gordon's alternating partial fraction, this loses no
 * digits to cancellation at T = double: the arithmetic performed is a SUBSET of
 * pfqn_gldsingle's, so the two agree to the last bit in every field.
 *
 * There is no gain on a station whose rates never settle, an infinite server
 * mu(m,n) = n being the usual case: it gets s_m = N and costs what it costs in
 * pfqn_gldsingle. The saving is over the OTHER stations, so a model carrying
 * one delay among M queues drops from O(M N^2) to O(N^2 + N sum_k s_k).
 *
 * Arithmetic: EXACT-CAPABLE, and the threshold scan is what keeps it so. It
 * compares rates with the field's own operator==, never a tolerance: at
 * T = Rational or T = Real<D> a tolerance has no meaning, and at T = double a
 * multiserver row repeats its tail exactly. A MISSED tie only costs time, since
 * the routine then behaves as pfqn_gldsingle; a FALSE tie would be a wrong
 * answer, which exact comparison cannot produce. The reference carries an
 * eps-relative tolerance instead, being confined to IEEE double.
 *
 * As in pfqn_gldsingle this port keeps only the linear recursion, the
 * reference's log-space branch being a range-management device for IEEE double,
 * and an infinite rate is accepted the same way: the term L/mu vanishes.
 */

#include <algorithm>
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
NcResult<T> pfqn_lldsingle(const Matrix<T>& L, int N, const Matrix<T>& mu) {
    if (!L.empty() && L.cols() != 1)
        throw InputError("pfqn_lldsingle: multiclass model detected, this routine is single class");
    if (N < 0) throw InputError("pfqn_lldsingle: negative population");

    const std::size_t M = L.empty() ? 0 : L.rows();
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    if (N == 0) return {one, 0.0};
    if (M == 0) return {zero, num_traits<T>::log_as_double(zero)};
    if (mu.rows() != M) throw InputError("pfqn_lldsingle: mu has the wrong station count");
    if (static_cast<int>(mu.cols()) < N)
        throw InputError("pfqn_lldsingle: mu has fewer rate columns than the population");

    const std::size_t Nu = static_cast<std::size_t>(N);

    // pfqn_gldsingle reads every rate column 1..N somewhere in its triangle and
    // rejects a zero there; the capped sweep reads a subset, so the rejection is
    // hoisted to keep the two routines refusing the same models.
    for (std::size_t m = 0; m < M; ++m)
        for (std::size_t t = 0; t < Nu; ++t)
            if (mu(m, t) == zero)
                throw NumericError(
                    "pfqn_lldsingle: a load-dependent rate is zero, the station cannot serve "
                    "and the normalizing constant diverges");

    // s[m]: the smallest offset past which row m of mu is constant, so that
    // mu(m,t) == mu(m,s[m]) for every t >= s[m]. Stored one-based.
    std::vector<std::size_t> s(M, Nu);
    for (std::size_t m = 0; m < M; ++m) {
        const T& tail = mu(m, Nu - 1);
        for (std::size_t n = Nu - 1; n >= 1; --n) {
            if (mu(m, n - 1) == tail)
                s[m] = n;
            else
                break;
        }
    }

    // gprev[n] = g(m-1, n, 1); cur[n * sm + (t-1)] = g(m, n, t)
    std::vector<T> gprev(Nu + 1, zero);
    gprev[0] = one;  // g(0, 0, 1) = 1, and g(0, n, 1) = 0 for n >= 1
    std::vector<T> cur;
    for (std::size_t m = 1; m <= M; ++m) {
        const std::size_t sm = s[m - 1];
        cur.assign((Nu + 1) * sm, zero);
        for (std::size_t t = 0; t < sm; ++t) cur[t] = one;  // g(m, 0, t) = 1
        for (std::size_t n = 1; n <= Nu; ++n) {
            // offsets above N-n+1 are never read back, exactly as in
            // pfqn_gldsingle, so the triangle is kept
            const std::size_t tmax = std::min(sm, Nu - n + 1);
            for (std::size_t t = 1; t <= tmax; ++t) {
                const std::size_t tsrc = std::min(t + 1, sm);
                cur[n * sm + (t - 1)] =
                    gprev[n] + L(m - 1, 0) * cur[(n - 1) * sm + (tsrc - 1)] / mu(m - 1, t - 1);
            }
        }
        for (std::size_t n = 0; n <= Nu; ++n) gprev[n] = cur[n * sm];
    }

    const T G = gprev[Nu];
    return {G, num_traits<T>::log_as_double(G)};
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_LLDSINGLE_H
