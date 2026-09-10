/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_HAREL_BOUNDS_H
#define LINE_API_PFQN_HAREL_BOUNDS_H

/**
 * Harel-Namn-Sturm throughput bounds for a single-class closed network.
 *
 * Templated port of jar/src/main/java/jline/api/pfqn/Pfqn_harel_bounds.java.
 * MATLAB has no counterpart, so the JAR is the reference.
 *
 * These are the SHARP bounds of the paper, distinct from the `sb` family
 * already in solver_ba_analyzer: `sb` uses only the first three power sums in
 * closed form, whereas this family evaluates the normalizing constant exactly
 * at small populations and extrapolates from it. Both cite Harel1999; they are
 * different results in it and neither subsumes the other.
 *
 * Write A_i = sum_j rho_j^i for the power sums of the relative utilizations.
 * Then
 *
 *   G(n)  = h_n(rho),  the complete homogeneous symmetric polynomial,
 *   TH(n) = G(n-1) / G(n),               the exact throughput at population n,
 *   LB    = N / (A_1 + (N-1) (A_N/A_1)^{1/(N-1)}),
 *   UB(n) = N / (A_1 + ((N-1)/(n-1)) (n/TH(n) - A_1)),   2 <= n <= N.
 *
 * G(n) IS the normalizing constant of the closed load-independent network at
 * population n, which is the oracle the port is tested against: G(n) computed
 * here must equal pfqn_ca on the same demands.
 *
 * TWO DELIBERATE NOTES ON THE PORT.
 *
 * First, the reference hardcodes G(0)..G(7) as expanded polynomials in the
 * power sums and REFUSES n > 7 with "G(n) polynomial not available". Those
 * expansions are the Newton-Girard recurrence
 *
 *   n G(n) = sum_{i=1..n} A_i G(n-i)
 *
 * unrolled by hand. The port evaluates the recurrence instead: it agrees term
 * for term with the reference at every n <= 7, needs no table, and is
 * EXACT-CAPABLE where the expanded form needs pow. The n <= 7 refusal on the
 * PUBLIC entry points is nevertheless kept, so the contract callers see is the
 * reference's; only the internal ceiling is gone. See _kb/03-api-layer.md.
 *
 * Second, the reference refuses a nonzero think time rather than folding it in,
 * because the bounds are derived for a network with no terminal population.
 * That refusal is reproduced: silently dropping Z would return a bound that
 * does not bound.
 *
 * Arithmetic: G and UB are EXACT-CAPABLE. LB needs an (N-1)-st root and is
 * TRANSCENDENTAL for N > 2; at N <= 2 the root is trivial and LB is exact too.
 */

#include <cstddef>
#include <string>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace pfqn {

/** Return value of pfqn_harel_bounds, mirroring Ret.pfqnHarelBounds. */
template <class T>
struct HarelBoundsResult {
    T LB;                ///< throughput lower bound at population N
    std::vector<T> UB;   ///< UB[n] for n = 2..maxUB; entries 0 and 1 are unset
    std::vector<T> TH;   ///< TH[n] = exact throughput at population n, n = 1..maxUB
    int N = 0;           ///< population the bounds are stated at
    std::size_t k = 0;   ///< number of stations
    int maxUB = 0;       ///< largest n at which an upper bound was formed
};

namespace detail {

/** Power sums A_i = sum_j rho_j^i, i = 1..maxPower; A[0] is unused. */
template <class T>
std::vector<T> harel_power_sums(const std::vector<T>& rho, int maxPower) {
    std::vector<T> A(static_cast<std::size_t>(maxPower) + 1, num_traits<T>::from_int(0));
    for (int i = 1; i <= maxPower; ++i) {
        T s = num_traits<T>::from_int(0);
        for (std::size_t j = 0; j < rho.size(); ++j) s += num_pow_int(rho[j], static_cast<unsigned>(i));
        A[static_cast<std::size_t>(i)] = s;
    }
    return A;
}

/** G(0..n) by the Newton-Girard recurrence n G(n) = sum_i A_i G(n-i). */
template <class T>
std::vector<T> harel_G(const std::vector<T>& A, int n) {
    if (static_cast<int>(A.size()) <= n)
        throw InputError("pfqn_harel_bounds: too few power sums for the requested population");
    std::vector<T> G(static_cast<std::size_t>(n) + 1, num_traits<T>::from_int(0));
    G[0] = num_traits<T>::from_int(1);
    for (int m = 1; m <= n; ++m) {
        T acc = num_traits<T>::from_int(0);
        for (int i = 1; i <= m; ++i)
            acc += T(A[static_cast<std::size_t>(i)] * G[static_cast<std::size_t>(m - i)]);
        G[static_cast<std::size_t>(m)] = T(acc / num_traits<T>::from_int(m));
    }
    return G;
}

/** The reference refuses a nonzero think time rather than folding it in. */
template <class T>
void harel_reject_thinktime(const T& Z, const std::string& who) {
    if (Z != num_traits<T>::from_int(0))
        throw InputError(who +
                         " is only valid for networks with zero think time; the provided think "
                         "time is nonzero");
}

/** Shared input screening of the loading vector. */
template <class T>
void harel_check_rho(const std::vector<T>& rho) {
    if (rho.empty())
        throw InputError("pfqn_harel_bounds: the loading vector must have at least one element");
    for (std::size_t i = 0; i < rho.size(); ++i)
        if (rho[i] <= num_traits<T>::from_int(0))
            throw InputError("pfqn_harel_bounds: all loading factors must be positive");
}

/** LB = N / (A1 + (N-1) (A_N/A_1)^{1/(N-1)}). */
template <class T>
T harel_lower_bound(const std::vector<T>& A, int N) {
    const T A1 = A[1];
    if (N == 1) return T(num_traits<T>::from_int(1) / A1);
    const T ratio = T(A[static_cast<std::size_t>(N)] / A1);
    // At N == 2 the exponent is one, so the root is the ratio itself and the
    // bound stays available in exact arithmetic.
    if (N == 2)
        return T(num_traits<T>::from_int(2) / T(A1 + ratio));
    if constexpr (!num_traits<T>::has_transcendental) {
        throw UnsupportedError(
            "pfqn_harel_lb needs an (N-1)-st root and is unavailable in exact arithmetic for "
            "N > 2");
    } else {
        using std::pow;
        const T root = pow(ratio, T(num_traits<T>::from_int(1) / num_traits<T>::from_int(N - 1)));
        return T(num_traits<T>::from_int(N) / T(A1 + num_traits<T>::from_int(N - 1) * root));
    }
}

/** UB(n) = N / (A1 + ((N-1)/(n-1)) (n/TH(n) - A1)). */
template <class T>
T harel_upper_from_th(const T& A1, int N, int n, const T& THn) {
    if (THn == num_traits<T>::from_int(0))
        throw NumericError("pfqn_harel_bounds: the throughput at the extrapolation point is zero");
    const T nOverTH = T(num_traits<T>::from_int(n) / THn);
    const T den = T(A1 + T(num_traits<T>::from_int(N - 1) / num_traits<T>::from_int(n - 1)) *
                             T(nOverTH - A1));
    if (den == num_traits<T>::from_int(0))
        throw NumericError("pfqn_harel_bounds: the upper-bound denominator vanishes");
    return T(num_traits<T>::from_int(N) / den);
}

}  // namespace detail

/**
 * Lower bound alone.
 *
 * @param rho (k) relative utilizations, all strictly positive
 * @param N   population, at least 1
 * @param Z   think time; must be zero
 */
template <class T>
T pfqn_harel_lb(const std::vector<T>& rho, int N, const T& Z) {
    detail::harel_reject_thinktime(Z, "pfqn_harel_lb");
    if (N < 1) throw InputError("pfqn_harel_lb: the population must be at least 1");
    detail::harel_check_rho(rho);
    const std::vector<T> A = detail::harel_power_sums(rho, N);
    return detail::harel_lower_bound(A, N);
}

/** Zero think time. */
template <class T>
T pfqn_harel_lb(const std::vector<T>& rho, int N) {
    return pfqn_harel_lb(rho, N, num_traits<T>::from_int(0));
}

/**
 * Upper bound extrapolated from the exact throughput at population n.
 *
 * @param n extrapolation point, 2 <= n <= min(N, 7)
 */
template <class T>
T pfqn_harel_ub(const std::vector<T>& rho, int N, int n, const T& Z) {
    detail::harel_reject_thinktime(Z, "pfqn_harel_ub");
    if (N < 1) throw InputError("pfqn_harel_ub: the population must be at least 1");
    if (n < 2) throw InputError("pfqn_harel_ub: the extrapolation point must be at least 2");
    if (n > N) throw InputError("pfqn_harel_ub: the extrapolation point cannot exceed N");
    // Kept from the reference, whose hardcoded G(n) table stops at 7.
    if (n > 7) throw InputError("pfqn_harel_ub: the extrapolation point cannot exceed 7");
    detail::harel_check_rho(rho);
    const std::vector<T> A = detail::harel_power_sums(rho, n);
    const std::vector<T> G = detail::harel_G(A, n);
    if (G[static_cast<std::size_t>(n)] == num_traits<T>::from_int(0))
        throw NumericError("pfqn_harel_ub: the normalizing constant vanishes");
    const T THn = T(G[static_cast<std::size_t>(n) - 1] / G[static_cast<std::size_t>(n)]);
    return detail::harel_upper_from_th(A[1], N, n, THn);
}

/** Zero think time. */
template <class T>
T pfqn_harel_ub(const std::vector<T>& rho, int N, int n) {
    return pfqn_harel_ub(rho, N, n, num_traits<T>::from_int(0));
}

/**
 * Both bounds, plus the exact throughputs the upper bounds extrapolate from.
 *
 * @param maxUB largest extrapolation point; defaults to min(N, 7) when <= 0
 */
template <class T>
HarelBoundsResult<T> pfqn_harel_bounds(const std::vector<T>& rho, int N, const T& Z, int maxUB) {
    detail::harel_reject_thinktime(Z, "pfqn_harel_bounds");
    if (N < 1) throw InputError("pfqn_harel_bounds: the population must be at least 1");
    detail::harel_check_rho(rho);
    const int effectiveMaxUB = maxUB > 0 ? maxUB : (N < 7 ? N : 7);
    if (effectiveMaxUB > 7)
        throw InputError("pfqn_harel_bounds: upper bounds are available only for n <= 7");
    if (effectiveMaxUB > N)
        throw InputError("pfqn_harel_bounds: the extrapolation point cannot exceed N");

    HarelBoundsResult<T> res;
    res.N = N;
    res.k = rho.size();
    res.maxUB = effectiveMaxUB;

    // The lower bound reads A up to N, the upper bounds only up to maxUB.
    const int maxPower = N > effectiveMaxUB ? N : effectiveMaxUB;
    const std::vector<T> A = detail::harel_power_sums(rho, maxPower);
    res.LB = detail::harel_lower_bound(A, N);

    const std::vector<T> G = detail::harel_G(A, effectiveMaxUB);
    const T zero = num_traits<T>::from_int(0);
    res.TH.assign(static_cast<std::size_t>(effectiveMaxUB) + 1, zero);
    res.UB.assign(static_cast<std::size_t>(effectiveMaxUB) + 1, zero);
    for (int n = 1; n <= effectiveMaxUB; ++n) {
        if (G[static_cast<std::size_t>(n)] == zero)
            throw NumericError("pfqn_harel_bounds: the normalizing constant vanishes");
        res.TH[static_cast<std::size_t>(n)] =
            T(G[static_cast<std::size_t>(n) - 1] / G[static_cast<std::size_t>(n)]);
    }
    for (int n = 2; n <= effectiveMaxUB; ++n)
        res.UB[static_cast<std::size_t>(n)] =
            detail::harel_upper_from_th(A[1], N, n, res.TH[static_cast<std::size_t>(n)]);
    return res;
}

/** Zero think time, default extrapolation ceiling min(N, 7). */
template <class T>
HarelBoundsResult<T> pfqn_harel_bounds(const std::vector<T>& rho, int N) {
    return pfqn_harel_bounds(rho, N, num_traits<T>::from_int(0), 0);
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_HAREL_BOUNDS_H
