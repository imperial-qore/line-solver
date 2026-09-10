/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_PFQN_LDBCMP_H
#define LINE_API_PFQN_PFQN_LDBCMP_H

/**
 * Anselmi-Cremonesi (2008) lower throughput bound for a closed single-class
 * BCMP network with load-dependent stations.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_ldbcmp.m. The bound (their eq. 15)
 * uses the fact that a closed BCMP network is, as N -> inf, equivalent to the
 * open network obtained by removing the bottleneck and injecting at rate
 * 1/Dmax; Algorithm 1 refines it to a monotone fixed point
 *
 *   X = (N - Qhat) / [ Dmax (b + N - Qhat) - b (Dmax X')^N Dmax ]
 *
 * where Qhat is the sum of the non-bottleneck open queue lengths, b the number
 * of bottleneck stations, and c(i) the Heffes load-dependence coefficient.
 *
 * APPLICABILITY. The bound requires N >= Qhat and every non-bottleneck
 * utilization below one. MATLAB returns NaN in both cases; the port reports it
 * through a flag on the result instead, since a NaN throughput propagates
 * silently while a flag has to be read.
 *
 * ARITHMETIC. The fixed point raises Dmax X to the integer power N, which is
 * num_pow_int and stays in the field, and everything else is an addition or a
 * division. The bound is therefore EXACT in rational arithmetic and is left
 * ungated -- but note that the iteration is a contraction, not a closed form,
 * so what is exact is each iterate, not the limit.
 */

#include <cstddef>
#include <limits>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace pfqn {

/** Return value of pfqn_ldbcmp, mirroring [Xlo, Rhi, Qhat]. */
template <class T>
struct LdBcmpBound {
    T Xlo;
    T Rhi;
    T Qhat;
    bool applicable;  ///< false where MATLAB returns NaN (N < Qhat, or rho >= 1)
};

/**
 * @param L (M) limiting demands, @param N population, @param Z think time
 * @param c (M) Heffes coefficients, empty for all fixed-rate stations
 * @param tol relative fixed-point tolerance (MATLAB default 1e-10)
 */
template <class T>
LdBcmpBound<T> pfqn_ldbcmp(const std::vector<T>& L, const T& N, const T& Z, const std::vector<T>& c,
                           const T& tol) {
    const std::size_t M = L.size();
    if (M == 0) throw InputError("pfqn_ldbcmp: empty demand vector");
    if (!c.empty() && c.size() != M) throw InputError("pfqn_ldbcmp: c has the wrong length");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);

    T Dm = L[0];
    for (const T& x : L)
        if (x > Dm) Dm = x;
    if (Dm <= zero) throw InputError("pfqn_ldbcmp: all demands are zero");
    // MATLAB's bottleneck test is abs(D - Dm) <= 1e-12*Dm; the tolerance is a
    // floating-point guard and is kept, so the two agree station for station.
    const T bt = T(num_traits<T>::from_double(1e-12) * Dm);
    long bmax = 0;
    std::vector<bool> isbott(M, false);
    for (std::size_t i = 0; i < M; ++i) {
        if (num_abs(T(L[i] - Dm)) <= bt) {
            isbott[i] = true;
            ++bmax;
        }
    }

    LdBcmpBound<T> r;
    r.applicable = true;
    const T lambda = T(one / Dm);
    T Qhat = zero;
    for (std::size_t i = 0; i < M; ++i) {
        if (isbott[i]) continue;
        const T rho = T(lambda * L[i]);
        if (rho >= one) {
            r.applicable = false;
            r.Xlo = zero;
            r.Rhi = zero;
            r.Qhat = zero;
            return r;
        }
        const T ci = c.empty() ? zero : c[i];
        Qhat += T(T(ci + one) * rho / T(one - rho));
    }
    Qhat += lambda * Z;
    r.Qhat = Qhat;
    if (T(N - Qhat) < zero) {
        r.applicable = false;
        r.Xlo = zero;
        r.Rhi = zero;
        return r;
    }

    const T a = T(N - Qhat);
    // Dmax X' is raised to the integer power N, so N must be integral here,
    // exactly as in MATLAB, where N indexes a closed population.
    const double Nd = num_traits<T>::to_double(N);
    if (Nd < 0.0 || Nd != std::floor(Nd))
        throw InputError("pfqn_ldbcmp: the population must be a non-negative integer");
    const unsigned Nu = static_cast<unsigned>(Nd);

    T Xprime = zero, Xlo = zero;
    for (int it = 0; it < 10000; ++it) {
        const T Xprev = Xlo;
        const T denom = T(Dm * T(num_traits<T>::from_int(bmax) + N - Qhat) -
                          num_traits<T>::from_int(bmax) * num_pow_int(T(Dm * Xprime), Nu) * Dm);
        if (denom == zero) throw NumericError("pfqn_ldbcmp: zero denominator in the fixed point");
        Xlo = T(a / denom);
        Xprime = Xlo;
        if (Xprev > zero && T(num_abs(T(Xprev - Xlo)) / Xprev) <= tol) break;
    }
    r.Xlo = Xlo;
    if (Xlo == zero) {
        // Exactly at the regime boundary N == Qhat the numerator N - Qhat is
        // zero, so the bound degenerates to X >= 0: still valid, and MATLAB
        // returns it with Rhi = N/0 = Inf rather than erroring. An exact
        // backend has no infinity to return there.
        if constexpr (num_traits<T>::has_transcendental) {
            r.Rhi = num_traits<T>::from_double(std::numeric_limits<double>::infinity());
            return r;
        } else {
            throw UnsupportedError(
                "pfqn_ldbcmp: at N == Qhat the throughput bound degenerates to zero and the "
                "response-time bound is infinite, which exact arithmetic cannot represent");
        }
    }
    r.Rhi = T(N / Xlo);
    return r;
}

template <class T>
LdBcmpBound<T> pfqn_ldbcmp(const std::vector<T>& L, const T& N, const T& Z) {
    return pfqn_ldbcmp(L, N, Z, std::vector<T>(), num_traits<T>::from_double(1e-10));
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_PFQN_LDBCMP_H
