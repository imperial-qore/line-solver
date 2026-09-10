/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_PFQN_SCB_H
#define LINE_API_PFQN_PFQN_SCB_H

/**
 * Dowdy-Carlson-Krantz-Tripathi (1992) single-class bounds of multi-class
 * queueing networks, J. ACM 39(1):188-213.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_scb.m, pfqn_scbgap.m,
 * pfqn_usumbound.m and pfqn_minclasses.m.
 *
 * SEMANTICS DIFFER FROM EVERY OTHER pfqn_* BOUND IN THIS TREE. aba/bjb/gb/...
 * bracket the exact solution OF THE GIVEN MODEL; pfqn_scb brackets the
 * multiclass system that the given single-class model aggregates. Its lower
 * side is therefore the EXACT single-class solution, not an approximation of
 * it, and mixing the family into an auto composite would compare two different
 * quantities.
 *
 * ARITHMETIC. pfqn_scb runs the exact single-class MVA recursion (additions,
 * multiplications and one division per population step) and then scales by a
 * rational factor, so it stays in the field and is left ungated. The three
 * combinatorial bounds are pure rational expressions in N, K and r.
 */

#include <algorithm>
#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace pfqn {

/** Return value of pfqn_scb, mirroring [Xlo, Xhi, Ulo, Uhi]. */
template <class T>
struct ScbBounds {
    T Xlo;               ///< lower bound on multiclass throughput X_R (= exact X_1)
    T Xhi;               ///< upper bound on X_R
    std::vector<T> Ulo;  ///< (K) lower bounds on the multiclass utilizations U_k,R
    std::vector<T> Uhi;  ///< (K) upper bounds on U_k,R
};

/**
 * Bracket on the throughput and the per-device utilizations of the UNKNOWN
 * multiclass system whose single-class counterpart has demands L at population N.
 *
 * Theorem 2 / Corollary 2: aggregating an R-class model into its single-class
 * counterpart can only understate performance, U_k,1 <= U_k,R and X_1 <= X_R,
 * and Corollary 1 makes the utilization ratio uniform, U_k,R/U_k,1 = X_R/X_1
 * for every k. Theorem 3 (their Expression 3) caps the relative throughput
 * error at (m-1)/(N+m-1), m = min(N,K), independently of the demands. The
 * single-server capacity U_k,R <= 1 caps the same ratio at 1/(X_1*max(L)),
 * tight on the paper's own worst case, so both are applied.
 *
 * @param L (K) demands of the queueing stations only; a delay station is not
 *          admitted, Theorem 3 resting on the delay-free balanced-network
 *          throughput N/((N+m-1)D)
 * @param N population (N >= 1)
 */
template <class T>
ScbBounds<T> pfqn_scb(const std::vector<T>& L, long N) {
    const std::size_t K = L.size();
    if (K == 0) throw InputError("pfqn_scb: requires at least one queueing station");
    if (N < 1) throw InputError("pfqn_scb: requires N >= 1");
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);

    // Exact single-class MVA at Z=0. This IS the lower bound (Theorem 2), so it
    // is computed exactly rather than bounded: a bounded X1 would not bracket X_R.
    std::vector<T> Q(K, zero), Rk(K, zero);
    T X1 = zero;
    for (long n = 1; n <= N; ++n) {
        T sumR = zero;
        for (std::size_t k = 0; k < K; ++k) {
            Rk[k] = T(L[k] * T(one + Q[k]));
            sumR += Rk[k];
        }
        X1 = T(num_traits<T>::from_int(n) / sumR);
        for (std::size_t k = 0; k < K; ++k) Q[k] = T(X1 * Rk[k]);
    }

    const long m = std::min<long>(N, static_cast<long>(K));
    T ratio = T(num_traits<T>::from_int(N + m - 1) / num_traits<T>::from_int(N));
    T Dmax = L[0];
    for (const T& x : L)
        if (x > Dmax) Dmax = x;
    if (T(X1 * Dmax) > zero) {
        // U_k,R <= 1 with the uniform ratio of Corollary 1. Tight at the worst case.
        const T cap = T(one / T(X1 * Dmax));
        if (cap < ratio) ratio = cap;
    }

    ScbBounds<T> b;
    b.Xlo = X1;
    b.Xhi = T(X1 * ratio);
    b.Ulo.resize(K);
    b.Uhi.resize(K);
    for (std::size_t k = 0; k < K; ++k) {
        b.Ulo[k] = T(X1 * L[k]);
        b.Uhi[k] = T(b.Ulo[k] * ratio);
    }
    return b;
}

/**
 * Demand-free bound on the relative throughput error incurred when r of the N
 * single-customer classes are merged into one class. With r = N this is the
 * full single-class aggregation error of their Theorem 3, at most 50%; with
 * r < N it is the partial-aggregation error of their Theorem 4. The bound never
 * reads the demands, so it can be attached as a certified error bar to any
 * result computed on merged chains.
 *
 * General case, dominating classes allowed (Expression 4, and with r = N
 * Expression 3): e = (min(r,K)-1)/(r+min(r,K)-1). Undominated case, every
 * customer placing the same total demand (Theorem 5 and its comment (3), which
 * lifts the N = R restriction): e = r(r-1)/(min(N,K)(2r-1)), valid for r <= K
 * only, smaller than the general case by the factor r/min(N,K) and equal to it
 * at r = K. THE DOMAIN IS NOT COSMETIC: Theorem 5 gives each of its R classes a
 * dedicated device, so r never exceeds K there, and comment (3) states the
 * generalization for r < K. Evaluated at r > K the expression climbs past the
 * general bound and past the 50% cap of Theorem 3, i.e. it stops being a bound,
 * so r > K is refused rather than returned.
 *
 * @param N total customers, one per class @param K devices
 * @param r classes merged into one (1 <= r <= N)
 * @param undominated true for the tighter Theorem-5 form, valid only when every
 *        customer's total device demand is equal, and only for r <= K
 */
template <class T>
T pfqn_scbgap(long N, long K, long r, bool undominated) {
    if (N < 1 || K < 1) throw InputError("pfqn_scbgap: requires N >= 1 and K >= 1");
    if (r < 1 || r > N) throw InputError("pfqn_scbgap: requires 1 <= r <= N");
    if (r == 1) return num_traits<T>::from_int(0);  // merging one class changes nothing
    if (undominated) {
        if (r > K)
            throw InputError(
                "pfqn_scbgap: the undominated (Theorem 5) form is defined for r <= K only; "
                "beyond it the expression exceeds the general bound and the 50% cap");
        // Theorem 5, comment (3)
        return T(num_traits<T>::from_int(r * (r - 1)) /
                 num_traits<T>::from_int(std::min<long>(N, K) * (2 * r - 1)));
    }
    const long m = std::min<long>(r, K);
    // Expression (4); r=N gives (3)
    return T(num_traits<T>::from_int(m - 1) / num_traits<T>::from_int(r + m - 1));
}

/** Full single-class aggregation: r = N, dominating classes allowed. */
template <class T>
T pfqn_scbgap(long N, long K) {
    return pfqn_scbgap<T>(N, K, N, false);
}

/**
 * Largest value the sum of device utilizations can take in any closed
 * product-form network with R classes, K devices and N customers (Theorem 6):
 * sum_k U_k,R <= (H-1) + (K-H+1)(N-H+1)/(K+N-2H+1), H = min(R,K). Demand-free
 * and nondecreasing in R, which is what makes it invertible into a lower bound
 * on the number of necessary classes; see pfqn_minclasses. The paper's worked
 * case is K = 2, N = 3, R = 1, giving 2N/(N+1) = 1.5.
 */
template <class T>
T pfqn_usumbound(long R, long K, long N) {
    if (N < 1 || K < 1) throw InputError("pfqn_usumbound: requires N >= 1 and K >= 1");
    if (R < 1 || R > N) throw InputError("pfqn_usumbound: requires 1 <= R <= N");
    const long H = std::min<long>(R, K);
    return T(num_traits<T>::from_int(H - 1) +
             T(num_traits<T>::from_int((K - H + 1) * (N - H + 1)) /
               num_traits<T>::from_int(K + N - 2 * H + 1)));
}

/**
 * Smallest number of customer classes R consistent with an observed sum of
 * device utilizations, by inverting the nondecreasing pfqn_usumbound. Only
 * measured quantities are needed -- the utilizations, the device count and the
 * population -- so the answer is available BEFORE any class-specific demand has
 * been characterized. An upper bound on R is meaningless and none is returned.
 *
 * The paper's example: K = 2, N = 3, Usum = 1.6 -> 2, a single class admitting
 * at most 2N/(N+1) = 1.5.
 *
 * @return least R in 1..N with pfqn_usumbound(R,K,N) >= Usum, or -1 where Usum
 *         exceeds min(N,K) and so is unattainable by ANY class structure. -1 is
 *         this port's integral encoding of the NaN MATLAB and Python return.
 */
template <class T>
long pfqn_minclasses(const T& Usum, long K, long N) {
    if (N < 1 || K < 1) throw InputError("pfqn_minclasses: requires N >= 1 and K >= 1");
    if (Usum < num_traits<T>::from_int(0))
        throw InputError("pfqn_minclasses: requires a nonnegative utilization sum");
    // The reference's relative slack, reproduced so a boundary case decides the
    // same way in both codebases; in exact arithmetic it only widens by 1e-12.
    const T one = num_traits<T>::from_int(1);
    const T scale = Usum > one ? Usum : one;
    const T tol = T(num_traits<T>::from_double(1e-12) * scale);
    for (long R = 1; R <= N; ++R)
        if (!(pfqn_usumbound<T>(R, K, N) < T(Usum - tol))) return R;
    return -1;
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_PFQN_SCB_H
