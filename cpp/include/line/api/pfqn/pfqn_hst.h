/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_HST_H
#define LINE_API_PFQN_HST_H

/**
 * Operational sensitivity of throughput to homogeneous-service-time (HST)
 * violations, and the constrained worst case (Suri 1983).
 *
 * Templated port of matlab/src/api/pfqn/pfqn_hst.m.
 *
 * A robustness certificate for a single-class closed product-form solution: how
 * far the predicted throughput can move when the HST assumption fails at one
 * station. That assumption states that the mean service time at station i does
 * not depend on the queue length there. Suri (1983) perturbs it to
 * S_i(n) = S_i (1 + a_n), one relative deviation per queue-length level n, and
 * shows (eq. 3.11) that to first order
 *
 *   (1/X0) dX0/da_n = c_n = P(n_i >= n+1)/u_i - P(n_i >= n),
 *
 * with u_i = L_i X0 the station utilization and the marginals taken from the
 * product-form solution, P(n_i >= n) = L_i^n G(N-n)/G(N). The naive certificate
 * |dX0/X0| <= (sum_n |c_n|) d follows from |a_n| <= d alone, and by Lemma 3.1
 * that total equals Q_i(N) - Q_i(N-1).
 *
 * That bound is loose because the deviations are not free: an operationally
 * consistent perturbation must leave the OBSERVED mean service time unchanged,
 * sum_n p_n a_n = 0 with p_n = P(n_i = n). The constrained problem (P1),
 *
 *   max |sum_n c_n a_n|  s.t.  |a_n| <= d,  sum_n p_n a_n = 0,
 *
 * is a one-constraint linear program, solved here exactly: its optimum sets
 * a_n = +/- d according to whether the ratio c_n / p_n exceeds a threshold, with
 * at most one fractional coordinate. On the paper's Figure 1 system it collapses
 * 0.831 d to 0.102 d.
 *
 * Reference: R. Suri, "Robustness of Queuing Network Formulas", JACM
 * 30(3):564-594, 1983 (eq. 3.11, Lemma 3.1, problem (P1)).
 *
 * Arithmetic: TRANSCENDENTAL. The marginals come from pfqn_rgf, whose recursion
 * is carried in the log domain.
 */

#include <algorithm>
#include <cstddef>
#include <limits>
#include <numeric>
#include <vector>

#include "line/api/pfqn/pfqn_rgf.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace pfqn {

/** Everything the HST certificate reports about one station. */
template <class T>
struct HstResult {
    std::size_t station;   ///< 0-based index of the station analysed
    T X;                   ///< product-form throughput
    T U;                   ///< utilization of that station
    T Q;                   ///< mean queue length there
    std::vector<T> Pgeq;   ///< P(n_i >= k), k = 0..N
    std::vector<T> p;      ///< P(n_i = k), k = 0..N
    std::vector<T> c;      ///< sensitivity coefficients c_k, k = 1..N (eq. 3.11)
    T total;               ///< sum_k |c_k|, the unconstrained certificate per unit d
    T worst;               ///< the (P1) optimum per unit d
    std::vector<T> astar;  ///< the worst-case deviation profile a_k / d, k = 1..N
};

/**
 * @param L   (M) service demands of the queueing stations
 * @param N   population, an integer of at least one job
 * @param Z   think time
 * @param ist 0-based station the HST perturbation is applied to
 */
template <class T>
HstResult<T> pfqn_hst(const std::vector<T>& L, int N, const T& Z, std::size_t ist) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_hst reads its marginals from pfqn_rgf and needs transcendental arithmetic");
    using std::exp;
    using std::log;
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    const std::size_t M = L.size();
    if (M == 0) throw InputError("pfqn_hst requires at least one queueing station");
    if (ist >= M)
        throw InputError("pfqn_hst: the station index is out of range for the supplied demands");
    if (N < 1) throw InputError("pfqn_hst requires an integer population of at least one job");
    if (L[ist] <= zero)
        throw InputError(
            "pfqn_hst: the requested station has zero demand, so its queue-length marginals are "
            "degenerate");

    const RgfResult<T> rgf = pfqn_rgf(L, N, Z);  // lg[k] = log G(k), k = 0..N
    const std::size_t Np = static_cast<std::size_t>(N);

    HstResult<T> res;
    res.station = ist;
    res.X = exp(T(rgf.lg[Np - 1] - rgf.lg[Np]));
    const T y = L[ist];
    const T logy = log(y);

    // P(n_i >= k) = y^k G(N-k)/G(N)
    res.Pgeq.assign(Np + 1, zero);
    for (std::size_t k = 0; k <= Np; ++k)
        res.Pgeq[k] = exp(T(num_traits<T>::from_int(static_cast<long>(k)) * logy + rgf.lg[Np - k] -
                            rgf.lg[Np]));
    res.p.assign(Np + 1, zero);
    for (std::size_t k = 0; k <= Np; ++k)
        res.p[k] = T(res.Pgeq[k] - (k + 1 <= Np ? res.Pgeq[k + 1] : zero));

    res.U = T(y * res.X);
    res.Q = zero;
    for (std::size_t k = 1; k <= Np; ++k) res.Q += res.Pgeq[k];

    // eq. (3.11): c_n = P(>= n+1)/u - P(>= n), n = 1..N
    res.c.assign(Np, zero);
    for (std::size_t n = 1; n <= Np; ++n) {
        const T Pn1 = (n + 1 <= Np) ? res.Pgeq[n + 1] : zero;
        res.c[n - 1] = T(T(Pn1 / res.U) - res.Pgeq[n]);
    }
    res.total = zero;
    for (std::size_t k = 0; k < Np; ++k) res.total += (res.c[k] < zero ? T(-res.c[k]) : res.c[k]);

    // (P1): one equality constraint plus a box. At the optimum
    // a_n = sign(c_n - lambda p_n) d, so sorting by the ratio c_n / p_n and
    // sweeping the split point enumerates every candidate lambda; the constraint
    // fixes the single fractional coordinate at the split.
    std::vector<T> pp(Np);
    for (std::size_t n = 0; n < Np; ++n) pp[n] = res.p[n + 1];
    std::vector<std::size_t> ord(Np);
    std::iota(ord.begin(), ord.end(), static_cast<std::size_t>(0));
    const T tiny = num_traits<T>::from_double(std::numeric_limits<double>::min());
    std::vector<T> ratio(Np);
    for (std::size_t n = 0; n < Np; ++n) ratio[n] = T(res.c[n] / (pp[n] > tiny ? pp[n] : tiny));
    std::stable_sort(ord.begin(), ord.end(),
                     [&ratio](std::size_t a, std::size_t b) { return ratio[b] < ratio[a]; });

    T best = zero;
    std::vector<T> astar(Np, zero);
    for (std::size_t k = 0; k <= Np; ++k) {
        std::vector<T> a(Np, T(-one));
        for (std::size_t t = 0; t < k; ++t) a[ord[t]] = one;
        for (std::size_t piv = 0; piv < Np; ++piv) {
            if (pp[piv] <= zero) continue;
            std::vector<T> aa = a;
            T rest = zero;
            for (std::size_t n = 0; n < Np; ++n) rest += pp[n] * aa[n];
            rest -= pp[piv] * aa[piv];
            const T v = T(-rest / pp[piv]);
            if (v < T(-one) || v > one) continue;
            aa[piv] = v;
            T obj = zero;
            for (std::size_t n = 0; n < Np; ++n) obj += res.c[n] * aa[n];
            const T mobj = obj < zero ? T(-obj) : obj;
            const T mbest = best < zero ? T(-best) : best;
            if (mobj > mbest) {
                best = obj;
                astar = aa;
            }
        }
    }
    if (best < zero) {  // the feasible set is symmetric
        best = T(-best);
        for (std::size_t n = 0; n < Np; ++n) astar[n] = T(-astar[n]);
    }
    res.worst = best;
    res.astar = astar;
    return res;
}

/** MATLAB default: the bottleneck station, argmax L. */
template <class T>
HstResult<T> pfqn_hst(const std::vector<T>& L, int N, const T& Z) {
    if (L.empty()) throw InputError("pfqn_hst requires at least one queueing station");
    std::size_t ist = 0;
    for (std::size_t i = 1; i < L.size(); ++i)
        if (L[i] > L[ist]) ist = i;
    return pfqn_hst(L, N, Z, ist);
}

/** MATLAB default: no think time and the bottleneck station. */
template <class T>
HstResult<T> pfqn_hst(const std::vector<T>& L, int N) {
    return pfqn_hst(L, N, num_traits<T>::from_int(0));
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_HST_H
