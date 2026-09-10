/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_PAS_IS_H
#define LINE_API_PFQN_PAS_IS_H

/**
 * Importance-sampling estimate of the normalizing constant of a single
 * communicating class of a cyclic two-station pass-and-swap (P&S) network with
 * swap graph H.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_pas_is.m together with its
 * placement-order helper matlab/src/api/pfqn/pas_placement.m, cross-checked
 * against jar/src/main/java/jline/api/pfqn/nc/Pfqn_pas_is.java.
 *
 * Model. Two order-independent stations (1 upstream, 2 downstream) hold all N
 * jobs. With a non-empty swap graph the ordered-state chain is reducible and
 * the recurrent communicating class D is the set of orderings that are
 * non-decreasing with respect to H (Comte and Dorsman, 2021). Writing Phi_m
 * for the balanced-fairness balance function,
 *
 *   G_C = sum_{c in D} sum_{k=0}^{ell} Phi_1(c_{1..k}) Phi_2(c_{ell..k+1}),
 *   Phi_m(q) = prod_{p=1}^{|q|} 1 / mu_m(n(q_{1..p})),  n(.) = prefix counts,
 *
 * which depends on the ordering only through the counts reached at each
 * position; that is the order-independence property.
 *
 * Auto-normalized IS. Orderings are drawn from D by placing, at each step, a
 * uniformly random placement-order-minimal present class. The SAME samples
 * feed numerator and denominator: with xi = 1 the estimate is G_C, and with
 * xi = (number of class-r jobs in the prefix) the ratio is the mean class-r
 * queue length at station 1. Auto-normalized IS is consistent but biased at
 * finite sample count for the ratio Q, while G_C itself is unbiased.
 *
 * Deviation from the reference, deliberate: 1/p(c) is accumulated as the
 * product of the integer branching factors rather than as exp(-sum log na).
 * Same number, no log/exp round trip.
 *
 * Arithmetic: INEXACT BY CONSTRUCTION. The output is a random variable, so no
 * arithmetic makes it exact; the gate is also required by lG = log(G).
 *
 * RNG contract: see pfqn_mc_common.h. Comparable to MATLAB only in
 * distribution, never stream for stream; reproducible within this port only
 * when the generator is passed in the same state.
 */

#include <cstddef>
#include <functional>
#include <vector>

#include "line/api/pfqn/pfqn_mc_common.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** Return value of pfqn_pas_is / pfqn_oi_is, mirroring [G, lG, Q]. */
template <class T>
struct PasIsResult {
    T G;           ///< estimate of the communicating-class normalizing constant
    double lG;     ///< log of the estimate
    Matrix<T> Q;   ///< (2 x R) mean per-class queue length, Q(1,:) = N - Q(0,:)
};

/**
 * The OI rank rate of a station as a function of the per-class COUNT vector:
 * the svcRateFun of an OI / P&S node. The argument is the (R) vector of job
 * counts of the prefix, exactly the `occ` row the MATLAB handles receive.
 * OI property P1 makes mu permutation-invariant, i.e. a function of the
 * counts; it is NOT in general a function of the support alone (an INF
 * station has mu(n) = sum_r n_r sigma_r).
 */
template <class T>
using OiRateFun = std::function<T(const std::vector<int>&)>;

/**
 * matlab/src/api/pfqn/pas_placement.m: transitive closure of the "must
 * precede" relation. P(i,j) is true iff class i must be placed before class j.
 * An empty or all-zero H yields the all-false closure, i.e. no constraint.
 */
inline Matrix<int> pas_placement(const Matrix<int>& H) {
    if (H.empty()) return Matrix<int>();
    const std::size_t R = H.rows();
    if (H.cols() != R) throw InputError("pas_placement: H must be square");
    Matrix<int> P(R, R, 0);
    for (std::size_t i = 0; i < R; ++i)
        for (std::size_t j = 0; j < R; ++j) P(i, j) = H(i, j) != 0 ? 1 : 0;
    for (std::size_t it = 0; it < R; ++it) {
        Matrix<int> Pn(R, R, 0);
        bool changed = false;
        for (std::size_t i = 0; i < R; ++i)
            for (std::size_t j = 0; j < R; ++j) {
                int v = P(i, j);
                if (!v)
                    for (std::size_t k = 0; k < R && !v; ++k)
                        if (P(i, k) && H(k, j) != 0) v = 1;
                Pn(i, j) = v;
                if (v != P(i, j)) changed = true;
            }
        P = Pn;
        if (!changed) break;
    }
    return P;
}

/**
 * @param N       (R) closed population vector
 * @param mu      the two OI rank-rate functions, station 1 then station 2
 * @param H       (R x R) swap-graph adjacency; empty or all zero for pure OI
 * @param samples number of importance samples
 * @param rng     explicit generator, advanced by the call
 * @param want_qlen estimate the per-class queue lengths as well as the
 *        constant. False estimates ONLY G: the prefix-count matrix is neither
 *        allocated nor written and its coefficients are not accumulated, and Q
 *        comes back zero. The ordering is drawn from the same stream either
 *        way, so G is unchanged to the last bit -- this is for the callers that
 *        want G(N - e_r) and read nothing else from it.
 */
template <class T>
PasIsResult<T> pfqn_pas_is(const std::vector<int>& N, const std::vector<OiRateFun<T>>& mu,
                           const Matrix<int>& H, std::size_t samples, McRng& rng,
                           bool want_qlen = true) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_pas_is requires transcendental arithmetic: it is a Monte Carlo estimator, "
                  "inexact by construction, and reports the log of its own estimate");

    const std::size_t R = N.size();
    if (mu.size() != 2)
        throw InputError(
            "pfqn_pas_is models a two-station pass-and-swap tandem: mu must have exactly two rate "
            "functions");
    for (int n : N)
        if (n < 0) throw InputError("pfqn_pas_is: negative population");
    if (!H.empty() && (H.rows() != R || H.cols() != R))
        throw InputError("pfqn_pas_is: H must be a (R x R) swap-graph adjacency matrix");
    if (samples == 0) throw InputError("pfqn_pas_is: at least one sample is required");

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    PasIsResult<T> res;
    res.Q = Matrix<T>(2, R, zero);
    long ell_l = 0;
    for (int n : N) ell_l += n;
    if (ell_l == 0) {
        res.G = one;
        res.lG = 0.0;
        return res;
    }
    const std::size_t ell = static_cast<std::size_t>(ell_l);

    const Matrix<int> P = pas_placement(H);

    // accum[0] carries xi = 1; accum[1 + r] carries xi = n_{1,r}. Without the
    // queue lengths only the first coefficient exists, and cnt1 -- an O(ell R)
    // write per sample that nothing else reads -- is not allocated at all.
    std::vector<T> accum(want_qlen ? R + 1 : 1, zero);
    std::vector<int> x(R), occ(R), occ2(R), avail(R);
    std::vector<std::size_t> c(ell);
    std::vector<T> Phi1(ell + 1), Phi2cut(ell + 1);
    Matrix<T> cnt1(want_qlen ? ell + 1 : 0, want_qlen ? R : 0, zero);

    for (std::size_t s = 0; s < samples; ++s) {
        // ---- draw a feasible ordering from D ------------------------------
        x = N;
        T invp = one;
        for (std::size_t p = 0; p < ell; ++p) {
            std::size_t na = 0;
            for (std::size_t j = 0; j < R; ++j) {
                if (x[j] <= 0) continue;
                bool blocked = false;
                if (!P.empty())
                    for (std::size_t i = 0; i < R && !blocked; ++i)
                        if (x[i] > 0 && P(i, j)) blocked = true;
                if (!blocked) avail[na++] = j;
            }
            if (na == 0)
                throw NumericError(
                    "pfqn_pas_is: swap graph induces no feasible ordering (cyclic placement "
                    "order)");
            const std::size_t pick = avail[mc_uniform_int(rng, na)];
            c[p] = pick;
            invp *= num_traits<T>::from_int(static_cast<long>(na));
            x[pick] -= 1;
        }

        // ---- prefix balance at station 1, and the per-class prefix counts ---
        occ.assign(R, 0);
        Phi1[0] = one;
        if (want_qlen)
            for (std::size_t r = 0; r < R; ++r) cnt1(0, r) = zero;
        T phi = one;
        for (std::size_t k = 0; k < ell; ++k) {
            const std::size_t cls = c[k];
            occ[cls] += 1;
            const T rate = mu[0](occ);
            if (rate == zero) throw NumericError("pfqn_pas_is: zero rank rate at station 1");
            phi /= rate;
            Phi1[k + 1] = phi;
            if (want_qlen)
                for (std::size_t r = 0; r < R; ++r)
                    cnt1(k + 1, r) = num_traits<T>::from_int(occ[r]);
        }

        // ---- reversed-suffix balance at station 2 --------------------------
        occ2.assign(R, 0);
        phi = one;
        Phi2cut[ell] = one;
        for (std::size_t k = ell; k >= 1; --k) {
            occ2[c[k - 1]] += 1;
            const T rate = mu[1](occ2);
            if (rate == zero) throw NumericError("pfqn_pas_is: zero rank rate at station 2");
            phi /= rate;
            Phi2cut[k - 1] = phi;
        }

        // ---- split convolution over every cut ------------------------------
        for (std::size_t k = 0; k <= ell; ++k) {
            const T w2 = k >= ell ? one : Phi2cut[k];
            const T w = Phi1[k] * w2;
            accum[0] += w * invp;
            if (want_qlen && k > 0)
                for (std::size_t r = 0; r < R; ++r) accum[1 + r] += w * cnt1(k, r) * invp;
        }
    }

    const T ns = num_traits<T>::from_int(static_cast<long>(samples));
    res.G = accum[0] / ns;
    res.lG = num_traits<T>::log_as_double(res.G);
    if (want_qlen)
        for (std::size_t r = 0; r < R; ++r) {
            const T q1 = res.G > zero ? (accum[1 + r] / ns) / res.G : zero;
            res.Q(0, r) = q1;
            res.Q(1, r) = num_traits<T>::from_int(N[r]) - q1;
        }
    return res;
}

/** Reference default of 1e4 samples. */
template <class T>
PasIsResult<T> pfqn_pas_is(const std::vector<int>& N, const std::vector<OiRateFun<T>>& mu,
                           const Matrix<int>& H, McRng& rng, bool want_qlen = true) {
    return pfqn_pas_is(N, mu, H, static_cast<std::size_t>(10000), rng, want_qlen);
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_PAS_IS_H
