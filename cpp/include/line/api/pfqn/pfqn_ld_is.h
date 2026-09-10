/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_LD_IS_H
#define LINE_API_PFQN_LD_IS_H

/**
 * Importance-sampling estimate of the normalizing constant of a closed
 * LOAD-DEPENDENT product-form network.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_ld_is.m, cross-checked against
 * jar/src/main/java/jline/api/pfqn/nc/Pfqn_ld_is.java.
 *
 * Identity. Every product-form station's balance function is the sum, over the
 * orderings q of a given per-class count vector n, of an ordered product of a
 * per-position factor,
 *
 *   F_i(n) = |n|!/prod_r(n_r!) prod_r L(i,r)^{n_r} / prod_{k=1}^{|n|} mu_i(k)
 *          = sum_{q: |q| = n} prod_{p=1}^{|n|} L(i,q_p) / mu_i(p),
 *
 * so with ell = sum(N) and a cut vector splitting an ordering c of all ell
 * jobs into S contiguous segments (one per station),
 *
 *   G(N) = sum_c sum_{cuts} prod_m prod_p L(m, seg_m(p)) / mu_m(p).
 *
 * The delay is the station with mu_Z(k) = k, a single server is mu_i(k) = 1
 * and a c-server queue is mu_i(k) = min(k,c).
 *
 * Estimator. An ordering c is drawn by placing a uniformly random present
 * class at each step; for the sampled c the inner sum over ALL cut vectors is
 * evaluated EXACTLY by the dynamic program A_0(0) = 1,
 * A_m(k) = sum_{j<=k} A_{m-1}(j) w_m(c_{j+1..k}), in O(S ell^2). The estimate
 * is the sample mean of A_S(ell)/p(c), which is unbiased for G(N).
 *
 * Deviation from the reference, deliberate: MATLAB accumulates log p(c) and
 * multiplies by exp(-logp). Here 1/p(c) is accumulated directly as the product
 * of the integer branching factors, which is the same number computed without
 * a log/exp round trip and is exact in every arithmetic. The two agree to
 * rounding in double.
 *
 * Arithmetic: INEXACT BY CONSTRUCTION. The estimator's output is a random
 * variable whose value depends on the drawn orderings, so no arithmetic makes
 * it exact; it is gated on has_transcendental because it reports lG = log(G)
 * and because the sampler itself needs a real-valued uniform stream.
 *
 * RNG contract: see pfqn_mc_common.h. Comparable to MATLAB only in
 * distribution, never stream for stream; reproducible within this port only
 * when the generator is passed in the same state.
 */

#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_ca.h"
#include "line/api/pfqn/pfqn_mc_common.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/**
 * @param L       (M x R) per-class demands at the M queueing stations
 * @param N       (R) closed population vector
 * @param Z       (R) aggregated think times; empty or all zero for no delay
 * @param mu      (M x k) load-dependent capacities, mu(i,k-1) with k jobs at
 *                station i; empty for the load-independent case mu = 1. A
 *                short row is extended with its last entry, as in the
 *                reference.
 * @param samples number of importance samples
 * @param rng     explicit generator, advanced by the call
 */
template <class T>
NcResult<T> pfqn_ld_is(const Matrix<T>& L, const std::vector<int>& N, const std::vector<T>& Z,
                       const Matrix<T>& mu, std::size_t samples, McRng& rng) {
    static_assert(num_traits<T>::has_transcendental,
                  "pfqn_ld_is requires transcendental arithmetic: it is a Monte Carlo estimator, "
                  "inexact by construction, and reports the log of its own estimate");

    const std::size_t M = L.empty() ? 0 : L.rows();
    const std::size_t R = N.size();
    if (!L.empty() && L.cols() != R)
        throw InputError("pfqn_ld_is: L must have as many columns as N has classes");
    if (!Z.empty() && Z.size() != R) throw InputError("pfqn_ld_is: Z has the wrong length");
    for (int n : N)
        if (n < 0) throw InputError("pfqn_ld_is: negative population");
    if (samples == 0) throw InputError("pfqn_ld_is: at least one sample is required");

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    long ell_l = 0;
    for (int n : N) ell_l += n;
    if (ell_l == 0) return {one, 0.0};
    const std::size_t ell = static_cast<std::size_t>(ell_l);

    // ---- station list: the M queues, plus the delay as mu_Z(k) = k ---------
    bool hasZ = false;
    for (std::size_t r = 0; r < Z.size(); ++r)
        if (Z[r] > zero) hasZ = true;
    const std::size_t S = M + (hasZ ? 1u : 0u);
    if (S == 0) throw InputError("pfqn_ld_is: no station carries any demand");

    Matrix<T> D(S, R, zero);
    Matrix<T> B(S, ell, one);
    for (std::size_t i = 0; i < M; ++i) {
        for (std::size_t r = 0; r < R; ++r) D(i, r) = L(i, r);
        if (mu.empty()) continue;
        if (mu.rows() != M) throw InputError("pfqn_ld_is: mu has the wrong station count");
        const std::size_t kmax = mu.cols() < ell ? mu.cols() : ell;
        for (std::size_t k = 0; k < kmax; ++k) B(i, k) = mu(i, k);
        for (std::size_t k = kmax; k < ell; ++k) B(i, k) = mu(i, mu.cols() - 1);
    }
    if (hasZ) {
        for (std::size_t r = 0; r < R; ++r) D(S - 1, r) = Z[r];
        for (std::size_t k = 0; k < ell; ++k) B(S - 1, k) = num_traits<T>::from_int(static_cast<long>(k) + 1);
    }
    for (std::size_t i = 0; i < S; ++i)
        for (std::size_t k = 0; k < ell; ++k)
            if (!(B(i, k) > zero))
                throw InputError("pfqn_ld_is: load-dependent capacities must be strictly positive");

    T acc = zero;
    std::vector<int> x(R);
    std::vector<std::size_t> c(ell);
    std::vector<std::size_t> avail(R);
    std::vector<T> A(ell + 1), Anew(ell + 1);

    for (std::size_t s = 0; s < samples; ++s) {
        // ---- draw an ordering, uniformly random present class at each step --
        x = N;
        T invp = one;  // 1/p(c) = product of the branching factors
        for (std::size_t p = 0; p < ell; ++p) {
            std::size_t na = 0;
            for (std::size_t r = 0; r < R; ++r)
                if (x[r] > 0) avail[na++] = r;
            const std::size_t pick = avail[mc_uniform_int(rng, na)];
            c[p] = pick;
            invp *= num_traits<T>::from_int(static_cast<long>(na));
            x[pick] -= 1;
        }

        // ---- exact inner sum over every cut vector, by dynamic programming --
        A.assign(ell + 1, zero);
        A[0] = one;
        for (std::size_t m = 0; m < S; ++m) {
            Anew.assign(ell + 1, zero);
            for (std::size_t j = 0; j <= ell; ++j) {
                if (A[j] == zero) continue;
                Anew[j] += A[j];  // empty segment at station m
                T w = one;
                for (std::size_t k = j; k < ell; ++k) {
                    w *= D(m, c[k]) / B(m, k - j);
                    if (w == zero) break;
                    Anew[k + 1] += A[j] * w;
                }
            }
            A.swap(Anew);
        }
        acc += A[ell] * invp;
    }

    const T G = acc / num_traits<T>::from_int(static_cast<long>(samples));
    return {G, num_traits<T>::log_as_double(G)};
}

/** Reference default of 1e4 samples. */
template <class T>
NcResult<T> pfqn_ld_is(const Matrix<T>& L, const std::vector<int>& N, const std::vector<T>& Z,
                       const Matrix<T>& mu, McRng& rng) {
    return pfqn_ld_is(L, N, Z, mu, static_cast<std::size_t>(10000), rng);
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_LD_IS_H
