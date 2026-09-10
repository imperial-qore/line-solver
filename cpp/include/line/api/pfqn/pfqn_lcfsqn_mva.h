/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_LCFSQN_MVA_H
#define LINE_API_PFQN_LCFSQN_MVA_H

/**
 * Exact mean value analysis of the two-station multiclass LCFS network of
 * Casale, QUESTA 2026 (station 1 LCFS, station 2 LCFS-PR).
 *
 * Templated port of matlab/src/api/pfqn/pfqn_lcfsqn_mva.m.
 *
 * The recursion over the population lattice carries, besides the queue lengths
 * Q and the throughputs T, the BACK PROBABILITIES B(s,r) that a class-r job
 * sits at the back of the queue at station s. Writing |n| for the total
 * population, n_k for the class-k population and A = prod_r alpha_r^{n_r},
 *
 *   Wnp = A [ 1 + Q_{n-e_k}(1,k)
 *               + sum_{r != k} (alpha_k/alpha_r) B_{n-e_k}(1,r)/B_{n-e_r}(1,k)
 *                              Q_{n-e_r}(1,k) ]
 *   Wpr = alpha_k^{|n|-1} beta_k [ 1 + Q_{n-e_k}(2,k)
 *               + sum_{r != k} (alpha_r/alpha_k) B_{n-e_k}(2,r)/B_{n-e_r}(2,k)
 *                              Q_{n-e_r}(2,k) ]
 *   B_n(1,k) = A n_k / (Wnp + Wpr),   B_n(2,k) = alpha_k^{|n|-1} beta_k n_k / (Wnp + Wpr)
 *   Q_n(s,k) = B_n(s,k) + sum_r B_n(s,r) Q_{n-e_r}(s,k)
 *   T_n(k)   = sum_r B_n(1,r) T_{n-e_r}(k) + (1/alpha_k) B_n(1,k) (1 - sum_r alpha_r T_{n-e_k}(r))
 *
 * with everything zero at n = 0. The reference special-cases |n| = 1; it is
 * not necessary, since the general step reduces to it (the r != k sum is empty
 * and every n - e_k term is zero), and the port therefore uses one uniform
 * recursion.
 *
 * Arithmetic: EXACT-CAPABLE. The reference is written entirely in a
 * scaled-log representation -- B is stored as a mantissa plus a log scale, and
 * every sum goes through a log-sum-exp helper -- and its own comments state
 * that this is a range device and "mathematically EXACT, no approximations are
 * made". That representation is dropped here: the recursion is evaluated
 * directly in T, which needs no scaling in an exact field and none in the
 * high-precision floats either.
 *
 * ONE SEMANTIC CONSEQUENCE of dropping it. The log-sum-exp helper cannot
 * represent a negative term, so the reference silently skips the throughput
 * contribution whenever (1 - sum_r alpha_r T_{n-e_k}(r)) is negative, and
 * skips any term whose factors are not strictly positive. Those guards are
 * artifacts of the representation, not of the model: the bracket is a
 * probability complement and is nonnegative on a consistent state. This port
 * includes every term unconditionally, which is what the recursion says.
 */

#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"
#include "line/util/population.h"

namespace line {
namespace pfqn {

template <class T>
struct LcfsMvaResult {
    std::vector<T> T_;  ///< (R) per-class throughput
    Matrix<T> Q;        ///< (2 x R) mean queue lengths
    Matrix<T> U;        ///< (2 x R) utilizations
    Matrix<T> B;        ///< (2 x R) back probabilities
};

/**
 * @param alpha (R) mean service times at the LCFS station
 * @param beta  (R) mean service times at the LCFS-PR station
 * @param N     (R) population per class
 */
template <class T>
LcfsMvaResult<T> pfqn_lcfsqn_mva(const std::vector<T>& alpha, const std::vector<T>& beta,
                                 const std::vector<int>& N) {
    const std::size_t R = alpha.size();
    if (beta.size() != R) throw InputError("pfqn_lcfsqn_mva: alpha and beta have different lengths");
    if (N.size() != R) throw InputError("pfqn_lcfsqn_mva: alpha and N have different lengths");

    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);

    LcfsMvaResult<T> res;
    res.T_.assign(R, zero);
    res.Q = Matrix<T>(2, R, zero);
    res.U = Matrix<T>(2, R, zero);
    res.B = Matrix<T>(2, R, zero);

    long K = 0;
    for (int v : N) {
        if (v < 0) throw InputError("pfqn_lcfsqn_mva: negative population");
        K += v;
    }
    if (K == 0) return res;
    for (std::size_t r = 0; r < R; ++r)
        if (!(alpha[r] > zero)) throw InputError("pfqn_lcfsqn_mva: a service time at the LCFS "
                                                 "station is not positive");

    const std::vector<std::size_t> prods = plane_sizes(N);
    const std::size_t total = population_count(N);
    // Per state: Q(2R), B(2R), Tp(R), laid out contiguously.
    std::vector<T> Q(total * 2 * R, zero), B(total * 2 * R, zero), Tp(total * R, zero);
    const auto QB = [&](std::vector<T>& v, std::size_t idx, std::size_t s, std::size_t r) -> T& {
        return v[idx * 2 * R + s * R + r];
    };

    std::vector<int> n(R, 0);
    bool more = next_pop(n, N);  // the origin stays all-zero
    while (more) {
        const std::size_t idx = pop_index(n, prods);
        int tot = 0;
        for (int v : n) tot += v;

        // A = prod_r alpha_r^{n_r}, shared by every class at this state.
        T A = one;
        for (std::size_t r = 0; r < R; ++r)
            A *= num_pow_int(alpha[r], static_cast<unsigned>(n[r]));

        // ---- back probabilities -------------------------------------------
        for (std::size_t k = 0; k < R; ++k) {
            if (n[k] == 0) continue;
            const std::size_t ik = idx - prods[k];
            T wnp = one + QB(Q, ik, 0, k);
            T wpr = one + QB(Q, ik, 1, k);
            for (std::size_t r = 0; r < R; ++r) {
                if (r == k || n[r] == 0) continue;
                const std::size_t ir = idx - prods[r];
                const T dnp = QB(B, ir, 0, k);
                const T dpr = QB(B, ir, 1, k);
                if (dnp == zero || dpr == zero)
                    throw NumericError("pfqn_lcfsqn_mva: a back probability vanished");
                wnp += alpha[k] / alpha[r] * (QB(B, ik, 0, r) / dnp) * QB(Q, ir, 0, k);
                wpr += alpha[r] / alpha[k] * (QB(B, ik, 1, r) / dpr) * QB(Q, ir, 1, k);
            }
            const T Apr = num_pow_int(alpha[k], static_cast<unsigned>(tot - 1)) * beta[k];
            const T W = A * wnp + Apr * wpr;
            if (W == zero) throw NumericError("pfqn_lcfsqn_mva: zero total waiting time");
            const T nk = num_traits<T>::from_int(n[k]);
            QB(B, idx, 0, k) = A * nk / W;
            QB(B, idx, 1, k) = Apr * nk / W;
        }

        // ---- queue lengths -------------------------------------------------
        for (std::size_t k = 0; k < R; ++k) {
            if (n[k] == 0) continue;
            for (std::size_t s = 0; s < 2; ++s) {
                T q = QB(B, idx, s, k);
                for (std::size_t r = 0; r < R; ++r) {
                    if (n[r] == 0) continue;
                    q += QB(B, idx, s, r) * QB(Q, idx - prods[r], s, k);
                }
                QB(Q, idx, s, k) = q;
            }
        }

        // ---- throughputs ---------------------------------------------------
        for (std::size_t k = 0; k < R; ++k) {
            if (n[k] == 0) continue;
            const std::size_t ik = idx - prods[k];
            T unp = zero;
            for (std::size_t r = 0; r < R; ++r) unp += alpha[r] * Tp[ik * R + r];
            T t = QB(B, idx, 0, k) * (one - unp) / alpha[k];
            for (std::size_t r = 0; r < R; ++r) {
                if (n[r] == 0) continue;
                t += QB(B, idx, 0, r) * Tp[(idx - prods[r]) * R + k];
            }
            Tp[idx * R + k] = t;
        }

        more = next_pop(n, N);
    }

    const std::size_t last = total - 1;
    for (std::size_t r = 0; r < R; ++r) {
        res.T_[r] = Tp[last * R + r];
        for (std::size_t s = 0; s < 2; ++s) {
            res.Q(s, r) = QB(Q, last, s, r);
            res.B(s, r) = QB(B, last, s, r);
        }
        res.U(0, r) = res.T_[r] * alpha[r];
        res.U(1, r) = res.T_[r] * beta[r];
    }
    return res;
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_LCFSQN_MVA_H
