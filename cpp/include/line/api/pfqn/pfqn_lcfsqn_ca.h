/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_LCFSQN_CA_H
#define LINE_API_PFQN_LCFSQN_CA_H

/**
 * Convolution algorithm for the two-station multiclass LCFS queueing network
 * of Casale, "A family of multiclass LCFS queueing networks with
 * order-dependent product-form solutions", QUESTA 2026.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_lcfsqn_ca.m.
 *
 * Station 1 is LCFS (non-preemptive) and station 2 is LCFS-PR
 * (preemptive-resume); alpha_r and beta_r are the class-r mean service times
 * at the two stations. The network is NOT product form in the classical sense
 * -- the balance function depends on the ORDER of the jobs in the LCFS queue,
 * not only on their counts -- so the constant obeys a coupled pair of
 * recursions over the population lattice rather than a single Buzen update:
 *
 *   V(0) = 1,   V(n) = prod_r alpha_r^{n_r} * sum_{r: n_r > 0} V(n - e_r)
 *   G(0) = 1,   G(n) = sum_{r: n_r > 0} alpha_r^{|n|-1} beta_r G(n - e_r) + V(n)
 *
 * V is the order-dependent auxiliary term contributed by the non-preemptive
 * station; it is returned alongside G because the mean-value routine
 * pfqn_lcfsqn_mva needs it.
 *
 * Arithmetic: EXACT-CAPABLE. Only integer powers, additions and
 * multiplications in the field of the inputs.
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
struct LcfsQnResult {
    T G;  ///< normalizing constant
    T V;  ///< order-dependent auxiliary term
};

/**
 * @param alpha (R) mean service times at the LCFS station
 * @param beta  (R) mean service times at the LCFS-PR station
 * @param N     (R) population per class
 */
template <class T>
LcfsQnResult<T> pfqn_lcfsqn_ca(const std::vector<T>& alpha, const std::vector<T>& beta,
                               const std::vector<int>& N) {
    const std::size_t R = alpha.size();
    if (beta.size() != R) throw InputError("pfqn_lcfsqn_ca: alpha and beta have different lengths");
    if (N.size() != R) throw InputError("pfqn_lcfsqn_ca: alpha and N have different lengths");

    const T one = num_traits<T>::from_int(1);
    const T zero = num_traits<T>::from_int(0);
    long K = 0;
    for (int v : N) {
        if (v < 0) throw InputError("pfqn_lcfsqn_ca: negative population");
        K += v;
    }
    if (K == 0) return {one, one};

    const std::vector<std::size_t> prods = plane_sizes(N);
    const std::size_t total = population_count(N);
    std::vector<T> G(total, zero), V(total, zero);

    std::vector<int> n(R, 0);
    bool more = true;
    while (more) {
        const std::size_t idx = pop_index(n, prods);
        int tot = 0;
        for (int v : n) tot += v;
        if (tot == 0) {
            G[idx] = one;
            V[idx] = one;
        } else {
            T gv = zero, vv = zero;
            for (std::size_t r = 0; r < R; ++r) {
                if (n[r] == 0) continue;
                const std::size_t idxr = idx - prods[r];
                vv += V[idxr];
                gv += num_pow_int(alpha[r], static_cast<unsigned>(tot - 1)) * beta[r] * G[idxr];
            }
            // V picks up the full prod_r alpha_r^{n_r} factor once per state.
            T apow = one;
            for (std::size_t r = 0; r < R; ++r)
                apow *= num_pow_int(alpha[r], static_cast<unsigned>(n[r]));
            V[idx] = apow * vv;
            G[idx] = gv + V[idx];
        }
        more = next_pop(n, N);
    }
    return {G[total - 1], V[total - 1]};
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_LCFSQN_CA_H
