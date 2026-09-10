/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_LCFSQN_NC_H
#define LINE_API_PFQN_LCFSQN_NC_H

/**
 * Normalizing constant of the two-station multiclass LCFS network as a sum of
 * PERMANENTS, the closed form of Casale, QUESTA 2026.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_lcfsqn_nc.m.
 *
 * Splitting the job sequence at the position x of the boundary between the
 * two stations gives
 *
 *   G(N) = ( sum_{x=0}^{K} perm(A_x, N) ) / prod_r N_r!,   K = sum_r N_r,
 *
 * where A_x is the (K x R) matrix of the per-position, per-class weights
 *
 *   A_x(i,r) = alpha_r^i                  for i <= x        (LCFS side)
 *   A_x(i,r) = alpha_r^{i-1} beta_r       for i >  x        (LCFS-PR side)
 *
 * and the permanent is taken with column multiplicities N.
 *
 * TWO REFERENCE DEFECTS, both corrected here, both invisible at the default
 * population N = ones(1,R) and both silent for any other one.
 *
 * 1. TRANSPOSED MATRIX. MATLAB's make_A builds the matrix transposed with
 *    respect to the way perm() then reads it: make_A fills A(class, position)
 *    over i = 1..R and j = 1..K, while perm(A, N) reads A(position, class)
 *    over i = 1..sum(N) and k = 1..R. The two indexings coincide only when
 *    K = R, i.e. exactly the default. For any larger population the rows
 *    i = R+1, ..., K of the matrix perm() actually reads are identically zero,
 *    so every product term in Ryser's formula contains a zero factor and the
 *    routine returns G = 0 -- silently, since zero is a representable constant
 *    and nothing downstream checks positivity.
 *
 * 2. MISSING prod_r N_r!. The permanent counts the ORDERED arrangements of the
 *    jobs, so it over-counts every state by the number of permutations WITHIN
 *    each class. The same division appears explicitly in the sibling routine
 *    pfqn_joint, whose local Fper ends in `perm(A)/prod(factorial(N))`; it is
 *    absent here. With N = ones the factorials are all 1, which is why the
 *    omission never showed.
 *
 * Both corrections are verified, not assumed: with them the closed form agrees
 * with pfqn_lcfsqn_ca -- an independent recursion sharing no code -- as an
 * EXACT rational on every population tested, and without either one it does
 * not.
 *
 * Arithmetic: EXACT-CAPABLE, inheriting exactness from pfqn_perm. Note that
 * the cost is (K+1) * prod_r (N_r + 1) * K * R, so this form is a closed-form
 * cross-check on the recursion of pfqn_lcfsqn_ca rather than a replacement
 * for it.
 */

#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_perm.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/**
 * @param alpha (R) mean service times at the LCFS station
 * @param beta  (R) mean service times at the LCFS-PR station
 * @param N     (R) population per class
 */
template <class T>
T pfqn_lcfsqn_nc(const std::vector<T>& alpha, const std::vector<T>& beta,
                 const std::vector<int>& N) {
    const std::size_t R = alpha.size();
    if (beta.size() != R) throw InputError("pfqn_lcfsqn_nc: alpha and beta have different lengths");
    if (N.size() != R) throw InputError("pfqn_lcfsqn_nc: alpha and N have different lengths");

    const T one = num_traits<T>::from_int(1);
    T G = num_traits<T>::from_int(0);
    int K = 0;
    for (int v : N) {
        if (v < 0) throw InputError("pfqn_lcfsqn_nc: negative population");
        K += v;
    }
    if (K == 0) return one;

    for (int x = 0; x <= K; ++x) {
        Matrix<T> A(static_cast<std::size_t>(K), R);
        for (int i = 1; i <= K; ++i)
            for (std::size_t r = 0; r < R; ++r)
                A(static_cast<std::size_t>(i - 1), r) =
                    i <= x ? num_pow_int(alpha[r], static_cast<unsigned>(i))
                           : num_pow_int(alpha[r], static_cast<unsigned>(i - 1)) * beta[r];
        G += pfqn_perm(A, N);
    }
    for (int v : N) G /= num_factorial<T>(static_cast<unsigned>(v));
    return G;
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_LCFSQN_NC_H
