/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_PFQN_GRNMOL_H
#define LINE_API_PFQN_PFQN_GRNMOL_H

/**
 * Normalizing constant by the closed-form Grundmann-Moeller rule.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_grnmol.m, which applies the rule
 * directly rather than through the successive-degree recursion of pfqn_cub:
 *
 *   G = [ (sum N + M - 1)! / prod_r N_r! ] * sum_{i=0}^{S} w_i H_i
 *   w_i = 2^{-2S} (-1)^i c_i^{2S+1} / (i! (i + c_i)!),   c_i = 2(S-i) + M
 *   H_i = sum_{|b| = S-i} prod_r ( ((2b+1)/c_i)' L(:,r) )^{N_r}
 *
 * with b ranging over the M-vectors of non-negative integers summing to S-i,
 * which is what the reference's call to matlab/util/sprod.m enumerates. The
 * point (2b+1)/c_i is barycentric by construction, since sum_m (2b_m+1) = c_i.
 *
 * REFERENCE DEFECT, reproduced as a rejection rather than as a wrong number.
 * pfqn_grnmol.m sets S = ceil(sum(N)-1)/2, which MATLAB parses as
 * (ceil(sum(N)-1))/2 and NOT as ceil((sum(N)-1)/2), the value pfqn_cub.m uses.
 * For an EVEN total population S is therefore a half-integer, c_i is a
 * half-integer, and the weight calls factorial(i + c_i) on a non-integer,
 * which MATLAB rejects outright ("N must be a matrix of non-negative
 * integers"). pfqn_grnmol is thus callable only for an ODD total population,
 * and the port throws InputError for an even one instead of silently choosing
 * one of the two readings of the expression. Use pfqn_cub for even
 * populations; it computes the same integral with the correct degree.
 *
 * ARITHMETIC. Every ingredient -- the barycentric points, the binomial-like
 * weights, the integer powers -- is rational, and the factorial prefactor is
 * formed as an exact factorial rather than through gammaln. The routine is
 * therefore EXACT in rational arithmetic and is deliberately left ungated:
 * at the reference degree it returns the exact normalizing constant, and
 * checking that against pfqn_ca is the sharpest test the rule admits.
 */

#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

namespace detail {

/**
 * Advance b to the next M-vector of non-negative integers with a fixed sum,
 * in the order matlab/util/sprod.m produces through multichoose. Returns
 * false once the enumeration is exhausted.
 */
inline bool next_composition(std::vector<long>& b, long total) {
    const std::size_t M = b.size();
    if (M < 2) return false;
    // Find the rightmost position that can be decremented with something to
    // its right to absorb the unit.
    for (std::size_t i = M - 1; i-- > 0;) {
        if (b[i] > 0) {
            b[i] -= 1;
            long rest = total;
            for (std::size_t k = 0; k <= i; ++k) rest -= b[k];
            for (std::size_t k = i + 1; k < M; ++k) b[k] = 0;
            b[i + 1] = rest;
            return true;
        }
    }
    return false;
}

}  // namespace detail

/**
 * @param L (M x R) demands, @param N (R) population with an ODD total
 * @return  the normalizing constant
 */
template <class T>
T pfqn_grnmol(const Matrix<T>& L, const std::vector<int>& N) {
    const std::size_t M = L.rows(), R = L.cols();
    if (N.size() != R) throw InputError("pfqn_grnmol: L and N disagree on the class count");
    if (M == 0) throw InputError("pfqn_grnmol: empty demand matrix");
    long Nt = 0;
    for (int v : N) {
        if (v < 0) throw InputError("pfqn_grnmol: negative population");
        Nt += v;
    }
    if (Nt == 0) return num_traits<T>::from_int(1);
    if (Nt % 2 == 0)
        throw InputError(
            "pfqn_grnmol: the reference is only callable for an odd total population "
            "(S = ceil(sum(N)-1)/2 is a half-integer otherwise); use pfqn_cub");
    const long S = (Nt - 1) / 2;
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);

    T G = zero;
    for (long i = 0; i <= S; ++i) {
        const long ci = 2 * (S - i) + static_cast<long>(M);
        T w = T(one / num_pow_int(num_traits<T>::from_int(2), static_cast<unsigned>(2 * S)));
        if (i % 2 == 1) w = T(-w);
        w *= num_pow_int(num_traits<T>::from_int(ci), static_cast<unsigned>(2 * S + 1));
        w /= num_factorial<T>(static_cast<unsigned>(i));
        w /= num_factorial<T>(static_cast<unsigned>(i + ci));

        const long tot = S - i;
        std::vector<long> b(M, 0);
        b[0] = tot;
        T H = zero;
        while (true) {
            T prod = one;
            for (std::size_t r = 0; r < R; ++r) {
                if (N[r] == 0) continue;
                T uL = zero;
                for (std::size_t m = 0; m < M; ++m)
                    uL += T(num_traits<T>::from_int(2 * b[m] + 1) / num_traits<T>::from_int(ci)) *
                          L(m, r);
                prod *= num_pow_int(uL, static_cast<unsigned>(N[r]));
            }
            H += prod;
            if (!detail::next_composition(b, tot)) break;
        }
        G += w * H;
    }

    G *= num_factorial<T>(static_cast<unsigned>(Nt + static_cast<long>(M) - 1));
    for (std::size_t r = 0; r < R; ++r) G /= num_factorial<T>(static_cast<unsigned>(N[r]));
    return G;
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_PFQN_GRNMOL_H
