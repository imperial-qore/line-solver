/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_PERM_H
#define LINE_API_PFQN_PERM_H

/**
 * Permanent of a matrix with repeated columns, by Ryser's formula.
 *
 * Templated port of matlab/src/util/perm.m, which the pfqn family uses in two
 * places: the joint queue-length probability pfqn_joint (through its local
 * Fper) and the LCFS normalizing constant pfqn_lcfsqn_nc.
 *
 * For an (n x R) matrix A whose column k is repeated m_k times, n = sum_k m_k,
 * the permanent of the expanded (n x n) matrix is
 *
 *   perm = (-1)^n sum_{0 <= f <= m} (-1)^{|f|} prod_k C(m_k, f_k)
 *                                   prod_{i=1}^{n} sum_k f_k A(i,k),
 *
 * which costs prod_k (m_k + 1) evaluations rather than the n! of the
 * definition. Collapsing the repeated columns is what makes the formula usable
 * here at all: the queueing applications have n jobs but only R distinct
 * classes, and R is small.
 *
 * Arithmetic: EXACT-CAPABLE. Additions, multiplications and integer binomials
 * only; the binomials are formed by the exact Pascal recurrence of num_nck,
 * not by the floating-point nck, so a large multiplicity does not lose
 * integrality.
 */

#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"
#include "line/util/population.h"

namespace line {
namespace pfqn {

/**
 * @param A (n x R) matrix, column k standing for m_k identical columns
 * @param m (R) column multiplicities, summing to n
 */
template <class T>
T pfqn_perm(const Matrix<T>& A, const std::vector<int>& m) {
    const std::size_t R = m.size();
    if (A.cols() != R) throw InputError("pfqn_perm: A and the multiplicities disagree in width");
    int n = 0;
    for (int v : m) {
        if (v < 0) throw InputError("pfqn_perm: negative multiplicity");
        n += v;
    }
    if (static_cast<int>(A.rows()) != n)
        throw InputError("pfqn_perm: A must have sum(m) rows");
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    if (n == 0) return one;

    T val = zero;
    std::vector<int> f(R, 0);
    bool more = true;
    while (more) {
        int fs = 0;
        for (int v : f) fs += v;
        T term = (fs % 2 == 0) ? one : -one;
        for (std::size_t k = 0; k < R; ++k) term *= num_nck<T>(m[k], f[k]);
        for (std::size_t i = 0; i < static_cast<std::size_t>(n); ++i) {
            T s = zero;
            for (std::size_t k = 0; k < R; ++k)
                if (f[k] != 0) s += num_traits<T>::from_int(f[k]) * A(i, k);
            term *= s;
        }
        val += term;
        more = next_pop(f, m);
    }
    return (n % 2 == 0) ? val : T(-val);
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_PERM_H
