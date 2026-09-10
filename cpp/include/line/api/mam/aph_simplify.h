/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_APH_SIMPLIFY_H
#define LINE_API_MAM_APH_SIMPLIFY_H

/**
 * Composition of two matrix-exponential distributions given in (alpha, T) form.
 *
 * Port of `matlab/lib/kpctoolbox/aph/aph_simplify.m`. The three patterns are
 * the three ways two activities can be composed in an LQN activity graph, which
 * is the caller this exists for (`updateMetricsMomentBased`):
 *
 *   Sequence  X1 + X2               the convolution
 *   Parallel  max(X1, X2)           both run, the composite ends with the last
 *   Branch    X1 w.p. p1, X2 w.p. p2  the mixture
 *
 * WHY (alpha, T) AND NOT A `Map`. A MAP pair (D0, D1) fixes the RESTART law as
 * well as the absorption law, and every one of these patterns is a statement
 * about a single passage: what happens after absorption is the caller's, not
 * the composition's. `aph_fit.h` converts in the other direction
 * (`aph_canonical`) precisely once, when a fitted APH has to become a process.
 *
 * THE ALPHA VECTORS ARE SUB-STOCHASTIC ON PURPOSE. `1 - alpha*e` is the mass
 * that skips the phase process entirely, i.e. an atom at zero, and each pattern
 * routes it explicitly -- in the sequence, the mass of X1 that starts already
 * absorbed enters X2's initial vector directly. Renormalizing alpha would drop
 * that atom and shorten the composite.
 */

#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/** A matrix-exponential law in (alpha, S) form: initial vector and subgenerator. */
template <class T>
struct AphPair {
    std::vector<T> alpha;  ///< (n) initial vector, possibly sub-stochastic
    Matrix<T> S;           ///< (n x n) subgenerator

    std::size_t order() const { return alpha.size(); }
};

/** The `pattern` argument of aph_simplify.m, by name. */
enum class AphPattern { Sequence = 1, Parallel = 2, Branch = 3, Loop = 4 };

namespace detail {

/** Exit-rate column of a subgenerator, -S*e, one entry per phase. */
template <class T>
std::vector<T> aph_exit_rates(const Matrix<T>& S) {
    const T zero = num_traits<T>::from_int(0);
    std::vector<T> out(S.rows(), zero);
    for (std::size_t i = 0; i < S.rows(); ++i) {
        T s = zero;
        for (std::size_t j = 0; j < S.cols(); ++j) s -= S(i, j);
        out[i] = s;
    }
    return out;
}

/** The atom at zero, 1 - alpha*e. */
template <class T>
T aph_zero_atom(const std::vector<T>& alpha) {
    T s = num_traits<T>::from_int(1);
    for (const T& a : alpha) s -= a;
    return s;
}

}  // namespace detail

/**
 * Compose two matrix-exponential laws, as aph_simplify.m does.
 *
 * `p1` and `p2` are read only by `Branch`, where they are the branch
 * probabilities; the other patterns take them as the reference does, i.e. they
 * are present in the signature and unused.
 */
template <class T>
AphPair<T> aph_simplify(const AphPair<T>& d1, const AphPair<T>& d2, const T& p1, const T& p2,
                        AphPattern pattern) {
    const T zero = num_traits<T>::from_int(0);
    const std::size_t n1 = d1.order(), n2 = d2.order();
    if (d1.S.rows() != n1 || d1.S.cols() != n1 || d2.S.rows() != n2 || d2.S.cols() != n2)
        throw InputError("aph_simplify: an (alpha, S) pair has mismatched dimensions");

    AphPair<T> out;
    switch (pattern) {
        case AphPattern::Sequence: {
            // alpha = [a1, (1 - a1*e) a2]; the second block is entered directly
            // by the mass of X1 that was already absorbed at time zero.
            const T atom1 = detail::aph_zero_atom(d1.alpha);
            const std::vector<T> exit1 = detail::aph_exit_rates(d1.S);
            out.alpha.assign(n1 + n2, zero);
            for (std::size_t i = 0; i < n1; ++i) out.alpha[i] = d1.alpha[i];
            for (std::size_t j = 0; j < n2; ++j) out.alpha[n1 + j] = T(atom1 * d2.alpha[j]);
            out.S = Matrix<T>(n1 + n2, n1 + n2, zero);
            for (std::size_t i = 0; i < n1; ++i) {
                for (std::size_t j = 0; j < n1; ++j) out.S(i, j) = d1.S(i, j);
                for (std::size_t j = 0; j < n2; ++j)
                    out.S(i, n1 + j) = T(exit1[i] * d2.alpha[j]);
            }
            for (std::size_t i = 0; i < n2; ++i)
                for (std::size_t j = 0; j < n2; ++j) out.S(n1 + i, n1 + j) = d2.S(i, j);
            return out;
        }
        case AphPattern::Parallel: {
            // States: the n1*n2 product block (both still running), then the n1
            // block (X2 already done) and the n2 block (X1 already done).
            const T atom1 = detail::aph_zero_atom(d1.alpha);
            const T atom2 = detail::aph_zero_atom(d2.alpha);
            const std::vector<T> exit1 = detail::aph_exit_rates(d1.S);
            const std::vector<T> exit2 = detail::aph_exit_rates(d2.S);
            const std::size_t np = n1 * n2, n = np + n1 + n2;
            out.alpha.assign(n, zero);
            for (std::size_t i = 0; i < n1; ++i)
                for (std::size_t j = 0; j < n2; ++j) out.alpha[i * n2 + j] = T(d1.alpha[i] * d2.alpha[j]);
            for (std::size_t i = 0; i < n1; ++i) out.alpha[np + i] = T(atom2 * d1.alpha[i]);
            for (std::size_t j = 0; j < n2; ++j) out.alpha[np + n1 + j] = T(atom1 * d2.alpha[j]);
            out.S = Matrix<T>(n, n, zero);
            // kron(T1, I2) + kron(I1, T2) on the product block
            for (std::size_t i = 0; i < n1; ++i)
                for (std::size_t j = 0; j < n2; ++j) {
                    const std::size_t r = i * n2 + j;
                    for (std::size_t k = 0; k < n1; ++k) out.S(r, k * n2 + j) += d1.S(i, k);
                    for (std::size_t k = 0; k < n2; ++k) out.S(r, i * n2 + k) += d2.S(j, k);
                    // kron(I1, -T2*e): X2 absorbs, X1 continues alone
                    out.S(r, np + i) = exit2[j];
                    // kron(-T1*e, I2): X1 absorbs, X2 continues alone
                    out.S(r, np + n1 + j) = exit1[i];
                }
            for (std::size_t i = 0; i < n1; ++i)
                for (std::size_t j = 0; j < n1; ++j) out.S(np + i, np + j) = d1.S(i, j);
            for (std::size_t i = 0; i < n2; ++i)
                for (std::size_t j = 0; j < n2; ++j)
                    out.S(np + n1 + i, np + n1 + j) = d2.S(i, j);
            return out;
        }
        case AphPattern::Branch: {
            out.alpha.assign(n1 + n2, zero);
            for (std::size_t i = 0; i < n1; ++i) out.alpha[i] = T(p1 * d1.alpha[i]);
            for (std::size_t j = 0; j < n2; ++j) out.alpha[n1 + j] = T(p2 * d2.alpha[j]);
            out.S = Matrix<T>(n1 + n2, n1 + n2, zero);
            for (std::size_t i = 0; i < n1; ++i)
                for (std::size_t j = 0; j < n1; ++j) out.S(i, j) = d1.S(i, j);
            for (std::size_t i = 0; i < n2; ++i)
                for (std::size_t j = 0; j < n2; ++j) out.S(n1 + i, n1 + j) = d2.S(i, j);
            return out;
        }
        default:
            break;
    }
    // The reference's fourth arm is commented out in aph_simplify.m and returns
    // nothing at all, so a caller asking for it in MATLAB gets an undefined
    // output rather than a loop composition. Refusing is the same statement.
    throw UnsupportedError(
        "aph_simplify: the loop pattern is not implemented in the reference either "
        "(aph_simplify.m leaves its fourth arm commented out)");
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_APH_SIMPLIFY_H
