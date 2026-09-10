/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MC_CTMC_PSEUDOSTOCHCOMP_H
#define LINE_API_MC_CTMC_PSEUDOSTOCHCOMP_H

/**
 * Pseudo stochastic complement of a CTMC partition.
 *
 * Templated port of jar/src/main/java/jline/api/mc/Ctmc_pseudostochcomp.java,
 * which has no MATLAB twin. The exact complement over a retained set I censors
 * the excursions through the complement Ic and needs the inverse of Q22; the
 * pseudo complement replaces that inverse by the single rank-one return law
 *   S = Q11 + Q12 1 y,   y = pi(Ic) Q21 / sum( pi(Ic) Q21 ),
 * i.e. every excursion into Ic is assumed to re-enter I through the stationary
 * re-entry distribution y, independently of where it left. This is exact when
 * Q22 is a single lumped state and is the Takahashi-style approximation
 * otherwise; it costs one stationary solve instead of one linear solve per
 * retained state.
 *
 * The default retained set, used when `keep` is empty, is the first
 * ceil(n/2) states, matching the reference.
 *
 * ARITHMETIC: field plus the stationary solve, so exact under Rational whenever
 * ctmc_solve is.
 */

#include <cstddef>
#include <vector>

#include "line/api/mc/ctmc_solve.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mc {

/** Blocks of the partition together with the pseudo complement over I. */
template <class T>
struct PseudoStochCompResult {
    Matrix<T> S;    ///< pseudo stochastic complement over the retained set
    Matrix<T> Q11;  ///< retained-to-retained block
    Matrix<T> Q12;  ///< retained-to-complement block
    Matrix<T> Q21;  ///< complement-to-retained block
    Matrix<T> Q22;  ///< complement-to-complement block
    Matrix<T> Tm;   ///< the rank-one return term S - Q11
};

/**
 * @param Q    (n x n) generator
 * @param keep 0-based indices of the retained set I; empty selects the first
 *             ceil(n/2) states
 */
template <class T>
PseudoStochCompResult<T> ctmc_pseudostochcomp(const Matrix<T>& Q,
                                              const std::vector<std::size_t>& keep) {
    const std::size_t n = Q.rows();
    if (Q.cols() != n) throw InputError("ctmc_pseudostochcomp: generator is not square");
    const T zero = num_traits<T>::from_int(0);

    std::vector<std::size_t> I = keep;
    if (I.empty())
        for (std::size_t i = 0; i < (n + 1) / 2; ++i) I.push_back(i);

    std::vector<bool> kept(n, false);
    for (std::size_t i = 0; i < I.size(); ++i) {
        if (I[i] >= n) throw InputError("ctmc_pseudostochcomp: a retained index is out of range");
        kept[I[i]] = true;
    }
    std::vector<std::size_t> Ic;
    for (std::size_t i = 0; i < n; ++i)
        if (!kept[i]) Ic.push_back(i);
    if (Ic.empty()) throw InputError("ctmc_pseudostochcomp: the complement set is empty");

    const std::size_t nk = I.size(), nd = Ic.size();
    PseudoStochCompResult<T> r;
    r.Q11 = Matrix<T>(nk, nk, zero);
    r.Q12 = Matrix<T>(nk, nd, zero);
    r.Q21 = Matrix<T>(nd, nk, zero);
    r.Q22 = Matrix<T>(nd, nd, zero);
    for (std::size_t a = 0; a < nk; ++a) {
        for (std::size_t b = 0; b < nk; ++b) r.Q11(a, b) = Q(I[a], I[b]);
        for (std::size_t b = 0; b < nd; ++b) r.Q12(a, b) = Q(I[a], Ic[b]);
    }
    for (std::size_t a = 0; a < nd; ++a) {
        for (std::size_t b = 0; b < nk; ++b) r.Q21(a, b) = Q(Ic[a], I[b]);
        for (std::size_t b = 0; b < nd; ++b) r.Q22(a, b) = Q(Ic[a], Ic[b]);
    }

    const std::vector<T> pie = ctmc_solve(Q);

    // y = pi(Ic) Q21, normalized to a probability over the re-entry states.
    std::vector<T> y(nk, zero);
    T sy = zero;
    for (std::size_t b = 0; b < nk; ++b) {
        T acc = zero;
        for (std::size_t a = 0; a < nd; ++a) acc += pie[Ic[a]] * r.Q21(a, b);
        y[b] = acc;
        sy += acc;
    }
    if (sy == zero) throw NumericError("ctmc_pseudostochcomp: no flow returns to the retained set");
    for (std::size_t b = 0; b < nk; ++b) y[b] = T(y[b] / sy);

    // Q12 1 y: the row sums of Q12 spread over the re-entry law.
    r.Tm = Matrix<T>(nk, nk, zero);
    for (std::size_t a = 0; a < nk; ++a) {
        T out = zero;
        for (std::size_t b = 0; b < nd; ++b) out += r.Q12(a, b);
        for (std::size_t b = 0; b < nk; ++b) r.Tm(a, b) = T(out * y[b]);
    }

    r.S = r.Q11;
    for (std::size_t a = 0; a < nk; ++a)
        for (std::size_t b = 0; b < nk; ++b) r.S(a, b) = T(r.S(a, b) + r.Tm(a, b));
    return r;
}

/** Default partition: the first ceil(n/2) states. */
template <class T>
PseudoStochCompResult<T> ctmc_pseudostochcomp(const Matrix<T>& Q) {
    return ctmc_pseudostochcomp(Q, std::vector<std::size_t>());
}

}  // namespace mc
}  // namespace line

#endif  // LINE_API_MC_CTMC_PSEUDOSTOCHCOMP_H
