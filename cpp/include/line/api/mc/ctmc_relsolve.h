/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MC_CTMC_RELSOLVE_H
#define LINE_API_MC_CTMC_RELSOLVE_H

/**
 * Equilibrium distribution relative to a reference state.
 *
 * Templated port of matlab/lib/kpctoolbox/mc/ctmc_relsolve.m. The balance
 * equations are closed with p(refstate) = 1 instead of sum(p) = 1, so the
 * result is NOT a probability vector: it is the stationary measure scaled so
 * that the reference state carries weight one. Dividing by its sum recovers
 * ctmc_solve. Keeping the unnormalized form is what makes the ratios usable
 * when the normalizing constant itself overflows.
 *
 * On a reducible generator the reference-state closure is meaningless across
 * components, so the reducible case falls back to solving each weakly connected
 * component and renormalizing globally, exactly as MATLAB does.
 */

#include <cstddef>
#include <vector>

#include "line/api/mc/ctmc_solve.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace mc {

/** Stationary measure scaled so that entry refstate equals one. */
template <class T>
std::vector<T> ctmc_relsolve(const Matrix<T>& Qin, std::size_t refstate) {
    const std::size_t n = Qin.rows();
    if (Qin.cols() != n) throw InputError("ctmc_relsolve: Q is not square");
    if (n == 0) throw InputError("ctmc_relsolve: Q is empty");
    if (refstate >= n) throw InputError("ctmc_relsolve: refstate out of range");
    const T zero = num_traits<T>::from_int(0);
    if (n == 1) return std::vector<T>(1, num_traits<T>::from_int(1));
    const Matrix<T> Q = ctmc_makeinfgen(Qin);
    const std::vector<std::vector<std::size_t>> comps = detail::weak_components(Q);
    if (comps.size() > 1) {
        std::vector<T> p(n, zero);
        for (std::size_t c = 0; c < comps.size(); ++c) {
            const Matrix<T> Qc = ctmc_makeinfgen(detail::submatrix(Q, comps[c]));
            const std::vector<T> pc = ctmc_solve(Qc);
            for (std::size_t k = 0; k < comps[c].size(); ++k) p[comps[c][k]] = pc[k];
        }
        T s = zero;
        for (std::size_t i = 0; i < n; ++i) s += p[i];
        for (std::size_t i = 0; i < n; ++i) p[i] = p[i] / s;
        return p;
    }
    bool allzero = true;
    for (std::size_t i = 0; i < n && allzero; ++i)
        for (std::size_t j = 0; j < n; ++j)
            if (Q(i, j) != zero) {
                allzero = false;
                break;
            }
    if (allzero)
        return std::vector<T>(n, num_traits<T>::from_int(1) /
                                     num_traits<T>::from_int(static_cast<int>(n)));
    // The last balance equation is redundant, so it is replaced by p(refstate) = 1.
    Matrix<T> A(n, n, zero);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j + 1 < n; ++j) A(j, i) = Q(i, j);
    A(n - 1, refstate) = num_traits<T>::from_int(1);
    std::vector<T> b(n, zero);
    b[n - 1] = num_traits<T>::from_int(1);
    return solve(A, b);
}

/** Reference state 1 in the MATLAB numbering, i.e. index 0 here. */
template <class T>
std::vector<T> ctmc_relsolve(const Matrix<T>& Q) {
    return ctmc_relsolve(Q, static_cast<std::size_t>(0));
}

}  // namespace mc
}  // namespace line

#endif
