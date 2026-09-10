/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_SCAT_H
#define LINE_API_PFQN_SCAT_H

/**
 * Neuse-Chandy SCAT (Self-Correcting Approximation Technique) approximate MVA.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_scat.m, cross-checked against
 * jar/src/main/java/jline/api/pfqn/mva/Pfqn_scat.java.
 *
 * SCAT shares the Linearizer fixed point: it carries the mean queue lengths at
 * the target population N and at the R reduced populations N - e_s, and
 * corrects the Bard-Schweitzer proportionality assumption with
 *
 *   Delta(i,r,s) = Q(i,r | N - e_s)/(N - e_s)_r - Q(i,r | N)/N_r,
 *
 * held fixed while an inner MVA fixed point is iterated. It differs from
 * Linearizer in that this correction is refreshed ONCE: SCAT stops after the
 * first pass, where Linearizer performs the fixed three passes of Chandy and
 * Neuse (1982), Sec. 4. Cost is therefore about one third of Linearizer's, and
 * accuracy sits between Bard-Schweitzer (the Delta == 0 special case, pfqn_bs)
 * and Linearizer. So this is pfqn_egflinearizer with alpha == 1 and one refresh
 * round, exactly as pfqn_linearizer is the same call with three.
 *
 * SCAT's second departure from Linearizer, fitting a probability mass function
 * centred on the mean queue length at queue-dependent centres instead of
 * propagating the MVA distribution recursion (Krzesinski and Greyling 1984,
 * Sec. 4), does not arise here: this entry point covers single-server and delay
 * stations only, exactly as pfqn_linearizer does. That mass function is
 * available separately as the "scat" marginal rule of pfqn_ab_amva.
 *
 * Arithmetic: inherited from pfqn_egflinearizer. At alpha == 1 the real power
 * degenerates to a rational operation, but the inner Core loop still stops on
 * enorm(Q_{k+1} - Q_k) < tol, so what comes back is the iterate the stopping
 * rule selected rather than the solution of a finite rational problem.
 *
 * Reference: D. Neuse, K. M. Chandy, "SCAT: A Heuristic Algorithm for Queueing
 * Network Models of Computing Systems", ACM SIGMETRICS Perform. Eval. Rev.
 * 10(3), 1981.
 */

#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_amva_common.h"
#include "line/api/pfqn/pfqn_egflinearizer.h"
#include "line/num/number.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/**
 * @param L       (M x R) service demands
 * @param N       (R) population per class
 * @param Z       (K x R) think times, summed over rows; may be empty
 * @param type    (M) scheduling discipline; carried, see pfqn_egflinearizer
 * @param tol     convergence tolerance
 * @param maxiter total inner-iteration budget
 * @param QN0     (M x R) warm start; may be empty
 */
template <class T>
LinearizerResult<T> pfqn_scat(const Matrix<T>& L, const std::vector<int>& N, const Matrix<T>& Z,
                              const std::vector<SchedStrategy>& type, double tol, int maxiter,
                              const Matrix<T>& QN0) {
    const std::vector<T> alpha(N.size(), num_traits<T>::from_int(1));
    // npasses == 1 is what separates SCAT from Linearizer: one Delta refresh, not three
    return pfqn_egflinearizer(L, N, Z, type, tol, maxiter, alpha, QN0, 1);
}

template <class T>
LinearizerResult<T> pfqn_scat(const Matrix<T>& L, const std::vector<int>& N, const Matrix<T>& Z) {
    return pfqn_scat(L, N, Z, std::vector<SchedStrategy>(), 1e-8, 1000, Matrix<T>());
}

template <class T>
LinearizerResult<T> pfqn_scat(const Matrix<T>& L, const std::vector<int>& N) {
    return pfqn_scat(L, N, Matrix<T>(), std::vector<SchedStrategy>(), 1e-8, 1000, Matrix<T>());
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_SCAT_H
