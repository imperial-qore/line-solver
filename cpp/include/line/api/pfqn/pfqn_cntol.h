/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_CNTOL_H
#define LINE_API_PFQN_CNTOL_H

/**
 * Chandy-Neuse population-scaled termination cutoff for approximate MVA.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_cntol.m. The cutoff
 * 1/(4000 + 16*sum(N)) is published in K. M. Chandy, D. Neuse, "Linearizer: A
 * Heuristic Algorithm for Queuing Network Models of Computing Systems",
 * Commun. ACM 25(2):126-134, 1982, p.129 and appendix. The iteration continues
 * while
 *
 *   max_{i,r} |Q^I(i,r) - Q^{I-1}(i,r)| / N_r > 1/(4000 + 16*|N|),
 *
 * |N| = sum(N). The paper motivates the scaling with |N|: at large populations
 * removing one job changes the queue lengths very little, so a fixed cutoff
 * would terminate the iteration prematurely. It also notes that the expression
 * stays below 0.00025 even at very small populations.
 *
 * The same expression is what LQNS uses as its termination test, set in the
 * SchweitzerCommon constructor of libmva/src/mva.cc; that code carries no
 * citation, and the paper above is its source.
 *
 * The cutoff is a stopping rule, not an algebraic quantity, so it is a double
 * in every arithmetic backend: it is compared against a double residual, and
 * making it exact would not make the fixed point it selects any more exact.
 *
 * Passing NaN as the tol argument of pfqn_bs / pfqn_egflinearizer selects BOTH
 * this cutoff and the normalized-maximum metric of the paper, which is the
 * published test; passing pfqn_cntol(N) as a plain number selects only the
 * cutoff, with those functions' own convergence metric. NaN is the sentinel
 * because it cannot collide with any legitimate tolerance and it is the one
 * form the MATLAB, Java and Python twins share (MATLAB and Python additionally
 * accept the string "cn").
 */

#include <cmath>
#include <vector>

#include "line/num/number.h"

namespace line {
namespace pfqn {

/** Termination cutoff at the given total population. */
inline double pfqn_cntol_total(double total_population) {
    return 1.0 / (4000.0 + 16.0 * total_population);
}

/** Termination cutoff at the given population vector. */
template <class T>
double pfqn_cntol(const std::vector<T>& N) {
    double total = 0.0;
    for (const T& n : N) total += num_traits<T>::to_double(n);
    return pfqn_cntol_total(total);
}

/** Termination cutoff at the given integer population vector. */
inline double pfqn_cntol(const std::vector<int>& N) {
    double total = 0.0;
    for (int n : N) total += static_cast<double>(n);
    return pfqn_cntol_total(total);
}

/** True when tol is the sentinel requesting the Chandy-Neuse test. */
inline bool is_cntol(double tol) { return std::isnan(tol); }

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_CNTOL_H
