/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_NCJD_H
#define LINE_API_PFQN_NCJD_H

/**
 * Joint-dependent name of pfqn_ncoi: the balance-function convolution of a
 * closed network whose station rates read the whole per-class occupancy vector.
 *
 * The two names denote the SAME routine because the balanced-fairness recursion
 * Phi_i(0) = 1, mu_i(n) Phi_i(n) = sum_{r: n_r>0} v_{i,r} Phi_i(n - e_r) never
 * inspects the structure of mu_i: it evaluates the callable at the full count
 * vector n. Order independence (mu_i constant on each support) is a modelling
 * restriction that buys insensitivity and a physical reading of Phi, not
 * something the convolution uses, so any joint-dependent scaling eta_i(n)
 * (sn.jdscaling) is admissible, with the proviso that the product form it
 * induces is the balanced-fair one matched to that rate.
 *
 * Use pfqn_ncoi when the model is genuinely order independent and the name
 * should say so, this one when the rate is a general joint dependence.
 * pfqn_clwjd is the transform route and is NOT merely a renaming: it needs the
 * rate to saturate at a finite cutoff.
 */

#include <vector>

#include "line/api/pfqn/pfqn_ncoi.h"

namespace line {
namespace pfqn {

/** @see pfqn_ncoi */
template <class T>
NcResult<T> pfqn_ncjd(const std::vector<T>& Z, const std::vector<int>& N,
                      const std::vector<OiRate<T>>& mu, const Matrix<T>& visits) {
    return pfqn_ncoi(Z, N, mu, visits);
}

/** @see pfqn_ncoi */
template <class T>
NcResult<T> pfqn_ncjd(const std::vector<T>& Z, const std::vector<int>& N,
                      const std::vector<OiRate<T>>& mu) {
    return pfqn_ncoi(Z, N, mu, Matrix<T>());
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_NCJD_H
