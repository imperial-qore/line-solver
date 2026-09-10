/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_MVAJD_H
#define LINE_API_PFQN_MVAJD_H

/**
 * Joint-dependent name of pfqn_mvaoi: the mean-value analysis of a closed
 * network whose station rates read the whole per-class occupancy vector.
 *
 * The two names denote the SAME routine because the recursion evaluates the rate
 * callable at a full occupancy vector, mu_i(s_i + e_r) with s_i the shift (the
 * occupancy already committed at the bottom of station i), and never inspects
 * the structure of mu_i. That is the "third form" of the Conditional MVA of
 * Casale (QUESTA 2009), a rate depending on the full per-class occupancy vector,
 * so any joint-dependent scaling eta_i(n) (sn.jdscaling) is admissible.
 *
 * Unlike the AMVA joint-dependence route, which evaluates eta at the MEAN
 * arrival-instant vector 1 + E[Q] and therefore collapses a support indicator to
 * 1, this routine evaluates the rate at exact integer occupancies and is exact
 * for the balanced-fair station, at the cost of walking
 * prod_r C(N_r+K+1,K+1) states with K joint-dependent stations.
 */

#include <functional>
#include <vector>

#include "line/api/pfqn/pfqn_mvaoi.h"

namespace line {
namespace pfqn {

/** @see pfqn_mvaoi */
template <class T>
MvaoiResult<T> pfqn_mvajd(const std::vector<T>& Z, const std::vector<int>& N,
                          const std::vector<std::function<T(const std::vector<int>&)>>& mu,
                          const Matrix<T>& Dli, const Matrix<T>& visits, bool want_soi) {
    return pfqn_mvaoi(Z, N, mu, Dli, visits, want_soi);
}

/** @see pfqn_mvaoi */
template <class T>
MvaoiResult<T> pfqn_mvajd(const std::vector<T>& Z, const std::vector<int>& N,
                          const std::vector<std::function<T(const std::vector<int>&)>>& mu,
                          const Matrix<T>& Dli, const Matrix<T>& visits) {
    return pfqn_mvaoi(Z, N, mu, Dli, visits);
}

/** @see pfqn_mvaoi */
template <class T>
MvaoiResult<T> pfqn_mvajd(const std::vector<T>& Z, const std::vector<int>& N,
                          const std::vector<std::function<T(const std::vector<int>&)>>& mu,
                          const Matrix<T>& Dli, bool want_soi) {
    return pfqn_mvaoi(Z, N, mu, Dli, want_soi);
}

/** @see pfqn_mvaoi */
template <class T>
MvaoiResult<T> pfqn_mvajd(const std::vector<T>& Z, const std::vector<int>& N,
                          const std::vector<std::function<T(const std::vector<int>&)>>& mu,
                          const Matrix<T>& Dli) {
    return pfqn_mvaoi(Z, N, mu, Dli);
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_MVAJD_H
