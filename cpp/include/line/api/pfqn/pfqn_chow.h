/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_CHOW_H
#define LINE_API_PFQN_CHOW_H

/**
 * @file
 * @ingroup api_pfqn
 * JMT-compatible Chow approximate MVA.
 *
 * JMT estimates every class's arrival-instant queue length at a station by
 * the aggregate queue length at the full population. This is the Bard
 * large-customer-population fixed point implemented by pfqn_lcp.
 */

#include <cstddef>
#include <vector>

#include "line/api/pfqn/pfqn_lcp.h"

namespace line {
namespace pfqn {

/** Legacy compatibility selector; JMT-compatible Chow ignores it. */
enum class ChowVariant { Forward, Backward };

/**
 * @brief JMT-compatible Chow approximate MVA.
 *
 * @param L       (M x R) demands
 * @param N       (R) populations
 * @param Z       (R) think times, empty for none
 * @param type    (M) per-station scheduling, empty for all PS
 * @param tol     convergence tolerance
 * @param maxiter iteration cap
 * @param QN0     queue lengths that warm-start the iteration; empty for a cold start
 * @param variant legacy compatibility selector; accepted but ignored
 */
template <class T>
AmvaResult<T> pfqn_chow(const Matrix<T>& L, const std::vector<T>& N, const std::vector<T>& Z,
                        const std::vector<AmvaSched>& type, double tol = 1e-6,
                        std::size_t maxiter = 1000, const Matrix<T>& QN0 = Matrix<T>(),
                        ChowVariant variant = ChowVariant::Forward) {
    (void)variant;
    return pfqn_lcp(L, N, Z, type, tol, maxiter, QN0);
}

template <class T>
AmvaResult<T> pfqn_chow(const Matrix<T>& L, const std::vector<T>& N, const std::vector<T>& Z) {
    return pfqn_chow(L, N, Z, std::vector<AmvaSched>());
}

template <class T>
AmvaResult<T> pfqn_chow(const Matrix<T>& L, const std::vector<T>& N) {
    return pfqn_chow(L, N, std::vector<T>(), std::vector<AmvaSched>());
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_CHOW_H
