/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FJ_RESPT_CLOSED_H
#define LINE_API_FJ_RESPT_CLOSED_H

/**
 * Varki bound on the residence time of a closed fork-join subnetwork.
 *
 * Templated port of matlab/src/api/fj/fj_respt_closed.m.
 *
 *   R_{P_K}(M) <= x [ H_K + A ]
 *
 * with A the mean number of jobs an arriving job finds at the subnetwork. In a
 * closed network made of the parallel subsystem alone every other job is
 * necessarily inside it, so A = M-1 and the bound is tight at K = 2.
 */

#include "line/api/fj/fj_harmonic.h"
#include "line/api/fj/fj_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace fj {

/** [R, exact] of fj_respt_closed. */
template <class T>
struct FJResptClosedResult {
    T R;
    bool exact;
};

/**
 * @param K number of parallel branches, K >= 1
 * @param x mean service time of each branch
 * @param M number of circulating jobs, M >= 1
 * @param A mean queue length seen on arrival
 * @return  the bound and whether it is known to be tight
 */
template <class T>
FJResptClosedResult<T> fj_respt_closed(unsigned K, const T& x, unsigned M, const T& A) {
    detail::require_positive_K(K, "fj_respt_closed");
    const T zero = num_traits<T>::from_int(0);
    if (!(x > zero)) throw InputError("fj_respt_closed: the mean service time must be positive");
    if (M < 1) throw InputError("fj_respt_closed: M must be a positive integer");
    if (A < zero)
        throw InputError("fj_respt_closed: the arrival-instant queue length must be non-negative");
    FJResptClosedResult<T> out;
    out.R = x * (fj_harmonic<T>(K) + A);
    out.exact = false;
    return out;
}

/**
 * The isolated parallel subsystem of Theorem 4.1, where A = M-1.
 *
 * @param K number of parallel branches, K >= 1
 * @param x mean service time of each branch
 * @param M number of circulating jobs, M >= 1
 * @return  the bound, flagged exact at K = 2
 */
template <class T>
FJResptClosedResult<T> fj_respt_closed(unsigned K, const T& x, unsigned M) {
    if (M < 1) throw InputError("fj_respt_closed: M must be a positive integer");
    FJResptClosedResult<T> out =
        fj_respt_closed<T>(K, x, M, num_traits<T>::from_int(static_cast<long>(M) - 1));
    out.exact = (K == 2);
    return out;
}

}  // namespace fj
}  // namespace line

#endif  // LINE_API_FJ_RESPT_CLOSED_H
