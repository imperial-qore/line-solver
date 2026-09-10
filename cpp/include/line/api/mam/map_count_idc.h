/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MAP_COUNT_IDC_H
#define LINE_API_MAM_MAP_COUNT_IDC_H

/**
 * Index of dispersion for counts (IDC) of a MAP at resolution t.
 *
 * Templated port of matlab/lib/kpctoolbox/map/map_count_idc.m. The IDC of the
 * counting process A(t) is the scaled variance-time curve
 *   I_a(t) = Var(A(t)) / E[A(t)],  t > 0,
 * interpolating between the interarrival SCV at t->0+ and map_idc at t->Inf
 * (Whitt-You, "A Robust Queueing Network Analyzer Based on Indices of
 * Dispersion", eq. 1). It is the traffic descriptor RQNA is built on.
 *
 * ARITHMETIC: transcendental. Var(A(t)) needs the matrix exponential
 * (map_count_var), so this refuses under Rational like its variance input.
 */

#include <cstddef>
#include <vector>

#include "line/api/mam/map_count_mean.h"
#include "line/api/mam/map_count_var.h"
#include "line/api/mam/map_moment.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace mam {

/**
 * @param m the MAP (D0, D1)
 * @param t window lengths (t > 0)
 * @return the IDC at each window length, in the order of t; 1 where the mean is
 *         zero (an orderly point process is locally Poisson as t -> 0)
 */
template <class T>
std::vector<T> map_count_idc(const Map<T>& m, const std::vector<T>& t) {
    static_assert(num_traits<T>::has_transcendental,
                  "map_count_idc requires transcendental arithmetic (matrix exponential)");
    const std::vector<T> mean = map_count_mean(m, t);
    const std::vector<T> var = map_count_var(m, t);
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    std::vector<T> out(t.size(), one);
    for (std::size_t k = 0; k < t.size(); ++k)
        if (mean[k] > zero) out[k] = T(var[k] / mean[k]);
    return out;
}

/** Scalar convenience: the IDC at a single window length t. */
template <class T>
T map_count_idc(const Map<T>& m, const T& t) {
    const std::vector<T> tv(1, t);
    return map_count_idc(m, tv)[0];
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MAP_COUNT_IDC_H
