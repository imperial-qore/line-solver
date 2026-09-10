/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_INFER_INFER_GET_QLEN_ARRIVAL_H
#define LINE_API_INFER_INFER_GET_QLEN_ARRIVAL_H

/**
 * Per-class queue lengths at arrival for the per-class sample format.
 *
 * Templated port of matlab/src/api/infer/infer_get_qlen_arrival.m. No JAR
 * counterpart.
 *
 * MATLAB reads its 6 x (K+1) cell array, taking data{3,k} as the arrival times
 * of class k IN MILLISECONDS and data{4,k} as the response times in seconds.
 * The cell container is a MATLAB storage detail with no meaning in C++, so the
 * port takes the two per-class sample vectors directly; the millisecond to
 * second conversion of the arrival times is kept, since it is what makes the
 * two vectors commensurable.
 *
 * The classes are concatenated in order, given identity job ids, handed to
 * infer_compute_ql_at_arrival, and split back per class.
 *
 * ARITHMETIC: one division and the event replay of
 * infer_compute_ql_at_arrival, so a finite field computation. Note that the
 * division by 1000 is exact in the exact instantiation and only there.
 */

#include <cstddef>
#include <vector>

#include "line/api/infer/infer_compute_ql_at_arrival.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace infer {

/**
 * @param at_ms (K) per-class arrival times, in milliseconds
 * @param rt    (K) per-class response times, in seconds
 * @return      (K) per-class (n_k x K) queue lengths at arrival
 */
template <class T>
std::vector<Matrix<T>> infer_get_qlen_arrival(const std::vector<std::vector<T>>& at_ms,
                                              const std::vector<std::vector<T>>& rt) {
    const std::size_t K = at_ms.size();
    if (rt.size() != K)
        throw InputError("infer_get_qlen_arrival: arrival and response time sets differ in size");

    const T thousand = num_traits<T>::from_int(1000);
    std::vector<T> at, rtall;
    std::vector<std::size_t> cls;
    std::vector<std::size_t> nobs(K, 0);
    for (std::size_t k = 0; k < K; ++k) {
        if (at_ms[k].size() != rt[k].size())
            throw InputError("infer_get_qlen_arrival: a class has mismatched sample counts");
        nobs[k] = at_ms[k].size();
        for (std::size_t i = 0; i < nobs[k]; ++i) {
            at.push_back(at_ms[k][i] / thousand);  // ms -> s
            rtall.push_back(rt[k][i]);
            cls.push_back(k);
        }
    }

    const std::size_t n = at.size();
    std::vector<long> jobid(n);
    for (std::size_t i = 0; i < n; ++i) jobid[i] = static_cast<long>(i) + 1;

    const Matrix<T> ql = infer_compute_ql_at_arrival(at, jobid, rtall, jobid, cls, K);

    std::vector<Matrix<T>> out;
    out.reserve(K);
    std::size_t counter = 0;
    for (std::size_t k = 0; k < K; ++k) {
        Matrix<T> qk(nobs[k], K, num_traits<T>::from_int(0));
        for (std::size_t i = 0; i < nobs[k]; ++i)
            for (std::size_t c = 0; c < K; ++c) qk(i, c) = ql(counter + i, c);
        counter += nobs[k];
        out.push_back(qk);
    }
    return out;
}

}  // namespace infer
}  // namespace line

#endif  // LINE_API_INFER_INFER_GET_QLEN_ARRIVAL_H
