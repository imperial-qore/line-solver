/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_EXPAND_H
#define LINE_API_PFQN_EXPAND_H

/**
 * Expand per-station metrics from a reduced model back to the original
 * station set.
 *
 * Templated port of matlab/src/api/pfqn/pfqn_expand.m. The inverse of
 * pfqn_unique: row i of each output is row mapping(i) of the corresponding
 * input, so every station of a replicated group receives the metrics computed
 * once for its representative.
 *
 * Arithmetic: EXACT-CAPABLE. The routine copies values and performs no
 * arithmetic at all, so it is exact in every arithmetic by construction and
 * carries no transcendental gate.
 */

#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace pfqn {

/** Return value of pfqn_expand, mirroring [QN_full,UN_full,CN_full]. */
template <class T>
struct ExpandResult {
    Matrix<T> QN;  ///< (M x R) queue lengths at the original stations
    Matrix<T> UN;  ///< (M x R) utilizations
    Matrix<T> CN;  ///< (M x R) residence times
};

/**
 * @param QN      (M' x R) reduced queue lengths
 * @param UN      (M' x R) reduced utilizations
 * @param CN      (M' x R) reduced residence times
 * @param mapping (M) 0-based unique-station index per original station
 */
template <class T>
ExpandResult<T> pfqn_expand(const Matrix<T>& QN, const Matrix<T>& UN, const Matrix<T>& CN,
                            const std::vector<std::size_t>& mapping) {
    const std::size_t R = QN.cols();
    const std::size_t M = mapping.size();
    if (UN.rows() != QN.rows() || CN.rows() != QN.rows() || UN.cols() != R || CN.cols() != R)
        throw InputError("pfqn_expand: the three metric matrices have different shapes");

    ExpandResult<T> res;
    res.QN = Matrix<T>(M, R);
    res.UN = Matrix<T>(M, R);
    res.CN = Matrix<T>(M, R);
    for (std::size_t i = 0; i < M; ++i) {
        const std::size_t u = mapping[i];
        if (u >= QN.rows()) throw InputError("pfqn_expand: mapping index out of range");
        for (std::size_t j = 0; j < R; ++j) {
            res.QN(i, j) = QN(u, j);
            res.UN(i, j) = UN(u, j);
            res.CN(i, j) = CN(u, j);
        }
    }
    return res;
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_EXPAND_H
