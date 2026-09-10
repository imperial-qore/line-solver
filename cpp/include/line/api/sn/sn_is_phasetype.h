/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SN_SN_IS_PHASETYPE_H
#define LINE_API_SN_SN_IS_PHASETYPE_H

/**
 * Port of matlab/src/api/sn/sn_is_phasetype.m.
 *
 * Whether a (D0, D1, ...) list is a valid phase-type / MAP representation:
 * D0 has non-negative off-diagonals, every further block is non-negative, and
 * the entry vector, when supplied, is non-negative. It does NOT check that the
 * rows sum to zero -- `refreshProcessRepresentations` calls it to decide
 * whether a fitted representation may be USED, and a representation that fails
 * only the row-sum test has a different problem.
 *
 * A representation that cannot be inspected -- an empty list, a shorter one
 * than (D0, D1), a non-square D0, a D0 carrying NaN -- returns TRUE, exactly as
 * the reference does: the predicate reports a DEFECT it can see, and "nothing
 * to see" is not a defect. A block after D0 that carries NaN is skipped for the
 * same reason.
 *
 * ARITHMETIC: field. Sign tests against GlobalConstants::Zero.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/lang/lang_types.h"
#include "line/util/matrix.h"

namespace line {
namespace api {

/**
 * @param maps the (D0, D1, ...) blocks
 * @param pie  the entry vector, empty when there is none to check
 */
template <class T>
bool sn_is_phasetype(const std::vector<Matrix<T>>& maps, const std::vector<T>& pie) {
    const double tol = lang::GlobalConstants::Zero;
    if (maps.size() < 2) return true;
    const Matrix<T>& D0 = maps[0];
    const std::size_t n = D0.rows();
    if (n == 0 || D0.cols() != n) return true;
    for (std::size_t a = 0; a < n; ++a)
        for (std::size_t b = 0; b < n; ++b) {
            const double v = num_traits<T>::to_double(D0(a, b));
            if (std::isnan(v)) return true;  // uninspectable, as the reference has it
        }
    for (std::size_t a = 0; a < n; ++a)
        for (std::size_t b = 0; b < n; ++b) {
            if (a == b) continue;
            if (num_traits<T>::to_double(D0(a, b)) < -tol) return false;
        }
    for (std::size_t k = 1; k < maps.size(); ++k) {
        const Matrix<T>& Dk = maps[k];
        bool hasnan = false;
        for (std::size_t a = 0; a < Dk.rows() && !hasnan; ++a)
            for (std::size_t b = 0; b < Dk.cols(); ++b)
                if (std::isnan(num_traits<T>::to_double(Dk(a, b)))) {
                    hasnan = true;
                    break;
                }
        if (hasnan) continue;
        for (std::size_t a = 0; a < Dk.rows(); ++a)
            for (std::size_t b = 0; b < Dk.cols(); ++b)
                if (num_traits<T>::to_double(Dk(a, b)) < -tol) return false;
    }
    for (std::size_t a = 0; a < pie.size(); ++a) {
        const double v = num_traits<T>::to_double(pie[a]);
        if (std::isnan(v)) return true;
        if (v < -tol) return false;
    }
    return true;
}

template <class T>
bool sn_is_phasetype(const std::vector<Matrix<T>>& maps) {
    return sn_is_phasetype(maps, std::vector<T>());
}

}  // namespace api
}  // namespace line

#endif  // LINE_API_SN_SN_IS_PHASETYPE_H
