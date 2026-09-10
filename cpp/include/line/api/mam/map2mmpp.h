/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MAP2MMPP_H
#define LINE_API_MAM_MAP2MMPP_H

/**
 * Read a MAP as an MMPP (matlab/lib/kpctoolbox/map/map2mmpp.m).
 *
 * Returns the modulating generator Q = D0 + D1 and the arrival-rate matrix
 * LAMBDA = D1. The representation is an MMPP only when D1 is diagonal; MATLAB
 * emits a warning in that case and returns anyway, whereas this port reports
 * the off-diagonal mass through the result struct so the caller decides.
 *
 * Pure additions and copies, exact at Rational: no transcendental gate.
 */

#include "line/api/mam/map_moment.h"
#include "line/num/number.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/** Result of map2mmpp. */
template <class T>
struct Map2mmppResult {
    Matrix<T> Q;        ///< modulating generator D0 + D1
    Matrix<T> LAMBDA;   ///< arrival rates D1
    T offdiag_norm;     ///< max |D1(i,j)| over i != j; zero iff a true MMPP
    bool is_mmpp;       ///< offdiag_norm == 0
};

/** Modulating generator and arrival-rate matrix of a MAP read as an MMPP. */
template <class T>
Map2mmppResult<T> map2mmpp(const Map<T>& m) {
    const T zero = num_traits<T>::from_int(0);
    Map2mmppResult<T> r;
    r.Q = map_infgen(m);
    r.LAMBDA = m.D1;
    r.offdiag_norm = zero;
    for (std::size_t i = 0; i < m.D1.rows(); ++i)
        for (std::size_t j = 0; j < m.D1.cols(); ++j) {
            if (i == j) continue;
            const T a = num_abs(T(m.D1(i, j)));
            if (a > r.offdiag_norm) r.offdiag_norm = a;
        }
    r.is_mmpp = (r.offdiag_norm == zero);
    return r;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MAP2MMPP_H
