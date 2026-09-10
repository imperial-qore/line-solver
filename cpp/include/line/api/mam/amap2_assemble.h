/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_AMAP2_ASSEMBLE_H
#define LINE_API_MAM_AMAP2_ASSEMBLE_H

/**
 * Assemble an AMAP(2) in one of the two canonical forms
 * (matlab/lib/m3a/m3a/amap2/amap2_assemble.m).
 *
 * l1, l2 are the MEAN holding times of the two phases; p1, p2 the two
 * branching probabilities. Form 1 is used for a non-negative autocorrelation
 * decay rate, form 2 for a negative one:
 *
 *   form 1: D0 = [ -1/l1  p1/l1 ; 0 -1/l2 ]  D1 = [ (1-p1)/l1  0      ]
 *                                                 [ (1-p2)/l2  p2/l2  ]
 *   form 2: D0 = [ -1/l1  p1/l1 ; 0 -1/l2 ]  D1 = [ 0          (1-p1)/l1 ]
 *                                                 [ (1-p2)/l2  p2/l2     ]
 *
 * Pure field arithmetic, exact at Rational: no transcendental gate.
 */

#include "line/api/mam/map_moment.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/** AMAP(2) in canonical form 1 (gamma >= 0) or 2 (gamma < 0). */
template <class T>
Map<T> amap2_assemble(const T& l1, const T& l2, const T& p1, const T& p2, int form) {
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    if (l1 == zero || l2 == zero) throw InputError("amap2_assemble: zero phase mean");
    Map<T> m;
    m.D0 = Matrix<T>(2, 2, zero);
    m.D1 = Matrix<T>(2, 2, zero);
    m.D0(0, 0) = -one / l1;
    m.D0(0, 1) = p1 / l1;
    m.D0(1, 1) = -one / l2;
    if (form == 1) {
        m.D1(0, 0) = (one - p1) / l1;
        m.D1(1, 0) = (one - p2) / l2;
        m.D1(1, 1) = p2 / l2;
    } else if (form == 2) {
        m.D1(0, 1) = (one - p1) / l1;
        m.D1(1, 0) = (one - p2) / l2;
        m.D1(1, 1) = p2 / l2;
    } else {
        throw InputError("amap2_assemble: form must be 1 (gamma >= 0) or 2 (gamma < 0)");
    }
    return m;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_AMAP2_ASSEMBLE_H
