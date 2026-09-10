/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_APH2_ASSEMBLE_H
#define LINE_API_MAM_APH2_ASSEMBLE_H

/**
 * Assemble an APH(2) from its canonical parameters
 * (matlab/lib/m3a/m3a/aph2/aph2_assemble.m).
 *
 * l1 and l2 are the MEAN holding times of the two phases (not rates) and p1
 * is the probability of continuing from phase 1 to phase 2:
 *
 *   D0 = [ -1/l1  p1/l1 ;  0  -1/l2 ]
 *   D1 = [ (1-p1)/l1  0 ;  1/l2  0  ]
 *
 * Pure field arithmetic -- four reciprocals and three products -- so this is
 * exact at Rational and is deliberately NOT gated on transcendental
 * arithmetic. The moment-matching step that produces l1, l2, p1 is the part
 * that needs a square root; see aph2_fitall.h.
 */

#include "line/api/mam/map_moment.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/** APH(2) with phase means l1, l2 and continuation probability p1. */
template <class T>
Map<T> aph2_assemble(const T& l1, const T& l2, const T& p1) {
    const T zero = num_traits<T>::from_int(0);
    const T one = num_traits<T>::from_int(1);
    if (l1 == zero || l2 == zero) throw InputError("aph2_assemble: zero phase mean");
    Map<T> m;
    m.D0 = Matrix<T>(2, 2, zero);
    m.D1 = Matrix<T>(2, 2, zero);
    m.D0(0, 0) = -one / l1;
    m.D0(0, 1) = p1 / l1;
    m.D0(1, 1) = -one / l2;
    m.D1(0, 0) = (one - p1) / l1;
    m.D1(1, 0) = one / l2;
    return m;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_APH2_ASSEMBLE_H
