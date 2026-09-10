/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MAP_MARK_H
#define LINE_API_MAM_MAP_MARK_H

/**
 * Mark the arrivals of a MAP with class probabilities
 * (matlab/lib/kpctoolbox/map/map_mark.m).
 *
 * The unmarked inter-arrival process is unchanged: Dc = prob(c) * D1, so the
 * class of an arrival is drawn independently of the phase. This is the
 * phase-independent special case of mmap_mark in mmap_lambda.h, which takes a
 * per-phase weight matrix instead.
 *
 * MATLAB warns and renormalizes when the probabilities do not sum to one;
 * this port renormalizes silently but rejects a non-positive total, since a
 * zero total has no meaningful normalization. Pure field arithmetic, exact at
 * Rational: no transcendental gate.
 */

#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/api/mam/mmap_lambda.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/** MMAP with the same inter-arrival process and arrivals marked by prob. */
template <class T>
Mmap<T> map_mark(const Map<T>& m, const std::vector<T>& prob) {
    const T zero = num_traits<T>::from_int(0);
    if (prob.empty()) throw InputError("map_mark: empty marking probabilities");
    T s = zero;
    for (const T& p : prob) {
        if (p < zero) throw InputError("map_mark: negative marking probability");
        s += p;
    }
    if (s == zero) throw InputError("map_mark: marking probabilities sum to zero");
    Mmap<T> out;
    out.D0 = m.D0;
    out.D1 = m.D1;
    for (std::size_t c = 0; c < prob.size(); ++c) {
        Matrix<T> Dc(m.D1.rows(), m.D1.cols());
        const T w = prob[c] / s;
        for (std::size_t i = 0; i < Dc.rows(); ++i)
            for (std::size_t j = 0; j < Dc.cols(); ++j) Dc(i, j) = w * m.D1(i, j);
        out.Dc.push_back(Dc);
    }
    return mmap_normalize(out);
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MAP_MARK_H
