/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MAP_PDF_H
#define LINE_API_MAM_MAP_PDF_H

/**
 * Probability density of the inter-arrival time of a MAP.
 *
 * Templated port of matlab/lib/kpctoolbox/map/map_pdf.m. Conditional on the
 * phase pie seen by an arrival the inter-arrival time is phase-type with
 * representation (pie, D0), so f(t) = pie exp(D0 t) (-D0) e. The JAR has no
 * counterpart of this function.
 *
 * ARITHMETIC: exp(D0 t) is a tolerance-controlled approximation, so this
 * requires transcendental arithmetic; pie itself is exact.
 */

#include <cstddef>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/expm.h"
#include "line/util/linalg.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/**
 * Probability density of the inter-arrival time at the given points.
 *
 * @param m    the MAP (D0, D1)
 * @param tset evaluation times, each >= 0
 * @return f(t) in the order of tset
 */
template <class T>
std::vector<T> map_pdf(const Map<T>& m, const std::vector<T>& tset) {
    static_assert(num_traits<T>::has_transcendental, "map_pdf requires transcendental arithmetic");
    const std::vector<T> pie = map_pie(m);
    const T zero = num_traits<T>::from_int(0);
    Matrix<T> negD0 = m.D0;
    for (std::size_t i = 0; i < negD0.rows(); ++i)
        for (std::size_t j = 0; j < negD0.cols(); ++j) negD0(i, j) = -negD0(i, j);
    const std::vector<T> e = ones<T>(m.order());
    const std::vector<T> negD0e = mulvec(negD0, e);

    std::vector<T> out;
    out.reserve(tset.size());
    for (std::size_t k = 0; k < tset.size(); ++k) {
        if (tset[k] < zero) throw InputError("map_pdf: negative evaluation point");
        std::vector<T> v = pie;
        if (!(tset[k] == zero)) v = vecmul(pie, expm(m.D0, tset[k]));
        T s = zero;
        for (std::size_t i = 0; i < v.size(); ++i) s += v[i] * negD0e[i];
        out.push_back(s);
    }
    return out;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MAP_PDF_H
