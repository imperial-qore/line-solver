/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_FES_MAP_LEVELS_H
#define LINE_API_FES_MAP_LEVELS_H

/**
 * Per-level processes of a load-dependent flow-equivalent server.
 *
 * Templated port of matlab/src/api/fes/fes_map_levels.m, mirrored by the JAR and
 * native Python.
 *
 * A flow-equivalent server is described by one MAP (F0^k, F1^k) per population
 * level k = 1..n. A single MAP is replicated over the levels and scaled by
 * min(k, mi), which reproduces a queue with mi servers and, for mi infinite, a
 * delay station serving at rate k mu. The scaling is exact for exponential
 * service and is the load-dependent rate approximation otherwise.
 *
 * ARITHMETIC: field operations only, exact at T = Rational.
 */

#include <cstddef>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace fes {

/** Replicate a load independent MAP over n levels, scaling level k by min(k, mi). */
template <class T>
std::vector<mam::Map<T>> fes_map_levels(const mam::Map<T>& map, std::size_t n, double mi = 1.0) {
    std::vector<mam::Map<T>> levels;
    levels.reserve(n);
    for (std::size_t k = 1; k <= n; ++k) {
        const double sd = (static_cast<double>(k) < mi) ? static_cast<double>(k) : mi;
        const T s = num_traits<T>::from_int(static_cast<long>(sd));
        mam::Map<T> lev = map;
        for (std::size_t i = 0; i < lev.D0.rows(); ++i)
            for (std::size_t j = 0; j < lev.D0.cols(); ++j) {
                lev.D0(i, j) = s * lev.D0(i, j);
                lev.D1(i, j) = s * lev.D1(i, j);
            }
        levels.push_back(lev);
    }
    return levels;
}

/** Validate a load dependent descriptor and return its first n levels. */
template <class T>
std::vector<mam::Map<T>> fes_map_levels(const std::vector<mam::Map<T>>& levels, std::size_t n) {
    if (levels.size() < n)
        throw InputError("fes_map_levels: the flow-equivalent server has fewer levels than required");
    const std::size_t mf = levels[0].order();
    for (std::size_t k = 0; k < n; ++k)
        if (levels[k].order() != mf)
            throw InputError("fes_map_levels: all levels must have the same number of phases");
    return std::vector<mam::Map<T>>(levels.begin(), levels.begin() + static_cast<long>(n));
}

}  // namespace fes
}  // namespace line

#endif  // LINE_API_FES_MAP_LEVELS_H
