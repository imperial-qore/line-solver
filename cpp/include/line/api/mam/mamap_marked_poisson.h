/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MAMAP_MARKED_POISSON_H
#define LINE_API_MAM_MAMAP_MARKED_POISSON_H

/**
 * The marked Poisson process every MAMAP fitter falls back to.
 *
 * Its own header only to break an include cycle: `mamap2m_fit.h` dispatches to
 * the sigma fitters in `mamap22_fit_fs.h` and `mamap22_fit_bs.h`, and those in
 * turn need this fallback. Keeping it here lets the dispatcher include the
 * fitters without the fitters including the dispatcher.
 *
 * A one-phase process carries no autocorrelation and no moment beyond the mean,
 * so this is what a fitter returns when the AMAP(2) has collapsed: the class
 * probabilities are the only descriptors left to honour.
 */

#include <cstddef>
#include <vector>

#include "line/api/mam/mmap_lambda.h"
#include "line/num/number.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {
namespace mamapdetail {

/** A marked Poisson process of the given mean and class law. */
template <class T>
Mmap<T> marked_poisson(const T& mean, const std::vector<T>& p) {
    const T one = num_traits<T>::from_int(1);
    Mmap<T> m;
    m.D0 = Matrix<T>(1, 1, T(-one / mean));
    m.D1 = Matrix<T>(1, 1, T(one / mean));
    m.Dc.assign(p.size(), Matrix<T>(1, 1, num_traits<T>::from_int(0)));
    for (std::size_t c = 0; c < p.size(); ++c) m.Dc[c](0, 0) = T(m.D1(0, 0) * p[c]);
    return m;
}

}  // namespace mamapdetail
}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MAMAP_MARKED_POISSON_H
