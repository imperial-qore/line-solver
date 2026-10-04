/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MMAP_BACKWARD_MOMENT_H
#define LINE_API_MAM_MMAP_BACKWARD_MOMENT_H

/**
 * @file
 * @ingroup api_mam
 * Class-conditional backward moments of an MMAP
 * (matlab/lib/m3a/m3a/mmap/mmap_backward_moment.m).
 *
 * Its own header, and not part of mmap_compress.h, because the fitting
 * families mmap_compress dispatches to (maph2m, mamap2m, mamap22) take backward
 * moments themselves: had they kept including mmap_compress.h, mmap_compress
 * could not include them back.
 */

#include <cstddef>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/api/mam/map_transform.h"
#include "line/api/mam/mmap_lambda.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/**
 * Class-conditional backward moments of an MMAP (mmap_backward_moment.m).
 *
 * @param orders the moment orders to compute
 * @param normalized true for B(c,k) with M_k = sum_c B(c,k) p_c, i.e. divided
 *        by the class probability p_c (the MATLAB default); false for the
 *        unnormalized form with M_k = sum_c B(c,k)
 * @param m the marked MAP whose backward moments are taken
 * @return B[c][h], the moment of order orders[h] for class c
 */
template <class T>
std::vector<std::vector<T>> mmap_backward_moment(const Mmap<T>& m,
                                                 const std::vector<unsigned>& orders,
                                                 bool normalized) {
    const std::size_t n = m.order();
    const std::size_t C = m.classes();
    const T zero = num_traits<T>::from_int(0);
    const std::vector<T> pie = map_pie(m.map());
    Matrix<T> negD0 = m.D0;
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) negD0(i, j) = -negD0(i, j);
    const Matrix<T> Minv = inverse(negD0);

    std::vector<std::vector<T>> B(C, std::vector<T>(orders.size(), zero));
    for (std::size_t c = 0; c < C; ++c) {
        T pa = num_traits<T>::from_int(1);
        if (normalized) {
            const std::vector<T> t = vecmul(vecmul(pie, Minv), m.Dc[c]);
            pa = zero;
            for (const T& v : t) pa += v;
            if (pa == zero)
                throw NumericError(
                    "mmap_backward_moment: class with zero arrival probability cannot be "
                    "normalized");
        }
        for (std::size_t h = 0; h < orders.size(); ++h) {
            const unsigned k = orders[h];
            const std::vector<T> t = vecmul(vecmul(pie, matpow(Minv, k + 1)), m.Dc[c]);
            T s = zero;
            for (const T& v : t) s += v;
            B[c][h] = num_factorial<T>(k) / pa * s;
        }
    }
    return B;
}

/** mmap_backward_moment with the MATLAB default, normalized. */
template <class T>
std::vector<std::vector<T>> mmap_backward_moment(const Mmap<T>& m,
                                                 const std::vector<unsigned>& orders) {
    return mmap_backward_moment(m, orders, true);
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MMAP_BACKWARD_MOMENT_H
