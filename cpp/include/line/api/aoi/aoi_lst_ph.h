/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_AOI_LST_PH_H
#define LINE_API_AOI_LST_PH_H

/**
 * Laplace-Stieltjes transform of a phase-type distribution PH(alpha, T).
 *
 * Templated port of matlab/src/api/aoi/aoi_lst_ph.m, cross-checked against
 * Aoi_lst.ph in jar/src/main/java/jline/api/aoi/Aoi_lst.java (identical).
 *
 *   H*(s) = alpha (s I - T)^{-1} t,   t = -T e
 *
 * MATLAB evaluates it as the linear solve alpha * ((sI - T) \ t) rather than
 * by forming the inverse, and this port does the same through the LU in
 * line/util/lu.h. A linear solve stays in the field, so a PH with rational
 * parameters has an exactly representable transform at every rational s --
 * which makes this the natural exact stand-in for aoi_lst_det, since an
 * Erlang-k with k -> inf approaches a constant while staying rational.
 */

#include <cstddef>
#include <vector>

#include "line/api/aoi/aoi_types.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace aoi {

/**
 * @param alpha initial probability row vector, length n
 * @param Tmat  sub-generator, n x n
 * @return      s -> alpha (sI - Tmat)^{-1} (-Tmat e)
 */
template <class T>
Lst<T> aoi_lst_ph(const std::vector<T>& alpha, const Matrix<T>& Tmat) {
    const std::size_t n = alpha.size();
    if (n == 0) throw InputError("aoi_lst_ph: alpha must be non-empty");
    if (Tmat.rows() != n || Tmat.cols() != n)
        throw InputError("aoi_lst_ph: the sub-generator must be square and match the length of alpha");

    // Exit rate vector t = -T e.
    std::vector<T> t = mulvec(Tmat, ones<T>(n));
    for (std::size_t i = 0; i < n; ++i) t[i] = -t[i];

    return [alpha, Tmat, t, n](const T& s) {
        Matrix<T> A(n, n, num_traits<T>::from_int(0));
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) A(i, j) = (i == j ? T(s - Tmat(i, j)) : T(-Tmat(i, j)));
        std::vector<T> x = t;
        const std::vector<std::size_t> piv = lu_factor(A);
        lu_solve(A, piv, x);
        T v = num_traits<T>::from_int(0);
        for (std::size_t i = 0; i < n; ++i) v += alpha[i] * x[i];
        return v;
    };
}

}  // namespace aoi
}  // namespace line

#endif  // LINE_API_AOI_LST_PH_H
