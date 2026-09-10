/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MOMENT_MOMENT_TENSOR_H
#define LINE_API_MOMENT_MOMENT_TENSOR_H

/**
 * Joint moment arrays and the mode products used by every joint conversion.
 *
 * Templated port of matlab/src/api/moment/moment_tensorsize.m,
 * moment_tensortrans.m and moment_jointtrans.m.
 *
 * The storage is COLUMN MAJOR, matching the MATLAB linear index exactly, so the
 * stride of dimension l is prod(sz(1:l-1)). Trailing singleton dimensions are
 * dropped on construction, which is the MATLAB convention that
 * moment_tensorsize enforces: a column vector of n+1 entries and an (n+1)x1
 * array of a degenerate second class are the same object, and both are read as
 * the d = 1 case.
 *
 * The mode argument here is 0-BASED, unlike the 1-based MATLAB argument.
 *
 * Reference:
 * A. Heindl and A. van de Liefvoort. Moment conversions for discrete
 * distributions. PMCCS, 2003.
 */

#include <vector>

#include "line/api/moment/moment_housematrix.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace moment {

/** Joint moment array of size (n_1+1)x...x(n_d+1), stored column major. */
template <class T>
struct MomentTensor {
    std::vector<std::size_t> sz;
    std::vector<T> data;

    MomentTensor() : sz(1, 0) {}

    /** A plain moment vector, i.e. the d = 1 case. */
    explicit MomentTensor(const std::vector<T>& v) : sz(1, v.size()), data(v) {}

    /** A zero-filled array of the given extents, with trailing singletons dropped. */
    explicit MomentTensor(const std::vector<std::size_t>& extents) : sz(extents) {
        normalize();
        data.assign(numel(), num_traits<T>::from_int(0));
    }

    std::size_t order() const { return sz.size(); }

    std::size_t numel() const {
        std::size_t n = 1;
        for (std::size_t l = 0; l < sz.size(); ++l) n *= sz[l];
        return n;
    }

    /** Stride of dimension l in the column-major layout. */
    std::size_t stride(std::size_t l) const {
        std::size_t s = 1;
        for (std::size_t k = 0; k < l; ++k) s *= sz[k];
        return s;
    }

    /** Linear index of the multi-index ord, whose entries are the moment orders. */
    std::size_t index(const std::vector<std::size_t>& ord) const {
        if (ord.size() != sz.size()) throw InputError("MomentTensor: wrong multi-index length");
        std::size_t ia = 0, s = 1;
        for (std::size_t l = 0; l < sz.size(); ++l) {
            if (ord[l] >= sz[l]) throw InputError("MomentTensor: multi-index out of range");
            ia += ord[l] * s;
            s *= sz[l];
        }
        return ia;
    }

    const T& at(const std::vector<std::size_t>& ord) const { return data[index(ord)]; }
    T& at(const std::vector<std::size_t>& ord) { return data[index(ord)]; }

    /** Drops the trailing singleton dimensions, keeping at least one. */
    void normalize() {
        while (sz.size() > 1 && sz.back() == 1) sz.pop_back();
        if (sz.empty()) sz.assign(1, 0);
    }
};

/** Extents with the trailing singleton dimensions removed. */
template <class T>
std::vector<std::size_t> moment_tensorsize(const MomentTensor<T>& A) {
    return A.sz;
}

/** Mode product: every fibre of A along dimension mode is replaced by Tm times that fibre. */
template <class T>
MomentTensor<T> moment_tensortrans(const MomentTensor<T>& A, const Matrix<T>& Tm,
                                   std::size_t mode) {
    if (mode >= A.order())
        throw InputError("moment_tensortrans: the mode must be a dimension index of A");
    const std::size_t n = A.sz[mode];
    if (Tm.cols() != n || Tm.rows() != n)
        throw InputError(
            "moment_tensortrans: T must be square with as many columns as the extent of the "
            "transformed dimension");
    const std::size_t inner = A.stride(mode);
    const std::size_t outer = A.numel() / (inner * n);
    MomentTensor<T> B(A.sz);
    std::vector<T> fibre(n);
    for (std::size_t o = 0; o < outer; ++o)
        for (std::size_t i = 0; i < inner; ++i) {
            const std::size_t base = o * inner * n + i;
            for (std::size_t k = 0; k < n; ++k) fibre[k] = A.data[base + k * inner];
            for (std::size_t r = 0; r < n; ++r) {
                T acc = num_traits<T>::from_int(0);
                for (std::size_t k = 0; k < n; ++k) acc += Tm(r, k) * fibre[k];
                B.data[base + r * inner] = acc;
            }
        }
    return B;
}

/** Applies the conversion matrix of one edge along every dimension of A. */
template <class T>
MomentTensor<T> moment_jointtrans(const MomentTensor<T>& A, MomentEdge edge) {
    MomentTensor<T> B = A;
    for (std::size_t mode = 0; mode < A.order(); ++mode)
        B = moment_tensortrans<T>(B, moment_housematrix<T>(edge, static_cast<int>(A.sz[mode]) - 1),
                                  mode);
    return B;
}

/** String overload matching the MATLAB call signature. */
template <class T>
MomentTensor<T> moment_jointtrans(const MomentTensor<T>& A, const std::string& edge) {
    return moment_jointtrans<T>(A, moment_edge_from_string(edge));
}

}  // namespace moment
}  // namespace line

#endif
