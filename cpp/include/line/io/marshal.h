/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_IO_MARSHAL_H
#define LINE_IO_MARSHAL_H

/**
 * Boundary marshalling for host bindings (MATLAB MEX, pybind11, the JSON CLI).
 *
 * Three things every binding needs and none of them should reinvent:
 *
 * 1. COLUMN-MAJOR INGEST AND EGRESS. MATLAB hands out column-major buffers and
 *    numpy's default is row-major; the port stores row-major. A transpose is
 *    unavoidable in one direction, so it is done once here, explicitly, rather
 *    than by every gateway. Note the honest cost: this is a copy, not a view.
 *    It is O(m n) against algorithms that are O(m n) at best and usually far
 *    worse (the convolution is O(prod(N+1) M R), the LU O(n^3)), so the copy is
 *    not the term that matters. Claiming zero copy here would be a lie that
 *    only pays off for a caller that never runs an algorithm.
 *
 * 2. EXACT VALUES ACROSS A DOUBLE BOUNDARY. Neither MATLAB nor numpy has a
 *    rational type. An exact result therefore crosses as the pair of decimal
 *    integer strings (numerator, denominator) plus a double approximation, so
 *    the host can rebuild it with sym(num)/sym(den) or fractions.Fraction and
 *    can see, rather than guess, that the double is a rounding of something
 *    exact.
 *
 * 3. ERROR IDENTITY. A host needs a stable identifier, not a prose message:
 *    mexErrMsgIdAndTxt takes an id, and a Python binding maps the id to an
 *    exception class. The mapping lives here so all bindings agree.
 */

#include <cstddef>
#include <string>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace io {

/**
 * Build a Matrix from a column-major host buffer (MATLAB's mxGetPr layout,
 * numpy's order='F'). Converts the element type through num_traits, so a
 * double buffer can feed an exact or high-precision instantiation.
 */
template <class T, class S>
Matrix<T> from_column_major(const S* data, std::size_t rows, std::size_t cols) {
    if (rows > 0 && cols > 0 && data == nullptr)
        throw InputError("from_column_major: null buffer with nonzero dimensions");
    Matrix<T> m(rows, cols);
    for (std::size_t j = 0; j < cols; ++j)
        for (std::size_t i = 0; i < rows; ++i)
            m(i, j) = num_traits<T>::from_double(static_cast<double>(data[j * rows + i]));
    return m;
}

/** Build a Matrix from a row-major host buffer (numpy's default order='C'). */
template <class T, class S>
Matrix<T> from_row_major(const S* data, std::size_t rows, std::size_t cols) {
    if (rows > 0 && cols > 0 && data == nullptr)
        throw InputError("from_row_major: null buffer with nonzero dimensions");
    Matrix<T> m(rows, cols);
    for (std::size_t i = 0; i < rows; ++i)
        for (std::size_t j = 0; j < cols; ++j)
            m(i, j) = num_traits<T>::from_double(static_cast<double>(data[i * cols + j]));
    return m;
}

/** Write a Matrix into a column-major host buffer, as doubles. */
template <class T>
void to_column_major(const Matrix<T>& m, double* out) {
    if (!m.empty() && out == nullptr) throw InputError("to_column_major: null destination");
    for (std::size_t j = 0; j < m.cols(); ++j)
        for (std::size_t i = 0; i < m.rows(); ++i)
            out[j * m.rows() + i] = num_traits<T>::to_double(m(i, j));
}

/** Write a vector into a host buffer, as doubles. */
template <class T>
void to_host(const std::vector<T>& v, double* out) {
    if (!v.empty() && out == nullptr) throw InputError("to_host: null destination");
    for (std::size_t i = 0; i < v.size(); ++i) out[i] = num_traits<T>::to_double(v[i]);
}

/** Numerator of an exact value as a decimal string; empty for inexact types. */
template <class T>
std::string exact_numerator(const T&);
template <class T>
std::string exact_denominator(const T&);

/**
 * An exact scalar as it crosses to a host without a rational type: the two
 * decimal integer strings, plus the double a host can use directly.
 */
struct ExactValue {
    std::string numerator;
    std::string denominator;
    double approx = 0.0;
    bool is_exact = false;  ///< false when the source arithmetic was inexact
};

/** Marshal any T; only the exact instantiation fills numerator/denominator. */
template <class T>
ExactValue marshal_scalar(const T& v) {
    ExactValue e;
    e.approx = num_traits<T>::to_double(v);
    e.is_exact = num_traits<T>::is_exact;
    if (num_traits<T>::is_exact) {
        e.numerator = exact_numerator(v);
        e.denominator = exact_denominator(v);
    }
    return e;
}

template <class T>
std::string exact_numerator(const T&) {
    return std::string();
}

template <>
inline std::string exact_numerator<Rational>(const Rational& v) {
    return num_traits<Rational>::numerator_str(v);
}

template <class T>
std::string exact_denominator(const T&) {
    return std::string();
}

template <>
inline std::string exact_denominator<Rational>(const Rational& v) {
    return num_traits<Rational>::denominator_str(v);
}

/**
 * Stable identifier for an error, for mexErrMsgIdAndTxt and for mapping to a
 * host exception class. The prose message stays in what(); hosts must key on
 * the identifier, which will not change with wording.
 */
inline const char* error_id(const Error& e) {
    if (dynamic_cast<const InputError*>(&e) != nullptr) return "line:input";
    if (dynamic_cast<const NumericError*>(&e) != nullptr) return "line:numeric";
    if (dynamic_cast<const UnsupportedError*>(&e) != nullptr) return "line:unsupported";
    return "line:error";
}

}  // namespace io
}  // namespace line

#endif  // LINE_IO_MARSHAL_H
