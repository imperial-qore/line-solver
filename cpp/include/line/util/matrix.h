/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_UTIL_MATRIX_H
#define LINE_UTIL_MATRIX_H

/**
 * Dense matrix and non-owning view.
 *
 * API functions take MatrixView<const T> for inputs and return Matrix<T>, so a
 * caller that already holds the data (a numpy array through pybind11, a JSON
 * buffer, another algorithm's output) passes it without a copy. Storage is
 * row-major with an explicit row stride so a view can also address a
 * submatrix or a column-major buffer transposed.
 */

#include <cstddef>
#include <initializer_list>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"

namespace line {

template <class T>
class MatrixView {
public:
    MatrixView() : data_(nullptr), rows_(0), cols_(0), stride_(0) {}
    MatrixView(T* data, std::size_t rows, std::size_t cols)
        : data_(data), rows_(rows), cols_(cols), stride_(cols) {}
    MatrixView(T* data, std::size_t rows, std::size_t cols, std::size_t stride)
        : data_(data), rows_(rows), cols_(cols), stride_(stride) {}

    std::size_t rows() const { return rows_; }
    std::size_t cols() const { return cols_; }
    std::size_t stride() const { return stride_; }
    bool empty() const { return rows_ == 0 || cols_ == 0; }
    T* data() const { return data_; }

    T& operator()(std::size_t i, std::size_t j) const { return data_[i * stride_ + j]; }
    T& at(std::size_t i, std::size_t j) const {
        if (i >= rows_ || j >= cols_) throw InputError("MatrixView index out of range");
        return data_[i * stride_ + j];
    }

private:
    T* data_;
    std::size_t rows_, cols_, stride_;
};

template <class T>
class Matrix {
public:
    Matrix() : rows_(0), cols_(0) {}
    Matrix(std::size_t rows, std::size_t cols) : rows_(rows), cols_(cols), v_(rows * cols) {}
    Matrix(std::size_t rows, std::size_t cols, const T& fill)
        : rows_(rows), cols_(cols), v_(rows * cols, fill) {}

    /** Row-major initializer, e.g. Matrix<double>{{1,2},{3,4}}. */
    Matrix(std::initializer_list<std::initializer_list<T>> rows) : rows_(rows.size()), cols_(0) {
        for (const auto& r : rows) cols_ = r.size() > cols_ ? r.size() : cols_;
        v_.assign(rows_ * cols_, T());
        std::size_t i = 0;
        for (const auto& r : rows) {
            std::size_t j = 0;
            for (const auto& x : r) v_[i * cols_ + (j++)] = x;
            ++i;
        }
    }

    static Matrix row(std::initializer_list<T> vals) {
        Matrix m(1, vals.size());
        std::size_t j = 0;
        for (const auto& x : vals) m(0, j++) = x;
        return m;
    }

    static Matrix col(std::initializer_list<T> vals) {
        Matrix m(vals.size(), 1);
        std::size_t i = 0;
        for (const auto& x : vals) m(i++, 0) = x;
        return m;
    }

    std::size_t rows() const { return rows_; }
    std::size_t cols() const { return cols_; }
    std::size_t size() const { return v_.size(); }
    bool empty() const { return v_.empty(); }

    T* data() { return v_.data(); }
    const T* data() const { return v_.data(); }

    T& operator()(std::size_t i, std::size_t j) { return v_[i * cols_ + j]; }
    const T& operator()(std::size_t i, std::size_t j) const { return v_[i * cols_ + j]; }

    /** Linear access, for vectors held as 1 x n or n x 1. */
    T& operator[](std::size_t k) { return v_[k]; }
    const T& operator[](std::size_t k) const { return v_[k]; }

    MatrixView<T> view() { return MatrixView<T>(v_.data(), rows_, cols_); }
    MatrixView<const T> view() const { return MatrixView<const T>(v_.data(), rows_, cols_); }
    operator MatrixView<const T>() const { return view(); }

    void fill(const T& x) { v_.assign(v_.size(), x); }

    Matrix transpose() const {
        Matrix r(cols_, rows_);
        for (std::size_t i = 0; i < rows_; ++i)
            for (std::size_t j = 0; j < cols_; ++j) r(j, i) = (*this)(i, j);
        return r;
    }

    T sum() const {
        T s = T();
        for (const auto& x : v_) s += x;
        return s;
    }

private:
    std::size_t rows_, cols_;
    std::vector<T> v_;
};

/** Deep copy of a view into an owning matrix, converting the element type. */
template <class T, class S>
Matrix<T> matrix_from(const MatrixView<S>& v) {
    Matrix<T> m(v.rows(), v.cols());
    for (std::size_t i = 0; i < v.rows(); ++i)
        for (std::size_t j = 0; j < v.cols(); ++j) m(i, j) = num_traits<T>::from_double(double(v(i, j)));
    return m;
}

}  // namespace line

#endif  // LINE_UTIL_MATRIX_H
