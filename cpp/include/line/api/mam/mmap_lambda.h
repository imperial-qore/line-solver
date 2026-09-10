/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MMAP_LAMBDA_H
#define LINE_API_MAM_MMAP_LAMBDA_H

/**
 * Marked MAP (MMAP) algebra: per-class rates, class probabilities,
 * superposition, normalization and scaling.
 *
 * Templated port of the M3A/kpctoolbox MMAP primitives (mmap_count_lambda.m,
 * mmap_lambda.m, mmap_pc.m, mmap_super.m, mmap_normalize.m, mmap_scale.m,
 * mmap_mark.m, mmap_hide.m, mmap_isfeasible.m).
 *
 * An MMAP is the tuple (D0, D1, D1^(1), ..., D1^(C)): D1 is the aggregate
 * arrival matrix and the per-class matrices partition it, sum_c D1^(c) = D1.
 * Every quantity here is rational in the entries, so the exact instantiation
 * carries the partition identity exactly; that identity is precisely what
 * marking and splitting operations break when they are wrong, and a rounded
 * check cannot tell a broken partition from an accumulation of error.
 *
 * Superposition uses the Kronecker sum, as in MATLAB's krons: for generators
 * A and B, krons(A,B) = kron(A, I) + kron(I, B), so the phase process of the
 * superposition is the product chain.
 */

#include <cstddef>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/api/mam/map_transform.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/linalg.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/** An MMAP: the underlying MAP plus the per-class arrival matrices. */
template <class T>
struct Mmap {
    Matrix<T> D0;
    Matrix<T> D1;
    std::vector<Matrix<T>> Dc;  ///< per-class matrices, sum_c Dc = D1

    std::size_t order() const { return D0.rows(); }
    std::size_t classes() const { return Dc.size(); }
    Map<T> map() const { return Map<T>{D0, D1}; }
};

/** Kronecker product. */
template <class T>
Matrix<T> kron(const Matrix<T>& A, const Matrix<T>& B) {
    Matrix<T> C(A.rows() * B.rows(), A.cols() * B.cols(), num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) {
            if (A(i, j) == num_traits<T>::from_int(0)) continue;
            for (std::size_t k = 0; k < B.rows(); ++k)
                for (std::size_t l = 0; l < B.cols(); ++l)
                    C(i * B.rows() + k, j * B.cols() + l) = A(i, j) * B(k, l);
        }
    return C;
}

/** Kronecker sum, MATLAB's krons: kron(A, I_nb) + kron(I_na, B). */
template <class T>
Matrix<T> krons(const Matrix<T>& A, const Matrix<T>& B) {
    if (A.rows() != A.cols() || B.rows() != B.cols())
        throw InputError("krons: both operands must be square");
    const Matrix<T> left = kron(A, eye<T>(B.rows()));
    const Matrix<T> right = kron(eye<T>(A.rows()), B);
    Matrix<T> C(left.rows(), left.cols());
    for (std::size_t i = 0; i < C.rows(); ++i)
        for (std::size_t j = 0; j < C.cols(); ++j) C(i, j) = left(i, j) + right(i, j);
    return C;
}

/**
 * Superposition of two MMAPs: the phase process is the product chain, and the
 * class list of the result is the concatenation of the two class lists
 * (mmap_super.m, 'default' option).
 */
template <class T>
Mmap<T> mmap_super(const Mmap<T>& a, const Mmap<T>& b) {
    const std::size_t na = a.order(), nb = b.order();
    const T zero = num_traits<T>::from_int(0);
    Mmap<T> s;
    s.D0 = krons(a.D0, b.D0);
    s.D1 = krons(a.D1, b.D1);
    const Matrix<T> zeroA(na, na, zero), zeroB(nb, nb, zero);
    for (std::size_t c = 0; c < a.classes(); ++c) s.Dc.push_back(krons(a.Dc[c], zeroB));
    for (std::size_t c = 0; c < b.classes(); ++c) s.Dc.push_back(krons(zeroA, b.Dc[c]));
    return s;
}

/** Per-class arrival rates, lambda_c = theta D1^(c) e. */
template <class T>
std::vector<T> mmap_count_lambda(const Mmap<T>& m) {
    const std::vector<T> theta = map_prob(m.map());
    std::vector<T> lk(m.classes(), num_traits<T>::from_int(0));
    for (std::size_t c = 0; c < m.classes(); ++c) {
        const std::vector<T> t = vecmul(theta, m.Dc[c]);
        for (const T& v : t) lk[c] += v;
    }
    return lk;
}

/** Alias kept for parity with the MATLAB name. */
template <class T>
std::vector<T> mmap_lambda(const Mmap<T>& m) {
    return mmap_count_lambda(m);
}

/** Class probabilities seen by an arriving job, pc = pie (-D0)^-1 D1^(c) e. */
template <class T>
std::vector<T> mmap_pc(const Mmap<T>& m) {
    Matrix<T> negD0 = m.D0;
    for (std::size_t i = 0; i < negD0.rows(); ++i)
        for (std::size_t j = 0; j < negD0.cols(); ++j) negD0(i, j) = -negD0(i, j);
    const Matrix<T> inv = inverse(negD0);
    const std::vector<T> pie = map_pie(m.map());
    std::vector<T> pc(m.classes(), num_traits<T>::from_int(0));
    for (std::size_t c = 0; c < m.classes(); ++c) {
        const std::vector<T> t = vecmul(vecmul(pie, inv), m.Dc[c]);
        for (const T& v : t) pc[c] += v;
    }
    return pc;
}

/** True when the per-class matrices partition D1 exactly and D0 is a generator. */
template <class T>
bool mmap_isfeasible(const Mmap<T>& m) {
    const T zero = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < m.order(); ++i) {
        for (std::size_t j = 0; j < m.order(); ++j) {
            T s = zero;
            for (std::size_t c = 0; c < m.classes(); ++c) {
                if (m.Dc[c](i, j) < zero) return false;
                s += m.Dc[c](i, j);
            }
            if (s != m.D1(i, j)) return false;
            if (i != j && m.D0(i, j) < zero) return false;
        }
        T row = zero;
        for (std::size_t j = 0; j < m.order(); ++j) row += m.D0(i, j) + m.D1(i, j);
        if (row != zero) return false;
    }
    return true;
}

/**
 * Feasibility of a marked MAP WITHIN A TOLERANCE, the semantics of
 * matlab/lib/m3a/m3a/mmap/mmap_isfeasible.m, whose second argument defaults to
 * 10^-map_feastol = 1e-8: every per-class matrix non-negative up to -tol, the
 * per-class matrices summing to D1 up to tol, and the underlying MAP feasible.
 *
 * The mmap_isfeasible above is the EXACT predicate: it
 * compares sums with == and rejects any negative entry however small. That is
 * the right check for an algebraically assembled MMAP, but no optimizer can
 * pass it -- the marking probabilities sum to one only to the tolerance of the
 * solve -- so the optimization-based fits report feasibility through this tolerant form,
 * which is what the reference actually applies to their output.
 */
template <class T>
bool mmap_isfeasible_tol(const Mmap<T>& m, const T& tol) {
    const T zero = num_traits<T>::from_int(0);
    if (!map_isfeasible(m.map(), tol)) return false;
    for (std::size_t c = 0; c < m.classes(); ++c)
        for (std::size_t i = 0; i < m.order(); ++i)
            for (std::size_t j = 0; j < m.order(); ++j)
                if (m.Dc[c](i, j) < -tol) return false;
    for (std::size_t i = 0; i < m.order(); ++i)
        for (std::size_t j = 0; j < m.order(); ++j) {
            T s = m.D1(i, j);
            for (std::size_t c = 0; c < m.classes(); ++c) s -= m.Dc[c](i, j);
            if (num_abs(T(s)) > tol) return false;
        }
    return true;
}

/**
 * Clamp negative off-diagonal and per-class entries to zero and rebuild D1 and
 * the diagonal of D0 from them (mmap_normalize.m).
 */
template <class T>
Mmap<T> mmap_normalize(const Mmap<T>& in) {
    const T zero = num_traits<T>::from_int(0);
    Mmap<T> m = in;
    const std::size_t n = m.order();
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            if (i != j && m.D0(i, j) < zero) m.D0(i, j) = zero;
    for (std::size_t c = 0; c < m.classes(); ++c)
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j)
                if (m.Dc[c](i, j) < zero) m.Dc[c](i, j) = zero;
    m.D1 = Matrix<T>(n, n, zero);
    for (std::size_t c = 0; c < m.classes(); ++c)
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) m.D1(i, j) += m.Dc[c](i, j);
    for (std::size_t i = 0; i < n; ++i) {
        m.D0(i, i) = zero;
        T s = zero;
        for (std::size_t j = 0; j < n; ++j) s += m.D0(i, j) + m.D1(i, j);
        m.D0(i, i) = -s;
    }
    return m;
}

/** Rescale time so that the mean inter-arrival time becomes M. */
template <class T>
Mmap<T> mmap_scale(const Mmap<T>& in, const T& M) {
    if (M == num_traits<T>::from_int(0)) throw InputError("mmap_scale: zero target mean");
    const T ratio = map_mean(in.map()) / M;
    Mmap<T> m = in;
    for (std::size_t i = 0; i < m.order(); ++i)
        for (std::size_t j = 0; j < m.order(); ++j) {
            m.D0(i, j) *= ratio;
            m.D1(i, j) *= ratio;
        }
    for (std::size_t c = 0; c < m.classes(); ++c)
        for (std::size_t i = 0; i < m.order(); ++i)
            for (std::size_t j = 0; j < m.order(); ++j) m.Dc[c](i, j) *= ratio;
    return m;
}

/**
 * Hide a subset of the marks (mmap_hide.m).
 *
 * The process is UNCHANGED -- D0 and the total D1 still describe the same point
 * process -- and only the observation of the hidden classes is removed, which is
 * why the reference renormalizes afterwards: with Dc zeroed for the hidden
 * classes, mmap_normalize rebuilds D1 as the sum of the SURVIVING marks and
 * moves the hidden arrivals into D0 as phase changes. The result is the MAP of
 * the visible class alone, embedded in the joint phase process.
 *
 * @param hide 0-based class indices to hide
 */
template <class T>
Mmap<T> mmap_hide(const Mmap<T>& in, const std::vector<std::size_t>& hide) {
    const T zero = num_traits<T>::from_int(0);
    Mmap<T> m = in;
    for (std::size_t k : hide) {
        if (k >= m.classes()) throw InputError("mmap_hide: class index out of range");
        m.Dc[k] = Matrix<T>(m.order(), m.order(), zero);
    }
    return mmap_normalize(m);
}

/** `mmap_hide(m, setdiff(1:K, keep))`: keep ONE mark, hide every other. */
template <class T>
Mmap<T> mmap_hide_but(const Mmap<T>& in, std::size_t keep) {
    std::vector<std::size_t> hide;
    for (std::size_t k = 0; k < in.classes(); ++k)
        if (k != keep) hide.push_back(k);
    return mmap_hide(in, hide);
}

/** Turn a MAP into a single-class MMAP (mmap_mark with one class). */
template <class T>
Mmap<T> mmap_mark(const Map<T>& base, const Matrix<T>& weights) {
    Mmap<T> m;
    m.D0 = base.D0;
    m.D1 = base.D1;
    const std::size_t C = weights.cols();
    if (C == 0) throw InputError("mmap_mark: no classes");
    for (std::size_t c = 0; c < C; ++c) {
        Matrix<T> Dc(base.D1.rows(), base.D1.cols());
        for (std::size_t i = 0; i < Dc.rows(); ++i)
            for (std::size_t j = 0; j < Dc.cols(); ++j) Dc(i, j) = base.D1(i, j) * weights(i, c);
        m.Dc.push_back(Dc);
    }
    return m;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MMAP_LAMBDA_H
