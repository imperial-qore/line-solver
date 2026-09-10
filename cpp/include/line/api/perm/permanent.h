/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PERM_PERMANENT_H
#define LINE_API_PERM_PERMANENT_H

/**
 * The PERMANENT of a matrix, exactly, by four algorithms.
 *
 * Port of python/line_solver/api/perm/exact.py. PYTHON-ONLY: no MATLAB or JAR
 * twin, so native Python is the reference.
 *
 * WHY A QUEUEING LIBRARY WANTS PERMANENTS. The normalizing constant of a closed
 * multiclass network with distinguishable jobs is a permanent of the demand
 * matrix with COLUMN MULTIPLICITIES given by the class populations -- see
 * `api/pfqn/pfqn_lcfsqn_nc.h`, which already relies on that identity. The
 * permanent looks like a determinant with every sign made positive, and that
 * single change removes the multilinear cancellation Gaussian elimination
 * lives on, which is why there is no polynomial algorithm and why four of them
 * are carried here rather than one.
 *
 * THE FOUR, and when each is the right choice:
 *
 *  - MULTIPLICITY: inclusion-exclusion over the DISTINCT columns, weighted by
 *    binomial coefficients. Its cost is set by the number of distinct columns,
 *    not by n, so on the matrices this library actually forms -- a station's
 *    demand repeated once per job of a class -- it is the only tractable one.
 *    This is the default for exactly that reason.
 *  - RYSER: the textbook 2^n subset sum, O(2^n n^2). No structure exploited,
 *    and the one to compare the others against.
 *  - RYSER-GRAY: the same sum walked in GRAY-CODE order, so consecutive subsets
 *    differ in one column and the row sums update in O(n) instead of O(n^2).
 *    Same value, one factor of n cheaper.
 *  - NAIVE: every permutation, n!. Unusable past about ten, and kept because it
 *    is the definition -- it is what the others are verified against.
 *
 * ALL FOUR RETURN THE SAME NUMBER, and the tests assert exactly that on random
 * matrices. A permanent has no cheap independent check, so agreement between an
 * O(n!) definition and an O(2^n) formula IS the verification.
 *
 * ARITHMETIC: field. Ryser's alternating sum cancels heavily, so on a large
 * matrix the exact instantiation is not a luxury.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace perm {

/** Which algorithm `permanent` should use. */
enum class PermMethod { Multiplicity = 0, Ryser, RyserGray, Naive };

namespace permdetail {

/** Binomial coefficient, exactly, at the sizes a permanent reaches. */
template <class T>
T binom(std::size_t n, std::size_t k) {
    if (k > n) return num_traits<T>::from_int(0);
    T v = num_traits<T>::from_int(1);
    for (std::size_t i = 0; i < k; ++i)
        v = T(v * num_traits<T>::from_int(static_cast<long>(n - i)) /
              num_traits<T>::from_int(static_cast<long>(i + 1)));
    return v;
}

/** The distinct columns of `m`, and how many times each occurs. */
template <class T>
void unique_columns(const Matrix<T>& m, Matrix<T>* uniq, std::vector<std::size_t>* mult) {
    const std::size_t n = m.rows(), c = m.cols();
    std::vector<std::size_t> rep;   // representative column of each group
    mult->clear();
    for (std::size_t j = 0; j < c; ++j) {
        bool found = false;
        for (std::size_t g = 0; g < rep.size() && !found; ++g) {
            bool same = true;
            for (std::size_t i = 0; i < n && same; ++i)
                if (m(i, j) != m(i, rep[g])) same = false;
            if (same) {
                ++(*mult)[g];
                found = true;
            }
        }
        if (!found) {
            rep.push_back(j);
            mult->push_back(1);
        }
    }
    *uniq = Matrix<T>(n, rep.size(), num_traits<T>::from_int(0));
    for (std::size_t g = 0; g < rep.size(); ++g)
        for (std::size_t i = 0; i < n; ++i) (*uniq)(i, g) = m(i, rep[g]);
}

/** The next vector of the box product 0 <= f_k <= mult_k, or false at the end. */
inline bool pprod_next(std::vector<std::size_t>& f, const std::vector<std::size_t>& mult) {
    for (std::size_t k = 0; k < f.size(); ++k) {
        if (f[k] < mult[k]) {
            ++f[k];
            return true;
        }
        f[k] = 0;
    }
    return false;
}

}  // namespace permdetail

/**
 * Inclusion-exclusion over the DISTINCT columns.
 *
 * With R distinct columns of multiplicities m_1..m_R the sum runs over the box
 * `0 <= f_k <= m_k` rather than over 2^n subsets, so a matrix whose columns
 * repeat -- which is what a class population produces -- costs
 * `prod (m_k + 1)` terms instead of `2^(sum m_k)`.
 */
template <class T>
T permanent_multiplicity(const Matrix<T>& m) {
    const std::size_t n = m.rows();
    if (n == 0) return num_traits<T>::from_int(1);
    if (m.cols() != n) throw InputError("permanent: the matrix must be square");

    Matrix<T> uniq;
    std::vector<std::size_t> mult;
    permdetail::unique_columns(m, &uniq, &mult);
    const std::size_t R = mult.size();

    T value = num_traits<T>::from_int(0);
    std::vector<std::size_t> f(R, 0);
    do {
        std::size_t fsum = 0;
        for (std::size_t k = 0; k < R; ++k) fsum += f[k];
        T term = num_traits<T>::from_int((fsum % 2 == 0) ? 1 : -1);
        for (std::size_t j = 0; j < R; ++j) term *= permdetail::binom<T>(mult[j], f[j]);
        for (std::size_t i = 0; i < n; ++i) {
            T s = num_traits<T>::from_int(0);
            for (std::size_t k = 0; k < R; ++k)
                s += num_traits<T>::from_int(static_cast<long>(f[k])) * uniq(i, k);
            term *= s;
        }
        value += term;
    } while (permdetail::pprod_next(f, mult));
    return T(num_traits<T>::from_int((n % 2 == 0) ? 1 : -1) * value);
}

/** Ryser's formula over explicit column subsets: O(2^n n^2). */
template <class T>
T permanent_ryser(const Matrix<T>& m) {
    const std::size_t n = m.rows();
    if (n == 0) return num_traits<T>::from_int(1);
    if (m.cols() != n) throw InputError("permanent: the matrix must be square");
    if (n > 30) throw InputError("permanent_ryser: 2^n subsets is not enumerable past n = 30");

    T total = num_traits<T>::from_int(0);
    const unsigned long long lim = 1ULL << n;
    for (unsigned long long sub = 0; sub < lim; ++sub) {
        std::size_t bits = 0;
        for (std::size_t j = 0; j < n; ++j)
            if (sub & (1ULL << j)) ++bits;
        T prod = num_traits<T>::from_int(1);
        for (std::size_t i = 0; i < n; ++i) {
            T rs = num_traits<T>::from_int(0);
            for (std::size_t j = 0; j < n; ++j)
                if (sub & (1ULL << j)) rs += m(i, j);
            prod *= rs;
        }
        total += num_traits<T>::from_int(((n - bits) % 2 == 0) ? 1 : -1) * prod;
    }
    return total;
}

/**
 * Ryser's formula in GRAY-CODE order: O(2^n n).
 *
 * Consecutive subsets differ in exactly one column, so the row sums are updated
 * rather than recomputed. That is where the factor of n goes.
 */
template <class T>
T permanent_ryser_gray(const Matrix<T>& m) {
    const std::size_t n = m.rows();
    if (n == 0) return num_traits<T>::from_int(1);
    if (m.cols() != n) throw InputError("permanent: the matrix must be square");
    if (n > 30) throw InputError("permanent_ryser_gray: 2^n steps is not enumerable past n = 30");

    std::vector<T> rowsum(n, num_traits<T>::from_int(0));
    std::vector<char> on(n, 0);
    T value = num_traits<T>::from_int(0);
    std::size_t bits = 0;

    const unsigned long long steps = (1ULL << n) - 1ULL;
    for (unsigned long long s = 1; s <= steps; ++s) {
        // The bit that changes between Gray codes s-1 and s is the index of the
        // lowest set bit of s.
        std::size_t j = 0;
        unsigned long long t = s;
        while ((t & 1ULL) == 0ULL) {
            t >>= 1;
            ++j;
        }
        on[j] = !on[j];
        const T sign = num_traits<T>::from_int(on[j] ? 1 : -1);
        for (std::size_t i = 0; i < n; ++i) rowsum[i] += sign * m(i, j);
        bits = on[j] ? bits + 1 : bits - 1;

        T prod = num_traits<T>::from_int(1);
        for (std::size_t i = 0; i < n; ++i) prod *= rowsum[i];
        value += num_traits<T>::from_int((bits % 2 == 0) ? 1 : -1) * prod;
    }
    return T(num_traits<T>::from_int((n % 2 == 0) ? 1 : -1) * value);
}

/** Every permutation: n!. The definition, and what the rest are checked against. */
template <class T>
T permanent_naive(const Matrix<T>& m) {
    const std::size_t n = m.rows();
    if (n == 0) return num_traits<T>::from_int(1);
    if (m.cols() != n) throw InputError("permanent: the matrix must be square");
    if (n > 12) throw InputError("permanent_naive: n! is not enumerable past n = 12");

    std::vector<std::size_t> p(n);
    for (std::size_t i = 0; i < n; ++i) p[i] = i;
    T value = num_traits<T>::from_int(0);
    do {
        T prod = num_traits<T>::from_int(1);
        for (std::size_t i = 0; i < n; ++i) prod *= m(i, p[i]);
        value += prod;
    } while (std::next_permutation(p.begin(), p.end()));
    return value;
}

/**
 * The permanent, by the chosen method.
 *
 * The default exploits repeated columns, which is the case this library
 * generates: a class of N jobs contributes N identical columns.
 */
template <class T>
T permanent(const Matrix<T>& m, PermMethod method = PermMethod::Multiplicity) {
    switch (method) {
        case PermMethod::Multiplicity: return permanent_multiplicity(m);
        case PermMethod::Ryser: return permanent_ryser(m);
        case PermMethod::RyserGray: return permanent_ryser_gray(m);
        case PermMethod::Naive: return permanent_naive(m);
    }
    throw InputError("permanent: unknown method");
}

/**
 * Round a matrix's entries onto a coarse lattice so repeated columns are found.
 *
 * NOT the reference's `preprocessing_ds`, despite what this comment used to
 * claim: that one is the Sinkhorn scaling to double stochasticity (python
 * `preprocessing_ds`, JAR `QueueingNetwork.preprocessingDS`) and returns a
 * rescaling factor alongside the matrix. This is an unrelated operation that
 * merely shared the name. The multiplicity algorithm keys on EXACT
 * column equality, so two demands that differ in the last bits are two distinct
 * columns and the saving is lost; snapping to a tolerance recovers it. This is
 * a deliberate perturbation of the input, not a numerical tidy-up, so it is a
 * separate call rather than something `permanent` does silently.
 */
template <class T>
Matrix<T> snap_to_lattice(const Matrix<T>& m, double tolerance = 0.001) {
    if (!(tolerance > 0.0)) throw InputError("snap_to_lattice: the tolerance must be positive");
    Matrix<T> out(m.rows(), m.cols(), num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < m.rows(); ++i)
        for (std::size_t j = 0; j < m.cols(); ++j) {
            const double v = num_traits<T>::to_double(m(i, j));
            out(i, j) = num_traits<T>::from_double(std::round(v / tolerance) * tolerance);
        }
    return out;
}

}  // namespace perm
}  // namespace line

#endif  // LINE_API_PERM_PERMANENT_H
