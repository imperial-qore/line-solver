/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MAP_MOMENT_EXTRA_H
#define LINE_API_MAM_MAP_MOMENT_EXTRA_H

/**
 * Three small MAP quantities the C++ tree had not carried: the factorial and
 * joint moments, and the MMAP generator.
 *
 * Port of `map_factorial_moment`, `map_joint_moment` (from
 * python/line_solver/api/mam/map_derivatives.py) and `mmap_infgen` (from
 * mmap_ops.py). PYTHON-ONLY: no MATLAB or JAR twin.
 *
 * ALL THREE REST ON ONE IDENTITY. The embedded chain of a MAP moves from one
 * arrival to the next with kernel `(-D0)^-1 D1`, and the time spent doing so
 * has, conditional on the phase, a phase-type law with generator D0. So
 *
 *     E[X^k]           = k! pie (-D0)^-k e,
 *     E[X_n^k X_n+1^l] = k! l! pie (-D0)^-k P (-D0)^-l e,   P = (-D0)^-1 D1,
 *
 * with `pie` the EMBEDDED at-arrivals law and `P` the embedded TRANSITION
 * KERNEL. The kernel is what separates the two intervals; dropping it entirely
 * would give the product of two marginals, i.e. the independent case.
 *
 * THE REFERENCE USES `D1` WHERE THE KERNEL BELONGS, and that is a defect, not a
 * convention. `P = (-D0)^-1 D1` carries one more resolvent than `D1` does, so
 * the reference's joint moment is short by exactly one factor of the mean and
 * is dimensionally wrong -- a moment of two times must scale as time squared.
 * A POISSON PROCESS SEES IT IMMEDIATELY: interarrivals are independent there,
 * so `E[X_n X_n+1]` must be `mean^2`, and at rate 1.5 that is 0.4444 while the
 * reference returns 0.6667, the mean itself. This port uses the kernel. Native
 * Python needs the same correction; it is recorded rather than made here,
 * because it changes a reference.
 *
 * `(-D0)^-1` IS APPLIED BY SOLVING, NOT BY INVERTING. The reference forms the
 * inverse explicitly and falls back to a pseudo-inverse when D0 is singular; a
 * singular D0 is a MAP whose phase process cannot leave some state, which is a
 * malformed input rather than something to smooth over, so this port lets the
 * factorization report it.
 *
 * ARITHMETIC: field.
 */

#include <cstddef>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/lu.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

namespace momextradetail {

/** Apply `(-D0)^-1` on the RIGHT of a row vector, i.e. solve `x (-D0) = v`. */
template <class T>
std::vector<T> right_solve_negD0(const Matrix<T>& D0, const std::vector<T>& v) {
    const std::size_t n = D0.rows();
    Matrix<T> A(n, n, num_traits<T>::from_int(0));
    // x (-D0) = v is (-D0)' x' = v', so the transpose is what gets factored.
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) A(i, j) = -D0(j, i);
    return solve(A, v);
}

}  // namespace momextradetail

/** k! pie (-D0)^-k e: the k-th factorial moment of the interarrival time. */
template <class T>
T map_factorial_moment(const Map<T>& m, std::size_t k) {
    const std::size_t n = m.D0.rows();
    if (n == 0 || m.D0.cols() != n) throw InputError("map_factorial_moment: D0 must be square");
    std::vector<T> x = map_pie(m);
    for (std::size_t i = 0; i < k; ++i) x = momextradetail::right_solve_negD0(m.D0, x);
    T s = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < n; ++i) s += x[i];
    T fact = num_traits<T>::from_int(1);
    for (std::size_t i = 2; i <= k; ++i) fact *= num_traits<T>::from_int(static_cast<long>(i));
    return T(fact * s);
}

/**
 * k! l! pie (-D0)^-k P (-D0)^-l e with P = (-D0)^-1 D1: the joint moment of
 * CONSECUTIVE interarrival times.
 *
 * The EMBEDDED KERNEL P sits between the two resolvents, not `D1`; see the
 * header for the measurement that separates the two.
 */
template <class T>
T map_joint_moment(const Map<T>& m, std::size_t k, std::size_t l) {
    const std::size_t n = m.D0.rows();
    if (n == 0 || m.D0.cols() != n) throw InputError("map_joint_moment: D0 must be square");
    if (m.D1.rows() != n || m.D1.cols() != n)
        throw InputError("map_joint_moment: D1 must match D0");

    std::vector<T> x = map_pie(m);
    for (std::size_t i = 0; i < k; ++i) x = momextradetail::right_solve_negD0(m.D0, x);
    // x P, with P = (-D0)^-1 D1: the resolvent FIRST, then the arrival. The
    // reference applies D1 alone and is short by one factor of the mean.
    x = momextradetail::right_solve_negD0(m.D0, x);
    std::vector<T> y(n, num_traits<T>::from_int(0));
    for (std::size_t j = 0; j < n; ++j) {
        T s = num_traits<T>::from_int(0);
        for (std::size_t i = 0; i < n; ++i) s += x[i] * m.D1(i, j);
        y[j] = s;
    }
    for (std::size_t i = 0; i < l; ++i) y = momextradetail::right_solve_negD0(m.D0, y);

    T s = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < n; ++i) s += y[i];
    T fact = num_traits<T>::from_int(1);
    for (std::size_t i = 2; i <= k; ++i) fact *= num_traits<T>::from_int(static_cast<long>(i));
    for (std::size_t i = 2; i <= l; ++i) fact *= num_traits<T>::from_int(static_cast<long>(i));
    return T(fact * s);
}

/** The generator of an MMAP: D0 plus every marked arrival matrix. */
template <class T>
Matrix<T> mmap_infgen(const Matrix<T>& D0, const std::vector<Matrix<T>>& Dk) {
    const std::size_t n = D0.rows();
    if (n == 0 || D0.cols() != n) throw InputError("mmap_infgen: D0 must be square");
    Matrix<T> Q = D0;
    for (std::size_t c = 0; c < Dk.size(); ++c) {
        if (Dk[c].rows() != n || Dk[c].cols() != n)
            throw InputError("mmap_infgen: every arrival matrix must match D0");
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) Q(i, j) += Dk[c](i, j);
    }
    return Q;
}

/**
 * The QBD blocks of a MAP/MAP/1 queue: backward, local and forward.
 *
 * The level is the queue length and the phase is the PAIR (arrival phase,
 * service phase), so every block is a Kronecker product with the identity of
 * the other process:
 *
 *   F = D1_arr (x) I     an arrival raises the level
 *   B = I (x) D1_srv     a completion lowers it
 *   L = D0_arr (x) I + I (x) D0_srv   both processes move, the level does not
 *
 * The ORDER of the factors is the state ordering and cannot be swapped
 * independently in the three: doing so in one alone transposes the phase index
 * and the chain silently describes a different queue.
 */
template <class T>
void qbd_blocks_mapmap1(const Matrix<T>& D0a, const Matrix<T>& D1a, const Matrix<T>& D0s,
                        const Matrix<T>& D1s, Matrix<T>* B, Matrix<T>* L, Matrix<T>* F) {
    const std::size_t na = D0a.rows(), ns = D0s.rows();
    if (na == 0 || ns == 0 || D0a.cols() != na || D1a.rows() != na || D1a.cols() != na ||
        D0s.cols() != ns || D1s.rows() != ns || D1s.cols() != ns)
        throw InputError("qbd_blocks_mapmap1: the two MAPs must be square and self-consistent");
    const T zero = num_traits<T>::from_int(0);
    const std::size_t m = na * ns;

    *F = Matrix<T>(m, m, zero);
    *B = Matrix<T>(m, m, zero);
    *L = Matrix<T>(m, m, zero);
    for (std::size_t i = 0; i < na; ++i)
        for (std::size_t j = 0; j < na; ++j)
            for (std::size_t k = 0; k < ns; ++k) {
                (*F)(i * ns + k, j * ns + k) += D1a(i, j);
                (*L)(i * ns + k, j * ns + k) += D0a(i, j);
            }
    for (std::size_t i = 0; i < na; ++i)
        for (std::size_t k = 0; k < ns; ++k)
            for (std::size_t l = 0; l < ns; ++l) {
                (*B)(i * ns + k, i * ns + l) += D1s(k, l);
                (*L)(i * ns + k, i * ns + l) += D0s(k, l);
            }
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MAP_MOMENT_EXTRA_H
