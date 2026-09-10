/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_AMVA_COMMON_H
#define LINE_API_PFQN_AMVA_COMMON_H

/**
 * Scaffolding shared by the approximate-MVA family. This header is NOT a port
 * of a MATLAB function: it collects the small utilities that
 * matlab/src/util/ provides to every AMVA routine (oner, pprod, sprod,
 * multichoose, enorm) plus the scheduling-discipline tag those routines branch
 * on, so that each pfqn_*.h below stays a 1:1 port of its own MATLAB file.
 *
 * Scheduling disciplines. The MATLAB routines take sn.sched-valued vectors but
 * only ever compare them against SchedStrategy.PS, SchedStrategy.FCFS and
 * SchedStrategy.INF, so the full enum is not part of the porting surface: the
 * three tags below are what the algorithms actually distinguish. Anything else
 * in a model is handled upstream, in the NetworkStruct layer, which is out of
 * scope here.
 */

#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"
#include "line/util/population.h"

namespace line {
namespace pfqn {

/** The three scheduling disciplines the AMVA and Schmidt recursions branch on. */
enum class SchedStrategy { PS, FCFS, INF };

/**
 * matlab/src/util/oner.m: decrement position r of N, with r given 1-based and
 * r == 0 meaning "leave N alone" (the s == 0 arm of every `for s=0:R` loop).
 * The result may go negative, exactly as in MATLAB; callers that cannot cope
 * with a negative population guard it themselves, as the MATLAB ones do.
 */
inline std::vector<int> oner(const std::vector<int>& N, std::size_t r) {
    std::vector<int> n(N);
    if (r >= 1) {
        if (r > n.size()) throw InputError("oner: class index out of range");
        n[r - 1] -= 1;
    }
    return n;
}

/** matlab/src/util/enorm.m: Frobenius norm of a matrix. */
template <class T>
T enorm(const Matrix<T>& A) {
    T s = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) s += A(i, j) * A(i, j);
    using std::sqrt;
    return sqrt(s);
}

/**
 * Frobenius norm of the difference of two equally shaped matrices, AS A DOUBLE.
 *
 * The only use of this quantity anywhere in the AMVA family is the stopping
 * test `enorm_diff(Q, Qlast) < tol`, and tol is a double. Returning a double
 * therefore loses nothing and buys the exact backend the whole family: sqrt is
 * not an operation of the rational field, so a T-valued version cannot be
 * instantiated there at all. The sum of squares is accumulated in T, exactly,
 * and only the final square root drops to double, so for T == double the
 * result is bit-identical to sqrt of the T-valued sum.
 */
template <class T>
double enorm_diff(const Matrix<T>& A, const Matrix<T>& B) {
    T s = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) {
            const T d = A(i, j) - B(i, j);
            s += d * d;
        }
    return std::sqrt(num_traits<T>::to_double(s));
}

/**
 * Sum the rows of a think-time matrix into a length-R vector, the `sum(Z,1)`
 * that every AMVA entry point performs on its Z argument. An empty Z gives
 * the all-zero vector.
 */
template <class T>
std::vector<T> sum_rows(const Matrix<T>& Z, std::size_t R) {
    std::vector<T> z(R, num_traits<T>::from_int(0));
    if (Z.empty()) return z;
    if (Z.cols() != R) throw InputError("think-time matrix has the wrong class count");
    for (std::size_t k = 0; k < Z.rows(); ++k)
        for (std::size_t r = 0; r < R; ++r) z[r] += Z(k, r);
    return z;
}

/**
 * matlab/src/util/multichoose.m and sprod.m, as an in-place odometer: the
 * compositions of `c` into R non-negative parts, i.e. the vectors n with
 * sum(n) == c. Start from first_composition and iterate until this returns
 * false. The enumeration order differs from MATLAB's recursive one; every use
 * site sums over the whole set, so only completeness matters.
 */
inline void first_composition(std::vector<int>& n, int c) {
    n.assign(n.size(), 0);
    if (!n.empty()) n[0] = c;
}

inline bool next_composition(std::vector<int>& n) {
    const long R = static_cast<long>(n.size());
    if (R < 2) return false;
    long i = R - 2;
    while (i >= 0 && n[i] == 0) --i;
    if (i < 0) return false;
    // Every entry strictly between i and the last is zero by the choice of i,
    // so the whole remainder is carried in the last slot.
    const int tail = n[R - 1];
    n[R - 1] = 0;
    n[i] -= 1;
    n[i + 1] = tail + 1;
    return true;
}

/** Multinomial coefficient sum(m)!/prod_i m_i!, exact in any arithmetic. */
template <class T>
T num_multinomial(const std::vector<int>& m) {
    int tot = 0;
    for (int v : m) tot += v;
    T r = num_factorial<T>(static_cast<unsigned>(tot));
    for (int v : m) r /= num_factorial<T>(static_cast<unsigned>(v));
    return r;
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_AMVA_COMMON_H
