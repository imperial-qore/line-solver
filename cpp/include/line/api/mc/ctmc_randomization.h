/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MC_CTMC_RANDOMIZATION_H
#define LINE_API_MC_CTMC_RANDOMIZATION_H

/**
 * Uniformization (randomization) of a CTMC: the embedded DTMC P = I + Q/q.
 *
 * Templated port of matlab/lib/kpctoolbox/mc/ctmc_randomization.m and
 * jar/src/main/java/jline/api/mc/Ctmc_randomization.java. The rate q must
 * dominate max_i |q_ii| or P is not stochastic; the result is then passed
 * through dtmc_makestochastic, which repairs the rounding of the division.
 *
 * REFERENCE DEFECT. The MATLAB default rate is max|Q| + rand, drawn from the
 * global unseeded stream, so two calls on the same generator return different
 * matrices and nothing downstream of it is reproducible. The default here is
 * the deterministic q = (21/20) max|Q|, which is the same rate ctmc_courtois
 * derives explicitly and which strictly exceeds max_i |q_ii| (so P stays
 * stochastic and aperiodic). Every quantity the callers in this port take from
 * P -- stationary vectors, SCC structure, aggregation matrices -- is invariant
 * to q, so the substitution changes no result, only its reproducibility.
 *
 * Every operation is a field operation, so this is exact at Rational.
 */

#include <cstddef>

#include "line/api/mc/dtmc_makestochastic.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace mc {

template <class T>
struct RandomizationResult {
    Matrix<T> P;  ///< uniformized stochastic matrix
    T q;          ///< rate actually used
};

/** Largest magnitude of any entry of Q; equals max_i |q_ii| for a generator. */
template <class T>
T ctmc_maxabs(const Matrix<T>& Q) {
    T m = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < Q.rows(); ++i)
        for (std::size_t j = 0; j < Q.cols(); ++j) {
            const T a = num_abs(T(Q(i, j)));
            if (a > m) m = a;
        }
    return m;
}

/**
 * @param Q generator
 * @param q uniformization rate, which must satisfy q >= max_i |q_ii|
 */
template <class T>
RandomizationResult<T> ctmc_randomization(const Matrix<T>& Q, const T& q) {
    const std::size_t n = Q.rows();
    if (Q.cols() != n) throw InputError("ctmc_randomization: generator is not square");
    if (!(q > num_traits<T>::from_int(0)))
        throw InputError("ctmc_randomization: the uniformization rate must be positive");
    Matrix<T> P(n, n);
    const T one = num_traits<T>::from_int(1);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            P(i, j) = Q(i, j) / q + (i == j ? one : num_traits<T>::from_int(0));
    RandomizationResult<T> r;
    r.P = dtmc_makestochastic(P);
    r.q = q;
    return r;
}

/** Deterministic default rate (21/20) max|Q|; see the defect note above. */
template <class T>
RandomizationResult<T> ctmc_randomization(const Matrix<T>& Q) {
    const T m = ctmc_maxabs(Q);
    // zero-transition uniformization rationale: see _kb/03-api-layer.md (cpp port notes: mc)
    const T q = (m == num_traits<T>::from_int(0)) ? num_traits<T>::from_int(1)
                                                  : T(m * num_traits<T>::from_rational(21, 20));
    return ctmc_randomization(Q, q);
}

}  // namespace mc
}  // namespace line

#endif  // LINE_API_MC_CTMC_RANDOMIZATION_H
