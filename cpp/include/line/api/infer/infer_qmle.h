/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_INFER_INFER_QMLE_H
#define LINE_API_INFER_INFER_QMLE_H

/**
 * Queue-length-based maximum-likelihood estimator of the service demands of a
 * closed queueing network.
 *
 * Templated port of matlab/src/api/infer/infer_qmle.m. The JAR has no
 * counterpart: jline/api/infer/ carries the LQN identification classes only,
 * so MATLAB is the sole reference.
 *
 * From the observed mean queue lengths Q, the populations N and the think
 * times Z, the demand of class j at station i is estimated by
 *
 *   D(i,j) = Q(i,j) / (N_j - sum_k Q(k,j)) * Z_j / (1 + sum_s Q(i,s) - Q(i,j)/N_j)
 *
 * i.e. the arrival-theorem residence time inverted for the demand, with the
 * think-time population N_j - sum_k Q(k,j) supplying the class throughput.
 *
 * ARITHMETIC: additions, multiplications and divisions of the inputs only, so
 * a finite field computation, exact in the exact instantiation.
 */

#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace infer {

/**
 * @param Q (M x R) mean queue lengths
 * @param N (R) population per class
 * @param Z (R) think time per class
 * @return  (M x R) estimated service demands
 */
template <class T>
Matrix<T> infer_qmle(const Matrix<T>& Q, const std::vector<T>& N, const std::vector<T>& Z) {
    const std::size_t M = Q.rows(), R = Q.cols();
    if (N.size() != R) throw InputError("infer_qmle: Q and N disagree on the class count");
    if (Z.size() != R) throw InputError("infer_qmle: Q and Z disagree on the class count");
    const T one = num_traits<T>::from_int(1);
    const T zero = num_traits<T>::from_int(0);

    // column sums of Q, the in-network population of each class
    std::vector<T> colsum(R, zero);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t j = 0; j < R; ++j) colsum[j] += Q(i, j);

    Matrix<T> D(M, R, zero);
    for (std::size_t i = 0; i < M; ++i) {
        T rowsum = zero;
        for (std::size_t s = 0; s < R; ++s) rowsum += Q(i, s);
        for (std::size_t j = 0; j < R; ++j) {
            const T think_pop = N[j] - colsum[j];
            if (think_pop == zero)
                throw NumericError("infer_qmle: a class has its whole population in the queues");
            if (N[j] == zero) throw InputError("infer_qmle: zero population");
            const T den = one + rowsum - Q(i, j) / N[j];
            if (den == zero) throw NumericError("infer_qmle: singular arrival-theorem denominator");
            D(i, j) = Q(i, j) / think_pop * Z[j] / den;
        }
    }
    return D;
}

}  // namespace infer
}  // namespace line

#endif  // LINE_API_INFER_INFER_QMLE_H
