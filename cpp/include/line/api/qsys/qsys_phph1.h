/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_QSYS_PHPH1_H
#define LINE_API_QSYS_QSYS_PHPH1_H

/**
 * The PH/PH/1 FCFS queue.
 *
 * ALGORITHM, AND HOW IT DIFFERS FROM THE MATLAB REFERENCE.
 * matlab/src/api/qsys/qsys_phph1.m converts the arrival PH to a MAP and then
 * calls BUTools' MMAPPH1FCFS, which is not transcribed here. This port computes
 * the SAME quantities from the port's
 * own QBD machinery. Both the arrival PH (alpha, T) and the service PH
 * (beta, S) become renewal MAPs,
 *
 *     arrival: D0 = T, D1 = (-T e) alpha,   service: D0 = S, D1 = (-S e) beta,
 *
 * and the resulting MAP/MAP/1 queue is solved as the level-independent QBD
 * described in qsys_mapmap1.h. Different algorithm, same quantity. Both
 * processes are renewal by construction, so the reference and the port model
 * the identical stochastic system; only the numerical route differs.
 *
 * MEASURED AGREEMENT (MATLAB R2025a, T = double). Metric order is
 * meanQueueLength / meanWaitingTime / meanSojournTime / utilization.
 *  - M/M/1 collapse, qsys_phph1(alpha = [1], T = [-2], beta = [1], S = [-3]):
 *    MATLAB 1.999999999999999 / 0.6666666666666656 / 0.9999999999999989 /
 *    0.6666666666666666; the port returns the textbook 2 / 0.666666666666668 /
 *    1 / 0.666666666666667. Relative differences 5e-16, 1.8e-15, 1.1e-15.
 *  - Erlang-2 arrivals alpha = [1 0], T = [-4 4; 0 -4] (lambda = 2) with
 *    Erlang-2 service beta = [1 0], S = [-6 6; 0 -6] (mean 1/3): MATLAB
 *    1.250000000000001 / 0.2916666666666675 / 0.6250000000000008 /
 *    0.6666666666666666; port 1.25 / 0.291666666666667 / 0.625000000000001 /
 *    0.666666666666667. Relative differences 8e-16, 1.7e-15, 1.1e-15.
 *  - Hyperexponential arrivals alpha = [0.7 0.3], T = diag(-5, -1)
 *    (mean 0.44, lambda = 25/11) with exponential service beta = [1],
 *    S = [-10/3]: MATLAB 3.273624198148769 / 1.140394647185455 /
 *    1.440394647185455 / 0.6818181818181819; port 3.27362419814877 /
 *    1.14039464718546 / 1.44039464718546 / 0.681818181818182. Relative
 *    differences below 4e-15 on every metric.
 *
 * The port reproduces the BUTools reference to machine precision throughout
 * this family; the residual is the reference's own rounding, not a truncation,
 * because these level distributions decay fast enough that the reference's
 * numQLProbs = 100 cutoff is not visible. The slower-decaying instances in
 * qsys_mapph1.h expose the reference's truncation floor instead.
 *
 * ARITHMETIC. Gated on num_traits<T>::has_transcendental, inherited from
 * qsys_mapmap1 and ultimately from the cyclic reduction that produces R; see
 * qbd_r.h. The two PH-to-MAP conversions are exact matrix algebra.
 */

#include <cstddef>
#include <vector>

#include "line/api/qsys/qsys_mapmap1.h"
#include "line/api/qsys/qsys_mapph1.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace qsys {

/**
 * PH/PH/1 by the exact QBD solution of the equivalent MAP/MAP/1 queue.
 *
 * @param alpha      PH arrival entry vector, length n
 * @param Tm         PH arrival sub-generator, n x n
 * @param beta       PH service entry vector, length m
 * @param S          PH service sub-generator, m x m
 * @param dist_size  how many entries of queueLengthDist to materialize
 */
template <class T>
MapMap1Result<T> qsys_phph1(const std::vector<T>& alpha, const Matrix<T>& Tm,
                            const std::vector<T>& beta, const Matrix<T>& S,
                            std::size_t dist_size) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_phph1 requires transcendental arithmetic");
    return qsys_mapmap1(detail::ph_to_map(alpha, Tm), detail::ph_to_map(beta, S), dist_size);
}

/** qsys_phph1 with 100 materialized levels, the reference's numQLProbs. */
template <class T>
MapMap1Result<T> qsys_phph1(const std::vector<T>& alpha, const Matrix<T>& Tm,
                            const std::vector<T>& beta, const Matrix<T>& S) {
    return qsys_phph1(alpha, Tm, beta, S, static_cast<std::size_t>(100));
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_QSYS_PHPH1_H
