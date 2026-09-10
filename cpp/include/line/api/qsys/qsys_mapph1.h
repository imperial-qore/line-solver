/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_QSYS_MAPPH1_H
#define LINE_API_QSYS_QSYS_MAPPH1_H

/**
 * The MAP/PH/1 FCFS queue.
 *
 * ALGORITHM, AND HOW IT DIFFERS FROM THE MATLAB REFERENCE.
 * matlab/src/api/qsys/qsys_mapph1.m computes these quantities by calling
 * BUTools' MMAPPH1FCFS, which is not transcribed here. This port computes the
 * SAME quantities from the port's own
 * QBD machinery: a PH service (sigma, S) is the renewal MAP
 *
 *     D0 = S,    D1 = (-S e) sigma,
 *
 * so MAP/PH/1 is literally a MAP/MAP/1 queue with a renewal service process,
 * and the problem is handed to qsys_mapmap1 (see qsys_mapmap1.h for the QBD
 * construction). Different algorithm, same quantity. Because the service MAP
 * built this way IS a renewal process, none of the caveats about discarded
 * service correlation in qsys_mapmap1.h apply to this function: the reference
 * and the port model the identical stochastic system here.
 *
 * MEASURED AGREEMENT (MATLAB R2025a, T = double). Metric order is
 * meanQueueLength / meanWaitingTime / meanSojournTime / utilization.
 *  - M/M/1 collapse, qsys_mapph1(D0 = [-2], D1 = [2], sigma = [1], S = [-3]):
 *    MATLAB 1.999999999999999 / 0.6666666666666656 / 0.9999999999999989 /
 *    0.6666666666666666; the port returns the textbook 2 / 0.666666666666668 /
 *    1 / 0.666666666666667. Relative differences 5e-16, 1.8e-15, 1.1e-15,
 *    1e-16.
 *  - Correlated MMPP2 arrival D0 = [-2.5 0.2; 0.1 -0.7], D1 = diag(2.3, 0.6)
 *    (lambda = 7/6) with Erlang-2 service sigma = [1 0], S = [-6 6; 0 -6]:
 *    MATLAB 0.8364637529649164 / 0.3836355977794533 / 0.7169689311127866 /
 *    0.3888888888888888; port 0.836463752964917 / 0.383635597779453 /
 *    0.716968931112786 / 0.388888888888889. Relative differences below 1e-15
 *    on every metric.
 *  - Erlang-2 arrival D0 = [-4 4; 0 -4], D1 = [0 0; 4 0] (lambda = 2) with
 *    hyperexponential service sigma = [0.6 0.4], S = diag(-8, -1.6): MATLAB
 *    2.133647429543165 / 0.7418237147715814 / 1.066823714771581 / 0.65; port
 *    2.13364742954318 / 0.741823714771588 / 1.06682371477159 / 0.65. Relative
 *    differences below 1e-15 on every metric.
 *
 * The port therefore reproduces the BUTools reference to machine precision on
 * this family. It deliberately does NOT reproduce LINE's own MATLAB
 * qbd_mapmap1 on the equivalent MAP pair, which returns 0.8364637492083693 and
 * 2.133647423495968 for the second and third cases: those are 4.5e-9 and
 * 2.8e-9 BELOW the values above because MATLAB's qbd_mapmap1 sums the level
 * distribution until the accumulated mass reaches 1 - 1e-10 and discards the
 * remaining tail times its level index, while this port sums the whole
 * geometric tail in closed form (see qbd_mapmap1.h).
 *
 * ARITHMETIC. Gated on num_traits<T>::has_transcendental, inherited from
 * qsys_mapmap1 and ultimately from the cyclic reduction that produces R; see
 * qbd_r.h. The PH-to-MAP conversion itself is exact matrix algebra.
 */

#include <cstddef>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/api/qsys/qsys_mapmap1.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace qsys {

namespace detail {

/**
 * The renewal MAP of a PH distribution: D0 = S, D1 = (-S e) beta. Exact matrix
 * algebra, so un-gated.
 */
template <class T>
mam::Map<T> ph_to_map(const std::vector<T>& beta, const Matrix<T>& S) {
    const std::size_t m = S.rows();
    if (S.cols() != m || beta.size() != m)
        throw InputError("ph_to_map: beta and S dimensions are inconsistent");
    mam::Map<T> out;
    out.D0 = S;
    out.D1 = Matrix<T>(m, m, num_traits<T>::from_int(0));
    for (std::size_t i = 0; i < m; ++i) {
        T exit = num_traits<T>::from_int(0);
        for (std::size_t j = 0; j < m; ++j) exit -= S(i, j);
        for (std::size_t j = 0; j < m; ++j) out.D1(i, j) = exit * beta[j];
    }
    return out;
}

}  // namespace detail

/**
 * MAP/PH/1 by the exact QBD solution of the equivalent MAP/MAP/1 queue.
 *
 * @param arrival    arrival MAP (D0, D1)
 * @param sigma      PH service entry vector, length m
 * @param S          PH service sub-generator, m x m
 * @param dist_size  how many entries of queueLengthDist to materialize
 */
template <class T>
MapMap1Result<T> qsys_mapph1(const mam::Map<T>& arrival, const std::vector<T>& sigma,
                             const Matrix<T>& S, std::size_t dist_size) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_mapph1 requires transcendental arithmetic");
    return qsys_mapmap1(arrival, detail::ph_to_map(sigma, S), dist_size);
}

/** qsys_mapph1 with 100 materialized levels, the reference's numQLProbs. */
template <class T>
MapMap1Result<T> qsys_mapph1(const mam::Map<T>& arrival, const std::vector<T>& sigma,
                             const Matrix<T>& S) {
    return qsys_mapph1(arrival, sigma, S, static_cast<std::size_t>(100));
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_QSYS_MAPPH1_H
