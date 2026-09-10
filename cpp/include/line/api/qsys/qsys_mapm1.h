/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_QSYS_MAPM1_H
#define LINE_API_QSYS_QSYS_MAPM1_H

/**
 * The MAP/M/1 FCFS queue, the single-server case of MAP/M/c.
 *
 * ALGORITHM, AND HOW IT DIFFERS FROM THE MATLAB REFERENCE.
 * matlab/src/api/qsys/qsys_mapm1.m is a one-line wrapper for qsys_mapmc with
 * c = 1, which in turn calls Q-MAM's Q_CT_MAP_M_C, which is not transcribed
 * here. This port keeps
 * the same delegation and calls this port's own qsys_mapmc (see qsys_mapmc.h),
 * which builds the quasi-birth-death process directly: at c = 1 the level
 * blocks are A0 = D1, A1 = D0 - mu I, A2 = mu I with a single boundary level,
 * and R comes from cyclic reduction.
 * Different algorithm, same quantity.
 *
 * At c = 1 this is also exactly a MAP/MAP/1 queue with an exponential service
 * MAP, so qsys_mapmap1(arrival, {[-mu], [mu]}) must return the same numbers by
 * a completely different route (cyclic reduction on the Kronecker-product QBD
 * rather than cyclic reduction on the arrival-phase QBD). That
 * cross-check is one of the assertions in the test file, and holds to 1e-13.
 *
 * MEASURED AGREEMENT (MATLAB R2025a, T = double). Metric order is
 * meanQueueLength / meanWaitingTime / meanSojournTime / utilization.
 *  - M/M/1 collapse, qsys_mapm1(D0 = [-2], D1 = [2], mu = 3): MATLAB
 *    1.999999994584355 / 0.6666666665442794 / 0.9999999998776128 /
 *    0.6666666666666666. This port returns the textbook 2 / 0.666666666666667
 *    / 1 / 0.666666666666667 to 1e-15, so the relative differences against
 *    MATLAB are 2.7e-9, 1.7e-10 and 1.2e-10. The residual is the reference's
 *    maxNumComp = 500 level truncation; the port sums the geometric tail in
 *    closed form and matches the exact formulas.
 *  - Correlated MMPP2 arrivals D0 = [-2.5 0.2; 0.1 -0.7], D1 = diag(2.3, 0.6)
 *    (lambda = 7/6), mu = 2: MATLAB 2.795385119831413 / 1.896044396731732 /
 *    2.396044396731732 / 0.5833333333333333; port 2.79538513100956 /
 *    1.89604439800819 / 2.39604439800819 / 0.583333333333333. Relative
 *    differences 4.0e-9, 6.7e-10, 5.3e-10.
 *  - Erlang-2 arrivals D0 = [-4 4; 0 -4], D1 = [0 0; 4 0] (lambda = 2),
 *    mu = 3: MATLAB 1.568729300352212 / 0.4510313188166638 /
 *    0.7843646521499972 / 0.6666666666666666; port 1.56872930440884 /
 *    0.451031318871088 / 0.784364652204421 / 0.666666666666667. Relative
 *    differences 2.6e-9, 1.2e-10, 6.9e-11.
 *
 * The cross-check against qsys_mapmap1 mentioned above was run on both
 * non-Poisson instances and agrees to 5e-16 on meanQueueLength and
 * meanWaitingTime, which identifies the 1e-9 residuals above as the
 * reference's truncation rather than this port's error.
 *
 * ARITHMETIC. Gated on num_traits<T>::has_transcendental, inherited from
 * qsys_mapmc: cyclic reduction drives R to a tolerance and never terminates in
 * a finite number of field operations. See qbd_r.h.
 */

#include <cstddef>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/api/qsys/qsys_mapmc.h"
#include "line/num/number.h"

namespace line {
namespace qsys {

/**
 * MAP/M/1 by the matrix-geometric solution, i.e. qsys_mapmc at c = 1.
 *
 * @param arrival    arrival MAP (D0, D1)
 * @param mu         exponential service rate
 * @param dist_size  how many entries of queueLengthDist to materialize
 */
template <class T>
MapMcResult<T> qsys_mapm1(const mam::Map<T>& arrival, const T& mu, std::size_t dist_size) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_mapm1 requires transcendental arithmetic");
    return qsys_mapmc(arrival, mu, 1u, dist_size);
}

/** qsys_mapm1 with 100 materialized levels. */
template <class T>
MapMcResult<T> qsys_mapm1(const mam::Map<T>& arrival, const T& mu) {
    return qsys_mapm1(arrival, mu, static_cast<std::size_t>(100));
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_QSYS_MAPM1_H
