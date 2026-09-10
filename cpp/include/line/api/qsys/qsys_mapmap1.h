/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_QSYS_MAPMAP1_H
#define LINE_API_QSYS_QSYS_MAPMAP1_H

/**
 * The MAP/MAP/1 FCFS queue: mean number in system, waiting time, sojourn time,
 * utilization and the queue-length distribution.
 *
 * ALGORITHM, AND HOW IT DIFFERS FROM THE MATLAB REFERENCE.
 * matlab/src/api/qsys/qsys_mapmap1.m obtains these quantities by calling
 * BUTools' MMAPPH1FCFS, which is not transcribed here; this port computes the
 * SAME quantities from the port's own
 * quasi-birth-death machinery (line/api/mam/qbd_mapmap1.h), which solves the
 * level-independent QBD whose level is the number in system and whose phase is
 * the pair (arrival phase, service phase). The algorithm is therefore
 * different, the quantity is the same.
 *
 * The reference does one further thing this port deliberately does not. Before
 * calling BUTools it replaces the service MAP (D0,D1) by the phase-type
 * distribution
 *
 *     sigma = theta D1 / (theta D1 e),   S = D0,   theta (D0+D1) = 0,
 *
 * i.e. by the stationary marginal of a service time. That embedded PH keeps the
 * service-time MARGINAL exactly but DISCARDS the serial correlation between
 * consecutive service times, so the reference solves a MAP/PH/1 queue in place
 * of the MAP/MAP/1 queue it names. The two coincide exactly when the service
 * MAP is a renewal process (D1 = (-D0 e) sigma, which includes every PH-renewal
 * service), and diverge otherwise. This port solves the genuine MAP/MAP/1 QBD,
 * so on a correlated service MAP it does not agree with the reference and is
 * not meant to: see the measured numbers below.
 *
 * MEASURED AGREEMENT (MATLAB R2025a, T = double). Metric order is
 * meanQueueLength / meanWaitingTime / meanSojournTime / utilization. Arrival
 * MAPs used below:
 *   P(2) = Poisson, D0 = [-2], D1 = [2]
 *   C    = correlated MMPP2, D0 = [-2.5 0.2; 0.1 -0.7], D1 = diag(2.3, 0.6),
 *          lambda = 7/6
 *
 * 1. M/M/1 collapse, qsys_mapmap1(P(2), service D0 = [-3], D1 = [3]). MATLAB
 *    1.999999999999999 / 0.6666666666666656 / 0.9999999999999989 /
 *    0.6666666666666666; port 2 / 0.666666666666668 / 1 / 0.666666666666667,
 *    the textbook values. Relative differences 5e-16, 1.8e-15, 1.1e-15.
 *
 * 2. Renewal service, where the reference's PH reduction is exact. Arrival C,
 *    service the Erlang-2 renewal MAP D0 = [-6 6; 0 -6], D1 = [0 0; 6 0].
 *    MATLAB 0.8364637529649164 / 0.3836355977794533 / 0.7169689311127866 /
 *    0.3888888888888888; port 0.836463752964917 / 0.383635597779453 /
 *    0.716968931112786 / 0.388888888888889. Relative differences below 1e-15
 *    on every metric: with a renewal service the two model the same system and
 *    the port reproduces the reference to machine precision.
 *
 * 3. Correlated service, where the reference's PH reduction is NOT exact.
 *    Arrival C, service MAP D0 = [-5 0.4; 0.2 -1.4], D1 = diag(4.6, 1.2)
 *    (lambda_s = 7/3, so rho = 1/2). MATLAB qsys_mapmap1 reports
 *    2.0305355320573 / 1.311887598906254 / 1.740459027477683 / 0.5; this port
 *    returns 2.54755753480255 / 1.75504931554505 / 2.18362074411647 / 0.5,
 *    relative differences of 2.03e-1, 2.53e-1 and 2.03e-1. The port's value is
 *    the correct MAP/MAP/1 answer: LINE's own MATLAB qbd_mapmap1 on the
 *    identical pair of MAPs returns QN = 2.547554867631855, UN =
 *    0.5000000000000002 and RN = 2.183618457970162, i.e. the port to within
 *    1.05e-6. That residual 1.05e-6 is MATLAB qbd_mapmap1's truncated level
 *    sum, which stops once the accumulated mass reaches 1 - 1e-10 and so
 *    discards the tail times its level index; this port sums the geometric
 *    tail in closed form (see qbd_mapmap1.h). The 20% gap against
 *    qsys_mapmap1.m is a different thing entirely: it is the service
 *    correlation the reference's PH reduction throws away.
 *
 * meanWaitingTime is meanSojournTime minus the mean service time 1/lambda_s,
 * and utilization is lambda_a/lambda_s, both matching the reference's own
 * definitions.
 *
 * ARITHMETIC. The function is gated on num_traits<T>::has_transcendental
 * because qbd_mapmap1 reaches R by cyclic reduction, a fixed-point iteration
 * driven to a tolerance that does not terminate in a finite number of field
 * operations; running it at exact rational arithmetic would produce a rational
 * with a denominator doubling every step and still not the exact R. See
 * qbd_r.h. Everything downstream of R here -- the mean, the level
 * probabilities, the Little's-law conversions -- is finite exact matrix
 * algebra and carries no additional error.
 */

#include <cstddef>
#include <vector>

#include "line/api/mam/map_moment.h"
#include "line/api/mam/qbd_mapmap1.h"
#include "line/api/mam/qbd_r.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace qsys {

/**
 * Return value of the MAP/MAP/1 family (qsys_mapmap1, qsys_mapph1, qsys_phph1),
 * carrying the same quantities as the MATLAB result struct.
 */
template <class T>
struct MapMap1Result {
    T meanQueueLength;                ///< E[N], number in system
    T meanWaitingTime;                ///< E[Wq], time in queue
    T meanSojournTime;                ///< E[W] = E[Wq] + E[S]
    T utilization;                    ///< rho = lambda_a / lambda_s
    std::vector<T> queueLengthDist;   ///< P(N = n), n = 0, 1, ...
};

/**
 * MAP/MAP/1 by the exact QBD solution.
 *
 * @param arrival    arrival MAP (D0, D1)
 * @param service    service MAP (D0, D1); its correlation IS honoured
 * @param dist_size  how many entries of queueLengthDist to materialize
 */
template <class T>
MapMap1Result<T> qsys_mapmap1(const mam::Map<T>& arrival, const mam::Map<T>& service,
                              std::size_t dist_size) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_mapmap1 requires transcendental arithmetic");
    if (dist_size == 0) throw InputError("qsys_mapmap1: dist_size must be positive");
    const T zero = num_traits<T>::from_int(0);
    const mam::QbdMapMap1Result<T> q = mam::qbd_mapmap1(arrival, service, T(zero), dist_size);

    const T lambda_a = mam::map_lambda(arrival);
    const T lambda_s = mam::map_lambda(service);
    const T mean_service = num_traits<T>::from_int(1) / lambda_s;

    MapMap1Result<T> r;
    r.meanQueueLength = q.QN;
    r.meanSojournTime = q.QN / lambda_a;
    r.meanWaitingTime = r.meanSojournTime - mean_service;
    r.utilization = lambda_a / lambda_s;
    r.queueLengthDist.assign(q.pqueue.rows(), zero);
    for (std::size_t k = 0; k < q.pqueue.rows(); ++k) {
        T s = zero;
        for (std::size_t j = 0; j < q.pqueue.cols(); ++j) s += q.pqueue(k, j);
        r.queueLengthDist[k] = s;
    }
    return r;
}

/** qsys_mapmap1 with 100 materialized levels, the reference's numQLProbs. */
template <class T>
MapMap1Result<T> qsys_mapmap1(const mam::Map<T>& arrival, const mam::Map<T>& service) {
    return qsys_mapmap1(arrival, service, static_cast<std::size_t>(100));
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_QSYS_MAPMAP1_H
