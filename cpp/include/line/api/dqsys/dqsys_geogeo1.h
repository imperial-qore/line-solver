/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_DQSYS_DQSYS_GEOGEO1_H
#define LINE_API_DQSYS_DQSYS_GEOGEO1_H

/**
 * Geo/Geo/1: the discrete-time single-server queue with geometric
 * interarrival and service times.
 *
 * Templated port of matlab/src/api/qsys/dqsys_geogeo1.m. Two timing conventions
 * are in use in the literature and both are supported, as in MATLAB:
 *
 *   LAS_DA (late arrival, delayed access): an arrival in a slot cannot be
 *          served in that slot. Empty probability 1 - rho, mean queue length
 *          a(1-a)/(s-a).
 *   EAS    (early arrival): the arrival is eligible immediately. Empty
 *          probability 1 - r with r = a(1-s)/(s(1-a)), mean queue length
 *          a(1-s)/(s-a).
 *
 * The distinction is not cosmetic: the two conventions give different empty
 * probabilities, different mean queue lengths and different mean service times
 * for the same (a, s). Daduna (LNCS 2046, Cor. 2.7) is the reference that ties
 * the LAS_DA form to the continuous-time Geo/Geo/1 term for term.
 *
 * Everything here is a rational function of a and s, so the exact
 * instantiation gives the stationary quantities with no rounding. That is
 * worth having in the slotted setting, where the interesting regime is
 * s - a small and the double evaluation of a(1-s)/(s(s-a)) loses digits
 * exactly there.
 */

#include <cstddef>
#include <string>

#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace dqsys {

enum class GeoConvention { LAS_DA, EAS };

template <class T>
struct GeoGeo1Result {
    GeoConvention convention = GeoConvention::LAS_DA;
    T arrivalProb;
    T serviceProb;
    T utilization;       ///< rho = a/s
    T throughput;        ///< a
    T emptyProb;
    T ratio;             ///< r = a(1-s)/(s(1-a))
    T meanQueueLength;
    T meanWaitingQueue;
    T meanSojournTime;
    T meanWaitingTime;
    T meanServiceTime;
};

/** Stationary queue-length pmf under the convention of the result. */
template <class T>
T dqsys_geogeo1_pmf(const GeoGeo1Result<T>& r, int n) {
    if (n < 0) throw InputError("dqsys_geogeo1_pmf: queue length must be nonnegative");
    const T one = num_traits<T>::from_int(1);
    if (r.convention == GeoConvention::EAS)
        return (one - r.ratio) * num_pow_int(r.ratio, static_cast<unsigned>(n));
    if (n == 0) return r.emptyProb;
    return r.emptyProb * (r.utilization / (one - r.arrivalProb)) *
           num_pow_int(r.ratio, static_cast<unsigned>(n - 1));
}

/**
 * @param a arrival probability per slot, in (0,1]
 * @param s service completion probability per slot, in (0,1]
 * @param convention slot-boundary convention (late arrival, early arrival)
 */
template <class T>
GeoGeo1Result<T> dqsys_geogeo1(const T& a, const T& s,
                              GeoConvention convention = GeoConvention::LAS_DA) {
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    if (a <= zero || a > one) throw InputError("dqsys_geogeo1: a must lie in (0,1]");
    if (s <= zero || s > one) throw InputError("dqsys_geogeo1: s must lie in (0,1]");
    if (a >= s) throw InputError("dqsys_geogeo1: the load a/s must be strictly less than 1");

    GeoGeo1Result<T> r;
    r.convention = convention;
    r.arrivalProb = a;
    r.serviceProb = s;
    r.utilization = a / s;
    r.throughput = a;
    r.ratio = a * (one - s) / (s * (one - a));
    r.meanWaitingTime = a * (one - s) / (s * (s - a));
    r.meanWaitingQueue = a * r.meanWaitingTime;

    if (convention == GeoConvention::LAS_DA) {
        r.emptyProb = one - r.utilization;
        r.meanQueueLength = a * (one - a) / (s - a);
        r.meanSojournTime = (one - a) / (s - a);
        r.meanServiceTime = one / s;
    } else {
        r.emptyProb = one - r.ratio;
        r.meanQueueLength = a * (one - s) / (s - a);
        r.meanSojournTime = (one - s) / (s - a);
        r.meanServiceTime = (one - s) / s;
    }
    return r;
}

}  // namespace dqsys
}  // namespace line

#endif  // LINE_API_DQSYS_DQSYS_GEOGEO1_H
