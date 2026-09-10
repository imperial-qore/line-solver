/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_QSYS_QSYS_MAPG1_H
#define LINE_API_QSYS_QSYS_MAPG1_H

/**
 * The MAP/G/1 FCFS queue, by moment-matching the general service time to a
 * phase-type distribution.
 *
 * Templated port of matlab/src/api/qsys/qsys_mapg1.m. The reference has two
 * parts: a service fit, written in the .m file itself, and the queue solution,
 * delegated to BUTools' MMAPPH1FCFS. This port transcribes the first and
 * replaces the second with the port's own QBD route, exactly as
 * qsys_mapph1.h already does -- MAP/PH/1 is a MAP/MAP/1 queue with a renewal
 * service process, and qsys_mapph1 reproduces MMAPPH1FCFS to 1e-15 on that
 * family (see the measured agreement table in qsys_mapph1.h).
 *
 * THE SERVICE FIT, branch by branch, is the reference's own fitServiceToPH:
 *
 *   3 or more moments : an acyclic PH matching (m1, m2, m3)
 *   exactly 2 moments : cv2 = m2/m1^2 - 1, and then
 *                       cv2 <= 0  -> Erlang-k, k = max(1, round(1/max(cv2, 0.01)))
 *                       cv2 <  1  -> Erlang-k, k = max(2, round(1/cv2))
 *                       cv2 == 1  -> exponential
 *                       cv2 >  1  -> balanced-means 2-phase hyperexponential,
 *                                    p = (1 + sqrt((cv2-1)/(cv2+1)))/2,
 *                                    rates 2p/m1 and 2(1-p)/m1
 *   1 moment          : exponential
 *
 * The Erlang branches match m1 exactly and m2 only through the rounded k, and
 * the hyperexponential branch matches both m1 and cv2 exactly (its cv2 is
 * 1/(2p(1-p)) - 1, which is the requested one for that p). Both are the
 * reference's approximations, reproduced rather than improved.
 *
 * THE 3-MOMENT BRANCH IS THE ONE SUBSTITUTION. The reference calls BUTools'
 * APHFrom3Moments, which this port does not transcribe (the same position
 * qsys_mapph1.h takes on MMAPPH1FCFS). It uses line::mam::aph_fit instead, the
 * m3a implementation of
 * the same Bobbio-Horvath-Telek canonical APH. MEASURED: on
 * (m1, m2, m3) = (1/3, 1/6, 1/9), the moments of Erlang(2) with mean 1/3, the
 * two agree on the order (3), on the sub-generator
 * ([-12 12 0; 0 -12 12; 0 0 -5]) and on the entry vector ([0.8 0 0.2]) to
 * 1e-14, so the queue results agree as well. They are not guaranteed to pick
 * the same representation everywhere -- the order search and the branch
 * conditions are written differently -- and where they differ the queue
 * results will differ too, because a queue depends on the whole service
 * distribution and not on its first three moments. A caller who needs the
 * reference's exact PH can pass it to qsys_mapph1 directly.
 *
 * WHAT IS NOT RETURNED. MMAPPH1FCFS also returns higher queue-length and
 * sojourn-time moments (ncMoms, stMoms beyond the first). The QBD route gives
 * the means and the queue-length distribution, and no higher moments are
 * fabricated from the truncated distribution: queueLengthMoments and
 * sojournTimeMoments simply have no counterpart here.
 *
 * REFERENCE DEFECTS in qsys_mapg1.m:
 *
 *  1. THE cv2 <= 0 BRANCH SILENTLY BUILDS AN ERLANG-100. Line 128 computes
 *     k = max(1, round(1/max(cv2, 0.01))), and for any cv2 <= 0 the max pins
 *     the denominator at 0.01, so k = 100 whatever the moments were. Nothing
 *     warns. Deterministic service therefore enters the QBD with a 100-phase
 *     service process (cv2 = 0.01, not 0), which is both expensive and a
 *     silent modelling decision. MATLAB reproduction, from matlab/:
 *         r = qsys_mapg1(-2, 2, [1/3, (1/3)^2])   % cv2 = 0 exactly
 *         r.meanQueueLength -> 1.33999999999909
 *     against the exact M/D/1 value 4/3 + ... (the Erlang-100 answer is
 *     1.3399999999 rather than the M/D/1 1.3333...). Reproduced here, since
 *     it is the reference's model choice, and pinned by a test that asserts
 *     the phase count is 100.
 *  2. m2 and m3 are read without any feasibility check: a moment set that is
 *     not PH-representable reaches the fit and fails there.
 *  3. The reference computes rho from the INPUT mean 1/serviceMoments(1) and
 *     not from the fitted PH, so a branch whose fit does not preserve m1 would
 *     report a utilization inconsistent with its own service process. Every
 *     branch does preserve m1, so this is latent rather than active; the port
 *     reproduces the reference's formula.
 *
 * ARITHMETIC. Gated on num_traits<T>::has_transcendental, inherited from
 * qsys_mapmap1 (the cyclic reduction behind R) and from aph_fit.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/mam/aph_fit.h"
#include "line/api/mam/map_moment.h"
#include "line/api/qsys/qsys_mapmap1.h"
#include "line/api/qsys/qsys_mapph1.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/matrix.h"

namespace line {
namespace qsys {

/** Which branch of fitServiceToPH was taken, for the caller and for tests. */
enum class MapG1ServiceFit {
    Exponential,       ///< one moment, or cv2 exactly 1
    Erlang,            ///< 0 < cv2 < 1, or the cv2 <= 0 branch
    Hyperexponential,  ///< cv2 > 1
    Acyclic            ///< three or more moments
};

/** Result of qsys_mapg1. */
template <class T>
struct MapG1Result {
    T meanQueueLength;               ///< E[N], number in system
    T meanWaitingTime;               ///< max(0, E[W] - m1), as in the reference
    T meanSojournTime;               ///< E[W]
    T utilization;                   ///< lambda m1
    std::vector<T> queueLengthDist;  ///< P(N = n), n = 0, 1, ...
    mam::Map<T> serviceFit;          ///< the fitted PH as its renewal MAP
    MapG1ServiceFit fitKind;         ///< which branch produced it
    std::size_t servicePhases;       ///< order of the fitted PH
};

namespace detail {

/** Exponential PH of the given mean, as its renewal MAP. */
template <class T>
mam::Map<T> mapg1_exponential(const T& mean) {
    if (!(mean > num_traits<T>::from_int(0)))
        throw InputError("qsys_mapg1: the mean service time must be positive");
    mam::Map<T> m;
    m.D0 = Matrix<T>(1, 1, T(-num_traits<T>::from_int(1) / mean));
    m.D1 = Matrix<T>(1, 1, T(num_traits<T>::from_int(1) / mean));
    return m;
}

/** Erlang-k PH of the given mean, as its renewal MAP (entry in phase 1). */
template <class T>
mam::Map<T> mapg1_erlang(const T& mean, std::size_t k) {
    if (k == 0) throw InputError("qsys_mapg1: an Erlang fit needs at least one phase");
    const T zero = num_traits<T>::from_int(0);
    const T mu = T(num_traits<T>::from_int(static_cast<long>(k)) / mean);
    mam::Map<T> m;
    m.D0 = Matrix<T>(k, k, zero);
    m.D1 = Matrix<T>(k, k, zero);
    for (std::size_t i = 0; i < k; ++i) {
        m.D0(i, i) = -mu;
        if (i + 1 < k) m.D0(i, i + 1) = mu;
    }
    m.D1(k - 1, 0) = mu;  // completion restarts in phase 1
    return m;
}

/**
 * Balanced-means two-phase hyperexponential with the given mean and squared
 * coefficient of variation, as its renewal MAP.
 */
template <class T>
mam::Map<T> mapg1_hyperexp2(const T& mean, const T& cv2) {
    using std::sqrt;
    const T zero = num_traits<T>::from_int(0), one = num_traits<T>::from_int(1);
    const T two = num_traits<T>::from_int(2);
    const T p = T(T(one + sqrt(T(T(cv2 - one) / T(cv2 + one)))) / two);
    const T l1 = T(two * p / mean), l2 = T(two * T(one - p) / mean);
    mam::Map<T> m;
    m.D0 = Matrix<T>(2, 2, zero);
    m.D1 = Matrix<T>(2, 2, zero);
    m.D0(0, 0) = -l1;
    m.D0(1, 1) = -l2;
    m.D1(0, 0) = l1 * p;
    m.D1(0, 1) = l1 * T(one - p);
    m.D1(1, 0) = l2 * p;
    m.D1(1, 1) = l2 * T(one - p);
    return m;
}

}  // namespace detail

/**
 * Fit a general service time to a PH, following qsys_mapg1.m's
 * fitServiceToPH. Exposed separately because the branch selection is the part
 * of the reference that is transcribed verbatim, and a caller may want the
 * fitted process without solving a queue.
 *
 * @param moments the first one, two or three raw moments of the service time
 * @param kind    out: which branch was taken
 */
template <class T>
mam::Map<T> qsys_mapg1_service_fit(const std::vector<T>& moments, MapG1ServiceFit& kind) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_mapg1_service_fit requires transcendental arithmetic");
    if (moments.empty()) throw InputError("qsys_mapg1: no service moments given");
    const T one = num_traits<T>::from_int(1);
    const T m1 = moments[0];

    if (moments.size() >= 3) {
        kind = MapG1ServiceFit::Acyclic;
        return mam::aph_fit(moments[0], moments[1], moments[2]).aph;
    }
    if (moments.size() == 2) {
        const T cv2 = T(T(moments[1] / T(m1 * m1)) - one);
        const T zero = num_traits<T>::from_int(0);
        if (cv2 <= zero) {
            // REFERENCE DEFECT 1: the floor at 0.01 makes this Erlang-100.
            const T floored = num_traits<T>::from_double(0.01);
            const T den = cv2 > floored ? cv2 : floored;  // max(cv2, 0.01)
            const double kd = std::round(1.0 / num_traits<T>::to_double(den));
            const std::size_t k = kd < 1.0 ? 1u : static_cast<std::size_t>(kd);
            kind = MapG1ServiceFit::Erlang;
            return detail::mapg1_erlang(m1, k);
        }
        if (cv2 < one) {
            const double kd = std::round(1.0 / num_traits<T>::to_double(cv2));
            const std::size_t k = kd < 2.0 ? 2u : static_cast<std::size_t>(kd);
            kind = MapG1ServiceFit::Erlang;
            return detail::mapg1_erlang(m1, k);
        }
        if (cv2 == one) {
            kind = MapG1ServiceFit::Exponential;
            return detail::mapg1_exponential(m1);
        }
        kind = MapG1ServiceFit::Hyperexponential;
        return detail::mapg1_hyperexp2(m1, cv2);
    }
    kind = MapG1ServiceFit::Exponential;
    return detail::mapg1_exponential(m1);
}

/**
 * The MAP/G/1 FCFS queue.
 *
 * @param arrival   arrival MAP (D0, D1)
 * @param moments   the first one, two or three raw moments of the service time
 * @param dist_size how many entries of queueLengthDist to materialize
 *                  (the reference's numQLProbs, default 100)
 */
template <class T>
MapG1Result<T> qsys_mapg1(const mam::Map<T>& arrival, const std::vector<T>& moments,
                          std::size_t dist_size) {
    static_assert(num_traits<T>::has_transcendental,
                  "qsys_mapg1 requires transcendental arithmetic");
    MapG1Result<T> r;
    r.serviceFit = qsys_mapg1_service_fit(moments, r.fitKind);
    r.servicePhases = r.serviceFit.order();

    const MapMap1Result<T> q = qsys_mapmap1(arrival, r.serviceFit, dist_size);
    r.meanQueueLength = q.meanQueueLength;
    r.meanSojournTime = q.meanSojournTime;
    r.queueLengthDist = q.queueLengthDist;

    // The reference recomputes both of these from the INPUT mean rather than
    // from the fitted process; see reference defect 3.
    const T zero = num_traits<T>::from_int(0);
    const T lambda = mam::map_lambda(arrival);
    r.utilization = T(lambda * moments[0]);
    const T wait = T(q.meanSojournTime - moments[0]);
    r.meanWaitingTime = wait < zero ? zero : wait;
    return r;
}

/** qsys_mapg1 with the reference's default of 100 materialized levels. */
template <class T>
MapG1Result<T> qsys_mapg1(const mam::Map<T>& arrival, const std::vector<T>& moments) {
    return qsys_mapg1(arrival, moments, static_cast<std::size_t>(100));
}

}  // namespace qsys
}  // namespace line

#endif  // LINE_API_QSYS_QSYS_MAPG1_H
