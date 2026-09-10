/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_MAP2_FIT_IDC_H
#define LINE_API_MAM_MAP2_FIT_IDC_H

/**
 * Fit a MAP(2) to three moments and an asymptotic index of dispersion.
 *
 * Templated port of matlab/src/api/mam/map2_fit_idc.m, mirrored by the JAR and
 * native Python.
 *
 * A MAP(2) has a geometrically decaying autocorrelation, so its index of
 * dispersion obeys
 *
 *   I = SCV + (SCV - 1) g2 / (1 - g2),
 *
 * as reported in Section 5.2.2 of Casale, Mi, Cherkasova and Smirni, IEEE Trans.
 * Soft. Eng. 37(5), 2011. Inverting it in closed form gives
 * g2 = (I - SCV)/(I - 1), and that decay rate is handed to `map2_fit`, the
 * explicit inverse characterization of Heindl, Horvath and Gross. A third moment
 * outside the feasible region is replaced by its lower limit (3/2) e2^2 / e1,
 * the largest heavy-tail decay a MAP(2) admits.
 *
 * THE EXPONENTIAL FALLBACK IS NOT MERELY A FEASIBILITY GUARD, and must not be
 * relaxed. When SCV <= 1 or I < SCV the reference returns an exponential,
 * because a flow-equivalent server whose service is exponential and load
 * dependent is EXACT for a product-form subnetwork by Norton's theorem, whereas
 * any MAP(2) fitted to the marginal inter-departure statistics is not: the
 * departure stream of a subnetwork is not independent of the rest of the model.
 * Fitting the sub-exponential SCV of a non-bursty aggregate was measured to cost
 * up to 2.2% of throughput on a three-station exponential network that the
 * exponential fallback reproduces exactly.
 *
 * ARITHMETIC: transcendental, inherited from map2_fit.
 */

#include "line/api/mam/map2_fit.h"
#include "line/api/mam/map_transform.h"
#include "line/lang/lang_types.h"
#include "line/num/number.h"
#include "line/util/matrix.h"

namespace line {
namespace mam {

/** The fitted process and which of the reference's five outcomes produced it. */
template <class T>
struct Map2FitIdcResult {
    Map<T> map;
    /**
     * 0 all four descriptors matched, 1 exponential because the burstiness is
     * not representable, 2 third moment clamped to its lower limit, 3 third
     * moment selected automatically, 4 the fit failed and an exponential is
     * returned.
     */
    int status = 0;
};

/**
 * @param e1 mean inter-arrival time
 * @param e2 second moment of the inter-arrival times
 * @param e3 third moment of the inter-arrival times
 * @param I  asymptotic index of dispersion
 */
template <class T>
Map2FitIdcResult<T> map2_fit_idc(const T& e1, const T& e2, const T& e3, const T& I) {
    static_assert(num_traits<T>::has_transcendental,
                  "map2_fit_idc inherits map2_fit's arithmetic");
    const T one = num_traits<T>::from_int(1);
    const T tol = num_traits<T>::from_double(lang::GlobalConstants::FineTol);
    const T scv = T((e2 - e1 * e1) / (e1 * e1));

    Map2FitIdcResult<T> out;
    if (!(scv > one + tol) || I < scv) {
        out.map = map_exponential_mean(e1);
        out.status = 1;
        return out;
    }

    const T g2 = T((I - scv) / (I - one));
    Map2FitResult<T> r = map2_fit(e1, e2, e3, g2);
    if (r.err == 0 && r.has_map) {
        out.map = r.map;
        out.status = 0;
        return out;
    }

    const T e3min = T(num_traits<T>::from_rational(3, 2) + num_traits<T>::from_double(1e-6)) *
                    T(e2 * e2 / e1);
    if (e3 < e3min) {
        r = map2_fit(e1, e2, e3min, g2);
        if (r.err == 0 && r.has_map) {
            out.map = r.map;
            out.status = 2;
            return out;
        }
    }

    // -1 is map2_fit's sentinel for "choose the third moment yourself".
    r = map2_fit(e1, e2, num_traits<T>::from_int(-1), g2);
    if (r.err == 0 && r.has_map) {
        out.map = r.map;
        out.status = 3;
        return out;
    }

    out.map = map_exponential_mean(e1);
    out.status = 4;
    return out;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_MAP2_FIT_IDC_H
