/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_AMAP2_FIT_GAMMA_H
#define LINE_API_MAM_AMAP2_FIT_GAMMA_H

/**
 * AMAP(2) fit of three moments and the autocorrelation decay rate
 * (matlab/lib/m3a/m3a/amap2/amap2_fit_gamma.m).
 *
 * A unit SCV short-circuits to a Poisson process, since the AMAP(2) canonical
 * forms are degenerate there. Otherwise every exact solution is enumerated
 * with amap2_fitall_gamma, normalized, and the first is returned.
 *
 * Not ported: the reference's approximate branch, which calls
 * amap2_adjust_gamma to relax (M2, M3, GAMMA) into the feasible region. All
 * four of its methods are driven by patternsearch / PSwarm / fmincon, so there
 * is no closed form to port. When no exact solution exists this port takes the
 * same final fallback as the in-tree MATLAB does when adjustment fails: a
 * Poisson process matching the mean, flagged through Amap2FitGammaResult.
 *
 * Gated on transcendental arithmetic through amap2_fitall_gamma.
 */

#include <vector>

#include "line/api/mam/amap2_fitall_gamma.h"
#include "line/api/mam/map_moment.h"
#include "line/api/mam/map_transform.h"
#include "line/num/number.h"

namespace line {
namespace mam {

/** Result of amap2_fit_gamma. */
template <class T>
struct Amap2FitGammaResult {
    Map<T> amap;                ///< the selected fit
    std::vector<Map<T>> amaps;  ///< every exact solution found
    bool poisson_fallback;      ///< true when no exact AMAP(2) exists (or SCV == 1)
};

/**
 * Fit an AMAP(2) to (M1, M2, M3, GAMMA). cvtol is the tolerance on
 * |M2 - 2 M1^2| below which the Poisson short-circuit fires (MATLAB 1e-6).
 */
template <class T>
Amap2FitGammaResult<T> amap2_fit_gamma(const T& M1, const T& M2, const T& M3, const T& GAMMA,
                                       const T& cvtol) {
    static_assert(num_traits<T>::has_transcendental,
                  "amap2_fit_gamma requires transcendental arithmetic");
    const T two = num_traits<T>::from_int(2);
    Amap2FitGammaResult<T> r;
    r.poisson_fallback = false;

    if (num_abs(T(M2 - two * M1 * M1)) < cvtol) {
        r.amap = map_exponential_mean(M1);
        r.amaps.push_back(r.amap);
        r.poisson_fallback = true;
        return r;
    }

    std::vector<Map<T>> all = amap2_fitall_gamma(M1, M2, M3, GAMMA);
    for (std::size_t j = 0; j < all.size(); ++j) all[j] = map_normalize(all[j]);
    r.amaps = all;

    if (r.amaps.empty()) {
        r.amap = map_exponential_mean(M1);
        r.amaps.push_back(r.amap);
        r.poisson_fallback = true;
        return r;
    }
    r.amap = r.amaps.front();
    return r;
}

/** amap2_fit_gamma with the MATLAB default cvtol = 1e-6. */
template <class T>
Amap2FitGammaResult<T> amap2_fit_gamma(const T& M1, const T& M2, const T& M3, const T& GAMMA) {
    return amap2_fit_gamma(M1, M2, M3, GAMMA, T(num_traits<T>::from_double(1e-6)));
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_AMAP2_FIT_GAMMA_H
