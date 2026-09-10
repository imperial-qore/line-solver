/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_MAM_APH2_FIT_H
#define LINE_API_MAM_APH2_FIT_H

/**
 * APH(2) fit of three moments, with a fallback to adjusted moments
 * (matlab/lib/m3a/m3a/aph2/aph2_fit.m).
 *
 * Exact fitting is attempted first with aph2_fitall; if that yields nothing
 * the moments are relaxed with aph2_adjust ('simple' method) and the fit is
 * repeated. The first solution is returned, together with all of them.
 *
 * Gated on transcendental arithmetic through aph2_fitall and aph2_adjust.
 */

#include <vector>

#include "line/api/mam/aph2_adjust.h"
#include "line/api/mam/aph2_fitall.h"
#include "line/api/mam/map_moment.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace mam {

/** Result of aph2_fit. */
template <class T>
struct Aph2FitResult {
    Map<T> aph;                ///< the selected fit, APHS.front()
    std::vector<Map<T>> aphs;  ///< every feasible form found
    bool adjusted;             ///< true when the moments had to be relaxed
};

/** Fit an APH(2) to (M1, M2, M3), relaxing the moments if necessary. */
template <class T>
Aph2FitResult<T> aph2_fit(const T& M1, const T& M2, const T& M3) {
    static_assert(num_traits<T>::has_transcendental, "aph2_fit requires transcendental arithmetic");
    Aph2FitResult<T> r;
    r.adjusted = false;
    r.aphs = aph2_fitall(M1, M2, M3);
    if (r.aphs.empty()) {
        const Aph2AdjustResult<T> adj = aph2_adjust(M1, M2, M3);
        r.adjusted = true;
        r.aphs = aph2_fitall(M1, adj.M2a, adj.M3a);
        if (r.aphs.empty()) throw NumericError("aph2_fit: feasibility could not be restored");
    }
    r.aph = r.aphs.front();
    return r;
}

}  // namespace mam
}  // namespace line

#endif  // LINE_API_MAM_APH2_FIT_H
