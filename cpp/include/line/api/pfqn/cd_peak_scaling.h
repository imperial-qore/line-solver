/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_PFQN_CD_PEAK_SCALING_H
#define LINE_API_PFQN_CD_PEAK_SCALING_H

/**
 * Peak of a class-dependence handle over the reachable population lattice.
 *
 * Templated port of matlab/src/api/pfqn/cd_peak_scaling.m. The handle beta(n)
 * returns either a scalar shared by every class or a per-class vector; the
 * peak is taken over both the lattice states 0 <= n <= NK and the classes,
 * because utilization is a per-station quantity and the whole station shares
 * one normalizer, exactly as max(lldscaling(ist,:)) does for the
 * load-dependent case. That single normalizer is what every solver divides by
 * when reporting U = T S / bmax at a station with limited class dependence.
 *
 * NON-FINITE VALUES are skipped, as in the reference (v = v(isfinite(v))): a
 * class-dependence handle may legitimately return Inf for an unreachable
 * composition, and letting that become the normalizer would zero every
 * utilization at the station. A state whose every entry is non-finite
 * contributes nothing, and if that is every state the peak stays at zero.
 *
 * ARITHMETIC. A maximum over evaluations of the caller's handle: no
 * transcendental function is applied, so the routine is EXACT in rational
 * arithmetic and is deliberately left ungated. Whether the result is exact
 * depends only on the handle.
 */

#include <cmath>
#include <cstddef>
#include <functional>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/population.h"

namespace line {
namespace pfqn {

/**
 * @param beta class-dependence handle, evaluated on a per-class count vector
 * @param NK   (R) per-class population bound
 */
template <class T>
T cd_peak_scaling(const std::function<std::vector<T>(const std::vector<int>&)>& beta,
                  const std::vector<int>& NK) {
    if (NK.empty()) throw InputError("cd_peak_scaling: empty population vector");
    for (int v : NK)
        if (v < 0) throw InputError("cd_peak_scaling: negative population");
    T bmax = num_traits<T>::from_int(0);
    std::vector<int> n(NK.size(), 0);
    bool more = true;
    while (more) {
        int tot = 0;
        for (int v : n) tot += v;
        if (tot > 0) {
            const std::vector<T> v = beta(n);
            for (const T& x : v) {
                if (!std::isfinite(num_traits<T>::to_double(x))) continue;
                if (x > bmax) bmax = x;
            }
        }
        more = next_pop(n, NK);
    }
    return bmax;
}

}  // namespace pfqn
}  // namespace line

#endif  // LINE_API_PFQN_CD_PEAK_SCALING_H
