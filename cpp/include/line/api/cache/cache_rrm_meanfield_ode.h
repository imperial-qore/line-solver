/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_CACHE_RRM_MEANFIELD_ODE_H
#define LINE_API_CACHE_RRM_MEANFIELD_ODE_H

/**
 * Mean-field drift of the RANDOM(m) multi-list cache.
 *
 * Templated port of matlab/src/api/cache/cache_rrm_meanfield_ode.m,
 * cross-checked against jar/src/main/java/jline/api/cache/
 * Cache_rrm_meanfield_ode.java.
 *
 * The state x(k,s) is the probability that item k occupies list s, with s = 0
 * meaning "not cached". A request for item k while it is in list s promotes it
 * to list s+1 and demotes a uniformly chosen occupant of list s+1; hence the
 * drift
 *
 *   dx(k,s)/dt = lambda(k) x(k,s-1)
 *              - sum_j lambda(j)/m(s) x(j,s-1) x(k,s)
 *              + sum_j lambda(j)/m(s+1) x(j,s) x(k,s+1) - lambda(k) x(k,s)
 *
 * for 1 <= s < h, without the last two terms at s = h (an item already in the
 * top list stays there), and dx(k,0)/dt = -sum_{s>=1} dx(k,s)/dt so that each
 * item's occupancies stay on the simplex.
 *
 * Bilinear in x with rational coefficients, so the drift is a field
 * expression: no transcendental requirement, and the exact instantiation
 * evaluates it without rounding. Note this is the right-hand side only. The
 * steady state (matlab's cache_rrm_meanfield script) needs a stiff ODE
 * integrator, which is out of scope for this header.
 *
 * The state vector is flat in MATLAB's column-major reshape order:
 * x[k + s*n] is x(k,s), k = 0..n-1 the item and s = 0..h the level.
 */

#include <cstddef>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace cache {

/**
 * @param x      (n*(h+1)) occupancies, x[k + s*n] = x(item k, level s)
 * @param lambda (n) per-item request rates
 * @param m      (h) list capacities
 * @return the drift, same layout as x
 */
template <class T>
std::vector<T> cache_rrm_meanfield_ode(const std::vector<T>& x, const std::vector<T>& lambda,
                                       const std::vector<int>& m) {
    const std::size_t n = lambda.size();
    const std::size_t h = m.size();
    if (x.size() != n * (h + 1))
        throw InputError("cache_rrm_meanfield_ode: state vector has the wrong length");

    const T zero = num_traits<T>::from_int(0);
    std::vector<T> dxdt(n * (h + 1), zero);

    for (std::size_t k = 0; k < n; ++k) {
        for (std::size_t s = 1; s <= h; ++s) {
            if (m[s - 1] == 0)
                throw InputError("cache_rrm_meanfield_ode: a list has zero capacity");
            const T ms = num_traits<T>::from_int(static_cast<long>(m[s - 1]));

            // outflow of item k from level s caused by promotions into level s
            T sum1 = zero;
            for (std::size_t j = 0; j < n; ++j)
                sum1 += lambda[j] / ms * x[j + (s - 1) * n] * x[k + s * n];

            // inflow from level s+1 (demotions), minus k's own promotion out
            T sum2 = zero;
            if (s < h) {
                if (m[s] == 0) throw InputError("cache_rrm_meanfield_ode: a list has zero capacity");
                const T ms1 = num_traits<T>::from_int(static_cast<long>(m[s]));
                for (std::size_t j = 0; j < n; ++j)
                    sum2 += lambda[j] / ms1 * x[j + s * n] * x[k + (s + 1) * n];
                sum2 -= lambda[k] * x[k + s * n];
            }

            dxdt[k + s * n] = lambda[k] * x[k + (s - 1) * n] - sum1 + sum2;
        }
        T acc = zero;
        for (std::size_t s = 1; s <= h; ++s) acc += dxdt[k + s * n];
        dxdt[k] = -acc;
    }
    return dxdt;
}

}  // namespace cache
}  // namespace line

#endif  // LINE_API_CACHE_RRM_MEANFIELD_ODE_H
