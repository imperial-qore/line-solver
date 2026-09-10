/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_CACHE_CACHE_RRM_MEANFIELD_H
#define LINE_API_CACHE_CACHE_RRM_MEANFIELD_H

/**
 * Steady state of the RANDOM(m) multi-list mean field.
 *
 * Port of matlab/src/api/cache/cache_rrm_meanfield.m, which is a SCRIPT rather
 * than a function: it fixes n = 7 items, m = [1 1 3], a fixed popularity
 * vector and a fixed initial condition, integrates the drift of
 * cache_rrm_meanfield_ode.m with ode23s over [0,10000] and prints the terminal
 * occupancy, the miss rate lambda*x(:,1) and the miss ratio. This header is
 * that computation with the data as arguments.
 *
 * The right-hand side was already ported (cache_rrm_meanfield_ode.h) and noted
 * there that the steady state needed a stiff integrator that the port did not
 * have. It now does: line/util/ode.h. The Jacobian is not written out here --
 * the drift is bilinear, so a numeric central-difference Jacobian is accurate
 * to about eps^(2/3) and the integrator forms it itself; that is the deliberate
 * choice, since unlike cache_miss_rmf nothing downstream needs an exact
 * Jacobian and an analytic one would be a second expression of the same
 * formula to keep in step with the first.
 *
 * The initial condition is the reference's: every item outside the cache,
 * x(k,0) = 1. Integrating from there to t = 1e4 is exactly what
 * cache_rrm_meanfield.m does.
 *
 * ARITHMETIC. Gated: the answer is the limit of a tolerance-driven
 * integration.
 */

#include <cstddef>
#include <vector>

#include "line/api/cache/cache_rrm_meanfield_ode.h"
#include "line/num/number.h"
#include "line/util/error.h"
#include "line/util/ode.h"

namespace line {
namespace cache {

/** Return value of cache_rrm_meanfield. */
template <class T>
struct CacheRrmMeanfieldResult {
    std::vector<T> x;  ///< (n*(h+1)) terminal occupancy, x[k + s*n] = x(item k, level s)
    T missrate;        ///< lambda . x(:,0)
    T missratio;       ///< missrate / sum(lambda)
};

/**
 * @param lambda (n) per-item request rates
 * @param m      (h) list capacities
 * @param tmax   integration horizon (the reference uses 1e4)
 */
template <class T>
CacheRrmMeanfieldResult<T> cache_rrm_meanfield(const std::vector<T>& lambda,
                                               const std::vector<int>& m, const T& tmax) {
    static_assert(num_traits<T>::has_transcendental,
                  "cache_rrm_meanfield requires transcendental arithmetic: the steady state is "
                  "the limit of a tolerance-driven integration of the mean-field drift");
    const std::size_t n = lambda.size();
    const std::size_t h = m.size();
    if (n == 0) throw InputError("cache_rrm_meanfield: no items");
    if (h == 0) throw InputError("cache_rrm_meanfield: no cache lists");

    const T zero = num_traits<T>::from_int(0);
    std::vector<T> x0(n * (h + 1), zero);
    for (std::size_t k = 0; k < n; ++k) x0[k] = num_traits<T>::from_int(1);

    const auto f = [&](const T& t, const std::vector<T>& x) {
        (void)t;
        return cache_rrm_meanfield_ode(x, lambda, m);
    };
    OdeOptions<T> opt;
    opt.rtol = num_traits<T>::from_double(1e-8);
    opt.atol = num_traits<T>::from_double(1e-10);
    opt.store_trajectory = false;

    CacheRrmMeanfieldResult<T> res;
    res.x = ode_rosenbrock4(f, T(num_traits<T>::from_int(0)), tmax, x0, opt).final_state();
    res.missrate = zero;
    T tot = zero;
    for (std::size_t k = 0; k < n; ++k) {
        res.missrate += lambda[k] * res.x[k];
        tot += lambda[k];
    }
    res.missratio = res.missrate / tot;
    return res;
}

/** cache_rrm_meanfield with the reference horizon tmax = 1e4. */
template <class T>
CacheRrmMeanfieldResult<T> cache_rrm_meanfield(const std::vector<T>& lambda,
                                               const std::vector<int>& m) {
    return cache_rrm_meanfield(lambda, m, T(num_traits<T>::from_int(10000)));
}

}  // namespace cache
}  // namespace line

#endif  // LINE_API_CACHE_CACHE_RRM_MEANFIELD_H
