/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_DA_DA_FPI_H
#define LINE_API_DA_DA_FPI_H

/**
 * Damped fixed-point iteration, the shared driver of the decomposition
 * algorithms.
 *
 * Templated port of matlab/src/api/da/da_fpi.m. The iteration is
 *   x_{k+1} = (1 - omega) x_ref + omega f(x_k, k)
 * stopping when the configured norm of the increment falls below iter_tol, or
 * after iter_max steps. MATLAB passes the iteration function as a handle
 * returning both the new iterate and the reference point it should be damped
 * against; the port takes a std::function with the same contract, so a caller
 * whose reference differs from its input (as in the Erlang fixed point) is
 * expressible without special-casing.
 *
 * THE INCREMENT NORM FOLLOWS MATLAB'S max(), which OMITS NaN and answers NaN
 * only when every entry is one. The first form of this port skipped NaN entries
 * but left `delta` at its initial 0, so an all-NaN increment -- a diverged
 * iterate -- reported CONVERGENCE, where MATLAB's `NaN < iter_tol` is false and
 * the loop continues to iter_max. Corrected here; a caller that wants to stop on
 * a divergence asks for it with `nanstop`.
 *
 * A tolerance-driven loop is inexact by construction, whatever the arithmetic:
 * the answer is the fixed point only to within iter_tol. The static_assert
 * records that, so nobody instantiates it at exact arithmetic expecting an
 * exact fixed point.
 */

#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <utility>
#include <vector>

#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace da {

/** Options mirroring the fields MATLAB reads off the options struct. */
struct FpiOptions {
    std::size_t iter_max = 10000;
    double iter_tol = 1e-8;
    double damping = 1.0;    ///< omega; 1 means no damping
    std::size_t miniter = 1; ///< iterations before the stopping test applies
    bool nanstop = false;    ///< stop when the increment norm is not finite
    /**
     * `config.da_norm`, the increment norm. MATLAB passes a function handle;
     * the two handles the reference actually installs are both RELATIVE --
     * `max(|xn-xr|./xr)` in solver_mam.m and `max(|xn-xr|./(xr+FineTol))` in
     * solver_mam_basic_mmap_inner.m -- so they are expressed here as a flag and
     * the denominator's offset rather than as a std::function, which would have
     * to be templated on T and would change the type of every existing caller's
     * options object.
     *
     * FALSE keeps the default absolute max-norm. TRUE with `relative_eps = 0`
     * reproduces the first handle, including its division by a zero reference
     * (the reference divides by `xr` unguarded, and an all-zero start therefore
     * yields NaN on the first sweep; `miniter` is what keeps that from stopping
     * the loop, exactly as in MATLAB).
     */
    bool relative_norm = false;
    double relative_eps = 0.0;
};

template <class T>
struct FpiResult {
    std::vector<T> x;
    std::size_t iterations = 0;
    bool converged = false;
};

/**
 * @param iterfun (x, iteration) -> (xnew, xref); xref is the point the damping
 *                and the increment norm are taken against
 * @param x0      initial iterate
 * @param options fixed-point options (tolerance, iteration cap, damping)
 */
template <class T>
FpiResult<T> da_fpi(
    const std::function<std::pair<std::vector<T>, std::vector<T>>(const std::vector<T>&, std::size_t)>&
        iterfun,
    const std::vector<T>& x0, const FpiOptions& options = FpiOptions()) {
    static_assert(num_traits<T>::has_transcendental,
                  "da_fpi requires transcendental arithmetic: it stops on a tolerance, so its "
                  "result is the fixed point only to within iter_tol whatever the arithmetic");
    if (x0.empty()) throw InputError("da_fpi: empty initial iterate");
    const T omega = num_traits<T>::from_double(options.damping);
    const T one = num_traits<T>::from_int(1);

    FpiResult<T> r;
    r.x = x0;
    for (std::size_t it = 1; it <= options.iter_max; ++it) {
        r.iterations = it;
        std::pair<std::vector<T>, std::vector<T>> step = iterfun(r.x, it);
        std::vector<T>& xnew = step.first;
        const std::vector<T>& xref = step.second;
        if (xnew.size() != r.x.size() || xref.size() != r.x.size())
            throw InputError("da_fpi: the iteration function changed the vector length");

        if (options.damping != 1.0)
            for (std::size_t i = 0; i < xnew.size(); ++i)
                xnew[i] = (one - omega) * xref[i] + omega * xnew[i];

        // MATLAB's max() OMITS NaN and answers NaN only when every entry is one,
        // which is what lets the relative norm survive a zero reference entry on
        // the first sweep instead of stopping the loop there.
        double delta = 0.0;
        bool anynum = false;
        for (std::size_t i = 0; i < xnew.size(); ++i) {
            double d = std::fabs(num_traits<T>::to_double(T(xnew[i] - xref[i])));
            if (options.relative_norm)
                d /= num_traits<T>::to_double(xref[i]) + options.relative_eps;
            if (std::isnan(d)) continue;
            if (!anynum || d > delta) delta = d;
            anynum = true;
        }
        if (!anynum && !xnew.empty()) delta = std::numeric_limits<double>::quiet_NaN();
        r.x = xnew;

        if (it >= options.miniter) {
            if (delta < options.iter_tol) {
                r.converged = true;
                break;
            }
            if (options.nanstop && !std::isfinite(delta)) break;
        }
    }
    return r;
}

}  // namespace da
}  // namespace line

#endif  // LINE_API_DA_DA_FPI_H
