/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SIM_SIM_VONNEUMANN_H
#define LINE_API_SIM_SIM_VONNEUMANN_H

/**
 * Von Neumann ratio test for randomness of a sequence.
 *
 * Port of matlab/src/api/sim/sim_vonneumann.m. The statistic is the ratio of the
 * mean square successive difference to the variance,
 *   ratio = sum_{i=1}^{b-1} (x_{i+1}-x_i)^2 / sum_{i=1}^{b} (x_i - xbar)^2,
 * with b the number of observations. Under the null hypothesis that x is i.i.d.
 * normal the ratio has mean 2 and variance 4(b-2)/((b-1)(b+1)), and (ratio-2)/sd
 * is asymptotically standard normal, so the two-sided p-value is 2(1-Phi(|z|)).
 * Serial correlation of either sign moves the ratio away from 2: positive
 * correlation shrinks the successive differences and pushes the ratio below 2,
 * negative correlation pushes it above.
 *
 * The null mean and variance above were confirmed by Monte Carlo over
 * b = 10, 16, 24, 32, 50 to within 0.3% in the reference.
 *
 * The test is TWO-SIDED and its rejection is used as a stopping rule by the
 * QUEST procedures, so alpha here is a stage significance (0.30 by default in
 * sim_fquest, decaying during warmup) and not the interval's coverage level.
 *
 * Reference: J. von Neumann, "Distribution of the Ratio of the Mean Square
 * Successive Difference to the Variance", Ann. Math. Statist. 12(4), 1941;
 * L. C. Young, "Randomness in Ordered Sequences", Ann. Math. Statist. 12, 1941.
 */

#include <cmath>
#include <cstddef>
#include <vector>

#include "line/api/sim/sim_dist.h"
#include "line/api/sim/sim_types.h"
#include "line/num/number.h"
#include "line/util/error.h"

namespace line {
namespace sim {

/**
 * Outcome of the von Neumann randomness test.
 *
 * The statistic stays in the working type T; the standardized value and the
 * p-value are doubles because they come out of a normal approximation, see the
 * arithmetic note in sim_types.h.
 */
template <class T>
struct VonNeumannResult {
    T ratio;              ///< The von Neumann ratio
    double zscore = 0.0;  ///< Standardized statistic (ratio-2)/sd
    double pvalue = 1.0;  ///< Two-sided p-value
    bool reject = false;  ///< True when pvalue < alpha, i.e. randomness is rejected
    std::size_t nobs = 0; ///< Number of observations b
};

/**
 * @param x     the sequence, at least 3 finite observations
 * @param alpha significance level in (0,1), 0.05 by default
 */
template <class T>
VonNeumannResult<T> sim_vonneumann(const std::vector<T>& x, double alpha = 0.05) {
    static_assert(num_traits<T>::has_transcendental,
                  "sim_vonneumann: the p-value is a normal tail, so exact arithmetic is refused");
    if (!(alpha > 0.0) || !(alpha < 1.0))
        throw InputError("sim_vonneumann: alpha must be a real scalar in (0,1)");

    const std::size_t b = x.size();
    if (b < 3)
        throw InputError("sim_vonneumann: at least 3 observations are required");
    for (std::size_t i = 0; i < b; ++i)
        if (!detail::num_isfinite(x[i]))
            throw InputError("sim_vonneumann: the sequence must be finite");

    T sum = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < b; ++i) sum += x[i];
    const T mean = sum / num_traits<T>::from_int(static_cast<long>(b));

    T den = num_traits<T>::from_int(0);
    for (std::size_t i = 0; i < b; ++i) {
        const T d = T(x[i] - mean);
        den += T(d * d);
    }
    if (!(den > num_traits<T>::from_int(0)))
        throw InputError("sim_vonneumann: the sequence is constant, the ratio is undefined");

    T num = num_traits<T>::from_int(0);
    for (std::size_t i = 1; i < b; ++i) {
        const T d = T(x[i] - x[i - 1]);
        num += T(d * d);
    }

    VonNeumannResult<T> r;
    r.nobs = b;
    r.ratio = T(num / den);
    const double bd = static_cast<double>(b);
    const double sd = std::sqrt(4.0 * (bd - 2.0) / ((bd - 1.0) * (bd + 1.0)));
    r.zscore = (num_traits<T>::to_double(r.ratio) - 2.0) / sd;
    r.pvalue = 2.0 * (1.0 - sim_normcdf(std::fabs(r.zscore)));
    r.reject = r.pvalue < alpha;
    return r;
}

}  // namespace sim
}  // namespace line

#endif  // LINE_API_SIM_SIM_VONNEUMANN_H
