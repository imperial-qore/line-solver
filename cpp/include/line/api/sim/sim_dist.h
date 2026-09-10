/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SIM_SIM_DIST_H
#define LINE_API_SIM_SIM_DIST_H

/**
 * Normal and Student t quantiles used by the output-analysis routines.
 *
 * Port of matlab/src/api/sim/sim_normcdf.m, sim_norminv.m and sim_tinv.m, which are
 * grouped here as the JAR groups them in jline.api.sim.SimDist and native Python
 * in line_solver.api.sim.dist. The four codebases reach the same values through
 * four different libraries -- MATLAB through erfc/erfinv/betaincinv to avoid a
 * toolbox dependency, the JAR through commons-math3, Python through SciPy and
 * this port through Boost.Math -- and agree to within 1e-10.
 *
 * These are deliberately plain doubles rather than templates on T. A quantile
 * order p and a significance level alpha are doubles at every call site in this
 * family, so the quantile they determine carries double information and no
 * more; see the arithmetic note in sim_types.h.
 *
 * The t quantile follows the MATLAB identity rather than a distribution object,
 * so the branch structure is comparable line by line with the reference:
 *   P(|T| > t) = betainc(nu/(nu+t^2), nu/2, 1/2),
 * inverted for the two-sided tail 2(1-p) and mapped back with
 *   t = sqrt(nu (1-z)/z),  z = betaincinv(2(1-p), nu/2, 1/2).
 * MATLAB's betaincinv(y, a, b) is Boost's ibeta_inv(a, b, y): both invert the
 * REGULARIZED incomplete beta, and getting the argument order wrong here is a
 * silently plausible wrong number rather than an error.
 */

#include <cmath>
#include <limits>

#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/erf.hpp>

#include "line/util/error.h"

namespace line {
namespace sim {

/**
 * Standard normal cumulative distribution function.
 *
 * @param z the argument
 * @return Phi(z)
 */
inline double sim_normcdf(double z) {
    return 0.5 * boost::math::erfc(-z / std::sqrt(2.0));
}

/**
 * Standard normal quantile function.
 *
 * @param p probability in [0,1]
 * @return Phi^{-1}(p), infinite at the endpoints
 */
inline double sim_norminv(double p) {
    if (!(p >= 0.0) || !(p <= 1.0))
        throw InputError("sim_norminv: the probability must lie in [0,1]");
    if (p == 0.0) return -std::numeric_limits<double>::infinity();
    if (p == 1.0) return std::numeric_limits<double>::infinity();
    return std::sqrt(2.0) * boost::math::erf_inv(2.0 * p - 1.0);
}

/**
 * Quantile function of Student's t distribution.
 *
 * @param p  probability in [0,1]
 * @param nu degrees of freedom, positive
 * @return the p-quantile of t with nu degrees of freedom
 */
inline double sim_tinv(double p, double nu) {
    if (!(nu > 0.0)) throw InputError("sim_tinv: nu must be a positive real scalar");
    if (!(p >= 0.0) || !(p <= 1.0))
        throw InputError("sim_tinv: the probability must lie in [0,1]");
    if (p == 0.5) return 0.0;
    if (p <= 0.0) return -std::numeric_limits<double>::infinity();
    if (p >= 1.0) return std::numeric_limits<double>::infinity();

    // reflect the lower half onto the upper half, the law is symmetric
    const bool flip = p < 0.5;
    const double pu = flip ? 1.0 - p : p;
    const double z = boost::math::ibeta_inv(0.5 * nu, 0.5, 2.0 * (1.0 - pu));
    const double t = std::sqrt(nu * (1.0 - z) / z);
    return flip ? -t : t;
}

}  // namespace sim
}  // namespace line

#endif  // LINE_API_SIM_SIM_DIST_H
