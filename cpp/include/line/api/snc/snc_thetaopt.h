/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_API_SNC_THETAOPT_H
#define LINE_API_SNC_THETAOPT_H

/**
 * Minimizes a Chernoff bound over the free parameter theta.
 *
 * Every bound in the snc domain holds for each theta > 0 for which the arrival
 * MGF is finite and the station is stable, so the reported bound is the infimum
 * over theta. The objective is evaluated on a logarithmic grid, non-finite
 * values (a diverging MGF, an unstable leftover rate) are discarded, and the
 * best grid point is refined by golden-section search in log10(theta).
 *
 * THE TWO-STAGE SEARCH IS NOT A CONVENIENCE: the feasible set is an interval
 * whose endpoints are not known in closed form once envelopes are composed, and
 * an unguarded local search steps into the infeasible region and terminates
 * there.
 *
 * Port of matlab/src/api/snc/snc_thetaopt.m, whose refinement is `fminbnd`
 * (golden section plus parabolic interpolation); the plain golden section here
 * reaches the same optimum on these smooth objectives, as the JAR port does.
 */

#include <algorithm>
#include <cmath>
#include <exception>
#include <functional>
#include <limits>
#include <vector>
#include "line/api/snc/snc_types.h"
#include "line/util/error.h"

namespace line {
namespace snc {

namespace detail {

/** Substituted for a non-finite objective, so the search can still compare it. */
constexpr double kSncInfeasible = 1e300;

inline double snc_safeval(const std::function<double(double)>& fun, double theta) {
    double v;
    try {
        v = fun(theta);
    } catch (const std::exception&) {
        return kSncInfeasible;
    }
    if (!std::isfinite(v)) return kSncInfeasible;
    return v;
}

}  // namespace detail

/**
 * @param fun      the objective, a function of theta
 * @param thetamax upper end of the search range
 * @return the minimum and its theta; (infinity, NaN) if nothing is feasible
 */
inline SncResult snc_thetaopt(const std::function<double(double)>& fun, double thetamax = 1e3) {
    if (thetamax <= 0) throw UnsupportedError("snc_thetaopt: thetamax must be positive");
    const int GRID = 600;
    const double loExp = -6.0, hiExp = std::log10(thetamax);
    std::vector<double> grid(GRID), fval(GRID);
    int imin = 0;
    for (int i = 0; i < GRID; ++i) {
        grid[i] = std::pow(10.0, loExp + (hiExp - loExp) * i / (GRID - 1.0));
        fval[i] = detail::snc_safeval(fun, grid[i]);
        if (fval[i] < fval[imin]) imin = i;
    }
    double val = fval[imin];
    if (val >= 1e299)
        return SncResult{std::numeric_limits<double>::infinity(),
                         std::numeric_limits<double>::quiet_NaN()};
    double theta = grid[imin];

    double lo = std::log10(grid[std::max(imin - 1, 0)]);
    double hi = std::log10(grid[std::min(imin + 1, GRID - 1)]);
    if (hi > lo) {
        const double invphi = (std::sqrt(5.0) - 1.0) / 2.0;
        double x1 = hi - invphi * (hi - lo), x2 = lo + invphi * (hi - lo);
        double f1 = detail::snc_safeval(fun, std::pow(10.0, x1));
        double f2 = detail::snc_safeval(fun, std::pow(10.0, x2));
        for (int it = 0; it < 200 && (hi - lo) > 1e-12; ++it) {
            if (f1 < f2) {
                hi = x2;
                x2 = x1;
                f2 = f1;
                x1 = hi - invphi * (hi - lo);
                f1 = detail::snc_safeval(fun, std::pow(10.0, x1));
            } else {
                lo = x1;
                x1 = x2;
                f1 = f2;
                x2 = lo + invphi * (hi - lo);
                f2 = detail::snc_safeval(fun, std::pow(10.0, x2));
            }
        }
        const double xopt = 0.5 * (lo + hi);
        const double vopt = detail::snc_safeval(fun, std::pow(10.0, xopt));
        if (vopt < val) {
            val = vopt;
            theta = std::pow(10.0, xopt);
        }
    }
    return SncResult{val, theta};
}

namespace detail {

/** Both envelopes at one theta, with `ok=false` when the composition is infeasible. */
struct SncPair {
    bool ok = false;
    Env a, s;
};

inline SncPair snc_pair(const Envelope& arv, const Envelope& srv, double theta) {
    SncPair p;
    p.a = arv(theta);
    p.s = srv(theta);
    p.ok = std::isfinite(p.a.sigma) && std::isfinite(p.a.rho) && std::isfinite(p.s.sigma) &&
           std::isfinite(p.s.rho) && p.s.rho > p.a.rho;
    return p;
}

}  // namespace detail

}  // namespace snc
}  // namespace line

#endif  // LINE_API_SNC_THETAOPT_H
