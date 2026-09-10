"""Stochastic network calculus: tail bounds, quantiles and mean bounds.

Native port of matlab/src/api/snc/snc_bound_*.m, snc_perc_*.m, snc_mean_*.m and
snc_thetaopt.m, cross-checked against jline.api.snc.

Every bound holds for each theta > 0 for which the arrival MGF is finite and the
station is stable, so what is reported is the infimum over theta, computed by
:func:`snc_thetaopt`. Arrival and service enter as CALLABLES of theta returning
``(sigma, rho)``; see :mod:`line_solver.api.snc.envelopes`.

Reference: M. Fidler, A. Rizk, "A Guide to the Stochastic Network Calculus",
IEEE Communications Surveys and Tutorials 17(1), 92-105, 2015, Sec. IV-B.
"""

import math

import numpy as np
from scipy.optimize import minimize_scalar

__all__ = [
    'snc_thetaopt',
    'snc_bound_backlog', 'snc_bound_delay',
    'snc_perc_backlog', 'snc_perc_delay',
    'snc_mean_delay', 'snc_mean_backlog',
]

_INFEASIBLE = 1e300


def snc_thetaopt(fun, thetamax=1e3):
    """Minimize a Chernoff bound over the free parameter theta.

    ``fun`` is evaluated on a logarithmic grid, non-finite values (a diverging
    MGF, an unstable leftover rate) are discarded, and the best grid point is
    refined by a bounded scalar minimization in log10(theta). The two-stage
    search is used because the feasible set is an interval whose endpoints are
    not known in closed form once envelopes are composed, and an unguarded local
    search steps into the infeasible region and terminates there.

    :param fun: callable theta -> scalar objective to minimize
    :param thetamax: upper end of the search range
    :return: ``(value, theta)``; ``(inf, nan)`` if no feasible theta exists
    """
    if thetamax <= 0:
        raise ValueError("snc_thetaopt: thetamax must be positive, got %g." % thetamax)

    def safeval(theta):
        try:
            v = float(fun(theta))
        except (ValueError, ZeroDivisionError, OverflowError):
            return _INFEASIBLE
        if not math.isfinite(v):
            return _INFEASIBLE
        return v

    grid = np.logspace(-6.0, math.log10(thetamax), 600)
    fval = np.array([safeval(t) for t in grid])
    imin = int(np.argmin(fval))
    val = float(fval[imin])
    if val >= 1e299:
        return float('inf'), float('nan')
    theta = float(grid[imin])

    lo = float(grid[max(imin - 1, 0)])
    hi = float(grid[min(imin + 1, grid.size - 1)])
    if hi > lo:
        res = minimize_scalar(lambda x: safeval(10.0 ** x),
                              bounds=(math.log10(lo), math.log10(hi)),
                              method='bounded', options={'xatol': 1e-10})
        if res.success and float(res.fun) < val:
            val = float(res.fun)
            theta = 10.0 ** float(res.x)
    return val, theta


def _pair(arv, srv, theta):
    """Both envelopes at one theta, or None when the composition is infeasible."""
    sA, rA = arv(theta)
    sS, rS = srv(theta)
    if not (math.isfinite(sA) and math.isfinite(sS)
            and math.isfinite(rA) and math.isfinite(rS)) or rS <= rA:
        return None
    return sA, rA, sS, rS


def snc_bound_backlog(arv, srv, b, thetamax=1e3):
    """Upper bound on ``P{B > b}``, minimized over theta and clipped at 1.

    ``P{B(t) > b} <= exp(-theta*(b-sigmaA-sigmaS))/(1-exp(-theta*(rhoS-rhoA)))``,
    the union bound over the start of the backlogged period summed as a
    geometric series on the unit-slot time axis. It is an upper bound on the
    tail, never an estimate of it: the decay rate is asymptotically exact and
    the prefactor is loose.

    :return: ``(eps, theta)``
    """
    if b < 0:
        raise ValueError("snc_bound_backlog: b must be nonnegative, got %g." % b)

    def obj(theta):
        p = _pair(arv, srv, theta)
        if p is None:
            return float('inf')
        sA, rA, sS, rS = p
        return math.exp(-theta * (b - sA - sS)) / (1.0 - math.exp(-theta * (rS - rA)))

    eps, theta = snc_thetaopt(obj, thetamax)
    if not math.isfinite(eps) or eps > 1.0:
        eps = 1.0  # no feasible theta, or the bound is vacuous at this level
    return eps, theta


def snc_bound_delay(arv, srv, d, thetamax=1e3):
    """Upper bound on ``P{D > d}``, minimized over theta and clipped at 1.

    ``P{D(t) > d} <= exp(-theta*(rhoS*d-sigmaA-sigmaS))/(1-exp(-theta*(rhoS-rhoA)))``,
    the horizontal rather than vertical deviation between the arrival and
    service envelopes.

    :return: ``(eps, theta)``
    """
    if d < 0:
        raise ValueError("snc_bound_delay: d must be nonnegative, got %g." % d)

    def obj(theta):
        p = _pair(arv, srv, theta)
        if p is None:
            return float('inf')
        sA, rA, sS, rS = p
        return math.exp(-theta * (rS * d - sA - sS)) / (1.0 - math.exp(-theta * (rS - rA)))

    eps, theta = snc_thetaopt(obj, thetamax)
    if not math.isfinite(eps) or eps > 1.0:
        eps = 1.0
    return eps, theta


def snc_perc_backlog(arv, srv, eps, thetamax=1e3):
    """Backlog quantile at violation probability ``eps``.

    Inverts :func:`snc_bound_backlog` in b at fixed theta and re-optimizes:
    ``b(theta) = sigmaA + sigmaS - log(eps*(1-exp(-theta*(rhoS-rhoA))))/theta``.
    The minimizing theta differs from the one of the forward bound at a given
    level, which is why the inversion is done in closed form.

    :return: ``(b, theta)``
    """
    if not (0.0 < eps < 1.0):
        raise ValueError("snc_perc_backlog: eps must lie in (0,1), got %g." % eps)

    def obj(theta):
        p = _pair(arv, srv, theta)
        if p is None:
            return float('inf')
        sA, rA, sS, rS = p
        return sA + sS - math.log(eps * (1.0 - math.exp(-theta * (rS - rA)))) / theta

    b, theta = snc_thetaopt(obj, thetamax)
    return max(b, 0.0), theta


def snc_perc_delay(arv, srv, eps, thetamax=1e3):
    """Delay quantile at violation probability ``eps``.

    ``d(theta) = (sigmaA + sigmaS - log(eps*(1-exp(-theta*(rhoS-rhoA))))/theta)/rhoS``,
    minimized over the feasible thetas. This is the deliverable of the family: a
    statistical delay guarantee, the quantity an SLA is written against, as
    opposed to the mean delay returned by the queueing-theoretic solvers.

    :return: ``(d, theta)``
    """
    if not (0.0 < eps < 1.0):
        raise ValueError("snc_perc_delay: eps must lie in (0,1), got %g." % eps)

    def obj(theta):
        p = _pair(arv, srv, theta)
        if p is None:
            return float('inf')
        sA, rA, sS, rS = p
        return (sA + sS - math.log(eps * (1.0 - math.exp(-theta * (rS - rA)))) / theta) / rS

    d, theta = snc_thetaopt(obj, thetamax)
    return max(d, 0.0), theta


def _mean_from_tail(logK, a):
    """Closed-form integral of ``min(1, K*exp(-a*x))`` over the positive axis."""
    if logK >= 0.0:
        return (logK + 1.0) / a
    return math.exp(logK) / a


def snc_mean_delay(arv, srv, thetamax=1e3):
    """Upper bound on ``E[D]``, from integrating the delay tail bound.

    For a nonnegative delay ``E[D] = int_0^inf P{D>d} dd``, and at fixed theta
    the clipped integral of ``K*exp(-theta*rhoS*d)`` is available in CLOSED
    FORM, so no quadrature is involved and the result is a bound rather than a
    bound plus a discretization error.

    IT IS A LOOSE MEAN BOUND AND THAT IS INHERENT: on the M/M/1 read in job
    units it returns 2.4x the exact ``1/(mu-lambda)`` at rho = 0.1 and 10.4x at
    rho = 0.95, because an integral over the whole axis is dominated by the
    prefactor rather than by the decay rate. Use :func:`snc_perc_delay` when the
    quantile is what matters.

    :return: ``(ED, theta)``
    """
    def obj(theta):
        p = _pair(arv, srv, theta)
        if p is None:
            return float('inf')
        sA, rA, sS, rS = p
        logK = theta * (sA + sS) - math.log(1.0 - math.exp(-theta * (rS - rA)))
        return _mean_from_tail(logK, theta * rS)

    return snc_thetaopt(obj, thetamax)


def snc_mean_backlog(arv, srv, thetamax=1e3):
    """Upper bound on ``E[B]``, from integrating the backlog tail bound.

    The counterpart of :func:`snc_mean_delay`, with ``a = theta``. The unit is
    the unit of the envelopes: jobs when the pair is
    :func:`snc_env_poisson` with :func:`snc_srv_exp`, units of work when it is
    :func:`snc_env_cpoisson` with :func:`snc_srv_rate`.

    :return: ``(EB, theta)``
    """
    def obj(theta):
        p = _pair(arv, srv, theta)
        if p is None:
            return float('inf')
        sA, rA, sS, rS = p
        logK = theta * (sA + sS) - math.log(1.0 - math.exp(-theta * (rS - rA)))
        return _mean_from_tail(logK, theta)

    return snc_thetaopt(obj, thetamax)
