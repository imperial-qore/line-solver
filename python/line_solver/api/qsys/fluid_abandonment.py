"""
Fluid model of a multiserver queue with customer abandonment.

Native Python twin of matlab/src/api/qsys/qsys_ggisgi_fluid.m: the steady state
of the G/GI/s+GI fluid model of W. Whitt (2006), Operations Research 54(1),
37-54, Theorem 3.1 and Corollary 3.2.
"""

from typing import Any, Callable, Dict, Optional, Sequence

import numpy as np


def _inv_ccdf(ccdf: Callable[[float], float], target: float, tol: float,
              max_time: Optional[float]) -> float:
    """
    Smallest w with F^c(w) = target, found by doubling then bisection. F^c is
    non-increasing, so the doubling either brackets the crossing or proves that
    the patience law never decays that far.
    """
    if ccdf(0.0) < target:
        raise ValueError('the patience ccdf is below 1/rho at t = 0, so it is not a ccdf')
    lo = 0.0
    if max_time is None:
        hi = 1.0
        while ccdf(hi) > target:
            hi *= 2.0
            if hi > 1e12:
                raise ValueError('the patience ccdf never falls to 1/rho, so the overloaded fluid '
                                 'model has no equilibrium: too little of the fluid is willing to '
                                 'abandon')
    else:
        hi = float(max_time)
        if ccdf(hi) > target:
            raise ValueError('the patience ccdf is still above 1/rho at maxTime')
    while hi - lo > tol * max(1.0, hi):
        mid = 0.5 * (lo + hi)
        if ccdf(mid) > target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def _integral(f: Callable[[float], float], a: float, b: float) -> float:
    """
    Composite Simpson rule on a fixed fine grid: the integrand is a ccdf, hence
    monotone and bounded, so a fixed grid is enough and is reproducible.
    """
    if b <= a:
        return 0.0
    n = 2000
    x = np.linspace(a, b, n + 1)
    y = np.array([f(float(xi)) for xi in x])
    w = np.ones(n + 1)
    w[1:-1:2] = 4.0
    w[2:-1:2] = 2.0
    return float((b - a) / (3.0 * n) * np.sum(w * y))


def qsys_ggisgi_fluid(lambda_val: float, mu: float, s: int,
                      patienceCcdf: Callable[[float], float],
                      servingCcdf: Optional[Callable[[float], float]] = None,
                      agePoints: Optional[Sequence[float]] = None,
                      tol: float = 1e-12,
                      maxTime: Optional[float] = None) -> Dict[str, Any]:
    """
    Steady state of the G/GI/s+GI fluid model.

    Scale the content by ``s`` and let ``s`` grow. Customers become quanta of
    fluid but their sojourns do not shrink, so the ages survive the limit: the
    state is the density ``b(x)`` of fluid in service of age ``x`` and the
    density ``q(x)`` of fluid waiting of age ``x``. With ``rho = lambda/(s*mu)``,

    * ``rho <= 1``: ``b(x) = rho G^c(x)``, ``q = 0``, no wait, no abandonment;
    * ``rho > 1``: ``b(x) = G^c(x)``, ``q(x) = rho F^c(x)`` on ``[0,w]``,

    the queue boundary ``w`` solving ``F^c(w) = 1/rho`` (eq. 3.6): fluid that
    survives its patience for ``w`` enters service, so the surviving fraction
    must equal the fraction ``1/rho`` the servers can absorb.

    Args:
        lambda_val: arrival rate
        mu: service rate of one server
        s: number of servers
        patienceCcdf: F^c(t) = P(patience > t)
        servingCcdf: G^c(x) = P(service > x), needed only for the in-service age
            density; defaults to the exponential of rate ``mu``
        agePoints: ages at which to return the two densities
        tol: bisection tolerance for w
        maxTime: largest age searched for w; the search grows automatically when
            this is None

    Returns:
        Dict with ``regime``, ``trafficIntensity``, ``offeredWait``,
        ``meanWait``, ``meanWaitServed``, ``meanWaitAbandon``, ``probAbandon``,
        ``meanQueueLength``, ``meanNumberInService``, ``meanNumber``,
        ``utilization``, ``throughput``, ``abandonRate``, and, when
        ``agePoints`` is given, ``serviceAgeDensity`` and ``queueAgeDensity``.

    References:
        W. Whitt (2006). Fluid models for multiserver queues with abandonments.
        Operations Research 54(1), 37-54.
    """
    if lambda_val <= 0:
        raise ValueError('The arrival rate lambda must be positive.')
    if mu <= 0:
        raise ValueError('The service rate mu must be positive.')
    if s < 1:
        raise ValueError('The number of servers s must be at least 1.')
    if not callable(patienceCcdf):
        raise ValueError('The patience ccdf must be a callable F^c(t) = P(T > t).')
    if servingCcdf is None:
        servingCcdf = lambda x, _mu=mu: float(np.exp(-_mu * x))

    rho = lambda_val / (s * mu)
    result: Dict[str, Any] = {'trafficIntensity': rho}

    if rho <= 1.0:
        # Underloaded and balanced, eq. (3.2): the queue is empty and the model
        # is the infinite-server fluid model.
        result['regime'] = 'balanced' if abs(rho - 1.0) <= np.finfo(float).eps else 'underloaded'
        w = 0.0
        mean_wait = 0.0
        prob_abandon = 0.0
        mean_wait_abandon = 0.0
    else:
        result['regime'] = 'overloaded'
        w = _inv_ccdf(patienceCcdf, 1.0 / rho, tol, maxTime)      # eq. (3.6)
        # Eq. (3.14): W = int_0^w F^c(t) dt = m_a F_e(w), over ALL fluid.
        mean_wait = _integral(patienceCcdf, 0.0, w)
        prob_abandon = 1.0 - 1.0 / rho
        # E[T | T <= w] = (W - w F^c(w)) / F(w) by parts, F^c(w) = 1/rho.
        mean_wait_abandon = (mean_wait - w / rho) / prob_abandon

    result['offeredWait'] = w
    result['meanWait'] = mean_wait
    result['meanWaitServed'] = w
    result['meanWaitAbandon'] = mean_wait_abandon
    result['probAbandon'] = prob_abandon
    result['meanQueueLength'] = lambda_val * mean_wait          # eq. (3.11), Little's law
    result['meanNumberInService'] = min(lambda_val / mu, float(s))
    result['meanNumber'] = result['meanNumberInService'] + result['meanQueueLength']
    result['utilization'] = min(rho, 1.0)
    result['throughput'] = min(lambda_val, s * mu)
    result['abandonRate'] = lambda_val - result['throughput']

    if agePoints is not None and len(agePoints) > 0:
        x = np.atleast_1d(np.asarray(agePoints, dtype=float))
        sigma = min(rho, 1.0)                                    # rate into service, per server
        result['agePoints'] = x
        result['serviceAgeDensity'] = sigma * np.array([servingCcdf(float(xi)) for xi in x])
        if rho > 1.0:
            result['queueAgeDensity'] = rho * np.array(
                [patienceCcdf(float(xi)) if xi <= w else 0.0 for xi in x])
        else:
            result['queueAgeDensity'] = np.zeros(x.size)

    return result
