"""
Gaussian approximations for heavily-loaded G/GI/n+GI queues.

Native Python twin of matlab/src/api/qsys/qsys_ggingi_tga.m, implementing the
truncated Gaussian approximation (TGA-G) of Y. Liu, W. Whitt and Y. Yu (2016),
Approximations for heavily-loaded G/GI/n+GI queues, Naval Research Logistics
63(3), 187-217.
"""

from math import erfc, exp, pi, sqrt
from typing import Any, Callable, Dict, Optional

import numpy as np


def _phi(x: float) -> float:
    """Standard normal density."""
    return exp(-x * x / 2.0) / sqrt(2.0 * pi)


def _Phi(x: float) -> float:
    """Standard normal cdf, through erfc so no statistics package is needed."""
    return erfc(-x / sqrt(2.0)) / 2.0


def _simpson(f: Callable[[float], float], a: float, b: float, n: int = 2000) -> float:
    """Composite Simpson rule on an even panel count."""
    if b <= a:
        return 0.0
    if n % 2 == 1:
        n += 1
    x = np.linspace(a, b, n + 1)
    y = np.array([float(f(float(xx))) for xx in x])
    w = np.ones(n + 1)
    w[1:-1:2] = 4.0
    w[2:-1:2] = 2.0
    return float((b - a) / (3.0 * n) * np.sum(w * y))


def _inv_ccdf(ccdf: Callable[[float], float], target: float, tol: float = 1e-12) -> float:
    """Smallest w with F^c(w) = target, by doubling then bisection."""
    lo, hi = 0.0, 1.0
    while ccdf(hi) > target:
        hi *= 2.0
        if hi > 1e12:
            raise ValueError('the patience ccdf never falls to 1/rho, so the overloaded model has '
                             'no fluid equilibrium')
    while hi - lo > tol * max(1.0, hi):
        mid = 0.5 * (lo + hi)
        if ccdf(mid) > target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def _truncated_moments(a: float):
    """
    Mean and variance of ``max(Z,-a)`` for a standard normal Z, the shape every
    truncated Gaussian measure below is built from.
    """
    m1 = _phi(a) - a * (1.0 - _Phi(a))
    m2 = _Phi(a) - a * _phi(a) + a * a * (1.0 - _Phi(a))
    return m1, max(0.0, m2 - m1 * m1)


def qsys_ggingi_tga(lambda_val: float, mu: float, n: int, ca: float, cs: float,
                    patienceCcdf: Callable[[float], float],
                    patiencePdf: Optional[Callable[[float], float]] = None,
                    serviceCcdf: Optional[Callable[[float], float]] = None) -> Dict[str, Any]:
    """
    Truncated Gaussian approximation (TGA-G) for the G/GI/n+GI queue.

    A general stationary arrival process of rate ``lambda_val`` and variability
    ``ca^2``, iid general service of mean ``1/mu`` and variability ``cs^2``, ``n``
    servers, unlimited waiting room, and iid general patience.

    THE APPROXIMATION IS A FLUID CENTRE PLUS A GAUSSIAN FLUCTUATION, TRUNCATED.
    In the efficiency-driven regime (``rho > 1`` held fixed as ``n`` grows) the
    fluid limit gives the centre -- all servers busy, waiting time
    ``w = F^-1(1-1/rho)``, queue ``Q = lambda int_0^w F^c`` -- and the
    many-server central limit theorem gives a NORMAL fluctuation of order
    ``sqrt(n)`` around it. Adding the two directly can produce negative queues
    and negative waits, so both are TRUNCATED at zero, which is what makes the
    formulas usable down to moderate overload; the paper reports good accuracy
    for ``rho > 1.02`` and abandonment rates below 2.

    Three independent sources of variability enter separately, which is what
    lets the exponential-service formula be generalized: the service law appears
    only as the factor ``(cs+1)rho`` in ``sigma_W^2`` (eq. 24), reducing to the
    exponential case at ``cs = 1``.

    An UNDERLOADED model (``rho <= 1``) has no queue in the limit; the number in
    system is then normal with the infinite-server variance, whose
    variance-to-mean ratio is the asymptotic peakedness of
    :func:`qsys_ggnm_diffusion`.

    Args:
        lambda_val: arrival rate
        mu: service rate of one server
        n: number of servers
        ca: coefficient of variation of the interarrival time
        cs: coefficient of variation of the service time
        patienceCcdf: F^c(x) = P(patience > x)
        patiencePdf: the patience density; differenced from the ccdf when absent
        serviceCcdf: G^c(x) = P(S > x), used only in the underloaded branch

    Returns:
        Dict with ``regime``, ``trafficIntensity``, ``fluidWait``,
        ``fluidQueueLength``, ``meanWait``, ``varWait``, ``meanQueueLength``,
        ``varQueueLength``, ``meanNumberInService``, ``meanNumber``,
        ``probDelay``, ``probAbandon``, ``sigmaW``, ``sigmaX``.

    References:
        Y. Liu, W. Whitt, Y. Yu (2016). Approximations for heavily-loaded
        G/GI/n+GI queues. Naval Research Logistics 63(3), 187-217.
    """
    if lambda_val <= 0 or mu <= 0:
        raise ValueError('The arrival and service rates must be positive.')
    n = int(round(n))
    if n < 1:
        raise ValueError('The number of servers n must be at least 1.')
    ca2, cs2 = ca ** 2, cs ** 2
    rho = lambda_val / (n * mu)
    lam_pn = lambda_val / n              # the per-server arrival rate of the reference's scaling

    if patiencePdf is None:
        h = 1e-6
        pdf = lambda x: max(0.0, (patienceCcdf(max(0.0, x - h)) - patienceCcdf(x + h)) / (2 * h))
    else:
        pdf = patiencePdf

    result: Dict[str, Any] = {'trafficIntensity': rho}
    if rho <= 1.0:
        # Underloaded: no queue in the limit, and the number in system is normal
        # with the infinite-server variance (eq. 10).
        if serviceCcdf is None:
            omega = 0.5
        else:
            num = _simpson(lambda x: float(serviceCcdf(x)) ** 2, 0.0,
                           _inv_ccdf(serviceCcdf, 1e-12))
            omega = num * mu
        varX = (lambda_val / mu) * (1.0 + (ca2 - 1.0) * omega)
        result.update({
            'regime': 'underloaded',
            'fluidWait': 0.0,
            'fluidQueueLength': 0.0,
            'meanWait': 0.0,
            'varWait': 0.0,
            'meanQueueLength': 0.0,
            'varQueueLength': 0.0,
            'meanNumberInService': lambda_val / mu,
            'meanNumber': lambda_val / mu,
            'probDelay': 0.0,
            'probAbandon': 0.0,
            'sigmaW': 0.0,
            'sigmaX': sqrt(varX),
        })
        return result

    # Overloaded: the fluid centre of Theorem 2.1(b).
    w = _inv_ccdf(patienceCcdf, 1.0 / rho)
    fw = float(pdf(w))
    if fw <= 0:
        raise ValueError('the patience density vanishes at the fluid waiting time, so the Gaussian '
                         'correction is undefined there')
    qPerServer = lam_pn * _simpson(lambda x: float(patienceCcdf(x)), 0.0, w)

    # Eq. (24): the service law enters only through the (cs+1)rho term, which is
    # 2rho for exponential service and recovers eq. (11) there.
    sigmaW2 = ((ca2 - 1.0) + (cs + 1.0) * rho) / (2.0 * mu * rho * rho * fw)
    sigmaX2 = mu * mu * sigmaW2 + lam_pn * _simpson(
        lambda x: float(patienceCcdf(x)) * (1.0 + (ca2 - 1.0) * float(patienceCcdf(x))), 0.0, w)
    sigmaW = sqrt(max(sigmaW2, 0.0))
    sigmaX = sqrt(max(sigmaX2, 0.0))

    aW = sqrt(n) * w / sigmaW                       # eq. (21)
    aX = sqrt(n) * qPerServer / sigmaX              # eq. (19)
    m1W, vW = _truncated_moments(aW)
    m1X, vX = _truncated_moments(aX)

    meanWait = w * (_Phi(aW) + _phi(aW) / aW)       # eq. (20) in moment form
    varWait = (sigmaW * sigmaW / n) * vW
    meanQueue = n * qPerServer * (_Phi(aX) + _phi(aX) / aX)
    varQueue = n * sigmaX * sigmaX * vX
    # E[B] = E[min(X_n, n)]: the servers are all busy but for the lower tail of
    # the Gaussian fluctuation.
    meanB = n - sqrt(n) * sigmaX * (_phi(aX) - aX * (1.0 - _Phi(aX)))
    probDelay = _Phi(aW)                            # eq. (22)
    # Eq. (23): a customer abandons when its patience falls short of its wait.
    probAbandon = _simpson(
        lambda x: (1.0 - _Phi(aW * (x / w - 1.0))) * float(pdf(x)), 0.0, max(w * 20.0, w + 20.0))
    probAbandon = min(max(probAbandon, 0.0), 1.0)

    result.update({
        'regime': 'overloaded',
        'fluidWait': w,
        'fluidQueueLength': n * qPerServer,
        'meanWait': meanWait,
        'varWait': varWait,
        'meanQueueLength': meanQueue,
        'varQueueLength': varQueue,
        'meanNumberInService': meanB,
        'meanNumber': meanB + meanQueue,
        'probDelay': probDelay,
        'probAbandon': probAbandon,
        'sigmaW': sigmaW,
        'sigmaX': sigmaX,
    })
    return result
