"""
Diffusion approximation for the G/GI/n/m queue.

Native Python twin of matlab/src/api/qsys/qsys_ggnm_diffusion.m, implementing
W. Whitt (2004), A diffusion approximation for the G/GI/n/m queue, Operations
Research 52(6), 922-941.
"""

from math import erfc, exp, log, pi, sqrt
from typing import Any, Callable, Dict, Optional

import numpy as np


def _phi(x: float) -> float:
    """Standard normal density."""
    return exp(-x * x / 2.0) / sqrt(2.0 * pi)


def _Phi(x: float) -> float:
    """Standard normal cdf, through erfc so no statistics package is needed."""
    return erfc(-x / sqrt(2.0)) / 2.0


def _peakedness_weight(serviceCcdf: Callable[[float], float], ES: float, tol: float,
                       panels: int) -> float:
    """
    omega_G = int G^c(x)^2 dx / int G^c(x) dx of eq. (1.7), by Simpson on a grid
    cut where the ccdf is negligible. The denominator is E[S], so only the
    numerator is actually integrated.
    """
    hi = 1.0
    while serviceCcdf(hi) > tol:
        hi *= 2.0
        if hi > 1e12:
            raise ValueError('the service ccdf does not decay, so its peakedness is undefined')
    x = np.linspace(0.0, hi, panels + 1)
    y = np.array([float(serviceCcdf(float(xx))) ** 2 for xx in x])
    w = np.ones(panels + 1)
    w[1:-1:2] = 4.0
    w[2:-1:2] = 2.0
    num = float(hi / (3.0 * panels) * np.sum(w * y))
    return num / ES


def qsys_ggnm_diffusion(lambda_val: float, mu: float, n: int, m: float, ca: float, cs: float,
                        serviceCcdf: Optional[Callable[[float], float]] = None,
                        tol: float = 1e-12, panels: int = 4000) -> Dict[str, Any]:
    """
    Diffusion approximation for the G/GI/n/m queue.

    A general arrival process characterized by its rate and its variability
    parameter ``ca^2``, iid general service times of mean ``1/mu`` and SCV
    ``cs^2``, ``n`` servers and ``m`` extra waiting spaces.

    THE APPROXIMATION IS ONE DIFFUSION WITH TWO REGIONS. Below the staffing
    level the queue behaves like an infinite-server system, whose limit is
    NORMAL with variance-to-mean ratio the ASYMPTOTIC PEAKEDNESS

        z = 1 + (ca^2 - 1) omega_G,   omega_G = int G^c(x)^2 dx / int G^c(x) dx

    (eqs. 1.6-1.7); above it the queue behaves like a single-server queue, whose
    limit is EXPONENTIAL with variability ``v = (ca^2 + cs^2)/2`` (eq. 3.7). The
    steady-state law is therefore a normal piece spliced to an exponential piece,
    and every measure below is an integral of that density (eq. 3.14).

    WHAT z SAYS. The service-time distribution enters the delay probability ONLY
    through ``omega_G``, which is 1 for deterministic service, 1/2 for
    exponential, and falls toward 0 as service gets more variable. So when
    ``ca^2 = 1`` the delay probability does not depend on the service law at all
    (``z = 1``), which is the long-standing M/GI/n approximation by M/M/n; away
    from ``ca^2 = 1`` it does, and this quantifies how much.

    The delay probability is ``alpha(beta/sqrt(z))`` with the Halfin-Whitt
    function ``alpha`` when ``m`` is infinite (eq. 3.10), so this generalizes
    :func:`qsys_mmk_qed`.

    Args:
        lambda_val: arrival rate
        mu: service rate of one server
        n: number of servers
        m: extra waiting spaces; ``float('inf')`` for an unbounded queue
        ca: coefficient of variation of the interarrival time
        cs: coefficient of variation of the service time
        serviceCcdf: G^c(x) = P(S > x); the exponential of rate ``mu`` by default
        tol: service-tail cut for the peakedness integral
        panels: Simpson panels for it

    Returns:
        Dict with ``beta`` (the QED server slack), ``gamma`` (the scaled waiting
        room), ``peakedness`` (z), ``variability`` (v), ``probDelay``,
        ``probBlock``, ``meanQueueLength`` (customers waiting),
        ``meanNumber`` (in system), ``meanWait``, ``utilization`` and
        ``throughput``.

    References:
        W. Whitt (2004). A diffusion approximation for the G/GI/n/m queue.
        Operations Research 52(6), 922-941.
    """
    if lambda_val <= 0 or mu <= 0:
        raise ValueError('The arrival and service rates must be positive.')
    n = int(round(n))
    if n < 1:
        raise ValueError('The number of servers n must be at least 1.')
    if m < 0:
        raise ValueError('The number of extra waiting spaces m must be non-negative.')

    ca2, cs2 = ca ** 2, cs ** 2
    ES = 1.0 / mu
    rho = lambda_val / (n * mu)
    beta = sqrt(n) * (1.0 - rho)                       # eq. (0.1)
    gamma = float('inf') if not np.isfinite(m) else m / sqrt(n)     # eq. (0.3)

    if serviceCcdf is None:
        omega = 0.5                                    # exponential service
    else:
        omega = _peakedness_weight(serviceCcdf, ES, tol, panels)
    z = 1.0 + (ca2 - 1.0) * omega                      # eq. (1.6), asymptotic peakedness
    if z <= 0:
        raise ValueError('the asymptotic peakedness came out non-positive; check ca and the '
                         'service ccdf')
    v = (ca2 + cs2) / 2.0                              # eq. (3.7) with the weight w = 1
    b = beta / sqrt(z)
    r = beta / v                                       # rate of the exponential piece

    # Mass on the exponential piece. The tail factor is 1 - exp(-r*gamma), which
    # is negative together with r when the queue is overloaded, so the ratio
    # stays positive and alpha stays in (0,1) on both sides of beta = 0.
    if np.isfinite(gamma):
        tail = -np.expm1(-r * gamma)
    else:
        tail = 1.0
    if abs(r) < 1e-14:
        # beta = 0: the exponential piece degenerates to a uniform on [0,gamma].
        if not np.isfinite(gamma):
            raise ValueError('with beta = 0 the queue needs a finite waiting room to be stable')
        alpha = 1.0 / (1.0 + _Phi(b) / (_phi(b) * gamma / sqrt(z)))
        meanAbove = gamma / 2.0
        densityAtTop = alpha / gamma
    else:
        alpha = 1.0 / (1.0 + b * _Phi(b) / (_phi(b) * tail))
        # Mean of the truncated exponential on [0,gamma] with rate r.
        if np.isfinite(gamma):
            e = exp(-r * gamma)
            meanAbove = (1.0 / r - (gamma + 1.0 / r) * e) / tail
            densityAtTop = alpha * r * e / tail
        else:
            meanAbove = 1.0 / r
            densityAtTop = 0.0

    # Mean of the normal piece, N(-beta, z) conditioned below 0.
    meanBelow = -beta - sqrt(z) * _phi(b) / _Phi(b)
    meanScaled = (1.0 - alpha) * meanBelow + alpha * meanAbove

    probDelay = alpha
    # Eq. (7.5): the loss rate of the diffusion at the upper boundary, divided by
    # the arrival rate, is the density there times v over sqrt(n).
    probBlock = densityAtTop * v / sqrt(n) if np.isfinite(gamma) else 0.0
    probBlock = min(max(probBlock, 0.0), 1.0)
    meanQueue = sqrt(n) * alpha * meanAbove
    meanNumber = n + sqrt(n) * meanScaled
    throughput = lambda_val * (1.0 - probBlock)
    return {
        'beta': beta,
        'gamma': gamma,
        'peakedness': z,
        'peakednessWeight': omega,
        'variability': v,
        'probDelay': probDelay,
        'probBlock': probBlock,
        'meanQueueLength': meanQueue,
        'meanNumber': meanNumber,
        'meanWait': meanQueue / throughput if throughput > 0 else 0.0,
        'utilization': min(rho, 1.0),
        'throughput': throughput,
        'trafficIntensity': rho,
    }
