"""
Extremal two-moment bounds for the GI/GI/1 queue.

Native Python twin of matlab/src/api/qsys/qsys_gig1_bnds_extremal.m: the tight
lower bound, the conjectured tight upper bound of Y. Chen and W. Whitt (2020),
Queueing Systems 94, 327-356, and the classical Kingman and Daley bounds it
improves on.
"""

from math import exp

import numpy as np
from scipy.special import gammaln
from typing import Any, Dict


def _delta(rho: float) -> float:
    """
    The D/M/1 root of eq. (3.5), ``delta = exp(-(1-delta)/rho)``, in (0,1).

    ``g(delta) = delta - exp(-(1-delta)/rho)`` is negative at 0 and positive just
    below 1, where the second root ``delta = 1`` sits, so bisection on [0,1)
    finds the wanted root without ever landing on the trivial one.
    """
    lo, hi = 0.0, 1.0 - 1e-15
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if mid - exp(-(1.0 - mid) / rho) < 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def _tight(rho: float, ca2: float, cs2: float, K: int, N: int) -> float:
    """
    Algorithm 1 of the reference: the mean waiting time of the extremal model,

        E[W(F0,Gu*)] = rho ca2 + rho^2 cs2/(2(1-rho)) + E[W(D(1/p),RS(D(rho),p))],

    the last term by Spitzer's identity ``sum_n E[Sn^+]/n`` with
    ``Sn = rho(NB(n,1-p)+n) - n/p``, evaluated from the negative binomial pmf.
    The pmf is stepped in n by the ratio ``P(n+1)/P(n) = ((n+k)/n)p`` rather than
    formed from factorials, which would overflow long before n reaches thousands.
    """
    p = 1.0 / (1.0 + ca2)
    total = rho * ca2 + rho * rho * cs2 / (2.0 * (1.0 - rho))
    n = np.arange(1, N + 1, dtype=float)
    # log P(NB(n,1-p)=k) = lgamma(n+k) - lgamma(k+1) - lgamma(n) + n log p + k log(1-p).
    # DIVERGENCE from Algorithm 1, which steps the pmf by the ratio
    # P(n+1)/P(n) = ((n+k)/n)p: the same numbers, but a K*N scalar recursion is
    # prohibitively slow in Python, while that recursion cannot be vectorized as
    # a cumulative product without overflowing before its (1-p)^k factor tames
    # it. The JAR and the C++ port keep the recursion, where loops are cheap.
    lognfac = gammaln(n)
    for k in range(1, K + 1):
        logp = (gammaln(n + k) - gammaln(k + 1.0) - lognfac + n * np.log(p)
                + k * np.log1p(-p))
        step = (n + k) * rho - n / p
        np.maximum(step, 0.0, out=step)
        nz = step > 0
        if not np.any(nz):
            continue
        total += float(np.sum(np.exp(logp[nz]) * step[nz] / n[nz]))
    return total


def qsys_gig1_bnds_extremal(lambda_val: float, mu: float, ca: float, cs: float,
                            K: int = 4000, N: int = 2000,
                            skipTight: bool = False) -> Dict[str, Any]:
    """
    Extremal two-moment bounds for the GI/GI/1 queue.

    Two moments do not determine ``E[W]``; they determine a SET of possible
    values, and the width of that set is the honest uncertainty in any
    two-moment approximation. The extremal laws attain its ends: the lower end
    with deterministic interarrival times and a three-point service law on
    multiples of that interval (closed form, eq. 2.12), the upper end
    asymptotically with TWO-POINT laws, an interarrival law with an atom at 0 and
    a service law whose upper atom runs to infinity as its probability vanishes.
    Making an interarrival time larger only empties the queue once; making a
    service time larger delays everyone behind it, which is why the two ends look
    so different.

    The upper end is reduced to a ``D(1/p)/RS(D(rho),p)/1`` model with
    ``p = 1/(1+ca^2)`` and evaluated by Spitzer's identity with the negative
    binomial pmf, so it is a truncated numerical limit rather than a formula; the
    closed-form companion (eq. 3.4) is within about 1% of it.

    Args:
        lambda_val: arrival rate
        mu: service rate
        ca: coefficient of variation of the interarrival time
        cs: coefficient of variation of the service time
        K: truncation of the negative binomial value
        N: truncation of the random-walk length
        skipTight: skip the O(K*N) tight bound and return the closed forms only

    Returns:
        Dict of TIMES IN QUEUE (add ``1/mu`` for response times):
        ``trafficIntensity``, ``lowerBound``, ``upperBound``,
        ``upperBoundClosed``, ``upperBoundDaley``, ``upperBoundKingman``,
        ``heavyTraffic``, ``delta``, ``relativeWidth``, ``tightComputed``.

    References:
        Y. Chen, W. Whitt (2020). Algorithms for the upper bound mean waiting
        time in the GI/GI/1 queue. Queueing Systems 94, 327-356.
    """
    if lambda_val <= 0 or mu <= 0:
        raise ValueError('The arrival and service rates must be positive.')
    rho = lambda_val / mu
    if rho >= 1:
        raise ValueError('The bounds require a stable queue, rho < 1.')
    ca2, cs2 = ca ** 2, cs ** 2

    # The reference sets E[U] = 1; every waiting time therefore carries the
    # factor 1/lambda, that time unit expressed in the caller's units.
    scale = 1.0 / lambda_val
    delta = _delta(rho)                                                     # (3.5)
    result: Dict[str, Any] = {
        'trafficIntensity': rho,
        'lowerBound': scale * rho * max((1 + cs2) * rho - 1, 0.0) / (2 * (1 - rho)),      # (2.12)
        'upperBoundKingman': scale * rho ** 2 * (ca2 / rho ** 2 + cs2) / (2 * (1 - rho)),  # (2.6)
        'upperBoundDaley': scale * rho ** 2 * ((2 - rho) * ca2 / rho + cs2) / (2 * (1 - rho)),
        'heavyTraffic': scale * rho ** 2 * (ca2 + cs2) / (2 * (1 - rho)),                  # (2.9)
        'delta': delta,
        'upperBoundClosed': scale * (2 * (1 - rho) * rho / (1 - delta) * ca2
                                     + rho ** 2 * cs2) / (2 * (1 - rho)),                  # (3.4)
    }
    if skipTight:
        result['upperBound'] = result['upperBoundClosed']
        result['tightComputed'] = False
    else:
        result['upperBound'] = scale * _tight(rho, ca2, cs2, K, N)
        result['tightComputed'] = True
    ub = result['upperBound']
    result['relativeWidth'] = (ub - result['lowerBound']) / ub if ub > 0 else 0.0
    return result
