"""
CTMC transient analysis via Fox-Glynn uniformization.

Poisson weights and truncation points are obtained by the Fox-Glynn algorithm,
with Jansen's correction to the right tail estimate. Neither exp(-q*t) nor
(q*t)^k/k! is ever formed, so the method is free of overflow and underflow and
needs no horizon splitting.

Key algorithms:
    ctmc_foxglynn: Transient probabilities by Fox-Glynn uniformization
    ctmc_foxglynn_weights: Truncation window and normalized Poisson weights
"""

import math
import numpy as np
from typing import Tuple

# Smallest rate for which the Fox-Glynn tail estimates, being asymptotic normal
# approximations, are applicable.
FOX_GLYNN_MIN_LAMBDA = 25.0


def _chernoff_exponent(lam: float, k: float) -> float:
    """
    Chernoff exponent lam*h(k/lam) with h(u) = u*log(u) - u + 1, so that
    exp(-e) dominates P{X >= k} for k > lam and P{X <= k} for k < lam under
    X ~ Poisson(lam). It certifies the Fox-Glynn estimates below, which are
    asymptotic and valid only above FOX_GLYNN_MIN_LAMBDA.
    """
    if k <= 0.0:
        return lam
    return lam - k + k * math.log(k / lam)


def _right_truncation_point(lam: float, tol: float) -> int:
    """
    Right truncation point R with P{X > R} <= tol/2 for X ~ Poisson(lam).

    The starting guess is the Fox-Glynn (1988) right tail estimate: with
    m = floor(lam) and shift s(k) = k*sqrt(2*lam) + 3/2, the tail
    P{X >= m + ceil(s(k))} is bounded by

        a(lam) * d(lam,k) * exp(-k^2/2) / (k*sqrt(2*pi)),
        a(lam) = (1 + 1/lam) * exp(1/16) * sqrt(2).

    Jansen's correction (2011), applied here as the factor

        d(lam,k) = 1 / (1 - exp(-(2/9)*s(k))),

    repairs the original statement, which drops this factor for the
    k-dependent shift and is therefore optimistic at moderate lam; d tends to
    one as lam grows, so the corrected bound agrees with Fox and Glynn's
    asymptotically. The guess is then tightened and, if necessary, grown until
    the Chernoff bound above holds, so the returned R is certified in any
    regime.
    """
    target = math.log(2.0 / tol)
    mode = int(math.floor(lam))
    r = mode
    if lam >= FOX_GLYNN_MIN_LAMBDA:
        a = (1.0 + 1.0 / lam) * math.exp(1.0 / 16.0) * math.sqrt(2.0)
        spread = math.sqrt(2.0 * lam)
        for k in range(1, 65):
            shift = k * spread + 1.5
            d = 1.0 / (1.0 - math.exp(-(2.0 / 9.0) * shift))
            bound = a * d * math.exp(-0.5 * k * k) / (k * math.sqrt(2.0 * math.pi))
            if bound <= 0.5 * tol:
                r = mode + int(math.ceil(shift))
                break
    while r > mode and _chernoff_exponent(lam, r) >= target:
        r -= 1
    while _chernoff_exponent(lam, r + 1.0) < target:
        r += 1
    return r


def _left_truncation_point(lam: float, tol: float) -> int:
    """
    Left truncation point L with P{X < L} <= tol/2 for X ~ Poisson(lam), zero
    when no truncation is admissible.

    The starting guess is the Fox-Glynn (1988) left tail estimate with
    b(lam) = (1 + 1/lam)*exp(1/(8*lam)) and shift k*sqrt(lam) + 3/2 below the
    mode. Unlike the right tail this one needs no Jansen factor, the left tail
    of a Poisson being lighter than its normal approximation. It is certified
    against the Chernoff bound in the same way.
    """
    target = math.log(2.0 / tol)
    mode = int(math.floor(lam))
    if _chernoff_exponent(lam, 0.0) < target:
        return 0
    l = 0
    if lam >= FOX_GLYNN_MIN_LAMBDA:
        b = (1.0 + 1.0 / lam) * math.exp(1.0 / (8.0 * lam))
        spread = math.sqrt(lam)
        for k in range(1, 65):
            bound = b * math.exp(-0.5 * k * k) / (k * math.sqrt(2.0 * math.pi))
            if bound <= 0.5 * tol:
                l = mode - int(math.floor(k * spread + 1.5))
                break
        l = max(l, 0)
    while l > 0 and _chernoff_exponent(lam, l - 1.0) < target:
        l -= 1
    while l < mode and _chernoff_exponent(lam, float(l)) >= target:
        l += 1
    return l


def ctmc_foxglynn_weights(lam: float, tol: float = 1e-12,
                          maxiter: int = -1) -> Tuple[int, int, np.ndarray]:
    """
    Fox-Glynn truncation window and normalized Poisson weights.

    Following Fox-Glynn, the weights are built by the two-sided recursion
    w[k-1] = w[k]*k/lam and w[k+1] = w[k]*lam/(k+1) anchored at the mode, so
    neither exp(-lam) nor lam^k/k! is ever evaluated and the overflow and
    underflow that limit the direct series cannot occur. Anchoring at
    w[mode] = 1 keeps the extreme weights near tol, far above the denormal
    threshold, making Fox and Glynn's rescaling of the mode weight unnecessary
    here. The normalizing sum is accumulated in increasing order of magnitude.

    Args:
        lam: Poisson rate, that is the uniformization constant times the horizon
        tol: Poisson tail-mass truncation tolerance
        maxiter: Cap on the right truncation point; nonpositive leaves it uncapped

    Returns:
        Tuple of (left truncation point, right truncation point, weights)
    """
    if tol <= 0.0:
        tol = 1e-12
    left = _left_truncation_point(lam, tol)
    right = _right_truncation_point(lam, tol)
    if maxiter > 0 and right > maxiter:
        right = maxiter
        left = min(left, right)

    n = right - left + 1
    w = np.zeros(n, dtype=np.float64)
    mode = min(max(int(math.floor(lam)), left), right)
    w[mode - left] = 1.0
    for k in range(mode, left, -1):
        w[k - 1 - left] = w[k - left] * k / lam
    for k in range(mode, right):
        w[k + 1 - left] = w[k - left] * lam / (k + 1)
    return left, right, w / np.sum(np.sort(w))


def ctmc_foxglynn(pi0: np.ndarray, Q: np.ndarray, t: float,
                  tol: float = 1e-12, maxiter: int = -1) -> np.ndarray:
    """
    Transient distribution of a CTMC by Fox-Glynn uniformization.

    Args:
        pi0: Initial probability distribution
        Q: Infinitesimal generator matrix
        t: Transient analysis period boundary [0,t]
        tol: Poisson tail-mass truncation tolerance
        maxiter: Maximum truncation depth; pass a nonpositive value to let the
            Fox-Glynn right truncation point size it

    Returns:
        Transient probability vector at time t
    """
    Q = np.asarray(Q, dtype=np.float64)
    pi0 = np.asarray(pi0, dtype=np.float64).flatten()

    q = 1.1 * float(np.max(np.abs(np.diag(Q))))
    lam = q * t
    if q <= 0.0 or lam <= 0.0:
        return pi0.copy()

    left, right, w = ctmc_foxglynn_weights(lam, tol, maxiter)

    Qs = np.eye(Q.shape[0]) + Q / q
    pi = np.zeros(Q.shape[0], dtype=np.float64)
    p = pi0.copy()
    for k in range(right + 1):
        if k >= left:
            pi += w[k - left] * p
        if k < right:
            p = p @ Qs
    return pi


__all__ = [
    'ctmc_foxglynn',
    'ctmc_foxglynn_weights',
]
