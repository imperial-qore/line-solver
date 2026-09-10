"""
Two-moment approximations for the maximum of iid random variables.

Native Python twin of matlab/src/api/qsys/qsys_maxima_twomoment.m, implementing
C. Crow, D. Goldberg and W. Whitt (2007), Two-moment approximations for maxima,
Operations Research 55(3), 532-548.
"""

from math import exp, log, sqrt
from typing import Any, Dict, Optional

import numpy as np

EULER = 0.5772156649015329


def qsys_maxima_twomoment(n: int, mean: float, cs2: float, q: Optional[float] = None,
                          exactFitted: bool = True) -> Dict[str, Any]:
    """
    Approximate the maximum of ``n`` iid non-negative variables from two moments.

    THE SHAPE OF THE ANSWER. For a law with an exponential-like tail the maximum
    of ``n`` samples grows like ``c~^2 (log n + ...)``: doubling ``n`` adds a
    constant, it does not scale the answer. What the two moments buy is the
    SLOPE ``c~^2`` of that logarithm and an offset ``eta``:

        x_n(q) = c~^2 [ log(n eta) - log log(1/q) ],  E[M_n] = c~^2 [log(n eta) + gamma]

    with, for ``cs2 >= 1``, ``c~^2 = cs2`` and ``eta = (cs2+1)/(2 cs2^2)`` from
    the H2 representative, and for ``cs2 < 1`` the shifted-exponential
    representative ``c~^2 = sqrt(cs2)``, ``eta = exp((1-sqrt(cs2))/sqrt(cs2))``.

    WHEN NOT TO USE IT. The extreme-value form needs ``n`` past a threshold
    ``n* ~ cs2/q``, because with a highly variable law most of the ``n`` samples
    come from the short component and only about ``n p`` of them can contend for
    the maximum. Measured against exact maxima, the closed form is within a few
    percent for ``n >= 100`` at ``cs2 = 4`` and ``16``, and useless at ``n = 10``
    for ``cs2 = 16`` -- which is exactly what ``n*`` predicts.

    AND WHEN TWO MOMENTS ARE NOT ENOUGH. Below ``cs2 = 1`` the maximum is
    genuinely family-dependent: an Erlang and a shifted exponential with the same
    two moments have maxima that differ by tens of percent and diverge as ``n``
    grows, because their tails decay at different rates. The paper's own caution.
    ``exactFitted`` therefore also returns the maximum computed exactly from the
    fitted representative, which is the reliable route it recommends.

    Args:
        n: the number of samples
        mean: the mean of the underlying law
        cs2: its squared coefficient of variation
        q: a quantile level in (0,1); the mean is returned when absent
        exactFitted: also compute the maximum exactly from the fitted
            representative distribution, by integrating ``1-F^n``

    Returns:
        Dict with ``value`` (the closed-form mean or quantile), ``slope``
        (``c~^2 mean``), ``eta``, ``threshold`` (n*), ``reliable`` (whether
        ``n >= n*``), ``family``, and, when requested, ``exactFittedValue``.

    References:
        C. Crow, D. Goldberg, W. Whitt (2007). Two-moment approximations for
        maxima. Operations Research 55(3), 532-548.
    """
    if n < 1:
        raise ValueError('At least one sample is required.')
    if mean <= 0:
        raise ValueError('The mean must be positive.')
    if cs2 <= 0:
        raise ValueError('The squared coefficient of variation must be positive.')
    if q is not None and not (0 < q < 1):
        raise ValueError('The quantile level must lie in (0,1).')

    if cs2 >= 1:
        ct = cs2
        eta = (cs2 + 1.0) / (2.0 * cs2 * cs2)
        family = 'H2'
    else:
        ct = sqrt(cs2)
        eta = exp((1.0 - sqrt(cs2)) / sqrt(cs2))
        family = 'shifted exponential'
    inner = log(n * eta) + (EULER if q is None else -log(log(1.0 / q)))
    value = mean * ct * inner
    qq = 0.5 if q is None else q
    threshold = cs2 / qq                      # eq. (4.22)
    out: Dict[str, Any] = {
        'value': value,
        'slope': mean * ct,
        'eta': eta,
        'threshold': threshold,
        'reliable': n >= threshold,
        'family': family,
    }

    if exactFitted:
        # The paper's other recommendation: fit the representative law, then
        # compute the maximum exactly from F^n rather than from its tail.
        if cs2 >= 1:
            # H2 with balanced means, matched on two moments.
            p1 = 0.5 * (1.0 + sqrt((cs2 - 1.0) / (cs2 + 1.0)))
            l1 = 2.0 * p1 / mean
            l2 = 2.0 * (1.0 - p1) / mean
            ccdf = lambda t: p1 * np.exp(-l1 * t) + (1.0 - p1) * np.exp(-l2 * t)
            hi = 40.0 * mean * max(cs2, 1.0)
        else:
            d = mean * (1.0 - sqrt(cs2))
            m = mean * sqrt(cs2)
            ccdf = lambda t: np.where(t <= d, 1.0, np.exp(-(t - d) / m))
            hi = d + 40.0 * m
        grid = np.linspace(0.0, hi, 200001)
        f = 1.0 - (1.0 - ccdf(grid)) ** n
        if q is None:
            out['exactFittedValue'] = float(np.trapezoid(f, grid))
        else:
            cdfn = (1.0 - ccdf(grid)) ** n
            idx = int(np.searchsorted(cdfn, q))
            out['exactFittedValue'] = float(grid[min(idx, grid.size - 1)])
    return out
