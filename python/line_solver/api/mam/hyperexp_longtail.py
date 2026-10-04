"""
Fitting a hyperexponential to a long-tail distribution.

Native Python twin of matlab/src/api/mam/hyperexp_fit_longtail.m, implementing
the recursive procedure of A. Feldmann and W. Whitt (1998), Fitting mixtures of
exponentials to long-tail distributions to analyze network performance models,
Performance Evaluation 31, 245-279, Section 4.
"""

from math import exp, log
from typing import Any, Callable, Dict, Optional, Sequence

import numpy as np


def _grid_error(fitted, ccdf, lo: float, hi: float) -> float:
    """
    Largest relative error on a log grid across the covered range. The fit is
    exact at the fitting arguments by construction, so this is what says whether
    it also holds BETWEEN them.
    """
    ts = np.exp(np.linspace(np.log(lo), np.log(hi), 200))
    worst = 0.0
    for t in ts:
        target = float(ccdf(float(t)))
        if target > 1e-300:
            worst = max(worst, abs(fitted(float(t)) - target) / target)
    return worst


def _fit_fixed_k(ccdf, k, c1, b, decade, points):
    """The recursion at a fixed component count; the public entry point retries."""
    return hyperexp_fit_longtail(ccdf, k, c1, b, decade, points)


def hyperexp_fit_longtail(ccdf: Callable[[float], float], k: Optional[int] = None,
                          c1: float = None,
                          b: float = 1.5, decade: float = 4.0,
                          points: Optional[Sequence[float]] = None) -> Dict[str, Any]:
    """
    Fit a hyperexponential to a long-tail distribution, recursively over time
    scales.

    WHY MOMENTS ARE THE WRONG HANDLE. A Pareto law with tail index below 2 has
    infinite variance, so no two- or three-moment fit exists at all; and even
    when the moments are finite, matching them says nothing about the several
    ORDERS OF MAGNITUDE of time scale over which a long-tail distribution
    actually acts. This procedure matches the CCDF ITSELF at points spread across
    those decades.

    THE RECURSION. Order the components so that ``lambda_1 < ... < lambda_k``.
    In the far tail only the slowest component survives, so ``(p_1, lambda_1)``
    can be fitted there alone, from the ccdf at ``c_1`` and ``b c_1``:

        lambda_1 = ln(F^c(c_1)/F^c(b c_1)) / ((b-1)c_1),
        p_1 = F^c(c_1) exp(lambda_1 c_1).

    Subtract that component from the ccdf and repeat one decade lower, and so on
    (eqs. 4.6-4.11). The last component takes whatever probability is left,
    ``p_k = 1 - sum_{j<k} p_j``, and its rate follows from the ccdf at ``c_k``
    (eqs. 4.12-4.14). This is Prony's method applied to a ccdf.

    Args:
        ccdf: F^c(t) = P(X > t) of the distribution to approximate
        k: number of exponential components; ``None`` takes one per decade
            between the 0.9 quantile and the 1e-6 quantile, which is the range
            the spacing ``decade`` can actually cover
        c1: the largest fitting argument; defaults to the point where the ccdf
            falls below 1e-6, which puts the slowest component in the real tail
        b: the within-scale spacing, 1 < b < c_i/c_{i+1}. The default pair
            (b, decade) = (1.5, 4) is not the paper's illustrative (2, 10): the
            algorithm is exact AT the fitting arguments and free between them,
            and measured on a Weibull(0.3) the tighter grid cuts the worst
            between-point error from about 54% to 12%, at the cost of more
            components. Pass (2, 10) for the paper's own figures
        decade: the ratio between successive fitting arguments,
            ``c_i = c_1 decade^-(i-1)``; recomputed automatically when c1 is None
            so that the k arguments span from the 0.9 quantile to the 1e-6 one
        points: explicit decreasing fitting arguments, overriding c1 and decade

    Returns:
        Dict with ``p`` (the mixing probabilities), ``lambda`` (the rates),
        ``points`` (the c_i used), ``mean`` of the fitted law, ``targetMean``
        of the original one (integrated over the covered range), ``coverage``
        (the interval the fit is constrained on) and
        ``maxRelError`` of the fitted ccdf at the fitting arguments. The last
        component matches only at ``c_k``, its weight being fixed by the total
        probability, so the error at ``b c_k`` is not zero by construction.

    References:
        A. Feldmann, W. Whitt (1998). Fitting mixtures of exponentials to
        long-tail distributions to analyze network performance models.
        Performance Evaluation 31, 245-279.
    """
    if not (b > 1):
        raise ValueError('The spacing b must exceed 1.')
    if decade <= b:
        raise ValueError('The decade ratio must exceed the spacing b, or the fitting arguments '
                         'would interleave.')

    def _quantile(prob: float) -> float:
        """Smallest t with F^c(t) <= prob, by doubling then bisection."""
        hi = 1.0
        while ccdf(hi) > prob:
            hi *= 2.0
            if hi > 1e15:
                raise ValueError('the ccdf does not decay, so there is no tail to fit')
        lo = 0.0
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if ccdf(mid) > prob:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)

    if k is None and points is None:
        # One component per decade between the body and the tail: that is what
        # the spacing decade = 10 buys, and asking for more components than
        # decades is exactly what breaks the recursion below.
        top = _quantile(1e-6) if c1 is None else c1
        body = _quantile(0.9)
        if body <= 0 or top <= body:
            raise ValueError('the ccdf gives no usable range of time scales to fit over')
        k = max(2, int(round(np.log(top / body) / np.log(decade))) + 1)
        # The recursion needs each component to dominate at its own scale. Near
        # the body of a law with a lot of mass there (a Pareto, say) that fails,
        # and the remaining probability runs out; back off one component at a
        # time until it holds. Only the AUTOMATIC count retries: an explicit k
        # that cannot be fitted is an error the caller asked for.
        for kk in range(k, 1, -1):
            try:
                return _fit_fixed_k(ccdf, kk, c1, b, decade, None)
            except ValueError:
                continue
        raise ValueError('no component count from %d down to 2 admits the recursion; the ccdf may '
                         'not be long-tailed enough for this scheme' % k)
    k = int(k)
    if k < 1:
        raise ValueError('At least one exponential component is required.')

    if points is not None:
        cs = np.asarray(points, dtype=float)
        if cs.size != k:
            raise ValueError('One fitting argument per component is required.')
        if np.any(np.diff(cs) >= 0):
            raise ValueError('The fitting arguments must be strictly decreasing.')
    else:
        if c1 is None:
            # Put the slowest component where the tail actually is.
            c1 = _quantile(1e-6)
        cs = np.array([c1 * decade ** (-i) for i in range(k)], dtype=float)

    p = np.zeros(k)
    lam = np.zeros(k)
    for i in range(k):
        ci = cs[i]
        # Eqs. (4.6)-(4.7): what the already-fitted, slower components leave.
        resid_c = float(ccdf(ci)) - float(np.sum(p[:i] * np.exp(-lam[:i] * ci)))
        resid_bc = float(ccdf(b * ci)) - float(np.sum(p[:i] * np.exp(-lam[:i] * b * ci)))
        if i < k - 1:
            if resid_c <= 0 or resid_bc <= 0 or resid_c <= resid_bc:
                raise ValueError(
                    'the residual ccdf is not positive and decreasing at fitting argument %g. The '
                    'recursion needs the fitting arguments well separated, c_i/c_(i+1) >> b, so '
                    'that only the slowest surviving component matters at each scale; widen '
                    'decade, lower k, or move c1 further into the tail' % ci)
            lam[i] = log(resid_c / resid_bc) / ((b - 1.0) * ci)      # eq. (4.10)
            p[i] = resid_c * exp(lam[i] * ci)                        # eq. (4.11)
        else:
            # Eqs. (4.12)-(4.14): the last component takes the rest of the mass.
            p[i] = 1.0 - float(np.sum(p[:i]))
            if p[i] <= 0:
                raise ValueError('the fitted components already carry all the probability, so the '
                                 'last one has none left; lower k or move c1 further into the tail')
            if resid_c <= 0:
                raise ValueError('the residual ccdf has gone non-positive at the last fitting '
                                 'argument; lower k or move c1 further into the tail')
            lam[i] = log(p[i] / resid_c) / ci                        # eq. (4.14)
        if lam[i] <= 0:
            raise ValueError('a non-positive rate came out of the fit at argument %g; the ccdf is '
                             'not decaying fast enough there for this many components' % ci)

    fitted = lambda t: float(np.sum(p * np.exp(-lam * t)))
    # The target mean, for the caller to compare against: the fit is only
    # constrained on [c_k, b c_1], and a mean lives wherever the body is, so a
    # k too small to reach the body shows up here and nowhere else.
    hi = cs[0] * b
    grid = np.linspace(0.0, hi, 20001)
    target_mean = float(np.trapezoid(np.array([float(ccdf(float(t))) for t in grid]), grid))
    errs = []
    for ci in cs:
        for t in (ci, b * ci):
            target = float(ccdf(t))
            if target > 0:
                errs.append(abs(fitted(t) - target) / target)
    return {
        'p': p,
        'lambda': lam,
        'points': cs,
        'mean': float(np.sum(p / lam)),
        'targetMean': target_mean,
        'coverage': (float(cs[-1]), float(cs[0] * b)),
        'maxRelError': max(errs) if errs else 0.0,
        'maxRelErrorGrid': _grid_error(fitted, ccdf, cs[-1], hi),
    }
