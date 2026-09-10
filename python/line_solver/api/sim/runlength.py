"""
Run-length planning for steady-state simulation.

Native Python twin of matlab/src/api/sim/sim_runlength.m, sim_asymvar_mm1.m and
sim_asymvar_ctmc.m, implementing the planning method of W. Whitt (1989),
Planning queueing simulations, Management Science 35(11), 1341-1366.
"""

from math import erfc, sqrt
from typing import Any, Dict, Optional, Sequence

import numpy as np


def _z(confidence: float) -> float:
    """Two-sided normal quantile for the given confidence, by bisection on erfc."""
    if not (0 < confidence < 1):
        raise ValueError('The confidence must lie in (0,1).')
    target = 1.0 - confidence
    lo, hi = 0.0, 40.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if erfc(mid / sqrt(2.0)) > target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def sim_asymvar_mm1(lambda_val: float, mu: float) -> Dict[str, float]:
    """
    Asymptotic variance of the M/M/1 number-in-system process.

    THE QUANTITY THAT MATTERS FOR PLANNING is not the variance of the process but
    its ASYMPTOTIC VARIANCE, ``sigma^2 = lim t Var(time-average over [0,t])``,
    which is twice the integral of the autocovariance. It is what says how long a
    run must be, because the time average of a positively correlated process
    converges at rate ``sigma^2/t``, not at ``Var(X)/t``.

    For M/M/1 with utilization ``rho``,

        E[N] = rho/(1-rho),  Var(N) = rho/(1-rho)^2,
        sigma^2 = 2 rho (1+rho) / (mu (1-rho)^4).

    The FOURTH power is the whole story: the variance of the process itself grows
    like ``(1-rho)^-2``, but the run length needed to average it away grows like
    ``(1-rho)^-4`` divided by the squared mean, i.e. like ``(1-rho)^-2``. A queue
    at rho = 0.9 needs about 100 times the run of one at rho = 0.

    Args:
        lambda_val: arrival rate
        mu: service rate

    Returns:
        Dict with ``mean``, ``variance``, ``asymptoticVariance`` and
        ``relaxationTime`` (``sigma^2/Var``, the correlation time scale).

    References:
        W. Whitt (1989). Planning queueing simulations. Management Science
        35(11), 1341-1366; J. Abate, W. Whitt (1988). The correlation functions
        of RBM and M/M/1. Stochastic Models 4(2), 315-359.
    """
    if lambda_val <= 0 or mu <= 0:
        raise ValueError('The arrival and service rates must be positive.')
    rho = lambda_val / mu
    if rho >= 1:
        raise ValueError('The queue must be stable, rho < 1.')
    mean = rho / (1.0 - rho)
    var = rho / (1.0 - rho) ** 2
    asym = 2.0 * rho * (1.0 + rho) / (mu * (1.0 - rho) ** 4)
    return {'mean': mean, 'variance': var, 'asymptoticVariance': asym,
            'relaxationTime': asym / var if var > 0 else 0.0}


def sim_asymvar_ctmc(A, f: Sequence[float], pi: Optional[Sequence[float]] = None) -> Dict[str, Any]:
    """
    Asymptotic variance of a reward on a continuous-time Markov chain.

    ``sigma^2 = 2 sum_x pi(x) g(x) d(x)`` where ``g = f - E_pi[f]`` and ``d``
    solves ``A d = -g`` with ``pi d = 0``: the DEVIATION vector, the accumulated
    future excess reward started from each state. This is the general form of the
    quantity :func:`sim_asymvar_mm1` gives in closed form for M/M/1, and it is
    what run-length planning needs for any model LINE can build a generator for.

    Args:
        A: the generator matrix, rows summing to zero
        f: the reward attached to each state
        pi: the stationary distribution; solved for when absent

    Returns:
        Dict with ``mean``, ``variance``, ``asymptoticVariance`` and
        ``deviation`` (the vector d).

    References:
        W. Whitt (1989). Planning queueing simulations. Management Science
        35(11), 1341-1366.
    """
    A = np.asarray(A, dtype=float)
    n = A.shape[0]
    if A.shape != (n, n):
        raise ValueError('The generator must be square.')
    if np.max(np.abs(A.sum(axis=1))) > 1e-8:
        raise ValueError('The generator rows must sum to zero.')
    f = np.asarray(f, dtype=float).reshape(n)
    if pi is None:
        M = np.vstack([A.T[:-1], np.ones(n)])
        b = np.zeros(n)
        b[-1] = 1.0
        pi = np.linalg.lstsq(M, b, rcond=None)[0]
    pi = np.asarray(pi, dtype=float).reshape(n)
    mean = float(pi @ f)
    g = f - mean
    # A d = -g pins d only up to a constant, so one equation of A is redundant
    # and one normalization replaces it. WHICH equation is dropped matters: the
    # rows of A are related by pi A = 0, so a row whose pi is tiny is only
    # nominally redundant and dropping it loses real information -- on a queue
    # truncated where pi has underflowed, that alone puts sigma^2 out by orders
    # of magnitude. Dropping the row with the LARGEST pi is the well-conditioned
    # choice.
    drop = int(np.argmax(pi))
    keep = [i for i in range(n) if i != drop]
    M2 = np.vstack([A[keep], pi])
    b2 = np.concatenate([-g[keep], [0.0]])
    d = np.linalg.solve(M2, b2)
    var = float(pi @ (g * g))
    return {'mean': mean, 'variance': var,
            'asymptoticVariance': 2.0 * float(pi @ (g * d)), 'deviation': d}


def sim_runlength(mean: float, asymptoticVariance: float, relPrecision: float = 0.05,
                  confidence: float = 0.95, runLength: Optional[float] = None) -> Dict[str, float]:
    """
    Run length needed for a steady-state estimate of a given relative precision.

    A time average over ``[0,t]`` has standard error ``sqrt(sigma^2/t)``, so a
    two-sided interval of half-width ``z sqrt(sigma^2/t)`` reaches relative
    precision ``eps`` when

        t* = (z/eps)^2 sigma^2 / mean^2.

    THE POINT OF THE FORMULA is that everything expensive is in
    ``sigma^2/mean^2``, the squared coefficient of variation of the TIME AVERAGE
    rather than of the process. Halving the tolerance quadruples the run.

    Args:
        mean: the steady-state mean being estimated
        asymptoticVariance: sigma^2 of that estimator
        relPrecision: the target half-width as a fraction of the mean
        confidence: the confidence level of the interval
        runLength: an actual run length, to report the precision it buys

    Returns:
        Dict with ``requiredRunLength``, ``z``, and, when ``runLength`` is given,
        ``halfWidth`` and ``achievedRelPrecision``.

    References:
        W. Whitt (1989). Planning queueing simulations. Management Science
        35(11), 1341-1366.
    """
    if mean == 0:
        raise ValueError('A relative precision is meaningless for a zero mean.')
    if asymptoticVariance < 0:
        raise ValueError('The asymptotic variance cannot be negative.')
    if relPrecision <= 0:
        raise ValueError('The relative precision must be positive.')
    z = _z(confidence)
    t = (z / relPrecision) ** 2 * asymptoticVariance / (mean * mean)
    out = {'requiredRunLength': t, 'z': z}
    if runLength is not None and runLength > 0:
        hw = z * sqrt(asymptoticVariance / runLength)
        out['halfWidth'] = hw
        out['achievedRelPrecision'] = hw / abs(mean)
    return out


def sim_runlength_plan(means, ciHalfWidth, samplesUsed: int,
                       relPrecision: float = 0.05,
                       confidence: float = 0.95) -> Dict[str, Any]:
    """
    How long a simulation run should have been, from the one it already did.

    A batch-means half-width H at confidence 1-alpha over a run of N samples
    pins the ASYMPTOTIC variance of the estimator,

        sigma^2 = (H/z)^2 N,   z = Phi^-1((1+confidence)/2),

    and that is the quantity a run length is planned from -- NOT the stationary
    variance, which on M/M/1 differs from it by a factor blowing up like
    ``(1-rho)^-2``. :func:`sim_runlength` then turns it into the sample count
    that reaches a requested RELATIVE precision.

    Args:
        means: one mean per (station, class)
        ciHalfWidth: the confidence-interval half-width of the same entries
        samplesUsed: the run length those half-widths came from
        relPrecision: the relative precision to plan for
        confidence: the level the half-widths were computed at

    Returns:
        Dict with ``relprecision``, ``confidence``, ``samplesUsed``,
        ``asymptoticVariance`` and ``requiredSamples``; an entry with a
        non-positive mean or half-width is left NaN, since there is nothing to
        plan from there.

    References:
        W. Whitt (1989). Planning queueing simulations. Management Science
        35(11), 1341-1366.

    See also:
        matlab/src/api/sim/sim_runlength_plan.m
    """
    means = np.asarray(means, dtype=float)
    ciHalfWidth = np.asarray(ciHalfWidth, dtype=float)
    if int(samplesUsed) <= 0:
        raise ValueError('sim_runlength_plan: the number of samples already used must be positive')
    from math import erf, sqrt
    # z = Phi^-1((1+confidence)/2), by bisection on erf so no statistics package
    # is needed here either.
    lo, hi = 0.0, 40.0
    target = float(confidence)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if erf(mid / sqrt(2.0)) < target:
            lo = mid
        else:
            hi = mid
    z = 0.5 * (lo + hi)
    asym = np.full(means.shape, np.nan)
    req = np.full(means.shape, np.nan)
    it = np.ndindex(means.shape)
    for idx in it:
        h = float(ciHalfWidth[idx]) if ciHalfWidth.shape == means.shape else float('nan')
        m = float(means[idx])
        if not np.isfinite(h) or h <= 0 or not np.isfinite(m) or m <= 0:
            continue
        av = (h / z) ** 2 * float(samplesUsed)
        asym[idx] = av
        req[idx] = sim_runlength(m, av, relPrecision=relPrecision,
                                 confidence=confidence)['requiredRunLength']
    return {'relprecision': relPrecision, 'confidence': confidence,
            'samplesUsed': int(samplesUsed), 'asymptoticVariance': asym,
            'requiredSamples': req}
