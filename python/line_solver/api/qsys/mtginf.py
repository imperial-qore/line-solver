"""
The Mt/G/infinity queue: exact time-varying analysis.

Native Python twin of matlab/src/api/qsys/qsys_mtginf.m, implementing the
"physics" of S. G. Eick, W. A. Massey and W. Whitt (1993), Operations Research
41(4), 731-742.
"""

from typing import Any, Callable, Dict, Optional, Sequence

import numpy as np


def _as_array(f: Callable, x: np.ndarray) -> np.ndarray:
    """
    Evaluate a user handle on a grid, accepting either an array-aware handle or
    a scalar one. Trying the array call first keeps the vectorized case fast,
    which matters because these grids have thousands of points.
    """
    try:
        y = np.asarray(f(x), dtype=float)
        if y.shape == x.shape:
            return y
        if y.size == 1:
            return np.full(x.shape, float(y))
    except Exception:
        pass
    return np.array([float(f(float(xi))) for xi in x])


def _simpson_grid(a: float, b: float, n: int):
    """Nodes and weights of the composite Simpson rule on an even panel count."""
    if n % 2 == 1:
        n += 1
    x = np.linspace(a, b, n + 1)
    w = np.ones(n + 1)
    w[1:-1:2] = 4.0
    w[2:-1:2] = 2.0
    return x, w * (b - a) / (3.0 * n)


def _tail_cut(ccdf: Callable[[float], float], tol: float, cap: float) -> float:
    """Smallest doubling point at which the service ccdf is below tol."""
    x = 1.0
    while float(ccdf(x)) > tol:
        x *= 2.0
        if x > cap:
            return cap
    return x


def _mean_curve(lambdaFun, serviceCcdf, t, startTime, cut, panels, unbounded) -> np.ndarray:
    """
    The Poisson mean m(t) alone, shared by the public entry point and by the
    finite-difference departure rate so that neither re-derives the other.
    """
    t = np.atleast_1d(np.asarray(t, dtype=float))
    if unbounded:
        xs, ws = _simpson_grid(0.0, cut, panels)
        gcs = _as_array(serviceCcdf, xs)
    out = np.zeros(t.size)
    for i, ti in enumerate(t):
        if unbounded:
            x, w, gc = xs, ws, gcs
        else:
            hi = min(cut, max(0.0, float(ti) - startTime))
            x, w = _simpson_grid(0.0, hi, panels)
            gc = _as_array(serviceCcdf, x)
        out[i] = float(np.sum(w * _as_array(lambdaFun, ti - x) * gc))
    return out


def qsys_mtginf(lambdaFun: Callable[[Any], Any], serviceCcdf: Callable[[Any], Any],
                ES: float, tvals: Sequence[float],
                startTime: float = -np.inf, ES2: Optional[float] = None,
                servicePdf: Optional[Callable[[Any], Any]] = None,
                tol: float = 1e-12, panels: int = 4000,
                maxAge: float = 1e12) -> Dict[str, Any]:
    """
    Exact time-varying analysis of the Mt/G/infinity queue.

    With a non-homogeneous Poisson arrival rate ``lambda(t)`` and iid service
    times ``S``, the number in system at time ``t`` is POISSON with mean

    .. math:: m(t) = E\\left[\\int_{t-S}^{t}\\lambda(u)du\\right]
              = E[S]\\,E[\\lambda(t-S_e)] = \\int_0^\\infty \\lambda(t-x)P(S>x)dx

    where ``S_e`` is the stationary-excess (equilibrium) law of ``S``, with
    density ``P(S>x)/E[S]``. This is exact, not an approximation: infinitely
    many servers mean customers never interact, so the model is a Poisson random
    measure and the whole distribution is known.

    THE PHYSICS. Writing the mean as ``E[S] E[lambda(t - S_e)]`` says the
    time-varying load is the stationary load ``E[S]lambda(t)`` subjected to a
    TIME LAG and a SPACE SHIFT: to first order
    ``m(t) ~ E[S] lambda(t - E[S_e])`` with ``E[S_e] = E[S^2]/(2E[S])``, so peak
    congestion lags peak arrival rate, and by more than the mean service time
    when the service law is variable. The pointwise stationary approximation
    ``E[S]lambda(t)`` is the zeroth-order term, which is why it misses the lag.

    Args:
        lambdaFun: the arrival rate, ideally array-aware; must accept arguments
            in the past when ``startTime`` is infinite
        serviceCcdf: G^c(x) = P(S > x)
        ES: the mean service time
        tvals: the times at which to evaluate
        startTime: time the system started empty; the default -inf assumes the
            arrival rate has been running forever
        ES2: the second moment of the service time, for the lag approximation
        servicePdf: the service density, used for the exact departure rate; when
            absent the departure rate comes from the flow balance
            ``m'(t) = lambda(t) - delta(t)`` by a central difference
        tol: service-tail cut for the age integral
        panels: Simpson panels for that integral
        maxAge: cap on the age integrated over

    Returns:
        Dict with ``times``, ``meanNumber`` (the Poisson mean m(t)),
        ``varNumber`` (equal to it), ``departureRate``, ``arrivalRate``,
        ``offeredLoadPSA`` (the pointwise stationary approximation
        ``E[S]lambda(t)``) and, when ``ES2`` is given, ``meanLag`` (``E[S_e]``)
        and ``lagApproximation`` (``E[S]lambda(t-E[S_e])``).

    References:
        S. G. Eick, W. A. Massey, W. Whitt (1993). The physics of the Mt/G/inf
        queue. Operations Research 41(4), 731-742.
    """
    if ES <= 0:
        raise ValueError('The mean service time ES must be positive.')
    t = np.atleast_1d(np.asarray(tvals, dtype=float))
    cut = _tail_cut(serviceCcdf, tol, maxAge)
    unbounded = not np.isfinite(startTime)

    # With an infinite past the age grid does not move with t, so the service
    # ccdf is evaluated once for every time point rather than once per point.
    if unbounded:
        xs, ws = _simpson_grid(0.0, cut, panels)
        gcs = _as_array(serviceCcdf, xs)
        pdfs = _as_array(servicePdf, xs) if servicePdf is not None else None

    mean = np.zeros(t.size)
    dep = np.zeros(t.size) if servicePdf is not None else None
    for i, ti in enumerate(t):
        if unbounded:
            x, w, gc = xs, ws, gcs
            pdf = pdfs
        else:
            hi = min(cut, max(0.0, float(ti) - startTime))
            x, w = _simpson_grid(0.0, hi, panels)
            gc = _as_array(serviceCcdf, x)
            pdf = _as_array(servicePdf, x) if servicePdf is not None else None
        lam = _as_array(lambdaFun, ti - x)
        # m(t) = int lambda(t-x) P(S>x) dx: the arrivals of age x still in service.
        mean[i] = float(np.sum(w * lam * gc))
        if pdf is not None:
            dep[i] = float(np.sum(w * lam * pdf))

    arrival = _as_array(lambdaFun, t)
    result: Dict[str, Any] = {
        'times': t,
        'meanNumber': mean,
        'varNumber': mean.copy(),          # Poisson: the variance is the mean
        'arrivalRate': arrival,
        'offeredLoadPSA': ES * arrival,
    }

    if dep is not None:
        result['departureRate'] = dep
    else:
        # Flow balance m'(t) = lambda(t) - delta(t), differentiated centrally.
        h = 1e-5 * max(1.0, float(np.max(np.abs(t))) if t.size else 1.0)
        up = _mean_curve(lambdaFun, serviceCcdf, t + h, startTime, cut, panels, unbounded)
        dn = _mean_curve(lambdaFun, serviceCcdf, t - h, startTime, cut, panels, unbounded)
        result['departureRate'] = arrival - (up - dn) / (2.0 * h)

    if ES2 is not None:
        lag = ES2 / (2.0 * ES)                       # E[S_e], the time lag
        result['meanLag'] = lag
        result['lagApproximation'] = ES * _as_array(lambdaFun, t - lag)

    return result
