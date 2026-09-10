"""Exact interval-valued MVA for single-class closed product-form networks.

Luthi and Haring, "Mean value analysis for queueing network models with
intervals as input parameters", Performance Evaluation 32(3):185-215, 1998.
Ported at parity from MATLAB pfqn_mva_interval.m.
"""

from typing import NamedTuple

import numpy as np

from .mva import pfqn_mva

__all__ = ['pfqn_mva_interval', 'PfqnMvaIntervalResult']

_ZERO_TOL = 1e-14


class PfqnMvaIntervalResult(NamedTuple):
    """Interval-valued mean performance measures. Every field is [lower, upper]."""
    X: np.ndarray      # throughput interval (2,)
    Q: np.ndarray      # per-station queue-length intervals (M, 2)
    U: np.ndarray      # per-station utilization enclosures (M, 2)
    R: np.ndarray      # per-station residence-time intervals (M, 2)
    Rtot: np.ndarray   # total response-time interval (2,)
    Qtot: np.ndarray   # total number of jobs at the stations (2,)


def _mva(L, n, z):
    """Ordinary single-class MVA at one corner: returns (X, Q, R)."""
    XN, _, QN, _, RN, _, _ = pfqn_mva(np.asarray(L, dtype=float).reshape(-1, 1),
                                      np.array([n], dtype=float),
                                      np.array([z], dtype=float))
    return float(np.ravel(XN)[0]), np.ravel(QN).astype(float), np.ravel(RN).astype(float)


def _endpoints(v, what):
    v = np.atleast_1d(np.asarray(v, dtype=float)).ravel()
    if v.size == 1:
        return float(v[0]), float(v[0])
    if v.size != 2:
        raise ValueError('%s must be given as a scalar or as a [lower upper] interval.' % what)
    return float(v[0]), float(v[1])


def pfqn_mva_interval(L, N, Z=0.0):
    """Exact hull of single-class MVA over an input box.

    Single-class MVA is monotone in every input: the throughput decreases in
    each demand and in the think time and increases in the population, the
    per-station queue length and residence time increase in the own demand and
    in the population and decrease in the other demands and in the think time,
    and the totals increase in every demand and in the population and decrease
    in the think time (Luthi and Haring 1998, Theorems 2-5, Table 1). By their
    Theorem 1 the exact range of a function monotone in each argument is
    attained at the endpoints of the input box, so each bound below is one
    ordinary MVA call at the corner that the sign pattern selects. This is the
    algorithm of their Fig. 2 and it costs 2*(m+2) MVA calls, m being the number
    of thick demand intervals; evaluating the MVA recursion in interval
    arithmetic instead would be a valid but far wider enclosure, since every
    input recurs at each step (the dependency problem, 14x too wide on the
    paper's own example).

    The returned interval is the exact hull of MVA over the input box, not a
    bound on the true network: it holds conditionally on the demands lying in
    the box, and says nothing about the accuracy of MVA itself. It must
    therefore not be composed with the brackets of SolverBA, which bracket the
    exact solution of a model whose demands are known.

    Delay stations are folded into Z, exactly as in pfqn_mva: a delay demand
    interval enters as a term of the think-time interval, and the hull of the
    sum is the sum of the hulls when the delays vary independently.
    Load-independent single-server queueing stations only, one class only; the
    monotonicity theorems cover no other case.

    Args:
        L: Service demand intervals (M, 2), column 0 lower, column 1 upper. An
            (M,) vector is read as a thin box.
        N: Population interval [nlo, nup], or a scalar for a thin population.
        Z: Think time interval [zlo, zup], or a scalar (default 0).

    Returns:
        PfqnMvaIntervalResult with fields X, Q, U, R, Rtot, Qtot.
    """
    L = np.asarray(L, dtype=float)
    if L.size == 0:
        raise ValueError('pfqn_mva_interval requires at least one queueing station.')
    if L.ndim == 1:
        L = np.column_stack([L, L])
    elif L.shape[1] == 1:
        L = np.column_stack([L[:, 0], L[:, 0]])
    elif L.shape[1] != 2:
        raise ValueError('pfqn_mva_interval is a single-class method: L must be M x 2, '
                         'one [lower upper] demand interval per station.')

    nlo, nup = _endpoints(N, 'population')
    zlo, zup = _endpoints(Z, 'think time')

    Llo = L[:, 0].copy()
    Lup = L[:, 1].copy()
    if np.any(L < 0) or zlo < 0:
        raise ValueError('demands and think times must be nonnegative.')
    if np.any(Llo > Lup) or nlo > nup or zlo > zup:
        raise ValueError('interval lower endpoints must not exceed the upper endpoints.')
    if nlo < 1:
        raise ValueError('pfqn_mva_interval requires a population interval with at least one '
                         'job; the monotonicity theorems assume n >= 1.')
    if abs(nlo - round(nlo)) > _ZERO_TOL or abs(nup - round(nup)) > _ZERO_TOL:
        raise ValueError('population endpoints must be integers; use pfqn_nintmva for a '
                         'nonintegral population.')
    nlo = int(round(nlo))
    nup = int(round(nup))

    M = L.shape[0]
    thick = Lup > Llo + _ZERO_TOL
    X = np.zeros(2)
    Q = np.zeros((M, 2))
    U = np.zeros((M, 2))
    R = np.zeros((M, 2))
    Rtot = np.zeros(2)
    Qtot = np.zeros(2)

    # S1: throughput upper bound, and the upper bounds of the stations whose
    # demand is thin (their own demand is fixed, so lowering the others
    # maximizes them).
    x1, q1, c1 = _mva(Llo, nup, zlo)
    X[1] = x1
    Q[~thick, 1] = q1[~thick]
    R[~thick, 1] = c1[~thick]

    # S2: the same quantities at the opposite corner, giving the lower bounds.
    x2, q2, c2 = _mva(Lup, nlo, zup)
    X[0] = x2
    Q[~thick, 0] = q2[~thick]
    R[~thick, 0] = c2[~thick]

    # S3/S4: the totals increase in the demands and the population and decrease
    # in the think time, so their corners differ from those of the throughput.
    _, q3, c3 = _mva(Llo, nlo, zup)
    Rtot[0] = c3.sum()
    Qtot[0] = q3.sum()
    _, q4, c4 = _mva(Lup, nup, zlo)
    Rtot[1] = c4.sum()
    Qtot[1] = q4.sum()

    # S5/S6: one pair of calls per thick station, its own demand at the endpoint
    # that maximizes (minimizes) it and the others at the opposite endpoint.
    for k in np.flatnonzero(thick):
        d = Llo.copy()
        d[k] = Lup[k]
        _, q5, c5 = _mva(d, nup, zlo)
        Q[k, 1] = q5[k]
        R[k, 1] = c5[k]
        d = Lup.copy()
        d[k] = Llo[k]
        _, q6, c6 = _mva(d, nlo, zup)
        Q[k, 0] = q6[k]
        R[k, 0] = c6[k]

    # U = X*D is not covered by the monotonicity table, so it is enclosed by the
    # product of the two intervals, intersected with the range of a single-server
    # utilization. Where the demand is thin the product is already exact.
    U[:, 0] = X[0] * Llo
    U[:, 1] = np.minimum(1.0, X[1] * Lup)

    return PfqnMvaIntervalResult(X=X, Q=Q, U=U, R=R, Rtot=Rtot, Qtot=Qtot)
