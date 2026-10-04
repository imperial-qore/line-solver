"""First-order LP relaxation of the achievable region of a multiclass open
Markovian queueing network (Bertsimas-Paschalidis-Tsitsiklis 1994).

Native port of matlab/src/api/npfqn/npfqn_bnd_bpt.m, cross-checked against
jar/src/main/java/jline/api/npfqn/Npfqn_bnd_bpt.java.
"""

from typing import Optional, Sequence

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix

__all__ = ['npfqn_bnd_bpt', 'NpfqnBndBptResult']


class NpfqnBndBptResult:
    """Outcome of one achievable-region LP solve.

    Attributes:
        zlb: lower bound on ``sum_r c[r]*x[r]``.
        x: the ``x`` block of the LP optimizer, length K. Only the objective
            value is a bound; an individual ``x[r]`` is a vertex coordinate,
            not a bound on class r, unless ``c`` is the r-th unit vector.
        lambda_: effective arrival rate of each class.
        rho: per-class utilization ``lambda_r/mu_r``.
        rho_station: per-station utilization.
        nvars: number of LP variables.
        nrows: number of LP rows.
    """

    __slots__ = ('zlb', 'x', 'lambda_', 'rho', 'rho_station', 'nvars', 'nrows')

    def __init__(self, zlb, x, lambda_, rho, rho_station, nvars, nrows):
        self.zlb = zlb
        self.x = x
        self.lambda_ = lambda_
        self.rho = rho
        self.rho_station = rho_station
        self.nvars = nvars
        self.nrows = nrows


def npfqn_bnd_bpt(lambda0: Sequence[float], mu: Sequence[float], P,
                  station_of: Sequence[int],
                  c: Optional[Sequence[float]] = None) -> NpfqnBndBptResult:
    """Lower bound on ``sum_r c[r]*x[r]`` over the achievable region.

    ``x[r]`` is the mean sojourn time of class ``r``; the bound is valid for
    EVERY non-idling scheduling policy.

    A "class" here is a buffer with its own exponential service rate and its
    own Markovian routing, so a station serving several customer types owns one
    class per type. The network is open: class ``r`` receives external Poisson
    arrivals at rate ``lambda0[r]`` and, on completing service, becomes class
    ``r'`` with probability ``P[r, r']`` or leaves with the row deficit.

    METHOD. Uniformize the chain and let ``R(t) = sum_r f(r) n_r(t)`` for an
    arbitrary vector ``f``. The steady-state balance of ``E[R^2]`` is an
    identity quadratic in ``f``; since it holds for every ``f``, the two sides'
    coefficient matrices agree entrywise. Diagonal entries give one equation
    per class, off-diagonal entries one per unordered pair, in the variables
    ``x_r``, ``I(r,l) = E[1{sigma(r) busy with r} n_l]`` and
    ``N(i,l) = E[1{station i idle} n_l]``. A third block states that the events
    "station i serves class r" and "station i idle" are mutually exclusive and
    exhaustive, so their terms sum to ``E[n_l] = lambda_l x_l``. Minimizing over
    this polyhedron is a relaxation, hence a lower bound.

    EXACT ON M/M/1: the LP reduces to ``mu*I11 - lambda^2*x = lambda`` and
    ``I11 + N11 = lambda*x`` with ``N11 >= 0``, whence ``x >= 1/(mu-lambda)``
    with equality.

    NOT INCLUDED, DELIBERATELY. The valid inequality ``I(r,r) >= rho_r`` would
    tighten the relaxation but is not part of the reference's characterization,
    and reproducing the reference's published bounds is the acceptance test.

    Args:
        lambda0: external Poisson arrival rate into each class (0 if none).
        mu: exponential service rate of each class.
        P: K x K routing, ``P[r, r']`` = P(class r becomes r' after service);
            row sums must not exceed 1.
        station_of: zero-based station index of each class.
        c: objective weights; ``None`` means all ones.

    Returns:
        NpfqnBndBptResult.

    Reference:
        D. Bertsimas, I. Paschalidis, J. Tsitsiklis (1994). Optimization of
        multiclass queueing networks: polyhedral and nonlinear characterizations
        of achievable performance. Annals of Applied Probability 4(1), 43-75.
        See also D. Bertsimas (1995), Queueing Systems 21, 337-389, Theorem 9,
        which restates the same characterization.
    """
    lambda0 = np.asarray(lambda0, dtype=float).ravel()
    mu = np.asarray(mu, dtype=float).ravel()
    station_of = np.asarray(station_of, dtype=int).ravel()
    P = np.asarray(P, dtype=float)
    K = lambda0.size
    if c is None:
        cost = np.ones(K)
    else:
        cost = np.asarray(c, dtype=float).ravel()
    if mu.size != K or station_of.size != K or cost.size != K:
        raise ValueError(
            "lambda0, mu, station_of and c must all have %d entries." % K)
    if P.shape != (K, K):
        raise ValueError("P must be %dx%d." % (K, K))
    if np.any(mu <= 0):
        raise ValueError("Every class needs a strictly positive service rate.")
    if np.any(P.sum(axis=1) > 1 + 1e-9):
        raise ValueError("The routing matrix has a row summing above one.")
    M = int(station_of.max()) + 1

    # ---- traffic equations, lambda = lambda0 + P' lambda ----
    lam = np.linalg.solve(np.eye(K) - P.T, lambda0)
    if np.any(lam < -1e-9):
        raise ValueError("The traffic equations have no nonnegative solution.")
    lam = np.maximum(lam, 0.0)
    rho = lam / mu
    rho_station = np.zeros(M)
    for r in range(K):
        rho_station[station_of[r]] += rho[r]
    bad = np.nonzero(rho_station >= 1 - 1e-12)[0]
    if bad.size:
        raise ValueError(
            "Station %d is saturated (rho=%.6g): no policy stabilizes the network."
            % (bad[0] + 1, rho_station[bad[0]]))

    # ---- variable layout: x(r) | I(r,l) | N(i,l) ----
    oI = K
    oN = K + K * K
    nv = K + K * K + M * K

    rows = []
    cols = []
    vals = []
    beq = []
    nrow = 0

    def emit(idx, val, rhs):
        nonlocal nrow
        rows.extend([nrow] * len(idx))
        cols.extend(idx)
        vals.extend(val)
        beq.append(rhs)
        nrow += 1

    # (a) diagonal equations, test function n_r^2
    for r in range(K):
        ci = [oI + r * K + r]
        cv = [2 * mu[r]]
        for w in range(K):
            if P[w, r] != 0:
                ci.append(oI + w * K + r)
                cv.append(-2 * mu[w] * P[w, r])
        ci.append(r)
        cv.append(-2 * lambda0[r] * lam[r])
        emit(ci, cv, 2 * lam[r] * (1 - P[r, r]))

    # (b) off-diagonal equations, test function n_r*n_s
    for r in range(1, K):
        for s in range(r):
            ci = [oI + r * K + s, oI + s * K + r]
            cv = [mu[r], mu[s]]
            for w in range(K):
                if P[w, r] != 0:
                    ci.append(oI + w * K + s)
                    cv.append(-mu[w] * P[w, r])
                if P[w, s] != 0:
                    ci.append(oI + w * K + r)
                    cv.append(-mu[w] * P[w, s])
            ci.append(s)
            cv.append(-lambda0[r] * lam[s])
            ci.append(r)
            cv.append(-lambda0[s] * lam[r])
            emit(ci, cv, -lam[r] * P[r, s] - lam[s] * P[s, r])

    # (c) exhaustiveness at each station
    for i in range(M):
        members = np.nonzero(station_of == i)[0]
        for l in range(K):
            ci = [oI + int(r) * K + l for r in members]
            cv = [1.0] * len(ci)
            ci.append(oN + i * K + l)
            cv.append(1.0)
            ci.append(l)
            cv.append(-lam[l])
            emit(ci, cv, 0.0)

    # Duplicate (row, column) entries are SUMMED by coo_matrix, which is the
    # accumulating semantics the reference's sparse() assembly relies on.
    Aeq = coo_matrix((vals, (rows, cols)), shape=(nrow, nv)).tocsr()

    f = np.zeros(nv)
    f[:K] = cost
    res = linprog(f, A_eq=Aeq, b_eq=np.asarray(beq, dtype=float),
                  bounds=(0, None), method='highs')
    if not res.success:
        raise ValueError(
            "The achievable-region LP did not solve to optimality (%s)." % res.message)

    return NpfqnBndBptResult(float(res.fun), np.asarray(res.x[:K], dtype=float),
                             lam, rho, rho_station, nv, nrow)
