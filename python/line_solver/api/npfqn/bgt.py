"""Piecewise-linear Lyapunov upper bound on the queue lengths of a multitype
Markovian queueing network (Bertsimas-Gamarnik-Tsitsiklis 2001).

Native port of matlab/src/api/npfqn/npfqn_bnd_bgt.m, cross-checked against
jar/src/main/java/jline/api/npfqn/Npfqn_bnd_bgt.java.
"""

from typing import List, Optional, Sequence

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix, vstack as sp_vstack

__all__ = ['npfqn_bnd_bgt', 'NpfqnBndBgtResult']


class NpfqnBndBgtResult:
    """Outcome of one GLP[dm] solve and of the Theorem 4 bound.

    Attributes:
        Qub: list of arrays, ``Qub[i][k]`` upper bounds ``E[Q(i,k)]``; ``inf``
            where the LP optimum leaves ``max_j L^j(i,k) = 0``.
        gamma: the drift certificate; strictly positive on success.
        Lmax: ``max`` over ``j`` and ``(i,k)`` of ``L``.
        L: the Lyapunov coefficients, shape ``(J, N)``.
        V: the per-station slack ``V_j``.
        B: the exception parameter of the smoothed Lyapunov function.
        U: the Theorem 4 bound on ``E[L^j'Q]``, the same for every ``j``.
        tail_ratio: geometric decay ratio of the tail bound.
        tail_step: step of the tail bound, ``2(Lmax+gamma/2)``.
        rho: per-class nominal load.
        rho_station: per-station nominal load.
        scale: the uniformization divisor applied to lambda and mu.
        class_type, class_stage, class_station: per-class index maps.
    """

    __slots__ = ('Qub', 'gamma', 'Lmax', 'L', 'V', 'B', 'U', 'tail_ratio', 'tail_step',
                 'rho', 'rho_station', 'scale', 'class_type', 'class_stage', 'class_station')

    def __init__(self, **kw):
        for k in self.__slots__:
            setattr(self, k, kw.get(k))


def npfqn_bnd_bgt(lambda_: Sequence[float], mu: Sequence[Sequence[float]],
                  sigma: Sequence[Sequence[int]],
                  J: Optional[int] = None) -> NpfqnBndBgtResult:
    """Upper bound the steady-state queue lengths of a multitype network.

    The bound is valid for EVERY work-conserving Markovian policy.

    MODEL. ``J`` single-server stations; ``I`` customer types; type ``i``
    arrives as a Poisson stream of rate ``lambda_[i]`` and passes through stages
    ``k = 0..len(mu[i])-1``, stage ``k`` being served at station ``sigma[i][k]``
    at exponential rate ``mu[i][k]``. Class ``(i,k)`` is the buffer of type
    ``i`` at stage ``k``; ``N = sum_i len(mu[i])`` is the number of classes.

    METHOD. Solve the Down-Meyn global-stability linear program GLP[dm], eq.
    (25)-(28) of the reference, in the piecewise-linear Lyapunov function
    ``Phi(x) = max_j L^j'x``::

        L^j(i,1) lambda_i + mu(i,k) (L^j(i,k+1) - L^j(i,k)) + V_j <= -gamma
                                                     for (i,k) in station j
        mu(i,k) (L^j(i,k+1) - L^j(i,k)) <= V_j       for (i,k) not in j
        (1/(J-1)) sum_{j' != j} L^j'(i,k) >= L^j(i,k)  for (i,k) not in j
        L, V, gamma >= 0

    with ``L^j(i,Ji+1) = 0``. A feasible solution with ``gamma > 0`` certifies
    that EVERY work-conserving policy is stable, and a smoothed ``Phi`` is then
    a Lyapunov function with drift ``gamma/4`` and an explicit exception
    parameter, giving the reference's Theorem 4 bound::

        E[L^j'Q] <= 16 N J^2 (J-1) (Lmax+gamma)^3 / gamma^2
                    + 8 (Lmax + gamma/2)^2 / gamma  =: U

    for every ``j``, whence ``E[Q(i,k)] <= U / max_j L^j(i,k)``.

    THE RATES ARE RESCALED so that ``sum_i lambda_i + sum_{i,k} mu(i,k) = 1``,
    the uniformization the reference imposes before Theorem 4. Queue lengths are
    counts and are unaffected by the time scale.

    NORMALIZATION, WHICH THE REFERENCE LEAVES OPEN. GLP[dm] is homogeneous and
    so is the bound, so this routine fixes ``L^j(i,k) <= 1`` and MAXIMIZES
    ``gamma``, then breaks ties among gamma-optimal solutions by maximizing
    ``sum L``: a degenerate optimum can otherwise zero some ``L^j(i,k)`` and
    report an infinite bound for a class for no reason.

    THE BOUND IS LOOSE, and knowingly so: the exception parameter carries
    ``(Lmax+gamma)^3/gamma^2`` and dominates as soon as ``J > 1``. What is sharp
    is the STABILITY CERTIFICATE ``gamma > 0`` and the geometric tail RATE.

    Reference:
        D. Bertsimas, D. Gamarnik, J. N. Tsitsiklis (2001). Performance of
        multiclass Markovian queueing networks via piecewise linear Lyapunov
        functions. Annals of Applied Probability 11(4), 1384-1428, Section 5.1
        (GLP[dm] of Down and Meyn 1997, and Theorem 4).
    """
    lam_in = np.asarray(lambda_, dtype=float).ravel()
    I = lam_in.size
    if len(mu) != I or len(sigma) != I:
        raise ValueError("mu and sigma must both have %d entries." % I)
    if J is None:
        J = 1 + max(int(np.max(np.asarray(sigma[i], dtype=int))) for i in range(I))

    # ---- flatten (i,k) into a class index ----
    class_type: List[int] = []
    class_stage: List[int] = []
    class_station: List[int] = []
    muc: List[float] = []
    first_of = np.zeros(I, dtype=int)
    for i in range(I):
        mi = np.asarray(mu[i], dtype=float).ravel()
        si = np.asarray(sigma[i], dtype=int).ravel()
        if mi.size != si.size:
            raise ValueError("mu[%d] and sigma[%d] have different lengths." % (i, i))
        if mi.size == 0:
            raise ValueError("Type %d has no stage." % (i + 1))
        if not lam_in[i] > 0:
            raise ValueError("Every type needs a strictly positive arrival rate.")
        first_of[i] = len(class_type)
        for k in range(mi.size):
            class_type.append(i)
            class_stage.append(k)
            class_station.append(int(si[k]))
            muc.append(float(mi[k]))
    N = len(muc)
    class_type = np.asarray(class_type, dtype=int)
    class_stage = np.asarray(class_stage, dtype=int)
    class_station = np.asarray(class_station, dtype=int)
    muc = np.asarray(muc, dtype=float)
    if np.any(muc <= 0):
        raise ValueError("Every stage needs a strictly positive service rate.")
    next_of = np.full(N, -1, dtype=int)
    for c in range(N - 1):
        if class_type[c + 1] == class_type[c]:
            next_of[c] = c + 1

    # ---- loads ----
    rho = lam_in[class_type] / muc
    rho_station = np.zeros(J)
    for c in range(N):
        rho_station[class_station[c]] += rho[c]
    bad = np.nonzero(rho_station >= 1)[0]
    if bad.size:
        raise ValueError("Station %d is saturated (rho=%.6g): the load condition of the "
                         "reference fails." % (bad[0] + 1, rho_station[bad[0]]))

    # ---- uniformization ----
    scale = float(lam_in.sum() + muc.sum())
    lam = lam_in / scale
    mus = muc / scale

    # ---- LP layout: L(j,c) -> j*N + c ; V(j) -> J*N + j ; gamma -> J*N+J ----
    oV = J * N
    ig = J * N + J
    nv = J * N + J + 1

    rows: List[int] = []
    cols: List[int] = []
    vals: List[float] = []
    nrow = 0

    def emit(idx, val):
        nonlocal nrow
        rows.extend([nrow] * len(idx))
        cols.extend(idx)
        vals.extend(val)
        nrow += 1

    for j in range(J):
        for c in range(N):
            if class_station[c] == j:
                ci = [j * N + int(first_of[class_type[c]]), j * N + c]
                cv = [lam[class_type[c]], -mus[c]]
                if next_of[c] >= 0:
                    ci.append(j * N + int(next_of[c]))
                    cv.append(mus[c])
                ci.extend([oV + j, ig])
                cv.extend([1.0, 1.0])
                emit(ci, cv)
            else:
                ci = [j * N + c]
                cv = [-mus[c]]
                if next_of[c] >= 0:
                    ci.append(j * N + int(next_of[c]))
                    cv.append(mus[c])
                ci.append(oV + j)
                cv.append(-1.0)
                emit(ci, cv)
                if J > 1:
                    ci = [j * N + c]
                    cv = [1.0]
                    for jp in range(J):
                        if jp != j:
                            ci.append(jp * N + c)
                            cv.append(-1.0 / (J - 1))
                    emit(ci, cv)

    # Duplicate (row, column) entries are SUMMED by coo_matrix, the accumulating
    # semantics the reference's sparse() assembly relies on.
    A = coo_matrix((vals, (rows, cols)), shape=(nrow, nv)).tocsr()
    b = np.zeros(nrow)

    bounds = [(0.0, 1.0)] * (J * N) + [(0.0, None)] * (J + 1)
    f = np.zeros(nv)
    f[ig] = -1.0
    res = linprog(f, A_ub=A, b_ub=b, bounds=bounds, method='highs')
    if not res.success:
        raise ValueError("GLP[dm] did not solve to optimality (%s)." % res.message)
    gamma = -float(res.fun)
    if not gamma > 0:
        raise ValueError("GLP[dm] has no solution with gamma > 0: this network is not certified "
                         "globally stable, so no finite piecewise-linear Lyapunov bound exists.")

    # Tie-break among gamma-optimal solutions: maximize sum L, so a degenerate
    # vertex does not report an infinite bound for a class it zeroed arbitrarily.
    pin = coo_matrix(([-1.0], ([0], [ig])), shape=(1, nv))
    A2 = sp_vstack([A, pin]).tocsr()
    b2 = np.concatenate([b, [-gamma]])
    f2 = np.zeros(nv)
    f2[:J * N] = -1.0
    x = np.asarray(res.x, dtype=float)
    res2 = linprog(f2, A_ub=A2, b_ub=b2, bounds=bounds, method='highs')
    if res2.success:
        x = np.asarray(res2.x, dtype=float)
        gamma = float(x[ig])

    L = x[:J * N].reshape(J, N)
    V = x[oV:oV + J]
    Lmax = float(L.max())

    # ---- Theorem 4 ----
    B = 16.0 * N * J * J * (J - 1) * (Lmax + gamma) ** 3 / gamma ** 2
    U = B + 8.0 * (Lmax + gamma / 2) ** 2 / gamma
    tail_step = 2 * (Lmax + gamma / 2)
    tail_ratio = (Lmax + gamma / 2) / (Lmax + 0.75 * gamma)

    best = L.max(axis=0)
    qub = np.full(N, np.inf)
    pos = best > 0
    qub[pos] = U / best[pos]

    Qub = [qub[class_type == i] for i in range(I)]
    return NpfqnBndBgtResult(Qub=Qub, gamma=gamma, Lmax=Lmax, L=L, V=V, B=B, U=U,
                             tail_ratio=tail_ratio, tail_step=tail_step, rho=rho,
                             rho_station=rho_station, scale=scale, class_type=class_type,
                             class_stage=class_stage, class_station=class_station)
