"""QD-AMVA: queue-dependent approximate mean value analysis.

Port of ``matlab/src/api/pfqn/pfqn_qdamva.m``, the queue-dependent AMVA of
Casale, Perez and Wang (IFIP PERFORMANCE 2015), on a closed multiclass
product-form network.

A Schweitzer/Bard core in which the class-r demand at station k is scaled by the
queue-dependence term ``g_k`` evaluated at the ARRIVAL-INSTANT total queue
length, ``g = pfqn_lldfun(1 + delta * Q.sum(axis=1), mu)``.

SETTING ``mu`` TO A CONSTANT ROW RECOVERS PLAIN SCHWEITZER AMVA ONLY FOR A
SINGLE CLASS. ``pfqn_lldfun`` does skip a constant row, so ``g == 1`` there, but
the residence time that remains is ``1 + delta * Q.sum(axis=1)`` with ONE
aggregate ``delta = (N.sum()-1)/N.sum()`` applied to the whole arrival-instant
queue, where Bard-Schweitzer shrinks the TAGGED class alone::

    1 + sum_{s != r} Q[k, s] + (N[r]-1)/N[r] * Q[k, r]

The two coincide iff ``K == 1``. Measured over 40 random three-class instances,
``pfqn_qdamva(L, N, Z, ones)`` departs from :func:`pfqn_bs` by up to 0.217 in
absolute queue length, and is the LESS accurate of the two on single-server
multiclass models (mean relative error on Q 0.069 against 0.056 at R = 3), the
aggregate delta buying nothing once ``g == 1``. This is the QD-AMVA closure, not
a defect of the port, but do not use the function as a Schweitzer oracle for
``K > 1``.

MU IS A DIMENSIONLESS RATE MULTIPLIER, NOT A RATE. ``mu[k, n]`` is the factor by
which station k serves faster when it holds n jobs. Two traps follow from
``pfqn_lldfun`` and are the reference's, not this port's:

- it SKIPS a station whose mu row is constant (its ``range(...) > 0`` gate), so
  a single-server station must be a row of ones and a c-server station
  ``minimum(1..smax, c)``. Passing a c-server station a constant row silently
  returns ``g = 1``, i.e. a single server.
- ``smax = mu.shape[1]`` must be at least ``ceil(sum(N))`` or the interpolation
  clamps the population and the top of the rate curve is never reached.

Delay stations are carried in ``Z``, not as rows of ``L``. Closed classes only:
an infinite ``N[r]`` is not supported.

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

from typing import Optional, Tuple

import numpy as np

from .utils import pfqn_lldfun

__all__ = ['pfqn_qdamva']


def pfqn_qdamva(L: np.ndarray,
                N: np.ndarray,
                Z: Optional[np.ndarray] = None,
                mu: Optional[np.ndarray] = None,
                Q0: Optional[np.ndarray] = None,
                tol: float = 1e-6,
                maxiter: int = 10000
                ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, int, np.ndarray]:
    """QD-AMVA on a closed multiclass product-form network.

    Args:
        L: (M x R) service demand matrix.
        N: (R,) population vector, finite.
        Z: (R,) think time vector; None means no think time.
        mu: (M x smax) queue-dependent rate multipliers; None means none.
        Q0: (M x R) initial guess; None means the reference's demand split.
        tol: convergence tolerance on the queue lengths.
        maxiter: maximum number of iterations.

    Returns:
        Q: (M x R) mean queue lengths.
        X: (R,) per-class throughputs.
        U: (M x R) per-class utilizations, carrying the g scaling.
        iter: number of iterations performed.
        R: (M x R) per-class residence times, ``Q = X * R``.
    """
    L = np.atleast_2d(np.asarray(L, dtype=float))
    M, K = L.shape
    N = np.asarray(N, dtype=float).reshape(-1)
    if N.size != K:
        raise ValueError('pfqn_qdamva: the population vector must have one entry per class')
    if Z is None:
        Z = np.zeros(K)
    Z = np.asarray(Z, dtype=float).reshape(-1)
    if Z.size != K:
        raise ValueError('pfqn_qdamva: the think-time vector must have one entry per class')
    if mu is not None:
        mu = np.atleast_2d(np.asarray(mu, dtype=float))

    X = np.zeros(K)
    U = np.zeros((M, K))
    R = np.zeros((M, K))
    it = 0

    Ntot = float(np.sum(N))
    if not Ntot > 0.0:
        # delta is undefined on an empty population, and the reference returns
        # the zero queue rather than dividing by it.
        return np.zeros((M, K)), X, U, it, R

    if Q0 is None:
        # Ltot = 0 for a class with no demand anywhere: L/Ltot is a NaN the
        # iteration never recovers from. Such a column arises routinely in a
        # layered fixed point, where a caller can start with no work at the
        # layer station, so the column is left at zero instead.
        Ltot = L.sum(axis=0)
        Q = np.zeros((M, K))
        nz = Ltot > 0
        if np.any(nz):
            Q[:, nz] = L[:, nz] / Ltot[nz] * N[nz]
    else:
        Q = np.array(Q0, dtype=float).reshape(M, K)

    delta = (Ntot - 1.0) / Ntot

    # Q*10 as the sentinel, as the reference notes, stalls on an all-zero seed:
    # the loop would exit before its first pass. Offset instead.
    Q_1 = Q + 10.0 * (1.0 + tol)
    while np.max(np.abs(Q - Q_1)) > tol and it < maxiter:
        it += 1
        Q_1 = Q.copy()

        # The arrival-instant total queue length, class independent.
        Ak = 1.0 + delta * Q.sum(axis=1)
        g = np.asarray(pfqn_lldfun(Ak, mu)).reshape(-1)

        for r in range(K):
            R[:, r] = L[:, r] * g * (1.0 + delta * Q.sum(axis=1))
            denom = Z[r] + float(np.sum(R[:, r]))
            X[r] = N[r] / denom if denom > 0.0 else 0.0
            Q[:, r] = X[r] * R[:, r]
            U[:, r] = L[:, r] * g * X[r]

    return Q, X, U, it, R
