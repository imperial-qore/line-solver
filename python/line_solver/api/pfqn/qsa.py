"""Queue-Shift Approximation (QSA) for closed product-form networks.

Port of ``matlab/src/api/pfqn/pfqn_qsa.m``.

Schweitzer, Serazzi and Broglia, "A Queue-Shift Approximation Technique for
Product-Form Queueing Networks", Tools'98, LNCS 1469, pp. 267-279. QSA
approximates the arrival-instant queue lengths through the absolute shift of
the aggregate queue length,

    Y_ri(K) = 1 + Q_i(K - e_r) - Q_i(K)      i in QC

in place of the fractional deviations of Linearizer, so the unknowns are one
per station rather than one per station-class. The core equation (13a) is
imposed at K, at every K - e_s and, in the three-level variant of eq. (16), at
every K - e_s - e_t with the affine extrapolation of eq. (15).
"""

from typing import Optional, Tuple

import numpy as np

from ...constants import SchedStrategy
from .mva import pfqn_bs

__all__ = ['pfqn_qsa']


def pfqn_qsa(L, N, Z=None, type=None, tol: float = 1e-10,
             maxiter: int = 100, levels: int = 3, QN0=None
             ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray,
                        np.ndarray, int]:
    """Queue-Shift Approximation for a closed product-form network.

    Args:
        L: Service demand matrix (M x R)
        N: Population vector (R,)
        Z: Think time vector (default 0)
        type: Scheduling strategy per station; ``SchedStrategy.INF`` marks a
            delay centre, whose demand enters the cycle time without a queueing
            term (the paper's DC set)
        tol: Residual tolerance of the Newton iteration (default 1e-10)
        maxiter: Maximum Newton iterations (default 100)
        levels: 2 for the two-level QSA of eq. (14), 3 for eq. (16)
        QN0: Warm start for the Bard-Schweitzer initialization (M x R)

    Returns:
        Tuple of (Q, U, W, C, X, totiter); W holds residence times and C cycle
        times, matching ``pfqn_linearizer``.
    """
    L = np.asarray(L, dtype=np.float64)
    N = np.asarray(N, dtype=np.float64).flatten()
    R = len(N)
    if L.ndim == 1:
        L = L.reshape(-1, 1) if R == 1 else L.reshape(1, -1)
    M = L.shape[0]

    Z = np.zeros(R) if Z is None else np.asarray(Z, dtype=np.float64)
    if Z.ndim > 1:
        Z = Z.sum(axis=0)
    Z = Z.flatten()

    if type is None:
        is_qc = np.ones(M, dtype=bool)
    else:
        tarr = np.asarray(type).flatten()[:M]
        is_qc = np.array([t != SchedStrategy.INF for t in tarr], dtype=bool)

    Q = np.zeros((M, R))
    U = np.zeros((M, R))
    W = np.zeros((M, R))
    C = np.zeros(R)
    X = np.zeros(R)

    if L.size == 0 or np.all(L.max(axis=0) == 0) or np.all(N <= 0):
        for r in range(R):
            if N[r] > 0 and Z[r] > 0:
                X[r] = N[r] / Z[r]
            U[:, r] = X[r] * L[:, r]
        return Q, U, W, C, X, 0

    # Populations touched by (16): K, every K - e_s, every K - e_s - e_t.
    pops = [N.copy()]
    s_idx = np.zeros(R, dtype=int) - 1
    p_idx = np.zeros((R, R), dtype=int) - 1
    for s in range(R):
        n = N.copy()
        n[s] -= 1.0
        if np.all(n >= 0):
            pops.append(n)
            s_idx[s] = len(pops) - 1
    if levels >= 3:
        for s in range(R):
            for t in range(s, R):
                n = N.copy()
                n[s] -= 1.0
                n[t] -= 1.0
                if np.all(n >= 0):
                    pops.append(n)
                    p_idx[s, t] = len(pops) - 1
                    p_idx[t, s] = p_idx[s, t]
    pops = np.array(pops)
    nP = pops.shape[0]

    # Bard-Schweitzer at every population supplies the Newton starting point.
    q = np.zeros((M, nP))
    for p in range(nP):
        q[:, p] = _aggbs(L, pops[p, :], Z, is_qc, QN0)
    qc = np.flatnonzero(is_qc)
    Ldc = L[~is_qc, :].sum(axis=0)

    x = q[np.ix_(qc, np.arange(nP))].reshape(-1).copy()
    n_unk = x.size
    F, adm = _resid(x, L, Z, pops, s_idx, p_idx, qc, Ldc, levels)
    fnrm = np.linalg.norm(F)
    totiter = 0
    for totiter in range(1, maxiter + 1):
        if fnrm < tol:
            totiter -= 1
            break
        J = np.zeros((n_unk, n_unk))
        for k in range(n_unk):
            h = 1e-7 * max(1.0, abs(x[k]))
            xp = x.copy()
            xp[k] += h
            Fp, _ = _resid(xp, L, Z, pops, s_idx, p_idx, qc, Ldc, levels)
            J[:, k] = (Fp - F) / h
        try:
            dx = np.linalg.solve(J, -F)
        except np.linalg.LinAlgError:
            dx = -np.linalg.pinv(J) @ F
        if not np.all(np.isfinite(dx)):
            dx = -np.linalg.pinv(J) @ F
        accepted = False
        lam = 1.0
        for _ in range(40):
            xn = x + lam * dx
            Fn, admn = _resid(xn, L, Z, pops, s_idx, p_idx, qc, Ldc, levels)
            if admn and np.linalg.norm(Fn) < fnrm:
                x, F, fnrm, adm = xn, Fn, np.linalg.norm(Fn), admn
                accepted = True
                break
            lam /= 2.0
        if not accepted:
            break

    # Disaggregate (13) at K into the per-class measures
    q[np.ix_(qc, np.arange(nP))] = x.reshape(len(qc), nP)
    Y0 = _shift(q, pops, s_idx, p_idx, -1, -1, levels)
    for r in range(R):
        if N[r] < 1:
            continue
        W[is_qc, r] = L[is_qc, r] * (q[is_qc, 0] + Y0[is_qc, r])
        W[~is_qc, r] = L[~is_qc, r]
        X[r] = N[r] / (Z[r] + W[:, r].sum())
        Q[:, r] = X[r] * W[:, r]
        U[:, r] = X[r] * L[:, r]
        C[r] = N[r] / X[r] - Z[r]
    return Q, U, W, C, X, totiter


def _resid(x, L, Z, pops, s_idx, p_idx, qc, Ldc, levels):
    """Residual of (13) imposed simultaneously at every population of (16).

    The second return value flags the side conditions of Remark 2 (non-negative
    queue lengths, positive cycle times).
    """
    M, R = L.shape
    nP = pops.shape[0]
    mq = len(qc)
    q = np.zeros((M, nP))
    q[np.ix_(qc, np.arange(nP))] = x.reshape(mq, nP)
    F = np.zeros((mq, nP))
    adm = bool(np.all(x >= 0))
    for p in range(nP):
        np_ = pops[p, :]
        s, t = _which(p, s_idx, p_idx)
        Y = _shift(q, pops, s_idx, p_idx, s, t, levels)
        A = q[qc, p][:, None] + Y[qc, :]   # 1 + Q_i(K - e_r) at the arrival instant
        acc = np.zeros(mq)
        for r in range(R):
            if np_[r] < 1:
                continue
            c = Z[r] + float(L[qc, r] @ A[:, r]) + Ldc[r]
            if not (c > 0) or not np.isfinite(c):
                adm = False
                c = np.finfo(float).eps
            acc += (np_[r] / c) * L[qc, r] * A[:, r]
        F[:, p] = q[qc, p] - acc
    return F.reshape(-1), adm


def _shift(q, pops, s_idx, p_idx, s, t, levels):
    """Shift matrix (M x R) of (16d)-(16e), or (15) when both s and t are set."""
    M = q.shape[0]
    R = pops.shape[1]
    Y = np.zeros((M, R))
    if s < 0:
        for r in range(R):
            if s_idx[r] >= 0:
                Y[:, r] = 1.0 + q[:, s_idx[r]] - q[:, 0]
    elif t < 0:
        if levels < 3:
            return _shift(q, pops, s_idx, p_idx, -1, -1, levels)  # (14)
        for r in range(R):
            if p_idx[s, r] >= 0 and pops[s_idx[s], r] >= 1:
                Y[:, r] = 1.0 + q[:, p_idx[s, r]] - q[:, s_idx[s]]
    else:
        Y = (_shift(q, pops, s_idx, p_idx, s, -1, levels)
             + _shift(q, pops, s_idx, p_idx, t, -1, levels)
             - _shift(q, pops, s_idx, p_idx, -1, -1, levels))
    return Y


def _which(p, s_idx, p_idx):
    """Decode a population index into the removed classes."""
    if p == 0:
        return -1, -1
    k = np.flatnonzero(s_idx == p)
    if k.size > 0:
        return int(k[0]), -1
    ss, tt = np.nonzero(p_idx == p)
    return int(ss[0]), int(tt[0])


def _aggbs(L, n, Z, is_qc, QN0: Optional[np.ndarray]):
    """Aggregate Bard-Schweitzer queue lengths at population n, with the
    delay-centre demands folded into the think time."""
    M, R = L.shape
    q = np.zeros(M)
    n = np.maximum(np.asarray(n, dtype=np.float64), 0.0)
    if np.all(n <= 0):
        return q
    Zeff = Z + L[~is_qc, :].sum(axis=0)
    if is_qc.any():
        if QN0 is None or np.size(QN0) == 0:
            Xb, QN, _, _, _ = pfqn_bs(L[is_qc, :], n, Zeff)
        else:
            QN0 = np.asarray(QN0, dtype=np.float64).reshape(M, R)
            Xb, QN, _, _, _ = pfqn_bs(L[is_qc, :], n, Zeff, 1e-6, 1000,
                                      QN0[is_qc, :])
        q[is_qc] = QN.sum(axis=1)
        Xb = np.asarray(Xb).flatten()
    else:
        Xb = np.zeros(R)
        for r in range(R):
            if n[r] >= 1 and Zeff[r] > 0:
                Xb[r] = n[r] / Zeff[r]
    for i in np.flatnonzero(~is_qc):
        q[i] = float(Xb @ L[i, :])
    return q
