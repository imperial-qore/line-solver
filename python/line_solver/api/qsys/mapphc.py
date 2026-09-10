"""
Exact MAP/PH/c FCFS queue.

THE STATE SPACE. With c identical servers the server identities carry no
information, so the service phases are held as a MULTISET: a configuration is
n = (n_1..n_ms) with sum(n) = k servers busy in phase i. There are
comb(ms+k-1, k) of them, the count of Asmussen and Moller (2001), against ms^k
for the ordered space. Levels 0..c-1 are the boundary (level = servers busy),
levels >= c repeat and carry the queue.

THE WAITING TIME. An arrival that finds j customers waiting ahead of it waits
for j+1 service completions, so Wq is the (j+1)-st event time of the
configuration MAP (Lc, Cdep) started at the arrival-epoch configuration.
Folding the matrix-geometric level distribution over j gives the LINEAR matrix
ODE G'(t) = G Lj + R G Cj with G(0) = (I-R)^-1 kron(D1,I)/lambda and
P(Wq > t) = pi_c G(t) e, so Wq is matrix-exponential. Its transform obeys the
generalized Sylvester equation g(sI-Lj) - R g Cj = G(0), and every moment reuses
that one operator with a different right-hand side.

References:
    S. Asmussen and J.R. Moller, "Calculation of the steady state waiting time
    distribution in GI/PH/c and MAP/PH/c queues", Queueing Systems 37(1):9-29,
    2001.
    D.P. Gaver, P.A. Jacobs, G. Latouche, "Finite birth-and-death models in
    randomly changing environments", Adv. Appl. Probab. 16:715-731, 1984.
    Original MATLAB: matlab/src/api/qsys/qsys_mapphc.m
"""

from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np
from scipy.linalg import expm

from ..mam.qbd import qbd_R_logred
from ..mam.ldqbd_mphc import ph_multisets


@dataclass
class MapPhcResult:
    """Result structure for the exact MAP/PH/c analysis."""
    meanQueueLength: float
    meanWaitingTime: float
    meanSojournTime: float
    utilization: float
    queueLengthDist: np.ndarray
    waitingTimeMoments: np.ndarray
    waitingTimeCCDF: Optional[np.ndarray]
    waitingTimePoints: Optional[np.ndarray]
    probWait: float
    phaseCount: int
    analyzer: str


def _find(rows: np.ndarray, key: np.ndarray) -> int:
    return int(np.flatnonzero((rows == key).all(axis=1))[0])


def _stat_left_null(G: np.ndarray) -> np.ndarray:
    """
    Left null vector of G normalized to sum one. The R-corrected level-c block
    has nonzero row sums, so the augmented system is used rather than a
    generator solve.
    """
    n = G.shape[0]
    B = np.hstack([G, np.ones((n, 1))])
    y = np.zeros(n + 1)
    y[n] = 1.0
    x, _, _, _ = np.linalg.lstsq(B.T, y, rcond=None)
    return x


def qsys_mapphc(D0, D1, alpha, S, c: int, max_num_comp: int = 500,
                num_w_moms: int = 3, w_points: Optional[Sequence[float]] = None) -> MapPhcResult:
    """
    Analyze a MAP/PH/c FCFS queue exactly.

    Args:
        D0, D1: arrival MAP of order ma
        alpha: PH service initial vector (ms,)
        S: PH service sub-generator (ms, ms)
        c: number of servers, identical so the service law is shared
        max_num_comp: cap on the queue length probabilities returned
        num_w_moms: how many waiting-time moments to return
        w_points: times at which to evaluate P(Wq > t)

    Returns:
        MapPhcResult
    """
    D0 = np.atleast_2d(np.asarray(D0, dtype=float))
    D1 = np.atleast_2d(np.asarray(D1, dtype=float))
    alpha = np.asarray(alpha, dtype=float).ravel()
    S = np.atleast_2d(np.asarray(S, dtype=float))
    c = int(c)
    ma = D0.shape[0]
    ms = S.shape[0]
    if D0.shape[1] != ma or D1.shape != (ma, ma):
        raise ValueError("D0 and D1 must be square and of equal order")
    if S.shape[1] != ms or alpha.size != ms:
        raise ValueError("alpha and S must have matching order")
    if c < 1:
        raise ValueError("The number of servers c must be at least one")

    s0 = -S @ np.ones(ms)

    theta = _stat_left_null(D0 + D1)
    lam = float(theta @ D1 @ np.ones(ma))
    mean_service = float(-alpha @ np.linalg.solve(S, np.ones(ms)))
    rho = lam * mean_service / c
    if rho >= 1.0:
        raise ValueError("The load %g of the system is not below one" % rho)

    cfg = [ph_multisets(ms, k) for k in range(c + 1)]

    Lcfg, Up, Dn = [], [], []
    for k in range(c + 1):
        Ck = cfg[k]
        nk = Ck.shape[0]
        Lk = np.zeros((nk, nk))
        for row in range(nk):
            n = Ck[row]
            for i in range(ms):
                if n[i] == 0:
                    continue
                for j in range(ms):
                    if j == i:
                        continue
                    m = n.copy(); m[i] -= 1; m[j] += 1
                    Lk[row, _find(Ck, m)] += n[i] * S[i, j]
                Lk[row, row] += n[i] * S[i, i]
        Lcfg.append(Lk)

        if k < c:
            Ck1 = cfg[k + 1]
            Uk = np.zeros((nk, Ck1.shape[0]))
            for row in range(nk):
                for j in range(ms):
                    m = Ck[row].copy(); m[j] += 1
                    Uk[row, _find(Ck1, m)] += alpha[j]
            Up.append(Uk)
        else:
            Up.append(None)

        if k > 0:
            Ckm = cfg[k - 1]
            Dk = np.zeros((nk, Ckm.shape[0]))
            for row in range(nk):
                n = Ck[row]
                for i in range(ms):
                    if n[i] == 0:
                        continue
                    m = n.copy(); m[i] -= 1
                    Dk[row, _find(Ckm, m)] += n[i] * s0[i]
            Dn.append(Dk)
        else:
            Dn.append(None)

    # Completion WITH an immediate restart: the repeating down block
    Cc = cfg[c]
    nc = Cc.shape[0]
    Cdep = np.zeros((nc, nc))
    for row in range(nc):
        n = Cc[row]
        for i in range(ms):
            if n[i] == 0:
                continue
            for j in range(ms):
                m = n.copy(); m[i] -= 1; m[j] += 1
                Cdep[row, _find(Cc, m)] += n[i] * s0[i] * alpha[j]

    Ima = np.eye(ma)
    Inc = np.eye(nc)
    A_up = np.kron(D1, Inc)
    A_loc = np.kron(D0, Inc) + np.kron(Ima, Lcfg[c])
    A_dn = np.kron(Ima, Cdep)
    R = qbd_R_logred(A_dn, A_loc, A_up)

    # Boundary levels 0..c, with the tail folded into level c through R
    sz = [ma * cfg[k].shape[0] for k in range(c + 1)]
    off = np.concatenate([[0], np.cumsum(sz)])
    tot = int(off[-1])
    Q = np.zeros((tot, tot))
    for k in range(c + 1):
        r0, r1 = off[k], off[k + 1]
        Ick = np.eye(cfg[k].shape[0])
        if k < c:
            Q[r0:r1, r0:r1] = np.kron(D0, Ick) + np.kron(Ima, Lcfg[k])
            Q[r0:r1, off[k + 1]:off[k + 2]] = np.kron(D1, Up[k])
        else:
            Q[r0:r1, r0:r1] = A_loc + R @ A_dn
        if k > 0:
            Q[r0:r1, off[k - 1]:off[k]] = np.kron(Ima, Dn[k])
    pi_vec = _stat_left_null(Q)

    n_op = nc * ma
    ImR_inv = np.linalg.inv(np.eye(n_op) - R)
    pi_c = pi_vec[off[c]:off[c + 1]]
    mass_boundary = float(np.sum(pi_vec[:off[c]]))
    mass_tail = float(pi_c @ ImR_inv @ np.ones(n_op))
    pi_vec = pi_vec / (mass_boundary + mass_tail)
    pi_c = pi_vec[off[c]:off[c + 1]]

    ql = [float(np.sum(pi_vec[off[k]:off[k + 1]])) for k in range(c)]
    tail = pi_c.copy()
    ql.append(float(np.sum(tail)))
    acc = float(np.sum(ql))
    while acc < 1 - 1e-12 and len(ql) < max_num_comp:
        tail = tail @ R
        ql.append(float(np.sum(tail)))
        acc += ql[-1]
    ql = np.asarray(ql)
    # E[N] in CLOSED FORM. max_num_comp caps the probabilities RETURNED, not the
    # mean: summing the truncated list loses the matrix-geometric tail, which at
    # rho -> 1 carries a first-order share of the mass. With pi_{c+j} = pi_c R^j,
    # sum_j (c+j) pi_c R^j e = pi_c [c (I-R)^-1 + R (I-R)^-2] e.
    u = ImR_inv @ np.ones(n_op)
    mean_ql = float(np.sum(np.arange(c) * ql[:c])
                    + pi_c @ (c * u + R @ (ImR_inv @ u)))

    # Waiting time
    Lj = np.kron(Ima, Lcfg[c])
    Cj = np.kron(Ima, Cdep)
    G0 = ImR_inv @ np.kron(D1, Inc) / lam
    prob_wait = float(pi_c @ G0 @ np.ones(n_op))

    Iop = np.eye(n_op)
    # X(-Lj) - R X Cj = rhs, vectorized column-major
    Kop = np.kron((-Lj).T, Iop) - np.kron(Cj.T, R)
    w_moms = np.zeros(max(num_w_moms, 1))
    g_prev = None
    for k in range(1, num_w_moms + 1):
        rhs = G0 if k == 1 else -(k - 1) * g_prev
        g = np.linalg.solve(Kop, rhs.reshape(-1, order='F')).reshape(n_op, n_op, order='F')
        w_moms[k - 1] = k * ((-1.0) ** (k - 1)) * float(pi_c @ g @ np.ones(n_op))
        g_prev = g
    mean_wt = float(w_moms[0]) if num_w_moms >= 1 else float('nan')

    w_ccdf = None
    pts = None
    if w_points is not None and len(np.atleast_1d(w_points)) > 0:
        pts = np.atleast_1d(np.asarray(w_points, dtype=float))
        Kt = np.kron(Lj.T, Iop) + np.kron(Cj.T, R)
        v0 = G0.reshape(-1, order='F')
        w_ccdf = np.zeros(pts.size)
        for it, t in enumerate(pts):
            Gt = (expm(Kt * t) @ v0).reshape(n_op, n_op, order='F')
            w_ccdf[it] = float(pi_c @ Gt @ np.ones(n_op))

    return MapPhcResult(
        meanQueueLength=mean_ql,
        meanWaitingTime=mean_wt,
        meanSojournTime=mean_wt + mean_service,
        utilization=rho,
        queueLengthDist=ql,
        waitingTimeMoments=w_moms,
        waitingTimeCCDF=w_ccdf,
        waitingTimePoints=pts,
        probWait=prob_wait,
        phaseCount=nc,
        analyzer="LINE:MAP/PH/%d" % c,
    )
