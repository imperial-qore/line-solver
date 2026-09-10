"""Discrete-time (slotted) matrix-analytic primitives and queues.

Everything here measures time in SLOTS and follows the late arrival system with
delayed access (LAS-DA): within a slot the service completion resolves first,
arrivals are appended at the end of the slot and cannot enter service before the
next one, and the level is read after both. This is the convention of the Q-MAM
discrete-time queues and of the LDES slotted engine, so the three codebases are
directly comparable.

A discrete phase-type law is a pair (alpha, A) with
P[X=k] = alpha @ A^(k-1) @ a, a = e - A @ e, k = 1,2,...
A batch arrival stream is a list [A_0, A_1, ...] where A_k carries the slots
delivering k events; a plain DMAP is the two-entry case.

MATLAB twins: dph_from_dist.m, dph_to_dmap.m, dmap_to_dph.m, dmap_is_renewal.m,
dmap_lambda.m, dmap_super.m, dmap_thin.m, dmap_compress.m,
dmap_compress_batch.m, mg1_dt_queue.m
"""

from typing import List, Optional, Sequence, Tuple

import numpy as np

from ...constants import ProcessType
from ..mc import dtmc_solve

__all__ = [
    'dph_from_dist', 'dph_to_dmap', 'dmap_to_dph', 'dmap_is_renewal',
    'dmap_moment', 'dmap_lambda', 'dmap_super', 'dmap_thin',
    'dmap_compress', 'dmap_compress_batch',
    'q_dt_map_map_1', 'q_dt_ph_ph_1', 'mg1_dt_queue',
]

_TOL = 1e-8


def dph_from_dist(proc_type, mean_slots: float, scv: float) -> Tuple[np.ndarray, np.ndarray]:
    """Exact discrete phase-type representation of a lattice-valued law.

    Geometric, Det and DiscreteUniform are represented EXACTLY, not
    moment-matched: a fitted surrogate would leave the lattice the caller relies
    on, so any other family raises instead.
    """
    if proc_type == ProcessType.GEOMETRIC:
        p = 1.0 / mean_slots
        if p > 1 + _TOL or p <= 0:
            raise ValueError(
                f"Geometric with mean {mean_slots} slots is outside the support {{1,2,...}}.")
        p = min(1.0, p)
        return np.array([[1.0]]), np.array([[1.0 - p]])

    if proc_type == ProcessType.DET:
        k = int(round(mean_slots))
        if abs(mean_slots - k) > _TOL * max(1.0, mean_slots) or k < 1:
            raise ValueError(
                f"Det of {mean_slots} slots is not a positive integral number of slots.")
        alpha = np.zeros((1, k))
        alpha[0, 0] = 1.0
        A = np.zeros((k, k))
        for i in range(k - 1):
            A[i, i + 1] = 1.0
        return alpha, A

    if proc_type == ProcessType.DUNIFORM:
        var_slots = scv * mean_slots ** 2
        width = np.sqrt(max(0.0, 12 * var_slots + 1)) - 1
        lo = int(round(mean_slots - width / 2))
        hi = int(round(mean_slots + width / 2))
        if lo < 1 or hi < lo:
            raise ValueError(
                f"DiscreteUniform spanning [{lo},{hi}] slots is outside the support {{1,2,...}}.")
        alpha = np.zeros((1, hi))
        alpha[0, 0] = 1.0
        A = np.zeros((hi, hi))
        for j in range(1, hi):
            # hazard of absorbing at step j, zero below the lower bound
            h = 0.0 if j < lo else 1.0 / (hi - j + 1)
            A[j - 1, j] = 1.0 - h
        return alpha, A

    raise ValueError(
        f"ProcessType {proc_type} has no exact discrete phase-type representation. "
        "The discrete-time path accepts Geometric, Det on the slot lattice, "
        "DiscreteUniform and DMAP.")


def dph_to_dmap(alpha: np.ndarray, A: np.ndarray) -> List[np.ndarray]:
    """Renewal DMAP [A, a @ alpha] of the DPH (alpha, A)."""
    A = np.asarray(A, dtype=float)
    alpha = np.asarray(alpha, dtype=float).reshape(1, -1)
    m = A.shape[0]
    a = np.ones((m, 1)) - A @ np.ones((m, 1))
    return [A.copy(), a @ alpha]


def dmap_is_renewal(D0: np.ndarray, D1: np.ndarray) -> bool:
    """True when D1 has rank one, i.e. the process renews at every event."""
    D1 = np.asarray(D1, dtype=float)
    if D1.shape[0] == 1:
        return True
    row_mass = D1.sum(axis=1)
    pivot = int(np.argmax(row_mass))
    if row_mass[pivot] <= 1e-14:
        return False
    rebuilt = np.outer(row_mass, D1[pivot, :] / row_mass[pivot])
    return bool(np.allclose(D1, rebuilt, atol=1e-8))


def dmap_to_dph(D0: np.ndarray, D1: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Discrete phase-type law underlying a renewal DMAP."""
    if not dmap_is_renewal(D0, D1):
        raise ValueError(
            "The DMAP does not renew at events, so it has no discrete phase-type form.")
    D1 = np.asarray(D1, dtype=float)
    row_mass = D1.sum(axis=1)
    pivot = int(np.argmax(row_mass))
    if row_mass[pivot] <= 1e-14:
        raise ValueError("The DMAP has no events: D1 is the zero matrix.")
    alpha = (D1[pivot, :] / row_mass[pivot]).reshape(1, -1)
    return alpha, np.asarray(D0, dtype=float)


def dmap_moment(D0: np.ndarray, D1: np.ndarray, orders: Sequence[int]) -> List[float]:
    """Raw moments of the interevent time of a DMAP, orders 1 to 3."""
    D0 = np.asarray(D0, dtype=float)
    D1 = np.asarray(D1, dtype=float)
    n = D0.shape[0]
    e = np.ones((n, 1))
    P = np.linalg.solve(np.eye(n) - D0, D1)
    al = np.asarray(dtmc_solve(P)).reshape(1, -1)
    A = np.linalg.inv(np.eye(n) - D0)
    Ae = A @ e
    m1 = float(al @ Ae)
    out = []
    for i in orders:
        if i == 1:
            out.append(m1)
        elif i == 2:
            out.append(float(2 * al @ A @ Ae) - m1)
        elif i == 3:
            out.append(float(6 * al @ A @ A @ Ae) - float(6 * al @ A @ Ae) + m1)
        else:
            raise ValueError("dmap_moment: raw moments of order > 3 not implemented")
    return out


def dmap_lambda(A: Sequence[np.ndarray]) -> float:
    """Mean number of EVENTS per slot, pi @ sum_k k*A_k @ e.

    A slot carrying a batch of two counts twice, which is what Little's law
    consumes downstream.
    """
    m = A[0].shape[0]
    P = np.zeros((m, m))
    for Ak in A:
        P = P + Ak
    pi = np.asarray(dtmc_solve(P)).reshape(1, -1)
    W = np.zeros((m, m))
    for k in range(1, len(A)):
        W = W + k * A[k]
    return float(np.asarray(pi @ W @ np.ones((m, 1))).reshape(-1)[0])


def dmap_super(A: Sequence[np.ndarray], B: Sequence[np.ndarray]) -> List[np.ndarray]:
    """Superposition, E_k = sum_{i+j=k} kron(A_i, B_j).

    NOT closed on DMAPs: two slotted streams fire in the same slot with positive
    probability, so the merged stream carries batches. Folding E_2 into E_1 would
    conserve neither the arrival rate nor the slot in which the work appears, so
    the batch dimension is kept and the station is solved as an M/G/1-type chain
    rather than a QBD.
    """
    p = len(A) - 1
    q = len(B) - 1
    E = []
    for k in range(p + q + 1):
        Ek = None
        for i in range(max(0, k - q), min(p, k) + 1):
            term = np.kron(A[i], B[k - i])
            Ek = term if Ek is None else Ek + term
        E.append(Ek)
    return E


def dmap_thin(A: Sequence[np.ndarray], p: float) -> List[np.ndarray]:
    """Bernoulli thinning, B_k = sum_{n>=k} C(n,k) p^k (1-p)^(n-k) A_n.

    The phase process is untouched, so this is exact for PROB/RAND routing.
    """
    if p < 0 or p > 1:
        raise ValueError(f"The routing probability must lie in [0,1], got {p}")
    from math import comb
    n = len(A) - 1
    B = []
    for k in range(n + 1):
        Bk = np.zeros_like(A[0], dtype=float)
        for j in range(k, n + 1):
            w = comb(j, k) * (p ** k) * ((1 - p) ** (j - k))
            if w > 0:
                Bk = Bk + w * A[j]
        B.append(Bk)
    # trailing zero batch levels carry no mass and only inflate the blocks
    while len(B) > 2 and np.max(np.abs(B[-1])) < 1e-14:
        B.pop()
    return B


def dmap_compress(DMAP: Sequence[np.ndarray], max_order: int) -> List[np.ndarray]:
    """Order reduction of a DMAP by matching three interevent moments.

    Correlation is NOT preserved, mirroring the continuous-time
    'mixture.order1' compression, and that is where the multi-station
    discrete-time path becomes approximate. Outside the DPH(2) region the
    fallback keeps the exact mean with a Geometric, so the arrival rate is
    conserved in every branch.
    """
    D0, D1 = DMAP[0], DMAP[1]
    if D0.shape[0] <= max_order:
        return [D0, D1]

    moms = dmap_moment(D0, D1, [1, 2, 3])
    try:
        from ...lib.thirdparty.butools.dph.canonical import DPH2From3Moments
        alpha, A = DPH2From3Moments(moms)
        cand = dph_to_dmap(np.asarray(alpha), np.asarray(A))
        if _is_feasible(cand):
            return cand
    except Exception:
        pass

    p = min(1.0, max(1e-14, 1.0 / moms[0]))
    return dph_to_dmap(np.array([[1.0]]), np.array([[1.0 - p]]))


def dmap_compress_batch(B: Sequence[np.ndarray], max_order: int) -> List[np.ndarray]:
    """Order reduction of a batch stream.

    Keeps the two features the downstream M/G/1-type solve consumes: the law of
    the time between NONEMPTY slots, matched to three moments, and the stationary
    batch-size distribution conditional on a nonempty slot, kept exactly. The
    event rate of the reduced stream equals the original one by construction.
    """
    if B[0].shape[0] <= max_order:
        return list(B)

    nb = len(B) - 1
    m = B[0].shape[0]
    Ptot = np.zeros((m, m))
    for Bk in B:
        Ptot = Ptot + Bk
    pi_phase = np.asarray(dtmc_solve(Ptot)).reshape(1, -1)
    e = np.ones((m, 1))
    qraw = np.array([float(pi_phase @ B[k] @ e) for k in range(1, nb + 1)])
    mass = float(qraw.sum())
    if mass <= 1e-14:
        raise ValueError("The batch stream carries no events, so it cannot be compressed.")

    marked = [B[0].copy(), Ptot - B[0]]
    marked_c = dmap_compress(marked, max_order)

    out = [marked_c[0]]
    for k in range(nb):
        out.append(marked_c[1] * (qraw[k] / mass))
    return out


def _is_feasible(DMAP: Sequence[np.ndarray]) -> bool:
    D0, D1 = np.asarray(DMAP[0]), np.asarray(DMAP[1])
    if D0.min() < -1e-10 or D1.min() < -1e-10:
        return False
    return bool(np.allclose((D0 + D1).sum(axis=1), 1.0, atol=1e-6))


def _qbd_r(Am1: np.ndarray, A0: np.ndarray, A1: np.ndarray) -> np.ndarray:
    """R matrix of a discrete-time QBD via BuTools cyclic reduction."""
    from ...lib.thirdparty.butools.mam.qbd import QBDFundamentalMatrices
    R = QBDFundamentalMatrices(np.matrix(Am1), np.matrix(A0), np.matrix(A1), matrices="R")
    return np.asarray(R, dtype=float)


def q_dt_map_map_1(C0, C1, D0, D1, max_num_comp: int = 1000) -> np.ndarray:
    """Queue length distribution of a discrete-time DMAP/DMAP/1/FCFS queue.

    Port of Q_DT_MAP_MAP_1.m of the Q-MAM library. The QBD blocks state the
    LAS-DA convention directly: A_1 = kron(C1,D0) says a job arriving at the end
    of a slot cannot be served within it, A_0 = kron(C0,D0) + kron(C1,D1) that a
    simultaneous arrival and completion leave the level unchanged, and the
    boundary B_k = kron(C_k, I) freezes the service phase while the system is
    empty.

    Only the queue length is returned. The waiting and sojourn pmfs the MATLAB
    routine also computes are deliberately not ported: no LINE caller consumes
    them, and the discrete-time solver path reads the queue length alone.
    """
    C0 = np.asarray(C0, dtype=float)
    C1 = np.asarray(C1, dtype=float)
    D0 = np.asarray(D0, dtype=float)
    D1 = np.asarray(D1, dtype=float)
    ma, ms = C0.shape[0], D0.shape[0]
    mtot = ma * ms

    pi_a = np.asarray(dtmc_solve(C0 + C1)).reshape(1, -1)
    avga = float(np.asarray(pi_a @ C1 @ np.ones((ma, 1))).reshape(-1)[0])
    pi_s = np.asarray(dtmc_solve(D0 + D1)).reshape(1, -1)
    avgs = float(np.asarray(pi_s @ D1 @ np.ones((ms, 1))).reshape(-1)[0])
    rho = avga / avgs
    if rho >= 1:
        raise ValueError(f"The load {rho} of the system exceeds one")

    Ims = np.eye(ms)
    Am1 = np.kron(C0, D1)
    A0 = np.kron(C0, D0) + np.kron(C1, D1)
    A1 = np.kron(C1, D0)
    Bm1 = np.kron(C0, D1)
    B0 = np.kron(C0, Ims)
    B1 = np.kron(C1, Ims)

    R = _qbd_r(Am1, A0, A1)

    # General boundary [B1; A0 + R*Am1]: the empty system has its own local
    # block, so levels 0 and 1 are solved together (QBD_pi.m else-branch)
    boundary = np.vstack([B1, A0 + R @ Am1])
    # rows: [level 0 ; level 1], columns: [to level 0, to level 1]
    joint = np.hstack([np.vstack([B0, Bm1]), boundary])
    pi01 = np.asarray(dtmc_solve(joint)).reshape(1, -1)

    temp = np.linalg.inv(np.eye(mtot) - R)
    pi0 = pi01[:, :mtot]
    pi1 = pi01[:, mtot:]
    normalizer = float(pi0.sum() + np.asarray(pi1 @ temp @ np.ones((mtot, 1))).reshape(-1)[0])
    pi0 = pi0 / normalizer
    pi1 = pi1 / normalizer

    levels = [pi0, pi1]
    summ = float(pi0.sum() + pi1.sum())
    cur = pi1
    it = 1
    while summ < 1 - 1e-10 and it < max_num_comp:
        cur = cur @ R
        levels.append(cur)
        summ += float(cur.sum())
        it += 1

    ql = np.array([float(lv.sum()) for lv in levels])
    total = ql.sum()
    if total > 0:
        ql = ql / total
    return ql


def q_dt_ph_ph_1(alpha, T, beta, S, max_num_comp: int = 1000) -> np.ndarray:
    """Queue length distribution of a discrete-time DPH/DPH/1/FCFS queue.

    A DPH renews at every event, so the streams are the DMAPs (T, t@alpha) and
    (S, s@beta) and the queue is the DMAP/DMAP/1 one. Routing through
    q_dt_map_map_1 keeps one statement of the LAS-DA convention in the codebase.
    """
    arv = dph_to_dmap(np.asarray(alpha), np.asarray(T))
    svc = dph_to_dmap(np.asarray(beta), np.asarray(S))
    return q_dt_map_map_1(arv[0], arv[1], svc[0], svc[1], max_num_comp)


def mg1_dt_queue(ARV: Sequence[np.ndarray], SVC: Sequence[np.ndarray],
                 max_num_comp: int = 500, want_departure: bool = False):
    """Discrete-time single-server queue with batch DMAP arrivals, DBMAP/DMAP/1.

    The chain is M/G/1-type because a slot may deliver a batch: with arrival
    matrices A_k and service pair (S0,S1),
    A^(-1) = kron(A_0,S1), A^(k) = kron(A_k,S0) + kron(A_{k+1},S1),
    B^(k) = kron(A_k,I), the boundary row holding the empty system where no
    service runs.

    Returns (QN, UN, TN, ql, dep) with dep the departure DMAP of the
    level-truncated chain, or None when not requested.
    """
    from ...lib.thirdparty.butools.mam.mg1gm1 import MG1FundamentalMatrix, MG1StationaryDistr

    S0 = np.asarray(SVC[0], dtype=float)
    S1 = np.asarray(SVC[1], dtype=float)
    ms = S0.shape[0]
    ma = ARV[0].shape[0]
    K = len(ARV) - 1
    Ims = np.eye(ms)

    lam = dmap_lambda(ARV)
    mu = dmap_lambda([S0, S1])
    if lam >= mu:
        raise ValueError(
            f"The discrete-time load {lam / mu} of the station is not below one "
            f"({lam} arrivals per slot against {mu} completions per busy slot).")

    Ablocks = [np.kron(ARV[0], S1)]
    for k in range(K + 1):
        blk = np.kron(ARV[k], S0)
        if k + 1 <= K:
            blk = blk + np.kron(ARV[k + 1], S1)
        Ablocks.append(blk)
    Bblocks = [np.kron(ARV[k], Ims) for k in range(K + 1)]
    # the boundary row must be as long as the repeating one for MG1StationaryDistr
    while len(Bblocks) < len(Ablocks):
        Bblocks.append(np.zeros_like(Bblocks[0]))

    Amat = [np.matrix(b) for b in Ablocks]
    Bmat = [np.matrix(b) for b in Bblocks]
    G = MG1FundamentalMatrix(Amat)
    pivec = np.asarray(MG1StationaryDistr(Amat, Bmat, G, max_num_comp)).reshape(1, -1)

    m = ma * ms
    nlev = pivec.shape[1] // m
    ql = np.array([float(pivec[0, i * m:(i + 1) * m].sum()) for i in range(nlev)])
    total = ql.sum()
    if total > 0:
        ql = ql / total

    QN = float(sum(i * ql[i] for i in range(nlev)))
    UN = 1.0 - float(ql[0])

    dep = None
    if want_departure:
        dep = _departure_process(ARV, S0, S1, ql, m, ms, K, Ims)
    return QN, UN, lam, ql, dep


def _departure_process(ARV, S0, S1, ql, m, ms, K, Ims) -> List[np.ndarray]:
    """Departure DMAP of the level-truncated chain.

    Levels 0..L with arrivals that would cross L held at L. The truncation is a
    level cut, not a rate change, so the departure rate stays within the mass
    left above L.
    """
    cum = np.cumsum(ql)
    idx = np.nonzero(cum > 1 - 1e-10)[0]
    L = int(idx[0]) if idx.size else len(ql) - 1
    L = max(1, L)

    nstates = (L + 1) * m
    D0 = np.zeros((nstates, nstates))
    D1 = np.zeros((nstates, nstates))
    for i in range(L + 1):
        rows = slice(i * m, (i + 1) * m)
        for k in range(K + 1):
            if i == 0:
                tgt = min(L, k)
                D0[rows, tgt * m:(tgt + 1) * m] += np.kron(ARV[k], Ims)
            else:
                tno = min(L, i + k)
                D0[rows, tno * m:(tno + 1) * m] += np.kron(ARV[k], S0)
                tdep = min(L, i - 1 + k)
                D1[rows, tdep * m:(tdep + 1) * m] += np.kron(ARV[k], S1)
    return [D0, D1]
