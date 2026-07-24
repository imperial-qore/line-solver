"""
Asymptotic Bounds for Product-Form Queueing Networks.

Implements various bounds on throughput and queue lengths for closed
queueing networks, including:
- Balanced Job Bounds (Zahorjan)
- Asymptotic Bounds (Zahorjan-Gittelsohn-Bryant)
- ZGSB Bounds (Zahorjan-Gittelsohn-Schweitzer-Bryant)
"""

import numpy as np
from typing import Union

def pfqn_xzabalow(
    L: np.ndarray,
    N: Union[int, float],
    Z: float
) -> float:
    """
    Lower ABA (asymptotic bound analysis) bound on throughput.

    Returns N / (Z + sum(L)*N), the ABA lower throughput bound for
    single-class closed queueing networks. This is NOT the classical
    Zahorjan-Balanced (balanced job bounds) lower bound; that one is
    pfqn_xzgsblow, which is tighter.

    Args:
        L: Service demand vector (M,).
        N: Population (scalar).
        Z: Think time.

    Returns:
        Lower bound on throughput.
    """
    L = np.asarray(L, dtype=float).ravel()
    M = len(L)

    # JIT dispatch for large problems

    Ltot = np.sum(L)
    return float(N / (Z + Ltot * N))

def pfqn_xzabaup(
    L: np.ndarray,
    N: Union[int, float],
    Z: float
) -> float:
    """
    Upper asymptotic bound on throughput (Zahorjan-Balanced).

    Provides a simple upper bound on system throughput for single-class
    closed queueing networks based on bottleneck analysis.

    Args:
        L: Service demand vector (M,).
        N: Population (scalar).
        Z: Think time.

    Returns:
        Upper bound on throughput.
    """
    L = np.asarray(L, dtype=float).ravel()
    M = len(L)

    # JIT dispatch for large problems

    return float(min(1.0 / np.max(L), N / (np.sum(L) + Z)))

def pfqn_qzgblow(
    L: np.ndarray,
    N: Union[int, float],
    Z: float,
    i: int
) -> float:
    """
    Lower asymptotic bound on queue length (Zahorjan-Gittelsohn-Bryant).

    Args:
        L: Service demand vector (M,).
        N: Population (scalar).
        Z: Think time.
        i: Station index (0-based).

    Returns:
        Lower bound on mean queue length at station i.
    """
    L = np.asarray(L, dtype=float).ravel()
    M = len(L)

    # JIT dispatch for large problems

    yi = N * L[i] / (Z + np.sum(L) + np.max(L) * N)

    if yi >= 1:
        return float(N)

    Qgb = yi / (1 - yi) - (yi ** (N + 1)) / (1 - yi)
    return float(max(0, Qgb))

def pfqn_qzgbup(
    L: np.ndarray,
    N: Union[int, float],
    Z: float,
    i: int
) -> float:
    """
    Upper asymptotic bound on queue length (Zahorjan-Gittelsohn-Bryant).

    Args:
        L: Service demand vector (M,).
        N: Population (scalar).
        Z: Think time.
        i: Station index (0-based).

    Returns:
        Upper bound on mean queue length at station i.
    """
    L = np.asarray(L, dtype=float).ravel()
    M = len(L)

    # JIT dispatch for large problems

    sigma = np.sum(L ** 2) / np.sum(L)

    # Compute upper bound on throughput at N-1
    if N > 1:
        X_N_minus_1 = pfqn_xzabaup(L, N - 1, Z)
    else:
        X_N_minus_1 = 0

    Yi = L[i] * min(
        1.0 / np.max(L),
        N / (Z + np.sum(L) + sigma * (N - 1 - Z * X_N_minus_1))
    )

    if Yi < 1:
        Qgb = Yi / (1 - Yi) - (Yi ** (N + 1)) / (1 - Yi)
        return float(max(0, Qgb))
    else:
        return float(N)

def pfqn_xzgsblow(
    L: np.ndarray,
    N: Union[int, float],
    Z: float
) -> float:
    """
    Lower asymptotic bound on throughput (Zahorjan-Gittelsohn-Schweitzer-Bryant).

    Provides a tighter lower bound than pfqn_xzabalow by accounting for
    queue length bounds.

    Args:
        L: Service demand vector (M,).
        N: Population (scalar).
        Z: Think time.

    Returns:
        Lower bound on throughput.
    """
    L = np.asarray(L, dtype=float).ravel()
    M = len(L)

    # JIT dispatch for large problems

    max_L = np.max(L)

    R = Z + np.sum(L) + max_L * (N - 1)

    for i in range(M):
        if L[i] < max_L:
            R = R + (L[i] - max_L) * pfqn_qzgblow(L, N - 1, Z, i)

    discriminant = R ** 2 - 4 * Z * max_L * (N - 1)
    if discriminant < 0:
        # Fall back to simple bound
        return pfqn_xzabalow(L, N, Z)

    X = 2 * N / (R + np.sqrt(discriminant))
    return float(X)

def pfqn_xzgsbup(
    L: np.ndarray,
    N: Union[int, float],
    Z: float
) -> float:
    """
    Upper asymptotic bound on throughput (Zahorjan-Gittelsohn-Schweitzer-Bryant).

    Provides a tighter upper bound than pfqn_xzabaup by accounting for
    queue length bounds.

    Args:
        L: Service demand vector (M,).
        N: Population (scalar).
        Z: Think time.

    Returns:
        Upper bound on throughput.
    """
    L = np.asarray(L, dtype=float).ravel()
    M = len(L)

    # JIT dispatch for large problems

    max_L = np.max(L)

    R = Z + np.sum(L) + max_L * (N - 1)

    for i in range(M):
        if L[i] < max_L:
            R = R + (L[i] - max_L) * pfqn_qzgbup(L, N - 1, Z, i)

    discriminant = R ** 2 - 4 * Z * max_L * N
    if discriminant < 0:
        # Fall back to simple bound
        return pfqn_xzabaup(L, N, Z)

    X = 2 * N / (R + np.sqrt(discriminant))
    return float(X)

# Majumdar-Woodside discipline codes for pfqn_mwrbb `sched`.
MWRBB_FIFO = 0
MWRBB_PS = 1
MWRBB_NPPRIO = 2   # non-preemptive priority
MWRBB_PPPRIO = 3   # preemptive priority
MWRBB_ABA = 4      # ABA full-contention (discipline-independent)


def _mwrbb_wrest(k, c, V, S, N, fup, fc, sched, prio, C):
    """Per-visit residence at station k for class c EXCLUDING the isolated
    higher-priority (1/f_c) term."""
    Vkc = V[k, c]
    Skc = S[k, c]
    d = sched[k]
    if d == MWRBB_FIFO:
        s = 0.0
        for m in range(C):
            pcm = 1.0 if fc * Vkc == 0 else min(1.0, (fup[m] * V[k, m]) / (fc * Vkc))
            s += N[m] * S[k, m] * pcm
        return s   # own service is the m=c term (= N_c S_kc)
    elif d == MWRBB_PS:
        dp = 0.0
        for m in range(C):
            Ncont = N[c] - 1 if m == c else N[m]
            term = Skc if fc * Vkc == 0 else min(Skc, (fup[m] * V[k, m] * S[k, m]) / (fc * Vkc))
            dp += Ncont * term
        return Skc + dp
    elif d == MWRBB_ABA:   # ABA full-contention (P_cm = 1)
        s = 0.0
        for m in range(C):
            s += N[m] * S[k, m]   # wait behind full service of all customers
        return s
    else:   # NPPRIO or PPPRIO
        dp = 0.0
        for m in range(C):
            if prio[m] == prio[c]:   # equal priority (includes c)
                Ncont = N[c] - 1 if m == c else N[m]
                pcm = 1.0 if fc * Vkc == 0 else min(1.0, (fup[m] * V[k, m]) / (fc * Vkc))
                dp += Ncont * S[k, m] * pcm
            # higher priority handled via Bh in the caller
        if d == MWRBB_NPPRIO:   # lower-priority water-filling (Lemma 4)
            L = [m for m in range(C) if prio[m] > prio[c]]
            L.sort(key=lambda m: S[k, m], reverse=True)
            budget = 1.0
            for l in L:
                capr = np.inf if fc * Vkc == 0 else (fup[l] * V[k, l]) / (fc * Vkc)
                al = 0.0 if N[l] <= 0 else min(budget / N[l], capr)
                if al < 0:
                    al = 0.0
                dp += N[l] * al * S[k, l]
                budget -= N[l] * al
                if budget < 0:
                    budget = 0.0
        return Skc + dp


def pfqn_mwrbb(
    V: np.ndarray,
    S: np.ndarray,
    N: np.ndarray,
    Z: np.ndarray = None,
    sched: np.ndarray = None,
    prio: np.ndarray = None,
):
    """
    Majumdar-Woodside robust box bounds on throughput for closed multiclass
    queueing networks with mixed scheduling disciplines.

    Computes distribution-insensitive (NBUE) upper and lower bounds on the
    per-class system throughput of a closed multiclass queueing network, per
    S. Majumdar and C.M. Woodside, "Robust bounds and throughput guarantees
    for closed multiclass queueing networks", Performance Evaluation 32 (1998)
    101-136. The upper bound intersects the no-contention bound (eq. 2) with
    the utilization-based bound (eq. 3) and is discipline-independent. The
    lower bound is the multiclass throughput guarantee of Theorem 2 (eq. 15):
    X_c >= N_c / (Z_c + sum_k V_kc (S_kc + d_kc+)), where d_kc+ depends on the
    discipline at station k -- FIFO (Theorem 1 / Lemma 1), processor sharing
    (Lemma 2), preemptive priority (Lemma 3), non-preemptive priority
    (Lemmas 4-5). The coupled inequalities are resolved by the interval-
    narrowing fixed point reproducing the BNR-Prolog robust box bounds; for a
    single FIFO class it reduces to the Muntz-Wong bounds. Only queueing
    stations are passed; Z aggregates the pure-delay stations.

    Args:
        V: (K, C) mean visits of class c at queueing station k.
        S: (K, C) mean service demand per visit of class c at station k.
        N: (C,) population of class c.
        Z: (C,) think time of class c (default zeros).
        sched: (K,) discipline code per station (0=FIFO, 1=PS,
            2=non-preemptive priority, 3=preemptive priority, 4=ABA
            full-contention discipline-independent); default all FIFO.
        prio: (C,) class priority, lower value = higher priority; default equal.

    Returns:
        Xlo: (C,) lower bound on class throughput (Theorem 2).
        Xup: (C,) upper bound on class throughput (eqs. 2-3).
        Wlo: (K, C) per-visit residence time consistent with the lower bound.
    """
    V = np.asarray(V, dtype=float)
    S = np.asarray(S, dtype=float)
    if V.ndim == 1:
        V = V.reshape(-1, 1)
        S = S.reshape(-1, 1)
    K, C = V.shape
    N = np.asarray(N, dtype=float).ravel()
    Z = np.zeros(C) if Z is None else np.asarray(Z, dtype=float).ravel()
    sched = np.zeros(K, dtype=int) if sched is None else np.asarray(sched).astype(int).ravel()
    prio = np.zeros(C) if prio is None else np.asarray(prio, dtype=float).ravel()

    # no-contention upper bound on the cycle rate f_c = X_c/N_c (eqs. 1-2)
    fup = np.zeros(C)
    flo = np.zeros(C)
    for c in range(C):
        fup[c] = 1.0 / (Z[c] + np.sum(V[:, c] * S[:, c]))

    maxiter = 20000
    tol = 1e-13
    for _ in range(maxiter):
        maxdelta = 0.0
        # utilization-based narrowing of the upper bounds (eq. 3)
        for c in range(C):
            cap = fup[c]
            for k in range(K):
                other = 0.0
                for m in range(C):
                    if m != c:
                        other += N[m] * V[k, m] * S[k, m] * flo[m]
                denomk = N[c] * V[k, c] * S[k, c]
                if denomk > 0:
                    cap = min(cap, (1.0 - other) / denomk)
            newfup = min(fup[c], max(cap, 0.0))
            maxdelta = max(maxdelta, abs(newfup - fup[c]))
            fup[c] = newfup
        # lower-bound narrowing (Theorem 2, eq. 15). Higher-priority delay
        # carries a 1/f_c factor, isolated as Bh: f_c = (1 - Bh) / DEN.
        for c in range(C):
            DEN = Z[c]
            Bh = 0.0
            fc = flo[c]
            for k in range(K):
                Vkc = V[k, c]
                if Vkc == 0:
                    continue
                if sched[k] in (MWRBB_PPPRIO, MWRBB_NPPRIO):
                    for m in range(C):
                        if prio[m] < prio[c]:
                            Bh += N[m] * fup[m] * V[k, m] * S[k, m]
                DEN += Vkc * _mwrbb_wrest(k, c, V, S, N, fup, fc, sched, prio, C)
            val = (1.0 - Bh) / DEN
            if val < 0:
                val = 0.0
            newflo = max(flo[c], val)
            maxdelta = max(maxdelta, abs(newflo - flo[c]))
            flo[c] = newflo
        if maxdelta < tol:
            break

    Xlo = N * flo
    Xup = N * fup
    Wlo = np.zeros((K, C))
    for c in range(C):
        fc = flo[c]
        for k in range(K):
            Vkc = V[k, c]
            if Vkc == 0:
                continue
            W = _mwrbb_wrest(k, c, V, S, N, fup, fc, sched, prio, C)
            if sched[k] in (MWRBB_NPPRIO, MWRBB_PPPRIO) and fc * Vkc > 0:
                for m in range(C):
                    if prio[m] < prio[c]:
                        W += N[m] * fup[m] * V[k, m] * S[k, m] / (fc * Vkc)
            Wlo[k, c] = W
    return Xlo, Xup, Wlo

__all__ = [
    'pfqn_xzabalow',
    'pfqn_xzabaup',
    'pfqn_qzgblow',
    'pfqn_qzgbup',
    'pfqn_xzgsblow',
    'pfqn_xzgsbup',
    'pfqn_mwrbb',
]
