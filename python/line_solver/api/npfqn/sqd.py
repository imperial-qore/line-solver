"""
NPFQN Smith Queue Decomposition (SQD) approximate MVA for Blocking-After-Service (BAS).

Solves a finite-buffer closed queueing network under Blocking-After-Service
(manufacturing/transfer) blocking directly from its NetworkStruct. The method is an
AMVA-style population recursion in which each finite-capacity station is described by
a load-dependent effective service rate calibrated from an M/M/1/K blocking
probability; downstream blocking is propagated through the effective routing between
service stations. Delay (INF/EXT) stations are treated as infinite-capacity pure-delay
nodes. Single-chain (chain-aggregated) demands only.

Originally contributed as SolverDBT by Avinash Bommareddy (Imperial College London
FYP, 2026); refactored here into an sn-based API function.
"""

import numpy as np
from dataclasses import dataclass
from typing import Optional

from ..sn import NetworkStruct, SchedStrategy, sn_get_demands_chain

_INITIAL_V1 = 692.192
_CALIBRATION_EPSILON = 0.05


@dataclass
class NpfqnSqdResult:
    """Per-station results of the BAS approximation (single chain).

    Attributes:
        X: (M,) per-station throughput
        Q: (M,) per-station queue length
        U: (M,) per-station utilization
        R: (M,) per-station residence time
    """
    X: np.ndarray
    Q: np.ndarray
    U: np.ndarray
    R: np.ndarray


def _is_delay(sched_value) -> bool:
    name = getattr(sched_value, 'name', str(sched_value))
    return name in ('INF', 'EXT')


def _mm1k_blocking(K: int, rho: float) -> float:
    """Steady-state blocking probability of an M/M/1/K queue at load rho."""
    if rho <= 1e-15:
        return 0.0
    if abs(rho - 1.0) < 1e-9:
        return 1.0 / (K + 1.0)
    return (1.0 - rho) * rho ** K / (1.0 - rho ** (K + 1))


def _compute_beta_gamma(K: int, p_block_down: float, mode: int):
    """Calibrate the (beta, gamma) shape parameters of the load-dependent effective
    service rate from the station capacity K and a downstream blocking probability."""
    if K <= 2 or K == np.iinfo(np.int64).max:
        return K, 1.0

    a = 2.0
    b = float(K)

    if mode == 2:        # blocking-aware
        pK = max(1e-6, min(p_block_down, 1.0 - 1e-6))
        Va = 1.0 - pK * (a - 1.0) / (b - 1.0)
        Vb = 1.0 - pK
    elif mode == 1:      # fixed heuristic
        Va = (b - a) / b
        Vb = _CALIBRATION_EPSILON
    else:                # mode 0: base
        return K, 1.0

    Va = min(Va, 0.999)
    Va = max(Va, 0.01)
    Vb = max(Vb, 1e-6)
    if Vb >= Va:
        return K, 1.0

    ln_va = np.log(Va)
    ln_vb = np.log(Vb)
    gamma = np.log(ln_va / ln_vb) / np.log((a - 1.0) / (b - 1.0))
    gamma = max(0.5, min(gamma, 10.0))
    beta = (a - 1.0) / (-ln_va) ** (1.0 / gamma)

    if not np.isfinite(gamma) or not np.isfinite(beta) or beta <= 0:
        return K, 1.0

    return beta, gamma


def _compute_effective_routing(sn: NetworkStruct, M: int, R: int, is_delay) -> np.ndarray:
    """Build station-to-station effective routing among service stations, collapsing
    pass-through delay (INF/EXT) stations into a single hop."""
    s2s = np.asarray(sn.stationToStateful).flatten()
    p = np.zeros((M, M))
    for i in range(M):
        if is_delay[i]:
            continue
        sf_i = int(s2s[i])
        for j in range(M):
            sf_j = int(s2s[j])
            p_ij = sn.rt[sf_i * R, sf_j * R]
            if p_ij <= 0:
                continue
            if not is_delay[j]:
                p[i, j] += p_ij
            else:
                for k in range(M):
                    if not is_delay[k]:
                        sf_k = int(s2s[k])
                        p_jk = sn.rt[sf_j * R, sf_k * R]
                        if p_jk > 0:
                            p[i, k] += p_ij * p_jk
    return p


def npfqn_sqd(
    sn: NetworkStruct,
    N: Optional[int] = None,
    calibration_mode: int = 0,
    server_blocking_time: bool = True,
    neighbor_mode: str = 'downstream',
    v1_policy: str = 'compound',
    initial_v1: Optional[np.ndarray] = None,
) -> NpfqnSqdResult:
    """Solve a BAS closed network from its NetworkStruct.

    Args:
        sn: chain-aggregated network structure
        N: total closed-class population (defaults to sn.nclosedjobs)
        calibration_mode: 0=base, 1=fixed heuristic, 2=blocking-aware
        server_blocking_time: add a manufacturing-blocking term to server time
        neighbor_mode: 'downstream' (routed) or 'ownserver' blocking aggregation
        v1_policy: 'compound' or 'fresh' load-dependent rate-scale update
        initial_v1: optional per-station initial V1 (None = _INITIAL_V1)

    Returns:
        NpfqnSqdResult with per-station X, Q, U, R.
    """
    if N is None:
        N = sn.nclosedjobs
    M = sn.nstations
    R = sn.nclasses
    INT_MAX = np.iinfo(np.int64).max

    demands = sn_get_demands_chain(sn)
    Vchain = demands.Vchain
    STchain = demands.STchain

    V = np.array([Vchain[i, 0] for i in range(M)], dtype=float)
    ST = np.array([STchain[i, 0] for i in range(M)], dtype=float)

    cap_arr = np.asarray(sn.cap).flatten()
    sched_dict = sn.sched if sn.sched else {}
    is_delay = np.zeros(M, dtype=bool)
    cap = np.zeros(M, dtype=np.int64)
    for i in range(M):
        is_delay[i] = _is_delay(sched_dict.get(i))
        if is_delay[i]:
            cap[i] = INT_MAX
        else:
            c = cap_arr[i]
            cap[i] = INT_MAX if (not np.isfinite(c) or c > 1e14) else int(c)

    p_eff = _compute_effective_routing(sn, M, R, is_delay)

    # Exact job-hole duality for a closed two-station cyclic BAS network. Its
    # throughput is symmetric about C/2 with C = K1 + K2 + 2 (each station holds up
    # to K_i jobs plus one blocked-from-upstream phantom slot), and the population-N
    # solution maps onto the accurate low-population dual at C-N:
    #   X_i(N)  = X_i(C-N)
    #   Q_i(N)  = (K_i + 1) - Q_other(C-N)     (reverse cycle = swap the two stations)
    #   U_i(N)  = X_i(N) * ST_i,   W_i(N) = Q_i(N) / X_i(N)
    # The AMVA recursion below is accurate up to the peak (N <= C/2) but over-predicts
    # past it because it does not capture the backward blocking that reduces throughput.
    # Validated against JMT to three digits (see FiniteBufferBASTest Tables 1-3).
    if (M == 2 and not is_delay[0] and not is_delay[1]
            and cap[0] < INT_MAX and cap[1] < INT_MAX
            and p_eff[0, 1] > 1.0 - 1e-9 and p_eff[1, 0] > 1.0 - 1e-9):
        C = int(cap[0]) + int(cap[1]) + 2
        Ndual = C - N
        if 2 * N > C and Ndual >= 1 and Ndual < N:
            dual = npfqn_sqd(sn, Ndual, calibration_mode, server_blocking_time,
                             neighbor_mode, v1_policy, initial_v1)
            XN = np.zeros(M)
            QN = np.zeros(M)
            UN = np.zeros(M)
            WN = np.zeros(M)
            for i in range(M):
                T_i = dual.X[i]
                Q_i = (int(cap[i]) + 1) - dual.Q[1 - i]
                U_i = min(1.0, T_i * ST[i])
                W_i = Q_i / T_i if T_i > 1e-15 else 0.0
                XN[i] = T_i
                QN[i] = Q_i
                UN[i] = U_i
                WN[i] = W_i
            return NpfqnSqdResult(X=XN, Q=QN, U=UN, R=WN)

    V1 = np.zeros(M)
    v1init = np.zeros(M)
    L_buf = np.zeros(M)
    L_svr = np.zeros(M)
    for i in range(M):
        v1init[i] = _INITIAL_V1 if initial_v1 is None else initial_v1[i]
        V1[i] = v1init[i]

    X = 0.0
    W_buf = np.zeros(M)
    W_svr = np.zeros(M)

    ownserver = (neighbor_mode == 'ownserver')
    fresh = (v1_policy == 'fresh')

    for pop in range(1, N + 1):

        # wait times
        for i in range(M):
            if is_delay[i]:
                W_buf[i] = 0.0
                W_svr[i] = ST[i]
                continue

            if ownserver:
                p_block_down = _mm1k_blocking(cap[i], X * V[i] * ST[i]) if cap[i] < INT_MAX else 0.0
            else:
                p_block_down = 0.0
                for j in range(M):
                    if (not is_delay[j]) and p_eff[i, j] > 0 and cap[j] < INT_MAX:
                        rho_j = X * V[j] * ST[j]
                        p_block_down += p_eff[i, j] * _mm1k_blocking(cap[j], rho_j)

            if ownserver:
                n = L_svr[i]
            else:
                n = 0.0
                for j in range(M):
                    n += p_eff[i, j] * L_svr[j]

            if n < 1e-10:
                mu_n = V1[i]
            else:
                beta, gamma = _compute_beta_gamma(int(cap[i]), p_block_down, calibration_mode)
                base = max(0.0, (n - 1.0) / beta)
                exp_arg = base ** gamma
                mu_n = n * V1[i] * np.exp(-exp_arg)   # Eq.13
            if mu_n < 1e-10:
                mu_n = 1e-10

            W_buf[i] = (1.0 / mu_n) * (1.0 + n)        # Eq.18
            W_svr[i] = ST[i] * (1.0 + L_svr[i])        # Eq.17

            if server_blocking_time:
                bt = 0.0
                for j in range(M):
                    if (not is_delay[j]) and p_eff[i, j] > 0 and cap[j] < INT_MAX:
                        rho_j = X * V[j] * ST[j]
                        p_bj = _mm1k_blocking(cap[j], rho_j)
                        denom = ST[i] + ST[j]
                        theta = ST[j] / denom if denom > 1e-15 else 0.0
                        bt += p_eff[i, j] * p_bj * ST[j] * theta
                W_svr[i] += bt

        # throughput
        sum_vw = float(np.sum(V * (W_buf + W_svr)))
        X = pop / sum_vw if sum_vw > 1e-15 else 0.0

        # queue lengths
        L_buf = X * V * W_buf
        L_svr = X * V * W_svr

        # adjust V1
        if pop < N:
            for i in range(M):
                if is_delay[i]:
                    continue
                if ownserver:
                    p_block = _mm1k_blocking(cap[i], X * V[i] * ST[i]) if cap[i] < INT_MAX else 0.0
                else:
                    p_block = 0.0
                    for j in range(M):
                        if (not is_delay[j]) and p_eff[i, j] > 0 and cap[j] < INT_MAX:
                            rho_j = X * V[j] * ST[j]
                            p_block += p_eff[i, j] * _mm1k_blocking(cap[j], rho_j)
                V1[i] = v1init[i] * (1.0 - p_block) if fresh else V1[i] * (1.0 - p_block)
                if V1[i] < 1e-10:
                    V1[i] = 1e-10

    XN = np.zeros(M)
    QN = np.zeros(M)
    UN = np.zeros(M)
    WN = np.zeros(M)
    for i in range(M):
        T_i = X * V[i]
        Q_i = L_buf[i] + L_svr[i]
        W_tot = W_buf[i] + W_svr[i]
        R_i = Q_i / T_i if T_i > 1e-15 else W_tot
        U_i = min(1.0, T_i * ST[i])
        XN[i] = T_i
        QN[i] = Q_i
        UN[i] = U_i
        WN[i] = R_i

    return NpfqnSqdResult(X=XN, Q=QN, U=UN, R=WN)


__all__ = ['npfqn_sqd', 'NpfqnSqdResult']
