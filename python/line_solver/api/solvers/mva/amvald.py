"""
Approximate MVA for Load-Dependent (AMVA-LD) networks.

This module implements the approximate MVA solver for load-dependent
queueing networks using QD-AMVA style corrections with under-relaxation.

References:
    Original MATLAB: matlab/src/solvers/MVA/solver_amvald.m
    Original MATLAB: matlab/src/solvers/MVA/solver_amvald_forward.m
"""

import numpy as np
from dataclasses import dataclass
from typing import Optional, Tuple, Dict
from enum import IntEnum

# Import the proper pfqn_lldfun from utils
from ...pfqn.utils import pfqn_lldfun


class SchedStrategy(IntEnum):
    """Scheduling strategies for queueing stations."""
    EXT = -1     # External arrival (source)
    INF = 0      # Infinite server (delay)
    FCFS = 1     # First-come first-served
    LCFS = 2     # Last-come first-served
    SIRO = 3     # Service in random order
    PS = 4       # Processor sharing
    LCFSPR = 5   # LCFS with preemptive resume
    HOL = 6      # Head of line (priority)
    DPS = 7      # Discriminatory processor sharing


@dataclass
class AmvaldOptions:
    """Options for AMVA-LD solver."""
    method: str = 'default'
    iter_tol: float = 1e-6
    iter_max: int = 1000
    tol: float = 1e-8
    verbose: bool = False

    @dataclass
    class Config:
        multiserver: str = 'default'
        np_priority: str = 'default'
        highvar: str = 'default'

    config: Config = None

    def __post_init__(self):
        if self.config is None:
            self.config = AmvaldOptions.Config()


@dataclass
class AmvaldResult:
    """Result from AMVA-LD solver."""
    Q: np.ndarray      # Queue lengths (M x K)
    U: np.ndarray      # Utilization (M x K)
    R: np.ndarray      # Response times (M x K)
    T: np.ndarray      # Throughputs (M x K)
    C: np.ndarray      # Cycle times (1 x K)
    X: np.ndarray      # System throughputs (1 x K)
    lG: float          # Log normalizing constant
    totiter: int       # Total iterations




def solver_amvald_forward(
    M: int, K: int,
    nservers: np.ndarray,
    lldscaling: Optional[np.ndarray],
    sched: Dict[int, int],
    gamma: np.ndarray,
    tau: np.ndarray,
    Qchain_in: np.ndarray,
    Xchain_in: np.ndarray,
    Uchain_in: np.ndarray,
    STchain_in: np.ndarray,
    Vchain_in: np.ndarray,
    Nchain_in: np.ndarray,
    options: AmvaldOptions
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Forward step of AMVA-LD: compute waiting times from current estimates.

    Args:
        M: Number of stations
        K: Number of chains
        nservers: Number of servers at each station (M,)
        lldscaling: Load-dependent scaling matrix (M x Nmax)
        sched: Scheduling strategy for each station
        gamma: Correction factors
        tau: Throughput differences
        Qchain_in: Current queue length estimates (M x K)
        Xchain_in: Current throughput estimates (K,)
        Uchain_in: Current utilization estimates (M x K)
        STchain_in: Service times (M x K)
        Vchain_in: Visit ratios (M x K)
        Nchain_in: Population per chain (K,)
        options: Solver options

    Returns:
        Wchain: Waiting times (M x K)
        STeff: Effective service times (M x K)
    """
    # Initialize gamma if empty
    if gamma.size == 0:
        gamma = np.zeros((M, K))

    Nt = np.sum(Nchain_in[np.isfinite(Nchain_in)])
    if np.all(np.isinf(Nchain_in)):
        delta = 1.0
    else:
        delta = (Nt - 1) / Nt if Nt > 0 else 0.0

    deltaclass = (Nchain_in - 1) / Nchain_in
    deltaclass[np.isinf(Nchain_in)] = 1.0
    deltaclass[Nchain_in == 0] = 0.0

    # Find class types
    ocl = np.where(np.isinf(Nchain_in))[0]  # open classes
    ccl = np.where(np.isfinite(Nchain_in) & (Nchain_in > 0))[0]  # closed classes
    nnzclasses = np.where(Nchain_in > 0)[0]  # non-zero classes

    # Compute arrival queue lengths
    interpTotArvlQlen = np.zeros(M)
    selfArvlQlenSeenByClosed = np.zeros((M, K))
    totArvlQlenSeenByClosed = np.zeros((M, K))
    stationaryQlen = np.zeros((M, K))

    for k in range(M):
        interpTotArvlQlen[k] = delta * np.sum(Qchain_in[k, nnzclasses])
        for r in nnzclasses:
            selfArvlQlenSeenByClosed[k, r] = deltaclass[r] * Qchain_in[k, r]
            totArvlQlenSeenByClosed[k, r] = (deltaclass[r] * Qchain_in[k, r] +
                                              np.sum(Qchain_in[k, nnzclasses]) - Qchain_in[k, r])
            stationaryQlen[k, r] = Qchain_in[k, r]

    # Initialize lldscaling if empty
    if lldscaling is None or lldscaling.size == 0:
        lldscaling = np.ones((M, max(1, int(np.ceil(Nt)))))

    # Compute LLD term (using proper cubic interpolation like MATLAB)
    lldterm = pfqn_lldfun(1 + interpTotArvlQlen, lldscaling)

    # Compute multi-server term
    msterm = pfqn_lldfun(1 + interpTotArvlQlen, None, nservers)

    # For FCFS-type stations, use Seidmann approximation
    for k in range(M):
        sched_k = sched.get(k, SchedStrategy.FCFS)
        if sched_k in [SchedStrategy.FCFS, SchedStrategy.SIRO, SchedStrategy.LCFSPR]:
            if nservers[k] > 0 and not np.isinf(nservers[k]):
                msterm[k] = 1.0 / nservers[k]

    # Compute effective service times
    Wchain = np.zeros((M, K))
    STeff = np.zeros((M, K))

    lldterm_2d = np.tile(lldterm.reshape(-1, 1), (1, K))

    for r in nnzclasses:
        for k in range(M):
            STeff[k, r] = STchain_in[k, r] * lldterm_2d[k, r] * msterm[k]

    # Compute waiting times based on scheduling strategy
    for ir, r in enumerate(nnzclasses):
        sd = np.setdiff1d(nnzclasses, [r])  # other classes

        for k in range(M):
            sched_k = sched.get(k, SchedStrategy.FCFS)

            if sched_k == SchedStrategy.INF:
                # Infinite server (delay)
                Wchain[k, r] = STeff[k, r]

            elif sched_k == SchedStrategy.PS:
                # Processor sharing
                if r in ocl:
                    totQlen = np.sum(Qchain_in[k, nnzclasses])
                    Wchain[k, r] = STeff[k, r] * (1 + totQlen)
                else:
                    Wchain[k, r] = STeff[k, r] * (1 + interpTotArvlQlen[k])
                    if len(ccl) > 0:
                        # Add linearizer correction
                        Wchain[k, r] += STeff[k, r] * ((Nt - 1) * gamma[r, k] if gamma.ndim == 2 and r < gamma.shape[0] and k < gamma.shape[1] else 0)

            elif sched_k in [SchedStrategy.FCFS, SchedStrategy.SIRO, SchedStrategy.LCFSPR]:
                # FCFS-type stations
                if STeff[k, r] > 0:
                    ns = nservers[k]

                    # Compute multiserver correction Bk
                    if ns > 1 and not np.isinf(ns):
                        rho_sum = np.sum(deltaclass[nnzclasses] * Xchain_in[nnzclasses] *
                                        Vchain_in[k, nnzclasses] * STeff[k, nnzclasses])
                        if rho_sum < 0.75:
                            Bk = deltaclass * Xchain_in * Vchain_in[k, :] * STeff[k, :] / ns
                        else:
                            Bk = (deltaclass * Xchain_in * Vchain_in[k, :] * STeff[k, :] / ns) ** ns
                        Bk = np.where(np.isfinite(Bk), Bk, 1.0)
                    else:
                        Bk = np.ones(K)

                    # Multi-server correction with serial think time
                    if ns > 1:
                        Wchain[k, r] = STeff[k, r] * (ns - 1)
                    else:
                        Wchain[k, r] = 0.0

                    # Add service time
                    Wchain[k, r] += STeff[k, r]

                    # Add queueing delay
                    if r in ocl:
                        Wchain[k, r] += (STeff[k, r] * deltaclass[r] * stationaryQlen[k, r] * Bk[r] +
                                        np.sum(STeff[k, sd] * Bk[sd] * stationaryQlen[k, sd]))
                    else:
                        Wchain[k, r] += (STeff[k, r] * selfArvlQlenSeenByClosed[k, r] * Bk[r] +
                                        np.sum(STeff[k, sd] * Bk[sd] * stationaryQlen[k, sd]))
            else:
                # Default: same as FCFS
                Wchain[k, r] = STeff[k, r] * (1 + selfArvlQlenSeenByClosed[k, r])

    return Wchain, STeff


def solver_amvald(
    sn,  # NetworkStruct
    Lchain: np.ndarray,
    STchain: np.ndarray,
    Vchain: np.ndarray,
    alpha: np.ndarray,
    Nchain: np.ndarray,
    SCVchain: np.ndarray,
    refstatchain: np.ndarray,
    options: Optional[AmvaldOptions] = None
) -> AmvaldResult:
    """
    Approximate MVA for Load-Dependent networks (AMVA-LD).

    Implements iterative approximate MVA with under-relaxation for
    load-dependent queueing networks.

    Args:
        sn: Network structure
        Lchain: Demands per chain (M x K)
        STchain: Service times per chain (M x K)
        Vchain: Visit ratios per chain (M x K)
        alpha: Class-to-chain mapping
        Nchain: Population per chain (K,)
        SCVchain: Squared coefficient of variation (M x K)
        refstatchain: Reference station per chain (K,)
        options: Solver options

    Returns:
        AmvaldResult with Q, U, R, T, C, X, lG, totiter
    """
    if options is None:
        options = AmvaldOptions()

    M = sn.nstations
    K = sn.nchains

    Nchain = np.asarray(Nchain).flatten()
    Nt = np.sum(Nchain[np.isfinite(Nchain)])
    delta = (Nt - 1) / Nt if Nt > 0 else 0.0
    deltaclass = (Nchain - 1) / Nchain
    deltaclass[np.isinf(Nchain)] = 1.0
    deltaclass[Nchain <= 0] = 0.0

    tol = options.iter_tol
    nservers = np.asarray(sn.nservers).flatten() if sn.nservers is not None else np.ones(M)
    lldscaling = sn.lldscaling if hasattr(sn, 'lldscaling') and sn.lldscaling is not None else None
    sched = sn.sched if sn.sched else {}

    # Initialize Q, X, U
    Uchain = np.zeros((M, K))
    Tchain = np.zeros((M, K))

    # Balanced initialization for queue lengths
    Qchain = np.ones((M, K))
    Qchain = Qchain / np.sum(Qchain, axis=0, keepdims=True) * Nchain.reshape(1, -1)
    Qchain[np.isinf(Qchain)] = 0.0
    Qchain[:, Nchain == 0] = 0.0

    # For open classes, zero out at reference station
    for r in range(K):
        if np.isinf(Nchain[r]):
            refstat = int(refstatchain[r, 0]) if refstatchain.ndim > 1 else int(refstatchain[r])
            if refstat < M:
                Qchain[refstat, r] = 0.0

    nnzclasses = np.where(Nchain > 0)[0]
    Xchain = 1.0 / np.sum(STchain, axis=0)
    Xchain[np.isinf(Xchain)] = 0.0

    # For open classes, use arrival rate
    for r in range(K):
        if np.isinf(Nchain[r]):
            refstat = int(refstatchain[r, 0]) if refstatchain.ndim > 1 else int(refstatchain[r])
            if refstat < M and STchain[refstat, r] > 0:
                Xchain[r] = 1.0 / STchain[refstat, r]

    # Initialize utilization
    for k in range(M):
        for r in nnzclasses:
            if np.isinf(nservers[k]):
                Uchain[k, r] = Vchain[k, r] * STchain[k, r] * Xchain[r]
            else:
                Uchain[k, r] = Vchain[k, r] * STchain[k, r] * Xchain[r] / nservers[k]

    # Initialize correction factors
    gamma = np.zeros((K, M))
    tau = np.zeros((K, K))

    # Main iteration loop
    omicron = 0.5  # under-relaxation parameter
    totiter = 0
    outer_iter = 0
    max_outer_iter = int(np.sqrt(options.iter_max))

    QchainOuter_1 = Qchain + np.inf

    while (outer_iter < 2 or np.max(np.abs(Qchain - QchainOuter_1)) > tol) and outer_iter < max_outer_iter:
        outer_iter += 1
        QchainOuter_1 = Qchain.copy()
        XchainOuter_1 = Xchain.copy()
        UchainOuter_1 = Uchain.copy()

        inner_iter = 0
        Qchain_1 = Qchain + np.inf

        while (inner_iter < 2 or np.max(np.abs(Qchain - Qchain_1)) > tol) and inner_iter <= max_outer_iter:
            inner_iter += 1

            Qchain_1 = Qchain.copy()
            Xchain_1 = Xchain.copy()
            Uchain_1 = Uchain.copy()

            # Forward step
            Wchain, STeff = solver_amvald_forward(
                M, K, nservers, lldscaling, sched, gamma, tau,
                Qchain_1, Xchain_1, Uchain_1, STchain, Vchain, Nchain, options
            )

            totiter += 1
            if totiter > options.iter_max:
                break

            # Update metrics
            Cchain_s = np.zeros(K)
            Rchain = np.zeros((M, K))

            for r in nnzclasses:
                if np.sum(Wchain[:, r]) == 0:
                    Xchain[r] = 0.0
                else:
                    if np.isinf(Nchain[r]):
                        Cchain_s[r] = np.dot(Vchain[:, r], Wchain[:, r])
                        # X remains constant for open classes
                    elif Nchain[r] == 0:
                        Xchain[r] = 0.0
                        Cchain_s[r] = 0.0
                    else:
                        Cchain_s[r] = np.dot(Vchain[:, r], Wchain[:, r])
                        Xchain[r] = omicron * Nchain[r] / Cchain_s[r] + (1 - omicron) * Xchain_1[r]

                for k in range(M):
                    Rchain[k, r] = Vchain[k, r] * Wchain[k, r]
                    Qchain[k, r] = omicron * Xchain[r] * Vchain[k, r] * Wchain[k, r] + (1 - omicron) * Qchain_1[k, r]
                    Tchain[k, r] = Xchain[r] * Vchain[k, r]
                    Uchain[k, r] = omicron * Vchain[k, r] * STeff[k, r] * Xchain[r] + (1 - omicron) * Uchain_1[k, r]

        # Update gamma for linearizer
        ccl = np.where(np.isfinite(Nchain) & (Nchain > 0))[0]
        for s in range(K):
            if np.isfinite(Nchain[s]) and Nchain[s] > 0:
                for k in range(M):
                    gamma[s, k] = np.sum(Qchain_1[k, :]) / (Nt - 1) - np.sum(QchainOuter_1[k, :]) / Nt if Nt > 1 else 0.0
                for r in nnzclasses:
                    tau[s, r] = Xchain_1[r] - XchainOuter_1[r]

    # Utilization capping
    for k in range(M):
        sched_k = sched.get(k, SchedStrategy.FCFS)
        if sched_k in [SchedStrategy.FCFS, SchedStrategy.SIRO, SchedStrategy.PS, SchedStrategy.LCFSPR]:
            U_sum = np.sum(Uchain[k, :])
            if U_sum > 1:
                for r in range(K):
                    if Vchain[k, r] * STeff[k, r] > 0:
                        denom = np.sum(Vchain[k, :] * STeff[k, :] * Xchain)
                        if denom > 0:
                            Uchain[k, r] = min(1.0, U_sum) * Vchain[k, r] * STeff[k, r] * Xchain[r] / denom

    # Compute final response times
    Rchain = np.zeros((M, K))
    for k in range(M):
        for r in range(K):
            if Tchain[k, r] > 0:
                Rchain[k, r] = Qchain[k, r] / Tchain[k, r]

    # Clean up non-finite values
    Xchain = np.where(np.isfinite(Xchain), Xchain, 0.0)
    Uchain = np.where(np.isfinite(Uchain), Uchain, 0.0)
    Rchain = np.where(np.isfinite(Rchain), Rchain, 0.0)

    # Zero out results for zero population chains
    for j in range(K):
        if Nchain[j] == 0:
            Xchain[j] = 0.0
            Uchain[:, j] = 0.0
            Qchain[:, j] = 0.0
            Rchain[:, j] = 0.0
            Tchain[:, j] = 0.0

    # Compute cycle times
    C_result = np.sum(Rchain, axis=0).reshape(1, -1)

    # Estimate normalizing constant
    ccl = np.where(np.isfinite(Nchain))[0]
    Nclosed = Nchain[ccl]
    Xclosed = Xchain[ccl]
    valid = Xclosed > options.tol
    if np.any(valid):
        lG = -np.sum(Nclosed[valid] * np.log(Xclosed[valid]))
    else:
        lG = np.nan

    return AmvaldResult(
        Q=Qchain,
        U=Uchain,
        R=Rchain,
        T=Tchain,
        C=C_result,
        X=Xchain.reshape(1, -1),
        lG=lG,
        totiter=totiter
    )


__all__ = [
    'solver_amvald',
    'solver_amvald_forward',
    'AmvaldOptions',
    'AmvaldResult',
    'SchedStrategy',
]
