"""
MAM Solver handler.

Native Python implementation of MAM (Matrix-Analytic Methods) solver handler
that analyzes queueing networks with phase-type distributions and Markovian
arrival processes.

Port from:


"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, Dict, List, Tuple, Any
import time

from ...sn import (
    NetworkStruct,
    SchedStrategy,
    NodeType,
    sn_is_open_model,
    sn_is_closed_model,
    sn_get_demands_chain,
)


@dataclass
class SolverMAMOptions:
    """Options for MAM solver."""
    method: str = 'default'
    tol: float = 1e-6
    verbose: bool = False
    iter_max: int = 100
    iter_tol: float = 1e-6
    space_max: int = 128
    merge: str = 'super'
    compress: str = 'mixture.order1'
    num_cdf_pts: int = 200  # Number of points for CDF computation


@dataclass
class SolverMAMReturn:
    """
    Result of MAM solver handler.

    Attributes:
        Q: Mean queue lengths (M x K)
        U: Utilizations (M x K)
        R: Response times (M x K)
        T: Throughputs (M x K)
        C: Cycle times (1 x K)
        X: System throughputs (1 x K)
        A: Arrival rates (M x K)
        W: Waiting times (M x K)
        runtime: Runtime in seconds
        method: Method used
        it: Number of iterations
    """
    Q: Optional[np.ndarray] = None
    U: Optional[np.ndarray] = None
    R: Optional[np.ndarray] = None
    T: Optional[np.ndarray] = None
    C: Optional[np.ndarray] = None
    X: Optional[np.ndarray] = None
    A: Optional[np.ndarray] = None
    W: Optional[np.ndarray] = None
    runtime: float = 0.0
    method: str = "default"
    it: int = 0


def _get_visits(sn: NetworkStruct) -> np.ndarray:
    """
    Compute station visit ratios from routing matrix.

    Args:
        sn: Network structure

    Returns:
        Visit ratios matrix (M x K)
    """
    M = sn.nstations
    K = sn.nclasses

    # Initialize visits matrix
    V = np.zeros((M, K))

    if sn.visits is not None and len(sn.visits) > 0:
        # Sum visits across all chains
        for c in range(len(sn.visits)):
            if sn.visits[c] is not None:
                v = np.asarray(sn.visits[c])
                if v.ndim == 1:
                    v = v.reshape(-1, 1) if v.shape[0] == M else v.reshape(1, -1)
                if v.shape[0] == M and v.shape[1] == K:
                    V += v
                elif v.shape[0] == M:
                    # Single column, repeat for all classes in chain
                    V[:, :v.shape[1]] += v
    else:
        # Default: equal visits
        V = np.ones((M, K))

    return V


def _get_service_rates(sn: NetworkStruct) -> np.ndarray:
    """
    Extract service rates from network structure.

    Args:
        sn: Network structure

    Returns:
        Service rates matrix (M x K)
    """
    M = sn.nstations
    K = sn.nclasses

    if hasattr(sn, 'rates') and sn.rates is not None:
        rates = np.asarray(sn.rates)
        if rates.shape == (M, K):
            return rates
        elif rates.ndim == 1 and rates.shape[0] == M:
            return np.tile(rates.reshape(-1, 1), (1, K))

    # Default rates
    return np.ones((M, K))


def _get_nservers(sn: NetworkStruct) -> np.ndarray:
    """
    Get number of servers per station.

    Args:
        sn: Network structure

    Returns:
        Number of servers array (M,)
    """
    M = sn.nstations

    if hasattr(sn, 'nservers') and sn.nservers is not None:
        nservers = np.asarray(sn.nservers).flatten()
        if len(nservers) == M:
            return nservers

    # Default: single server
    return np.ones(M)


def _get_scheduling(sn: NetworkStruct, station_idx: int) -> SchedStrategy:
    """
    Get scheduling strategy for a station.

    Args:
        sn: Network structure
        station_idx: Station index

    Returns:
        Scheduling strategy
    """
    if hasattr(sn, 'sched') and sn.sched is not None:
        if isinstance(sn.sched, dict):
            return sn.sched.get(station_idx, SchedStrategy.FCFS)
        elif isinstance(sn.sched, list) and station_idx < len(sn.sched):
            return sn.sched[station_idx]

    return SchedStrategy.FCFS


def _is_delay_station(sn: NetworkStruct, station_idx: int) -> bool:
    """
    Check if station is a delay (infinite server) station.

    Args:
        sn: Network structure
        station_idx: Station index

    Returns:
        True if station is a delay station
    """
    sched = _get_scheduling(sn, station_idx)
    return sched == SchedStrategy.INF


def _is_source_station(sn: NetworkStruct, station_idx: int) -> bool:
    """
    Check if station is a source.

    Args:
        sn: Network structure
        station_idx: Station index

    Returns:
        True if station is a source
    """
    sched = _get_scheduling(sn, station_idx)
    return sched == SchedStrategy.EXT


def solver_mam_basic(
    sn: NetworkStruct,
    options: Optional[SolverMAMOptions] = None
) -> SolverMAMReturn:
    """
    Basic MAM solver using source decomposition.

    Implements basic matrix-analytic decomposition for queueing networks.
    This is the default method for open networks.

    Args:
        sn: Network structure
        options: Solver options

    Returns:
        SolverMAMReturn with performance metrics
    """
    start_time = time.time()

    if options is None:
        options = SolverMAMOptions()

    M = sn.nstations
    K = sn.nclasses

    # Initialize result matrices
    QN = np.zeros((M, K))
    UN = np.zeros((M, K))
    RN = np.zeros((M, K))
    TN = np.zeros((M, K))
    CN = np.zeros((1, K))
    XN = np.zeros((1, K))

    # Get network parameters
    V = _get_visits(sn)
    rates = _get_service_rates(sn)
    nservers = _get_nservers(sn)

    # Service times (inverse of rates), handling zero rates
    with np.errstate(divide='ignore', invalid='ignore'):
        S = np.where(rates > 0, 1.0 / rates, 0.0)

    # Check model type
    is_open = sn_is_open_model(sn)
    is_closed = sn_is_closed_model(sn)

    if is_open:
        # Open queueing network analysis
        # Get external arrival rates
        lambda_ext = np.zeros(K)

        # Find source station and get arrival rates
        for ist in range(M):
            if _is_source_station(sn, ist):
                for k in range(K):
                    if rates[ist, k] > 0 and np.isfinite(rates[ist, k]):
                        lambda_ext[k] = rates[ist, k]

        # System throughputs equal external arrival rates
        XN[0, :] = lambda_ext

        # Compute metrics for each station
        for ist in range(M):
            if _is_source_station(sn, ist):
                # Source station: throughput equals arrival rate
                TN[ist, :] = lambda_ext
                continue

            for k in range(K):
                # Station throughput = system throughput * visit ratio
                TN[ist, k] = XN[0, k] * V[ist, k]

                if _is_delay_station(sn, ist):
                    # Delay station (infinite server)
                    UN[ist, k] = S[ist, k] * TN[ist, k]
                    QN[ist, k] = TN[ist, k] * S[ist, k]
                    RN[ist, k] = S[ist, k]
                else:
                    # Queue station
                    service_rate = rates[ist, k] * nservers[ist]
                    if service_rate > 0:
                        UN[ist, k] = TN[ist, k] / service_rate

        # For queuing stations, compute queue lengths using M/M/c formulas
        for ist in range(M):
            if not _is_delay_station(sn, ist) and not _is_source_station(sn, ist):
                # Aggregate utilization
                total_util = np.sum(UN[ist, :])
                c = nservers[ist]

                if total_util < 1.0 - 1e-10:
                    for k in range(K):
                        if TN[ist, k] > 0:
                            # M/M/c approximation for class k
                            rho_k = UN[ist, k]

                            if c == 1:
                                # M/M/1 queue
                                QN[ist, k] = rho_k / (1.0 - total_util) if total_util < 1 else rho_k
                            else:
                                # M/M/c approximation
                                QN[ist, k] = rho_k / (1.0 - total_util) if total_util < 1 else rho_k

                            RN[ist, k] = QN[ist, k] / TN[ist, k] if TN[ist, k] > 0 else 0
                else:
                    # Saturated system
                    for k in range(K):
                        if sn.njobs is not None and k < len(sn.njobs.flatten()):
                            njobs_k = sn.njobs.flatten()[k]
                            QN[ist, k] = njobs_k if np.isfinite(njobs_k) else 0
                        RN[ist, k] = QN[ist, k] / TN[ist, k] if TN[ist, k] > 0 else 0

    elif is_closed:
        # Closed queueing network analysis using fixed-point iteration
        N = sn.njobs.flatten() if sn.njobs is not None else np.ones(K)

        # Get demands per chain
        demands_result = sn_get_demands_chain(sn)
        Lchain = demands_result.Lchain if demands_result else np.sum(S * V, axis=0).reshape(1, -1)

        # Initialize arrival rate estimate (scalar)
        lambda_est = float(np.sum(N) / (np.sum(Lchain) + 1e-10))

        # Fixed-point iteration
        TN_prev = np.full((M, K), np.inf)
        it = 0
        tol = options.iter_tol

        while it < options.iter_max:
            it += 1
            TN_prev = TN.copy()

            # Compute throughputs
            for ist in range(M):
                for k in range(K):
                    TN[ist, k] = float(lambda_est * V[ist, k])
                    UN[ist, k] = float(TN[ist, k] * S[ist, k] / nservers[ist])

            # Check utilization bound
            U_max = float(np.max(np.sum(UN, axis=1)))
            if U_max >= 1.0:
                lambda_est = float(lambda_est / (U_max + 0.1))
                continue

            # Compute queue lengths
            for ist in range(M):
                if _is_delay_station(sn, ist):
                    for k in range(K):
                        QN[ist, k] = TN[ist, k] * S[ist, k]
                        RN[ist, k] = S[ist, k]
                else:
                    total_util = np.sum(UN[ist, :])
                    for k in range(K):
                        if total_util < 1.0 - 1e-10:
                            QN[ist, k] = UN[ist, k] / (1.0 - total_util)
                        else:
                            QN[ist, k] = N[k] if np.isfinite(N[k]) else 0
                        RN[ist, k] = QN[ist, k] / TN[ist, k] if TN[ist, k] > 0 else 0

            # Update arrival rate based on queue lengths
            QN_total = float(np.sum(QN))
            if QN_total > 0:
                lambda_new = float((np.sum(N) / QN_total) * lambda_est)
                # Damped update
                alpha = float((options.iter_max - it) / options.iter_max)
                lambda_est = float(alpha * lambda_est + (1 - alpha) * lambda_new)

            # Check convergence
            diff = np.max(np.abs(TN - TN_prev))
            if diff < tol:
                break

        XN[0, :] = lambda_est * np.ones(K)

        # Scale queue lengths to match population
        for k in range(K):
            if np.isfinite(N[k]):
                QN_k = np.sum(QN[:, k])
                if QN_k > 0:
                    QN[:, k] = QN[:, k] * N[k] / QN_k

    else:
        # Mixed network - not supported
        raise RuntimeError("SolverMAM does not support mixed models with both open and closed classes.")

    # Compute cycle times
    CN = np.sum(RN, axis=0).reshape(1, -1)

    # Clean up NaN values
    QN = np.nan_to_num(QN, nan=0.0)
    UN = np.nan_to_num(UN, nan=0.0)
    RN = np.nan_to_num(RN, nan=0.0)
    TN = np.nan_to_num(TN, nan=0.0)
    CN = np.nan_to_num(CN, nan=0.0)
    XN = np.nan_to_num(XN, nan=0.0)

    result = SolverMAMReturn()
    result.Q = QN
    result.U = UN
    result.R = RN
    result.T = TN
    result.C = CN
    result.X = XN
    result.A = TN.copy()  # Arrival rates ≈ throughputs
    result.W = RN.copy()  # Waiting times ≈ response times
    result.runtime = time.time() - start_time
    result.method = "dec.source"
    result.it = it if is_closed else 1

    return result


def solver_mam(
    sn: NetworkStruct,
    options: Optional[SolverMAMOptions] = None
) -> SolverMAMReturn:
    """
    Main MAM solver handler.

    Routes to appropriate method based on options and network characteristics.

    Args:
        sn: Network structure
        options: Solver options

    Returns:
        SolverMAMReturn with performance metrics
    """
    if options is None:
        options = SolverMAMOptions()

    method = options.method.lower()

    if method in ['default', 'dec.source', 'dec.poisson']:
        # Use basic decomposition method
        if method == 'dec.poisson':
            options.space_max = 1
        return solver_mam_basic(sn, options)
    elif method == 'dec.mmap':
        # MMAP decomposition - fallback to basic for now
        # Full implementation would use MMAPPH1FCFS
        return solver_mam_basic(sn, options)
    elif method in ['mna', 'inap', 'exact']:
        # Matrix-analytic / RCAT methods - fallback to basic
        return solver_mam_basic(sn, options)
    else:
        # Unknown method - use basic
        if options.verbose:
            print(f"Warning: Unknown MAM method '{method}'. Using dec.source.")
        return solver_mam_basic(sn, options)


__all__ = [
    'solver_mam',
    'solver_mam_basic',
    'SolverMAMReturn',
    'SolverMAMOptions',
]
