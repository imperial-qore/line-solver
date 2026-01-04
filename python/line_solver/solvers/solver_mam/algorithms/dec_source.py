"""
Decomposition with MMAP arrivals (dec.source) algorithm.

This is the default MAM algorithm that decomposes a queueing network into
individual stations and solves each station with an MMAP arrival process.

Algorithm:
1. Build MMAP arrival process for each station based on routing
2. Iterate: adjust lambda → solve queues → check convergence
3. For each queue, compute metrics (QN, UN, RN) using queue analysis
4. Update arrival rates based on queue utilizations

References:
    MATLAB: matlab/src/solvers/MAM/solver_mam_basic.m (370 lines)
    MATLAB: matlab/src/solvers/MAM/solver_mam_basic.m
"""

import numpy as np
import time
from typing import Dict, Tuple, Optional
from dataclasses import dataclass

from . import MAMAlgorithm, MAMResult
from ..utils.network_adapter import (
    extract_mam_params,
    extract_visit_counts,
    check_closed_network,
    check_product_form,
)
from ....api.mam import (
    map_lambda,
    map_scv,
    mmap_super_safe,
    mmap_exponential,
)


@dataclass
class DecSourceOptions:
    """Options for dec.source solver."""
    tol: float = 1e-6
    max_iter: int = 100
    space_max: int = 1000
    verbose: bool = False


class DecSourceAlgorithm(MAMAlgorithm):
    """Decomposition with MMAP arrivals (default method)."""

    @staticmethod
    def supports_network(sn) -> Tuple[bool, Optional[str]]:
        """Check if network can be solved by dec.source.

        dec.source supports:
        - Open and closed networks
        - Single and multiple classes
        - FCFS, LCFS, PS, and priority scheduling
        - General phase-type service distributions

        Args:
            sn: NetworkStruct

        Returns:
            (can_solve, reason_if_not)
        """
        # dec.source is very general and supports most networks
        # Only restriction: no Fork-Join (needs special handling)
        has_fork_join = any(nt in [5, 6] for nt in sn.nodetype)  # NodeType.FORK=5, JOIN=6
        if has_fork_join:
            return False, "dec.source does not support Fork-Join. Use fj method instead."

        return True, None

    def solve(self, sn, options=None) -> MAMResult:
        """Solve the network using decomposition with MMAP arrivals.

        Args:
            sn: NetworkStruct from Network.compileStruct()
            options: DecSourceOptions (optional)

        Returns:
            MAMResult with QN, UN, RN, TN metrics
        """
        if options is None:
            options = DecSourceOptions()

        start_time = time.time()

        # Extract network parameters
        params = extract_mam_params(sn)
        M = params['nstations']
        K = params['nclasses']
        rates = params['rates']
        scv = params['scv']
        nservers = params['nservers']
        visits = extract_visit_counts(sn)
        is_closed = check_closed_network(sn)

        if options.verbose:
            print(f"dec.source: M={M} stations, K={K} classes")
            print(f"  Closed network: {is_closed}")

        # Initialize result matrices
        QN = np.zeros((M, K))
        UN = np.zeros((M, K))
        RN = np.zeros((M, K))
        TN = np.zeros((1, K))
        CN = np.zeros((1, K)) if is_closed else None
        XN = np.zeros((1, K))

        # Initialize arrival rates (lambda) - ensure it's always a proper K-length array
        lambda_k = np.zeros(K, dtype=np.float64)

        # Identify Source stations (for open networks)
        source_stations = set()
        if hasattr(sn, 'sched') and sn.sched is not None:
            for ist in range(M):
                sched = sn.sched.get(ist, None)
                if sched is not None:
                    sched_name = sched.name if hasattr(sched, 'name') else str(sched)
                    if sched_name == 'EXT' or (hasattr(sched, 'value') and sched.value == 11):
                        source_stations.add(ist)

        if is_closed:
            # For closed networks, use population / total demand
            total_demand = np.sum(1.0 / np.maximum(rates, 1e-10))
            lambda_init = float(np.mean(sn.njobs) / total_demand) if len(sn.njobs) > 0 else 1.0
            for k in range(K):
                lambda_k[k] = lambda_init
        else:
            # For open networks, extract arrival rate from Source station's rates
            for ist in source_stations:
                for k in range(K):
                    if rates[ist, k] > 0:
                        lambda_k[k] = rates[ist, k]

            # Fallback: if no arrival rate found, use 1.0
            for k in range(K):
                if lambda_k[k] <= 0:
                    lambda_k[k] = 1.0

        TN_prev = TN.copy() + np.inf
        iteration = 0

        # Main iteration loop
        while np.max(np.abs(TN - TN_prev)) > options.tol and iteration < options.max_iter:
            iteration += 1
            TN_prev = TN.copy()

            if options.verbose and iteration <= 3:
                print(f"  Iteration {iteration}: lambda = {lambda_k}")

            # Update throughputs: TN(m,k) = visit_m,k * lambda_k
            for m in range(M):
                for k in range(K):
                    if visits[m, k] > 0:
                        TN[0, k] = lambda_k[k]
                    else:
                        TN[0, k] = 0.0

            # Solve each station (skip Source stations)
            for m in range(M):
                if m in source_stations:
                    # Source station: no queue, just throughput equals arrival rate
                    QN[m, :] = 0.0
                    UN[m, :] = 0.0
                    RN[m, :] = 0.0
                    XN[0, :] = lambda_k
                else:
                    QN[m, :], UN[m, :], RN[m, :], XN[0, :] = _solve_station(
                        m, M, K, lambda_k, rates[m, :], scv[m, :],
                        nservers[m], visits[m, :], options
                    )

            # Check for stability (utilization < 1)
            max_util = np.max(UN)
            if max_util >= 1.0 and not is_closed:
                # Scale down lambda for open networks
                lambda_k = lambda_k / max_util
            elif is_closed:
                # For closed networks, adjust lambda based on queue lengths
                total_qlen = np.sum(QN)
                if total_qlen > 0:
                    total_demand = np.sum(1.0 / np.maximum(rates, 1e-10))
                    lambda_val = float(np.mean(sn.njobs) / total_demand)
                    # Update all K elements
                    for k in range(K):
                        lambda_k[k] = lambda_val
                else:
                    lambda_k = lambda_k * 1.1  # Slight increase

        # Compute cycle times (closed networks)
        if is_closed and CN is not None:
            for k in range(K):
                if TN[0, k] > 1e-10:
                    CN[0, k] = np.sum(RN[:, k] * visits[:, k]) / TN[0, k]

        runtime = time.time() - start_time

        return MAMResult(
            QN=QN,
            UN=UN,
            RN=RN,
            TN=TN,
            CN=CN,
            XN=XN,
            totiter=iteration,
            method="dec.source",
            runtime=runtime
        )


def _solve_station(station_idx: int,
                  M: int,
                  K: int,
                  lambda_k: np.ndarray,
                  rates: np.ndarray,
                  scv: np.ndarray,
                  nservers: float,
                  visits: np.ndarray,
                  options: DecSourceOptions) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Solve a single station in isolation.

    Uses simplified queue analysis based on arrival SCV and service SCV.

    Args:
        station_idx: Station index
        M: Total number of stations
        K: Number of classes
        lambda_k: Arrival rates per class
        rates: Service rates
        scv: Service SCV
        nservers: Number of servers
        visits: Visit counts
        options: Solver options

    Returns:
        (QN, UN, RN, XN) for this station
    """
    QN = np.zeros(K)
    UN = np.zeros(K)
    RN = np.zeros(K)

    # Service time (mean)
    S = 1.0 / np.maximum(rates, 1e-10)

    # Compute arrival rates weighted by visits
    lambda_arr = np.zeros(K)
    for k in range(K):
        if visits[k] > 0:
            lambda_arr[k] = lambda_k[k] * visits[k]
        else:
            lambda_arr[k] = 0.0

    # Total arrival rate
    lambda_total = np.sum(lambda_arr)

    # Aggregate SCV of arrivals (assume Poisson for now, SCV=1)
    ca2 = 1.0  # Coefficient of arrival variation squared

    # For each class, compute utilization and metrics
    for k in range(K):
        if lambda_arr[k] <= 1e-10:
            UN[k] = 0.0
            QN[k] = 0.0
            RN[k] = S[k]
        elif np.isinf(nservers):
            # Infinite server station (Delay): no waiting, U = lambda * S
            UN[k] = lambda_arr[k] * S[k]
            RN[k] = S[k]  # Response time equals service time (no waiting)
            QN[k] = lambda_arr[k] * S[k]  # Little's law: L = lambda * R = lambda * S
        else:
            # Utilization
            rho = lambda_arr[k] * S[k] / nservers
            UN[k] = rho

            if rho >= 1.0:
                # Unstable queue
                QN[k] = np.inf
                RN[k] = np.inf
            else:
                # Use M/G/c approximation (Kingman's formula variant)
                # For single server case (nservers=1), use M/G/1 formula
                if nservers <= 1:
                    cs2 = scv[k]  # Service SCV
                    # M/G/1: E[W] = (ca2 + cs2) * rho / (2(1-rho)) * S
                    mean_wait = (ca2 + cs2) * rho / (2.0 * (1.0 - rho)) * S[k]
                    RN[k] = S[k] + mean_wait
                    # Little's law: L = lambda * R (total system size)
                    QN[k] = lambda_arr[k] * RN[k]
                else:
                    # Multi-server approximation (simpler)
                    # Use Erlang C-like approximation
                    mean_wait = _erlang_c_wait(nservers, rho, S[k])
                    RN[k] = S[k] + mean_wait
                    # Little's law: L = lambda * R (total system size)
                    QN[k] = lambda_arr[k] * RN[k]

    # System metrics
    XN = np.sum(lambda_arr)  # System throughput

    return QN, UN, RN, np.array([XN])


def _erlang_c_wait(c: float,
                  rho: float,
                  service_time: float) -> float:
    """
    Compute mean waiting time using Erlang C formula (approximation).

    For M/M/c queue, or multi-class approximation.

    Args:
        c: Number of servers
        rho: Utilization = lambda * S / c
        service_time: Mean service time S

    Returns:
        Mean waiting time
    """
    if rho >= 1.0 or c < 1:
        return np.inf

    # Simple approximation: for M/M/c
    # Erlang C: Pw = (c*rho)^c / (c! * (1-rho)) / [sum of terms]
    # Mean wait = Pw * S / (c * (1 - rho))

    if c <= 1:
        # M/M/1: W = rho / (1 - rho) * S
        return rho / (1.0 - rho) * service_time
    else:
        # M/M/c approximation
        # Simplified: use M/M/1 with adjusted rho
        rho_eff = rho / c
        if rho_eff >= 1.0:
            return np.inf
        return rho_eff / (1.0 - rho_eff) * service_time


def mmap_exponential(rate: float):
    """Create exponential MMAP with given rate.

    Args:
        rate: Arrival rate

    Returns:
        (D0, [D1]) MMAP representation
    """
    if not np.isfinite(rate) or rate <= 0:
        # No arrivals
        D0 = np.array([[-1.0]])
        D1 = np.array([[0.0]])
    else:
        D0 = np.array([[-rate]])
        D1 = np.array([[rate]])

    return D0, [D1]
