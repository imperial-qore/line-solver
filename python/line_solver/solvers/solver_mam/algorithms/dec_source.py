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
from ....api.solvers.mam.handler import (
    solver_mam_basic as handler_solver_mam_basic,
    SolverMAMOptions as HandlerOptions,
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

        Uses the handler implementation that matches MATLAB solver_mam_basic.m

        Args:
            sn: NetworkStruct from Network.compileStruct()
            options: DecSourceOptions or SolverMAMOptions (optional)

        Returns:
            MAMResult with QN, UN, RN, TN metrics
        """
        start_time = time.time()

        # Convert options to handler options
        handler_opts = HandlerOptions()
        if options is not None:
            if hasattr(options, 'tol'):
                handler_opts.tol = options.tol
            if hasattr(options, 'max_iter'):
                handler_opts.iter_max = options.max_iter
            if hasattr(options, 'verbose'):
                handler_opts.verbose = options.verbose
            if hasattr(options, 'space_max'):
                handler_opts.space_max = options.space_max

        # Call the handler implementation (matches MATLAB solver_mam_basic.m)
        handler_result = handler_solver_mam_basic(sn, handler_opts)

        runtime = time.time() - start_time

        # Convert handler result to MAMResult
        M = handler_result.Q.shape[0]
        K = handler_result.Q.shape[1]

        # TN should be (M x K) for station-class throughputs
        TN = handler_result.T
        if TN.ndim == 1:
            TN = TN.reshape(1, -1)

        return MAMResult(
            QN=handler_result.Q,
            UN=handler_result.U,
            RN=handler_result.R,
            TN=TN,
            CN=handler_result.C,
            XN=handler_result.X,
            totiter=handler_result.it,
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
