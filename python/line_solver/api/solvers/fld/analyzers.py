"""
FLD Solver analyzers.

Native Python implementation of FLD solver analyzers that orchestrate
method selection and provide the main entry point for fluid analysis.

Port from:



"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any, Tuple
import time

from ...sn import (
    NetworkStruct,
    SchedStrategy,
    NodeType,
)
# Import standalone FluFluQueue using direct file import to avoid butools package dependency
# This is necessary because the queues __init__.py imports modules that require external butools
import importlib.util as _imputil
import os as _os
_fluflu_path = _os.path.join(_os.path.dirname(__file__), '..', '..', 'butools', 'queues', 'flufluqueue.py')
_fluflu_spec = _imputil.spec_from_file_location('flufluqueue', _fluflu_path)
_fluflu_module = _imputil.module_from_spec(_fluflu_spec)
_fluflu_spec.loader.exec_module(_fluflu_module)
FluFluQueue = _fluflu_module.FluFluQueue
FluFluResultData = _fluflu_module.FluFluResult
from .handler import solver_fld, SolverFLDOptions, SolverFLDReturn


@dataclass
class FLDResult:
    """
    Result of FLD solver analysis.

    Attributes:
        QN: Mean queue lengths (M x K)
        UN: Utilizations (M x K)
        RN: Response times (M x K)
        TN: Throughputs (M x K)
        CN: Cycle times (1 x K)
        XN: System throughputs (1 x K)
        AN: Arrival rates (M x K)
        WN: Waiting times (M x K)
        QNt: Transient queue lengths
        UNt: Transient utilizations
        TNt: Transient throughputs
        t: Time vector
        odeStateVec: Final ODE state vector
        iter: Number of iterations
        runtime: Runtime in seconds
        method: Method used
    """
    QN: Optional[np.ndarray] = None
    UN: Optional[np.ndarray] = None
    RN: Optional[np.ndarray] = None
    TN: Optional[np.ndarray] = None
    CN: Optional[np.ndarray] = None
    XN: Optional[np.ndarray] = None
    AN: Optional[np.ndarray] = None
    WN: Optional[np.ndarray] = None
    QNt: Optional[List] = None
    UNt: Optional[List] = None
    TNt: Optional[List] = None
    t: Optional[np.ndarray] = None
    odeStateVec: Optional[np.ndarray] = None
    iter: int = 0
    runtime: float = 0.0
    method: str = ""


def matrix_method_analyzer(
    sn: NetworkStruct,
    options: Optional[SolverFLDOptions] = None
) -> SolverFLDReturn:
    """
    Matrix method analyzer for fluid analysis.

    Implements the matrix-based fluid analysis as per Ruuskanen et al., PEVA 151 (2021).
    Uses ODE integration with the transition rate matrix W to compute transient and
    steady-state performance metrics.

    Args:
        sn: Network structure
        options: Solver options

    Returns:
        SolverFLDReturn with performance metrics
    """
    if options is None:
        options = SolverFLDOptions(method='matrix')
    else:
        options.method = 'matrix'

    return solver_fld(sn, options)


def closing_method_analyzer(
    sn: NetworkStruct,
    options: Optional[SolverFLDOptions] = None
) -> SolverFLDReturn:
    """
    Closing method analyzer for fluid analysis.

    Uses closing approximations and state-dependent methods for networks
    that cannot be handled by the matrix method (e.g., DPS, open classes).

    Implementation based on Section 7.2 of Ruuskanen et al., PEVA 151 (2021),
    with adaptations for state-dependent scheduling.

    The closing method:
    1. Computes throughputs at each station assuming full network load
    2. Applies state-dependent approximations for scheduling strategies
    3. Solves for eventual visit probabilities
    4. Computes response times and utilizations

    Args:
        sn: Network structure
        options: Solver options

    Returns:
        SolverFLDReturn with performance metrics

    References:
        MATLAB: matlab/src/solvers/FLD/solver_fluid_closing.m
    """
    start_time = time.time()

    if options is None:
        options = SolverFLDOptions(method='closing')
    else:
        options.method = 'closing'

    M = sn.nstations
    K = sn.nclasses

    # Initialize result matrices
    Q = np.zeros((M, K))
    U = np.zeros((M, K))
    R = np.zeros((M, K))
    T = np.zeros((M, K))

    # Get scheduling strategies and visit ratios
    sched_dict = sn.sched if sn.sched else {}
    rates = sn.rates if sn.rates is not None else np.ones((M, K))
    visits = sn.visits if sn.visits else {0: np.ones((M, K))}

    # Compute service times
    S = 1.0 / (rates + 1e-10)

    # Get aggregate visit ratios
    V_agg = np.zeros((M, K))
    for c in visits:
        if isinstance(visits[c], np.ndarray):
            V_agg += visits[c]

    # Compute aggregate throughputs
    total_lambda = 0
    for k in range(K):
        if sn.njobs[k] == np.inf:  # Open class
            total_lambda += np.sum(rates[:, k])

    # Main iteration to compute throughputs
    for m in range(M):
        station_sched = sched_dict.get(m, SchedStrategy.FCFS)

        for k in range(K):
            # Compute throughput at station m for class k
            if V_agg[m, k] > 0:
                if station_sched == SchedStrategy.INF:
                    # Delay station: throughput = arrival rate
                    T[m, k] = total_lambda if k < len(rates) else 0
                    Q[m, k] = T[m, k] * S[m, k] * V_agg[m, k]
                    U[m, k] = Q[m, k]
                    R[m, k] = S[m, k] * V_agg[m, k]

                elif station_sched == SchedStrategy.PS:
                    # Processor Sharing: exact geometric approximation
                    # Simplified computation for state-dependent PS
                    T[m, k] = total_lambda * V_agg[m, k] / np.sum(V_agg[m, :])
                    rho_m = np.sum(T[m, :]) * S[m, 0] if S[m, 0] > 0 else 0

                    if rho_m < 1.0:
                        # Use approximation for queue length
                        U[m, k] = T[m, k] * S[m, k]
                        Q[m, k] = U[m, k] / (1.0 - rho_m + 1e-10)
                        R[m, k] = Q[m, k] / (T[m, k] + 1e-10)
                    else:
                        Q[m, k] = sn.njobs[k] if np.isfinite(sn.njobs[k]) else 1
                        R[m, k] = Q[m, k] / (T[m, k] + 1e-10)
                        U[m, k] = 1.0

                elif station_sched == SchedStrategy.DPS:
                    # Discriminatory Processor Sharing: weighted by schedparam
                    T[m, k] = total_lambda * V_agg[m, k] / np.sum(V_agg[m, :])

                    # Get weight for this class (default 1.0)
                    w_k = 1.0
                    if hasattr(sn, 'schedparam') and sn.schedparam is not None:
                        if m < sn.schedparam.shape[0] and k < sn.schedparam.shape[1]:
                            w_k = sn.schedparam[m, k]

                    # Compute weighted total for all classes at this station
                    weighted_total = 0.0
                    for kk in range(K):
                        w_kk = 1.0
                        if hasattr(sn, 'schedparam') and sn.schedparam is not None:
                            if m < sn.schedparam.shape[0] and kk < sn.schedparam.shape[1]:
                                w_kk = sn.schedparam[m, kk]
                        # Weight by expected queue length contribution
                        weighted_total += w_kk * V_agg[m, kk]

                    # Effective service time scaled by weight ratio
                    if weighted_total > 0:
                        weight_share = w_k * V_agg[m, k] / weighted_total
                    else:
                        weight_share = 1.0 / K

                    S_eff = S[m, k] / w_k if w_k > 0 else S[m, k]
                    rho_m = np.sum(T[m, :] * S[m, :])

                    if rho_m < 1.0:
                        U[m, k] = T[m, k] * S[m, k]
                        Q[m, k] = U[m, k] / (1.0 - rho_m + 1e-10)
                        R[m, k] = S_eff + (Q[m, k] - U[m, k]) / (T[m, k] + 1e-10) * weight_share
                    else:
                        Q[m, k] = sn.njobs[k] if np.isfinite(sn.njobs[k]) else 1
                        R[m, k] = Q[m, k] / (T[m, k] + 1e-10)
                        U[m, k] = 1.0

                else:  # FCFS, etc.
                    # General approximation for other scheduling strategies
                    # Use MVA-like approximation
                    T[m, k] = total_lambda * V_agg[m, k] / np.sum(V_agg[m, :])
                    mu_m = rates[m, k] if k < rates.shape[1] else 1.0
                    U[m, k] = T[m, k] * S[m, k] / (sn.nservers[m] if m < len(sn.nservers) else 1)

                    if U[m, k] < 1.0:
                        # Use diffusion approximation for FCFS
                        rho_m = np.sum(T[m, :]) / (sn.nservers[m] if m < len(sn.nservers) else 1)
                        Q[m, k] = U[m, k] + (U[m, k] * rho_m) / (1.0 - rho_m + 1e-10)
                        R[m, k] = Q[m, k] / (T[m, k] + 1e-10)
                    else:
                        Q[m, k] = sn.njobs[k] if np.isfinite(sn.njobs[k]) else 1
                        R[m, k] = Q[m, k] / (T[m, k] + 1e-10)
                        U[m, k] = 1.0

    # Compute cycle times
    C = np.sum(R, axis=0, keepdims=True)

    # Compute system throughputs
    X = np.zeros((1, K))
    for k in range(K):
        if np.sum(C[:, k]) > 0:
            if np.isfinite(sn.njobs[k]):
                X[0, k] = sn.njobs[k] / (np.sum(C[:, k]) + 1e-10)
            else:
                X[0, k] = np.sum(T[:, k])

    # Clean up numerical issues
    Q = np.abs(Q)
    U = np.minimum(U, 1.0)
    R = np.abs(R)
    T = np.abs(T)
    Q[np.isnan(Q)] = 0
    U[np.isnan(U)] = 0
    R[np.isnan(R)] = 0
    T[np.isnan(T)] = 0

    result = SolverFLDReturn(
        Q=Q,
        U=U,
        R=R,
        T=T,
        C=C,
        X=X,
        Qt=[Q],
        Ut=[U],
        Tt=[T],
        t=np.array([0.0]),
        odeStateVec=None,
        runtime=time.time() - start_time,
        method='closing',
        it=1
    )

    return result


def _has_dps_scheduling(sn: NetworkStruct) -> bool:
    """
    Check if network has DPS (Discriminatory Processor Sharing) scheduling.

    Args:
        sn: Network structure

    Returns:
        True if any station has DPS scheduling
    """
    sched_dict = sn.sched if sn.sched else {}
    for i in range(sn.nstations):
        station_sched = sched_dict.get(i)
        if station_sched == SchedStrategy.DPS:
            return True
    return False


def _has_open_classes(sn: NetworkStruct) -> bool:
    """
    Check if network has open classes (sources).

    Args:
        sn: Network structure

    Returns:
        True if any node is a Source
    """
    for i in range(len(sn.nodetype)):
        if sn.nodetype[i] == NodeType.SOURCE:
            return True
    return False


def solver_fld_analyzer(
    sn: NetworkStruct,
    options: Optional[SolverFLDOptions] = None
) -> FLDResult:
    """
    FLD Analyzer - main entry point for fluid analysis.

    Selects appropriate fluid method based on network characteristics
    and solver options, then performs the analysis.

    Supported methods:
        - 'matrix': Matrix-based fluid analysis (Ruuskanen et al.)
        - 'closing': Closing approximation for state-dependent networks
        - 'statedep': State-dependent fluid analysis
        - 'mfq': Markovian Fluid Queue for single-queue open systems
        - 'default': Automatic method selection

    Args:
        sn: Network structure
        options: Solver options (method, tolerance, verbosity)

    Returns:
        FLDResult with all performance metrics

    Raises:
        RuntimeError: For unsupported configurations
    """
    start_time = time.time()

    if options is None:
        options = SolverFLDOptions()

    method = options.method.lower()

    result = FLDResult()
    ret = None

    # Handle MFQ method separately (returns FLDResult directly)
    if method == 'mfq':
        return mfq_analyzer(sn, options)

    # Determine actual method
    if method == 'default':
        has_open = _has_open_classes(sn)
        has_dps = _has_dps_scheduling(sn)

        if has_open or has_dps:
            # Use closing method for open classes or DPS
            method = 'closing'
            ret = closing_method_analyzer(sn, options)
        else:
            # Use matrix method for closed networks
            method = 'matrix'
            ret = matrix_method_analyzer(sn, options)

    elif method == 'matrix':
        if _has_dps_scheduling(sn):
            if options.verbose:
                print("Warning: Matrix method does not support DPS scheduling. Using closing method.")
            method = 'closing'
            ret = closing_method_analyzer(sn, options)
        else:
            ret = matrix_method_analyzer(sn, options)

    elif method in ['closing', 'statedep']:
        ret = closing_method_analyzer(sn, options)

    else:
        # Unsupported method - try matrix as fallback
        if options.verbose:
            print(f"Warning: Unsupported FLD method '{method}'. Using matrix method.")
        method = 'matrix'
        ret = matrix_method_analyzer(sn, options)

    # Copy results from handler return
    if ret is not None:
        result.QN = ret.Q
        result.UN = ret.U
        result.RN = ret.R
        result.TN = ret.T
        result.CN = ret.C
        result.XN = ret.X
        result.AN = ret.T.copy() if ret.T is not None else None  # Arrival rates ≈ throughputs
        result.WN = ret.R.copy() if ret.R is not None else None  # Waiting times ≈ response times
        result.QNt = ret.Qt
        result.UNt = ret.Ut
        result.TNt = ret.Tt
        result.t = ret.t
        result.odeStateVec = ret.odeStateVec
        result.iter = ret.it

    result.runtime = time.time() - start_time
    result.method = method

    return result


def mfq_analyzer(
    sn: NetworkStruct,
    options: Optional[SolverFLDOptions] = None
) -> FLDResult:
    """
    MFQ (Markovian Fluid Queue) analyzer for single-queue open systems.

    Uses BUTools FluFluQueue function to compute exact steady-state queue length
    and sojourn time moments for single-queue networks with phase-type arrivals
    and service distributions.

    Applicability:
    - Single-queue open system: Source -> Queue -> Sink
    - Single-server (c=1) or infinite-server (c=Inf)
    - Exponential or phase-type distributions

    Args:
        sn: Network structure
        options: Solver options

    Returns:
        FLDResult with performance metrics

    Raises:
        RuntimeError: If topology is not suitable for MFQ
    """
    start_time = time.time()

    if options is None:
        options = SolverFLDOptions(method='mfq')

    M = sn.nstations
    K = sn.nclasses

    result = FLDResult()
    result.QN = np.zeros((M, K))
    result.UN = np.zeros((M, K))
    result.RN = np.zeros((M, K))
    result.TN = np.zeros((M, K))
    result.WN = np.zeros((M, K))
    result.AN = np.zeros((M, K))

    # Validate topology
    topology_info = _validate_mfq_topology(sn)
    if not topology_info['is_valid']:
        raise RuntimeError(f"MFQ requires single-queue topology: {topology_info['error_msg']}")

    source_station = topology_info['source_station']
    queue_station = topology_info['queue_station']

    # Process each open class
    for k in range(K):
        njobs_k = sn.njobs[0, k] if hasattr(sn.njobs, 'shape') else sn.njobs[k]
        if njobs_k == float('inf') or njobs_k == np.inf:
            # Open class
            class_result = _analyze_open_class_mfq(sn, k, source_station, queue_station, options)

            result.QN[queue_station, k] = class_result['mean_queue_length']
            result.RN[queue_station, k] = class_result['mean_response_time']
            result.TN[queue_station, k] = class_result['throughput']
            result.UN[queue_station, k] = class_result['utilization']

            # Throughput at source equals queue throughput
            result.TN[source_station, k] = class_result['throughput']

    # Create dummy state vector for compatibility
    result.odeStateVec = np.zeros((1, 1))
    result.runtime = time.time() - start_time
    result.method = 'mfq'

    return result


def _validate_mfq_topology(sn: NetworkStruct) -> dict:
    """
    Validates that the network has a single-queue topology suitable for MFQ.
    """
    info = {
        'is_valid': False,
        'error_msg': '',
        'source_node': -1,
        'queue_node': -1,
        'sink_node': -1,
        'source_station': -1,
        'queue_station': -1
    }

    # Check for open classes
    has_open_class = False
    for k in range(sn.nclasses):
        njobs_k = sn.njobs[0, k] if hasattr(sn.njobs, 'shape') else sn.njobs[k]
        if njobs_k == float('inf') or njobs_k == np.inf:
            has_open_class = True
            break

    if not has_open_class:
        info['error_msg'] = "Not an open model - all classes are closed"
        return info

    # Find Source, Queue, and Sink nodes
    source_node = -1
    queue_node = -1
    sink_node = -1
    source_count = 0
    queue_count = 0
    sink_count = 0

    for i in range(sn.nnodes):
        node_type = sn.nodetype[i]
        if node_type == NodeType.SOURCE:
            source_node = i
            source_count += 1
        elif node_type == NodeType.QUEUE:
            queue_node = i
            queue_count += 1
        elif node_type == NodeType.SINK:
            sink_node = i
            sink_count += 1

    # Validate counts
    if source_count == 0:
        info['error_msg'] = "No source node found"
        return info
    if source_count > 1:
        info['error_msg'] = f"Multiple source nodes found ({source_count})"
        return info
    if sink_count == 0:
        info['error_msg'] = "No sink node found"
        return info
    if sink_count > 1:
        info['error_msg'] = f"Multiple sink nodes found ({sink_count})"
        return info
    if queue_count == 0:
        info['error_msg'] = "No queue node found"
        return info
    if queue_count > 1:
        info['error_msg'] = f"Multiple queue nodes found ({queue_count}) - MFQ supports single queue only"
        return info

    # Get station indices
    source_station = int(sn.nodeToStation[source_node])
    queue_station = int(sn.nodeToStation[queue_node])

    # Check server count (single-server or infinite-server)
    c = sn.nservers[queue_station, 0] if hasattr(sn.nservers, 'shape') else sn.nservers[queue_station]
    if c > 1 and c < float('inf'):
        info['error_msg'] = f"Multi-server queue (c={c}) not supported - MFQ requires c=1 or c=Inf"
        return info

    info['is_valid'] = True
    info['source_node'] = source_node
    info['queue_node'] = queue_node
    info['sink_node'] = sink_node
    info['source_station'] = source_station
    info['queue_station'] = queue_station
    return info


def _analyze_open_class_mfq(
    sn: NetworkStruct,
    class_idx: int,
    source_station: int,
    queue_station: int,
    options: SolverFLDOptions
) -> dict:
    """
    Analyzes an open class using FluFluQueue.
    """
    # Extract arrival rate
    lambda_rate = sn.rates[source_station, class_idx]
    if np.isnan(lambda_rate) or lambda_rate <= 0:
        raise RuntimeError(f"No valid arrival rate for class {class_idx} at source")

    # Extract service rate (mean)
    mu_rate = sn.rates[queue_station, class_idx]
    if np.isnan(mu_rate) or mu_rate <= 0:
        raise RuntimeError(f"No valid service rate for class {class_idx} at queue")

    # Check stability
    rho = lambda_rate / mu_rate
    if rho >= 1.0:
        # Unstable system
        return {
            'mean_queue_length': float('inf'),
            'mean_response_time': float('inf'),
            'throughput': lambda_rate,
            'utilization': 1.0
        }

    # Extract process matrices
    arrival_params = _extract_arrival_process(sn, source_station, class_idx, lambda_rate)
    service_params = _extract_service_process(sn, queue_station, class_idx, mu_rate)

    # Check if simple M/M/1 case
    is_simple_exponential = arrival_params['is_simple'] and service_params['is_simple']

    if is_simple_exponential:
        # Use analytical M/M/1 formulas
        mean_L = rho / (1 - rho)
        mean_W = 1.0 / (mu_rate - lambda_rate)
        fl_moms = [mean_L]
        st_moms = [mean_W]
    else:
        # Use FluFluQueue for general case
        prec = max(options.tol, 1e-14)
        srv0stop = True  # Work-conserving behavior

        flu_result = FluFluQueue(
            Qin=arrival_params['Q'],
            Rin=arrival_params['R'],
            Qout=service_params['Q'],
            Rout=service_params['R'],
            srv0stop=srv0stop,
            numFluidMoments=2,
            numSojournMoments=2,
            prec=prec
        )

        if flu_result.fluidMoments is None:
            raise RuntimeError("FluFluQueue did not return fluid moments")
        if flu_result.sojournMoments is None:
            raise RuntimeError("FluFluQueue did not return sojourn moments")

        fl_moms = flu_result.fluidMoments
        st_moms = flu_result.sojournMoments

    # Map moments to LINE metrics
    mean_queue_length = fl_moms[0]
    mean_response_time = st_moms[0]

    # Throughput from Little's Law: X = L / W
    FINE_TOL = 1e-10
    if mean_response_time > FINE_TOL:
        throughput = mean_queue_length / mean_response_time
    else:
        throughput = lambda_rate

    # Utilization: U = lambda * E[S]
    utilization = min(1.0, lambda_rate / mu_rate)

    return {
        'mean_queue_length': mean_queue_length,
        'mean_response_time': mean_response_time,
        'throughput': throughput,
        'utilization': utilization
    }


def _extract_arrival_process(
    sn: NetworkStruct,
    station_idx: int,
    class_idx: int,
    lambda_rate: float
) -> dict:
    """
    Extracts arrival process parameters (Q, R matrices) from network structure.
    """
    # Check if proc data is available
    station = sn.stations[station_idx] if hasattr(sn, 'stations') else station_idx
    jobclass = sn.jobclasses[class_idx] if hasattr(sn, 'jobclasses') else class_idx

    proc = None
    if hasattr(sn, 'proc') and sn.proc is not None:
        if station in sn.proc and jobclass in sn.proc[station]:
            proc = sn.proc[station][jobclass]

    if proc is None or len(proc) == 0:
        # Simple Poisson arrival (1-state)
        return {
            'Q': np.zeros((1, 1)),
            'R': np.array([[lambda_rate]]),
            'is_simple': True
        }

    # Check number of phases
    n_phases = proc[0].shape[0]
    if n_phases == 1:
        # Single state: Poisson arrival
        return {
            'Q': np.zeros((1, 1)),
            'R': np.array([[lambda_rate]]),
            'is_simple': True
        }

    # Multi-phase: Extract D0, D1 and convert to (Q, R) format
    # D0 + D1 = Q (generator), sum of D1 row = R (fluid rate)
    D0 = np.asarray(proc[0])
    D1 = np.asarray(proc[1])

    Q = D0 + D1
    R = np.zeros((n_phases, n_phases))
    for i in range(n_phases):
        R[i, i] = np.sum(D1[i, :])

    return {'Q': Q, 'R': R, 'is_simple': False}


def _extract_service_process(
    sn: NetworkStruct,
    station_idx: int,
    class_idx: int,
    mu_rate: float
) -> dict:
    """
    Extracts service process parameters (Q, R matrices) from network structure.
    """
    # Check if proc data is available
    station = sn.stations[station_idx] if hasattr(sn, 'stations') else station_idx
    jobclass = sn.jobclasses[class_idx] if hasattr(sn, 'jobclasses') else class_idx

    proc = None
    if hasattr(sn, 'proc') and sn.proc is not None:
        if station in sn.proc and jobclass in sn.proc[station]:
            proc = sn.proc[station][jobclass]

    if proc is None or len(proc) == 0:
        # Simple exponential service (1-phase)
        return {
            'Q': np.zeros((1, 1)),
            'R': np.array([[mu_rate]]),
            'is_simple': True
        }

    # Check number of phases
    n_phases = proc[0].shape[0]
    if n_phases == 1:
        # Single phase: exponential service
        return {
            'Q': np.zeros((1, 1)),
            'R': np.array([[mu_rate]]),
            'is_simple': True
        }

    # Multi-phase: Extract D0, D1 and convert to (Q, R) format
    D0 = np.asarray(proc[0])
    D1 = np.asarray(proc[1])

    Q = D0 + D1
    R = np.zeros((n_phases, n_phases))
    for i in range(n_phases):
        R[i, i] = np.sum(D1[i, :])

    return {'Q': Q, 'R': R, 'is_simple': False}


__all__ = [
    'FLDResult',
    'solver_fld_analyzer',
    'matrix_method_analyzer',
    'closing_method_analyzer',
    'mfq_analyzer',
]
