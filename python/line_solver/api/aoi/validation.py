"""
Topology validation and parameter extraction for AoI solver.

This module provides validation functions to check if a network is compatible
with AoI analysis and utilities to extract parameters for the bufferless or
single-buffer AoI solvers.
"""

import numpy as np
from typing import Tuple, Dict, Any, Optional
from ...api.sn import NetworkStruct, NodeType, SchedStrategy


def aoi_is_aoi(sn: NetworkStruct) -> Tuple[bool, Dict[str, Any]]:
    """
    Validate if network topology is compatible with AoI analysis.

    Age of Information analysis is supported for single-queue systems with
    specific properties:
    - **Bufferless (capacity=1)**: PH/PH/1/1 or PH/PH/1/1* (preemptive)
    - **Single-buffer (capacity=2)**: M/PH/1/2 or M/PH/1/2* (replacement)

    Requirements:
    - Exactly one open job class (no closed chains)
    - Single queue station (source → queue → sink)
    - Single server (nservers=1)
    - Supported scheduling: FCFS, LCFS, LCFSPR
    - Capacity: 1 or 2 (validated by queue node)
    - If capacity=2: arrivals must be exponential (Poisson)

    Parameters
    ----------
    sn : NetworkStruct
        Network structure from `Network.getStruct()`

    Returns
    -------
    is_valid : bool
        True if network is valid AoI topology
    info : dict
        Validation info with keys:
        - 'is_valid': bool
        - 'systemType': 'bufferless' or 'singlebuffer' (if valid)
        - 'queue_idx': index of queue station (if valid)
        - 'capacity': queue capacity (if valid)
        - 'reason': error message (if not valid)
        - 'scheduling': scheduling strategy name
    """
    info = {'is_valid': False}

    # Check: Single open class
    if sn.nclasses != 1:
        info['reason'] = f"AoI requires single class, found {sn.nclasses}"
        return False, info

    # Check: Open class (no closed jobs)
    # Note: njobs array exists for all networks - open classes have njobs = inf or 0
    # closed classes have finite njobs > 0
    has_closed_jobs = sn.nclosedjobs > 0
    if not has_closed_jobs and hasattr(sn, 'njobs') and sn.njobs is not None:
        njobs_arr = np.asarray(sn.njobs)
        # Open classes have njobs = inf, closed classes have finite positive njobs
        has_closed_jobs = np.any(np.isfinite(njobs_arr) & (njobs_arr > 0))
    if has_closed_jobs:
        info['reason'] = "AoI requires open network (no closed classes)"
        return False, info

    # Check: Single service station (exclude Source/Sink from count)
    # Count only Queue/Delay nodes, not Source/Sink nodes
    # NodeType: 0=Source, 1=Sink, 2=Queue, 3=Delay, etc.
    service_stations = []
    for i, nt in enumerate(sn.nodetype):
        if nt in [2, 3]:  # Queue or Delay
            service_stations.append(i)

    if len(service_stations) != 1:
        info['reason'] = f"AoI requires single queue, found {len(service_stations)} service stations"
        return False, info

    queue_node_idx = service_stations[0]

    # Convert node index to station index
    # Try both snake_case and camelCase attribute names
    if hasattr(sn, 'nodeToStation'):
        queue_station_idx = sn.nodeToStation[queue_node_idx]
    elif hasattr(sn, 'node_to_station'):
        queue_station_idx = sn.node_to_station[queue_node_idx]
    else:
        queue_station_idx = 0

    # Fallback: if nodeToStation returned -1, try direct indexing
    if queue_station_idx < 0:
        # Assume first station that's not a source
        for i in range(sn.nstations):
            if i < len(sn.isstation) and sn.isstation[i] and sn.nodetype[i] != 0:
                queue_station_idx = i
                break
        if queue_station_idx < 0:
            queue_station_idx = 0

    # Check: Single server (must have at least nservers array)
    if sn.nservers is None or len(sn.nservers) == 0:
        # No server info, assume single server (default for M/M/1)
        pass
    elif queue_station_idx < len(sn.nservers) and sn.nservers[queue_station_idx] != 1:
        info['reason'] = f"AoI requires single server, found {sn.nservers[queue_station_idx]}"
        return False, info

    # Check: Queue capacity
    if not hasattr(sn, 'cap') and not hasattr(sn, 'capacity'):
        capacity = 1  # Default: bufferless
    else:
        cap_array = sn.cap if hasattr(sn, 'cap') else sn.capacity
        if cap_array is None:
            capacity = 1
        else:
            capacity = cap_array[queue_station_idx] if isinstance(cap_array, np.ndarray) and queue_station_idx < len(cap_array) else cap_array
            # Handle infinite capacity (unbounded)
            if np.isinf(capacity):
                capacity = 1  # Default to bufferless for inf capacity
            else:
                capacity = int(capacity) if capacity is not None else 1

    if capacity not in [1, 2]:
        info['reason'] = f"AoI supports capacity 1 (bufferless) or 2 (single-buffer), found {capacity}"
        return False, info

    # Determine system type
    system_type = 'bufferless' if capacity == 1 else 'singlebuffer'

    # Check: Scheduling strategy
    sched_id = sn.schedid[queue_station_idx] if hasattr(sn, 'schedid') and queue_station_idx < len(sn.schedid) else SchedStrategy.FCFS
    sched_name = _get_sched_name(sched_id)

    supported_scheds = {SchedStrategy.FCFS, SchedStrategy.LCFS, SchedStrategy.LCFSPR}
    if sched_id not in supported_scheds:
        info['reason'] = f"AoI supports FCFS/LCFS/LCFSPR scheduling, found {sched_name}"
        return False, info

    # For single-buffer: verify exponential arrivals
    if capacity == 2:
        if not _is_exponential_arrival(sn):
            info['reason'] = "AoI single-buffer requires exponential (Poisson) arrivals"
            return False, info

    # All checks passed
    info['is_valid'] = True
    info['systemType'] = system_type
    info['queue_idx'] = queue_station_idx  # Return station index, not node index
    info['queue_node_idx'] = queue_node_idx  # Also store node index for reference
    info['capacity'] = capacity
    info['scheduling'] = sched_name

    return True, info


def aoi_extract_params(sn: NetworkStruct, aoi_info: Dict[str, Any], options: Any) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """
    Extract parameters for AoI solver from network structure.

    For bufferless (capacity=1) systems:
        Returns: {tau, T, sigma, S, p}

    For single-buffer (capacity=2) systems:
        Returns: {lambda_rate, sigma, S, r}

    Parameters
    ----------
    sn : NetworkStruct
        Network structure
    aoi_info : dict
        Validation info from aoi_is_aoi() with 'systemType' and 'queue_idx'
    options : SolverFLDOptions
        Solver options (may contain aoi_preemption override)

    Returns
    -------
    aoi_params : dict
        Parameters for the selected AoI solver
    solver_config : dict
        Configuration metadata: {'systemType', 'queue_idx', 'preemption_type', ...}
    """
    from .conversion import aoi_dist2ph

    queue_idx = aoi_info['queue_idx']
    system_type = aoi_info['systemType']
    sched_id = sn.schedid[queue_idx] if hasattr(sn, 'schedid') else SchedStrategy.FCFS

    # Extract arrival process (if bufferless)
    if system_type == 'bufferless':
        # Get arrival process
        arrival_proc = None
        if sn.proc is not None:
            arrival_proc = sn.proc.get((-1, 0)) if isinstance(sn.proc, dict) else None

        if arrival_proc is None:
            # Fallback: use lambda_arr if available
            lambda_arr = sn.lambda_arr[0] if hasattr(sn, 'lambda_arr') and sn.lambda_arr is not None else 0.5
            # Create exponential arrival (default)
            arrival_proc = type('ExpDist', (), {'mean': 1.0 / lambda_arr})()

        # Convert to PH representation
        tau, T = aoi_dist2ph(arrival_proc)
    else:
        # Single-buffer: exponential arrival rate
        lambda_arr = sn.lambda_arr[0] if hasattr(sn, 'lambda_arr') and sn.lambda_arr is not None else 0.5
        tau = None  # Not used for M/PH/1/2
        T = None

    # Extract service process
    service_proc = None
    if sn.proc is not None:
        service_proc = sn.proc.get((queue_idx, 0)) if isinstance(sn.proc, dict) else None

    if service_proc is None:
        # Fallback: use rates if available
        mu = sn.rates[queue_idx, 0] if hasattr(sn, 'rates') and sn.rates is not None else 1.0
        # Create exponential service (default)
        service_proc = type('ExpDist', (), {'mean': 1.0 / mu})()

    # Convert service to PH representation
    sigma, S = aoi_dist2ph(service_proc)

    # Determine preemption/replacement probability
    if system_type == 'bufferless':
        # p: preemption probability (0=FCFS, 1=LCFSPR)
        if hasattr(options, 'aoi_preemption') and options.aoi_preemption is not None:
            p = options.aoi_preemption
        else:
            # Auto-detect from scheduling strategy
            p = 1.0 if sched_id == SchedStrategy.LCFSPR else 0.0

        aoi_params = {
            'tau': tau,
            'T': T,
            'sigma': sigma,
            'S': S,
            'p': p,
        }
        solver_config = {
            'systemType': 'bufferless',
            'queue_idx': queue_idx,
            'preemption_type': 'LCFSPR' if p > 0.5 else 'FCFS',
        }
    else:
        # Single-buffer: M/PH/1/2
        # r: replacement probability (0=FCFS, 1=replacement)
        if hasattr(options, 'aoi_preemption') and options.aoi_preemption is not None:
            r = options.aoi_preemption
        else:
            # Auto-detect: LCFSPR → replacement
            r = 1.0 if sched_id == SchedStrategy.LCFSPR else 0.0

        aoi_params = {
            'lambda_rate': lambda_arr,
            'sigma': sigma,
            'S': S,
            'r': r,
        }
        solver_config = {
            'systemType': 'singlebuffer',
            'queue_idx': queue_idx,
            'replacement_type': 'replacement' if r > 0.5 else 'FCFS',
        }

    return aoi_params, solver_config


def _get_sched_name(sched_id: int) -> str:
    """Get human-readable scheduling strategy name."""
    sched_names = {
        SchedStrategy.FCFS: 'FCFS',
        SchedStrategy.LCFS: 'LCFS',
        SchedStrategy.PS: 'PS',
        SchedStrategy.DPS: 'DPS',
        SchedStrategy.SIRO: 'SIRO',
        SchedStrategy.INF: 'INF',
        SchedStrategy.LCFSPR: 'LCFSPR',
    }
    return sched_names.get(sched_id, f'Unknown({sched_id})')


def _is_exponential_arrival(sn: NetworkStruct) -> bool:
    """
    Check if arrival process is exponential (Poisson).

    For single-buffer systems, only exponential arrivals are supported.
    """
    # Try to get arrival process
    arrival_proc = sn.proc.get((-1, 0))  # Source → Queue (class 0)
    if arrival_proc is None:
        # No process defined, assume exponential
        return True

    # Check if it's exponential
    proc_name = type(arrival_proc).__name__
    return proc_name in ['Exp', 'Exponential'] or (hasattr(arrival_proc, 'cv') and abs(arrival_proc.cv - 1.0) < 0.01)
