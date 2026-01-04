"""
Network adapter for MAM solvers.

Extracts MAM-specific parameters from NetworkStruct and provides
utilities for converting network specifications to MAM-friendly formats.

This adapter handles:
- Station type identification (Delay, Queue, Fork, Join)
- Service distribution extraction
- Visit count computation
- Routing matrix handling
- Chain identification
"""

import numpy as np
from typing import List, Dict, Tuple, Optional
from ....api.sn.network_struct import NetworkStruct, NodeType, SchedStrategy


def extract_mam_params(sn: NetworkStruct) -> Dict:
    """
    Extract MAM-specific parameters from NetworkStruct.

    Args:
        sn: Compiled NetworkStruct from Network model

    Returns:
        Dictionary containing:
            - stations: List of station indices (excluding Sources/Sinks)
            - nstations: Number of stations
            - nclasses: Number of job classes
            - rates: (M, K) service rates matrix
            - scv: (M, K) service SCV matrix
            - visits: (M, K) visit counts matrix
            - routing: Routing matrix information
            - sched: (M,) scheduling strategies
            - node_types: List of NodeType for each station
            - is_delay: (M,) boolean mask for delay nodes
            - is_queue: (M,) boolean mask for queue nodes
            - is_fork: (M,) boolean mask for fork nodes
            - is_join: (M,) boolean mask for join nodes
    """
    M = sn.nstations
    K = sn.nclasses

    # Extract basic dimensions
    params = {
        'nstations': M,
        'nclasses': K,
        'nchains': sn.nchains,
        'nclosedjobs': sn.nclosedjobs,
    }

    # Extract service parameters
    rates = np.asarray(sn.rates, dtype=np.float64)
    scv = np.asarray(sn.scv, dtype=np.float64)

    # Ensure correct shapes
    if rates.ndim == 1:
        rates = rates.reshape(-1, 1)
    if scv.ndim == 1:
        scv = scv.reshape(-1, 1)

    params['rates'] = rates
    params['scv'] = scv

    # Extract visit counts
    try:
        if hasattr(sn, 'visits') and sn.visits is not None:
            visits = np.asarray(sn.visits, dtype=np.float64)
            if visits.ndim == 1:
                visits = visits.reshape(-1, 1)
            if visits.shape[0] != M or visits.shape[1] != K:
                visits = np.ones((M, K), dtype=np.float64)
        else:
            # Fallback: use unit visits
            visits = np.ones((M, K), dtype=np.float64)
    except:
        visits = np.ones((M, K), dtype=np.float64)

    params['visits'] = visits

    # Extract number of servers
    nservers = np.asarray(sn.nservers, dtype=np.float64).flatten()
    if len(nservers) < M:
        nservers = np.pad(nservers, (0, M - len(nservers)), constant_values=1)
    params['nservers'] = nservers[:M]

    # Identify station types
    node_types = []
    is_delay = np.zeros(M, dtype=bool)
    is_queue = np.zeros(M, dtype=bool)
    is_fork = np.zeros(M, dtype=bool)
    is_join = np.zeros(M, dtype=bool)

    for m in range(M):
        # Find the node index for this station
        node_idx = sn.stationToNode[m] if len(sn.stationToNode) > m else m
        if node_idx >= len(sn.nodetype):
            node_type = NodeType.QUEUE
        else:
            node_type = sn.nodetype[node_idx]

        node_types.append(node_type)

        if node_type == NodeType.DELAY:
            is_delay[m] = True
        elif node_type == NodeType.QUEUE:
            is_queue[m] = True
        elif node_type == NodeType.FORK:
            is_fork[m] = True
        elif node_type == NodeType.JOIN:
            is_join[m] = True

    params['node_types'] = node_types
    params['is_delay'] = is_delay
    params['is_queue'] = is_queue
    params['is_fork'] = is_fork
    params['is_join'] = is_join

    # Extract scheduling strategies
    sched = np.zeros(M, dtype=int)
    for m in range(M):
        node_idx = sn.stationToNode[m] if len(sn.stationToNode) > m else m
        if node_idx in sn.sched:
            strategy = sn.sched[node_idx]
            # Convert SchedStrategy enum to int if needed
            sched[m] = strategy.value if hasattr(strategy, 'value') else int(strategy)
        else:
            # Default: FCFS
            sched[m] = SchedStrategy.FCFS.value if hasattr(SchedStrategy.FCFS, 'value') else int(SchedStrategy.FCFS)

    params['sched'] = sched

    # Extract routing information
    if sn.rt is not None:
        params['routing'] = np.asarray(sn.rt, dtype=np.float64)
    else:
        # Fallback: identity routing
        params['routing'] = np.eye(M, dtype=np.float64)

    # Extract population (for closed networks)
    if sn.njobs is not None:
        params['njobs'] = np.asarray(sn.njobs, dtype=np.float64).flatten()
    else:
        params['njobs'] = np.inf * np.ones(K, dtype=np.float64)

    # Extract node names
    params['nodenames'] = sn.nodenames if sn.nodenames else [f'Node_{m}' for m in range(M)]
    params['classnames'] = sn.classnames if sn.classnames else [f'Class_{k}' for k in range(K)]

    # Extract reference station and class for each chain
    if len(sn.refstat) > 0:
        params['refstat'] = np.asarray(sn.refstat, dtype=int).flatten()
    else:
        params['refstat'] = np.zeros(K, dtype=int)

    return params


def identify_station_types(sn: NetworkStruct) -> Dict[str, List[int]]:
    """
    Identify and classify stations by type.

    Args:
        sn: Compiled NetworkStruct

    Returns:
        Dictionary with keys:
            - 'delays': Indices of delay nodes
            - 'queues': Indices of queue nodes
            - 'forks': Indices of fork nodes
            - 'joins': Indices of join nodes
    """
    classification = {
        'delays': [],
        'queues': [],
        'forks': [],
        'joins': [],
        'sources': [],
        'sinks': [],
    }

    M = sn.nstations
    for m in range(M):
        node_idx = sn.stationToNode[m] if len(sn.stationToNode) > m else m
        if node_idx >= len(sn.nodetype):
            node_type = NodeType.QUEUE
        else:
            node_type = sn.nodetype[node_idx]

        if node_type == NodeType.DELAY:
            classification['delays'].append(m)
        elif node_type == NodeType.QUEUE:
            classification['queues'].append(m)
        elif node_type == NodeType.FORK:
            classification['forks'].append(m)
        elif node_type == NodeType.JOIN:
            classification['joins'].append(m)
        elif node_type == NodeType.SOURCE:
            classification['sources'].append(m)
        elif node_type == NodeType.SINK:
            classification['sinks'].append(m)

    return classification


def build_routing_matrix(sn: NetworkStruct) -> np.ndarray:
    """
    Extract or construct routing matrix from NetworkStruct.

    Args:
        sn: Compiled NetworkStruct

    Returns:
        (M, M) routing matrix where element [i,j] is the probability
        that a job departs from station i to station j
    """
    M = sn.nstations

    if sn.rt is not None:
        # Use existing routing matrix, handling possible multi-dimensional case
        rt = np.asarray(sn.rt, dtype=np.float64)
        if rt.ndim == 3:
            # Multi-class routing: collapse to aggregate
            rt = np.mean(rt, axis=2)
        elif rt.ndim == 1:
            # Flatten to 2D
            side = int(np.sqrt(len(rt)))
            if side * side == len(rt):
                rt = rt.reshape(side, side)
        return rt.astype(np.float64)
    else:
        # Default: self-loop or connection matrix based approach
        if sn.connmatrix is not None:
            # Use connection matrix as basis
            conn = np.asarray(sn.connmatrix, dtype=np.float64)
            # Normalize rows to get routing probabilities
            rowsums = conn.sum(axis=1, keepdims=True)
            rowsums[rowsums == 0] = 1
            rt = conn / rowsums
            return rt
        else:
            # Fallback: identity (self-loops only)
            return np.eye(M, dtype=np.float64)


def get_service_distribution(sn: NetworkStruct, station_idx: int, class_idx: int) -> Dict:
    """
    Get service distribution parameters for a station-class pair.

    Args:
        sn: Compiled NetworkStruct
        station_idx: Station index (0-based)
        class_idx: Class index (0-based)

    Returns:
        Dictionary with keys:
            - 'rate': Service rate (lambda = 1/mean)
            - 'mean': Mean service time
            - 'scv': Squared coefficient of variation
            - 'type': Distribution type ('M' for exponential, 'G' for general, etc.)
    """
    M = sn.nstations
    K = sn.nclasses

    if station_idx >= M or class_idx >= K:
        raise ValueError(f"Invalid station/class indices: ({station_idx}, {class_idx})")

    # Extract rate and SCV
    rates = np.asarray(sn.rates, dtype=np.float64)
    scv_vals = np.asarray(sn.scv, dtype=np.float64)

    # Handle different array shapes
    if rates.ndim == 1:
        if station_idx < len(rates):
            rate = rates[station_idx]
        else:
            rate = 1.0
    elif rates.ndim >= 2:
        if rates.shape[0] > station_idx and rates.shape[1] > class_idx:
            rate = rates[station_idx, class_idx]
        elif rates.shape[0] > station_idx:
            rate = rates[station_idx, 0]
        else:
            rate = 1.0
    else:
        rate = 1.0

    if scv_vals.ndim == 1:
        if station_idx < len(scv_vals):
            scv = scv_vals[station_idx]
        else:
            scv = 1.0
    elif scv_vals.ndim >= 2:
        if scv_vals.shape[0] > station_idx and scv_vals.shape[1] > class_idx:
            scv = scv_vals[station_idx, class_idx]
        elif scv_vals.shape[0] > station_idx:
            scv = scv_vals[station_idx, 0]
        else:
            scv = 1.0
    else:
        scv = 1.0

    mean = 1.0 / max(rate, 1e-10)

    # Classify distribution type
    if abs(scv - 1.0) < 1e-6:
        dist_type = 'M'  # Exponential (Markovian)
    elif scv < 1.0:
        dist_type = 'E'  # Erlang-like (SCV < 1)
    elif scv > 1.0:
        dist_type = 'H'  # Hyperexponential-like (SCV > 1)
    else:
        dist_type = 'G'  # General

    return {
        'rate': float(rate),
        'mean': float(mean),
        'scv': float(scv),
        'type': dist_type,
    }


def extract_visit_counts(sn: NetworkStruct) -> np.ndarray:
    """
    Extract visit counts matrix from NetworkStruct.

    Args:
        sn: Compiled NetworkStruct

    Returns:
        (M, K) visit counts matrix
    """
    M = sn.nstations
    K = sn.nclasses

    # Try to use visits from network
    try:
        if hasattr(sn, 'visits') and sn.visits is not None:
            visits = np.asarray(sn.visits, dtype=np.float64)
            if visits.ndim == 1:
                visits = visits.reshape(-1, 1)
            if visits.shape[0] < M:
                # Pad with unit visits
                visits = np.pad(visits, ((0, M - visits.shape[0]), (0, 0)), constant_values=1)
            if visits.shape[1] < K:
                # Pad with unit visits
                visits = np.pad(visits, ((0, 0), (0, K - visits.shape[1])), constant_values=1)
            return visits[:M, :K]
        else:
            # Fallback: unit visits
            return np.ones((M, K), dtype=np.float64)
    except:
        # Fallback: unit visits on any error
        return np.ones((M, K), dtype=np.float64)


def check_singleclass_network(sn: NetworkStruct) -> bool:
    """
    Check if network has only a single job class.

    Args:
        sn: Compiled NetworkStruct

    Returns:
        True if network has exactly one job class
    """
    return sn.nclasses == 1


def check_closed_network(sn: NetworkStruct) -> bool:
    """
    Check if network is closed (no sources or sinks with infinite arrivals).

    Args:
        sn: Compiled NetworkStruct

    Returns:
        True if network is closed
    """
    if sn.njobs is None or len(sn.njobs) == 0:
        return False
    return all(np.isfinite(n) and n > 0 for n in sn.njobs)


def check_product_form(sn: NetworkStruct) -> Tuple[bool, Optional[str]]:
    """
    Check if network satisfies product-form assumptions.

    Args:
        sn: Compiled NetworkStruct

    Returns:
        Tuple of (is_product_form, reason_if_not)
    """
    # Product-form networks require:
    # - FCFS queues with exponential service
    # - Delay nodes
    # - Probabilistic routing
    # - No priorities or special scheduling

    for m in range(sn.nstations):
        node_idx = sn.stationToNode[m] if len(sn.stationToNode) > m else m
        if node_idx in sn.sched:
            sched = sn.sched[node_idx]
            if sched not in [SchedStrategy.FCFS, SchedStrategy.INF]:
                return False, f"Non-product-form scheduling at station {m}: {sched}"

    # Check for exponential service (SCV ≈ 1)
    scv = np.asarray(sn.scv, dtype=np.float64)
    if np.any(np.abs(scv - 1.0) > 0.1):  # Allow small tolerance
        return False, "Non-exponential service distributions"

    return True, None


def is_fork_join_network(sn: NetworkStruct) -> bool:
    """
    Check if network has Fork-Join topology.

    Args:
        sn: Compiled NetworkStruct

    Returns:
        True if network appears to have Fork-Join structure
    """
    # Simple check: presence of at least one Fork and one Join
    has_fork = any(nt == NodeType.FORK for nt in sn.nodetype)
    has_join = any(nt == NodeType.JOIN for nt in sn.nodetype)
    return has_fork and has_join
