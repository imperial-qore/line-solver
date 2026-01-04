"""
CTMC Solver handler.

Native Python implementation of CTMC (Continuous-Time Markov Chain) solver
handler that analyzes queueing networks through exact state-space enumeration.

The CTMC solver builds the complete state space and infinitesimal generator
matrix, then solves for steady-state probabilities to compute performance metrics.

Port from:

"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, Dict, List, Tuple, Any
from itertools import product
import time

import warnings
from ...sn import (
    NetworkStruct,
    SchedStrategy,
    NodeType,
    sn_is_open_model,
    sn_is_closed_model,
    sn_has_open_classes,
)
from ...mc import ctmc_solve, ctmc_makeinfgen


@dataclass
class SolverCTMCOptions:
    """Options for CTMC solver."""
    method: str = 'default'
    tol: float = 1e-6
    verbose: bool = False
    cutoff: int = 10  # Cutoff for open class populations
    hide_immediate: bool = True  # Hide immediate transitions
    state_space_gen: str = 'default'  # 'default', 'full', 'reachable'


@dataclass
class SolverCTMCReturn:
    """
    Result of CTMC solver handler.

    Attributes:
        Q: Mean queue lengths (M x K)
        U: Utilizations (M x K)
        R: Response times (M x K)
        T: Throughputs (M x K)
        C: Cycle times (1 x K)
        X: System throughputs (1 x K)
        pi: Steady-state distribution
        infgen: Infinitesimal generator matrix
        space: State space matrix
        runtime: Runtime in seconds
        method: Method used
    """
    Q: Optional[np.ndarray] = None
    U: Optional[np.ndarray] = None
    R: Optional[np.ndarray] = None
    T: Optional[np.ndarray] = None
    C: Optional[np.ndarray] = None
    X: Optional[np.ndarray] = None
    pi: Optional[np.ndarray] = None
    infgen: Optional[np.ndarray] = None
    space: Optional[np.ndarray] = None
    runtime: float = 0.0
    method: str = "default"


def _enumerate_state_space(
    sn: NetworkStruct,
    cutoff = 10
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Enumerate the state space for a queueing network.

    For closed networks, enumerates all valid job distributions.
    For open networks, uses cutoff to bound the state space.

    Args:
        sn: Network structure
        cutoff: Maximum jobs per station for open networks.
                Can be an int (same cutoff for all), or a matrix (M x K)
                with per-station, per-class cutoffs.

    Returns:
        Tuple of (state_space, state_space_aggr)
    """
    M = sn.nstations
    K = sn.nclasses
    # Use sn_has_open_classes instead of sn_is_open_model to properly handle
    # mixed networks (which have both open and closed classes)
    is_open = sn_has_open_classes(sn)

    # Get population constraints
    if sn.njobs is not None:
        N = sn.njobs.flatten()
    else:
        N = np.ones(K)

    # Check if model contains Cache or Router nodes - these require reduced cutoff
    # This matches MATLAB's behavior where cache_replc_routing uses cutoff=1
    has_cache_or_router = False
    if hasattr(sn, 'nodetype') and sn.nodetype is not None:
        for nt in sn.nodetype:
            nt_val = int(nt.value) if hasattr(nt, 'value') else int(nt)
            if nt_val == 6 or nt_val == 7:  # CACHE or ROUTER
                has_cache_or_router = True
                break

    # For models with Cache/Router nodes, limit cutoff to 1 to prevent state explosion
    # This matches MATLAB's approach in cache_replc_routing.m which uses cutoff=1
    if has_cache_or_router and (np.isscalar(cutoff) or np.atleast_2d(cutoff).size == 1):
        effective_cutoff = min(cutoff, 1) if np.isscalar(cutoff) else min(np.atleast_2d(cutoff).flat[0], 1)
        cutoff = effective_cutoff

    # Handle cutoff as scalar or matrix
    cutoff_arr = np.atleast_2d(cutoff)
    is_matrix_cutoff = cutoff_arr.shape[0] > 1 or cutoff_arr.shape[1] > 1

    # Identify Cache and Router stations - these have immediate processing (capacity 1)
    # This matches MATLAB's spaceGeneratorNodes.m behavior
    # Note: We need to map station indices to node indices since nodetype is indexed by node
    cache_router_stations = set()
    if hasattr(sn, 'nodetype') and sn.nodetype is not None:
        # Get station to node mapping if available
        station_to_node = None
        if hasattr(sn, 'stationToNode') and sn.stationToNode is not None:
            station_to_node = sn.stationToNode

        for ist in range(M):
            # Get the node index for this station
            if station_to_node is not None and ist < len(station_to_node):
                node_idx = int(station_to_node[ist])
            else:
                node_idx = ist  # Fallback: assume station index = node index

            if node_idx >= 0 and node_idx < len(sn.nodetype):
                nt = sn.nodetype[node_idx]
                # NodeType.CACHE = 6, NodeType.ROUTER = 7
                nt_val = int(nt.value) if hasattr(nt, 'value') else int(nt)
                if nt_val == 6 or nt_val == 7:  # CACHE or ROUTER
                    cache_router_stations.add(ist)

    def get_cutoff(ist: int, k: int) -> int:
        """Get cutoff for station ist and class k."""
        # Cache and Router nodes have capacity 1 (immediate processing)
        # This prevents state space explosion in models with these nodes
        if ist in cache_router_stations:
            return 1
        if is_matrix_cutoff:
            # Matrix cutoff: index by station and class
            # Handle different matrix layouts
            if cutoff_arr.shape[0] >= M and cutoff_arr.shape[1] >= K:
                return int(cutoff_arr[ist, k])
            elif cutoff_arr.shape[0] >= K and cutoff_arr.shape[1] >= M:
                # Transposed: (K x M) layout
                return int(cutoff_arr[k, ist])
            else:
                # Fallback: use first valid value or default
                return int(cutoff_arr.flat[0]) if cutoff_arr.size > 0 else 10
        else:
            return int(cutoff_arr.flat[0])

    states = []

    # Helper function for closed class distributions
    def _enumerate_distributions(n_jobs: int, n_stations: int) -> List[List[int]]:
        """Enumerate all ways to distribute n_jobs among n_stations."""
        if n_stations == 1:
            return [[n_jobs]]
        if n_jobs == 0:
            return [[0] * n_stations]

        result = []
        for i in range(n_jobs + 1):
            for rest in _enumerate_distributions(n_jobs - i, n_stations - 1):
                result.append([i] + rest)
        return result

    # Determine which classes are open vs closed
    open_classes = [k for k in range(K) if np.isinf(N[k])]
    closed_classes = [k for k in range(K) if np.isfinite(N[k]) and N[k] > 0]

    # Identify Source and Sink stations (closed class jobs should NOT be there)
    source_station = -1
    sink_station = -1
    if hasattr(sn, 'sched') and sn.sched is not None:
        for ist, sched_val in sn.sched.items():
            if sched_val == SchedStrategy.EXT or (isinstance(sched_val, int) and sched_val == 16):
                source_station = ist
    # Also check nodetype for Sink
    if hasattr(sn, 'nodetype') and sn.nodetype is not None:
        for ist in range(M):
            if ist < len(sn.nodetype) and sn.nodetype[ist] == NodeType.SINK:
                sink_station = ist

    # Stations where closed class jobs can reside (exclude Source/Sink)
    closed_valid_stations = [ist for ist in range(M) if ist != source_station and ist != sink_station]
    n_closed_stations = len(closed_valid_stations)

    if is_open and closed_classes:
        # Mixed network: combine closed class distributions with open class enumeration
        # For closed classes: enumerate valid distributions (conserving population)
        #                     only among valid (non-Source/Sink) stations
        # For open classes: enumerate 0 to cutoff at each station (except Source for arrivals)

        # Generate closed class distributions among valid stations only
        closed_class_dists = []
        for k in closed_classes:
            n_k = int(N[k])
            # Distribute among valid stations only
            class_dists = _enumerate_distributions(n_k, n_closed_stations)
            closed_class_dists.append((k, class_dists))

        # Generate open class ranges (per non-Source station for queue occupancy)
        # Open class jobs can be at any station except Source (which generates arrivals)
        open_class_ranges = []
        for k in open_classes:
            station_ranges = []
            for ist in range(M):
                if ist == source_station:
                    # Source doesn't hold jobs - only range is 0
                    station_ranges.append(range(1))  # Just [0]
                else:
                    c = get_cutoff(ist, k)
                    station_ranges.append(range(c + 1))
            open_class_ranges.append((k, station_ranges))

        # Combine: iterate over closed distributions × open combinations
        from itertools import product as iter_product

        # Get all combinations of closed class distributions
        if closed_class_dists:
            closed_combos = list(iter_product(*[dists for _, dists in closed_class_dists]))
        else:
            closed_combos = [()]

        # Get all combinations of open class values at each station
        if open_class_ranges:
            # Flatten: for each open class, product over stations
            open_station_products = []
            for k, station_ranges in open_class_ranges:
                open_station_products.append(list(iter_product(*station_ranges)))
            open_combos = list(iter_product(*open_station_products))
        else:
            open_combos = [()]

        # Build states from combinations
        for closed_combo in closed_combos:
            for open_combo in open_combos:
                # Build state vector: [n_11, n_12, ..., n_1K, n_21, ..., n_MK]
                state = [0] * (M * K)
                # Fill in closed classes (mapping from valid station indices to full indices)
                for idx, (k, _) in enumerate(closed_class_dists):
                    dist = closed_combo[idx]
                    for valid_idx, ist in enumerate(closed_valid_stations):
                        state[ist * K + k] = dist[valid_idx]
                # Fill in open classes
                for idx, (k, _) in enumerate(open_class_ranges):
                    station_vals = open_combo[idx]
                    for ist in range(M):
                        state[ist * K + k] = station_vals[ist]
                states.append(state)

    elif is_open:
        # Pure open network: enumerate all (station, class) combinations up to cutoff
        # State has M*K dimensions: jobs at each (station, class) pair
        ranges = []
        for ist in range(M):
            for k in range(K):
                c = get_cutoff(ist, k)
                ranges.append(range(c + 1))

        for state in product(*ranges):
            states.append(list(state))
    else:
        # Pure closed network: enumerate valid distributions
        # For multiclass, we need to combine distributions
        all_class_dists = []
        for k in range(K):
            n_k = int(N[k]) if np.isfinite(N[k]) else 0
            class_dists = _enumerate_distributions(n_k, M)
            all_class_dists.append(class_dists)

        # Combine all class distributions
        for combo in product(*all_class_dists):
            # combo is tuple of distributions, one per class
            # Convert to state vector: [n_11, n_12, ..., n_1K, n_21, ..., n_MK]
            state = []
            for ist in range(M):
                for k in range(K):
                    state.append(combo[k][ist])
            states.append(state)

    if not states:
        states = [[0] * (M * K)]

    state_space = np.array(states, dtype=np.float64)

    # Aggregated state space: sum over classes at each station
    state_space_aggr = np.zeros((len(states), M))
    for i, state in enumerate(states):
        for ist in range(M):
            for k in range(K):
                state_space_aggr[i, ist] += state[ist * K + k]

    return state_space, state_space_aggr


def _build_state_index_map(state_space: np.ndarray) -> dict:
    """
    Build a hash map from state tuples to indices for O(1) lookup.

    Args:
        state_space: Enumerated state space (n_states x state_dim)

    Returns:
        Dictionary mapping state tuples to their indices
    """
    state_map = {}
    for i, state in enumerate(state_space):
        # Convert to tuple of integers for hashing
        state_key = tuple(int(x) for x in state)
        state_map[state_key] = i
    return state_map


def _find_state_index_fast(state_map: dict, state: np.ndarray) -> int:
    """
    Find index of state using hash map (O(1) lookup).

    Args:
        state_map: Hash map from state tuples to indices
        state: State vector to find

    Returns:
        Index of state, or -1 if not found
    """
    state_key = tuple(int(x) for x in state)
    return state_map.get(state_key, -1)


def _build_generator(
    sn: NetworkStruct,
    state_space: np.ndarray,
    options: SolverCTMCOptions
) -> np.ndarray:
    """
    Build the infinitesimal generator matrix for the queueing network.

    Args:
        sn: Network structure
        state_space: Enumerated state space
        options: Solver options

    Returns:
        Infinitesimal generator matrix Q
    """
    M = sn.nstations
    K = sn.nclasses
    n_states = state_space.shape[0]

    Q = np.zeros((n_states, n_states))

    # Handle cutoff as scalar or matrix for state validation
    cutoff = options.cutoff
    cutoff_arr = np.atleast_2d(cutoff)
    is_matrix_cutoff = cutoff_arr.shape[0] > 1 or cutoff_arr.shape[1] > 1

    def is_within_cutoff(new_state: np.ndarray) -> bool:
        """Check if state is within cutoff bounds."""
        if is_matrix_cutoff:
            # Reshape state to (M, K) and compare with cutoff matrix
            state_matrix = new_state.reshape(M, K)
            # Handle different cutoff matrix layouts
            if cutoff_arr.shape[0] >= M and cutoff_arr.shape[1] >= K:
                return np.all(state_matrix <= cutoff_arr[:M, :K])
            elif cutoff_arr.shape[0] >= K and cutoff_arr.shape[1] >= M:
                return np.all(state_matrix <= cutoff_arr[:K, :M].T)
            else:
                # Fallback: use scalar comparison
                return np.all(new_state <= cutoff_arr.flat[0])
        else:
            return np.all(new_state <= cutoff_arr.flat[0])

    # Build hash map for O(1) state lookup instead of O(n) linear search
    state_map = _build_state_index_map(state_space)

    # Get service rates
    if hasattr(sn, 'rates') and sn.rates is not None:
        rates = np.asarray(sn.rates)
    else:
        rates = np.ones((M, K))

    # Get routing probabilities
    if hasattr(sn, 'rt') and sn.rt is not None:
        P = np.asarray(sn.rt)
    else:
        # Default: uniform routing
        P = np.ones((M * K, M * K)) / (M * K)

    # Get number of servers
    if hasattr(sn, 'nservers') and sn.nservers is not None:
        nservers = np.asarray(sn.nservers).flatten()
    else:
        nservers = np.ones(M)

    # Get load-dependent scaling
    lldscaling = None
    if hasattr(sn, 'lldscaling') and sn.lldscaling is not None:
        lldscaling = np.asarray(sn.lldscaling)

    # Track which stations are infinite servers
    inf_server_stations = set()
    if hasattr(sn, 'sched') and sn.sched is not None:
        for ist, sched_val in sn.sched.items():
            if sched_val == SchedStrategy.INF:
                inf_server_stations.add(ist)

    def get_load_scaling(ist: int, total_jobs: int) -> float:
        """Get service rate scaling factor for station ist with total_jobs present."""
        if total_jobs <= 0:
            return 0.0  # No jobs, no service
        # For infinite servers, each job gets its own server
        if ist in inf_server_stations:
            return float(total_jobs)
        if lldscaling is not None and ist < lldscaling.shape[0]:
            # lldscaling[ist, n-1] gives scaling when n jobs present
            idx = total_jobs - 1
            if lldscaling.ndim > 1 and idx < lldscaling.shape[1]:
                return lldscaling[ist, idx]
            elif lldscaling.ndim > 1:
                return lldscaling[ist, -1]  # Use last value if beyond range
        # Fall back to min(n, c) behavior
        c = nservers[ist] if ist < len(nservers) else 1
        return min(total_jobs, c)

    # Identify Source station (EXT scheduling) for open network arrivals
    # Use sn_has_open_classes to properly handle mixed networks
    is_open = sn_has_open_classes(sn)
    source_station = -1
    if is_open and hasattr(sn, 'sched'):
        for ist, sched_val in sn.sched.items():
            # Check for EXT scheduling (Source station)
            if sched_val == SchedStrategy.EXT or (isinstance(sched_val, int) and sched_val == 11):
                source_station = ist
                break

    # Build transitions
    for s, state in enumerate(state_space):
        # Reshape state to (M, K)
        state_matrix = state.reshape(M, K)

        # === External arrivals from Source ===
        # For open networks, Source generates arrivals that enter the next station
        if source_station >= 0:
            for k in range(K):
                # Arrival rate at source
                arrival_rate = rates[source_station, k] if source_station < rates.shape[0] and k < rates.shape[1] else 0

                if arrival_rate > 0:
                    # Find routing from source to next station
                    src_idx = source_station * K + k
                    for jst in range(M):
                        if jst == source_station:
                            continue  # Skip self-routing for source
                        for r in range(K):
                            dst_idx = jst * K + r
                            if src_idx < P.shape[0] and dst_idx < P.shape[1]:
                                p = P[src_idx, dst_idx]
                            else:
                                p = 0

                            if p > 0:
                                # Create new state with arrival at destination
                                new_state = state.copy()
                                new_state[jst * K + r] += 1  # Arrival to jst, class r

                                # Check if within cutoff
                                if is_within_cutoff(new_state):
                                    ns = _find_state_index_fast(state_map, new_state)
                                    if ns >= 0:
                                        Q[s, ns] += arrival_rate * p

        # === Service completions and internal routing ===
        for ist in range(M):
            # Skip source station - it doesn't hold jobs for service
            if ist == source_station:
                continue

            for k in range(K):
                n_ik = state_matrix[ist, k]
                if n_ik <= 0:
                    continue

                # Service completion at station ist, class k
                # Rate depends on scheduling discipline
                if hasattr(sn, 'sched') and sn.sched is not None:
                    sched = sn.sched.get(ist, SchedStrategy.FCFS)
                else:
                    sched = SchedStrategy.FCFS

                mu = rates[ist, k] if ist < rates.shape[0] and k < rates.shape[1] else 1.0

                if sched == SchedStrategy.INF:
                    # Infinite server: rate = n * mu
                    rate = n_ik * mu
                elif sched == SchedStrategy.DPS:
                    # Discriminatory Processor Sharing: weighted by schedparam
                    total_jobs = np.sum(state_matrix[ist, :])
                    if total_jobs > 0:
                        # Get weight for this class (default 1.0)
                        w_k = 1.0
                        if hasattr(sn, 'schedparam') and sn.schedparam is not None:
                            if ist < sn.schedparam.shape[0] and k < sn.schedparam.shape[1]:
                                w_k = sn.schedparam[ist, k]
                        # Compute weighted share: sum of w_j * n_j for all classes
                        weighted_total = 0.0
                        for kk in range(K):
                            n_kk = state_matrix[ist, kk]
                            w_kk = 1.0
                            if hasattr(sn, 'schedparam') and sn.schedparam is not None:
                                if ist < sn.schedparam.shape[0] and kk < sn.schedparam.shape[1]:
                                    w_kk = sn.schedparam[ist, kk]
                            weighted_total += w_kk * n_kk
                        # Rate = (w_k * n_k / weighted_total) * scaling * mu
                        scaling = get_load_scaling(ist, int(total_jobs))
                        if weighted_total > 0:
                            rate = (w_k * n_ik / weighted_total) * scaling * mu
                        else:
                            rate = 0
                    else:
                        rate = 0
                else:
                    # FCFS, PS, etc.: rate = scaling * mu
                    total_jobs = np.sum(state_matrix[ist, :])
                    scaling = get_load_scaling(ist, int(total_jobs))
                    rate = scaling * mu * (n_ik / total_jobs) if total_jobs > 0 else 0

                if rate <= 0:
                    continue

                # Calculate routing probabilities
                src_idx = ist * K + k
                total_routing = 0.0

                # Find destination states based on routing
                for jst in range(M):
                    if jst == source_station:
                        continue  # Can't route to source
                    for r in range(K):
                        # Routing probability from (ist, k) to (jst, r)
                        dst_idx = jst * K + r

                        if src_idx < P.shape[0] and dst_idx < P.shape[1]:
                            p = P[src_idx, dst_idx]
                        else:
                            p = 0

                        if p <= 0:
                            continue

                        total_routing += p

                        # Create new state
                        new_state = state.copy()
                        new_state[ist * K + k] -= 1  # Departure from ist, class k
                        new_state[jst * K + r] += 1  # Arrival to jst, class r

                        # Check if new state is valid
                        if np.any(new_state < 0):
                            continue

                        # Find new state index using O(1) hash lookup
                        ns = _find_state_index_fast(state_map, new_state)
                        if ns >= 0:
                            Q[s, ns] += rate * p

                # === Departures to sink (for open networks) ===
                # If routing doesn't sum to 1, remaining probability exits to sink
                if is_open and total_routing < 1.0:
                    exit_prob = 1.0 - total_routing
                    if exit_prob > 0:
                        # Create new state with job leaving the system
                        new_state = state.copy()
                        new_state[ist * K + k] -= 1

                        if np.all(new_state >= 0):
                            ns = _find_state_index_fast(state_map, new_state)
                            if ns >= 0:
                                Q[s, ns] += rate * exit_prob

    # Make valid generator (set diagonal)
    Q = ctmc_makeinfgen(Q)

    return Q


def _compute_metrics_from_distribution(
    sn: NetworkStruct,
    pi: np.ndarray,
    state_space: np.ndarray
) -> Dict[str, np.ndarray]:
    """
    Compute performance metrics from steady-state distribution.

    Args:
        sn: Network structure
        pi: Steady-state probability distribution
        state_space: State space matrix

    Returns:
        Dictionary with Q, U, R, T matrices
    """
    M = sn.nstations
    K = sn.nclasses

    # Initialize metrics
    QN = np.zeros((M, K))
    UN = np.zeros((M, K))
    RN = np.zeros((M, K))
    TN = np.zeros((M, K))

    # Get service rates
    if hasattr(sn, 'rates') and sn.rates is not None:
        rates = np.asarray(sn.rates)
    else:
        rates = np.ones((M, K))

    # Get number of servers
    if hasattr(sn, 'nservers') and sn.nservers is not None:
        nservers = np.asarray(sn.nservers).flatten()
    else:
        nservers = np.ones(M)

    # Get load-dependent scaling
    lldscaling = None
    if hasattr(sn, 'lldscaling') and sn.lldscaling is not None:
        lldscaling = np.asarray(sn.lldscaling)

    # Track which stations are infinite servers
    inf_server_stations = set()
    if hasattr(sn, 'sched') and sn.sched is not None:
        for ist, sched_val in sn.sched.items():
            if sched_val == SchedStrategy.INF:
                inf_server_stations.add(ist)

    def get_load_scaling(ist: int, total_jobs: int) -> float:
        """Get service rate scaling factor for station ist with total_jobs present."""
        if total_jobs <= 0:
            return 0.0  # No jobs, no service
        # For infinite servers, each job gets its own server
        if ist in inf_server_stations:
            return float(total_jobs)
        if lldscaling is not None and ist < lldscaling.shape[0]:
            idx = total_jobs - 1
            if lldscaling.ndim > 1 and idx < lldscaling.shape[1]:
                return lldscaling[ist, idx]
            elif lldscaling.ndim > 1:
                return lldscaling[ist, -1]
        c = nservers[ist] if ist < len(nservers) else 1
        return min(total_jobs, c)

    def get_max_scaling(ist: int) -> float:
        """Get maximum scaling factor for station (for utilization normalization)."""
        if lldscaling is not None and ist < lldscaling.shape[0]:
            if lldscaling.ndim > 1:
                return np.max(lldscaling[ist, :])
            return lldscaling[ist]
        c = nservers[ist] if ist < len(nservers) else 1
        return c

    # Compute expected queue lengths E[n_ik]
    for s, state in enumerate(state_space):
        state_matrix = state.reshape(M, K)
        for ist in range(M):
            for k in range(K):
                QN[ist, k] += pi[s] * state_matrix[ist, k]

    # Compute throughputs from expected service completions
    # T_i,k = E[service rate for class k at station i]
    for s, state in enumerate(state_space):
        state_matrix = state.reshape(M, K)
        for ist in range(M):
            total_at_station = np.sum(state_matrix[ist, :])
            if total_at_station <= 0:
                continue

            scaling = get_load_scaling(ist, int(total_at_station))

            if hasattr(sn, 'sched') and sn.sched is not None:
                sched = sn.sched.get(ist, SchedStrategy.FCFS)
            else:
                sched = SchedStrategy.FCFS

            for k in range(K):
                n_k = state_matrix[ist, k]
                if n_k <= 0:
                    continue

                mu = rates[ist, k] if ist < rates.shape[0] and k < rates.shape[1] else 1.0

                if sched == SchedStrategy.INF:
                    # Infinite server: each job served independently
                    TN[ist, k] += pi[s] * n_k * mu
                else:
                    # PS, FCFS, etc.: share of capacity proportional to number of jobs
                    TN[ist, k] += pi[s] * scaling * mu * (n_k / total_at_station)

    # Compute expected capacity (E[scaling]) at each station
    # E[scaling] = sum over states of pi(s) * lldscaling(station, n)
    expected_scaling = np.zeros(M)
    for s, state in enumerate(state_space):
        state_matrix = state.reshape(M, K)
        for ist in range(M):
            total_at_station = np.sum(state_matrix[ist, :])
            scaling = get_load_scaling(ist, int(total_at_station))
            expected_scaling[ist] += pi[s] * scaling
    # Ensure no division by zero
    expected_scaling = np.maximum(expected_scaling, 1e-10)

    # Compute utilizations using MATLAB's formulas
    # For INF servers: UN = QN
    # For LLD PS/DPS/etc: UN = E[n_k / n_total] (expected fraction of jobs that are class k)
    for ist in range(M):
        if hasattr(sn, 'sched') and sn.sched is not None:
            sched = sn.sched.get(ist, SchedStrategy.FCFS)
        else:
            sched = SchedStrategy.FCFS

        if sched == SchedStrategy.INF:
            # For infinite servers, utilization = queue length
            for k in range(K):
                UN[ist, k] = QN[ist, k]
        elif sched in [SchedStrategy.PS, SchedStrategy.DPS, SchedStrategy.GPS]:
            # For PS/DPS with load-dependence: UN = E[n_k / n_total]
            # This is already computed below via state iteration
            pass
        else:
            # For FCFS and others: UN = T / (c * mu) where c is number of servers
            c = nservers[ist] if ist < len(nservers) else 1
            for k in range(K):
                mu = rates[ist, k] if ist < rates.shape[0] and k < rates.shape[1] else 1.0
                if mu > 0 and c > 0:
                    UN[ist, k] = TN[ist, k] / (c * mu)

    # For PS/DPS with LLD, compute UN = E[n_k / n_total]
    for s, state in enumerate(state_space):
        state_matrix = state.reshape(M, K)
        for ist in range(M):
            if hasattr(sn, 'sched') and sn.sched is not None:
                sched = sn.sched.get(ist, SchedStrategy.FCFS)
            else:
                sched = SchedStrategy.FCFS

            if sched in [SchedStrategy.PS, SchedStrategy.DPS, SchedStrategy.GPS]:
                n_total = np.sum(state_matrix[ist, :])
                if n_total > 0:
                    for k in range(K):
                        n_k = state_matrix[ist, k]
                        # With equal schedparam weights: n_k / n_total
                        UN[ist, k] += pi[s] * n_k / n_total

    # Compute response times: R = Q / T (Little's law)
    for ist in range(M):
        for k in range(K):
            if TN[ist, k] > 0:
                RN[ist, k] = QN[ist, k] / TN[ist, k]
            else:
                RN[ist, k] = 0

    # Handle Source station metrics for open networks
    # Source doesn't hold jobs - QLen represents cutoff buffer which should be 0 for reporting
    # Throughput should be the arrival rate
    if hasattr(sn, 'sched') and sn.sched is not None:
        for ist, sched_val in sn.sched.items():
            # Check for EXT scheduling (Source station)
            is_source = (sched_val == SchedStrategy.EXT or
                        (isinstance(sched_val, int) and sched_val == 11))
            if is_source and ist < M:
                for k in range(K):
                    # Source doesn't hold jobs in the traditional sense
                    QN[ist, k] = 0.0
                    UN[ist, k] = 0.0
                    RN[ist, k] = 0.0
                    # Throughput is the arrival rate
                    if ist < rates.shape[0] and k < rates.shape[1]:
                        TN[ist, k] = rates[ist, k]

    return {
        'Q': QN,
        'U': UN,
        'R': RN,
        'T': TN
    }


def solver_ctmc_basic(
    sn: NetworkStruct,
    options: Optional[SolverCTMCOptions] = None
) -> SolverCTMCReturn:
    """
    Basic CTMC solver using state-space enumeration.

    Enumerates all valid states, builds the infinitesimal generator,
    and solves for steady-state distribution.

    Args:
        sn: Network structure
        options: Solver options

    Returns:
        SolverCTMCReturn with performance metrics
    """
    start_time = time.time()

    if options is None:
        options = SolverCTMCOptions()

    M = sn.nstations
    K = sn.nclasses

    # Enumerate state space
    state_space, state_space_aggr = _enumerate_state_space(sn, options.cutoff)

    # Mandatory truncation warning for open/mixed models
    if sn_has_open_classes(sn):
        print(f"CTMC solver using state space cutoff = {options.cutoff} for open/mixed model.")
        warnings.warn(
            "State space truncation may cause inaccurate results. "
            "Consider varying cutoff to assess sensitivity.",
            UserWarning
        )

    if options.verbose:
        print(f"CTMC state space size: {len(state_space)}")

    # Build generator matrix
    Q = _build_generator(sn, state_space, options)

    # Solve for steady-state distribution
    pi = ctmc_solve(Q)

    # Compute performance metrics
    metrics = _compute_metrics_from_distribution(sn, pi, state_space)

    QN = metrics['Q']
    UN = metrics['U']
    RN = metrics['R']
    TN = metrics['T']

    # Compute cycle times and system throughput
    CN = np.sum(RN, axis=0).reshape(1, -1)
    XN = np.zeros((1, K))
    for k in range(K):
        ref_stat = int(sn.refstat[k]) if hasattr(sn, 'refstat') and k < len(sn.refstat) else 0
        if ref_stat < M:
            XN[0, k] = TN[ref_stat, k]

    # Clean up NaN values
    QN = np.nan_to_num(QN, nan=0.0)
    UN = np.nan_to_num(UN, nan=0.0)
    RN = np.nan_to_num(RN, nan=0.0)
    TN = np.nan_to_num(TN, nan=0.0)
    CN = np.nan_to_num(CN, nan=0.0)
    XN = np.nan_to_num(XN, nan=0.0)

    result = SolverCTMCReturn()
    result.Q = QN
    result.U = UN
    result.R = RN
    result.T = TN
    result.C = CN
    result.X = XN
    result.pi = pi
    result.infgen = Q
    result.space = state_space
    result.runtime = time.time() - start_time
    result.method = "basic"

    return result


def solver_ctmc(
    sn: NetworkStruct,
    options: Optional[SolverCTMCOptions] = None
) -> SolverCTMCReturn:
    """
    Main CTMC solver handler.

    Routes to appropriate method based on options and network characteristics.

    Args:
        sn: Network structure
        options: Solver options

    Returns:
        SolverCTMCReturn with performance metrics
    """
    if options is None:
        options = SolverCTMCOptions()

    method = options.method.lower()

    if method in ['default', 'basic']:
        return solver_ctmc_basic(sn, options)
    else:
        # Unknown method - use basic
        if options.verbose:
            print(f"Warning: Unknown CTMC method '{method}'. Using basic.")
        return solver_ctmc_basic(sn, options)


__all__ = [
    'solver_ctmc',
    'solver_ctmc_basic',
    'SolverCTMCReturn',
    'SolverCTMCOptions',
]
