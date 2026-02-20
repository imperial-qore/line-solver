"""
SSA Solver handler.

Native Python implementation of SSA (Stochastic Simulation Algorithm) solver
handler that analyzes queueing networks through discrete-event simulation.

The SSA solver uses Gillespie's algorithm to simulate sample paths and
estimate steady-state performance metrics.

Port from:

"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, Dict, List, Tuple, Any
import time
from collections import deque

from ...sn import (
    NetworkStruct,
    SchedStrategy,
    NodeType,
    RoutingStrategy,
    sn_is_open_model,
    sn_is_closed_model,
    sn_has_open_classes,
    sn_get_arvr_from_tput,
)
from ....lang.base import ReplacementStrategy


# =============================================================================
# Cache State Management
# =============================================================================

@dataclass
class CacheState:
    """
    State of a single cache node during simulation.

    Uses a deque to represent the cache where index 0 is the most
    recently used (front) and index -1 is least recently used (back).
    Items are stored as 1-indexed (to match MATLAB).

    Supports multiple replacement policies:
    - LRU: Move to front on hit, evict from back on miss
    - FIFO: Don't move on hit, evict from back (oldest) on miss
    - RR: Don't move on hit, evict random item on miss
    - SFIFO: Segmented FIFO (treated as FIFO for simplicity)
    """
    node_idx: int  # Node index in the network
    capacity: int  # Total cache capacity
    items: deque  # Items currently in cache (LRU order, front = MRU)
    num_items: int  # Total number of possible items (catalog size)
    pread: np.ndarray  # Item access probabilities (Zipf distribution)
    hitclass: np.ndarray  # hitclass[k] = output class for hits from class k
    missclass: np.ndarray  # missclass[k] = output class for misses from class k
    replacement: int = ReplacementStrategy.LRU  # Replacement policy

    # Statistics
    hits: int = 0
    misses: int = 0

    def lookup_and_update(self, item: int, rng: np.random.Generator = None) -> bool:
        """
        Look up an item in the cache and update state using the configured
        replacement policy.

        Args:
            item: Item to look up (1-indexed)
            rng: Random number generator for RR policy (uses global if None)

        Returns:
            True if cache hit, False if cache miss
        """
        if item in self.items:
            # Cache HIT
            self.hits += 1
            if self.replacement == ReplacementStrategy.LRU:
                # LRU: move item to front (MRU position)
                self.items.remove(item)
                self.items.appendleft(item)
            # FIFO, RR, SFIFO: don't move item on hit
            return True
        else:
            # Cache MISS - add item, evict if necessary
            self.misses += 1
            if len(self.items) >= self.capacity:
                # Need to evict an item
                if self.replacement == ReplacementStrategy.RR:
                    # RR: evict random item
                    items_list = list(self.items)
                    if rng is not None:
                        evict_idx = rng.integers(0, len(items_list))
                    else:
                        evict_idx = np.random.randint(0, len(items_list))
                    evicted = items_list[evict_idx]
                    self.items.remove(evicted)
                else:
                    # LRU, FIFO, SFIFO: evict from back (oldest/LRU)
                    self.items.pop()
            # Add new item to front
            self.items.appendleft(item)
            return False

    def sample_item(self, rng: np.random.Generator = None) -> int:
        """
        Sample which item is being requested according to the access distribution.

        Args:
            rng: Random number generator (uses global if None)

        Returns:
            Item index (1-indexed)
        """
        if rng is not None:
            r = rng.random()
        else:
            r = np.random.random()

        # Cumulative distribution sampling
        cumsum = np.cumsum(self.pread)
        for i, c in enumerate(cumsum):
            if r <= c:
                return i + 1  # 1-indexed items
        return len(self.pread)  # Return last item if numerical issues


# =============================================================================
# Detailed Queue State for Priority Scheduling (HOL, FCFS, etc.)
# =============================================================================

@dataclass
class DetailedQueueState:
    """
    Detailed state of a queueing station for accurate simulation of
    buffer-order-dependent scheduling strategies like HOL, FCFS, LCFS.

    Tracks which jobs are in service vs in buffer separately, enabling
    accurate modeling of non-preemptive priority scheduling.
    """
    station_idx: int  # Station index
    nservers: int  # Number of servers
    nclasses: int  # Number of job classes
    sched_strategy: str  # Scheduling strategy name ('HOL', 'FCFS', etc.)
    classprio: np.ndarray  # Priority of each class (lower = higher priority)

    # State tracking
    in_service: np.ndarray = field(default=None)  # in_service[class] = jobs from class in service
    buffer: List = field(default_factory=list)  # Ordered list of class indices in buffer (FIFO order)

    def __post_init__(self):
        if self.in_service is None:
            self.in_service = np.zeros(self.nclasses, dtype=int)

    def total_jobs(self) -> int:
        """Total jobs at station (in service + in buffer)."""
        return int(np.sum(self.in_service)) + len(self.buffer)

    def jobs_in_service(self) -> int:
        """Total jobs currently being served."""
        return int(np.sum(self.in_service))

    def jobs_in_buffer(self) -> int:
        """Total jobs waiting in buffer."""
        return len(self.buffer)

    def jobs_per_class(self) -> np.ndarray:
        """Get total jobs per class (service + buffer)."""
        result = self.in_service.copy()
        for class_idx in self.buffer:
            result[class_idx] += 1
        return result

    def process_arrival(self, class_idx: int) -> bool:
        """
        Process a job arrival.

        Args:
            class_idx: Class index of arriving job

        Returns:
            True if job entered service, False if queued in buffer
        """
        if self.jobs_in_service() < self.nservers:
            # Server available - job enters service immediately
            self.in_service[class_idx] += 1
            return True
        else:
            # All servers busy - job goes to buffer
            self.buffer.append(class_idx)
            return False

    def process_departure(self, class_idx: int) -> Optional[int]:
        """
        Process a job departure from service.

        Args:
            class_idx: Class index of departing job

        Returns:
            Class index of job that enters service from buffer, or None if buffer empty
        """
        if self.in_service[class_idx] <= 0:
            return None

        # Remove job from service
        self.in_service[class_idx] -= 1

        # Select next job from buffer based on scheduling strategy
        if len(self.buffer) == 0:
            return None

        next_job_class = self._select_from_buffer()
        if next_job_class is not None:
            self.in_service[next_job_class] += 1

        return next_job_class

    def _select_from_buffer(self) -> Optional[int]:
        """
        Select next job from buffer based on scheduling strategy.

        Returns:
            Class index of selected job, or None if buffer empty
        """
        if len(self.buffer) == 0:
            return None

        if self.sched_strategy == 'HOL':
            # HOL: Select highest priority (lowest priority value), FIFO within priority
            # Find the highest priority among jobs in buffer
            min_prio = float('inf')
            for class_idx in self.buffer:
                prio = self.classprio[class_idx] if class_idx < len(self.classprio) else 0
                if prio < min_prio:
                    min_prio = prio

            # Find the first (FIFO) job with that priority
            for i, class_idx in enumerate(self.buffer):
                prio = self.classprio[class_idx] if class_idx < len(self.classprio) else 0
                if prio == min_prio:
                    self.buffer.pop(i)
                    return class_idx

        elif self.sched_strategy == 'LCFS':
            # LCFS: Take the last job (most recent arrival)
            return self.buffer.pop()

        elif self.sched_strategy == 'LCFSPRIO':
            # LCFS with priority: highest priority first, LCFS within priority
            min_prio = float('inf')
            for class_idx in self.buffer:
                prio = self.classprio[class_idx] if class_idx < len(self.classprio) else 0
                if prio < min_prio:
                    min_prio = prio

            # Find the last (LCFS) job with that priority
            for i in range(len(self.buffer) - 1, -1, -1):
                class_idx = self.buffer[i]
                prio = self.classprio[class_idx] if class_idx < len(self.classprio) else 0
                if prio == min_prio:
                    self.buffer.pop(i)
                    return class_idx

        else:
            # Default FCFS: Take the first job
            return self.buffer.pop(0)

        return None

    def get_departure_rate(self, class_idx: int, service_rate: float) -> float:
        """
        Get departure rate for a specific class.

        For detailed state, rate = service_rate * jobs_in_service[class]

        Args:
            class_idx: Class index
            service_rate: Service rate (mu) for this class

        Returns:
            Total departure rate for this class
        """
        return service_rate * self.in_service[class_idx]

    def copy(self) -> 'DetailedQueueState':
        """Create a deep copy of this state."""
        new_state = DetailedQueueState(
            station_idx=self.station_idx,
            nservers=self.nservers,
            nclasses=self.nclasses,
            sched_strategy=self.sched_strategy,
            classprio=self.classprio.copy() if self.classprio is not None else None,
            in_service=self.in_service.copy(),
            buffer=list(self.buffer)
        )
        return new_state


def _needs_detailed_state(sched_strategy) -> bool:
    """
    Check if a scheduling strategy requires detailed state tracking.

    Args:
        sched_strategy: Scheduling strategy (enum or string)

    Returns:
        True if detailed state is needed for accurate simulation
    """
    sched_name = sched_strategy.name if hasattr(sched_strategy, 'name') else str(sched_strategy)
    detailed_scheds = {'HOL', 'FCFS', 'LCFS', 'LCFSPRIO', 'SIRO', 'SEPT', 'LEPT'}
    return sched_name in detailed_scheds


def _init_detailed_queue_states(sn: NetworkStruct) -> Dict[int, DetailedQueueState]:
    """
    Initialize detailed queue states for stations that need them.

    Args:
        sn: Network structure

    Returns:
        Dictionary mapping station index to DetailedQueueState
    """
    detailed_states = {}

    M = sn.nstations
    K = sn.nclasses
    classprio = getattr(sn, 'classprio', None)
    if classprio is not None:
        classprio = classprio.flatten() if classprio.ndim > 1 else classprio
    else:
        classprio = np.zeros(K)

    nservers = getattr(sn, 'nservers', None)
    if nservers is None:
        nservers = np.ones(M)
    else:
        nservers = np.atleast_1d(nservers).flatten()

    for ist in range(M):
        # Get scheduling strategy (sn.sched is a dict mapping station index to strategy)
        if hasattr(sn, 'sched') and sn.sched is not None:
            sched = sn.sched.get(ist, None)
        else:
            sched = None
        if sched is None:
            continue

        sched_name = sched.name if hasattr(sched, 'name') else str(sched)

        if _needs_detailed_state(sched):
            c = int(nservers[ist]) if ist < len(nservers) else 1
            if np.isinf(c):
                c = 10000  # Large number for "infinite" servers

            detailed_states[ist] = DetailedQueueState(
                station_idx=ist,
                nservers=c,
                nclasses=K,
                sched_strategy=sched_name,
                classprio=classprio.copy()
            )

    return detailed_states


def _sync_detailed_states_from_aggregate(
    detailed_states: Dict[int, DetailedQueueState],
    state: np.ndarray,
    sn: NetworkStruct
) -> None:
    """
    Synchronize detailed queue states from the aggregate state array.

    This is critical for correct initialization: the detailed states must
    reflect the initial job distribution from the aggregate state.

    For each station with detailed state tracking, we populate jobs into
    service slots first (up to nservers), then remaining jobs go to buffer.
    For HOL scheduling, we process classes in priority order.

    Args:
        detailed_states: Dictionary of DetailedQueueState objects to populate
        state: Aggregate state array [M x K] with job counts
        sn: Network structure
    """
    if not detailed_states:
        return

    K = sn.nclasses
    classprio = getattr(sn, 'classprio', None)
    if classprio is not None:
        classprio = classprio.flatten() if classprio.ndim > 1 else classprio
    else:
        classprio = np.zeros(K)

    for ist, ds in detailed_states.items():
        # Reset state
        ds.in_service = np.zeros(ds.nclasses, dtype=int)
        ds.buffer = []

        # Get jobs at this station from aggregate state
        jobs_per_class = state[ist, :].astype(int)
        total_jobs = int(np.sum(jobs_per_class))

        if total_jobs == 0:
            continue

        # For HOL scheduling, allocate servers by priority order
        # Lower priority value = higher priority
        if ds.sched_strategy == 'HOL':
            # Get classes sorted by priority (highest priority first)
            class_order = sorted(range(K), key=lambda k: classprio[k] if k < len(classprio) else 0)

            servers_used = 0
            for k in class_order:
                n_k = jobs_per_class[k]
                if n_k <= 0:
                    continue

                # How many of class k can enter service?
                can_serve = min(n_k, ds.nservers - servers_used)
                ds.in_service[k] = can_serve
                servers_used += can_serve

                # Remaining jobs go to buffer
                for _ in range(n_k - can_serve):
                    ds.buffer.append(k)
        else:
            # For FCFS/SIRO/etc: allocate servers proportionally or FIFO
            # Simple approach: fill servers in class order, rest to buffer
            servers_used = 0
            for k in range(K):
                n_k = jobs_per_class[k]
                if n_k <= 0:
                    continue

                can_serve = min(n_k, ds.nservers - servers_used)
                ds.in_service[k] = can_serve
                servers_used += can_serve

                for _ in range(n_k - can_serve):
                    ds.buffer.append(k)


def _init_cache_states(
    sn: NetworkStruct,
    model = None
) -> Dict[int, CacheState]:
    """
    Initialize cache states for all cache nodes in the network.

    Args:
        sn: Network structure
        model: Optional model object with cache node details

    Returns:
        Dictionary mapping node index to CacheState
    """
    cache_states = {}

    if not hasattr(sn, 'nodetype') or sn.nodetype is None:
        return cache_states

    K = sn.nclasses

    for ind in range(sn.nnodes):
        if ind >= len(sn.nodetype):
            continue

        if sn.nodetype[ind] != NodeType.CACHE:
            continue

        # Get cache parameters from nodeparam
        if not hasattr(sn, 'nodeparam') or sn.nodeparam is None:
            continue

        if ind not in sn.nodeparam:
            continue

        cache_param = sn.nodeparam[ind]

        # Extract parameters
        nitems = getattr(cache_param, 'nitems', 100)
        itemcap = getattr(cache_param, 'itemcap', np.array([10]))
        if hasattr(itemcap, '__iter__'):
            capacity = int(np.sum(itemcap))
        else:
            capacity = int(itemcap)

        # Get replacement strategy (default to LRU)
        replacement = ReplacementStrategy.LRU
        if hasattr(cache_param, 'replacestrat'):
            rs = cache_param.replacestrat
            if rs is not None:
                if isinstance(rs, int):
                    replacement = rs
                elif hasattr(rs, 'value'):
                    replacement = rs.value
                else:
                    replacement = int(rs)

        # Get hit/miss class mappings
        hitclass = np.zeros(K, dtype=int)
        missclass = np.zeros(K, dtype=int)

        if hasattr(cache_param, 'hitclass'):
            hc = np.atleast_1d(cache_param.hitclass).flatten()
            for k in range(min(len(hc), K)):
                hitclass[k] = int(hc[k])

        if hasattr(cache_param, 'missclass'):
            mc = np.atleast_1d(cache_param.missclass).flatten()
            for k in range(min(len(mc), K)):
                missclass[k] = int(mc[k])

        # Get access probabilities (pread)
        pread = np.ones(nitems) / nitems  # Default: uniform

        if hasattr(cache_param, 'pread') and cache_param.pread is not None:
            # pread can be a dict or a list
            pread_data = cache_param.pread
            if isinstance(pread_data, dict):
                # Dict mapping class index to probability list
                for class_k, probs in pread_data.items():
                    if probs is not None:
                        pread = np.asarray(probs).flatten()
                        pread_sum = np.sum(pread)
                        if pread_sum > 0:
                            pread = pread / pread_sum
                        break
            elif isinstance(pread_data, list):
                # List of probability arrays, one per class
                for probs in pread_data:
                    if probs is not None:
                        probs_arr = np.asarray(probs).flatten()
                        if len(probs_arr) > 0 and np.sum(probs_arr) > 0:
                            pread = probs_arr
                            pread_sum = np.sum(pread)
                            if pread_sum > 0:
                                pread = pread / pread_sum
                            break

        # Try to get from model's cache node if available
        if model is not None and hasattr(model, '_nodes'):
            try:
                cache_node = model._nodes[ind]
                if hasattr(cache_node, '_read_process') and cache_node._read_process:
                    for jobclass, dist in cache_node._read_process.items():
                        if dist is not None:
                            # Zipf distribution
                            if hasattr(dist, '_s') and hasattr(dist, '_H'):
                                k = np.arange(1, nitems + 1)
                                pread = (1.0 / k ** dist._s) / dist._H
                                break
                            # DiscreteSampler with probs
                            elif hasattr(dist, '_probs'):
                                probs = np.asarray(dist._probs).flatten()
                                if len(probs) >= nitems:
                                    pread = probs[:nitems]
                                else:
                                    pread = np.zeros(nitems)
                                    pread[:len(probs)] = probs
                                pread_sum = np.sum(pread)
                                if pread_sum > 0:
                                    pread = pread / pread_sum
                                break
            except Exception:
                pass

        # Create cache state (initially empty)
        cache_states[ind] = CacheState(
            node_idx=ind,
            capacity=capacity,
            items=deque(maxlen=capacity),
            num_items=nitems,
            pread=pread,
            hitclass=hitclass,
            missclass=missclass,
            replacement=replacement,
            hits=0,
            misses=0
        )

    return cache_states


def _get_cache_station_mapping(sn: NetworkStruct) -> Dict[int, int]:
    """
    Get mapping from cache node indices to station indices.

    Returns:
        Dict mapping cache node index to station index
    """
    cache_to_station = {}

    if not hasattr(sn, 'nodetype') or sn.nodetype is None:
        return cache_to_station

    if not hasattr(sn, 'nodeToStation') or sn.nodeToStation is None:
        return cache_to_station

    for ind in range(sn.nnodes):
        if ind >= len(sn.nodetype):
            continue

        if sn.nodetype[ind] == NodeType.CACHE:
            ist = sn.nodeToStation[ind]
            if ist >= 0:
                cache_to_station[ind] = int(ist)

    return cache_to_station


@dataclass
class SolverSSAOptions:
    """Options for SSA solver."""
    method: str = 'default'
    tol: float = 1e-6
    verbose: bool = False
    samples: int = 10000  # Number of simulation samples/events
    timespan: Tuple[float, float] = (0.0, float('inf'))  # Simulation time window
    seed: int = 0  # Random seed for reproducibility
    cutoff: float = np.inf  # Match MATLAB default: infinity means no truncation
    confidence_level: float = 0.95  # Confidence level for CI


@dataclass
class SolverSSAReturn:
    """
    Result of SSA solver handler.

    Attributes:
        Q: Mean queue lengths (M x K)
        U: Utilizations (M x K)
        R: Response times (M x K)
        T: Throughputs (M x K)
        C: Cycle times (1 x K)
        X: System throughputs (1 x K)
        Q_ci: Queue length confidence intervals
        U_ci: Utilization confidence intervals
        R_ci: Response time confidence intervals
        T_ci: Throughput confidence intervals
        total_time: Total simulated time
        runtime: Runtime in seconds
        method: Method used
        samples: Number of samples collected
        A: Arrival rates (M x K)
    """
    Q: Optional[np.ndarray] = None
    U: Optional[np.ndarray] = None
    R: Optional[np.ndarray] = None
    T: Optional[np.ndarray] = None
    A: Optional[np.ndarray] = None
    C: Optional[np.ndarray] = None
    X: Optional[np.ndarray] = None
    Q_ci: Optional[np.ndarray] = None
    U_ci: Optional[np.ndarray] = None
    R_ci: Optional[np.ndarray] = None
    T_ci: Optional[np.ndarray] = None
    total_time: float = 0.0
    runtime: float = 0.0
    method: str = "default"
    samples: int = 0


def _init_state(
    sn: NetworkStruct,
    cutoff: int = 10
) -> np.ndarray:
    """
    Initialize the simulation state.

    For closed networks, places all jobs at reference stations.
    For open networks, starts with empty queues.

    Args:
        sn: Network structure
        cutoff: Population cutoff for open classes

    Returns:
        Initial state vector (M x K)
    """
    M = sn.nstations
    K = sn.nclasses
    state = np.zeros((M, K))

    if sn.njobs is not None:
        N = sn.njobs.flatten()
    else:
        N = np.zeros(K)

    for k in range(K):
        if np.isfinite(N[k]) and N[k] > 0:
            # Closed class: place jobs at reference station
            ref_stat = int(sn.refstat[k]) if hasattr(sn, 'refstat') and k < len(sn.refstat) else 0
            if ref_stat < M:
                state[ref_stat, k] = N[k]

    return state


def _init_rrobin_state(
    sn: NetworkStruct
) -> Dict[Tuple[int, int], Tuple[List[Tuple[int, int]], int]]:
    """
    Initialize round-robin routing state.

    For each (station, class) pair that uses RROBIN routing, this function
    identifies the list of valid destinations and initializes the counter.

    Args:
        sn: Network structure

    Returns:
        Dict mapping (station, class) to (destinations_list, current_index)
        where destinations_list is a list of (dest_station, dest_class) tuples
    """
    rrobin_state = {}
    M = sn.nstations
    K = sn.nclasses

    # Check if routing strategies are defined
    if not hasattr(sn, 'routing') or sn.routing is None or sn.routing.size == 0:
        return rrobin_state

    # Get the routing probability matrix
    if hasattr(sn, 'rt') and sn.rt is not None:
        P = np.asarray(sn.rt)
    else:
        return rrobin_state

    # Get node-to-station mapping
    nodeToStation = None
    if hasattr(sn, 'nodeToStation') and sn.nodeToStation is not None:
        nodeToStation = np.asarray(sn.nodeToStation).flatten()

    routing = np.asarray(sn.routing)

    # Iterate over all nodes and classes to find RROBIN routing
    N = sn.nnodes if hasattr(sn, 'nnodes') else M
    for ind in range(N):
        for k in range(K):
            # Check routing strategy for this node-class pair
            if ind < routing.shape[0] and k < routing.shape[1]:
                strategy = routing[ind, k]
                # Check for RROBIN (value 3) - handle both int and enum
                is_rrobin = (strategy == RoutingStrategy.RROBIN or
                            strategy == 3 or
                            (hasattr(strategy, 'value') and strategy.value == 3))

                if is_rrobin:
                    # Find the station index for this node
                    if nodeToStation is not None and ind < len(nodeToStation):
                        ist = int(nodeToStation[ind])
                        if ist < 0:
                            continue  # Not a station
                    else:
                        ist = ind  # Assume node index = station index

                    # Find all valid destinations from the routing matrix
                    src_idx = ist * K + k
                    destinations = []

                    for jst in range(M):
                        for r in range(K):
                            dst_idx = jst * K + r
                            if src_idx < P.shape[0] and dst_idx < P.shape[1]:
                                p = P[src_idx, dst_idx]
                                if p > 0:
                                    destinations.append((jst, r))

                    if destinations:
                        rrobin_state[(ist, k)] = (destinations, 0)

    return rrobin_state


def _get_rrobin_destination(
    rrobin_state: Dict[Tuple[int, int], Tuple[List[Tuple[int, int]], int]],
    ist: int,
    k: int
) -> Optional[Tuple[int, int]]:
    """
    Get the next destination for round-robin routing.

    Args:
        rrobin_state: Current round-robin state
        ist: Source station index
        k: Source class index

    Returns:
        (dest_station, dest_class) or None if not using RROBIN
    """
    key = (ist, k)
    if key not in rrobin_state:
        return None

    destinations, idx = rrobin_state[key]
    if not destinations:
        return None

    return destinations[idx]


def _advance_rrobin(
    rrobin_state: Dict[Tuple[int, int], Tuple[List[Tuple[int, int]], int]],
    ist: int,
    k: int
) -> None:
    """
    Advance the round-robin counter after a routing decision.

    Args:
        rrobin_state: Current round-robin state (modified in-place)
        ist: Source station index
        k: Source class index
    """
    key = (ist, k)
    if key not in rrobin_state:
        return

    destinations, idx = rrobin_state[key]
    if destinations:
        new_idx = (idx + 1) % len(destinations)
        rrobin_state[key] = (destinations, new_idx)


def _resolve_routing_through_nonstation_nodes(
    sn: NetworkStruct,
    start_node: int,
    start_class: int,
    P_nodes: np.ndarray,
    max_hops: int = 10
) -> List[Tuple[int, int, float]]:
    """
    Follow routing through non-station nodes (ClassSwitch, Router) until reaching stations.

    Args:
        sn: Network structure
        start_node: Starting node index (non-station)
        start_class: Starting class index
        P_nodes: Node-indexed routing probability matrix (rtnodes)
        max_hops: Maximum number of hops to prevent infinite loops

    Returns:
        List of (station_index, class_index, probability) tuples for final destinations
    """
    I = sn.nnodes
    K = sn.nclasses

    nodeToStation = None
    if hasattr(sn, 'nodeToStation') and sn.nodeToStation is not None:
        nodeToStation = np.asarray(sn.nodeToStation).flatten().astype(int)
    else:
        nodeToStation = np.arange(I, dtype=int)

    isstation = None
    if hasattr(sn, 'isstation') and sn.isstation is not None:
        isstation = np.asarray(sn.isstation).flatten()

    # Use BFS to resolve all paths through non-station nodes
    # Each entry: (node, class, accumulated_probability)
    queue = [(start_node, start_class, 1.0)]
    results = []
    visited_hops = 0

    while queue and visited_hops < max_hops * len(queue):
        visited_hops += 1
        current_node, current_class, prob = queue.pop(0)

        if prob < 1e-10:
            continue

        src_idx = current_node * K + current_class

        # Find all destinations from this node/class
        found_any = False
        for jnd in range(I):
            for r in range(K):
                dst_idx = jnd * K + r
                if src_idx < P_nodes.shape[0] and dst_idx < P_nodes.shape[1]:
                    p = P_nodes[src_idx, dst_idx]
                    if p > 0:
                        found_any = True
                        new_prob = prob * p

                        # Check if destination is a station
                        if isstation is not None and jnd < len(isstation) and isstation[jnd]:
                            jst = nodeToStation[jnd]
                            if jst >= 0:
                                results.append((jst, r, new_prob))
                        elif hasattr(sn, 'nodetype') and sn.nodetype is not None and jnd < len(sn.nodetype):
                            nt = sn.nodetype[jnd]
                            nt_val = int(nt.value) if hasattr(nt, 'value') else int(nt)

                            if nt_val == 1:  # NodeType.SINK = 1
                                # Job exits to Sink
                                results.append((-2, r, new_prob))
                            elif nt_val == 6:  # NodeType.CACHE = 6
                                # Reached a Cache node - encode specially
                                results.append((-(100 + jnd), r, new_prob))
                            elif nt_val in (7, 8):  # NodeType.ROUTER = 7, NodeType.CLASSSWITCH = 8
                                # Continue through this non-station node
                                queue.append((jnd, r, new_prob))
                            else:
                                # Other non-station types - treat as sink
                                results.append((-2, r, new_prob))
                        else:
                            # Unknown node type - treat as sink
                            results.append((-2, r, new_prob))

        if not found_any:
            # No outgoing routing - job exits system
            results.append((-2, current_class, prob))

    return results


def _get_enabled_transitions(
    sn: NetworkStruct,
    state: np.ndarray,
    rrobin_state: Optional[Dict] = None,
    detailed_states: Optional[Dict[int, DetailedQueueState]] = None,
    cutoff: int = 10
) -> List[Tuple[int, int, int, int, float]]:
    """
    Find all enabled transitions in current state.

    Returns list of (source_station, source_class, dest_station, dest_class, rate)
    for each enabled transition.

    Args:
        sn: Network structure
        state: Current state (M x K)
        rrobin_state: Round-robin routing state (optional)
        detailed_states: Detailed queue states for HOL/FCFS stations (optional)
        cutoff: Population cutoff for open classes (default: 10)

    Returns:
        List of enabled transitions with their rates
    """
    M = sn.nstations
    K = sn.nclasses
    transitions = []

    # Compute class capacity limits for open classes
    # For open classes (njobs = inf), capacity is limited by cutoff
    # For closed classes, capacity is the total chain population
    classcap = np.full((M, K), np.inf)
    if sn.njobs is not None:
        N = sn.njobs.flatten()
        for k in range(K):
            if np.isinf(N[k]):
                # Open class: apply cutoff
                classcap[:, k] = cutoff
            else:
                # Closed class: capacity is total chain population
                if hasattr(sn, 'chains') and sn.chains is not None:
                    chains = np.asarray(sn.chains)
                    for c in range(chains.shape[0]):
                        if k < chains.shape[1] and chains[c, k] > 0:
                            # Find all classes in this chain and sum their populations
                            chain_pop = 0
                            for kk in range(K):
                                if kk < chains.shape[1] and chains[c, kk] > 0:
                                    if np.isfinite(N[kk]):
                                        chain_pop += N[kk]
                            classcap[:, k] = chain_pop
                            break
                else:
                    classcap[:, k] = N[k]

    # Get service rates
    if hasattr(sn, 'rates') and sn.rates is not None:
        rates = np.asarray(sn.rates)
    else:
        rates = np.ones((M, K))

    # Get routing probabilities - use rtnodes (node-indexed) to properly handle
    # routing to non-station nodes like Sink
    I = sn.nnodes if hasattr(sn, 'nnodes') else M  # Number of nodes
    if hasattr(sn, 'rtnodes') and sn.rtnodes is not None:
        P = np.asarray(sn.rtnodes)
        use_node_routing = True
    elif hasattr(sn, 'rt') and sn.rt is not None:
        P = np.asarray(sn.rt)
        use_node_routing = False
    else:
        # Default: uniform routing to next station
        P = np.zeros((M * K, M * K))
        for i in range(M):
            for k in range(K):
                next_i = (i + 1) % M
                P[i * K + k, next_i * K + k] = 1.0
        use_node_routing = False

    # Station to node mapping (for rtnodes lookup)
    stationToNode = None
    nodeToStation = None
    isstation = None
    if use_node_routing:
        if hasattr(sn, 'stationToNode') and sn.stationToNode is not None:
            stationToNode = np.asarray(sn.stationToNode).flatten().astype(int)
        else:
            stationToNode = np.arange(M, dtype=int)
        if hasattr(sn, 'nodeToStation') and sn.nodeToStation is not None:
            nodeToStation = np.asarray(sn.nodeToStation).flatten().astype(int)
        else:
            nodeToStation = np.arange(I, dtype=int)
        if hasattr(sn, 'isstation') and sn.isstation is not None:
            isstation = np.asarray(sn.isstation).flatten()

    # Get number of servers
    if hasattr(sn, 'nservers') and sn.nservers is not None:
        nservers = np.asarray(sn.nservers).flatten()
    else:
        nservers = np.ones(M)

    # External arrivals (for open/mixed networks)
    # Use sn_has_open_classes instead of sn_is_open_model to properly handle
    # mixed networks (which have both open and closed classes)
    is_open = sn_has_open_classes(sn)
    if is_open:
        for k in range(K):
            # Check if class k has external arrivals
            if sn.njobs is not None:
                N = sn.njobs.flatten()
                if np.isinf(N[k]):
                    # Find arrival rate
                    for ist in range(M):
                        if hasattr(sn, 'sched') and sn.sched is not None:
                            sched = sn.sched.get(ist, SchedStrategy.FCFS)
                        else:
                            sched = SchedStrategy.FCFS

                        # Compare by name to handle enum type mismatches
                        sched_name = sched.name if hasattr(sched, 'name') else str(sched)
                        is_ext = sched_name == 'EXT' or sched == SchedStrategy.EXT

                        if is_ext:
                            arr_rate = rates[ist, k] if ist < rates.shape[0] and k < rates.shape[1] else 0
                            if arr_rate > 0:
                                # External arrival - check all destination nodes and classes
                                # (class switching may occur, e.g., Cache hit/miss)
                                if use_node_routing:
                                    # Use node-indexed routing (rtnodes)
                                    ind = stationToNode[ist]
                                    src_idx = ind * K + k
                                    for jnd in range(I):
                                        for dst_k in range(K):
                                            dst_idx = jnd * K + dst_k
                                            if src_idx < P.shape[0] and dst_idx < P.shape[1]:
                                                p = P[src_idx, dst_idx]
                                                if p > 0:
                                                    # Convert destination node to station or handle Cache
                                                    if isstation is not None and jnd < len(isstation) and isstation[jnd]:
                                                        jst = nodeToStation[jnd]
                                                        # Check capacity before adding arrival transition
                                                        if state[jst, dst_k] < classcap[jst, dst_k]:
                                                            transitions.append((-1, k, jst, dst_k, arr_rate * p))
                                                    elif hasattr(sn, 'nodetype') and sn.nodetype is not None and jnd < len(sn.nodetype):
                                                        # Check destination node type
                                                        nt = sn.nodetype[jnd]
                                                        nt_val = int(nt.value) if hasattr(nt, 'value') else int(nt)
                                                        if nt_val == 6:  # NodeType.CACHE = 6
                                                            # Encode Cache node as special destination: -(100 + node_idx)
                                                            transitions.append((-1, k, -(100 + jnd), dst_k, arr_rate * p))
                                                        elif nt_val in (7, 8):  # NodeType.ROUTER = 7, NodeType.CLASSSWITCH = 8
                                                            # Follow through non-station routing nodes to find final destination
                                                            resolved = _resolve_routing_through_nonstation_nodes(sn, jnd, dst_k, P)
                                                            for final_st, final_class, final_prob in resolved:
                                                                # Check capacity for final destination
                                                                if final_st >= 0 and final_st < M:
                                                                    if state[final_st, final_class] < classcap[final_st, final_class]:
                                                                        transitions.append((-1, k, final_st, final_class, arr_rate * p * final_prob))
                                                                else:
                                                                    # Non-station destination (sink, etc.)
                                                                    transitions.append((-1, k, final_st, final_class, arr_rate * p * final_prob))
                                                    # Other non-station destinations (like immediate routing) are handled elsewhere
                                else:
                                    # Use station-indexed routing (rt)
                                    src_idx = ist * K + k
                                    for jst in range(M):
                                        for dst_k in range(K):
                                            dst_idx = jst * K + dst_k
                                            if src_idx < P.shape[0] and dst_idx < P.shape[1]:
                                                p = P[src_idx, dst_idx]
                                                if p > 0:
                                                    # Check capacity before adding arrival transition
                                                    if state[jst, dst_k] < classcap[jst, dst_k]:
                                                        transitions.append((-1, k, jst, dst_k, arr_rate * p))

    # Service completions
    for ist in range(M):
        total_at_station = np.sum(state[ist, :])
        if total_at_station <= 0:
            continue

        if hasattr(sn, 'sched') and sn.sched is not None:
            sched = sn.sched.get(ist, SchedStrategy.FCFS)
        else:
            sched = SchedStrategy.FCFS

        # Compare by name to handle enum type mismatches
        sched_name = sched.name if hasattr(sched, 'name') else str(sched)

        # Skip source stations
        if sched_name == 'EXT' or sched == SchedStrategy.EXT:
            continue

        # Get class priorities for priority-based scheduling
        classprio = None
        if hasattr(sn, 'classprio') and sn.classprio is not None:
            classprio = sn.classprio.flatten() if sn.classprio.ndim > 1 else sn.classprio

        # For priority scheduling, find the highest priority class(es) present
        # In LINE, lower priority value = higher priority
        is_priority_sched = sched_name in ('HOL', 'PSPRIO', 'DPSPRIO', 'GPSPRIO',
                                            'FCFSPRIO', 'LCFSPRIO', 'FCFSPRPRIO',
                                            'LCFSPRPRIO', 'FCFSPIPRIO', 'LCFSPIPRIO')

        min_prio_present = None
        if is_priority_sched and classprio is not None:
            # Find the minimum (highest) priority among classes with jobs at this station
            prio_values = []
            for kk in range(K):
                if state[ist, kk] > 0:
                    prio_values.append(classprio[kk] if kk < len(classprio) else 0)
            if prio_values:
                min_prio_present = min(prio_values)

        for k in range(K):
            n_ik = state[ist, k]
            if n_ik <= 0:
                continue

            mu = rates[ist, k] if ist < rates.shape[0] and k < rates.shape[1] else 1.0

            if sched_name == 'INF':
                # Infinite server: rate = n * mu
                rate = n_ik * mu
            elif sched_name == 'HOL':
                # HOL (Head of Line) is NON-PREEMPTIVE priority: jobs in service complete
                # regardless of priority. Priority only affects which job enters service NEXT.
                #
                # If detailed state is available, use exact departure rate based on jobs IN SERVICE.
                # Otherwise, use approximate allocation in aggregate state.
                if detailed_states is not None and ist in detailed_states:
                    # Use detailed state for accurate HOL simulation
                    ds = detailed_states[ist]
                    rate = ds.get_departure_rate(k, mu)
                    if rate <= 0:
                        continue
                else:
                    # Fallback: approximate by allocating servers to classes in priority order
                    c = nservers[ist] if ist < len(nservers) else 1
                    if np.isinf(c):
                        c = total_at_station

                    if total_at_station <= c:
                        # All jobs fit in servers - each job departs at its own rate
                        rate = n_ik * mu
                    else:
                        # More jobs than servers - allocate servers by priority
                        class_prio = classprio[k] if classprio is not None and k < len(classprio) else 0

                        # Count jobs in higher priority classes (lower priority value = higher priority)
                        jobs_higher_prio = 0
                        if classprio is not None:
                            for kk in range(K):
                                if state[ist, kk] > 0:
                                    kk_prio = classprio[kk] if kk < len(classprio) else 0
                                    if kk_prio < class_prio:
                                        jobs_higher_prio += state[ist, kk]

                        # Servers available for this class = c - servers used by higher priority
                        servers_available = max(0, c - jobs_higher_prio)
                        # This class uses min(its jobs, available servers)
                        servers_for_class = min(n_ik, servers_available)

                        if servers_for_class <= 0:
                            continue

                        # Departure rate = servers serving this class * mu
                        rate = servers_for_class * mu
            elif sched_name in ('PSPRIO', 'DPSPRIO', 'GPSPRIO'):
                # Priority PS: if ni <= S, all get service; else only highest priority
                c = nservers[ist] if ist < len(nservers) else 1
                if np.isinf(c):
                    c = total_at_station

                if total_at_station <= c:
                    # All jobs can be served - regular PS
                    active_servers = min(total_at_station, c)
                    rate = active_servers * mu * (n_ik / total_at_station)
                else:
                    # More jobs than servers - only highest priority gets service
                    class_prio = classprio[k] if classprio is not None and k < len(classprio) else 0
                    if min_prio_present is not None and class_prio > min_prio_present:
                        # Lower priority - no service
                        continue
                    # Count jobs only in highest priority group
                    n_prio = 0
                    if classprio is not None and min_prio_present is not None:
                        for kk in range(K):
                            if state[ist, kk] > 0:
                                kk_prio = classprio[kk] if kk < len(classprio) else 0
                                if kk_prio == min_prio_present:
                                    n_prio += state[ist, kk]
                    else:
                        n_prio = total_at_station
                    active_servers = min(n_prio, c)
                    rate = active_servers * mu * (n_ik / n_prio) if n_prio > 0 else 0
            elif sched_name == 'FCFS':
                # FCFS: only jobs in service can complete
                # If detailed state is available, use exact departure rate.
                # Otherwise, fall back to PS approximation.
                if detailed_states is not None and ist in detailed_states:
                    ds = detailed_states[ist]
                    rate = ds.get_departure_rate(k, mu)
                    if rate <= 0:
                        continue
                else:
                    # Fallback: PS approximation
                    c = nservers[ist] if ist < len(nservers) else 1
                    if np.isinf(c):
                        c = total_at_station
                    active_servers = min(total_at_station, c)
                    rate = active_servers * mu * (n_ik / total_at_station)
            else:
                # PS, SIRO, etc.: rate = min(n, c) * mu * fraction
                c = nservers[ist] if ist < len(nservers) else 1
                if np.isinf(c):
                    c = total_at_station
                active_servers = min(total_at_station, c)
                rate = active_servers * mu * (n_ik / total_at_station)

            if rate <= 0:
                continue

            # Check if this station-class uses round-robin routing
            rr_dest = None
            if rrobin_state is not None:
                rr_dest = _get_rrobin_destination(rrobin_state, ist, k)

            if rr_dest is not None:
                # Round-robin: deterministic destination
                jst, r = rr_dest
                transitions.append((ist, k, jst, r, rate))
            else:
                # Probabilistic routing: find destinations from routing matrix
                total_routing_prob = 0.0

                if use_node_routing:
                    # Use node-indexed routing (rtnodes) to properly handle
                    # routing to non-station nodes like Sink
                    ind = stationToNode[ist]
                    src_idx = ind * K + k

                    for jnd in range(I):
                        for r in range(K):
                            dst_idx = jnd * K + r

                            if src_idx < P.shape[0] and dst_idx < P.shape[1]:
                                p = P[src_idx, dst_idx]
                            else:
                                p = 0

                            if p > 0:
                                total_routing_prob += p
                                # Convert destination node to station or handle Cache
                                if isstation is not None and jnd < len(isstation) and isstation[jnd]:
                                    jst = nodeToStation[jnd]
                                    # Check capacity at destination station for the destination class
                                    # For open classes with cutoff, block arrivals when at capacity
                                    if state[jst, r] < classcap[jst, r]:
                                        transitions.append((ist, k, jst, r, rate * p))
                                elif hasattr(sn, 'nodetype') and sn.nodetype is not None and jnd < len(sn.nodetype):
                                    # Check if destination is a Cache node (non-station)
                                    nt = sn.nodetype[jnd]
                                    nt_val = int(nt.value) if hasattr(nt, 'value') else int(nt)
                                    if nt_val == 6:  # NodeType.CACHE = 6
                                        # Encode Cache node as special destination: -(100 + node_idx)
                                        transitions.append((ist, k, -(100 + jnd), r, rate * p))
                                    elif nt_val == 1:  # NodeType.SINK = 1
                                        # Job exits to Sink
                                        transitions.append((ist, k, -2, r, rate * p))
                                    elif nt_val in (7, 8):  # NodeType.ROUTER = 7, NodeType.CLASSSWITCH = 8
                                        # Follow through non-station routing nodes to find final destination
                                        resolved = _resolve_routing_through_nonstation_nodes(sn, jnd, r, P)
                                        for final_st, final_class, final_prob in resolved:
                                            # Check capacity for final destination
                                            if final_st >= 0 and final_st < M:
                                                if state[final_st, final_class] < classcap[final_st, final_class]:
                                                    transitions.append((ist, k, final_st, final_class, rate * p * final_prob))
                                            else:
                                                # Non-station destination (sink, etc.)
                                                transitions.append((ist, k, final_st, final_class, rate * p * final_prob))
                                    else:
                                        # Other non-station types - treat as sink
                                        transitions.append((ist, k, -2, r, rate * p))
                                else:
                                    # Non-station destination (like Sink) - job leaves system
                                    transitions.append((ist, k, -2, r, rate * p))
                else:
                    # Use station-indexed routing (rt)
                    src_idx = ist * K + k

                    for jst in range(M):
                        for r in range(K):
                            dst_idx = jst * K + r

                            if src_idx < P.shape[0] and dst_idx < P.shape[1]:
                                p = P[src_idx, dst_idx]
                            else:
                                p = 0

                            if p > 0:
                                total_routing_prob += p
                                # Check capacity at destination station for the destination class
                                if state[jst, r] < classcap[jst, r]:
                                    transitions.append((ist, k, jst, r, rate * p))

                    # For open networks: if routing probs sum to < 1, remaining goes to sink
                    # Use destination -2 to indicate departure from system
                    if total_routing_prob < 1.0 - 1e-6:
                        sink_prob = 1.0 - total_routing_prob
                        transitions.append((ist, k, -2, k, rate * sink_prob))

    return transitions


def _fire_transition(
    state: np.ndarray,
    transition: Tuple[int, int, int, int, float]
) -> np.ndarray:
    """
    Fire a transition and return new state.

    Args:
        state: Current state
        transition: (src_station, src_class, dst_station, dst_class, rate)
            src_st = -1 means external arrival
            dst_st = -2 means departure to sink

    Returns:
        New state after transition
    """
    new_state = state.copy()
    src_st, src_k, dst_st, dst_k, _ = transition

    if src_st >= 0:
        new_state[src_st, src_k] -= 1

    # dst_st = -2 means departure to sink (job leaves system)
    if dst_st >= 0:
        new_state[dst_st, dst_k] += 1

    return new_state


def _has_cache_nodes(sn: NetworkStruct) -> bool:
    """Check if the network has any cache nodes."""
    if not hasattr(sn, 'nodetype') or sn.nodetype is None:
        return False
    return any(nt == NodeType.CACHE for nt in sn.nodetype)


def _build_station_to_cache_map(sn: NetworkStruct) -> Dict[int, int]:
    """
    Build mapping from station index to cache node index.

    Returns:
        Dict mapping station index to cache node index
    """
    station_to_cache = {}

    if not hasattr(sn, 'nodetype') or sn.nodetype is None:
        return station_to_cache

    if not hasattr(sn, 'stationToNode') or sn.stationToNode is None:
        return station_to_cache

    for ist in range(sn.nstations):
        ind = int(sn.stationToNode[ist])
        if ind >= 0 and ind < len(sn.nodetype):
            if sn.nodetype[ind] == NodeType.CACHE:
                station_to_cache[ist] = ind

    return station_to_cache


def solver_ssa_with_cache(
    sn: NetworkStruct,
    options: Optional[SolverSSAOptions] = None,
    model = None
) -> SolverSSAReturn:
    """
    SSA solver with proper cache simulation.

    This version properly simulates cache state and hit/miss decisions
    using LRU replacement policy and Zipf access distribution.

    Args:
        sn: Network structure
        options: Solver options
        model: Optional model object for accessing cache node details

    Returns:
        SolverSSAReturn with performance metrics
    """
    start_time = time.time()

    if options is None:
        options = SolverSSAOptions()

    # Set random seed
    if options.seed > 0:
        np.random.seed(options.seed)

    M = sn.nstations
    K = sn.nclasses

    # Initialize cache states
    cache_states = _init_cache_states(sn, model)

    # Build station to cache node mapping
    station_to_cache = _build_station_to_cache_map(sn)

    # Build cache station set for quick lookup
    cache_stations = set(station_to_cache.keys())

    # Initialize state
    state = _init_state(sn, options.cutoff)

    # Initialize round-robin routing state
    rrobin_state = _init_rrobin_state(sn)

    # Initialize detailed queue states for HOL/FCFS stations
    detailed_states = _init_detailed_queue_states(sn)
    # Synchronize detailed states with initial aggregate state
    _sync_detailed_states_from_aggregate(detailed_states, state, sn)

    # Statistics accumulators
    total_time = 0.0
    time_weighted_queue = np.zeros((M, K))
    time_weighted_util = np.zeros((M, K))
    departure_counts = np.zeros((M, K))
    arrival_counts = np.zeros((M, K))

    # Get number of servers for utilization calculation
    if hasattr(sn, 'nservers') and sn.nservers is not None:
        nservers = np.asarray(sn.nservers).flatten()
    else:
        nservers = np.ones(M)

    # Get routing probabilities
    if hasattr(sn, 'rt') and sn.rt is not None:
        P = np.asarray(sn.rt)
    else:
        P = np.zeros((M * K, M * K))

    samples_collected = 0
    max_time = options.timespan[1] if np.isfinite(options.timespan[1]) else 1e6

    # Print progress helper
    last_printed = [0]  # Use list to allow modification in nested function

    def print_progress():
        if options.verbose and samples_collected > 0 and samples_collected % 100 == 0:
            if last_printed[0] == 0:
                print(f"\nSSA samples: {samples_collected:6d}", end='', flush=True)
            else:
                # Use backspaces to overwrite the number
                print(f"\b\b\b\b\b\b{samples_collected:6d}", end='', flush=True)
            last_printed[0] = samples_collected

    while samples_collected < options.samples and total_time < max_time:
        # Find enabled transitions (excluding cache-specific routing for now)
        transitions = _get_enabled_transitions(sn, state, rrobin_state, detailed_states, options.cutoff)

        if not transitions:
            if options.verbose:
                print(f"\nSSA: No enabled transitions at sample {samples_collected}")
            break

        # Total rate
        total_rate = sum(t[4] for t in transitions)

        if total_rate <= 0:
            break

        # Sample time to next event (exponential)
        dt = np.random.exponential(1.0 / total_rate)

        # Update time-weighted statistics before transition
        time_weighted_queue += dt * state

        # Note: Throughput is computed from departure_counts / total_time
        # Do NOT accumulate all enabled transition rates - only count actual departures

        # Utilization: fraction of servers busy
        for ist in range(M):
            total_at_station = np.sum(state[ist, :])
            c = nservers[ist] if ist < len(nservers) else 1
            if np.isinf(c):
                c = total_at_station if total_at_station > 0 else 1

            if hasattr(sn, 'sched') and sn.sched is not None:
                sched = sn.sched.get(ist, SchedStrategy.FCFS)
            else:
                sched = SchedStrategy.FCFS

            sched_name = sched.name if hasattr(sched, 'name') else str(sched)

            if sched_name == 'INF':
                for k in range(K):
                    time_weighted_util[ist, k] += dt * state[ist, k]
            elif sched_name in ('PSPRIO', 'DPSPRIO', 'GPSPRIO'):
                # Priority PS: if ni <= c, all get service; else only highest priority
                busy_servers = min(total_at_station, c)
                if total_at_station <= c:
                    # All jobs can be served - regular PS
                    for k in range(K):
                        if total_at_station > 0:
                            time_weighted_util[ist, k] += dt * busy_servers * (state[ist, k] / total_at_station) / c
                else:
                    # More jobs than servers - only highest priority gets credit
                    classprio = getattr(sn, 'classprio', None)
                    if classprio is not None:
                        # Find minimum (highest) priority among classes with jobs
                        min_prio_present = None
                        for kk in range(K):
                            if state[ist, kk] > 0:
                                kk_prio = classprio[kk] if kk < len(classprio) else 0
                                if min_prio_present is None or kk_prio < min_prio_present:
                                    min_prio_present = kk_prio
                        # Count jobs only in highest priority group
                        n_prio = 0
                        for kk in range(K):
                            if state[ist, kk] > 0:
                                kk_prio = classprio[kk] if kk < len(classprio) else 0
                                if kk_prio == min_prio_present:
                                    n_prio += state[ist, kk]
                        # Credit utilization only to highest priority classes
                        active_servers = min(n_prio, c)
                        for k in range(K):
                            if state[ist, k] > 0:
                                k_prio = classprio[k] if k < len(classprio) else 0
                                if k_prio == min_prio_present and n_prio > 0:
                                    time_weighted_util[ist, k] += dt * active_servers * (state[ist, k] / n_prio) / c
                    else:
                        # No priority info - fall back to regular PS
                        for k in range(K):
                            if total_at_station > 0:
                                time_weighted_util[ist, k] += dt * busy_servers * (state[ist, k] / total_at_station) / c
            elif sched_name == 'HOL':
                # HOL (Head of Line) is NON-PREEMPTIVE: servers allocated to classes in priority order
                classprio = getattr(sn, 'classprio', None)
                if classprio is not None:
                    classprio = classprio.flatten() if classprio.ndim > 1 else classprio
                    if total_at_station <= c:
                        # All jobs fit in servers - each job gets utilization credit
                        for k in range(K):
                            if state[ist, k] > 0:
                                time_weighted_util[ist, k] += dt * state[ist, k] / c
                    else:
                        # More jobs than servers - allocate servers by priority
                        for k in range(K):
                            if state[ist, k] > 0:
                                k_prio = classprio[k] if k < len(classprio) else 0
                                # Count jobs in higher priority classes
                                jobs_higher_prio = 0
                                for kk in range(K):
                                    if state[ist, kk] > 0:
                                        kk_prio = classprio[kk] if kk < len(classprio) else 0
                                        if kk_prio < k_prio:
                                            jobs_higher_prio += state[ist, kk]
                                # Servers available for this class
                                servers_available = max(0, c - jobs_higher_prio)
                                servers_for_class = min(state[ist, k], servers_available)
                                if servers_for_class > 0:
                                    time_weighted_util[ist, k] += dt * servers_for_class / c
                else:
                    # No priority info - fall back to FCFS-like behavior
                    busy_servers = min(total_at_station, c)
                    for k in range(K):
                        if total_at_station > 0:
                            time_weighted_util[ist, k] += dt * busy_servers * (state[ist, k] / total_at_station) / c
            elif sched_name == 'FCFS':
                # FCFS: use detailed state if available for accurate utilization
                if detailed_states is not None and ist in detailed_states:
                    ds = detailed_states[ist]
                    for k in range(K):
                        # Utilization = jobs in service / number of servers
                        time_weighted_util[ist, k] += dt * ds.in_service[k] / c
                else:
                    # Fallback: PS approximation
                    busy_servers = min(total_at_station, c)
                    for k in range(K):
                        if total_at_station > 0:
                            time_weighted_util[ist, k] += dt * busy_servers * (state[ist, k] / total_at_station) / c
            else:
                busy_servers = min(total_at_station, c)
                for k in range(K):
                    if total_at_station > 0:
                        time_weighted_util[ist, k] += dt * busy_servers * (state[ist, k] / total_at_station) / c

        total_time += dt

        # Select which transition fires
        rand = np.random.random() * total_rate
        cumsum = 0.0
        selected = None
        for t in transitions:
            cumsum += t[4]
            if cumsum >= rand:
                selected = t
                break

        if selected is None:
            selected = transitions[-1]

        src_st, src_k, dst_st, dst_k, _ = selected

        # Check if this transition involves arrival at a Cache node
        # In cache networks, Cache is a non-station stateful node.
        # Cache destinations are encoded as: dst_st = -(100 + cache_node_idx)
        # When a job arrives at a Cache node, we need to:
        # 1. Immediately process the READ event (sample item, check cache)
        # 2. Route to appropriate hit/miss class destination

        if dst_st <= -100:
            # Job is arriving at a Cache node (encoded as -(100 + node_idx))
            cache_node_idx = -(dst_st + 100)
            cache_state = cache_states.get(cache_node_idx)

            if cache_state is not None:
                # Process cache READ event
                # Sample which item the job is requesting
                item = cache_state.sample_item()

                # Check cache and update state (LRU)
                is_hit = cache_state.lookup_and_update(item)

                # Determine output class based on hit/miss
                if is_hit:
                    output_class = int(cache_state.hitclass[dst_k])
                else:
                    output_class = int(cache_state.missclass[dst_k])

                # Find destination for the output class from cache
                # Jobs leave cache immediately and go to next station
                # Use node-indexed routing (rtnodes) to find downstream destination
                # May need to traverse through multiple non-station nodes (e.g., Router)
                I = sn.nnodes if hasattr(sn, 'nnodes') else M
                if hasattr(sn, 'rtnodes') and sn.rtnodes is not None:
                    P_nodes = np.asarray(sn.rtnodes)
                else:
                    P_nodes = None

                if hasattr(sn, 'nodeToStation') and sn.nodeToStation is not None:
                    node_to_st = np.asarray(sn.nodeToStation).flatten()
                else:
                    node_to_st = np.arange(I, dtype=int)

                found_dest = False

                # Resolve routing through non-station nodes until we reach a station or Sink
                current_node = cache_node_idx
                current_class = output_class
                max_hops = 10  # Prevent infinite loops

                for _ in range(max_hops):
                    if P_nodes is None:
                        break

                    # Look for downstream destination from current node/class
                    src_idx = current_node * K + current_class

                    # Collect all possible destinations with their probabilities
                    candidates = []
                    total_prob = 0.0
                    for jnd in range(I):
                        for r in range(K):
                            dst_idx = jnd * K + r
                            if src_idx < P_nodes.shape[0] and dst_idx < P_nodes.shape[1]:
                                p = P_nodes[src_idx, dst_idx]
                                if p > 0:
                                    candidates.append((jnd, r, p))
                                    total_prob += p

                    if not candidates:
                        # No downstream routing found - job exits system
                        next_node = -1
                        next_class = -1
                        break

                    # Sample one destination based on probabilities
                    rand = np.random.random() * total_prob
                    cumsum = 0.0
                    next_node = -1
                    next_class = -1
                    for jnd, r, p in candidates:
                        cumsum += p
                        if rand <= cumsum:
                            next_node = jnd
                            next_class = r
                            break
                    # Fallback to last candidate if numerical issues
                    if next_node < 0:
                        next_node, next_class, _ = candidates[-1]

                    if next_node < 0:
                        break

                    # Check if next node is a station
                    if next_node < len(node_to_st) and node_to_st[next_node] >= 0:
                        jst = int(node_to_st[next_node])
                        # Record departure from source
                        if src_st >= 0:
                            departure_counts[src_st, src_k] += 1

                        # Record arrival at final destination
                        arrival_counts[jst, next_class] += 1

                        # Fire the modified transition
                        new_state = state.copy()
                        if src_st >= 0:
                            new_state[src_st, src_k] -= 1
                        new_state[jst, next_class] += 1
                        state = new_state

                        found_dest = True
                        break
                    else:
                        # Destination is another non-station node (e.g., Router)
                        # Check if it's a Sink (exits system)
                        if hasattr(sn, 'nodetype') and sn.nodetype is not None:
                            nt = sn.nodetype[next_node] if next_node < len(sn.nodetype) else None
                            if nt is not None:
                                nt_val = int(nt.value) if hasattr(nt, 'value') else int(nt)
                                if nt_val == 1:  # NodeType.SINK = 1
                                    # Job exits to Sink
                                    if src_st >= 0:
                                        departure_counts[src_st, src_k] += 1
                                    new_state = state.copy()
                                    if src_st >= 0:
                                        new_state[src_st, src_k] -= 1
                                    state = new_state
                                    found_dest = True
                                    break
                        # Continue routing through this non-station node
                        current_node = next_node
                        current_class = next_class

                if not found_dest:
                    # No station destination found - job exits system
                    if src_st >= 0:
                        departure_counts[src_st, src_k] += 1
                    new_state = state.copy()
                    if src_st >= 0:
                        new_state[src_st, src_k] -= 1
                    state = new_state
            else:
                # No cache state found, fall back to normal transition
                if src_st >= 0:
                    departure_counts[src_st, src_k] += 1
                # For Cache without state, just process as sink
                new_state = state.copy()
                if src_st >= 0:
                    new_state[src_st, src_k] -= 1
                state = new_state
        else:
            # Normal transition (not involving cache)
            if src_st >= 0:
                departure_counts[src_st, src_k] += 1
            if dst_st >= 0:
                arrival_counts[dst_st, dst_k] += 1
            state = _fire_transition(state, selected)

        # Update detailed queue states for arrivals/departures at HOL/FCFS stations
        if detailed_states:
            # Process departure from source station
            if src_st >= 0 and src_st in detailed_states:
                detailed_states[src_st].process_departure(src_k)
            # Process arrival at destination station
            if dst_st >= 0 and dst_st in detailed_states:
                detailed_states[dst_st].process_arrival(dst_k)

        # Advance round-robin counter if applicable
        if src_st >= 0 and rrobin_state:
            _advance_rrobin(rrobin_state, src_st, src_k)

        samples_collected += 1
        print_progress()

    if options.verbose:
        print()  # Newline after progress

    # Compute average metrics
    if total_time > 0:
        QN = time_weighted_queue / total_time
        UN = time_weighted_util / total_time
        # Compute throughput from actual departure counts divided by total time
        # This matches MATLAB's approach of counting actual events
        TN = departure_counts / total_time
    else:
        QN = np.zeros((M, K))
        UN = np.zeros((M, K))
        TN = np.zeros((M, K))

    # Set Source station throughput to arrival rate
    # External arrivals don't have departures counted (src_st=-1 in transitions),
    # so we need to set TN[Source] = arrival_rate explicitly
    if hasattr(sn, 'rates') and sn.rates is not None:
        rates = np.asarray(sn.rates)
        for ist in range(M):
            if hasattr(sn, 'sched') and sn.sched is not None:
                sched = sn.sched.get(ist, SchedStrategy.FCFS)
            else:
                sched = SchedStrategy.FCFS
            sched_name = sched.name if hasattr(sched, 'name') else str(sched)
            if sched_name == 'EXT' or sched == SchedStrategy.EXT:
                # Source station: throughput = arrival rate
                for k in range(K):
                    if ist < rates.shape[0] and k < rates.shape[1]:
                        TN[ist, k] = rates[ist, k]

    # Response times via Little's law: R = Q / T
    RN = np.zeros((M, K))
    for ist in range(M):
        for k in range(K):
            if TN[ist, k] > 0:
                RN[ist, k] = QN[ist, k] / TN[ist, k]

    # Cycle times and system throughput
    CN = np.sum(RN, axis=0).reshape(1, -1)
    XN = np.zeros((1, K))
    for k in range(K):
        ref_stat = int(sn.refstat[k]) if hasattr(sn, 'refstat') and k < len(sn.refstat) else 0
        if ref_stat < M:
            XN[0, k] = TN[ref_stat, k]

    # Store cache hit/miss statistics in result for debugging
    # Also update sn.nodeparam with actual hit/miss probabilities for arrival rate computation
    cache_stats = {}
    for node_idx, cs in cache_states.items():
        total_accesses = cs.hits + cs.misses
        if total_accesses > 0:
            hit_rate = cs.hits / total_accesses
            miss_rate = cs.misses / total_accesses
            cache_stats[node_idx] = {
                'hits': cs.hits,
                'misses': cs.misses,
                'hit_rate': hit_rate,
                'miss_rate': miss_rate
            }
            # Update sn.nodeparam with actual hit/miss probabilities
            # This is needed for sn_get_arvr_from_tput to compute correct arrival rates
            if hasattr(sn, 'nodeparam') and sn.nodeparam is not None and node_idx in sn.nodeparam:
                nodeparam = sn.nodeparam[node_idx]
                if nodeparam is not None and hasattr(nodeparam, 'hitclass'):
                    # Get the number of classes that use this cache
                    hitclass = getattr(nodeparam, 'hitclass', None)
                    if hitclass is not None:
                        nclasses = len(np.atleast_1d(hitclass))
                        # Create arrays of hit/miss probabilities for each originating class
                        # For now, use the same hit/miss rate for all classes
                        nodeparam.actualhitprob = np.full(nclasses, hit_rate)
                        nodeparam.actualmissprob = np.full(nclasses, miss_rate)

    # Clean up NaN values
    QN = np.nan_to_num(QN, nan=0.0)
    UN = np.nan_to_num(UN, nan=0.0)
    RN = np.nan_to_num(RN, nan=0.0)
    TN = np.nan_to_num(TN, nan=0.0)
    CN = np.nan_to_num(CN, nan=0.0)
    XN = np.nan_to_num(XN, nan=0.0)

    # Compute arrival rates from throughputs using routing matrix
    AN = sn_get_arvr_from_tput(sn, TN, None)

    result = SolverSSAReturn()
    result.Q = QN
    result.U = UN
    result.R = RN
    result.T = TN
    result.A = AN
    result.C = CN
    result.X = XN
    result.total_time = total_time
    result.runtime = time.time() - start_time
    result.method = "serial"
    result.samples = samples_collected

    # Store cache statistics for debugging
    result.cache_stats = cache_stats  # type: ignore

    return result


def solver_ssa_basic(
    sn: NetworkStruct,
    options: Optional[SolverSSAOptions] = None,
    model = None
) -> SolverSSAReturn:
    """
    Basic SSA solver using Gillespie algorithm.

    Simulates the queueing network and collects statistics for
    performance metric estimation.

    If the network contains cache nodes, delegates to solver_ssa_with_cache
    for proper cache state simulation.

    Args:
        sn: Network structure
        options: Solver options
        model: Optional model object for cache details

    Returns:
        SolverSSAReturn with performance metrics
    """
    # Check if network has cache nodes
    if _has_cache_nodes(sn):
        return solver_ssa_with_cache(sn, options, model)

    start_time = time.time()

    if options is None:
        options = SolverSSAOptions()

    # Set random seed
    if options.seed > 0:
        np.random.seed(options.seed)

    M = sn.nstations
    K = sn.nclasses

    # Initialize state
    state = _init_state(sn, options.cutoff)

    # Initialize round-robin routing state
    rrobin_state = _init_rrobin_state(sn)

    # Initialize detailed queue states for HOL/FCFS stations
    detailed_states = _init_detailed_queue_states(sn)

    # Synchronize detailed states with initial aggregate state
    # This is critical: detailed states must reflect the initial job distribution
    _sync_detailed_states_from_aggregate(detailed_states, state, sn)

    # Statistics accumulators
    total_time = 0.0
    time_weighted_queue = np.zeros((M, K))
    time_weighted_util = np.zeros((M, K))
    departure_counts = np.zeros((M, K))
    arrival_counts = np.zeros((M, K))

    # Get number of servers for utilization calculation
    if hasattr(sn, 'nservers') and sn.nservers is not None:
        nservers = np.asarray(sn.nservers).flatten()
    else:
        nservers = np.ones(M)

    samples_collected = 0
    max_time = options.timespan[1] if np.isfinite(options.timespan[1]) else 1e6

    while samples_collected < options.samples and total_time < max_time:
        # Find enabled transitions
        transitions = _get_enabled_transitions(sn, state, rrobin_state, detailed_states, options.cutoff)

        if not transitions:
            # Deadlock or absorbing state
            if options.verbose:
                print(f"SSA: No enabled transitions at sample {samples_collected}")
            break

        # Total rate
        total_rate = sum(t[4] for t in transitions)

        if total_rate <= 0:
            break

        # Sample time to next event (exponential)
        dt = np.random.exponential(1.0 / total_rate)

        # Update time-weighted statistics before transition
        time_weighted_queue += dt * state

        # Note: Throughput is computed from departure_counts / total_time
        # Do NOT accumulate all enabled transition rates - only count actual departures

        # Utilization: fraction of servers busy
        for ist in range(M):
            total_at_station = np.sum(state[ist, :])
            c = nservers[ist] if ist < len(nservers) else 1
            if np.isinf(c):
                c = total_at_station if total_at_station > 0 else 1

            if hasattr(sn, 'sched') and sn.sched is not None:
                sched = sn.sched.get(ist, SchedStrategy.FCFS)
            else:
                sched = SchedStrategy.FCFS

            # Compare by name to handle enum type mismatches
            sched_name = sched.name if hasattr(sched, 'name') else str(sched)

            if sched_name == 'INF':
                # Delay station: utilization = jobs in service
                for k in range(K):
                    time_weighted_util[ist, k] += dt * state[ist, k]
            elif sched_name in ('PSPRIO', 'DPSPRIO', 'GPSPRIO'):
                # Priority PS: if ni <= c, all get service; else only highest priority
                busy_servers = min(total_at_station, c)
                if total_at_station <= c:
                    # All jobs can be served - regular PS
                    for k in range(K):
                        if total_at_station > 0:
                            time_weighted_util[ist, k] += dt * busy_servers * (state[ist, k] / total_at_station) / c
                else:
                    # More jobs than servers - only highest priority gets credit
                    classprio = getattr(sn, 'classprio', None)
                    if classprio is not None:
                        # Find minimum (highest) priority among classes with jobs
                        min_prio_present = None
                        for kk in range(K):
                            if state[ist, kk] > 0:
                                kk_prio = classprio[kk] if kk < len(classprio) else 0
                                if min_prio_present is None or kk_prio < min_prio_present:
                                    min_prio_present = kk_prio
                        # Count jobs only in highest priority group
                        n_prio = 0
                        for kk in range(K):
                            if state[ist, kk] > 0:
                                kk_prio = classprio[kk] if kk < len(classprio) else 0
                                if kk_prio == min_prio_present:
                                    n_prio += state[ist, kk]
                        # Credit utilization only to highest priority classes
                        active_servers = min(n_prio, c)
                        for k in range(K):
                            if state[ist, k] > 0:
                                k_prio = classprio[k] if k < len(classprio) else 0
                                if k_prio == min_prio_present and n_prio > 0:
                                    time_weighted_util[ist, k] += dt * active_servers * (state[ist, k] / n_prio) / c
                    else:
                        # No priority info - fall back to regular PS
                        for k in range(K):
                            if total_at_station > 0:
                                time_weighted_util[ist, k] += dt * busy_servers * (state[ist, k] / total_at_station) / c
            elif sched_name == 'HOL':
                # HOL (Head of Line) is NON-PREEMPTIVE: use detailed state if available
                if detailed_states is not None and ist in detailed_states:
                    # Use detailed state for accurate utilization
                    ds = detailed_states[ist]
                    for k in range(K):
                        # Utilization = jobs in service / number of servers
                        time_weighted_util[ist, k] += dt * ds.in_service[k] / c
                else:
                    # Fallback: allocate servers to classes in priority order
                    classprio = getattr(sn, 'classprio', None)
                    if classprio is not None:
                        classprio = classprio.flatten() if classprio.ndim > 1 else classprio
                        if total_at_station <= c:
                            # All jobs fit in servers - each job gets utilization credit
                            for k in range(K):
                                if state[ist, k] > 0:
                                    time_weighted_util[ist, k] += dt * state[ist, k] / c
                        else:
                            # More jobs than servers - allocate servers by priority
                            for k in range(K):
                                if state[ist, k] > 0:
                                    k_prio = classprio[k] if k < len(classprio) else 0
                                    # Count jobs in higher priority classes
                                    jobs_higher_prio = 0
                                    for kk in range(K):
                                        if state[ist, kk] > 0:
                                            kk_prio = classprio[kk] if kk < len(classprio) else 0
                                            if kk_prio < k_prio:
                                                jobs_higher_prio += state[ist, kk]
                                    # Servers available for this class
                                    servers_available = max(0, c - jobs_higher_prio)
                                    servers_for_class = min(state[ist, k], servers_available)
                                    if servers_for_class > 0:
                                        time_weighted_util[ist, k] += dt * servers_for_class / c
                    else:
                        # No priority info - fall back to FCFS-like behavior
                        busy_servers = min(total_at_station, c)
                        for k in range(K):
                            if total_at_station > 0:
                                time_weighted_util[ist, k] += dt * busy_servers * (state[ist, k] / total_at_station) / c
            elif sched_name == 'FCFS':
                # FCFS: use detailed state if available for accurate utilization
                if detailed_states is not None and ist in detailed_states:
                    ds = detailed_states[ist]
                    for k in range(K):
                        # Utilization = jobs in service / number of servers
                        time_weighted_util[ist, k] += dt * ds.in_service[k] / c
                else:
                    # Fallback: PS approximation
                    busy_servers = min(total_at_station, c)
                    for k in range(K):
                        if total_at_station > 0:
                            time_weighted_util[ist, k] += dt * busy_servers * (state[ist, k] / total_at_station) / c
            else:
                # Queue: fraction of servers busy
                busy_servers = min(total_at_station, c)
                for k in range(K):
                    if total_at_station > 0:
                        time_weighted_util[ist, k] += dt * busy_servers * (state[ist, k] / total_at_station) / c

        total_time += dt

        # Select which transition fires
        rand = np.random.random() * total_rate
        cumsum = 0.0
        selected = None
        for t in transitions:
            cumsum += t[4]
            if cumsum >= rand:
                selected = t
                break

        if selected is None:
            selected = transitions[-1]

        # Record departure and arrival
        src_st, src_k, dst_st, dst_k, _ = selected
        if src_st >= 0:
            departure_counts[src_st, src_k] += 1
        if dst_st >= 0:
            arrival_counts[dst_st, dst_k] += 1

        # Fire transition
        state = _fire_transition(state, selected)

        # Update detailed queue states for arrivals/departures at HOL/FCFS stations
        if detailed_states:
            # Process departure from source station
            if src_st >= 0 and src_st in detailed_states:
                detailed_states[src_st].process_departure(src_k)
            # Process arrival at destination station
            if dst_st >= 0 and dst_st in detailed_states:
                detailed_states[dst_st].process_arrival(dst_k)

            # DEBUG: Verify detailed state invariant (total jobs match aggregate)
            for ist, ds in detailed_states.items():
                ds_total = int(np.sum(ds.in_service)) + len(ds.buffer)
                agg_total = int(np.sum(state[ist, :]))
                if ds_total != agg_total:
                    # State mismatch detected - this is a bug!
                    verbose_val = options.verbose.value if hasattr(options.verbose, 'value') else options.verbose
                    if verbose_val >= 1:
                        print(f"WARNING: State mismatch at station {ist}: detailed={ds_total}, aggregate={agg_total}")
                        print(f"  in_service={ds.in_service}, buffer={ds.buffer}")
                        print(f"  aggregate={state[ist, :]}")
                    # Try to resync - add missing jobs to buffer
                    if agg_total > ds_total:
                        # Aggregate has more jobs - add extras to buffer
                        for k in range(K):
                            detailed_count_k = ds.in_service[k] + sum(1 for x in ds.buffer if x == k)
                            agg_count_k = int(state[ist, k])
                            while detailed_count_k < agg_count_k:
                                ds.buffer.append(k)
                                detailed_count_k += 1

        # Advance round-robin counter if applicable
        if src_st >= 0 and rrobin_state:
            _advance_rrobin(rrobin_state, src_st, src_k)

        samples_collected += 1

    # Compute average metrics
    if total_time > 0:
        QN = time_weighted_queue / total_time
        UN = time_weighted_util / total_time
        # Compute throughput from actual departure counts divided by total time
        # This matches MATLAB's approach of counting actual events
        TN = departure_counts / total_time
    else:
        QN = np.zeros((M, K))
        UN = np.zeros((M, K))
        TN = np.zeros((M, K))

    # Set Source station throughput to arrival rate
    # External arrivals don't have departures counted (src_st=-1 in transitions),
    # so we need to set TN[Source] = arrival_rate explicitly
    if hasattr(sn, 'rates') and sn.rates is not None:
        rates = np.asarray(sn.rates)
        for ist in range(M):
            if hasattr(sn, 'sched') and sn.sched is not None:
                sched = sn.sched.get(ist, SchedStrategy.FCFS)
            else:
                sched = SchedStrategy.FCFS
            sched_name = sched.name if hasattr(sched, 'name') else str(sched)
            if sched_name == 'EXT' or sched == SchedStrategy.EXT:
                # Source station: throughput = arrival rate
                for k in range(K):
                    if ist < rates.shape[0] and k < rates.shape[1]:
                        TN[ist, k] = rates[ist, k]

    # Response times via Little's law: R = Q / T
    RN = np.zeros((M, K))
    for ist in range(M):
        for k in range(K):
            if TN[ist, k] > 0:
                RN[ist, k] = QN[ist, k] / TN[ist, k]

    # Cycle times and system throughput
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

    # Compute arrival rates from throughputs using routing matrix
    AN = sn_get_arvr_from_tput(sn, TN, None)

    result = SolverSSAReturn()
    result.Q = QN
    result.U = UN
    result.R = RN
    result.T = TN
    result.A = AN
    result.C = CN
    result.X = XN
    result.total_time = total_time
    result.runtime = time.time() - start_time
    result.method = "serial"
    result.samples = samples_collected

    return result


def _run_ssa_replica(args: Tuple) -> Dict:
    """
    Run a single SSA simulation replica for parallel execution.

    This function is designed to be picklable for multiprocessing.

    Args:
        args: Tuple of (sn, samples, seed, cutoff, timespan)

    Returns:
        Dict with raw metrics from this replica
    """
    sn, samples, seed, cutoff, timespan = args

    # Set seed for this replica
    np.random.seed(seed)

    M = sn.nstations
    K = sn.nclasses

    # Initialize state
    state = _init_state(sn, cutoff)

    # Statistics accumulators
    Q_accum = np.zeros((M, K))
    service_completions = np.zeros((M, K))
    arrivals = np.zeros((M, K))
    busy_time = np.zeros((M, K))
    response_time_sum = np.zeros((M, K))
    response_count = np.zeros((M, K))

    current_time = timespan[0] if timespan else 0.0
    max_time = timespan[1] if timespan and len(timespan) > 1 else float('inf')
    last_event_time = current_time

    # Get rates
    rates = _compute_rates(sn, state)
    total_rate = np.sum(rates)

    for _ in range(samples):
        if total_rate <= 0:
            break

        if current_time >= max_time:
            break

        # Time to next event
        dt = np.random.exponential(1.0 / total_rate) if total_rate > 0 else float('inf')

        if current_time + dt > max_time:
            dt = max_time - current_time

        # Accumulate time-weighted statistics
        Q_accum += state * dt
        for i in range(M):
            for r in range(K):
                if state[i, r] > 0:
                    busy_time[i, r] += dt

        current_time += dt

        if current_time >= max_time:
            break

        # Select and execute event
        flat_rates = rates.flatten()
        probs = flat_rates / total_rate
        event_idx = np.random.choice(len(probs), p=probs)
        station = event_idx // K
        job_class = event_idx % K

        # Departure from station
        if state[station, job_class] > 0:
            state[station, job_class] -= 1
            service_completions[station, job_class] += 1

            # Route to next station
            next_station = _get_next_station(sn, station, job_class)
            if next_station >= 0 and next_station < M:
                state[next_station, job_class] += 1
                arrivals[next_station, job_class] += 1

        # Update rates
        rates = _compute_rates(sn, state)
        total_rate = np.sum(rates)

        last_event_time = current_time

    # Compute metrics from this replica
    total_time = last_event_time - (timespan[0] if timespan else 0.0)
    if total_time <= 0:
        total_time = 1.0

    Q_mean = Q_accum / total_time
    U_mean = busy_time / total_time
    T_mean = service_completions / total_time

    # Response times via Little's Law: R = Q / T
    R_mean = np.zeros((M, K))
    for i in range(M):
        for r in range(K):
            if T_mean[i, r] > 0:
                R_mean[i, r] = Q_mean[i, r] / T_mean[i, r]

    return {
        'Q': Q_mean,
        'U': U_mean,
        'R': R_mean,
        'T': T_mean,
        'total_time': total_time,
        'completions': service_completions.sum()
    }


def solver_ssa_parallel(
    sn: NetworkStruct,
    options: Optional[SolverSSAOptions] = None,
    n_workers: int = None,
    n_replicas: int = 4
) -> SolverSSAReturn:
    """
    Parallel SSA solver using multiprocessing.

    Runs multiple independent simulation replicas in parallel and
    aggregates results with confidence intervals.

    Args:
        sn: Network structure
        options: Solver options
        n_workers: Number of worker processes (default: CPU count)
        n_replicas: Number of independent replicas (default: 4)

    Returns:
        SolverSSAReturn with aggregated performance metrics and CIs
    """
    import multiprocessing as mp
    from scipy import stats

    start_time = time.time()

    if options is None:
        options = SolverSSAOptions()

    if n_workers is None:
        n_workers = min(mp.cpu_count(), n_replicas)

    # Prepare arguments for each replica
    base_seed = options.seed if options.seed > 0 else int(time.time())
    samples_per_replica = max(1000, options.samples // n_replicas)

    replica_args = []
    for i in range(n_replicas):
        seed = base_seed + i * 12345
        replica_args.append((sn, samples_per_replica, seed, options.cutoff, options.timespan))

    # Run replicas in parallel
    if options.verbose:
        print(f"Running {n_replicas} SSA replicas with {n_workers} workers...")

    try:
        with mp.Pool(processes=n_workers) as pool:
            results = pool.map(_run_ssa_replica, replica_args)
    except Exception as e:
        # Fallback to serial if multiprocessing fails
        if options.verbose:
            print(f"Parallel execution failed ({e}), falling back to serial.")
        results = [_run_ssa_replica(args) for args in replica_args]

    # Aggregate results
    M = sn.nstations
    K = sn.nclasses

    Q_all = np.array([r['Q'] for r in results])
    U_all = np.array([r['U'] for r in results])
    R_all = np.array([r['R'] for r in results])
    T_all = np.array([r['T'] for r in results])

    # Mean across replicas
    Q_mean = np.mean(Q_all, axis=0)
    U_mean = np.mean(U_all, axis=0)
    R_mean = np.mean(R_all, axis=0)
    T_mean = np.mean(T_all, axis=0)

    # Confidence intervals
    alpha = 1.0 - options.confidence_level
    t_critical = stats.t.ppf(1 - alpha / 2, n_replicas - 1) if n_replicas > 1 else 1.96

    Q_std = np.std(Q_all, axis=0, ddof=1) if n_replicas > 1 else np.zeros((M, K))
    U_std = np.std(U_all, axis=0, ddof=1) if n_replicas > 1 else np.zeros((M, K))
    R_std = np.std(R_all, axis=0, ddof=1) if n_replicas > 1 else np.zeros((M, K))
    T_std = np.std(T_all, axis=0, ddof=1) if n_replicas > 1 else np.zeros((M, K))

    se_factor = t_critical / np.sqrt(n_replicas)

    Q_ci = Q_std * se_factor
    U_ci = U_std * se_factor
    R_ci = R_std * se_factor
    T_ci = T_std * se_factor

    # System throughput
    X = np.sum(T_mean, axis=0, keepdims=True)

    # Cycle times (for closed networks)
    C = np.zeros((1, K))
    if sn.njobs is not None:
        N = sn.njobs.flatten()
        for r in range(K):
            if X[0, r] > 0 and np.isfinite(N[r]) and N[r] > 0:
                C[0, r] = N[r] / X[0, r]

    total_time = sum(r['total_time'] for r in results) / n_replicas
    runtime = time.time() - start_time

    # Compute arrival rates from throughputs using routing matrix
    A_mean = sn_get_arvr_from_tput(sn, T_mean, None)

    if options.verbose:
        print(f"Parallel SSA completed in {runtime:.3f}s ({n_replicas} replicas)")

    return SolverSSAReturn(
        Q=Q_mean,
        U=U_mean,
        R=R_mean,
        T=T_mean,
        A=A_mean,
        C=C,
        X=X,
        Q_ci=Q_ci,
        U_ci=U_ci,
        R_ci=R_ci,
        T_ci=T_ci,
        total_time=total_time,
        runtime=runtime,
        method='parallel',
        samples=options.samples
    )


def solver_ssa(
    sn: NetworkStruct,
    options: Optional[SolverSSAOptions] = None,
    model = None
) -> SolverSSAReturn:
    """
    Main SSA solver handler.

    Routes to appropriate simulation method based on options.

    Args:
        sn: Network structure
        options: Solver options
        model: Optional model object for cache details

    Returns:
        SolverSSAReturn with performance metrics
    """
    if options is None:
        options = SolverSSAOptions()

    method = options.method.lower()

    if method in ['default', 'serial', 'ssa']:
        return solver_ssa_basic(sn, options, model)
    elif method in ['parallel', 'para', 'ssa.parallel']:
        return solver_ssa_parallel(sn, options)
    elif method == 'nrm':
        # Next Reaction Method - same algorithm, different selection
        return solver_ssa_basic(sn, options, model)
    else:
        # Unknown method
        if options.verbose:
            print(f"Warning: Unknown SSA method '{method}'. Using serial.")
        return solver_ssa_basic(sn, options, model)


__all__ = [
    'solver_ssa',
    'solver_ssa_basic',
    'solver_ssa_parallel',
    'solver_ssa_with_cache',
    'SolverSSAReturn',
    'SolverSSAOptions',
]
