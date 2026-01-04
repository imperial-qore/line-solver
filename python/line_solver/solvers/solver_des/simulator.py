"""
Core SimPy-based simulation engine for DES solver.

This module contains the SimPySimulator class that manages the discrete
event simulation using SimPy, including network structure parsing,
event processing, and statistics collection.
"""

import simpy
import numpy as np
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Any, Callable, Tuple
from enum import IntEnum

from .des_options import DESOptions, DESResult
from .scheduling.base import Customer
from .statistics.warmup import MSER5TransientDetector


class NodeType(IntEnum):
    """Node types in the queueing network."""
    SOURCE = 0
    SINK = 1
    QUEUE = 2
    DELAY = 3
    FORK = 4
    JOIN = 5
    ROUTER = 6
    CLASSSWITCH = 7
    CACHE = 8
    PLACE = 9
    TRANSITION = 10
    LOGGER = 11


class RoutingStrategy(IntEnum):
    """
    Routing strategies that determine how jobs are dispatched to downstream stations.

    Attributes:
        DISABLED: Routing disabled (jobs are dropped)
        RAND: Random uniform selection among destinations
        PROB: Probabilistic routing using pre-defined probabilities
        RROBIN: Round-robin cycling through destinations
        WRROBIN: Weighted round-robin with probability-based weights
        JSQ: Join Shortest Queue - route to least loaded destination
        KCHOICES: Power of K choices - sample K, pick shortest queue
        FIRING: Used for transitions in Petri net models
    """
    DISABLED = -1
    RAND = 0
    PROB = 1
    RROBIN = 2
    WRROBIN = 3
    JSQ = 4
    KCHOICES = 5
    FIRING = 6

    @staticmethod
    def from_string(obj: Any) -> 'RoutingStrategy':
        """Convert a string or JPype object to RoutingStrategy enum."""
        obj_str = str(obj).upper()
        if 'DISABLED' in obj_str:
            return RoutingStrategy.DISABLED
        elif 'RAND' in obj_str:
            return RoutingStrategy.RAND
        elif 'WRROBIN' in obj_str:
            return RoutingStrategy.WRROBIN
        elif 'RROBIN' in obj_str:
            return RoutingStrategy.RROBIN
        elif 'JSQ' in obj_str:
            return RoutingStrategy.JSQ
        elif 'KCHOICES' in obj_str:
            return RoutingStrategy.KCHOICES
        elif 'FIRING' in obj_str:
            return RoutingStrategy.FIRING
        elif 'PROB' in obj_str:
            return RoutingStrategy.PROB
        else:
            return RoutingStrategy.PROB  # Default to probabilistic


class ServerState(IntEnum):
    """Server operational states for cold start/teardown modeling."""
    OFF = 0          # Server is off (not available)
    SETUP = 1        # Server is setting up (cold start)
    ACTIVE = 2       # Server is active (can serve customers)
    DELAYOFF = 3     # Server is scheduled to turn off after delay


class HeteroSchedPolicy(IntEnum):
    """
    Scheduling policies for heterogeneous multiserver queues (native Python).

    Defines how jobs are assigned to server types in heterogeneous
    multiserver queues where different server types may have different
    speeds and class compatibilities.

    Attributes:
        ORDER: Order-based assignment (servers assigned by type order)
        ALIS: Assign to Longest Idle Server (round-robin among types)
        ALFS: Assign to Least Flexible Server (prefer specialized servers)
        FAIRNESS: Fairness-based assignment (round-robin for balanced utilization)
        FSF: Fastest Server First (prefer faster server types by rate)
        RAIS: Random Assignment with Idle Selection
    """
    ORDER = 0
    ALIS = 1
    ALFS = 2
    FAIRNESS = 3
    FSF = 4
    RAIS = 5

    @staticmethod
    def from_string(obj: Any) -> 'HeteroSchedPolicy':
        """Convert a string or JPype object to HeteroSchedPolicy enum."""
        obj_str = str(obj).upper()
        if obj_str == 'ORDER':
            return HeteroSchedPolicy.ORDER
        elif obj_str == 'ALIS':
            return HeteroSchedPolicy.ALIS
        elif obj_str == 'ALFS':
            return HeteroSchedPolicy.ALFS
        elif obj_str == 'FAIRNESS':
            return HeteroSchedPolicy.FAIRNESS
        elif obj_str == 'FSF':
            return HeteroSchedPolicy.FSF
        elif obj_str == 'RAIS':
            return HeteroSchedPolicy.RAIS
        else:
            return HeteroSchedPolicy.ORDER  # Default


@dataclass
class ServerSelection:
    """Result of server selection for heterogeneous queues."""
    server_id: int  # Global server ID (-1 if none available)
    server_type_id: int  # Server type ID (-1 if homogeneous)


@dataclass
class SimulatorConfig:
    """
    Configuration derived from NetworkStruct.

    Contains classified node indices and extracted parameters
    needed for simulation.
    """
    num_classes: int
    num_stations: int
    num_nodes: int

    # Node classification (lists of node indices)
    source_nodes: List[int] = field(default_factory=list)
    source_stations: List[int] = field(default_factory=list)
    service_nodes: List[int] = field(default_factory=list)
    service_stations: List[int] = field(default_factory=list)
    is_delay_node: List[bool] = field(default_factory=list)
    sink_nodes: List[int] = field(default_factory=list)
    fork_nodes: List[int] = field(default_factory=list)
    join_nodes: List[int] = field(default_factory=list)
    router_nodes: List[int] = field(default_factory=list)
    class_switch_nodes: List[int] = field(default_factory=list)

    # SPN nodes
    place_nodes: List[int] = field(default_factory=list)
    transition_nodes: List[int] = field(default_factory=list)

    # Parameters per service node [svc_idx]
    num_servers: Dict[int, int] = field(default_factory=dict)
    buffer_capacities: Dict[int, int] = field(default_factory=dict)
    sched_strategies: Dict[int, Any] = field(default_factory=dict)

    # Rates [svc_idx][class_id] or [src_idx][class_id]
    lambdas: Optional[np.ndarray] = None  # Arrival rates
    mus: Optional[np.ndarray] = None  # Service rates

    # Class properties
    is_open_class: List[bool] = field(default_factory=list)
    is_closed_class: List[bool] = field(default_factory=list)
    closed_class_population: List[int] = field(default_factory=list)
    reference_station: List[int] = field(default_factory=list)
    class_priorities: List[int] = field(default_factory=list)

    # Routing matrix (node*class x node*class)
    routing_matrix: Optional[np.ndarray] = None

    # Node-to-station mapping
    node_to_station: Optional[np.ndarray] = None

    # ==================== Heterogeneous Server Support ====================
    # These fields support heterogeneous multiserver queues where different
    # server types can have different service rates and class compatibility.

    # Number of server types per service node: [svc_idx] -> count (0 if homogeneous)
    num_server_types: Dict[int, int] = field(default_factory=dict)

    # Servers per type: servers_per_type[svc_idx][type_id] -> count
    servers_per_type: Dict[int, List[int]] = field(default_factory=dict)

    # Server-class compatibility: server_compat[svc_idx][type_id][class_id] -> bool
    server_compat: Dict[int, List[List[bool]]] = field(default_factory=dict)

    # Heterogeneous service rates: hetero_mus[svc_idx][type_id][class_id] -> rate
    hetero_mus: Dict[int, List[List[float]]] = field(default_factory=dict)

    # Heterogeneous scheduling policy per service node
    hetero_sched_policy: Dict[int, HeteroSchedPolicy] = field(default_factory=dict)

    # Server ID to type ID mapping: server_to_type[svc_idx][server_id] -> type_id
    server_to_type: Dict[int, List[int]] = field(default_factory=dict)

    # ==================== Routing Strategy Support ====================
    # These fields support non-probabilistic routing strategies (RAND, RROBIN, etc.)

    # Routing strategy per node per class: [node_idx][class_idx] -> RoutingStrategy
    node_routing_strategies: Dict[int, Dict[int, RoutingStrategy]] = field(default_factory=dict)

    # K parameter for Power-of-K-Choices routing
    kchoices_k: int = 2

    # ==================== Load-Dependent Service Support ====================
    # These fields support state-dependent service rates where the service rate
    # depends on the number of jobs at the station.

    # Load-dependent scaling factors: lld_scaling[svc_idx][n-1] = scaling factor when n jobs present
    lld_scaling: Dict[int, np.ndarray] = field(default_factory=dict)

    # Flags for which stations have load-dependent service
    is_load_dependent: Dict[int, bool] = field(default_factory=dict)

    # ==================== Cache Node Support ====================
    # These fields support cache nodes with LRU/FIFO/RR replacement policies

    # Cache nodes list
    cache_nodes: List[int] = field(default_factory=list)

    # Cache configuration per node: cache_config[node_idx] -> CacheConfig
    cache_config: Dict[int, 'CacheConfig'] = field(default_factory=dict)


class ReplacementPolicy(IntEnum):
    """Cache replacement policies."""
    LRU = 0    # Least Recently Used
    FIFO = 1   # First In First Out
    RR = 2     # Random Replacement


@dataclass
class CacheConfig:
    """
    Configuration for a cache node.

    Attributes:
        node_idx: Node index
        num_items: Number of cacheable items
        capacity: Cache capacity (number of items that can be stored)
        replacement_policy: LRU, FIFO, or RR
        hit_class: Class to switch to on cache hit (per input class)
        miss_class: Class to switch to on cache miss (per input class)
        access_probs: Probability of accessing each item (for Zipf-like patterns)
    """
    node_idx: int = 0
    num_items: int = 100
    capacity: int = 10
    replacement_policy: ReplacementPolicy = ReplacementPolicy.LRU
    hit_class: Optional[Dict[int, int]] = None
    miss_class: Optional[Dict[int, int]] = None
    access_probs: Optional[np.ndarray] = None


@dataclass
class CacheState:
    """
    Runtime state for a cache node.

    Attributes:
        cache_contents: List of item IDs currently in cache (ordered by access time for LRU)
        total_hits: Total cache hits per class
        total_misses: Total cache misses per class
    """
    cache_contents: List[int] = field(default_factory=list)
    total_hits: Dict[int, int] = field(default_factory=dict)
    total_misses: Dict[int, int] = field(default_factory=dict)


@dataclass
class StationStats:
    """Statistics for a single service station."""
    num_classes: int

    # Time-weighted queue length
    total_queue_time: np.ndarray = field(init=False)
    last_queue_update_time: np.ndarray = field(init=False)
    current_queue_length: np.ndarray = field(init=False)

    # Time-weighted utilization
    total_busy_time: np.ndarray = field(init=False)
    last_busy_update_time: np.ndarray = field(init=False)
    current_busy_servers: np.ndarray = field(init=False)

    # Completions and response times
    completed_customers: np.ndarray = field(init=False)
    arrived_customers: np.ndarray = field(init=False)
    dropped_customers: np.ndarray = field(init=False)
    response_time_sum: np.ndarray = field(init=False)
    response_time_count: np.ndarray = field(init=False)

    # Individual response time samples for CDF computation
    response_time_samples: Dict[int, List[float]] = field(init=False)

    def __post_init__(self):
        K = self.num_classes
        self.total_queue_time = np.zeros(K)
        self.last_queue_update_time = np.zeros(K)
        self.current_queue_length = np.zeros(K, dtype=int)
        self.total_busy_time = np.zeros(K)
        self.last_busy_update_time = np.zeros(K)
        self.current_busy_servers = np.zeros(K, dtype=int)
        self.completed_customers = np.zeros(K, dtype=int)
        self.arrived_customers = np.zeros(K, dtype=int)
        self.dropped_customers = np.zeros(K, dtype=int)
        self.response_time_sum = np.zeros(K)
        self.response_time_count = np.zeros(K, dtype=int)
        # Initialize response time sample lists for each class
        self.response_time_samples = {k: [] for k in range(K)}

    def update_queue(self, class_id: int, current_time: float, delta: int = 0):
        """Update time-weighted queue length."""
        elapsed = current_time - self.last_queue_update_time[class_id]
        if elapsed > 0:
            self.total_queue_time[class_id] += elapsed * self.current_queue_length[class_id]
        self.current_queue_length[class_id] += delta
        self.last_queue_update_time[class_id] = current_time

    def update_busy(self, class_id: int, current_time: float, delta: int = 0):
        """Update time-weighted busy servers."""
        elapsed = current_time - self.last_busy_update_time[class_id]
        if elapsed > 0:
            self.total_busy_time[class_id] += elapsed * self.current_busy_servers[class_id]
        self.current_busy_servers[class_id] += delta
        self.last_busy_update_time[class_id] = current_time

    def record_arrival(self, class_id: int):
        """Record an arrival."""
        self.arrived_customers[class_id] += 1

    def record_drop(self, class_id: int):
        """Record a dropped customer."""
        self.dropped_customers[class_id] += 1

    def record_completion(self, class_id: int, response_time: float):
        """Record a service completion."""
        self.completed_customers[class_id] += 1
        self.response_time_sum[class_id] += response_time
        self.response_time_count[class_id] += 1
        # Store individual sample for CDF computation
        self.response_time_samples[class_id].append(response_time)

    def get_avg_queue_length(self, class_id: int, total_time: float) -> float:
        """Get average queue length for class."""
        if total_time <= 0:
            return 0.0
        return self.total_queue_time[class_id] / total_time

    def get_utilization(self, class_id: int, total_time: float, num_servers: int) -> float:
        """Get utilization for class."""
        if total_time <= 0 or num_servers <= 0:
            return 0.0
        return self.total_busy_time[class_id] / (total_time * num_servers)

    def get_throughput(self, class_id: int, total_time: float) -> float:
        """Get throughput for class."""
        if total_time <= 0:
            return 0.0
        return self.completed_customers[class_id] / total_time

    def get_avg_response_time(self, class_id: int) -> float:
        """Get average response time for class."""
        if self.response_time_count[class_id] <= 0:
            return 0.0
        return self.response_time_sum[class_id] / self.response_time_count[class_id]


@dataclass
class InServiceFiring:
    """
    Tracks an in-progress firing with phase information.

    Used for phase-level race semantics in SPN transitions.
    """
    firing_id: int
    trans_idx: int
    mode: Any  # TransitionModeInfo
    current_phase: int
    next_event_time: float
    consumed_tokens: Dict[int, np.ndarray]  # place_idx -> tokens consumed


@dataclass
class TransitionModeInfo:
    """
    Information about a transition mode in an SPN.

    Mirrors TransitionModeInfo from Kotlin implementation.
    """
    mode_idx: int
    mode_name: str
    timing_strategy: str  # 'TIMED' or 'IMMEDIATE'
    priority: int
    weight: float
    num_servers: int  # Max concurrent firings (use float('inf') for unlimited)
    enabling_conditions: np.ndarray  # [num_places, num_classes] - required tokens
    inhibiting_conditions: np.ndarray  # [num_places, num_classes] - max tokens threshold
    firing_outcomes: np.ndarray  # [num_nodes, num_classes] - tokens produced
    firing_distribution: Optional[Any] = None  # Raw distribution object
    # Distribution parameters for proper sampling
    firing_proc_id: str = 'EXP'  # EXP, ERLANG, HYPEREXP, PH, COXIAN, etc.
    firing_mean: float = 1.0
    firing_phases: int = 1
    firing_pie: Optional[np.ndarray] = None  # Initial phase probabilities
    firing_d0: Optional[np.ndarray] = None  # Sub-generator matrix
    firing_d1: Optional[np.ndarray] = None  # Transition matrix


class SimPySimulator:
    """
    Core DES engine using SimPy.

    Manages the simulation environment, processes, and statistics collection.
    Architecture for discrete event simulation.

    Attributes:
        env: SimPy simulation environment
        config: Parsed network configuration
        stats: Per-station statistics
        options: DES solver options
    """

    def __init__(
        self,
        sn: Any,
        options: DESOptions,
        init_sol: Optional[np.ndarray] = None,
        model: Optional[Any] = None,
    ):
        """
        Initialize simulator from network structure.

        Args:
            sn: NetworkStruct (from model.getStruct())
            options: DESOptions configuration
            init_sol: Optional initial state matrix [1 x (M*K)]
            model: Optional original model for accessing node-level info
        """
        self.sn = sn
        self.options = options
        self.init_sol = init_sol
        self.model = model  # Store model for accessing Transition distributions

        # SimPy environment
        self.env = simpy.Environment()

        # Random number generator
        seed = options.seed if options.seed > 0 else None
        self.rng = np.random.default_rng(seed)

        # Initialize SPN state before parsing (will be populated during parsing)
        self.transition_modes: Dict[int, List[TransitionModeInfo]] = {}
        self.place_node_to_idx: Dict[int, int] = {}
        self.transition_node_to_idx: Dict[int, int] = {}

        # Parse network structure
        self.config = self._parse_network_structure()

        # Per-station statistics
        self.stats: Dict[int, StationStats] = {}
        for svc_idx in range(len(self.config.service_nodes)):
            self.stats[svc_idx] = StationStats(self.config.num_classes)

        # Phase 3-4: Scheduler infrastructure (job-based and preemptive scheduling)
        from .scheduling.factory import create_scheduler
        self.schedulers: Dict[int, Any] = {}  # [svc_idx] -> SchedulingStrategy instance

        # Phase 3: Initialize schedulers for each service node
        self._initialize_schedulers()

        # Phase 5: Initialize Fork/Join manager
        from .nodes.fork_join import ForkJoinManager
        self.fj_manager = ForkJoinManager()
        self.fork_tracking: Dict[int, Any] = {}  # [parent_job_id] -> ForkJobInfo
        self._parse_fork_join_params()

        # Initialize Join buffer statistics
        self.join_buffer_stats: Dict[int, Any] = {}
        for join_idx, join_node_idx in enumerate(self.config.join_nodes):
            self.join_buffer_stats[join_idx] = StationStats(self.config.num_classes)

        # System-level statistics
        self.system_completed = np.zeros(self.config.num_classes, dtype=int)
        self.system_response_time_sum = np.zeros(self.config.num_classes)

        # Phase 4: Process tracking for preemption
        self.active_processes: Dict[int, Dict[int, Any]] = {}  # [svc_idx][server_id] -> SimPy process
        self.active_customers: Dict[int, Dict[int, Customer]] = {}  # [svc_idx][server_id] -> Customer

        # Queue state per service node
        self.wait_queues: Dict[int, List[Customer]] = {}
        self.server_busy: Dict[int, List[bool]] = {}
        for svc_idx in range(len(self.config.service_nodes)):
            self.wait_queues[svc_idx] = []
            num_servers = self.config.num_servers.get(svc_idx, 1)
            if num_servers < 1000000:  # Finite servers
                self.server_busy[svc_idx] = [False] * num_servers
            else:
                self.server_busy[svc_idx] = []  # Infinite servers

            # Initialize preemption tracking
            self.active_processes[svc_idx] = {}
            self.active_customers[svc_idx] = {}

        # Phase 6b: Server state management for cold start/teardown
        # Track server state for each server (OFF, SETUP, ACTIVE, DELAYOFF)
        self.server_states: Dict[int, Dict[int, ServerState]] = {}  # [svc_idx][server_id] -> ServerState
        self.server_setup_times: Dict[int, float] = {}  # [svc_idx] -> setup time distribution
        self.server_delayoff_times: Dict[int, float] = {}  # [svc_idx] -> delay-off time distribution
        self.delayoff_events: Dict[Tuple[int, int], Any] = {}  # [(svc_idx, server_id)] -> SimPy event

        for svc_idx in range(len(self.config.service_nodes)):
            num_servers = self.config.num_servers.get(svc_idx, 1)
            if num_servers < 1000000:  # Finite servers only
                self.server_states[svc_idx] = {}
                for server_id in range(num_servers):
                    # All servers start in ACTIVE state (ready to serve)
                    self.server_states[svc_idx][server_id] = ServerState.ACTIVE
            # Setup and delayoff times would be configured per service node (deferred to extension)

        # ==================== Heterogeneous Server State ====================
        # Runtime state for heterogeneous multiserver queues

        # Busy count per server type: busy_count_per_type[svc_idx][type_id] -> count
        self.busy_count_per_type: Dict[int, List[int]] = {}

        # Server type order for ALIS/FAIRNESS round-robin: [svc_idx] -> list of type indices
        self.server_type_order: Dict[int, List[int]] = {}

        # ALFS order (sorted by flexibility - least compatible classes first)
        self.alfs_order: Dict[int, List[int]] = {}

        # Initialize heterogeneous server state from config
        self._initialize_heterogeneous_servers()

        # ==================== Routing Strategy State ====================
        # Round-robin counters for RROBIN/WRROBIN routing
        self.round_robin_counters: Dict[int, int] = {}  # [node_idx] -> counter
        for i in range(self.config.num_nodes):
            self.round_robin_counters[i] = 0

        # K parameter for Power-of-K-Choices routing
        self.kchoices_k = self.config.kchoices_k

        # ==================== Cache State ====================
        # Cache state per node: cache_states[node_idx] -> CacheState
        self.cache_states: Dict[int, CacheState] = {}
        for node_idx in self.config.cache_nodes:
            self.cache_states[node_idx] = CacheState()
            # Initialize hit/miss counters for all classes
            for k in range(self.config.num_classes):
                self.cache_states[node_idx].total_hits[k] = 0
                self.cache_states[node_idx].total_misses[k] = 0

        # Simulation state
        self.warmup_done = False
        self.warmup_end_time = 0.0
        self.total_event_count = 0
        self.next_job_id = 0

        # MSER-5 warmup detection
        self.mser_enabled = (options.tranfilter.lower() == "mser5")
        self.mser_batch_size = options.mserbatch
        self.mser_observation_interval = 0  # Set during simulation
        self.last_mser_event_count = 0
        self.mser_observations: Dict[int, Dict[int, List[float]]] = {}  # [svc_idx][class_id] -> list
        self.mser_observation_times: List[float] = []
        self.last_mser_sample_time = 0.0

        # Initialize MSER observation storage for each service node/class
        for svc_idx in range(len(self.config.service_nodes)):
            self.mser_observations[svc_idx] = {}
            for k in range(self.config.num_classes):
                self.mser_observations[svc_idx][k] = []

        # Track last queue time for interval-based averaging
        self.last_mser_queue_time: Dict[int, np.ndarray] = {}
        for svc_idx in range(len(self.config.service_nodes)):
            self.last_mser_queue_time[svc_idx] = np.zeros(self.config.num_classes)

        # SPN state (place tokens and transition in-service counts)
        self.place_tokens: Dict[int, np.ndarray] = {}
        self.transition_in_service: Dict[int, Dict[int, int]] = {}
        self.place_token_time: Dict[int, np.ndarray] = {}
        self.last_place_update: Dict[int, float] = {}
        self.place_completions: Dict[int, np.ndarray] = {}
        # Track in-service tokens per transition per class (for closed network QLen)
        self.transition_in_service_time: Dict[int, np.ndarray] = {}
        self.transition_in_service_count: Dict[int, np.ndarray] = {}  # Current count per class
        self.last_transition_update: Dict[int, float] = {}

        # Transient analysis state
        self.is_transient_mode = options.is_transient_mode()
        self.transient_interval = options.transient_interval
        self.transient_max_samples = options.transient_max_samples
        self.transient_end_time = float('inf')
        if options.timespan is not None and len(options.timespan) >= 2:
            self.transient_end_time = float(options.timespan[1])

        # Transient data storage: list of (time, metrics) tuples
        self.transient_times: List[float] = []
        # transient_metrics[svc_idx][class_id] = list of (queue_length, utilization, throughput) tuples
        self.transient_queue_lengths: Dict[int, Dict[int, List[float]]] = {}
        self.transient_utilizations: Dict[int, Dict[int, List[float]]] = {}
        self.transient_throughputs: Dict[int, Dict[int, List[float]]] = {}
        self.last_transient_sample_time = 0.0

        # Initialize transient storage for each service node/class
        if self.is_transient_mode:
            for svc_idx in range(len(self.config.service_nodes)):
                self.transient_queue_lengths[svc_idx] = {}
                self.transient_utilizations[svc_idx] = {}
                self.transient_throughputs[svc_idx] = {}
                for k in range(self.config.num_classes):
                    self.transient_queue_lengths[svc_idx][k] = []
                    self.transient_utilizations[svc_idx][k] = []
                    self.transient_throughputs[svc_idx][k] = []

    def _parse_network_structure(self) -> SimulatorConfig:
        """Parse NetworkStruct into SimulatorConfig."""
        sn = self.sn

        # Handle both native Python sn and JPype-wrapped sn
        if hasattr(sn, 'nclasses'):
            num_classes = int(sn.nclasses)
        else:
            num_classes = 1

        if hasattr(sn, 'nstations'):
            num_stations = int(sn.nstations)
        else:
            num_stations = 1

        if hasattr(sn, 'nnodes'):
            num_nodes = int(sn.nnodes)
        else:
            num_nodes = num_stations

        config = SimulatorConfig(
            num_classes=num_classes,
            num_stations=num_stations,
            num_nodes=num_nodes,
        )

        # Get node types
        nodetype = None
        if hasattr(sn, 'nodetype'):
            nodetype = sn.nodetype

        # Classify nodes by type
        if nodetype is not None:
            for i in range(num_nodes):
                node_type = self._get_node_type(nodetype, i)

                if node_type == NodeType.SOURCE:
                    config.source_nodes.append(i)
                    if hasattr(sn, 'nodeToStation'):
                        station_idx = self._get_int_value(sn.nodeToStation, i)
                        config.source_stations.append(station_idx if station_idx is not None else i)
                elif node_type == NodeType.QUEUE:
                    config.service_nodes.append(i)
                    config.is_delay_node.append(False)
                    if hasattr(sn, 'nodeToStation'):
                        station_idx = self._get_int_value(sn.nodeToStation, i)
                        config.service_stations.append(station_idx if station_idx is not None else i)
                elif node_type == NodeType.DELAY:
                    config.service_nodes.append(i)
                    config.is_delay_node.append(True)
                    if hasattr(sn, 'nodeToStation'):
                        station_idx = self._get_int_value(sn.nodeToStation, i)
                        config.service_stations.append(station_idx if station_idx is not None else i)
                elif node_type == NodeType.SINK:
                    config.sink_nodes.append(i)
                elif node_type == NodeType.FORK:
                    config.fork_nodes.append(i)
                elif node_type == NodeType.JOIN:
                    config.join_nodes.append(i)
                elif node_type == NodeType.ROUTER:
                    config.router_nodes.append(i)
                elif node_type == NodeType.CACHE:
                    config.cache_nodes.append(i)
                elif node_type == NodeType.CLASSSWITCH:
                    config.class_switch_nodes.append(i)
                elif node_type == NodeType.PLACE:
                    config.place_nodes.append(i)
                elif node_type == NodeType.TRANSITION:
                    config.transition_nodes.append(i)

        # Parse class properties
        for k in range(num_classes):
            # Open vs closed class
            is_open = True
            if hasattr(sn, 'njobs') and sn.njobs is not None:
                njobs_k = self._get_value(sn.njobs, k)
                is_open = np.isinf(njobs_k) if njobs_k is not None else True

            config.is_open_class.append(is_open)
            config.is_closed_class.append(not is_open)

            # Closed class population
            pop = 0
            if not is_open and hasattr(sn, 'njobs'):
                njobs_k = self._get_value(sn.njobs, k)
                pop = int(njobs_k) if njobs_k is not None and np.isfinite(njobs_k) else 0
            config.closed_class_population.append(pop)

            # Reference station
            ref = 0
            if hasattr(sn, 'refstat') and sn.refstat is not None:
                ref_k = self._get_value(sn.refstat, k)
                ref = int(ref_k) if ref_k is not None else 0
            config.reference_station.append(ref)

            # Class priority
            prio = 0
            if hasattr(sn, 'classprio') and sn.classprio is not None:
                prio_k = self._get_value(sn.classprio, k)
                prio = int(prio_k) if prio_k is not None else 0
            config.class_priorities.append(prio)

        # Extract rates
        num_sources = len(config.source_nodes)
        num_service = len(config.service_nodes)

        # Arrival rates [src_idx, class_id]
        config.lambdas = np.zeros((num_sources, num_classes))
        if hasattr(sn, 'rates') and sn.rates is not None:
            for src_idx, station_idx in enumerate(config.source_stations):
                for k in range(num_classes):
                    rate = self._get_matrix_value(sn.rates, station_idx, k)
                    if rate is not None and rate > 0 and np.isfinite(rate):
                        config.lambdas[src_idx, k] = rate

        # Service rates [svc_idx, class_id]
        config.mus = np.zeros((num_service, num_classes))
        if hasattr(sn, 'rates') and sn.rates is not None:
            for svc_idx, station_idx in enumerate(config.service_stations):
                for k in range(num_classes):
                    rate = self._get_matrix_value(sn.rates, station_idx, k)
                    if rate is not None and rate > 0 and np.isfinite(rate):
                        config.mus[svc_idx, k] = rate

        # Extract number of servers
        for svc_idx in range(num_service):
            if config.is_delay_node[svc_idx]:
                config.num_servers[svc_idx] = 10000000  # Effectively infinite
            elif hasattr(sn, 'nservers') and sn.nservers is not None:
                station_idx = config.service_stations[svc_idx]
                nservers = self._get_value(sn.nservers, station_idx)
                if nservers is None or (np.isscalar(nservers) and np.isinf(nservers)):
                    config.num_servers[svc_idx] = 10000000  # Effectively infinite
                else:
                    config.num_servers[svc_idx] = int(nservers) if nservers else 1
            else:
                config.num_servers[svc_idx] = 1

        # Extract buffer capacities
        for svc_idx in range(num_service):
            station_idx = config.service_stations[svc_idx]
            if hasattr(sn, 'cap') and sn.cap is not None:
                cap = self._get_value(sn.cap, station_idx)
                if cap is not None and cap > 0 and np.isfinite(cap):
                    config.buffer_capacities[svc_idx] = int(cap)
                else:
                    config.buffer_capacities[svc_idx] = 10000000
            else:
                config.buffer_capacities[svc_idx] = 10000000

        # Extract routing matrix
        if hasattr(sn, 'rtnodes') and sn.rtnodes is not None:
            config.routing_matrix = np.asarray(sn.rtnodes)

        # Node-to-station mapping
        if hasattr(sn, 'nodeToStation') and sn.nodeToStation is not None:
            config.node_to_station = np.asarray(sn.nodeToStation)

        # Parse routing strategies from sn.routing
        self._parse_routing_strategies(sn, config)

        # Parse load-dependent service rates
        self._parse_load_dependent_service(sn, config)

        # Parse heterogeneous server configuration
        self._parse_heterogeneous_servers(sn, config)

        # Parse cache node configuration
        if config.cache_nodes:
            self._parse_cache_nodes(sn, config)

        # Parse SPN-specific parameters if places/transitions exist
        if config.place_nodes or config.transition_nodes:
            self._parse_spn_parameters(sn, config)

        return config

    def _parse_routing_strategies(self, sn: Any, config: SimulatorConfig) -> None:
        """
        Parse routing strategies from NetworkStruct.

        Extracts per-node, per-class routing strategies from sn.routing.
        This supports RAND, RROBIN, WRROBIN, JSQ, KCHOICES strategies
        in addition to the default probabilistic (PROB) routing.

        Args:
            sn: NetworkStruct (from model.getStruct())
            config: SimulatorConfig to populate
        """
        num_nodes = config.num_nodes
        num_classes = config.num_classes

        # Check if routing map exists
        if not hasattr(sn, 'routing') or sn.routing is None:
            return

        # sn.routing is a routing mapping
        # We need to iterate over nodes and jobclasses
        for node_idx in range(num_nodes):
            config.node_routing_strategies[node_idx] = {}

            # Try to get routing for this node
            node = None
            if hasattr(sn, 'nodes') and sn.nodes is not None:
                try:
                    node = sn.nodes[node_idx]
                except (IndexError, TypeError):
                    continue

            if node is None:
                continue

            try:
                # Try to get the node's routing map
                node_routing = None
                if hasattr(sn.routing, 'get'):
                    node_routing = sn.routing.get(node)
                elif hasattr(sn.routing, 'containsKey'):
                    if sn.routing.containsKey(node):
                        node_routing = sn.routing.get(node)

                if node_routing is None:
                    continue

                # For each job class, get the routing strategy
                for class_idx in range(num_classes):
                    jobclass = None
                    if hasattr(sn, 'jobclasses') and sn.jobclasses is not None:
                        try:
                            jobclass = sn.jobclasses[class_idx]
                        except (IndexError, TypeError):
                            continue

                    if jobclass is None:
                        continue

                    try:
                        strategy = None
                        if hasattr(node_routing, 'get'):
                            strategy = node_routing.get(jobclass)
                        elif hasattr(node_routing, 'containsKey'):
                            if node_routing.containsKey(jobclass):
                                strategy = node_routing.get(jobclass)

                        if strategy is not None:
                            config.node_routing_strategies[node_idx][class_idx] = \
                                RoutingStrategy.from_string(str(strategy))
                    except (AttributeError, TypeError):
                        pass

            except (AttributeError, TypeError):
                continue

    def _parse_load_dependent_service(self, sn: Any, config: SimulatorConfig) -> None:
        """
        Parse load-dependent service rates from NetworkStruct.

        Reads sn.lldscaling matrix and sets up per-station scaling factors.
        lldscaling[station, n-1] gives the scaling factor when n jobs are at that station.
        The scaling factor multiplies the service rate.

        Args:
            sn: NetworkStruct (from model.getStruct())
            config: SimulatorConfig to populate
        """
        num_service = len(config.service_nodes)

        # Initialize all stations as not load-dependent
        for svc_idx in range(num_service):
            config.is_load_dependent[svc_idx] = False

        # Check if lldscaling exists
        if not hasattr(sn, 'lldscaling') or sn.lldscaling is None:
            return

        try:
            # Convert to numpy array
            lldscaling = None
            if hasattr(sn.lldscaling, 'isEmpty'):
                if sn.lldscaling.isEmpty():
                    return
                # Convert to numpy
                lldscaling = self._java_matrix_to_numpy(sn.lldscaling)
            elif hasattr(sn.lldscaling, '__len__'):
                lldscaling = np.asarray(sn.lldscaling)
            else:
                lldscaling = np.asarray(sn.lldscaling)

            if lldscaling is None or lldscaling.size == 0:
                return

            # lldscaling shape should be [num_stations, max_population]
            # lldscaling[station, n-1] = scaling factor when n jobs present
            for svc_idx in range(num_service):
                station_idx = config.service_stations[svc_idx]
                if station_idx >= lldscaling.shape[0]:
                    continue

                max_n = lldscaling.shape[1] if lldscaling.ndim > 1 else 1

                # Extract scaling factors for this station
                if lldscaling.ndim > 1:
                    station_scaling = lldscaling[station_idx, :]
                else:
                    station_scaling = np.array([lldscaling[station_idx]])

                # Check if this station has any non-trivial scaling (not all 1s)
                has_scaling = False
                for n in range(max_n):
                    scale = station_scaling[n] if n < len(station_scaling) else 1.0
                    if scale != 1.0 and scale > 0:
                        has_scaling = True
                        break

                if has_scaling:
                    config.is_load_dependent[svc_idx] = True
                    # Store scaling factors (replace 0s with 1s)
                    scaling = np.array([
                        scale if scale > 0 else 1.0
                        for scale in station_scaling
                    ])
                    config.lld_scaling[svc_idx] = scaling

        except Exception:
            # Silently ignore parsing errors
            pass

    def _parse_cache_nodes(self, sn: Any, config: SimulatorConfig) -> None:
        """
        Parse cache node configuration from NetworkStruct.

        Extracts cache parameters from sn.nodeparam for each cache node:
        - Number of items (nitems)
        - Item capacity (itemcap)
        - Replacement strategy (replacestrat)
        - Hit/miss class mappings (hitclass, missclass)
        - Access probabilities (pread)

        Args:
            sn: NetworkStruct (from model.getStruct())
            config: SimulatorConfig to populate
        """
        num_classes = config.num_classes

        # Check if nodeparam exists
        if not hasattr(sn, 'nodeparam') or sn.nodeparam is None:
            return

        # Get nodes list for nodeparam lookup
        nodes = None
        if hasattr(sn, 'nodes') and sn.nodes is not None:
            nodes = sn.nodes

        for cache_node_idx in config.cache_nodes:
            try:
                # Get cache node parameters
                cache_param = None
                if nodes is not None:
                    cache_node = None
                    try:
                        if hasattr(nodes, '__getitem__'):
                            cache_node = nodes[cache_node_idx]
                        elif hasattr(nodes, 'get'):
                            cache_node = nodes.get(cache_node_idx)
                    except (IndexError, TypeError, KeyError):
                        pass

                    if cache_node is not None:
                        if hasattr(sn.nodeparam, 'get'):
                            cache_param = sn.nodeparam.get(cache_node)
                        elif hasattr(sn.nodeparam, '__getitem__'):
                            cache_param = sn.nodeparam[cache_node]

                if cache_param is None:
                    # Default cache config with minimal settings
                    cache_config = CacheConfig(
                        node_idx=cache_node_idx,
                        num_items=100,
                        capacity=10,
                        replacement_policy=ReplacementPolicy.LRU,
                        hit_class={k: k for k in range(num_classes)},
                        miss_class={k: k for k in range(num_classes)},
                        access_probs=None
                    )
                    config.cache_config[cache_node_idx] = cache_config
                    continue

                # Extract number of items
                num_items = 100  # default
                if hasattr(cache_param, 'nitems'):
                    num_items = int(cache_param.nitems) if cache_param.nitems else 100

                # Extract capacity (total across all levels)
                capacity = 10  # default
                if hasattr(cache_param, 'itemcap') and cache_param.itemcap is not None:
                    try:
                        if hasattr(cache_param.itemcap, 'sumRows'):
                            capacity = int(cache_param.itemcap.sumRows())
                        elif hasattr(cache_param.itemcap, '__iter__'):
                            capacity = int(sum(cache_param.itemcap))
                        else:
                            capacity = int(cache_param.itemcap)
                    except (TypeError, ValueError):
                        capacity = 10

                # Extract replacement strategy
                replacement = ReplacementPolicy.LRU  # default
                if hasattr(cache_param, 'replacestrat') and cache_param.replacestrat is not None:
                    strat_str = str(cache_param.replacestrat).upper()
                    if 'FIFO' in strat_str:
                        replacement = ReplacementPolicy.FIFO
                    elif 'RR' in strat_str or 'RANDOM' in strat_str:
                        replacement = ReplacementPolicy.RR
                    else:
                        replacement = ReplacementPolicy.LRU

                # Extract hit/miss class mappings
                hit_class = {k: k for k in range(num_classes)}  # default: same class
                miss_class = {k: k for k in range(num_classes)}  # default: same class

                if hasattr(cache_param, 'hitclass') and cache_param.hitclass is not None:
                    try:
                        hit_matrix = cache_param.hitclass
                        for k in range(num_classes):
                            if hasattr(hit_matrix, 'get'):
                                val = hit_matrix.get(0, k)
                                if val >= 0:
                                    hit_class[k] = int(val)
                            elif hasattr(hit_matrix, '__getitem__'):
                                val = hit_matrix[0, k] if hasattr(hit_matrix, 'shape') else hit_matrix[k]
                                if val >= 0:
                                    hit_class[k] = int(val)
                    except (IndexError, TypeError, KeyError):
                        pass

                if hasattr(cache_param, 'missclass') and cache_param.missclass is not None:
                    try:
                        miss_matrix = cache_param.missclass
                        for k in range(num_classes):
                            if hasattr(miss_matrix, 'get'):
                                val = miss_matrix.get(0, k)
                                if val >= 0:
                                    miss_class[k] = int(val)
                            elif hasattr(miss_matrix, '__getitem__'):
                                val = miss_matrix[0, k] if hasattr(miss_matrix, 'shape') else miss_matrix[k]
                                if val >= 0:
                                    miss_class[k] = int(val)
                    except (IndexError, TypeError, KeyError):
                        pass

                # Extract access probabilities
                access_probs = None
                if hasattr(cache_param, 'pread') and cache_param.pread is not None:
                    try:
                        pread = cache_param.pread
                        # pread is a mapping of priority read values
                        if hasattr(pread, 'get'):
                            # Get first server's probabilities
                            probs_list = pread.get(0)
                            if probs_list is not None:
                                access_probs = np.array(list(probs_list))
                        elif hasattr(pread, '__iter__'):
                            access_probs = np.array(list(pread))
                    except (TypeError, ValueError):
                        access_probs = None

                cache_config = CacheConfig(
                    node_idx=cache_node_idx,
                    num_items=num_items,
                    capacity=capacity,
                    replacement_policy=replacement,
                    hit_class=hit_class,
                    miss_class=miss_class,
                    access_probs=access_probs
                )
                config.cache_config[cache_node_idx] = cache_config

            except Exception:
                # Fallback to default config
                cache_config = CacheConfig(
                    node_idx=cache_node_idx,
                    num_items=100,
                    capacity=10,
                    replacement_policy=ReplacementPolicy.LRU,
                    hit_class={k: k for k in range(num_classes)},
                    miss_class={k: k for k in range(num_classes)},
                    access_probs=None
                )
                config.cache_config[cache_node_idx] = cache_config

    def _parse_heterogeneous_servers(self, sn: Any, config: SimulatorConfig) -> None:
        """
        Parse heterogeneous server configuration from NetworkStruct.

        Extracts server types, compatibility matrices, and per-type service rates
        for heterogeneous multiserver queues.

        Args:
            sn: NetworkStruct (from model.getStruct())
            config: SimulatorConfig to populate
        """
        num_service = len(config.service_nodes)
        num_classes = config.num_classes

        # Check if heterogeneous server configuration exists
        if not hasattr(sn, 'nservertypes') or sn.nservertypes is None:
            return

        for svc_idx in range(num_service):
            if config.is_delay_node[svc_idx]:
                continue  # Skip delay nodes

            station_idx = config.service_stations[svc_idx]

            # Get number of server types for this station
            n_types = 0
            try:
                if hasattr(sn.nservertypes, 'get'):
                    n_types = int(sn.nservertypes.get(station_idx))
                elif hasattr(sn.nservertypes, '__getitem__'):
                    n_types = int(sn.nservertypes[station_idx])
                else:
                    n_types = int(self._get_value(sn.nservertypes, station_idx) or 0)
            except (TypeError, IndexError, ValueError):
                n_types = 0

            if n_types <= 0:
                continue

            config.num_server_types[svc_idx] = n_types

            # Get station object for map lookups
            station = None
            if hasattr(sn, 'stations') and sn.stations is not None:
                try:
                    station = sn.stations[station_idx]
                except (IndexError, TypeError):
                    pass

            # Get servers per type
            servers_per_type = [1] * n_types
            if hasattr(sn, 'serverspertype') and sn.serverspertype is not None and station is not None:
                try:
                    servers_matrix = sn.serverspertype.get(station)
                    if servers_matrix is not None:
                        for t in range(n_types):
                            try:
                                val = servers_matrix.get(t) if hasattr(servers_matrix, 'get') else servers_matrix[t]
                                servers_per_type[t] = int(val) if val else 1
                            except (IndexError, TypeError):
                                pass
                except (AttributeError, TypeError):
                    pass
            config.servers_per_type[svc_idx] = servers_per_type

            # Update total server count
            total_servers = sum(servers_per_type)
            config.num_servers[svc_idx] = total_servers

            # Build server ID to type ID mapping
            server_to_type = []
            for type_id in range(n_types):
                for _ in range(servers_per_type[type_id]):
                    server_to_type.append(type_id)
            config.server_to_type[svc_idx] = server_to_type

            # Get server-class compatibility matrix
            compat = [[True] * num_classes for _ in range(n_types)]
            if hasattr(sn, 'servercompat') and sn.servercompat is not None and station is not None:
                try:
                    compat_matrix = sn.servercompat.get(station)
                    if compat_matrix is not None:
                        for t in range(n_types):
                            for k in range(num_classes):
                                try:
                                    if hasattr(compat_matrix, 'get'):
                                        val = compat_matrix.get(t, k)
                                    else:
                                        val = compat_matrix[t, k] if hasattr(compat_matrix, '__getitem__') else True
                                    compat[t][k] = bool(val) if val is not None else True
                                except (IndexError, TypeError):
                                    pass
                except (AttributeError, TypeError):
                    pass
            config.server_compat[svc_idx] = compat

            # Get heterogeneous service rates
            hetero_mus = [[0.0] * num_classes for _ in range(n_types)]
            if hasattr(sn, 'heterorates') and sn.heterorates is not None and station is not None:
                try:
                    station_rates = sn.heterorates.get(station)
                    if station_rates is not None:
                        for type_id in range(n_types):
                            type_rates = station_rates.get(type_id)
                            if type_rates is not None:
                                for class_id in range(num_classes):
                                    try:
                                        rate = type_rates.get(class_id)
                                        if rate is not None and rate > 0 and np.isfinite(rate):
                                            hetero_mus[type_id][class_id] = float(rate)
                                    except (TypeError, KeyError):
                                        pass
                except (AttributeError, TypeError):
                    pass
            config.hetero_mus[svc_idx] = hetero_mus

            # Get heterogeneous scheduling policy
            policy = HeteroSchedPolicy.ORDER
            if hasattr(sn, 'heteroschedpolicy') and sn.heteroschedpolicy is not None and station is not None:
                try:
                    station_policy = sn.heteroschedpolicy.get(station)
                    if station_policy is not None:
                        policy = HeteroSchedPolicy.from_string(str(station_policy))
                except (AttributeError, TypeError):
                    pass
            config.hetero_sched_policy[svc_idx] = policy

    def _initialize_heterogeneous_servers(self) -> None:
        """
        Initialize heterogeneous server runtime state.

        Sets up busy counters, type ordering, and ALFS priority lists
        based on the parsed configuration.
        """
        for svc_idx in range(len(self.config.service_nodes)):
            n_types = self.config.num_server_types.get(svc_idx, 0)
            if n_types <= 0:
                continue

            # Initialize busy count per type
            self.busy_count_per_type[svc_idx] = [0] * n_types

            # Initialize server type order for ALIS/FAIRNESS (round-robin)
            self.server_type_order[svc_idx] = list(range(n_types))

            # Initialize ALFS order (sorted by number of compatible classes, ascending)
            # Least flexible servers (fewer compatible classes) are preferred
            compat = self.config.server_compat.get(svc_idx, [])
            if compat:
                flexibility = []
                for type_id in range(n_types):
                    num_compatible = sum(1 for c in compat[type_id] if c)
                    flexibility.append((type_id, num_compatible))
                # Sort by flexibility (ascending)
                flexibility.sort(key=lambda x: x[1])
                self.alfs_order[svc_idx] = [t[0] for t in flexibility]
            else:
                self.alfs_order[svc_idx] = list(range(n_types))

    def _initialize_schedulers(self) -> None:
        """
        Phase 3: Initialize schedulers for each service node.

        Creates appropriate scheduler instances based on scheduling strategies
        from the network configuration. Supports all 32+ scheduling disciplines.
        """
        from .scheduling.factory import create_scheduler

        num_service = len(self.config.service_nodes)

        for svc_idx in range(num_service):
            station_idx = self.config.service_stations[svc_idx]

            # Get scheduling strategy for this station (default: FCFS)
            sched_strategy = None
            if hasattr(self.sn, 'sched') and self.sn.sched is not None:
                sched_value = self._get_value(self.sn.sched, station_idx)
                sched_strategy = sched_value

            if sched_strategy is None:
                sched_strategy = 'FCFS'  # Default scheduling discipline

            # Extract service rates per class (for SEPT/LEPT schedulers)
            service_rates = []
            for k in range(self.config.num_classes):
                rate = self.config.mus[svc_idx, k] if svc_idx < self.config.mus.shape[0] else 0.0
                service_rates.append(rate)

            # Extract class weights from schedparam (for DPS/GPS schedulers)
            class_weights = None
            if hasattr(self.sn, 'schedparam') and self.sn.schedparam is not None:
                if station_idx < self.sn.schedparam.shape[0]:
                    class_weights = []
                    for k in range(self.config.num_classes):
                        if k < self.sn.schedparam.shape[1]:
                            class_weights.append(self.sn.schedparam[station_idx, k])
                        else:
                            class_weights.append(1.0)

            # Get number of servers for this service node
            num_servers = self.config.num_servers.get(svc_idx, 1)

            # Create scheduler instance
            try:
                scheduler = create_scheduler(
                    strategy=sched_strategy,
                    num_classes=self.config.num_classes,
                    num_servers=num_servers,
                    service_rates=service_rates,
                    class_weights=class_weights,
                    class_priorities=self.config.class_priorities
                )
                self.schedulers[svc_idx] = scheduler
            except Exception as e:
                # Fallback to FCFS if scheduler creation fails
                import sys
                print(f"Warning: Failed to create scheduler {sched_strategy} for station {station_idx}: {e}", file=sys.stderr)
                scheduler = create_scheduler(
                    strategy='FCFS',
                    num_classes=self.config.num_classes,
                    num_servers=num_servers,
                    service_rates=service_rates,
                    class_priorities=self.config.class_priorities
                )
                self.schedulers[svc_idx] = scheduler

    def _parse_fork_join_params(self) -> None:
        """
        Phase 5: Parse Fork/Join parameters and initialize ForkJoinManager.

        Extracts fork fan-out and join quorum parameters from the network structure
        and creates ForkNode and JoinNode instances in the manager.
        """
        # Create Fork nodes in the manager
        for fork_idx, fork_node_idx in enumerate(self.config.fork_nodes):
            # Extract fan-out parameter (default: 1)
            fan_out = 1
            if hasattr(self.sn, 'nodeparam') and self.sn.nodeparam is not None:
                try:
                    param = self.sn.nodeparam.get(self.sn.nodes[fork_node_idx]) if hasattr(self.sn, 'nodes') else None
                    if param and hasattr(param, 'fanOut'):
                        fan_out_val = param.fanOut
                        if not np.isnan(fan_out_val):
                            fan_out = int(fan_out_val)
                except Exception:
                    pass

            # Find all output destinations from routing matrix
            output_nodes = []
            if self.config.routing_matrix is not None and fork_node_idx < self.config.routing_matrix.shape[0]:
                for dest_node_idx in range(self.config.routing_matrix.shape[1]):
                    prob = self.config.routing_matrix[fork_node_idx, dest_node_idx]
                    if prob > 0:
                        output_nodes.append(dest_node_idx)

            # Find corresponding join node (scan all join nodes for matching reference)
            join_node_idx = None
            for j_idx, j_node in enumerate(self.config.join_nodes):
                # Simple heuristic: first join node after this fork
                if j_node > fork_node_idx:
                    join_node_idx = j_node
                    break

            # Register fork in manager
            self.fj_manager.add_fork(
                node_id=fork_node_idx,
                fan_out=fan_out,
                output_nodes=output_nodes,
                join_node_id=join_node_idx
            )

        # Create Join nodes in the manager
        for join_idx, join_node_idx in enumerate(self.config.join_nodes):
            # Extract quorum parameter (default: 0 = all tasks required)
            quorum = 0
            if hasattr(self.sn, 'nodeparam') and self.sn.nodeparam is not None:
                try:
                    param = self.sn.nodeparam.get(self.sn.nodes[join_node_idx]) if hasattr(self.sn, 'nodes') else None
                    if param and hasattr(param, 'quorum'):
                        quorum_val = param.quorum
                        if not np.isnan(quorum_val):
                            quorum = int(quorum_val)
                except Exception:
                    pass

            # Find corresponding fork node (scan all fork nodes)
            fork_node_idx = None
            for f_idx, f_node in enumerate(self.config.fork_nodes):
                if f_node < join_node_idx:
                    fork_node_idx = f_node

            # Register join in manager
            self.fj_manager.add_join(
                node_id=join_node_idx,
                fork_node_id=fork_node_idx,
                quorum=quorum
            )

    def _create_service_generator(self, svc_idx: int) -> Callable[[int], float]:
        """
        Phase 3: Create a service time generator for a service node.

        Returns a function that takes a class_id and returns a service time sample.
        Used by schedulers (especially SJF/LJF) to pre-sample service times.

        Args:
            svc_idx: Service node index

        Returns:
            Function (class_id: int) -> float that returns service time sample
        """
        def service_gen(class_id: int) -> float:
            """Generate service time for given class at this station."""
            if svc_idx >= len(self.config.mus):
                return 0.0

            service_rate = self.config.mus[svc_idx, class_id]
            if service_rate > 0:
                return self.rng.exponential(1.0 / service_rate)
            return 0.0

        return service_gen

    def _parse_spn_parameters(self, sn: Any, config: SimulatorConfig) -> None:
        """
        Extract SPN transition modes and parameters from NetworkStruct.

        Parses transition parameters including modes, timing strategies,
        priorities, weights, server counts, and arc conditions.
        """
        num_places = len(config.place_nodes)
        num_classes = config.num_classes
        num_nodes = config.num_nodes

        # Initialize transition modes dictionary
        self.transition_modes: Dict[int, List[TransitionModeInfo]] = {}

        # Initialize place-related mappings
        self.place_node_to_idx: Dict[int, int] = {
            node_idx: place_idx for place_idx, node_idx in enumerate(config.place_nodes)
        }
        self.transition_node_to_idx: Dict[int, int] = {
            node_idx: trans_idx for trans_idx, node_idx in enumerate(config.transition_nodes)
        }

        # Parse each transition
        for trans_idx, node_idx in enumerate(config.transition_nodes):
            modes = []

            # Try to get node parameters
            param = None
            if hasattr(sn, 'nodeparam') and sn.nodeparam is not None:
                try:
                    param = sn.nodeparam.get(node_idx) if hasattr(sn.nodeparam, 'get') else sn.nodeparam[node_idx]
                except (KeyError, IndexError, TypeError):
                    param = None

            if param is None:
                # Create a default single mode with Exp(1)
                mode_info = TransitionModeInfo(
                    mode_idx=0,
                    mode_name="Mode0",
                    timing_strategy='TIMED',
                    priority=0,
                    weight=1.0,
                    num_servers=1,
                    enabling_conditions=np.zeros((num_places, num_classes)),
                    inhibiting_conditions=np.full((num_places, num_classes), float('inf')),
                    firing_outcomes=np.zeros((num_nodes, num_classes)),
                    firing_distribution=None
                )
                modes.append(mode_info)
            else:
                # Get number of modes
                num_modes = getattr(param, 'nmodes', 1) if hasattr(param, 'nmodes') else 1
                if num_modes is None or num_modes < 1:
                    num_modes = 1

                for mode_idx in range(num_modes):
                    # Extract mode name
                    mode_name = f"Mode{mode_idx}"
                    if hasattr(param, 'modenames') and param.modenames is not None:
                        try:
                            mode_name = str(param.modenames[mode_idx])
                        except (IndexError, TypeError):
                            pass

                    # Extract timing strategy
                    timing_strategy = 'TIMED'
                    if hasattr(param, 'timingstrategies') and param.timingstrategies is not None:
                        try:
                            ts = param.timingstrategies[mode_idx]
                            ts_name = str(ts.name() if callable(getattr(ts, 'name', None)) else ts).upper()
                            if 'IMMEDIATE' in ts_name:
                                timing_strategy = 'IMMEDIATE'
                        except (IndexError, TypeError):
                            pass

                    # Extract priority
                    priority = 0
                    if hasattr(param, 'firingprio') and param.firingprio is not None:
                        try:
                            priority = int(param.firingprio[mode_idx])
                        except (IndexError, TypeError, ValueError):
                            pass

                    # Extract weight
                    weight = 1.0
                    if hasattr(param, 'fireweight') and param.fireweight is not None:
                        try:
                            weight = float(param.fireweight[mode_idx])
                        except (IndexError, TypeError, ValueError):
                            pass

                    # Extract number of servers
                    num_servers = 1
                    if hasattr(param, 'nmodeservers') and param.nmodeservers is not None:
                        try:
                            nms = param.nmodeservers
                            if hasattr(nms, 'getNumRows') and hasattr(nms, 'get'):
                                # Extract value at position mode_idx
                                cols = int(nms.getNumCols())
                                if mode_idx < cols:
                                    ns = float(nms.get(0, mode_idx))
                                    num_servers = int(ns) if not np.isinf(ns) else float('inf')
                            else:
                                ns = nms[mode_idx]
                                num_servers = int(ns) if ns is not None and not np.isinf(ns) else float('inf')
                        except (IndexError, TypeError, ValueError):
                            pass

                    # Extract enabling conditions (required tokens from places)
                    enabling = np.zeros((num_places, num_classes))
                    if hasattr(param, 'enabling') and param.enabling is not None:
                        try:
                            enabling = self._extract_condition_matrix(
                                param.enabling, mode_idx, num_places, num_classes, config.place_nodes
                            )
                        except Exception:
                            pass

                    # Extract inhibiting conditions (max tokens threshold)
                    # In Petri net semantics: "inhibit if tokens >= threshold"
                    # Convention: 0, NaN, Infinity mean "no inhibition" (use infinity)
                    inhibiting = np.full((num_places, num_classes), float('inf'))
                    if hasattr(param, 'inhibiting') and param.inhibiting is not None:
                        try:
                            raw_inhibiting = self._extract_condition_matrix(
                                param.inhibiting, mode_idx, num_places, num_classes, config.place_nodes
                            )
                            # Treat 0, NaN, Inf, and negative values as "no inhibition"
                            for p in range(num_places):
                                for c in range(num_classes):
                                    val = raw_inhibiting[p, c]
                                    if np.isnan(val) or np.isinf(val) or val <= 0:
                                        inhibiting[p, c] = float('inf')
                                    else:
                                        inhibiting[p, c] = val
                        except Exception:
                            pass

                    # Extract firing outcomes (tokens produced)
                    firing = np.zeros((num_nodes, num_classes))
                    if hasattr(param, 'firing') and param.firing is not None:
                        try:
                            firing = self._extract_outcome_matrix(
                                param.firing, mode_idx, num_nodes, num_classes
                            )
                        except Exception:
                            pass

                    # Extract distribution parameters
                    firing_dist = None
                    firing_proc_id = 'EXP'
                    firing_mean = 1.0
                    firing_phases = 1
                    firing_pie = None
                    firing_d0 = None
                    firing_d1 = None

                    # Get process ID (distribution type)
                    if hasattr(param, 'firingprocid') and param.firingprocid is not None:
                        try:
                            procids = list(param.firingprocid.values())
                            if mode_idx < len(procids):
                                proc_id_obj = procids[mode_idx]
                                if proc_id_obj is not None:
                                    firing_proc_id = str(proc_id_obj.name()) if hasattr(proc_id_obj, 'name') else str(proc_id_obj)
                        except Exception:
                            pass

                    # Get number of phases
                    if hasattr(param, 'firingphases') and param.firingphases is not None:
                        try:
                            fp = param.firingphases
                            if hasattr(fp, 'getNumRows') and hasattr(fp, 'get'):
                                # Extract value at position mode_idx
                                cols = int(fp.getNumCols())
                                if mode_idx < cols:
                                    firing_phases = int(fp.get(0, mode_idx))
                            else:
                                phases_arr = np.asarray(fp).flatten()
                                if mode_idx < len(phases_arr):
                                    firing_phases = int(phases_arr[mode_idx])
                        except Exception:
                            pass

                    # Get initial phase probabilities
                    if hasattr(param, 'firingpie') and param.firingpie is not None:
                        try:
                            pies = list(param.firingpie.values())
                            if mode_idx < len(pies):
                                pie_obj = pies[mode_idx]
                                if pie_obj is not None:
                                    if hasattr(pie_obj, 'getNumRows') and hasattr(pie_obj, 'get'):
                                        rows = int(pie_obj.getNumRows())
                                        cols = int(pie_obj.getNumCols())
                                        # Pie is usually a row vector (1 x num_phases)
                                        if rows == 1:
                                            firing_pie = np.array([float(pie_obj.get(0, c)) for c in range(cols)])
                                        else:
                                            firing_pie = np.array([float(pie_obj.get(r, 0)) for r in range(rows)])
                                    else:
                                        firing_pie = np.asarray(pie_obj).flatten()
                        except Exception:
                            pass

                    # Get D0 and D1 matrices from firingproc
                    if hasattr(param, 'firingproc') and param.firingproc is not None:
                        try:
                            procs = list(param.firingproc.values())
                            if mode_idx < len(procs):
                                firing_dist = procs[mode_idx]
                                if firing_dist is not None:
                                    # firingproc is a MatrixCell containing [D0, D1]
                                    d0_obj = firing_dist.get(0) if hasattr(firing_dist, 'get') else None
                                    d1_obj = firing_dist.get(1) if hasattr(firing_dist, 'get') else None
                                    if d0_obj is not None:
                                        firing_d0 = self._java_matrix_to_numpy(d0_obj)
                                    if d1_obj is not None:
                                        firing_d1 = self._java_matrix_to_numpy(d1_obj)
                        except Exception:
                            pass

                    # Calculate mean from D0 matrix and pie
                    if firing_d0 is not None:
                        try:
                            diag = np.diag(firing_d0)
                            rates = -diag  # Exit rates (positive)
                            if len(rates) > 0 and np.all(rates > 0):
                                if firing_proc_id.upper() == 'HYPEREXP' or firing_proc_id.upper() == 'HEXP':
                                    # HyperExp: weighted mean = sum(pie[i] / rate[i])
                                    if firing_pie is not None and len(firing_pie) == len(rates):
                                        firing_mean = np.sum(firing_pie / rates)
                                    else:
                                        # Fallback: use uniform weights
                                        firing_mean = np.mean(1.0 / rates)
                                else:
                                    # Erlang/PH: sequential phases, mean = sum(1/rate[i])
                                    firing_mean = np.sum(1.0 / rates)
                        except Exception:
                            pass

                    # Fallback: if distribution info not found in nodeparam,
                    # try to get it from the model's Transition nodes directly
                    if firing_proc_id == 'EXP' and firing_phases == 1 and self.model is not None:
                        try:
                            transition_node = self._get_transition_node(node_idx)
                            if transition_node is not None:
                                dist_info = self._extract_transition_distribution(
                                    transition_node, mode_idx
                                )
                                if dist_info:
                                    firing_proc_id = dist_info.get('proc_id', firing_proc_id)
                                    firing_mean = dist_info.get('mean', firing_mean)
                                    firing_phases = dist_info.get('phases', firing_phases)
                        except Exception:
                            pass

                    mode_info = TransitionModeInfo(
                        mode_idx=mode_idx,
                        mode_name=mode_name,
                        timing_strategy=timing_strategy,
                        priority=priority,
                        weight=weight,
                        num_servers=num_servers,
                        enabling_conditions=enabling,
                        inhibiting_conditions=inhibiting,
                        firing_outcomes=firing,
                        firing_distribution=firing_dist,
                        firing_proc_id=firing_proc_id,
                        firing_mean=firing_mean,
                        firing_phases=firing_phases,
                        firing_pie=firing_pie,
                        firing_d0=firing_d0,
                        firing_d1=firing_d1
                    )
                    modes.append(mode_info)

            self.transition_modes[trans_idx] = modes

    def _get_transition_node(self, node_idx: int) -> Optional[Any]:
        """Get Transition node from the model by node index."""
        if self.model is None:
            return None
        try:
            nodes = self.model.getNodes()
            if node_idx < len(nodes):
                node = nodes[node_idx]
                if 'Transition' in type(node).__name__:
                    return node
        except Exception:
            pass
        return None

    def _extract_transition_distribution(
        self, transition_node: Any, mode_idx: int
    ) -> Optional[Dict[str, Any]]:
        """
        Extract distribution info from a Transition node for a specific mode.

        Returns dict with proc_id, mean, phases, or None if extraction fails.
        """
        try:
            modes = transition_node.getModes()
            if mode_idx >= len(modes):
                return None

            mode = modes[mode_idx]
            dist = transition_node.getDistribution(mode)
            if dist is None:
                return None

            result = {}

            # Get distribution name/type
            if hasattr(dist, 'getName'):
                name = dist.getName().upper()
                if 'ERLANG' in name:
                    result['proc_id'] = 'ERLANG'
                elif 'HYPEREXP' in name or 'HEXP' in name:
                    result['proc_id'] = 'HYPEREXP'
                elif 'EXP' in name:
                    result['proc_id'] = 'EXP'
                elif 'PH' in name or 'APH' in name:
                    result['proc_id'] = 'PH'
                elif 'COXIAN' in name:
                    result['proc_id'] = 'COXIAN'
                else:
                    result['proc_id'] = name

            # Get mean
            if hasattr(dist, 'getMean'):
                result['mean'] = float(dist.getMean())

            # Get number of phases
            if hasattr(dist, 'getNumberOfPhases'):
                result['phases'] = int(dist.getNumberOfPhases())
            elif hasattr(dist, 'getOrder'):
                result['phases'] = int(dist.getOrder())

            return result if result else None
        except Exception:
            return None

    def _extract_condition_matrix(
        self, cond: Any, mode_idx: int, num_places: int, num_classes: int, place_nodes: List[int]
    ) -> np.ndarray:
        """Extract enabling/inhibiting condition matrix for a mode."""
        result = np.zeros((num_places, num_classes))
        try:
            # Handle various formats of condition storage
            if hasattr(cond, 'get'):
                mode_cond = cond.get(mode_idx)
            elif hasattr(cond, '__getitem__'):
                mode_cond = cond[mode_idx]
            else:
                return result

            if mode_cond is not None:
                # Handle Java Matrix objects
                if hasattr(mode_cond, 'getNumRows') and hasattr(mode_cond, 'get'):
                    num_rows = int(mode_cond.getNumRows())
                    num_cols = int(mode_cond.getNumCols())
                    for place_idx, node_idx in enumerate(place_nodes):
                        if node_idx < num_rows:
                            for class_idx in range(min(num_classes, num_cols)):
                                result[place_idx, class_idx] = float(mode_cond.get(node_idx, class_idx))
                else:
                    arr = np.asarray(mode_cond)
                    for place_idx, node_idx in enumerate(place_nodes):
                        if node_idx < arr.shape[0]:
                            for class_idx in range(min(num_classes, arr.shape[1] if arr.ndim > 1 else 1)):
                                if arr.ndim > 1:
                                    result[place_idx, class_idx] = arr[node_idx, class_idx]
                                else:
                                    result[place_idx, class_idx] = arr[node_idx]
        except Exception:
            pass
        return result

    def _extract_outcome_matrix(
        self, outcome: Any, mode_idx: int, num_nodes: int, num_classes: int
    ) -> np.ndarray:
        """Extract firing outcome matrix for a mode."""
        result = np.zeros((num_nodes, num_classes))
        try:
            if hasattr(outcome, 'get'):
                mode_outcome = outcome.get(mode_idx)
            elif hasattr(outcome, '__getitem__'):
                mode_outcome = outcome[mode_idx]
            else:
                return result

            if mode_outcome is not None:
                # Handle Java Matrix objects
                if hasattr(mode_outcome, 'getNumRows') and hasattr(mode_outcome, 'get'):
                    num_rows = int(mode_outcome.getNumRows())
                    num_cols = int(mode_outcome.getNumCols())
                    rows = min(num_nodes, num_rows)
                    cols = min(num_classes, num_cols)
                    for r in range(rows):
                        for c in range(cols):
                            result[r, c] = float(mode_outcome.get(r, c))
                else:
                    arr = np.asarray(mode_outcome)
                    rows = min(num_nodes, arr.shape[0])
                    cols = min(num_classes, arr.shape[1] if arr.ndim > 1 else 1)
                    if arr.ndim > 1:
                        result[:rows, :cols] = arr[:rows, :cols]
                    else:
                        result[:rows, 0] = arr[:rows]
        except Exception:
            pass
        return result

    def _get_node_type(self, nodetype: Any, idx: int) -> NodeType:
        """Get node type for a given index."""
        try:
            nt = nodetype[idx]

            # Handle numeric node types (from mock or native structs)
            if isinstance(nt, (int, np.integer)):
                # Map numeric values to NodeType enum
                # Java NodeType ordinals: SOURCE=0, DELAY=1, QUEUE=2, CLASSSWITCH=3, FORK=4,
                # JOIN=5, SINK=6, ROUTER=7, CACHE=8, LOGGER=9, PLACE=10, TRANSITION=11
                # Python NodeType values defined above:
                # SOURCE=0, SINK=1, QUEUE=2, DELAY=3, FORK=4, JOIN=5, ROUTER=6, etc.
                node_type_map = {
                    0: NodeType.SOURCE,
                    1: NodeType.SINK,
                    2: NodeType.QUEUE,
                    3: NodeType.DELAY,
                    4: NodeType.FORK,
                    5: NodeType.JOIN,
                    6: NodeType.ROUTER,
                    7: NodeType.CLASSSWITCH,
                    8: NodeType.CACHE,
                    9: NodeType.PLACE,
                    10: NodeType.TRANSITION,
                    11: NodeType.LOGGER,
                }
                return node_type_map.get(int(nt), NodeType.QUEUE)

            # Handle both Python enum (name is property) and JPype enum (name() is method)
            if hasattr(nt, 'name'):
                name_attr = nt.name
                # If it's a method (JPype), call it; if it's a property (Python), use as-is
                if callable(name_attr):
                    name = str(name_attr()).upper()
                else:
                    name = str(name_attr).upper()
            else:
                name = str(nt).upper()

            if 'SOURCE' in name:
                return NodeType.SOURCE
            elif 'SINK' in name:
                return NodeType.SINK
            elif 'QUEUE' in name:
                return NodeType.QUEUE
            elif 'DELAY' in name:
                return NodeType.DELAY
            elif 'FORK' in name:
                return NodeType.FORK
            elif 'JOIN' in name:
                return NodeType.JOIN
            elif 'ROUTER' in name:
                return NodeType.ROUTER
            elif 'CLASSSWITCH' in name or 'CLASS_SWITCH' in name:
                return NodeType.CLASSSWITCH
            elif 'PLACE' in name:
                return NodeType.PLACE
            elif 'TRANSITION' in name:
                return NodeType.TRANSITION
            else:
                return NodeType.QUEUE
        except Exception:
            return NodeType.QUEUE

    def _java_matrix_to_numpy(self, mat: Any) -> Optional[np.ndarray]:
        """Convert a Java Matrix object to numpy array."""
        try:
            if mat is None:
                return None
            if hasattr(mat, 'getNumRows') and hasattr(mat, 'get'):
                rows = int(mat.getNumRows())
                cols = int(mat.getNumCols())
                result = np.zeros((rows, cols))
                for r in range(rows):
                    for c in range(cols):
                        result[r, c] = float(mat.get(r, c))
                return result
            else:
                return np.asarray(mat)
        except Exception:
            return None

    def _get_value(self, arr: Any, idx: int) -> Optional[float]:
        """Get value from array (handles various array types)."""
        try:
            if hasattr(arr, 'get'):
                val = arr.get(idx)
            elif hasattr(arr, '__getitem__'):
                val = arr[idx]
            else:
                return None
            # Handle numpy arrays - extract scalar to avoid deprecation warning
            if hasattr(val, 'item'):
                return float(val.item())
            return float(val)
        except Exception:
            return None

    def _get_int_value(self, arr: Any, idx: int) -> Optional[int]:
        """Get integer value from array (handles various array types including JPype and 2D arrays)."""
        try:
            if hasattr(arr, 'get'):
                val = arr.get(idx)
            elif hasattr(arr, 'ndim') and arr.ndim == 2 and arr.shape[0] == 1:
                # Handle (1, N) arrays - flatten to access element at column idx
                val = arr[0, idx]
            elif hasattr(arr, '__getitem__'):
                val = arr[idx]
            else:
                return None
            # Handle numpy arrays - extract scalar to avoid deprecation warning
            if hasattr(val, 'item'):
                return int(val.item())
            return int(val)
        except Exception:
            return None

    def _get_matrix_value(self, mat: Any, row: int, col: int) -> Optional[float]:
        """Get value from matrix (handles various matrix types)."""
        try:
            if hasattr(mat, 'get'):
                val = mat.get(row, col)
            elif hasattr(mat, '__getitem__'):
                val = mat[row, col]
            else:
                return None
            # Handle numpy arrays - extract scalar to avoid deprecation warning
            if hasattr(val, 'item'):
                return float(val.item())
            return float(val)
        except Exception:
            return None

    def _get_next_job_id(self) -> int:
        """Get next unique job ID."""
        job_id = self.next_job_id
        self.next_job_id += 1
        return job_id

    def simulate(self, max_events: int) -> None:
        """
        Run simulation (steady-state or transient).

        In steady-state mode:
            - Runs until max_events service completions
            - Uses MSER-5 warmup detection

        In transient mode (when timespan or transient=True):
            - Runs until end_time or max_events (whichever first)
            - No warmup period (starts collecting from time 0)
            - Collects time-series data at regular intervals

        Args:
            max_events: Maximum number of service completions
        """
        # Set up MSER observation interval (~1000 observations like Java)
        # Only enabled for steady-state mode
        if self.mser_enabled and not self.is_transient_mode:
            target_observations = 1000
            self.mser_observation_interval = max(1, max_events // target_observations)
            self.last_mser_sample_time = 0.0
            self.last_mser_event_count = 0
            # Initialize last queue time snapshots
            for svc_idx in range(len(self.config.service_nodes)):
                self.last_mser_queue_time[svc_idx] = np.zeros(self.config.num_classes)

        # Start arrival processes for open classes
        for src_idx, src_node in enumerate(self.config.source_nodes):
            for k in range(self.config.num_classes):
                if self.config.is_open_class[k]:
                    arrival_rate = self.config.lambdas[src_idx, k]
                    if arrival_rate > 0:
                        self.env.process(self._arrival_process(src_idx, k))

        # Initialize closed class populations
        self._init_closed_class_populations()

        # Initialize SPN state if places/transitions exist
        if self.config.place_nodes or self.config.transition_nodes:
            self._init_spn_state()

        # Run simulation
        if self.is_transient_mode:
            # Transient mode: run until end_time or max_events
            self._run_transient_simulation(max_events)
        else:
            # Steady-state mode: run until max_events
            while self.total_event_count < max_events:
                try:
                    self.env.step()
                except simpy.core.EmptySchedule:
                    break

            # Apply MSER-5 warmup detection if enabled
            if self.mser_enabled:
                self._apply_mser5_truncation()

    def _run_transient_simulation(self, max_events: int) -> None:
        """
        Run simulation in transient mode.

        Collects time-series data without warmup period.
        Stops when end_time is reached or max_events completed.

        Args:
            max_events: Maximum number of service completions
        """
        # Sample initial state at time 0
        self._collect_transient_sample()

        # Run simulation
        while self.total_event_count < max_events:
            try:
                self.env.step()

                # Check time-based termination
                if self.env.now >= self.transient_end_time:
                    break

                # Collect transient sample at regular intervals
                if self.env.now - self.last_transient_sample_time >= self.transient_interval:
                    if len(self.transient_times) < self.transient_max_samples:
                        self._collect_transient_sample()

            except simpy.core.EmptySchedule:
                break

        # Collect final sample
        if len(self.transient_times) < self.transient_max_samples:
            self._collect_transient_sample()

    def _collect_transient_sample(self) -> None:
        """
        Collect a transient sample of current system state.

        Records queue lengths, utilizations, and throughputs at the current time
        for all service nodes and classes.
        """
        current_time = self.env.now
        self.transient_times.append(current_time)
        self.last_transient_sample_time = current_time

        num_service = len(self.config.service_nodes)

        for svc_idx in range(num_service):
            stats = self.stats[svc_idx]
            num_servers = self.config.num_servers.get(svc_idx, 1)
            is_delay = self.config.is_delay_node[svc_idx] if svc_idx < len(self.config.is_delay_node) else False

            for k in range(self.config.num_classes):
                # Current queue length
                queue_length = float(stats.current_queue_length[k])

                # Compute instantaneous utilization
                # (busy time since last sample / elapsed time)
                if current_time > 0:
                    # Update time-weighted stats first
                    stats.update_queue(k, current_time, delta=0)

                    # Compute utilization
                    if is_delay:
                        util = 0.0  # Delay nodes have infinite servers
                    else:
                        total_busy_time = stats.total_busy_time[k]
                        util = total_busy_time / (current_time * num_servers)
                        util = min(1.0, max(0.0, util))
                else:
                    util = 0.0

                # Compute instantaneous throughput
                if current_time > 0:
                    throughput = stats.completed_customers[k] / current_time
                else:
                    throughput = 0.0

                # Store values
                self.transient_queue_lengths[svc_idx][k].append(queue_length)
                self.transient_utilizations[svc_idx][k].append(util)
                self.transient_throughputs[svc_idx][k].append(throughput)

    def _init_closed_class_populations(self) -> None:
        """Initialize closed class job populations at reference stations."""
        for k in range(self.config.num_classes):
            if self.config.is_closed_class[k]:
                pop = self.config.closed_class_population[k]
                ref_station = self.config.reference_station[k]

                # Find service node index for reference station
                svc_idx = None
                for idx, station in enumerate(self.config.service_stations):
                    if station == ref_station:
                        svc_idx = idx
                        break

                if svc_idx is not None and pop > 0:
                    # Place initial jobs at reference station
                    for _ in range(pop):
                        customer = Customer(
                            class_id=k,
                            priority=self.config.class_priorities[k],
                            system_arrival_time=0.0,
                            queue_arrival_time=0.0,
                            job_id=self._get_next_job_id(),
                        )
                        self._arrive_at_service_node(svc_idx, customer)

    def _init_spn_state(self) -> None:
        """Initialize SPN place tokens and transition server counts."""
        num_classes = self.config.num_classes

        # Initialize place tokens
        for place_idx, node_idx in enumerate(self.config.place_nodes):
            # Get initial token count from state if available
            initial_tokens = np.zeros(num_classes)
            if hasattr(self.sn, 'state') and self.sn.state is not None:
                try:
                    # state is indexed by stateful index, not node index
                    # Convert node_idx to stateful_idx using nodeToStateful
                    stateful_idx = None
                    if hasattr(self.sn, 'nodeToStateful') and self.sn.nodeToStateful is not None:
                        n2sf = np.asarray(self.sn.nodeToStateful).flatten()
                        if node_idx < len(n2sf):
                            stateful_idx = int(n2sf[node_idx])

                    if stateful_idx is not None and stateful_idx >= 0:
                        # sn.state is an array of objects, each element is the state for that stateful node
                        state_arr = self.sn.state
                        if stateful_idx < len(state_arr):
                            node_state = state_arr[stateful_idx]
                            if node_state is not None:
                                node_state = np.asarray(node_state).flatten()
                                # node_state contains token counts per class
                                for k in range(min(num_classes, len(node_state))):
                                    initial_tokens[k] = node_state[k]
                except Exception:
                    pass

            self.place_tokens[place_idx] = initial_tokens.astype(int)
            self.place_token_time[place_idx] = np.zeros(num_classes)
            self.last_place_update[place_idx] = 0.0
            self.place_completions[place_idx] = np.zeros(num_classes, dtype=int)

        # Initialize transition server counts and in-service tracking
        num_classes = self.config.num_classes
        for trans_idx in range(len(self.config.transition_nodes)):
            modes = self.transition_modes.get(trans_idx, [])
            self.transition_in_service[trans_idx] = {m.mode_idx: 0 for m in modes}
            self.transition_in_service_time[trans_idx] = np.zeros(num_classes)
            self.transition_in_service_count[trans_idx] = np.zeros(num_classes, dtype=int)
            self.last_transition_update[trans_idx] = 0.0

        # Initialize phase-level tracking for SPN race semantics
        self.in_service_firings: List[InServiceFiring] = []
        self.next_firing_id = 0
        # Track phase completions for throughput (MATLAB SSA compatibility)
        self.spn_phase_completions: Dict[int, np.ndarray] = {}
        for place_idx in range(len(self.config.place_nodes)):
            self.spn_phase_completions[place_idx] = np.zeros(num_classes, dtype=int)

        # Fire initially enabled transitions using phase-level race
        self._start_phase_level_simulation()

    def _start_phase_level_simulation(self) -> None:
        """
        Start SPN simulation using winner-selection race semantics.

        This approach samples complete firing times for all enabled modes,
        picks the winner (minimum time), and fires only the winner.
        This matches the Kotlin DES implementation.
        """
        if not self.config.transition_nodes:
            return

        # Fire initially enabled transitions
        # This will sample firing times, pick winners, and schedule completions.
        # When completions occur, _transition_completion calls _check_and_fire_transitions
        # to trigger more firings (cascade).
        self._check_and_fire_transitions()

    def _phase_level_race_process(self):
        """
        SimPy process for phase-level race simulation.

        This process continuously:
        1. Starts new firings for enabled modes
        2. Samples phase completion times for all in-service firings
        3. Advances the minimum time, completing or advancing phases
        """
        num_places = len(self.config.place_nodes)
        num_classes = self.config.num_classes

        while True:
            # Start new firings for enabled modes
            self._start_enabled_firings()

            # If no firings in progress, wait for tokens (shouldn't happen in closed network)
            if not self.in_service_firings:
                # Check if simulation should end
                if self.total_event_count >= self.options.samples:
                    break
                # Wait a small time and retry
                yield self.env.timeout(0.001)
                continue

            # Sample phase completion times for all in-service firings
            current_time = self.env.now
            min_time = float('inf')
            winner_idx = -1

            for i, firing in enumerate(self.in_service_firings):
                # Calculate remaining time to next phase completion
                if firing.next_event_time <= current_time:
                    # Need to sample new phase time
                    phase_time = self._sample_phase_time(firing.mode, firing.current_phase)
                    firing.next_event_time = current_time + phase_time

                if firing.next_event_time < min_time:
                    min_time = firing.next_event_time
                    winner_idx = i

            if winner_idx < 0:
                break

            # Wait until the winning phase completes
            wait_time = min_time - current_time
            if wait_time > 0:
                yield self.env.timeout(wait_time)

            # Process the winning phase completion
            winner = self.in_service_firings[winner_idx]
            self._process_phase_completion(winner, winner_idx)

            # Increment event count
            self.total_event_count += 1

            # Check termination condition
            if self.total_event_count >= self.options.samples:
                break

    def _start_enabled_firings(self) -> None:
        """
        Start firing events using phase-level race semantics.

        Multi-server support: keeps firing until no more tokens/servers available.
        Phase-level racing: all modes in a race group race at the phase level.
        This matches MATLAB SSA's CTMC-based throughput calculation.
        """
        num_places = len(self.config.place_nodes)
        num_classes = self.config.num_classes
        current_time = self.env.now

        # Keep firing until no more firings are possible
        fired_any = True
        while fired_any:
            fired_any = False

            for trans_idx in range(len(self.config.transition_nodes)):
                modes = self.transition_modes.get(trans_idx, [])
                if not modes:
                    continue

                # Find enabled modes (respects server availability and token counts)
                enabled = [m for m in modes if self._is_mode_enabled_with_servers(trans_idx, m)]
                if not enabled:
                    continue

                # Priority filtering
                max_prio = max(m.priority for m in enabled)
                high_prio = [m for m in enabled if m.priority == max_prio]

                # Start a race group: consume tokens once, all modes race at phase level
                first_mode = high_prio[0]

                # Update place statistics before consuming tokens
                for place_idx in range(num_places):
                    elapsed = current_time - self.last_place_update.get(place_idx, 0.0)
                    if elapsed > 0:
                        tokens = self.place_tokens.get(place_idx, np.zeros(num_classes))
                        self.place_token_time[place_idx] = self.place_token_time.get(
                            place_idx, np.zeros(num_classes)
                        ) + tokens * elapsed
                        self.last_place_update[place_idx] = current_time

                # Consume input tokens (once for all racing modes)
                consumed_tokens: Dict[int, np.ndarray] = {}
                for place_idx in range(num_places):
                    consumed = np.zeros(num_classes)
                    for class_idx in range(num_classes):
                        required = int(first_mode.enabling_conditions[place_idx, class_idx])
                        if required > 0:
                            self.place_tokens[place_idx][class_idx] -= required
                            consumed[class_idx] = required
                            if place_idx not in self.place_completions:
                                self.place_completions[place_idx] = np.zeros(num_classes, dtype=int)
                            self.place_completions[place_idx][class_idx] += required
                            if trans_idx not in self.transition_in_service_count:
                                self.transition_in_service_count[trans_idx] = np.zeros(num_classes, dtype=int)
                            self.transition_in_service_count[trans_idx][class_idx] += required
                    if np.any(consumed > 0):
                        consumed_tokens[place_idx] = consumed

                # Use a common race_id for all modes in this race group
                race_id = self.next_firing_id

                # Start ALL modes racing in this race group
                for mode in high_prio:
                    # Increment in-service counter per mode
                    if trans_idx not in self.transition_in_service:
                        self.transition_in_service[trans_idx] = {}
                    if mode.mode_idx not in self.transition_in_service[trans_idx]:
                        self.transition_in_service[trans_idx][mode.mode_idx] = 0
                    self.transition_in_service[trans_idx][mode.mode_idx] += 1

                    # Sample initial phase time for this mode
                    initial_phase_time = self._sample_phase_time(mode, 0)
                    firing = InServiceFiring(
                        firing_id=race_id,  # Same race_id for all modes in this race
                        trans_idx=trans_idx,
                        mode=mode,
                        current_phase=0,
                        next_event_time=current_time + initial_phase_time,
                        consumed_tokens=consumed_tokens if mode == first_mode else {}
                    )
                    self.in_service_firings.append(firing)

                self.next_firing_id += 1
                fired_any = True

    def _is_mode_enabled_with_servers(self, trans_idx: int, mode) -> bool:
        """Check if a mode is enabled (has tokens AND has free servers)."""
        # Check server availability for this specific mode
        current_in_service = 0
        if trans_idx in self.transition_in_service:
            current_in_service = self.transition_in_service[trans_idx].get(mode.mode_idx, 0)
        if current_in_service >= mode.num_servers:
            return False

        # Check token availability
        return self._is_mode_enabled(trans_idx, mode)

    def _sample_phase_time(self, mode: TransitionModeInfo, phase: int) -> float:
        """
        Sample time for a single phase transition.

        For Erlang(k, λ), each phase has rate λ = k/mean.
        For Exp, there's only one phase with rate 1/mean.
        For HyperExp, uses the selected branch rate.
        """
        proc_id = mode.firing_proc_id.upper() if mode.firing_proc_id else 'EXP'

        if proc_id == 'IMMEDIATE':
            return 0.0

        if proc_id == 'EXP':
            return self.rng.exponential(mode.firing_mean)

        elif proc_id == 'ERLANG':
            # For Erlang(mean, k), each phase has rate k/mean
            k = mode.firing_phases if mode.firing_phases > 0 else 1
            rate = k / mode.firing_mean if mode.firing_mean > 0 else 1.0
            return self.rng.exponential(1.0 / rate)

        elif proc_id == 'HYPEREXP' or proc_id == 'HEXP':
            # HyperExp is memoryless - sample from mixture
            if mode.firing_d0 is not None:
                rates = -np.diag(mode.firing_d0)
                if len(rates) >= 2 and np.all(rates > 0):
                    # Get mixture probability
                    if mode.firing_pie is not None and len(mode.firing_pie) > 0:
                        p = mode.firing_pie[0]
                    elif mode.firing_d1 is not None and mode.firing_d1.size > 0:
                        lambda1 = rates[0]
                        p = mode.firing_d1[0, 0] / lambda1 if lambda1 > 0 else 0.5
                    else:
                        p = 0.5
                    p = np.clip(p, 0.0, 1.0)
                    # Choose phase and sample
                    branch = 0 if self.rng.random() < p else 1
                    return self.rng.exponential(1.0 / rates[branch])
            return self.rng.exponential(mode.firing_mean)

        elif proc_id == 'PH' or proc_id == 'APH' or proc_id == 'COXIAN':
            # For general PH, use D0 to get phase rate
            if mode.firing_d0 is not None and phase < mode.firing_d0.shape[0]:
                exit_rate = -mode.firing_d0[phase, phase]
                if exit_rate > 0:
                    return self.rng.exponential(1.0 / exit_rate)
            return self.rng.exponential(mode.firing_mean)

        else:
            return self.rng.exponential(mode.firing_mean)

    def _process_phase_completion(self, firing: InServiceFiring, firing_idx: int) -> None:
        """
        Process completion of a phase in a firing.

        Either advances to next phase or completes the firing.
        """
        mode = firing.mode
        proc_id = mode.firing_proc_id.upper() if mode.firing_proc_id else 'EXP'
        num_phases = mode.firing_phases if mode.firing_phases > 0 else 1
        current_time = self.env.now
        num_classes = self.config.num_classes

        # Track phase completion for throughput (MATLAB SSA compatibility)
        # Each phase completion contributes to throughput, not just firing completions
        for place_idx in range(len(self.config.place_nodes)):
            if place_idx not in self.spn_phase_completions:
                self.spn_phase_completions[place_idx] = np.zeros(num_classes, dtype=int)
            # Increment for the class that consumed the token
            for class_idx in range(num_classes):
                # Use the first mode's enabling conditions to find which class is involved
                if mode.enabling_conditions[place_idx, class_idx] > 0:
                    self.spn_phase_completions[place_idx][class_idx] += 1
                    break  # Only count once per phase completion

        # Check if this is the last phase (absorption)
        is_last_phase = False
        if proc_id == 'EXP' or proc_id == 'IMMEDIATE' or proc_id == 'HYPEREXP' or proc_id == 'HEXP':
            is_last_phase = True  # Single phase distributions
        elif proc_id == 'ERLANG':
            is_last_phase = (firing.current_phase >= num_phases - 1)
        elif proc_id == 'PH' or proc_id == 'APH' or proc_id == 'COXIAN':
            # Check D0 for absorption probability
            if mode.firing_d0 is not None:
                phase = firing.current_phase
                if phase < mode.firing_d0.shape[0]:
                    exit_rate = -mode.firing_d0[phase, phase]
                    # Sum of transition rates to other phases
                    trans_sum = 0.0
                    for j in range(mode.firing_d0.shape[0]):
                        if j != phase and mode.firing_d0[phase, j] > 0:
                            trans_sum += mode.firing_d0[phase, j]
                    # Probability of absorption
                    p_absorb = 1.0 - trans_sum / exit_rate if exit_rate > 0 else 1.0
                    is_last_phase = self.rng.random() < p_absorb
            else:
                is_last_phase = True
        else:
            is_last_phase = True

        if is_last_phase:
            # Complete the firing
            self._complete_firing(firing, firing_idx)
        else:
            # Advance to next phase
            firing.current_phase += 1
            firing.next_event_time = current_time + self._sample_phase_time(mode, firing.current_phase)

    def _complete_firing(self, firing: InServiceFiring, firing_idx: int) -> None:
        """
        Complete a firing: produce output tokens and cancel all other racing firings.

        In SPN race semantics, when one mode completes, all other racing modes
        are canceled and must restart the race from the beginning.
        """
        trans_idx = firing.trans_idx
        mode = firing.mode
        num_places = len(self.config.place_nodes)
        num_classes = self.config.num_classes
        current_time = self.env.now

        # Update place statistics before producing tokens
        for place_idx in range(num_places):
            elapsed = current_time - self.last_place_update.get(place_idx, 0.0)
            if elapsed > 0:
                tokens = self.place_tokens.get(place_idx, np.zeros(num_classes))
                self.place_token_time[place_idx] = self.place_token_time.get(
                    place_idx, np.zeros(num_classes)
                ) + tokens * elapsed
                self.last_place_update[place_idx] = current_time

        # Update in-service time statistics
        if trans_idx in self.transition_in_service_count:
            if hasattr(self, 'last_transition_update'):
                trans_elapsed = current_time - self.last_transition_update.get(trans_idx, 0.0)
                if trans_elapsed > 0:
                    self.transition_in_service_time[trans_idx] += (
                        self.transition_in_service_count[trans_idx] * trans_elapsed
                    )
                    self.last_transition_update[trans_idx] = current_time

        # Decrement in-service counter for the completed firing
        self.transition_in_service[trans_idx][mode.mode_idx] -= 1

        # Cancel other racing firings in the SAME race group (same race_id/firing_id)
        # In phase-level racing, when one mode completes, all others in the race are canceled
        race_id = firing.firing_id
        to_remove = []
        for i, other_firing in enumerate(self.in_service_firings):
            if other_firing.firing_id == race_id:
                to_remove.append(i)
                if other_firing is not firing:
                    # Decrement in-service counter for canceled loser
                    other_trans = other_firing.trans_idx
                    other_mode = other_firing.mode
                    if other_trans in self.transition_in_service:
                        if other_mode.mode_idx in self.transition_in_service[other_trans]:
                            self.transition_in_service[other_trans][other_mode.mode_idx] -= 1

        # Remove firings from this race group
        for i in reversed(to_remove):
            del self.in_service_firings[i]

        # Decrement in-service count based on consumed tokens
        for place_idx, consumed in firing.consumed_tokens.items():
            for class_idx in range(num_classes):
                if consumed[class_idx] > 0:
                    if trans_idx in self.transition_in_service_count:
                        self.transition_in_service_count[trans_idx][class_idx] -= int(consumed[class_idx])
                        if self.transition_in_service_count[trans_idx][class_idx] < 0:
                            self.transition_in_service_count[trans_idx][class_idx] = 0

        # Produce output tokens
        num_nodes = self.config.num_nodes
        for node_idx in range(num_nodes):
            for class_idx in range(num_classes):
                produced = int(mode.firing_outcomes[node_idx, class_idx])
                if produced > 0:
                    place_idx = self.place_node_to_idx.get(node_idx)
                    if place_idx is not None:
                        self.place_tokens[place_idx][class_idx] += produced

    # ==================== SPN Core Methods (Legacy) ====================

    def _check_and_fire_transitions(self) -> None:
        """
        Fire all enabled transitions using race semantics.

        Loops until no more transitions can fire, supporting multi-server
        concurrent firing as in the Kotlin implementation.
        """
        if not self.config.transition_nodes:
            return

        fired_any = True
        while fired_any:
            fired_any = False
            for trans_idx in range(len(self.config.transition_nodes)):
                modes = self.transition_modes.get(trans_idx, [])
                if not modes:
                    continue

                # Find enabled modes
                enabled = [m for m in modes if self._is_mode_enabled(trans_idx, m)]
                if not enabled:
                    continue

                # Priority filtering - select highest priority modes
                max_prio = max(m.priority for m in enabled)
                high_prio = [m for m in enabled if m.priority == max_prio]

                # Race semantics: pick winner
                winner, firing_time = self._select_mode_by_race(trans_idx, high_prio)
                if winner is not None:
                    self._fire_transition(trans_idx, winner, firing_time)
                    fired_any = True

    def _is_mode_enabled(self, trans_idx: int, mode: TransitionModeInfo) -> bool:
        """
        Check if a transition mode can fire.

        Checks server availability, enabling conditions (required tokens),
        and inhibiting conditions (max token threshold).
        """
        # Server availability check
        in_service = self.transition_in_service.get(trans_idx, {}).get(mode.mode_idx, 0)
        if in_service >= mode.num_servers:
            return False

        num_places = len(self.config.place_nodes)
        num_classes = self.config.num_classes

        # Enabling conditions: required tokens must be present
        for place_idx in range(num_places):
            for class_idx in range(num_classes):
                required = mode.enabling_conditions[place_idx, class_idx]
                available = self.place_tokens.get(place_idx, np.zeros(num_classes))[class_idx]
                if required > 0 and available < required:
                    return False

        # Inhibiting conditions: current tokens must be below threshold
        for place_idx in range(num_places):
            for class_idx in range(num_classes):
                threshold = mode.inhibiting_conditions[place_idx, class_idx]
                current = self.place_tokens.get(place_idx, np.zeros(num_classes))[class_idx]
                if threshold < float('inf') and current >= threshold:
                    return False

        return True

    def _select_mode_by_race(
        self, trans_idx: int, modes: List[TransitionModeInfo]
    ) -> Tuple[Optional[TransitionModeInfo], float]:
        """
        Select winning mode by sampling firing times (race semantics).

        For immediate transitions, uses negative weight for tie-breaking.
        For timed transitions, samples firing time from distribution.
        """
        if not modes:
            return None, 0.0

        if len(modes) == 1:
            mode = modes[0]
            if mode.timing_strategy == 'IMMEDIATE':
                return mode, 0.0
            else:
                return mode, self._sample_firing_time(trans_idx, mode)

        min_time = float('inf')
        winner = None

        for mode in modes:
            if mode.timing_strategy == 'IMMEDIATE':
                # Use negative weight for tie-breaking (larger weight = more negative = wins)
                t = -mode.weight
            else:
                t = self._sample_firing_time(trans_idx, mode)

            if t < min_time:
                min_time = t
                winner = mode

        # Convert negative times (from immediate) to 0
        actual_time = 0.0 if min_time < 0 else min_time
        return winner, actual_time

    def _fire_transition(
        self, trans_idx: int, mode: TransitionModeInfo, firing_time: float
    ) -> None:
        """
        Fire a transition: consume input tokens and schedule completion.
        """
        current_time = self.env.now
        num_places = len(self.config.place_nodes)
        num_classes = self.config.num_classes

        # Update place statistics before consuming tokens
        for place_idx in range(num_places):
            elapsed = current_time - self.last_place_update.get(place_idx, 0.0)
            if elapsed > 0:
                tokens = self.place_tokens.get(place_idx, np.zeros(num_classes))
                self.place_token_time[place_idx] = self.place_token_time.get(
                    place_idx, np.zeros(num_classes)
                ) + tokens * elapsed
                self.last_place_update[place_idx] = current_time

        # Update in-service time before changing counts
        if trans_idx in self.transition_in_service_count:
            elapsed = current_time - self.last_place_update.get(-1, 0.0)  # Use a dummy key for transition time
            if elapsed > 0 and hasattr(self, 'last_transition_update'):
                trans_elapsed = current_time - self.last_transition_update.get(trans_idx, 0.0)
                self.transition_in_service_time[trans_idx] += self.transition_in_service_count[trans_idx] * trans_elapsed

        # Consume input tokens (from enabling conditions)
        for place_idx in range(num_places):
            for class_idx in range(num_classes):
                required = int(mode.enabling_conditions[place_idx, class_idx])
                if required > 0:
                    self.place_tokens[place_idx][class_idx] -= required
                    # Track completions
                    if place_idx not in self.place_completions:
                        self.place_completions[place_idx] = np.zeros(num_classes, dtype=int)
                    self.place_completions[place_idx][class_idx] += required
                    # Track in-service per class
                    if trans_idx not in self.transition_in_service_count:
                        self.transition_in_service_count[trans_idx] = np.zeros(num_classes, dtype=int)
                    self.transition_in_service_count[trans_idx][class_idx] += required

        # Update last transition update time
        if not hasattr(self, 'last_transition_update'):
            self.last_transition_update: Dict[int, float] = {}
        self.last_transition_update[trans_idx] = current_time

        # Increment in-service counter
        if trans_idx not in self.transition_in_service:
            self.transition_in_service[trans_idx] = {}
        if mode.mode_idx not in self.transition_in_service[trans_idx]:
            self.transition_in_service[trans_idx][mode.mode_idx] = 0
        self.transition_in_service[trans_idx][mode.mode_idx] += 1

        # Schedule completion event
        self.env.process(self._transition_completion(trans_idx, mode, firing_time))

    def _transition_completion(
        self, trans_idx: int, mode: TransitionModeInfo, delay: float
    ):
        """
        SimPy process for transition firing completion.

        Waits for firing delay, then produces output tokens and triggers cascade.
        """
        # Wait for firing time
        yield self.env.timeout(delay)

        current_time = self.env.now
        num_places = len(self.config.place_nodes)
        num_classes = self.config.num_classes

        # Update in-service time statistics before decrementing
        if trans_idx in self.transition_in_service_count:
            if hasattr(self, 'last_transition_update'):
                trans_elapsed = current_time - self.last_transition_update.get(trans_idx, 0.0)
                if trans_elapsed > 0:
                    self.transition_in_service_time[trans_idx] += self.transition_in_service_count[trans_idx] * trans_elapsed
                    self.last_transition_update[trans_idx] = current_time

        # Decrement in-service counter
        self.transition_in_service[trans_idx][mode.mode_idx] -= 1

        # Update place statistics before producing tokens
        for place_idx in range(num_places):
            elapsed = current_time - self.last_place_update.get(place_idx, 0.0)
            if elapsed > 0:
                tokens = self.place_tokens.get(place_idx, np.zeros(num_classes))
                self.place_token_time[place_idx] = self.place_token_time.get(
                    place_idx, np.zeros(num_classes)
                ) + tokens * elapsed
                self.last_place_update[place_idx] = current_time

        # Produce output tokens and decrement in-service count per class
        num_nodes = self.config.num_nodes
        for node_idx in range(num_nodes):
            for class_idx in range(num_classes):
                produced = int(mode.firing_outcomes[node_idx, class_idx])
                if produced > 0:
                    # Check if destination is a place
                    place_idx = self.place_node_to_idx.get(node_idx)
                    if place_idx is not None:
                        self.place_tokens[place_idx][class_idx] += produced
                    # Decrement in-service count for this class
                    if trans_idx in self.transition_in_service_count:
                        self.transition_in_service_count[trans_idx][class_idx] -= produced
                        # Ensure non-negative
                        if self.transition_in_service_count[trans_idx][class_idx] < 0:
                            self.transition_in_service_count[trans_idx][class_idx] = 0

        # Trigger cascade - check if more transitions can fire
        self._check_and_fire_transitions()

        # Increment event count
        self.total_event_count += 1

    def _sample_firing_time(self, trans_idx: int, mode: TransitionModeInfo) -> float:
        """
        Sample firing time from the mode's distribution.

        Supports EXP, ERLANG, HYPEREXP, and general PH distributions.
        Falls back to Exp(1) if no distribution is specified.
        """
        proc_id = mode.firing_proc_id.upper() if mode.firing_proc_id else 'EXP'

        try:
            if proc_id == 'EXP' or proc_id == 'IMMEDIATE':
                # Exponential distribution: sample Exp(mean)
                return self.rng.exponential(mode.firing_mean)

            elif proc_id == 'ERLANG':
                # Erlang(mean, k): sum of k independent Exp(k/mean) random variables
                k = mode.firing_phases
                if k <= 0:
                    k = 1
                rate = k / mode.firing_mean if mode.firing_mean > 0 else 1.0
                return sum(self.rng.exponential(1.0 / rate) for _ in range(k))

            elif proc_id == 'HYPEREXP' or proc_id == 'HEXP':
                # HyperExponential: mixture of exponentials
                # Extract rates from D0 diagonal and compute probabilities from D1
                # This matches the Kotlin implementation which computes p = D1[0,0] / lambda1
                if mode.firing_d0 is not None:
                    rates = -np.diag(mode.firing_d0)
                    if len(rates) >= 2 and np.all(rates > 0):
                        lambda1 = rates[0]
                        lambda2 = rates[1]
                        # Extract probability from D1 (as in Kotlin: p = D1[0,0] / lambda1)
                        if mode.firing_d1 is not None and mode.firing_d1.size > 0:
                            p = mode.firing_d1[0, 0] / lambda1 if lambda1 > 0 else 0.5
                        elif mode.firing_pie is not None and len(mode.firing_pie) > 0:
                            # Fallback to pie if D1 not available
                            p = mode.firing_pie[0]
                        else:
                            p = 0.5
                        # Clamp probability to valid range
                        p = np.clip(p, 0.0, 1.0)
                        probs = np.array([p, 1.0 - p])
                        # Choose phase according to probabilities
                        phase = self.rng.choice(2, p=probs)
                        return self.rng.exponential(1.0 / rates[phase])
                # Fallback: use mean
                return self.rng.exponential(mode.firing_mean)

            elif proc_id == 'PH' or proc_id == 'APH' or proc_id == 'COXIAN':
                # Phase-type distribution: simulate absorption time
                return self._sample_ph_distribution(mode.firing_d0, mode.firing_pie)

            elif proc_id == 'DET' or proc_id == 'DETERMINISTIC':
                # Deterministic: fixed time
                return mode.firing_mean

            else:
                # Unknown type: use exponential with the mean
                return self.rng.exponential(mode.firing_mean)

        except Exception:
            # Fallback: Exp(mean)
            return self.rng.exponential(mode.firing_mean if mode.firing_mean > 0 else 1.0)

    def _sample_ph_distribution(
        self, d0: Optional[np.ndarray], pie: Optional[np.ndarray]
    ) -> float:
        """
        Sample from a phase-type distribution using the D0 matrix and initial probabilities.

        The PH distribution is the time until absorption in a Markov chain.
        D0 is the sub-generator (diagonal = -exit rates, off-diagonal = transition rates).
        """
        if d0 is None:
            return self.rng.exponential(1.0)

        n = d0.shape[0]
        if n == 0:
            return self.rng.exponential(1.0)

        # Normalize initial probabilities
        if pie is not None:
            probs = pie.flatten()
            probs = probs / np.sum(probs) if np.sum(probs) > 0 else np.ones(n) / n
        else:
            probs = np.ones(n) / n

        # Start in initial phase
        current_phase = self.rng.choice(n, p=probs)
        total_time = 0.0

        while True:
            # Exit rate from current phase
            exit_rate = -d0[current_phase, current_phase]
            if exit_rate <= 0:
                break

            # Time in current phase
            time_in_phase = self.rng.exponential(1.0 / exit_rate)
            total_time += time_in_phase

            # Transition probabilities to other phases
            # Probability of going to phase j: d0[current, j] / exit_rate (for j != current)
            # Probability of absorption: 1 - sum(transition probs)
            trans_probs = np.zeros(n + 1)  # +1 for absorption state
            for j in range(n):
                if j != current_phase and d0[current_phase, j] > 0:
                    trans_probs[j] = d0[current_phase, j] / exit_rate
            trans_probs[n] = 1.0 - np.sum(trans_probs[:n])  # Absorption
            trans_probs = np.clip(trans_probs, 0, 1)
            trans_probs = trans_probs / np.sum(trans_probs)  # Normalize

            # Choose next state
            next_state = self.rng.choice(n + 1, p=trans_probs)
            if next_state == n:  # Absorbed
                break
            current_phase = next_state

            # Safety limit to prevent infinite loops
            if total_time > 1e6:
                break

        return total_time

    # ==================== End SPN Core Methods ====================

    def _arrival_process(self, src_idx: int, class_id: int):
        """
        SimPy process for external arrivals from a source.

        Replaces Java's ExternalArrival event.
        """
        while True:
            # Generate inter-arrival time
            arrival_rate = self.config.lambdas[src_idx, class_id]
            if arrival_rate <= 0:
                break

            interarrival_time = self.rng.exponential(1.0 / arrival_rate)
            yield self.env.timeout(interarrival_time)

            # Create customer
            customer = Customer(
                class_id=class_id,
                priority=self.config.class_priorities[class_id],
                system_arrival_time=self.env.now,
                queue_arrival_time=self.env.now,
                random_rank=self.rng.random(),
                job_id=self._get_next_job_id(),
            )

            # Route to first destination
            src_node = self.config.source_nodes[src_idx]
            dest_node, dest_class = self._route_from_node(src_node, class_id)

            if dest_node is not None:
                customer.class_id = dest_class
                self._arrive_at_node(dest_node, customer)

    def _arrive_at_node(self, node_idx: int, customer: Customer) -> None:
        """Handle customer arrival at any node type."""
        if node_idx in self.config.sink_nodes:
            self._depart_system(customer)
        elif node_idx in self.config.service_nodes:
            svc_idx = self.config.service_nodes.index(node_idx)
            self._arrive_at_service_node(svc_idx, customer)
        elif node_idx in self.config.router_nodes:
            # Pass through router
            dest_node, dest_class = self._route_from_node(node_idx, customer.class_id)
            if dest_node is not None:
                customer.class_id = dest_class
                self._arrive_at_node(dest_node, customer)
        elif node_idx in self.config.cache_nodes:
            # Process cache access (hit or miss) and switch class accordingly
            new_class = self._process_cache_access(node_idx, customer.class_id)
            customer.class_id = new_class
            # Route to next destination
            dest_node, dest_class = self._route_from_node(node_idx, new_class)
            if dest_node is not None:
                customer.class_id = dest_class
                self._arrive_at_node(dest_node, customer)
        elif node_idx in self.config.class_switch_nodes:
            # Pass through class switch
            dest_node, dest_class = self._route_from_node(node_idx, customer.class_id)
            if dest_node is not None:
                customer.class_id = dest_class
                self._arrive_at_node(dest_node, customer)
        elif node_idx in self.config.place_nodes:
            # Deposit token in place and check transitions
            place_idx = self.place_node_to_idx.get(node_idx)
            if place_idx is not None:
                self.place_tokens[place_idx][customer.class_id] += 1
                self._check_and_fire_transitions()
        elif node_idx in self.config.fork_nodes:
            # Phase 5: Handle Fork node (parent disappears, children created)
            self._handle_fork_arrival(node_idx, customer)
            # Parent absorbed - no further routing
        elif node_idx in self.config.join_nodes:
            # Phase 5: Handle Join node (forked children synchronize)
            self._handle_join_arrival(node_idx, customer)
            # Join handles routing after sync

    def _arrive_at_service_node(self, svc_idx: int, customer: Customer) -> None:
        """Handle arrival at a service node (Queue or Delay)."""
        class_id = customer.class_id
        current_time = self.env.now

        # Update queue arrival time
        customer.queue_arrival_time = current_time

        # Record arrival
        self.stats[svc_idx].record_arrival(class_id)

        # Check capacity
        current_length = (
            sum(self.stats[svc_idx].current_queue_length) +
            sum(self.stats[svc_idx].current_busy_servers)
        )
        buffer_cap = self.config.buffer_capacities.get(svc_idx, 10000000)

        if current_length >= buffer_cap:
            # Drop customer
            self.stats[svc_idx].record_drop(class_id)
            return

        # Phase 6: Check balking (customer refuses to join based on queue length)
        if self._should_balk(svc_idx, class_id, current_length):
            self.stats[svc_idx].record_drop(class_id)
            return

        # Update queue stats
        self.stats[svc_idx].update_queue(class_id, current_time, delta=1)

        # Check for available server
        is_delay = self.config.is_delay_node[svc_idx] if svc_idx < len(self.config.is_delay_node) else False

        if is_delay:
            # Delay node: infinite servers, start service immediately
            self.env.process(self._service_process(svc_idx, -1, customer))
        else:
            # Phase 6: Check for heterogeneous servers
            n_types = self.config.num_server_types.get(svc_idx, 0)

            if n_types > 0:
                # Heterogeneous queue: use class-aware server selection
                selection = self._find_available_server_for_class(svc_idx, class_id)
                if selection.server_id >= 0:
                    # Server found - start service with heterogeneous rate
                    self.env.process(self._service_process(
                        svc_idx, selection.server_id, customer, selection.server_type_id
                    ))
                else:
                    # No compatible server available - queue the customer
                    scheduler = self.schedulers.get(svc_idx)
                    if scheduler:
                        service_gen = self._create_service_generator(svc_idx)
                        scheduler.arrive(customer, current_time, service_gen)
                    else:
                        if svc_idx not in self.wait_queues:
                            self.wait_queues[svc_idx] = []
                        self.wait_queues[svc_idx].append(customer)
            else:
                # Phase 3: Homogeneous queue - use scheduler for queue management
                scheduler = self.schedulers[svc_idx]
                service_gen = self._create_service_generator(svc_idx)

                # Delegate arrival to scheduler
                accepted, service_info = scheduler.arrive(customer, current_time, service_gen)

                if not accepted:
                    # Customer rejected (shouldn't happen with current schedulers)
                    self.stats[svc_idx].record_drop(class_id)
                    self.stats[svc_idx].update_queue(class_id, current_time, delta=-1)
                    return

                if service_info is not None:
                    # Phase 4: Check if preemptive scheduler returned preemption info
                    if scheduler.is_preemptive() and isinstance(service_info, tuple) and len(service_info) == 3:
                        server_id, new_customer, preempted_info = service_info

                        if preempted_info is not None:
                            # Handle preemption: interrupt the active process
                            self._handle_preemption(svc_idx, server_id, preempted_info)

                        # Start new customer on the server
                        self.env.process(self._service_process(svc_idx, server_id, new_customer))
                    else:
                        # Non-preemptive or standard arrival
                        server_id, customer_to_serve = service_info if isinstance(service_info, tuple) else (service_info, customer)
                        self.env.process(self._service_process(svc_idx, server_id, customer_to_serve))
                # else: customer is queued inside scheduler (waiting for service)

    def _find_available_server(self, svc_idx: int) -> int:
        """Find an available server. Returns -1 if none."""
        if svc_idx not in self.server_busy:
            return -1
        for i, busy in enumerate(self.server_busy[svc_idx]):
            if not busy:
                return i
        return -1

    def _find_available_server_for_class(self, svc_idx: int, class_id: int) -> ServerSelection:
        """
        Find an available server that can serve the given job class.

        For heterogeneous queues, applies the configured scheduling policy
        (FSF, ALIS, ALFS, FAIRNESS, ORDER, RAIS) to select among compatible
        server types with available capacity.

        Args:
            svc_idx: Service node index
            class_id: Job class ID

        Returns:
            ServerSelection with server_id and server_type_id (-1 if none available)
        """
        n_types = self.config.num_server_types.get(svc_idx, 0)

        # Homogeneous case: use simple server selection
        if n_types <= 0:
            server_id = self._find_available_server(svc_idx)
            return ServerSelection(server_id=server_id, server_type_id=-1)

        # Heterogeneous case: find compatible server type with capacity
        compat = self.config.server_compat.get(svc_idx, [])
        servers_per_type = self.config.servers_per_type.get(svc_idx, [])
        busy_per_type = self.busy_count_per_type.get(svc_idx, [])
        server_to_type = self.config.server_to_type.get(svc_idx, [])
        policy = self.config.hetero_sched_policy.get(svc_idx, HeteroSchedPolicy.ORDER)

        # Get compatible types with available capacity
        available_types = []
        for type_id in range(n_types):
            if type_id < len(compat) and class_id < len(compat[type_id]):
                if compat[type_id][class_id]:  # Compatible with this class
                    if type_id < len(servers_per_type) and type_id < len(busy_per_type):
                        if busy_per_type[type_id] < servers_per_type[type_id]:  # Has capacity
                            available_types.append(type_id)

        if not available_types:
            return ServerSelection(server_id=-1, server_type_id=-1)

        # Apply scheduling policy to select server type
        selected_type_id = self._select_server_type(svc_idx, class_id, available_types, policy)

        if selected_type_id < 0:
            return ServerSelection(server_id=-1, server_type_id=-1)

        # Find a free server of the selected type
        for server_id, busy in enumerate(self.server_busy.get(svc_idx, [])):
            if not busy:
                if server_id < len(server_to_type) and server_to_type[server_id] == selected_type_id:
                    return ServerSelection(server_id=server_id, server_type_id=selected_type_id)

        return ServerSelection(server_id=-1, server_type_id=-1)

    def _select_server_type(
        self, svc_idx: int, class_id: int, available_types: List[int], policy: HeteroSchedPolicy
    ) -> int:
        """
        Select a server type according to the heterogeneous scheduling policy.

        Args:
            svc_idx: Service node index
            class_id: Job class ID
            available_types: List of type IDs with compatible capacity
            policy: Heterogeneous scheduling policy

        Returns:
            Selected server type ID, or -1 if none
        """
        if not available_types:
            return -1

        if policy == HeteroSchedPolicy.ORDER:
            # Order-based: select first available type (by type ID order)
            return min(available_types)

        elif policy == HeteroSchedPolicy.FSF:
            # Fastest Server First: select type with highest service rate for this class
            hetero_mus = self.config.hetero_mus.get(svc_idx, [])
            best_type = available_types[0]
            best_rate = 0.0
            for type_id in available_types:
                if type_id < len(hetero_mus) and class_id < len(hetero_mus[type_id]):
                    rate = hetero_mus[type_id][class_id]
                    if rate > best_rate:
                        best_rate = rate
                        best_type = type_id
            return best_type

        elif policy == HeteroSchedPolicy.ALIS:
            # Assign to Longest Idle Server: round-robin among available types
            type_order = self.server_type_order.get(svc_idx, list(range(len(available_types))))
            for type_id in type_order:
                if type_id in available_types:
                    # Move this type to end of order (round-robin)
                    type_order.remove(type_id)
                    type_order.append(type_id)
                    return type_id
            return available_types[0]

        elif policy == HeteroSchedPolicy.FAIRNESS:
            # Fairness: round-robin for balanced utilization (same as ALIS)
            type_order = self.server_type_order.get(svc_idx, list(range(len(available_types))))
            for type_id in type_order:
                if type_id in available_types:
                    type_order.remove(type_id)
                    type_order.append(type_id)
                    return type_id
            return available_types[0]

        elif policy == HeteroSchedPolicy.ALFS:
            # Assign to Least Flexible Server: prefer specialized servers
            alfs_order = self.alfs_order.get(svc_idx, list(range(len(available_types))))
            for type_id in alfs_order:
                if type_id in available_types:
                    return type_id
            return available_types[0]

        elif policy == HeteroSchedPolicy.RAIS:
            # Random Assignment with Idle Selection
            import random
            return random.choice(available_types)

        else:
            # Default to ORDER
            return min(available_types)

    def _mark_server_busy(self, svc_idx: int, server_id: int, server_type_id: int) -> None:
        """
        Mark a server as busy and update per-type busy counts.

        Args:
            svc_idx: Service node index
            server_id: Global server ID
            server_type_id: Server type ID (-1 for homogeneous)
        """
        if svc_idx in self.server_busy and server_id < len(self.server_busy[svc_idx]):
            self.server_busy[svc_idx][server_id] = True

        if server_type_id >= 0 and svc_idx in self.busy_count_per_type:
            if server_type_id < len(self.busy_count_per_type[svc_idx]):
                self.busy_count_per_type[svc_idx][server_type_id] += 1

    def _mark_server_free(self, svc_idx: int, server_id: int, server_type_id: int) -> None:
        """
        Mark a server as free and update per-type busy counts.

        Args:
            svc_idx: Service node index
            server_id: Global server ID
            server_type_id: Server type ID (-1 for homogeneous)
        """
        if svc_idx in self.server_busy and server_id < len(self.server_busy[svc_idx]):
            self.server_busy[svc_idx][server_id] = False

        if server_type_id >= 0 and svc_idx in self.busy_count_per_type:
            if server_type_id < len(self.busy_count_per_type[svc_idx]):
                self.busy_count_per_type[svc_idx][server_type_id] = max(
                    0, self.busy_count_per_type[svc_idx][server_type_id] - 1
                )

    def _generate_hetero_service_time(self, svc_idx: int, class_id: int, server_type_id: int) -> float:
        """
        Generate service time using heterogeneous service rate.

        For heterogeneous queues, uses the per-(type, class) service rate.
        Falls back to homogeneous rate if hetero rate not available.
        Applies load-dependent scaling if configured.

        Args:
            svc_idx: Service node index
            class_id: Job class ID
            server_type_id: Server type ID (-1 for homogeneous)

        Returns:
            Service time (exponentially distributed, scaled for load dependence)
        """
        service_time = 0.0

        # Try heterogeneous rate first
        if server_type_id >= 0:
            hetero_mus = self.config.hetero_mus.get(svc_idx, [])
            if server_type_id < len(hetero_mus) and class_id < len(hetero_mus[server_type_id]):
                rate = hetero_mus[server_type_id][class_id]
                if rate > 0 and np.isfinite(rate):
                    service_time = self.rng.exponential(1.0 / rate)

        # Fall back to homogeneous rate if no hetero rate
        if service_time == 0.0:
            if self.config.mus is not None and svc_idx < self.config.mus.shape[0]:
                rate = self.config.mus[svc_idx, class_id]
                if rate > 0 and np.isfinite(rate):
                    service_time = self.rng.exponential(1.0 / rate)

        # Apply load-dependent scaling if this station has load dependence
        if self.config.is_load_dependent.get(svc_idx, False) and service_time > 0:
            scaling_factor = self._get_load_dependent_scaling(svc_idx)
            if scaling_factor > 0:
                # Scaling factor multiplies the rate, so divide service time
                service_time /= scaling_factor

        return service_time

    def _get_load_dependent_scaling(self, svc_idx: int) -> float:
        """
        Get the load-dependent scaling factor for a station.

        The scaling factor is based on the current number of jobs at the station.
        lldscaling[station, n-1] gives the scaling factor when n jobs are present.

        Args:
            svc_idx: Service node index

        Returns:
            Scaling factor (1.0 if not load-dependent or out of range)
        """
        if svc_idx not in self.config.lld_scaling:
            return 1.0

        # Get current total jobs at this station (queue + in service)
        total_jobs = 0
        if svc_idx in self.stats:
            total_jobs = (
                int(sum(self.stats[svc_idx].current_queue_length)) +
                int(sum(self.stats[svc_idx].current_busy_servers))
            )

        if total_jobs <= 0:
            return 1.0  # First job gets full rate

        scaling_array = self.config.lld_scaling[svc_idx]
        # lldscaling[n-1] gives scaling when n jobs present
        scaling_idx = min(total_jobs - 1, len(scaling_array) - 1)
        return scaling_array[scaling_idx]

    def _handle_preemption(self, svc_idx: int, server_id: int, preempted_info: Any) -> None:
        """
        Phase 4: Handle preemption of a job currently in service.

        Interrupts the SimPy process for the preempted customer.
        The scheduler has already updated its internal state (saved preemption
        record, decremented busy servers, queued the preempted job, etc.).

        Args:
            svc_idx: Service node index
            server_id: Server ID where preemption occurred
            preempted_info: Tuple of (preempted_customer, preemption_record)
        """
        if server_id in self.active_processes.get(svc_idx, {}):
            process = self.active_processes[svc_idx][server_id]
            if process and process.is_alive:
                try:
                    process.interrupt()  # Raises simpy.Interrupt in the process
                except Exception:
                    pass  # Process may have already completed

            # Clean up process tracking
            del self.active_processes[svc_idx][server_id]

        if server_id in self.active_customers.get(svc_idx, {}):
            del self.active_customers[svc_idx][server_id]

    def _handle_fork_arrival(self, fork_node_idx: int, customer: Customer) -> None:
        """
        Phase 5: Handle parent job arrival at Fork node.

        Creates forked children and routes them to parallel branches.
        Parent job is absorbed (ceases to exist) at the Fork.

        Args:
            fork_node_idx: Fork node index
            customer: Parent customer (job to be forked)
        """
        from .nodes.fork_join import ForkChild

        fork = self.fj_manager.get_fork(fork_node_idx)
        if fork is None:
            return

        current_time = self.env.now

        # Fork the customer into parallel children
        children = fork.fork(customer, current_time)

        # Route each child to its destination
        for child in children:
            child_customer = child.customer
            child_customer.parent_job_id = child.parent_job_id

            # Route to destination
            self._arrive_at_node(child.destination_node, child_customer)

    def _handle_join_arrival(self, join_node_idx: int, customer: Customer) -> None:
        """
        Phase 5: Handle forked child arrival at Join node.

        Synchronizes forked children and routes parent when quorum is met.

        Args:
            join_node_idx: Join node index
            customer: Forked child customer
        """
        from .nodes.fork_join import ForkChild

        join = self.fj_manager.get_join(join_node_idx)
        if join is None:
            return

        current_time = self.env.now
        join_idx = self.config.join_nodes.index(join_node_idx) if join_node_idx in self.config.join_nodes else 0

        # Update Join buffer queue length stats (increment)
        if self.warmup_done:
            self.join_buffer_stats[join_idx].update_queue(customer.class_id, current_time, delta=1)
            self.join_buffer_stats[join_idx].record_arrival(customer.class_id)

        # Create ForkChild from arriving customer
        child = ForkChild(
            child_id=customer.customer_id,
            parent_job_id=customer.parent_job_id,
            task_index=0,
            destination_node=join_node_idx,
            customer=customer
        )

        # Check synchronization condition
        parent_customer = join.arrive(child, current_time)

        if parent_customer is not None:
            # Synchronization complete - all required tasks have arrived!

            # Update Join buffer stats (decrement for all forked tasks)
            if self.warmup_done:
                # Get fork info to know how many tasks to decrement
                fork_info = None
                for fork_node_idx in self.config.fork_nodes:
                    fork = self.fj_manager.get_fork(fork_node_idx)
                    if fork and customer.parent_job_id in fork.fork_jobs:
                        fork_info = fork.fork_jobs[customer.parent_job_id]
                        break

                if fork_info:
                    # Record response time for each forked task
                    for task_idx in range(fork_info.total_tasks):
                        self.join_buffer_stats[join_idx].record_completion(customer.class_id, 0.0)

                    # Decrement queue length by total tasks
                    self.join_buffer_stats[join_idx].update_queue(
                        customer.class_id, current_time, delta=-fork_info.total_tasks
                    )

            # Route parent job from Join to next destination
            dest_node, dest_class = self._route_from_node(join_node_idx, parent_customer.class_id)
            if dest_node is not None:
                parent_customer.class_id = dest_class
                self._arrive_at_node(dest_node, parent_customer)

        # If parent_customer is None, forked job absorbed at Join (waiting for sync)

    def _can_server_serve(self, svc_idx: int, server_id: int) -> bool:
        """
        Phase 6b: Check if server can serve customers (must be in ACTIVE state).

        Args:
            svc_idx: Service node index
            server_id: Server ID

        Returns:
            True if server is ACTIVE, False otherwise
        """
        if svc_idx not in self.server_states:
            return True  # No state tracking (infinite servers)

        state = self.server_states[svc_idx].get(server_id, ServerState.ACTIVE)
        return state == ServerState.ACTIVE

    def _transition_server_to_setup(self, svc_idx: int, server_id: int) -> None:
        """
        Phase 6b: Transition server from OFF to SETUP (cold start).

        Args:
            svc_idx: Service node index
            server_id: Server ID
        """
        if svc_idx not in self.server_states or server_id not in self.server_states[svc_idx]:
            return

        if self.server_states[svc_idx][server_id] == ServerState.OFF:
            self.server_states[svc_idx][server_id] = ServerState.SETUP
            # Setup time would be sampled from distribution (deferred to extension)
            # For now, setup is instantaneous
            self._transition_server_to_active(svc_idx, server_id)

    def _transition_server_to_active(self, svc_idx: int, server_id: int) -> None:
        """
        Phase 6b: Transition server from SETUP to ACTIVE.

        Args:
            svc_idx: Service node index
            server_id: Server ID
        """
        if svc_idx not in self.server_states or server_id not in self.server_states[svc_idx]:
            return

        if self.server_states[svc_idx][server_id] == ServerState.SETUP:
            self.server_states[svc_idx][server_id] = ServerState.ACTIVE
            # Server is now ready to serve customers

    def _initiate_server_delayoff(self, svc_idx: int, server_id: int) -> None:
        """
        Phase 6b: Initiate server shutdown after delay (for FaaS cold shutdown).

        Transitions server from ACTIVE to DELAYOFF, then OFF after delay.

        Args:
            svc_idx: Service node index
            server_id: Server ID
        """
        if svc_idx not in self.server_states or server_id not in self.server_states[svc_idx]:
            return

        if self.server_states[svc_idx][server_id] != ServerState.ACTIVE:
            return  # Can only delayoff from ACTIVE state

        self.server_states[svc_idx][server_id] = ServerState.DELAYOFF

        # Schedule server shutdown after delay (deferred to extension with configurable delays)
        # For now, shutdown is immediate
        self._shutdown_server(svc_idx, server_id)

    def _cancel_server_delayoff(self, svc_idx: int, server_id: int) -> None:
        """
        Phase 6b: Cancel scheduled server shutdown (customer arrival during DELAYOFF).

        Transitions server from DELAYOFF back to ACTIVE.

        Args:
            svc_idx: Service node index
            server_id: Server ID
        """
        if svc_idx not in self.server_states or server_id not in self.server_states[svc_idx]:
            return

        if self.server_states[svc_idx][server_id] == ServerState.DELAYOFF:
            self.server_states[svc_idx][server_id] = ServerState.ACTIVE
            # Cancel the scheduled shutdown event if it exists
            key = (svc_idx, server_id)
            if key in self.delayoff_events:
                event = self.delayoff_events[key]
                if not event.processed:
                    event.cancel()
                del self.delayoff_events[key]

    def _shutdown_server(self, svc_idx: int, server_id: int) -> None:
        """
        Phase 6b: Shut down server (transition to OFF).

        Args:
            svc_idx: Service node index
            server_id: Server ID
        """
        if svc_idx not in self.server_states or server_id not in self.server_states[svc_idx]:
            return

        self.server_states[svc_idx][server_id] = ServerState.OFF

    def _should_balk(self, svc_idx: int, class_id: int, queue_length: int) -> bool:
        """
        Phase 6: Check if customer should balk (refuse to join queue).

        Default: No balking. Can be overridden with balking parameters.

        Args:
            svc_idx: Service node index
            class_id: Customer class
            queue_length: Current queue length

        Returns:
            True if customer should balk, False otherwise
        """
        # Default: no balking (balking parameters not currently exposed from network config)
        return False

    def _reneging_process(self, svc_idx: int, customer: Customer) -> None:
        """
        Phase 7a: SimPy process for customer reneging (abandonment).

        Removes customer from queue if patience time expires before service starts.

        Args:
            svc_idx: Service node index
            customer: Customer who may renege
        """
        import simpy

        try:
            # Wait for patience time
            yield self.env.timeout(customer.patience_time)

            # Check if still waiting (not in service)
            if svc_idx in self.wait_queues:
                wait_queue = self.wait_queues[svc_idx]
                if customer in wait_queue:
                    # Customer reneges - remove from queue
                    customer.has_reneged = True
                    wait_queue.remove(customer)

                    # Update statistics
                    self.stats[svc_idx].update_queue(customer.class_id, self.env.now, delta=-1)
                    self.stats[svc_idx].record_drop(customer.class_id)

                    # Attempt retrial if configured
                    if customer.patience_time < float('inf'):
                        # Simple retrial: schedule re-attempt after delay
                        retrial_delay = self.rng.exponential(1.0)  # Default: exponential with rate 1
                        self.env.process(self._retrial_process(svc_idx, customer, retrial_delay))

        except simpy.Interrupt:
            # Customer started service before patience expired
            pass

    def _retrial_process(self, svc_idx: int, customer: Customer, delay: float) -> None:
        """
        Phase 7b: SimPy process for customer retrials.

        Schedules a re-attempt after initial rejection or reneging.

        Args:
            svc_idx: Service node index
            customer: Customer attempting retrial
            delay: Delay before retry
        """
        # Wait for retrial delay
        yield self.env.timeout(delay)

        # Re-attempt arrival
        self._arrive_at_service_node(svc_idx, customer)

    def _service_process(
        self, svc_idx: int, server_id: int, customer: Customer, server_type_id: int = -1
    ):
        """
        SimPy process for serving a customer.

        Phase 3: Supports pre-sampled service times for job-based scheduling.
        Phase 4: Supports preemption via simpy.Interrupt.
        Phase 6: Supports heterogeneous servers with per-type service rates.

        Replaces Java's ServiceCompletion/DelayCompletion events.

        Args:
            svc_idx: Service node index
            server_id: Global server ID (-1 for delay nodes)
            customer: Customer being served
            server_type_id: Server type ID for heterogeneous queues (-1 for homogeneous)
        """
        import simpy

        class_id = customer.class_id
        current_time = self.env.now
        scheduler = self.schedulers.get(svc_idx)

        # Mark server busy with heterogeneous support
        if server_id >= 0:
            self._mark_server_busy(svc_idx, server_id, server_type_id)

        # Phase 4: Track process reference for preemption
        if scheduler and scheduler.is_preemptive():
            self.active_processes[svc_idx][server_id] = self.env.active_process
            self.active_customers[svc_idx][server_id] = customer

        # Update busy stats
        self.stats[svc_idx].update_busy(class_id, current_time, delta=1)

        # Phase 3: Use pre-sampled service time if available
        if customer.service_time >= 0:
            service_time = customer.service_time
        else:
            # Phase 6: Use heterogeneous service time generation
            service_time = self._generate_hetero_service_time(svc_idx, class_id, server_type_id)

        try:
            # Wait for service completion
            service_start = self.env.now
            yield self.env.timeout(service_time)

            # === Normal completion path ===
            current_time = self.env.now

            # Mark server free with heterogeneous support
            if server_id >= 0:
                self._mark_server_free(svc_idx, server_id, server_type_id)

            # Update stats
            self.stats[svc_idx].update_queue(class_id, current_time, delta=-1)
            self.stats[svc_idx].update_busy(class_id, current_time, delta=-1)

            # Record completion
            response_time = current_time - customer.queue_arrival_time
            self.stats[svc_idx].record_completion(class_id, response_time)
            self.total_event_count += 1

            # Check if MSER sample should be collected
            if self.mser_enabled:
                events_since_last = self.total_event_count - self.last_mser_event_count
                if events_since_last >= self.mser_observation_interval:
                    self._collect_mser_sample()

            # Phase 4: Clean up process tracking
            if scheduler and scheduler.is_preemptive():
                if server_id in self.active_processes.get(svc_idx, {}):
                    del self.active_processes[svc_idx][server_id]
                if server_id in self.active_customers.get(svc_idx, {}):
                    del self.active_customers[svc_idx][server_id]

            # Route customer to next destination
            node_idx = self.config.service_nodes[svc_idx]
            dest_node, dest_class = self._route_from_node(node_idx, class_id)

            if dest_node is not None:
                customer.class_id = dest_class
                self._arrive_at_node(dest_node, customer)

            # Phase 3: Use scheduler to get next customer
            # For heterogeneous queues, we need to find compatible server for next customer
            if scheduler:
                next_customer = scheduler.on_departure(customer, current_time, server_id)
                if next_customer is not None:
                    # For heterogeneous: check if this server type is compatible with next customer
                    next_type_id = server_type_id
                    if server_type_id >= 0:
                        compat = self.config.server_compat.get(svc_idx, [])
                        if (server_type_id < len(compat) and
                            next_customer.class_id < len(compat[server_type_id]) and
                            not compat[server_type_id][next_customer.class_id]):
                            # Server not compatible with next customer, try to find another
                            next_selection = self._find_available_server_for_class(
                                svc_idx, next_customer.class_id
                            )
                            if next_selection.server_id >= 0:
                                self.env.process(self._service_process(
                                    svc_idx, next_selection.server_id, next_customer,
                                    next_selection.server_type_id
                                ))
                            # else: customer stays in queue
                        else:
                            self.env.process(self._service_process(
                                svc_idx, server_id, next_customer, server_type_id
                            ))
                    else:
                        self.env.process(self._service_process(
                            svc_idx, server_id, next_customer, server_type_id
                        ))
            # Fallback: Check manual wait queue (backward compatibility)
            elif server_id >= 0 and self.wait_queues.get(svc_idx):
                next_customer = self.wait_queues[svc_idx].pop(0)
                self.env.process(self._service_process(
                    svc_idx, server_id, next_customer, server_type_id
                ))

        except simpy.Interrupt:
            # === Preemption path ===
            # Job was preempted before completion
            elapsed = self.env.now - service_start

            # Mark server free with heterogeneous support
            if server_id >= 0:
                self._mark_server_free(svc_idx, server_id, server_type_id)

            # Update busy stats (no longer in service)
            self.stats[svc_idx].update_busy(class_id, self.env.now, delta=-1)

            # Phase 4: Clean up process tracking
            if server_id in self.active_processes.get(svc_idx, {}):
                del self.active_processes[svc_idx][server_id]
            if server_id in self.active_customers.get(svc_idx, {}):
                del self.active_customers[svc_idx][server_id]

            # NOTE: NO completion recorded, NO routing, NO next customer
            # The preempting job already started service in the arrival handler
            # The scheduler has saved preemption state for later resume

    def _depart_system(self, customer: Customer) -> None:
        """Handle customer departure from system (arrival at Sink)."""
        class_id = customer.class_id
        system_response_time = self.env.now - customer.system_arrival_time

        self.system_completed[class_id] += 1
        self.system_response_time_sum[class_id] += system_response_time

    def _collect_mser_sample(self) -> None:
        """
        Collect MSER-5 observation sample.

        Records time-weighted average queue length over the interval since
        the last sample. This mirrors Java's collectMSERSample() method.
        """
        current_time = self.env.now
        interval_duration = current_time - self.last_mser_sample_time

        if interval_duration <= 0:
            return

        # Update all queue stats before sampling
        for svc_idx in range(len(self.config.service_nodes)):
            for k in range(self.config.num_classes):
                self.stats[svc_idx].update_queue(k, current_time)

        self.mser_observation_times.append(current_time)

        # Record time-weighted average queue length over this interval
        for svc_idx in range(len(self.config.service_nodes)):
            for k in range(self.config.num_classes):
                # Compute interval queue time
                current_queue_time = self.stats[svc_idx].total_queue_time[k]
                last_queue_time = self.last_mser_queue_time[svc_idx][k]
                interval_queue_time = current_queue_time - last_queue_time

                # Average queue length over this interval
                avg_queue_length = interval_queue_time / interval_duration
                self.mser_observations[svc_idx][k].append(avg_queue_length)

                # Update snapshot
                self.last_mser_queue_time[svc_idx][k] = current_queue_time

        self.last_mser_sample_time = current_time
        self.last_mser_event_count = self.total_event_count

    def _apply_mser5_truncation(self) -> None:
        """
        Apply MSER-5 to determine warmup truncation point.

        Uses the aggregated queue length observations to find the optimal
        truncation point that minimizes variance.
        """
        # Aggregate observations across all queues/classes for truncation detection
        all_observations: List[float] = []

        # Collect all observations in time order
        num_obs = len(self.mser_observation_times)
        for obs_idx in range(num_obs):
            total_qlen = 0.0
            count = 0
            for svc_idx in range(len(self.config.service_nodes)):
                for k in range(self.config.num_classes):
                    if obs_idx < len(self.mser_observations[svc_idx][k]):
                        total_qlen += self.mser_observations[svc_idx][k][obs_idx]
                        count += 1
            if count > 0:
                all_observations.append(total_qlen / count)

        if len(all_observations) < 20:  # Not enough data
            self.warmup_end_time = 0.0
            return

        # Apply MSER-5 detector
        detector = MSER5TransientDetector(batch_size=self.mser_batch_size, min_batches=10)
        truncation_idx = detector.detect_warmup(all_observations)

        if truncation_idx > 0 and truncation_idx < len(self.mser_observation_times):
            self.warmup_end_time = self.mser_observation_times[truncation_idx]
        else:
            self.warmup_end_time = 0.0

    def _route_from_node(self, node_idx: int, class_id: int) -> Tuple[Optional[int], int]:
        """
        Determine next destination from routing matrix using the appropriate routing strategy.

        Supports PROB (probabilistic), RAND (random uniform), RROBIN (round-robin),
        WRROBIN (weighted round-robin), JSQ (join shortest queue), and KCHOICES
        (power of K choices) routing strategies.

        Returns:
            Tuple of (destination_node_idx, destination_class_id)
        """
        if self.config.routing_matrix is None:
            return (None, class_id)

        num_classes = self.config.num_classes
        num_nodes = self.config.num_nodes

        # Row in routing matrix: node_idx * num_classes + class_id
        from_idx = node_idx * num_classes + class_id

        # Get probabilities for all destinations
        probs = []
        destinations = []

        for to_node in range(num_nodes):
            for to_class in range(num_classes):
                to_idx = to_node * num_classes + to_class
                if (from_idx < self.config.routing_matrix.shape[0] and
                    to_idx < self.config.routing_matrix.shape[1]):
                    p = self.config.routing_matrix[from_idx, to_idx]
                    if p > 0:
                        probs.append(p)
                        destinations.append((to_node, to_class))

        if not destinations:
            return (None, class_id)

        # If only one destination, return it directly
        if len(destinations) == 1:
            return destinations[0]

        # Get routing strategy for this node and class
        strategy = RoutingStrategy.PROB  # Default to probabilistic
        if node_idx in self.config.node_routing_strategies:
            class_strategies = self.config.node_routing_strategies[node_idx]
            if class_id in class_strategies:
                strategy = class_strategies[class_id]

        # Dispatch to appropriate routing strategy
        if strategy == RoutingStrategy.RAND:
            return self._select_random_destination(destinations)
        elif strategy == RoutingStrategy.RROBIN:
            return self._select_round_robin_destination(node_idx, destinations)
        elif strategy == RoutingStrategy.WRROBIN:
            return self._select_weighted_round_robin_destination(node_idx, destinations, probs)
        elif strategy == RoutingStrategy.JSQ:
            return self._select_jsq_destination(destinations)
        elif strategy == RoutingStrategy.KCHOICES:
            return self._select_kchoices_destination(destinations)
        else:
            # PROB (probabilistic) - default behavior
            return self._select_probabilistic_destination(destinations, probs)

    def _select_probabilistic_destination(
        self, destinations: List[Tuple[int, int]], probs: List[float]
    ) -> Tuple[int, int]:
        """
        Select destination using probabilistic routing (PROB strategy).

        Uses the routing matrix probabilities to sample a destination.
        """
        probs_arr = np.array(probs)
        probs_arr = probs_arr / probs_arr.sum()
        choice = self.rng.choice(len(destinations), p=probs_arr)
        return destinations[choice]

    def _select_random_destination(
        self, destinations: List[Tuple[int, int]]
    ) -> Tuple[int, int]:
        """
        Select destination using random uniform routing (RAND strategy).

        Each destination has equal probability of selection.
        """
        idx = int(self.rng.random() * len(destinations))
        idx = min(idx, len(destinations) - 1)  # Ensure valid index
        return destinations[idx]

    def _select_round_robin_destination(
        self, node_idx: int, destinations: List[Tuple[int, int]]
    ) -> Tuple[int, int]:
        """
        Select destination using round-robin routing (RROBIN strategy).

        Cycles through destinations in order.
        """
        counter = self.round_robin_counters.get(node_idx, 0)
        idx = counter % len(destinations)
        self.round_robin_counters[node_idx] = counter + 1
        return destinations[idx]

    def _select_weighted_round_robin_destination(
        self, node_idx: int, destinations: List[Tuple[int, int]], weights: List[float]
    ) -> Tuple[int, int]:
        """
        Select destination using weighted round-robin routing (WRROBIN strategy).

        Uses routing probabilities as weights for the round-robin schedule.
        The counter cycles through a virtual sequence where each destination
        appears proportionally to its weight.
        """
        # Normalize weights to integers (multiply by scale factor)
        scale = 100
        int_weights = [max(1, int(w * scale)) for w in weights]
        total_weight = sum(int_weights)

        counter = self.round_robin_counters.get(node_idx, 0)
        position = counter % total_weight
        self.round_robin_counters[node_idx] = counter + 1

        # Find destination corresponding to this position
        cumulative = 0
        for i, dest in enumerate(destinations):
            cumulative += int_weights[i]
            if position < cumulative:
                return dest

        return destinations[-1]

    def _select_jsq_destination(
        self, destinations: List[Tuple[int, int]]
    ) -> Tuple[int, int]:
        """
        Select destination using Join Shortest Queue routing (JSQ strategy).

        Routes to the destination node with the fewest jobs (queue + in service).
        """
        min_queue = float('inf')
        min_dest = destinations[0]

        for dest_node, dest_class in destinations:
            # Check if destination is a service node (Queue or Delay)
            if dest_node < 0 or dest_node >= self.config.num_nodes:
                continue

            # Get station index for this node
            station_idx = -1
            if self.config.node_to_station is not None and dest_node < len(self.config.node_to_station):
                station_idx = int(self.config.node_to_station[dest_node])

            if station_idx < 0:
                continue

            # Find service node index for this station
            svc_idx = -1
            for idx, st_idx in enumerate(self.config.service_stations):
                if st_idx == station_idx:
                    svc_idx = idx
                    break

            if svc_idx >= 0 and svc_idx in self.stats:
                # Sum queue lengths across all classes at this station
                total_queue = sum(self.stats[svc_idx].current_queue_length)
                if total_queue < min_queue:
                    min_queue = total_queue
                    min_dest = (dest_node, dest_class)

        return min_dest

    def _select_kchoices_destination(
        self, destinations: List[Tuple[int, int]]
    ) -> Tuple[int, int]:
        """
        Select destination using Power of K Choices routing (KCHOICES strategy).

        Randomly selects K candidates and routes to the one with shortest queue.
        """
        k = self.kchoices_k

        if len(destinations) <= k:
            # If K >= number of destinations, use JSQ
            return self._select_jsq_destination(destinations)

        # Randomly select K candidates
        indices = self.rng.choice(len(destinations), size=k, replace=False)
        candidates = [destinations[i] for i in indices]

        # Select the one with shortest queue among candidates
        return self._select_jsq_destination(candidates)

    def _process_cache_access(self, node_idx: int, current_class: int) -> int:
        """
        Process a cache access for a job, returning the new class after hit/miss determination.

        This method:
        1. Samples which item is being accessed (using access probabilities or uniform)
        2. Checks if item is in cache (hit) or not (miss)
        3. Updates cache state based on replacement policy
        4. Returns appropriate class (hit_class or miss_class)

        Args:
            node_idx: Cache node index
            current_class: Current job class

        Returns:
            New class after cache access (hit_class or miss_class)
        """
        cache_config = self.config.cache_config.get(node_idx)
        if cache_config is None:
            return current_class

        cache_state = self.cache_states.get(node_idx)
        if cache_state is None:
            return current_class

        # Initialize hit/miss counters if needed
        if current_class not in cache_state.total_hits:
            cache_state.total_hits[current_class] = 0
        if current_class not in cache_state.total_misses:
            cache_state.total_misses[current_class] = 0

        # Sample which item is being accessed
        if cache_config.access_probs is not None and len(cache_config.access_probs) > 0:
            # Use specified access probabilities
            item_idx = int(self.rng.choice(len(cache_config.access_probs), p=cache_config.access_probs))
        else:
            # Uniform random access
            item_idx = int(self.rng.integers(0, cache_config.num_items))

        # Check if item is in cache
        if item_idx in cache_state.cache_contents:
            # Cache hit
            cache_state.total_hits[current_class] += 1
            self._handle_cache_hit(cache_config, cache_state, item_idx)
            new_class = cache_config.hit_class.get(current_class, current_class)
        else:
            # Cache miss
            cache_state.total_misses[current_class] += 1
            self._handle_cache_miss(cache_config, cache_state, item_idx)
            new_class = cache_config.miss_class.get(current_class, current_class)

        return new_class

    def _handle_cache_hit(self, cache_config: CacheConfig, cache_state: CacheState, item_idx: int) -> None:
        """
        Handle a cache hit - move item according to replacement policy.

        Args:
            cache_config: Cache configuration
            cache_state: Cache state
            item_idx: Item that was accessed (hit)
        """
        policy = cache_config.replacement_policy

        if policy == ReplacementPolicy.LRU:
            # Move to front (most recently used position)
            cache_state.cache_contents.remove(item_idx)
            cache_state.cache_contents.insert(0, item_idx)
        elif policy == ReplacementPolicy.FIFO:
            # FIFO: don't move on hit, item stays where it is
            pass
        elif policy == ReplacementPolicy.RR:
            # Random replacement: don't move on hit
            pass

    def _handle_cache_miss(self, cache_config: CacheConfig, cache_state: CacheState, item_idx: int) -> None:
        """
        Handle a cache miss - insert item according to replacement policy.

        Args:
            cache_config: Cache configuration
            cache_state: Cache state
            item_idx: Item that was accessed (miss, now needs to be inserted)
        """
        policy = cache_config.replacement_policy
        capacity = cache_config.capacity

        if policy == ReplacementPolicy.LRU:
            # Insert at front (most recently used position)
            cache_state.cache_contents.insert(0, item_idx)
            # Evict if over capacity (LRU: remove from end)
            while len(cache_state.cache_contents) > capacity:
                cache_state.cache_contents.pop()

        elif policy == ReplacementPolicy.FIFO:
            # Insert at end (FIFO order)
            cache_state.cache_contents.append(item_idx)
            # Evict oldest (FIFO: remove from front)
            while len(cache_state.cache_contents) > capacity:
                cache_state.cache_contents.pop(0)

        elif policy == ReplacementPolicy.RR:
            # Insert at front
            cache_state.cache_contents.insert(0, item_idx)
            # Random eviction if over capacity
            while len(cache_state.cache_contents) > capacity:
                evict_idx = int(self.rng.integers(0, len(cache_state.cache_contents)))
                cache_state.cache_contents.pop(evict_idx)

    def get_result(self) -> DESResult:
        """Compile and return simulation results."""
        result = DESResult()

        num_stations = self.config.num_stations
        num_classes = self.config.num_classes
        num_service = len(self.config.service_nodes)
        sim_time = self.env.now

        # Determine truncation index for MSER observations
        truncation_obs_idx = 0
        if self.mser_enabled and self.warmup_end_time > 0:
            # Find observation index corresponding to warmup_end_time
            for i, obs_time in enumerate(self.mser_observation_times):
                if obs_time >= self.warmup_end_time:
                    truncation_obs_idx = i
                    break

        # Compute effective simulation time (post-warmup)
        effective_sim_time = sim_time - self.warmup_end_time if self.warmup_end_time > 0 else sim_time

        # Initialize result matrices
        result.QN = np.zeros((num_stations, num_classes))
        result.UN = np.zeros((num_stations, num_classes))
        result.RN = np.zeros((num_stations, num_classes))
        result.TN = np.zeros((num_stations, num_classes))
        result.AN = np.zeros((num_stations, num_classes))

        # Fill in metrics for service nodes
        for svc_idx in range(num_service):
            station_idx = self.config.service_stations[svc_idx]
            num_servers = self.config.num_servers.get(svc_idx, 1)
            is_delay = self.config.is_delay_node[svc_idx] if svc_idx < len(self.config.is_delay_node) else False
            stats = self.stats[svc_idx]

            # Finalize time-weighted stats
            for k in range(num_classes):
                stats.update_queue(k, sim_time)
                stats.update_busy(k, sim_time)

            for k in range(num_classes):
                # For queue length, use MSER observations if available (more accurate post-warmup)
                if (self.mser_enabled and
                    svc_idx in self.mser_observations and
                    k in self.mser_observations[svc_idx] and
                    len(self.mser_observations[svc_idx][k]) > truncation_obs_idx + 5):
                    # Use average of post-warmup MSER observations
                    post_warmup_obs = self.mser_observations[svc_idx][k][truncation_obs_idx:]
                    qlen = np.mean(post_warmup_obs) if post_warmup_obs else 0.0
                else:
                    qlen = stats.get_avg_queue_length(k, sim_time)

                result.QN[station_idx, k] = qlen

                # For Delay (infinite server) nodes, utilization = queue length
                # since all customers are always in service
                if is_delay:
                    result.UN[station_idx, k] = qlen
                else:
                    result.UN[station_idx, k] = stats.get_utilization(k, sim_time, min(num_servers, 1000000))
                result.RN[station_idx, k] = stats.get_avg_response_time(k)
                result.TN[station_idx, k] = stats.get_throughput(k, sim_time)
                result.AN[station_idx, k] = stats.arrived_customers[k] / sim_time if sim_time > 0 else 0.0

        # Fill in metrics for SPN places (if any)
        for place_idx, node_idx in enumerate(self.config.place_nodes):
            # Map node_idx to station_idx (places may have different station mapping)
            station_idx = None
            if self.config.node_to_station is not None:
                try:
                    n2s = self.config.node_to_station
                    if n2s.ndim == 1 and node_idx < len(n2s):
                        station_idx = int(n2s[node_idx])
                    elif n2s.ndim == 2 and node_idx < n2s.shape[0]:
                        station_idx = int(n2s[node_idx, 0])
                    else:
                        station_idx = node_idx
                except (IndexError, TypeError):
                    station_idx = node_idx
            else:
                # Use node_idx as station_idx if no mapping
                station_idx = node_idx

            if station_idx is not None and station_idx < num_stations:
                # Finalize place statistics
                current_time = sim_time
                elapsed = current_time - self.last_place_update.get(place_idx, 0.0)
                if elapsed > 0:
                    tokens = self.place_tokens.get(place_idx, np.zeros(num_classes))
                    self.place_token_time[place_idx] = self.place_token_time.get(
                        place_idx, np.zeros(num_classes)
                    ) + tokens * elapsed
                    self.last_place_update[place_idx] = current_time

                # Finalize in-service time for all transitions
                for trans_idx in range(len(self.config.transition_nodes)):
                    if trans_idx in self.transition_in_service_count:
                        if hasattr(self, 'last_transition_update'):
                            trans_elapsed = current_time - self.last_transition_update.get(trans_idx, 0.0)
                            if trans_elapsed > 0:
                                self.transition_in_service_time[trans_idx] += (
                                    self.transition_in_service_count[trans_idx] * trans_elapsed
                                )
                                self.last_transition_update[trans_idx] = current_time

                for k in range(num_classes):
                    # Queue length = average tokens in place only
                    # For SPNs, QLen at a place is tokens AT the place, not in-service tokens
                    # This matches MATLAB SSA behavior
                    place_qlen = 0.0
                    if effective_sim_time > 0:
                        place_qlen = self.place_token_time.get(
                            place_idx, np.zeros(num_classes)
                        )[k] / effective_sim_time

                    result.QN[station_idx, k] = place_qlen

                    # Throughput = firing completions per unit time
                    # For SPNs, count token completions (firings that consumed/produced tokens)
                    if sim_time > 0:
                        result.TN[station_idx, k] = self.place_completions.get(
                            place_idx, np.zeros(num_classes, dtype=int)
                        )[k] / sim_time

        # System throughput
        result.XN = np.zeros((1, num_classes))
        for k in range(num_classes):
            result.XN[0, k] = self.system_completed[k] / sim_time if sim_time > 0 else 0.0

        # Compute confidence intervals if enabled
        if self.options.confint > 0 and self.mser_enabled:
            self._compute_confidence_intervals(
                result, num_stations, num_classes, num_service, truncation_obs_idx
            )

        # Populate transient results if in transient mode
        if self.is_transient_mode and self.transient_times:
            self._populate_transient_results(result, num_stations, num_classes, num_service)

        # Collect response time samples for CDF computation
        result.resp_time_samples = {}
        result.pass_time_samples = {}
        for svc_idx in range(num_service):
            station_idx = self.config.service_stations[svc_idx]
            result.resp_time_samples[station_idx] = {}
            result.pass_time_samples[station_idx] = {}
            for k in range(num_classes):
                # Get samples collected during simulation
                samples = self.stats[svc_idx].response_time_samples.get(k, [])
                result.resp_time_samples[station_idx][k] = samples
                # For now, passage time samples are the same as response time samples
                result.pass_time_samples[station_idx][k] = samples

        # Metadata
        result.total_events = self.total_event_count
        result.simulation_time = sim_time
        result.warmup_end_time = self.warmup_end_time

        return result

    def _populate_transient_results(
        self, result: DESResult, num_stations: int, num_classes: int, num_service: int
    ) -> None:
        """
        Populate transient analysis results (QNt, UNt, TNt, t).

        Creates time-indexed arrays of performance metrics collected during simulation.

        Args:
            result: DESResult object to populate
            num_stations: Number of stations
            num_classes: Number of job classes
            num_service: Number of service nodes
        """
        num_time_points = len(self.transient_times)

        # Time vector
        result.t = np.array(self.transient_times)

        # Initialize transient matrices [time_points x stations x classes]
        result.QNt = np.zeros((num_time_points, num_stations, num_classes))
        result.UNt = np.zeros((num_time_points, num_stations, num_classes))
        result.TNt = np.zeros((num_time_points, num_stations, num_classes))

        # Fill in transient data for service nodes
        for svc_idx in range(num_service):
            station_idx = self.config.service_stations[svc_idx]

            for k in range(num_classes):
                if svc_idx in self.transient_queue_lengths and k in self.transient_queue_lengths[svc_idx]:
                    q_data = self.transient_queue_lengths[svc_idx][k]
                    u_data = self.transient_utilizations[svc_idx][k]
                    t_data = self.transient_throughputs[svc_idx][k]

                    # Ensure data length matches time points
                    data_len = min(len(q_data), num_time_points)
                    for t_idx in range(data_len):
                        result.QNt[t_idx, station_idx, k] = q_data[t_idx]
                        result.UNt[t_idx, station_idx, k] = u_data[t_idx]
                        result.TNt[t_idx, station_idx, k] = t_data[t_idx]

    def _compute_confidence_intervals(
        self, result: DESResult, num_stations: int, num_classes: int,
        num_service: int, truncation_obs_idx: int
    ) -> None:
        """
        Compute confidence intervals for performance metrics using Overlapping Batch Means (OBM).

        Populates result.QNCI, result.UNCI, result.RNCI, result.TNCI with
        half-widths of confidence intervals.

        Args:
            result: DESResult object to populate
            num_stations: Number of stations
            num_classes: Number of job classes
            num_service: Number of service nodes
            truncation_obs_idx: Index where warmup ends in observations
        """
        from scipy import stats as scipy_stats

        confint_level = self.options.confint
        alpha = 1.0 - confint_level

        # Initialize CI matrices
        result.QNCI = np.zeros((num_stations, num_classes))
        result.UNCI = np.zeros((num_stations, num_classes))
        result.RNCI = np.zeros((num_stations, num_classes))
        result.TNCI = np.zeros((num_stations, num_classes))

        # Minimum observations needed for valid CI
        min_obs = self.options.ciminobs
        min_batches = self.options.ciminbatch

        for svc_idx in range(num_service):
            station_idx = self.config.service_stations[svc_idx]

            for k in range(num_classes):
                # Queue length CI from MSER observations
                if (svc_idx in self.mser_observations and
                    k in self.mser_observations[svc_idx]):
                    post_warmup = self.mser_observations[svc_idx][k][truncation_obs_idx:]
                    if len(post_warmup) >= min_obs:
                        ci_hw = self._compute_obm_ci(post_warmup, alpha, min_batches)
                        result.QNCI[station_idx, k] = ci_hw

                # For other metrics, compute CI from stats if we have response time samples
                stats = self.stats[svc_idx]
                if hasattr(stats, 'response_times') and k in stats.response_times:
                    resp_times = stats.response_times[k]
                    if len(resp_times) >= min_obs:
                        ci_hw = self._compute_obm_ci(resp_times, alpha, min_batches)
                        result.RNCI[station_idx, k] = ci_hw

    def _compute_obm_ci(
        self, observations: List[float], alpha: float, min_batches: int
    ) -> float:
        """
        Compute confidence interval half-width using Overlapping Batch Means.

        Args:
            observations: Post-warmup observations
            alpha: Significance level (e.g., 0.05 for 95% CI)
            min_batches: Minimum number of batches

        Returns:
            Half-width of confidence interval
        """
        from scipy import stats as scipy_stats

        n = len(observations)
        if n < min_batches:
            return 0.0

        # Batch size for OBM (with 50% overlap)
        overlap = self.options.obmoverlap
        batch_size = max(n // (min_batches * 2), 1)  # Ensure reasonable batch size
        step = max(1, int(batch_size * (1 - overlap)))  # Step with overlap

        # Compute batch means
        batch_means = []
        i = 0
        while i + batch_size <= n:
            batch = observations[i:i + batch_size]
            batch_means.append(np.mean(batch))
            i += step

        if len(batch_means) < 2:
            return 0.0

        # Compute sample variance of batch means (adjusted for overlap)
        batch_mean_var = np.var(batch_means, ddof=1)

        # Adjust variance for overlapping batches
        # The adjustment factor depends on the overlap fraction
        # For 50% overlap, the correction factor is approximately 4/3
        overlap_correction = 1.0 / (1.0 - overlap) if overlap < 1.0 else 1.0
        adjusted_var = batch_mean_var * overlap_correction

        # Degrees of freedom
        dof = len(batch_means) - 1

        # t critical value
        t_crit = scipy_stats.t.ppf(1 - alpha / 2, dof)

        # Standard error
        se = np.sqrt(adjusted_var / len(batch_means))

        # Half-width
        return t_crit * se
