"""
Node implementations for LINE native Python.

This module provides concrete node types: Queue, Delay, Source, Sink, Fork, Join.
Ported from MATLAB implementation in matlab/src/lang/nodes/
"""

from typing import Optional, Dict, Union
from enum import IntEnum
import numpy as np

from .base import (
    Node, Station, StatefulNode, JobClass, NodeType, SchedStrategy,
    SchedStrategyType, JoinStrategy, DropStrategy, ReplacementStrategy,
    RoutingStrategy
)
from ..distributions.continuous import Immediate


class Queue(Station):
    """
    Queue station for customer service.

    A Queue represents a queueing station with one or more servers
    and configurable scheduling strategy.
    """

    def __init__(self, model, name: str, sched_strategy: SchedStrategy = SchedStrategy.FCFS):
        """
        Initialize a Queue node.

        Args:
            model: Network instance
            name: Queue name
            sched_strategy: Scheduling strategy (default: FCFS)

        Raises:
            ValueError: If scheduling strategy is invalid
        """
        super().__init__(NodeType.QUEUE, name)
        self.set_model(model)
        model.add_node(self)

        self._sched_strategy = sched_strategy
        self._sched_param = {}  # Dict[JobClass, scheduling_parameter]
        self._service_process = {}  # Dict[JobClass, Distribution]
        self._server_types = None  # For heterogeneous servers

        # Initialize scheduling policy type
        self._set_sched_policy()

    def _set_sched_policy(self) -> None:
        """Determine if scheduling strategy is preemptive or non-preemptive."""
        preemptive_strategies = {
            SchedStrategy.LCFSPR,
            SchedStrategy.FCFSPR,
            SchedStrategy.EDF
        }
        self._sched_policy = (
            SchedStrategyType.PR if self._sched_strategy in preemptive_strategies
            else SchedStrategyType.NP
        )

    def get_sched_strategy(self) -> SchedStrategy:
        """Get scheduling strategy."""
        return self._sched_strategy

    def set_sched_strategy(self, strategy: SchedStrategy) -> None:
        """
        Set scheduling strategy.

        Args:
            strategy: Scheduling strategy (FCFS, LCFS, PS, etc.)
        """
        self._sched_strategy = strategy
        self._set_sched_policy()
        self._invalidate_java()

    def get_sched_policy(self) -> SchedStrategyType:
        """Get scheduling policy (preemptive or non-preemptive)."""
        return self._sched_policy

    def set_service(self, jobclass: JobClass, distribution, weight: float = None) -> None:
        """
        Set service time distribution for a job class.

        Args:
            jobclass: Job class
            distribution: Service distribution (Exp, Erlang, HyperExp, etc.)
                         or Workflow for activity-based distributions
            weight: Optional weight/scale factor for the distribution
        """
        # Store the distribution (and optional weight for future use)
        self._service_process[jobclass] = distribution
        if weight is not None:
            # Store weight metadata for analysis
            if not hasattr(self, '_service_weights'):
                self._service_weights = {}
            self._service_weights[jobclass] = weight
            # Also set scheduling parameter for DPS/GPS scheduling
            self._sched_param[jobclass] = weight
        self._invalidate_java()
        # Invalidate model's cached struct since service distribution changed
        if self._model is not None:
            self._model._sn = None

    def get_service(self, jobclass: JobClass):
        """
        Get service distribution for a job class.

        Args:
            jobclass: Job class

        Returns:
            Service distribution or None if not set
        """
        return self._service_process.get(jobclass)

    def set_strategy_param(self, jobclass: JobClass, param) -> None:
        """
        Set scheduling parameter for a job class.

        For PS/DPS/GPS: weight or priority parameter
        For other strategies: may be unused

        Args:
            jobclass: Job class
            param: Scheduling parameter (weight, priority, etc.)
        """
        self._sched_param[jobclass] = param
        self._invalidate_java()

    def get_strategy_param(self, jobclass: JobClass):
        """Get scheduling parameter for a job class."""
        return self._sched_param.get(jobclass, 1.0)

    def is_service_defined(self, jobclass: JobClass) -> bool:
        """Check if service distribution is defined for a class."""
        return jobclass in self._service_process

    def is_service_disabled(self, jobclass: JobClass) -> bool:
        """Check if service is explicitly disabled for a class."""
        return self._service_process.get(jobclass) is None

    def get_service_rates(self) -> Dict[JobClass, float]:
        """
        Get mean service rates for all classes.

        Returns:
            Dict mapping classes to mean service rates
        """
        rates = {}
        for jobclass, dist in self._service_process.items():
            if dist is not None and hasattr(dist, 'getMean'):
                rates[jobclass] = dist.getMean()
        return rates

    # PascalCase aliases for MATLAB compatibility
    getSchedStrategy = get_sched_strategy
    setSchedStrategy = set_sched_strategy
    getSchedPolicy = get_sched_policy
    setService = set_service
    getService = get_service
    setStrategyParam = set_strategy_param
    getStrategyParam = get_strategy_param
    getServiceRates = get_service_rates

    # =========================================================================
    # HETEROGENEOUS SERVER METHODS
    # =========================================================================

    def is_heterogeneous(self) -> bool:
        """
        Check if this queue has heterogeneous servers.

        Returns:
            True if server types have been defined
        """
        return self._server_types is not None and len(self._server_types) > 0

    def add_server_type(self, server_type) -> None:
        """
        Add a server type to this queue.

        Args:
            server_type: ServerType object to add

        Raises:
            ValueError: If server type is already in queue
        """
        if self._server_types is None:
            self._server_types = []

        if server_type in self._server_types:
            raise ValueError(f"Server type '{server_type.name}' already in queue")

        server_type.set_id(len(self._server_types))
        server_type.set_parent_queue(self)
        self._server_types.append(server_type)
        self._invalidate_java()

    def get_server_types(self):
        """Get the list of server types."""
        return self._server_types if self._server_types else []

    def set_hetero_sched_policy(self, policy) -> None:
        """
        Set the heterogeneous scheduling policy.

        Args:
            policy: HeteroSchedPolicy value
        """
        self._hetero_sched_policy = policy
        self._invalidate_java()

    def get_hetero_sched_policy(self):
        """Get the heterogeneous scheduling policy."""
        return getattr(self, '_hetero_sched_policy', None)

    def set_hetero_service(self, jobclass, server_type, distribution) -> None:
        """
        Set service distribution for a specific job class and server type.

        Args:
            jobclass: Job class
            server_type: ServerType object
            distribution: Service distribution
        """
        if not hasattr(self, '_hetero_service'):
            self._hetero_service = {}

        key = (jobclass, server_type)
        self._hetero_service[key] = distribution
        self._invalidate_java()

    def get_hetero_service(self, jobclass, server_type):
        """
        Get service distribution for a specific job class and server type.

        Args:
            jobclass: Job class
            server_type: ServerType object

        Returns:
            Service distribution or None
        """
        if not hasattr(self, '_hetero_service'):
            return None
        return self._hetero_service.get((jobclass, server_type))

    def get_total_num_of_servers(self) -> int:
        """
        Get total number of servers across all server types.

        Returns:
            Total server count
        """
        if not self._server_types:
            return self.number_of_servers
        return sum(st.get_num_of_servers() for st in self._server_types)

    # PascalCase aliases for heterogeneous methods
    isHeterogeneous = is_heterogeneous
    addServerType = add_server_type
    getServerTypes = get_server_types
    setHeteroSchedPolicy = set_hetero_sched_policy
    getHeteroSchedPolicy = get_hetero_sched_policy
    setHeteroService = set_hetero_service
    getHeteroService = get_hetero_service
    getTotalNumOfServers = get_total_num_of_servers

    # =========================================================================
    # POLLING AND SWITCHOVER METHODS
    # =========================================================================

    def set_polling_type(self, polling_type, k=None) -> None:
        """
        Set the polling type for this queue.

        Args:
            polling_type: PollingType enum value (EXHAUSTIVE, GATED, KLIMITED)
            k: For KLIMITED polling, the maximum number of jobs to serve (default: 1)
        """
        self._polling_type = polling_type
        self._polling_k = k if k is not None else 1
        self._invalidate_java()

        # Set default Immediate switchover for all classes (matching MATLAB behavior)
        # This is required for JMT to properly simulate polling systems
        model = self.get_model()
        if model is not None:
            classes = model.get_classes()
            for jobclass in classes:
                self.set_switchover(jobclass, Immediate())

    def get_polling_type(self):
        """Get the polling type for this queue."""
        return getattr(self, '_polling_type', None)

    def set_switchover(self, from_class, to_class_or_dist, distribution=None) -> None:
        """
        Set switchover time between job classes.

        Can be called as:
        - set_switchover(from_class, to_class, distribution) - for class-to-class switchover
        - set_switchover(jobclass, distribution) - for single class switchover time

        Args:
            from_class: Source job class
            to_class_or_dist: Target job class, OR distribution if only 2 args
            distribution: Switchover time distribution (optional if 2 args)
        """
        if not hasattr(self, '_switchover'):
            self._switchover = {}

        if distribution is None:
            # Called as set_switchover(class, distribution)
            self._switchover[from_class] = to_class_or_dist
        else:
            # Called as set_switchover(from_class, to_class, distribution)
            self._switchover[(from_class, to_class_or_dist)] = distribution
        self._invalidate_java()

    def get_switchover(self, from_class, to_class):
        """Get switchover time distribution between classes."""
        if not hasattr(self, '_switchover'):
            return None
        return self._switchover.get((from_class, to_class))

    # =========================================================================
    # ROUTING METHODS
    # =========================================================================

    def set_prob_routing(self, jobclass, destination, prob: float) -> None:
        """
        Set probabilistic routing to a destination node.

        Args:
            jobclass: Job class
            destination: Destination node
            prob: Routing probability (0 to 1)
        """
        if not hasattr(self, '_prob_routing'):
            self._prob_routing = {}
        if jobclass not in self._prob_routing:
            self._prob_routing[jobclass] = {}
        self._prob_routing[jobclass][destination] = prob
        self._invalidate_java()

    def get_prob_routing(self, jobclass):
        """Get probabilistic routing for a job class."""
        if not hasattr(self, '_prob_routing'):
            return {}
        return self._prob_routing.get(jobclass, {})

    def set_state(self, state) -> None:
        """
        Set initial state for this node.

        Args:
            state: Array of initial job counts per class
        """
        self._state = np.asarray(state)
        self._state_explicitly_set = True
        self._invalidate_java()
        # Invalidate the model struct so _refresh_state() runs again
        model = self.get_model() if hasattr(self, 'get_model') else None
        if model is not None and hasattr(model, '_has_struct'):
            model._has_struct = False
            model._sn = None

    def get_state(self):
        """Get the current state of this node."""
        return getattr(self, '_state', None)

    def set_routing_weight(self, jobclass: JobClass, destination, weight: float) -> None:
        """
        Set routing weight for weighted round-robin routing.

        Args:
            jobclass: Job class
            destination: Destination node
            weight: Routing weight for this destination
        """
        if not hasattr(self, '_routing_weights'):
            self._routing_weights = {}
        if jobclass not in self._routing_weights:
            self._routing_weights[jobclass] = {}
        self._routing_weights[jobclass][destination] = weight
        self._invalidate_java()

    def get_routing_weight(self, jobclass: JobClass, destination) -> float:
        """Get routing weight for a job class and destination."""
        if not hasattr(self, '_routing_weights'):
            return 1.0
        if jobclass not in self._routing_weights:
            return 1.0
        return self._routing_weights[jobclass].get(destination, 1.0)

    # PascalCase aliases
    setPollingType = set_polling_type
    getPollingType = get_polling_type
    setSwitchover = set_switchover
    getSwitchover = get_switchover
    setProbRouting = set_prob_routing
    getProbRouting = get_prob_routing
    setState = set_state
    getState = get_state
    setRoutingWeight = set_routing_weight
    getRoutingWeight = get_routing_weight


class Delay(Queue):
    """
    Delay node (infinite server).

    A Delay represents an infinite-capacity queue (think time station in closed models).
    All jobs immediately begin service without queueing.
    """

    def __init__(self, model, name: str):
        """
        Initialize a Delay node.

        Args:
            model: Network instance
            name: Delay node name
        """
        super().__init__(model, name, SchedStrategy.INF)
        # Override node_type to DELAY (Queue sets it to QUEUE)
        self._node_type = NodeType.DELAY
        self._number_of_servers = np.inf

    def set_number_of_servers(self, value: int) -> None:
        """
        Set number of servers (always infinity for Delay).

        Args:
            value: Ignored (always set to infinity)

        Raises:
            ValueError: If not infinity
        """
        if not np.isinf(value):
            raise ValueError(f"[{self.name}] Delay node must have infinite servers")
        self._number_of_servers = np.inf
        self._invalidate_java()


class Source(Station):
    """
    Source node for external arrivals.

    A Source represents the entry point for open-class jobs.
    Each network can have at most one Source node.
    """

    def __init__(self, model, name: str):
        """
        Initialize a Source node.

        Args:
            model: Network instance
            name: Source node name
        """
        super().__init__(NodeType.SOURCE, name)
        self.set_model(model)
        model.add_node(self)

        self._arrival_process = {}  # Dict[JobClass, Distribution]
        self._number_of_servers = np.inf

    def set_arrival(self, jobclass: JobClass, distribution) -> None:
        """
        Set arrival distribution for a job class.

        Args:
            jobclass: Job class
            distribution: Arrival distribution (Exp, Erlang, etc.)
        """
        self._arrival_process[jobclass] = distribution
        self._invalidate_java()

    def get_arrival(self, jobclass: JobClass):
        """
        Get arrival distribution for a job class.

        Args:
            jobclass: Job class

        Returns:
            Arrival distribution or None if not set
        """
        return self._arrival_process.get(jobclass)

    def get_arrival_rates(self) -> Dict[JobClass, float]:
        """
        Get arrival rates (1/mean) for all classes.

        Returns:
            Dict mapping classes to arrival rates
        """
        rates = {}
        for jobclass, dist in self._arrival_process.items():
            if dist is not None and hasattr(dist, 'getMean'):
                rates[jobclass] = 1.0 / dist.getMean()
        return rates

    # PascalCase aliases for MATLAB compatibility
    setArrival = set_arrival
    getArrival = get_arrival
    getArrivalRates = get_arrival_rates


class Sink(Node):
    """
    Sink node for job departure.

    A Sink is the exit point for open-class jobs.
    Each network can have at most one Sink node.
    """

    def __init__(self, model, name: str):
        """
        Initialize a Sink node.

        Args:
            model: Network instance
            name: Sink node name
        """
        super().__init__(NodeType.SINK, name)
        self.set_model(model)
        model.add_node(self)


class Fork(Node):
    """
    Fork node for parallel processing.

    A Fork splits an incoming job into multiple parallel tasks
    that must be synchronized at a corresponding Join node.
    """

    def __init__(self, model, name: str):
        """
        Initialize a Fork node.

        Args:
            model: Network instance
            name: Fork node name
        """
        super().__init__(NodeType.FORK, name)
        self.set_model(model)
        model.add_node(self)

        self._tasks_per_link = None
        self._capacity = np.inf

    def set_tasks_per_link(self, ntasks: np.ndarray) -> None:
        """
        Set number of tasks spawned per outgoing link (experimental).

        Args:
            ntasks: Array with task counts for each link
        """
        self._tasks_per_link = np.array(ntasks)
        self._invalidate_java()

    def get_tasks_per_link(self) -> Optional[np.ndarray]:
        """Get tasks per link."""
        return self._tasks_per_link

    @property
    def capacity(self) -> float:
        """Get capacity."""
        return self._capacity

    @capacity.setter
    def capacity(self, value: float) -> None:
        """Set capacity."""
        self._capacity = value
        self._invalidate_java()


class Join(Station):
    """
    Join node for parallel processing synchronization.

    A Join waits for all parallel tasks from a corresponding Fork
    before releasing the job downstream.
    """

    def __init__(self, model, name: str, fork: Optional[Fork] = None):
        """
        Initialize a Join node.

        Args:
            model: Network instance
            name: Join node name
            fork: Associated Fork node (optional)
        """
        super().__init__(NodeType.JOIN, name)
        self.set_model(model)
        model.add_node(self)

        self._fork = fork
        self._join_strategy = {}  # Dict[JobClass, JoinStrategy]
        self._required = {}  # Dict[JobClass, required_count] (-1 = all)
        self._number_of_servers = np.inf

    def set_fork(self, fork: Fork) -> None:
        """
        Set associated Fork node.

        Args:
            fork: Fork node to pair with
        """
        self._fork = fork
        self._invalidate_java()

    def get_fork(self) -> Optional[Fork]:
        """Get associated Fork node."""
        return self._fork

    def set_strategy(self, jobclass: JobClass, strategy: JoinStrategy) -> None:
        """
        Set join strategy for a job class.

        Args:
            jobclass: Job class
            strategy: Join strategy (STD, QUORUM, etc.)
        """
        self._join_strategy[jobclass] = strategy
        self._invalidate_java()

    def get_strategy(self, jobclass: JobClass) -> JoinStrategy:
        """Get join strategy for a job class."""
        return self._join_strategy.get(jobclass, JoinStrategy.STD)

    def set_required(self, jobclass: JobClass, nrequired: int) -> None:
        """
        Set required number of tasks to wait for (quorum join).

        Args:
            jobclass: Job class
            nrequired: Number of required tasks (-1 = all)
        """
        self._required[jobclass] = nrequired
        self._invalidate_java()

    def get_required(self, jobclass: JobClass) -> int:
        """Get required task count."""
        return self._required.get(jobclass, -1)


class Router(Node):
    """
    Router node for routing decisions.

    A Router is a node that routes jobs without service delay.
    Can be used for probabilistic routing and routing-dependent decisions.
    """

    def __init__(self, model, name: str):
        """
        Initialize a Router node.

        Args:
            model: Network instance
            name: Router node name
        """
        super().__init__(NodeType.ROUTER, name)
        self.set_model(model)
        model.add_node(self)


class ClassSwitch(Node):
    """
    ClassSwitch node for job class switching.

    A ClassSwitch allows jobs to change class without service delay.
    Useful for modeling class-dependent routing and scheduling.
    """

    def __init__(self, model, name: str, cs_matrix=None):
        """
        Initialize a ClassSwitch node.

        Args:
            model: Network instance
            name: ClassSwitch node name
            cs_matrix: Optional K×K class switching probability matrix where
                      element (i,j) is the probability that a job in class i
                      switches to class j. If not provided, must be set later
                      using set_class_switching_matrix().
        """
        super().__init__(NodeType.CLASSSWITCH, name)
        self.set_model(model)
        model.add_node(self)

        # Set switching matrix if provided
        if cs_matrix is not None:
            self._switch_matrix = np.asarray(cs_matrix)
        else:
            self._switch_matrix = None  # Class switching probabilities

    def init_class_switch_matrix(self) -> np.ndarray:
        """
        Initialize and return a class switching matrix.

        Creates a K×K matrix of zeros where K is the number of classes.
        Use this to create a template that can be filled with switching
        probabilities before calling set_class_switching_matrix().

        Returns:
            np.ndarray: K×K matrix initialized to zeros

        Example:
            >>> csmatrix = cs_node.init_class_switch_matrix()
            >>> csmatrix[0, 1] = 0.3  # 30% switch from class 0 to class 1
            >>> csmatrix[0, 0] = 0.7  # 70% stay in class 0
            >>> csmatrix[1, 0] = 1.0  # 100% switch from class 1 to class 0
            >>> cs_node.set_class_switching_matrix(csmatrix)
        """
        K = self._model.get_number_of_classes()
        return np.zeros((K, K))

    initClassSwitchMatrix = init_class_switch_matrix
    init_class_switching_matrix = init_class_switch_matrix  # Alternative alias

    def set_class_switching_matrix(self, cs_matrix: np.ndarray):
        """
        Set the class switching probability matrix.

        Args:
            cs_matrix: K×K matrix where element (i,j) is the probability
                      that a job in class i switches to class j.
                      Each row should sum to 1.0.

        Example:
            >>> csmatrix = cs_node.init_class_switch_matrix()
            >>> csmatrix[0, 0] = 0.7
            >>> csmatrix[0, 1] = 0.3
            >>> cs_node.set_class_switching_matrix(csmatrix)
        """
        self._switch_matrix = np.asarray(cs_matrix)

    setClassSwitchingMatrix = set_class_switching_matrix

    def get_class_switching_matrix(self) -> Optional[np.ndarray]:
        """
        Get the current class switching matrix.

        Returns:
            np.ndarray or None: The K×K switching matrix, or None if not set
        """
        return self._switch_matrix

    getClassSwitchingMatrix = get_class_switching_matrix

    def set_switch_probability(self, from_class, to_class, probability: float):
        """
        Set the probability of switching from one class to another.

        This is a convenience method that handles indexing automatically.

        Args:
            from_class: Source job class (object or 0-based index)
            to_class: Target job class (object or 0-based index)
            probability: Switching probability (0.0 to 1.0)

        Example:
            >>> cs_node.set_switch_probability(class1, class2, 0.3)
            >>> cs_node.set_switch_probability(class1, class1, 0.7)
        """
        if self._switch_matrix is None:
            K = self._model.get_number_of_classes()
            self._switch_matrix = np.zeros((K, K))

        # Handle class objects or indices
        if hasattr(from_class, 'get_index0'):
            from_idx = from_class.get_index0()
        elif hasattr(from_class, 'get_index'):
            from_idx = from_class.get_index() - 1  # Convert 1-based to 0-based
        else:
            from_idx = int(from_class)

        if hasattr(to_class, 'get_index0'):
            to_idx = to_class.get_index0()
        elif hasattr(to_class, 'get_index'):
            to_idx = to_class.get_index() - 1  # Convert 1-based to 0-based
        else:
            to_idx = int(to_class)

        self._switch_matrix[from_idx, to_idx] = probability

    setSwitchProbability = set_switch_probability


class Cache(Station):
    """
    Multi-level cache node with configurable replacement strategy.

    A Cache node models content caching behavior where arriving jobs request
    items from a catalog. Items may be cached (hit) or fetched (miss), with
    jobs potentially switching classes based on hit/miss outcomes.

    The cache supports multiple levels (hierarchy) with different capacities
    and uses configurable replacement policies (LRU, FIFO, RR, SFIFO).

    Cache analysis can be performed using algorithms from api.cache:
    - cache_mva: Mean Value Analysis for hit/miss probabilities
    - cache_erec: Exact recursive computation
    - cache_spm: Singular Perturbation Method (approximate)

    Args:
        model: Network instance
        name: Cache node name
        num_items: Number of items in the catalog (n)
        item_level_cap: Capacity of each cache level (int or array)
        replacement_strategy: Cache replacement policy (LRU, FIFO, RR, SFIFO)

    Example:
        >>> model = Network('CacheModel')
        >>> cache = Cache(model, 'MyCache', num_items=100,
        ...               item_level_cap=10, replacement_strategy=ReplacementStrategy.LRU)
        >>> job_class = ClosedClass(model, 'Request', 5, delay)
        >>> hit_class = ClosedClass(model, 'Hit', 0, delay)
        >>> miss_class = ClosedClass(model, 'Miss', 0, delay)
        >>> cache.set_read(job_class, Zipf(1.2, 100))  # Zipf popularity
        >>> cache.set_hit_class(job_class, hit_class)
        >>> cache.set_miss_class(job_class, miss_class)

    Reference:
        - Cache algorithms: Che, H. et al. "Hierarchical Web Caching Systems"
        - TTL approximation: Fofack et al. "Analysis of TTL-based Cache Networks"
    """

    def __init__(
        self,
        model,
        name: str,
        num_items: int,
        item_level_cap: Union[int, np.ndarray, list],
        replacement_strategy: ReplacementStrategy = ReplacementStrategy.LRU
    ):
        """
        Initialize a Cache node.

        Args:
            model: Network instance
            name: Cache node name
            num_items: Total number of items in catalog (n >= 1)
            item_level_cap: Capacity per cache level. Can be:
                - int: Single-level cache with given capacity
                - array/list: Multi-level cache with capacity per level
            replacement_strategy: Replacement policy (default: LRU)

        Raises:
            ValueError: If num_items < 1 or item_level_cap invalid
        """
        super().__init__(NodeType.CACHE, name)
        self.set_model(model)
        model.add_node(self)

        # Validate parameters
        if num_items < 1:
            raise ValueError(f"num_items must be >= 1, got {num_items}")

        self._num_items = int(num_items)

        # Convert item_level_cap to array
        if isinstance(item_level_cap, (int, float)):
            self._item_level_cap = np.array([int(item_level_cap)])
        else:
            self._item_level_cap = np.array(item_level_cap, dtype=int).ravel()

        if np.any(self._item_level_cap < 1):
            raise ValueError("All cache level capacities must be >= 1")

        self._num_levels = len(self._item_level_cap)
        self._replacement_strategy = replacement_strategy

        # Cache nodes have infinite servers (instant service)
        self._number_of_servers = np.inf

        # Scheduling (caches are non-preemptive, use INF scheduling like Delay)
        self._sched_strategy = SchedStrategy.INF
        self._sched_policy = SchedStrategyType.NP

        # Hit/miss class mappings: input_class -> output_class
        self._hit_class = {}   # Dict[JobClass, JobClass]
        self._miss_class = {}  # Dict[JobClass, JobClass]

        # Read (popularity) distributions per class
        self._read_process = {}  # Dict[JobClass, Distribution]

        # Access probabilities matrix (optional)
        self._access_prob = None  # (num_items, num_classes)

        # Access cost matrix (optional)
        self._access_cost = None  # (num_items, num_classes)

        # Result storage (populated by solver)
        self._actual_hit_prob = None   # Actual hit probabilities
        self._actual_miss_prob = None  # Actual miss probabilities

    def is_station(self) -> bool:
        """
        Check if node is a station.

        Cache nodes are NOT stations in LINE's semantics, matching MATLAB behavior.
        They are stateful nodes that immediately process arriving jobs (class switching
        based on cache hit/miss) but don't queue or serve jobs like stations.

        Returns:
            False: Cache is a stateful non-station node
        """
        return False

    # =========================================================================
    # Properties
    # =========================================================================

    @property
    def num_items(self) -> int:
        """Get the number of items in the catalog."""
        return self._num_items

    @property
    def num_levels(self) -> int:
        """Get the number of cache levels."""
        return self._num_levels

    @property
    def item_level_cap(self) -> np.ndarray:
        """Get the capacity array for each cache level."""
        return self._item_level_cap.copy()

    @property
    def total_capacity(self) -> int:
        """Get the total cache capacity across all levels."""
        return int(np.sum(self._item_level_cap))

    @property
    def replacement_strategy(self) -> ReplacementStrategy:
        """Get the cache replacement strategy."""
        return self._replacement_strategy

    @replacement_strategy.setter
    def replacement_strategy(self, value: ReplacementStrategy):
        """Set the cache replacement strategy."""
        self._replacement_strategy = value
        self._invalidate_java()

    # =========================================================================
    # Hit/Miss Class Configuration
    # =========================================================================

    def set_hit_class(self, input_class: JobClass, output_class: JobClass) -> None:
        """
        Set the output class for cache hits.

        When a job of input_class experiences a cache hit, it transitions
        to output_class for further processing.

        Args:
            input_class: Incoming job class making the request
            output_class: Resulting job class after cache hit
        """
        self._hit_class[input_class] = output_class
        self._invalidate_java()

    def get_hit_class(self, input_class: JobClass) -> Optional[JobClass]:
        """Get the output class for cache hits."""
        return self._hit_class.get(input_class)

    def set_miss_class(self, input_class: JobClass, output_class: JobClass) -> None:
        """
        Set the output class for cache misses.

        When a job of input_class experiences a cache miss, it transitions
        to output_class for further processing.

        Args:
            input_class: Incoming job class making the request
            output_class: Resulting job class after cache miss
        """
        self._miss_class[input_class] = output_class
        self._invalidate_java()

    def get_miss_class(self, input_class: JobClass) -> Optional[JobClass]:
        """Get the output class for cache misses."""
        return self._miss_class.get(input_class)

    # =========================================================================
    # Read (Popularity) Distribution Configuration
    # =========================================================================

    def set_read(self, jobclass: JobClass, distribution) -> None:
        """
        Set the read (popularity) distribution for a job class.

        The distribution determines which items are requested by jobs
        of the given class. Typically a discrete distribution like Zipf.

        Args:
            jobclass: Job class making requests
            distribution: Popularity distribution (e.g., Zipf, Uniform)
        """
        self._read_process[jobclass] = distribution
        self._invalidate_java()

    def get_read(self, jobclass: JobClass):
        """Get the read distribution for a job class."""
        return self._read_process.get(jobclass)

    # =========================================================================
    # Access Probabilities (Alternative to Distribution)
    # =========================================================================

    def set_access_prob(self, access_prob: np.ndarray) -> None:
        """
        Set the access probability matrix directly.

        Alternative to set_read() for specifying item access probabilities
        as a matrix.

        Args:
            access_prob: (num_items, num_classes) probability matrix.
                         Each column should sum to 1.
        """
        self._access_prob = np.asarray(access_prob)
        self._invalidate_java()

    def get_access_prob(self) -> Optional[np.ndarray]:
        """Get the access probability matrix."""
        return self._access_prob

    # =========================================================================
    # Result Methods (Populated by Solver)
    # =========================================================================

    def set_result_hit_prob(self, hit_prob: np.ndarray) -> None:
        """Set the computed hit probability (called by solver)."""
        self._actual_hit_prob = np.asarray(hit_prob)

    def set_result_miss_prob(self, miss_prob: np.ndarray) -> None:
        """Set the computed miss probability (called by solver)."""
        self._actual_miss_prob = np.asarray(miss_prob)

    def get_hit_ratio(self) -> Optional[np.ndarray]:
        """
        Get the hit ratio (probability) for each class.

        Returns:
            Hit probability array, or None if not computed yet.
        """
        return self._actual_hit_prob

    def get_miss_ratio(self) -> Optional[np.ndarray]:
        """
        Get the miss ratio (probability) for each class.

        Returns:
            Miss probability array, or None if not computed yet.
        """
        return self._actual_miss_prob

    # =========================================================================
    # Utility Methods
    # =========================================================================

    def get_gamma_matrix(self, nclasses: int = 1) -> np.ndarray:
        """
        Build the gamma (access factor) matrix for cache analysis.

        The gamma matrix has shape (num_items, num_levels) and contains
        the cache access intensity factors used by cache_mva/cache_erec.

        Args:
            nclasses: Number of job classes

        Returns:
            Gamma matrix for cache analysis algorithms
        """
        n = self._num_items
        h = self._num_levels

        gamma = np.zeros((n, h))

        # If access_prob is set, use it to derive gamma
        if self._access_prob is not None:
            # Aggregate across classes
            access_agg = np.sum(self._access_prob, axis=1) if self._access_prob.ndim > 1 else self._access_prob
            # Normalize to get access intensities
            total = np.sum(access_agg)
            if total > 0:
                gamma[:, 0] = access_agg / total
        elif self._read_process:
            # Extract probabilities from read distributions (e.g., Zipf)
            access_probs = np.zeros(n)
            for jobclass, dist in self._read_process.items():
                if dist is not None:
                    # Check if it's a Zipf distribution
                    if hasattr(dist, '_s') and hasattr(dist, '_H'):
                        # Zipf distribution: P(k) = (1/k^s) / H
                        k = np.arange(1, n + 1)
                        probs = (1.0 / k ** dist._s) / dist._H
                        access_probs += probs
                    # Check if it's a DiscreteSampler with probabilities
                    elif hasattr(dist, '_probs'):
                        probs = np.asarray(dist._probs)
                        if len(probs) >= n:
                            access_probs += probs[:n]
                        else:
                            access_probs[:len(probs)] += probs
                    # Check if it has evalPMF method (generic discrete distribution)
                    elif hasattr(dist, 'evalPMF'):
                        for i in range(n):
                            access_probs[i] += dist.evalPMF(i + 1)  # 1-indexed items
            # Normalize
            total = np.sum(access_probs)
            if total > 0:
                gamma[:, 0] = access_probs / total
            else:
                gamma[:, 0] = 1.0 / n
        else:
            # Default: uniform access
            gamma[:, 0] = 1.0 / n

        # For multi-level caches, propagate to other levels
        for level in range(1, h):
            gamma[:, level] = gamma[:, level - 1]

        return gamma

    def get_capacity_vector(self) -> np.ndarray:
        """
        Get the cache capacity vector.

        Returns:
            Array of cache level capacities
        """
        return self._item_level_cap.copy()

    # =========================================================================
    # Aliases (MATLAB compatibility)
    # =========================================================================

    setHitClass = set_hit_class
    getHitClass = get_hit_class
    setMissClass = set_miss_class
    getMissClass = get_miss_class
    setRead = set_read
    getRead = get_read
    setAccessProb = set_access_prob
    getAccessProb = get_access_prob
    setResultHitProb = set_result_hit_prob
    setResultMissProb = set_result_miss_prob
    getHitRatio = get_hit_ratio
    getMissRatio = get_miss_ratio
    hit_ratio = get_hit_ratio
    miss_ratio = get_miss_ratio

    # Property accessors (MATLAB style)
    def getNumItems(self) -> int:
        """Get number of items (MATLAB compatibility)."""
        return self.num_items

    def getNumLevels(self) -> int:
        """Get number of cache levels (MATLAB compatibility)."""
        return self.num_levels

    def getCapacity(self) -> int:
        """Get total cache capacity (MATLAB compatibility)."""
        return self.total_capacity

    def getReplacementStrategy(self) -> ReplacementStrategy:
        """Get replacement strategy (MATLAB compatibility)."""
        return self.replacement_strategy

    def getItemLevelCap(self) -> np.ndarray:
        """Get per-level capacity array (MATLAB compatibility)."""
        return self.item_level_cap


class TimingStrategy(IntEnum):
    """
    Timing strategies for Petri net transitions.
    """
    TIMED = 0       # Timed transition with delay
    IMMEDIATE = 1   # Immediate transition (zero delay)


class Mode:
    """
    A firing mode for a Petri net transition.

    Modes define how a transition can fire, including:
    - Enabling conditions (required tokens from input places)
    - Inhibiting conditions (blocking tokens in places)
    - Firing outcomes (tokens produced to output places)
    - Timing distribution for the firing delay
    """

    def __init__(self, transition: 'Transition', name: str, index: int):
        """
        Initialize a Mode.

        Args:
            transition: Parent Transition object
            name: Mode name
            index: 1-based index of this mode
        """
        self._transition = transition
        self._name = name
        self._index = index  # 1-based index

    @property
    def name(self) -> str:
        """Get mode name."""
        return self._name

    @property
    def index(self) -> int:
        """Get 1-based index (MATLAB compatibility)."""
        return self._index

    def __index__(self) -> int:
        """Allow Mode to be used as an index (0-based for Python arrays)."""
        return self._index - 1

    def __int__(self) -> int:
        """Convert to int (1-based index)."""
        return self._index

    def __repr__(self) -> str:
        return f"Mode('{self._name}', index={self._index})"


class Place(Station):
    """
    Place node for Stochastic Petri Nets.

    A Place represents a location where tokens (jobs) can accumulate.
    Places have infinite capacity (like Delay nodes) and infinite servers.

    In queueing network terms, a Place is similar to a delay station
    but tokens are processed when they move to Transitions.
    """

    def __init__(self, model, name: str):
        """
        Initialize a Place node.

        Args:
            model: Network instance
            name: Place name
        """
        super().__init__(NodeType.PLACE, name)
        self.set_model(model)
        model.add_node(self)

        # Places have infinite capacity and servers (like delay nodes)
        self._number_of_servers = np.inf
        self._capacity = np.inf
        self._sched_strategy = SchedStrategy.INF
        self._class_cap = {}  # Per-class capacity
        self._state = None  # Token state per class

    def set_class_capacity(self, jobclass: JobClass, capacity: float) -> None:
        """
        Set per-class capacity limit.

        Args:
            jobclass: Job class
            capacity: Capacity for this class
        """
        self._class_cap[jobclass] = capacity
        self._invalidate_java()

    setClassCapacity = set_class_capacity

    def set_state(self, state) -> None:
        """
        Set initial token state for this place.

        Args:
            state: Array of initial token counts per class
        """
        self._state = np.asarray(state)
        self._state_explicitly_set = True
        self._invalidate_java()

    setState = set_state

    def get_state(self):
        """Get the current state of this place."""
        return self._state

    getState = get_state


class Transition(StatefulNode):
    """
    Transition node for Stochastic Petri Nets.

    A Transition consumes tokens from input Places and produces tokens
    to output Places according to defined firing modes.

    Each Transition can have multiple firing modes, each with:
    - Enabling conditions: required tokens from input places
    - Inhibiting conditions: blocking conditions from places
    - Firing outcomes: tokens produced to output places
    - Timing distribution: delay before firing
    - Priority and weight for conflict resolution
    """

    def __init__(self, model, name: str):
        """
        Initialize a Transition node.

        Args:
            model: Network instance
            name: Transition name
        """
        super().__init__(NodeType.TRANSITION, name)
        self.set_model(model)
        model.add_node(self)

        self._capacity = np.inf

        # Mode-indexed data structures (lists, 0-based internal indexing)
        self._mode_names = []  # List of mode names
        self._modes = []  # List of Mode objects
        self._enabling_conditions = []  # List of (nnodes x nclasses) matrices
        self._inhibiting_conditions = []  # List of (nnodes x nclasses) matrices
        self._firing_outcomes = []  # List of (nnodes x nclasses) matrices
        self._number_of_servers = []  # List of server counts per mode
        self._timing_strategies = []  # List of TimingStrategy per mode
        self._distributions = []  # List of distributions per mode
        self._firing_priorities = []  # List of priorities per mode
        self._firing_weights = []  # List of weights per mode

    def add_mode(self, mode_name: str) -> Mode:
        """
        Add a new firing mode to this transition.

        Args:
            mode_name: Name for the new mode

        Returns:
            Mode object representing the new mode
        """
        nclasses = self._model.get_number_of_classes()
        nnodes = self._model.get_number_of_nodes()

        self._mode_names.append(mode_name)

        # Initialize matrices for this mode
        self._enabling_conditions.append(np.zeros((nnodes, nclasses)))
        self._inhibiting_conditions.append(np.full((nnodes, nclasses), np.inf))
        self._firing_outcomes.append(np.zeros((nnodes, nclasses)))

        # Initialize mode parameters
        self._number_of_servers.append(1)
        self._timing_strategies.append(TimingStrategy.TIMED)
        self._distributions.append(None)  # Will be set by setDistribution
        self._firing_priorities.append(1.0)
        self._firing_weights.append(1.0)

        # Create and store Mode object
        mode = Mode(self, mode_name, len(self._modes) + 1)  # 1-based index
        self._modes.append(mode)

        self._invalidate_java()
        return mode

    addMode = add_mode

    def get_number_of_modes(self) -> int:
        """Get the number of firing modes."""
        return len(self._mode_names)

    getNumberOfModes = get_number_of_modes

    def get_modes(self) -> list:
        """Get list of Mode objects."""
        return self._modes

    getModes = get_modes

    def _get_mode_index(self, mode) -> int:
        """
        Get 0-based index from a mode identifier.

        Args:
            mode: Mode object, 1-based index, or 0-based index

        Returns:
            0-based index for internal arrays
        """
        if isinstance(mode, Mode):
            return mode._index - 1  # Mode._index is 1-based
        elif hasattr(mode, '__index__'):
            # Could be a Mode or int-like object
            idx = mode.__index__()
            # If the index is >= len(modes), assume it's 1-based
            if idx >= len(self._modes):
                return idx - 1  # 1-based to 0-based
            return idx
        else:
            # Assume 1-based integer index (MATLAB convention)
            return int(mode) - 1

    def set_enabling_conditions(self, mode, jobclass: JobClass,
                                 input_node, enabling_condition: int) -> None:
        """
        Set enabling conditions for a mode.

        Args:
            mode: Mode object or index (1-based)
            jobclass: Job class
            input_node: Input Place node
            enabling_condition: Number of tokens required
        """
        mode_idx = self._get_mode_index(mode)

        # Get node index
        if hasattr(input_node, 'get_index0'):
            node_idx = input_node.get_index0()
        elif hasattr(input_node, 'get_index'):
            node_idx = input_node.get_index() - 1
        else:
            node_idx = self._model.get_node_index(input_node.name) - 1

        # Get class index
        if hasattr(jobclass, 'get_index0'):
            class_idx = jobclass.get_index0()
        elif hasattr(jobclass, 'get_index'):
            class_idx = jobclass.get_index() - 1
        else:
            class_idx = int(jobclass)

        self._enabling_conditions[mode_idx][node_idx, class_idx] = enabling_condition
        self._invalidate_java()

    setEnablingConditions = set_enabling_conditions

    def set_inhibiting_conditions(self, mode, jobclass: JobClass,
                                   input_node, inhibiting_condition: int) -> None:
        """
        Set inhibiting conditions for a mode.

        The transition cannot fire if the place has >= inhibiting_condition tokens.

        Args:
            mode: Mode object or index (1-based)
            jobclass: Job class
            input_node: Input Place node
            inhibiting_condition: Number of tokens that inhibit firing
        """
        mode_idx = self._get_mode_index(mode)

        # Get node index
        if hasattr(input_node, 'get_index0'):
            node_idx = input_node.get_index0()
        elif hasattr(input_node, 'get_index'):
            node_idx = input_node.get_index() - 1
        else:
            node_idx = self._model.get_node_index(input_node.name) - 1

        # Get class index
        if hasattr(jobclass, 'get_index0'):
            class_idx = jobclass.get_index0()
        elif hasattr(jobclass, 'get_index'):
            class_idx = jobclass.get_index() - 1
        else:
            class_idx = int(jobclass)

        self._inhibiting_conditions[mode_idx][node_idx, class_idx] = inhibiting_condition
        self._invalidate_java()

    setInhibitingConditions = set_inhibiting_conditions

    def set_firing_outcome(self, mode, jobclass: JobClass,
                           output_node, firing_outcome: int) -> None:
        """
        Set firing outcome for a mode.

        Args:
            mode: Mode object or index (1-based)
            jobclass: Job class
            output_node: Output Place or Sink node
            firing_outcome: Number of tokens produced
        """
        mode_idx = self._get_mode_index(mode)

        # Get node index
        if hasattr(output_node, 'get_index0'):
            node_idx = output_node.get_index0()
        elif hasattr(output_node, 'get_index'):
            node_idx = output_node.get_index() - 1
        else:
            node_idx = self._model.get_node_index(output_node.name) - 1

        # Get class index
        if hasattr(jobclass, 'get_index0'):
            class_idx = jobclass.get_index0()
        elif hasattr(jobclass, 'get_index'):
            class_idx = jobclass.get_index() - 1
        else:
            class_idx = int(jobclass)

        self._firing_outcomes[mode_idx][node_idx, class_idx] = firing_outcome
        self._invalidate_java()

    setFiringOutcome = set_firing_outcome

    def set_distribution(self, mode, distribution) -> None:
        """
        Set firing time distribution for a mode.

        Args:
            mode: Mode object or index (1-based)
            distribution: Timing distribution (Exp, Erlang, etc.)
        """
        mode_idx = self._get_mode_index(mode)
        self._distributions[mode_idx] = distribution
        self._invalidate_java()

    setDistribution = set_distribution

    def get_distribution(self, mode):
        """Get firing time distribution for a mode."""
        mode_idx = self._get_mode_index(mode)
        return self._distributions[mode_idx]

    getDistribution = get_distribution

    def set_number_of_servers(self, mode, num_servers: int) -> None:
        """
        Set number of servers for a mode.

        Args:
            mode: Mode object or index (1-based)
            num_servers: Number of servers (or GlobalConstants.MaxInt for infinite)
        """
        mode_idx = self._get_mode_index(mode)
        self._number_of_servers[mode_idx] = num_servers
        self._invalidate_java()

    setNumberOfServers = set_number_of_servers

    def set_timing_strategy(self, mode, timing_strategy: TimingStrategy) -> None:
        """
        Set timing strategy for a mode.

        Args:
            mode: Mode object or index (1-based)
            timing_strategy: TIMED or IMMEDIATE
        """
        mode_idx = self._get_mode_index(mode)
        self._timing_strategies[mode_idx] = timing_strategy
        self._invalidate_java()

    setTimingStrategy = set_timing_strategy

    def set_firing_priorities(self, mode, priority: float) -> None:
        """
        Set firing priority for a mode.

        Args:
            mode: Mode object or index (1-based)
            priority: Priority value (higher = more priority)
        """
        mode_idx = self._get_mode_index(mode)
        self._firing_priorities[mode_idx] = priority
        self._invalidate_java()

    setFiringPriorities = set_firing_priorities

    def set_firing_weights(self, mode, weight: float) -> None:
        """
        Set firing weight for a mode.

        Args:
            mode: Mode object or index (1-based)
            weight: Weight for probabilistic selection among enabled modes
        """
        mode_idx = self._get_mode_index(mode)
        self._firing_weights[mode_idx] = weight
        self._invalidate_java()

    setFiringWeights = set_firing_weights

    def set_mode_names(self, mode, mode_name: str) -> None:
        """
        Set name for a mode.

        Args:
            mode: Mode object or index (1-based)
            mode_name: New name for the mode
        """
        mode_idx = self._get_mode_index(mode)
        self._mode_names[mode_idx] = mode_name
        self._invalidate_java()

    setModeNames = set_mode_names

    def get_service_rates(self):
        """
        Get service rates for all modes.

        Returns:
            Tuple of (map, mu, phi) lists for PH-type service rates
        """
        nmodes = self.get_number_of_modes()
        map_list = [None] * nmodes
        mu_list = [None] * nmodes
        phi_list = [None] * nmodes

        for m in range(nmodes):
            dist = self._distributions[m]
            if dist is None:
                map_list[m] = [[np.nan], [np.nan]]
                mu_list[m] = np.nan
                phi_list[m] = np.nan
            elif hasattr(dist, 'getProcess'):
                map_list[m] = dist.getProcess()
                mu_list[m] = dist.getMu() if hasattr(dist, 'getMu') else None
                phi_list[m] = dist.getPhi() if hasattr(dist, 'getPhi') else None
            elif hasattr(dist, 'getMean'):
                rate = 1.0 / dist.getMean()
                map_list[m] = [[-rate], [rate]]
                mu_list[m] = [rate]
                phi_list[m] = [1]
            else:
                map_list[m] = [[np.nan], [np.nan]]
                mu_list[m] = np.nan
                phi_list[m] = np.nan

        return map_list, mu_list, phi_list

    getServiceRates = get_service_rates


class Logger(Node):
    """
    Logger node for recording job passage.

    A Logger node records arrival and departure timestamps for jobs
    passing through it. Used internally by getCdfRespT to collect
    response time samples via transient simulation.

    Ported from MATLAB's Logger class in matlab/src/lang/nodes/Logger.m
    """

    def __init__(self, model, name: str, log_file_name: str):
        """
        Initialize a Logger node.

        Args:
            model: Network instance
            name: Logger name
            log_file_name: Full path to the log file
        """
        import os

        super().__init__(NodeType.LOGGER, name)
        self.set_model(model)
        model.add_node(self)

        # Parse file name and path
        dirname, basename = os.path.split(log_file_name)
        self._file_name = basename
        self._file_path = dirname if dirname else model.get_log_path()

        # Logging options (matching MATLAB defaults)
        self._want_start_time = False
        self._want_logger_name = False
        self._want_timestamp = True
        self._want_job_id = True
        self._want_job_class = True
        self._want_time_same_class = False
        self._want_time_any_class = False

        # Scheduling (Logger uses FCFS, non-preemptive)
        self._sched_policy = SchedStrategyType.NP
        self._sched_strategy = SchedStrategy.FCFS
        self._capacity = float('inf')

    @property
    def file_name(self) -> str:
        """Get the log file name."""
        return self._file_name

    @property
    def file_path(self) -> str:
        """Get the log file path."""
        return self._file_path

    # Getters for logging options
    def get_start_time(self) -> bool:
        return self._want_start_time

    def get_logger_name(self) -> bool:
        return self._want_logger_name

    def get_timestamp(self) -> bool:
        return self._want_timestamp

    def get_job_id(self) -> bool:
        return self._want_job_id

    def get_job_class(self) -> bool:
        return self._want_job_class

    def get_time_same_class(self) -> bool:
        return self._want_time_same_class

    def get_time_any_class(self) -> bool:
        return self._want_time_any_class

    # Setters for logging options
    def set_start_time(self, value: bool) -> None:
        self._want_start_time = value

    def set_logger_name(self, value: bool) -> None:
        self._want_logger_name = value

    def set_timestamp(self, value: bool) -> None:
        self._want_timestamp = value

    def set_job_id(self, value: bool) -> None:
        self._want_job_id = value

    def set_job_class(self, value: bool) -> None:
        self._want_job_class = value

    def set_time_same_class(self, value: bool) -> None:
        self._want_time_same_class = value

    def set_time_any_class(self, value: bool) -> None:
        self._want_time_any_class = value

    def set_prob_routing(self, jobclass: JobClass, destination, probability: float) -> None:
        """
        Set probabilistic routing to a destination.

        Args:
            jobclass: Job class
            destination: Destination node
            probability: Routing probability
        """
        self.set_routing(jobclass, RoutingStrategy.PROB, destination, probability)

    # CamelCase aliases
    getStartTime = get_start_time
    getLoggerName = get_logger_name
    getTimestamp = get_timestamp
    getJobID = get_job_id
    getJobClass = get_job_class
    getTimeSameClass = get_time_same_class
    getTimeAnyClass = get_time_any_class
    setStartTime = set_start_time
    setLoggerName = set_logger_name
    setTimestamp = set_timestamp
    setJobID = set_job_id
    setJobClass = set_job_class
    setTimeSameClass = set_time_same_class
    setTimeAnyClass = set_time_any_class
    setProbRouting = set_prob_routing
    fileName = property(lambda self: self._file_name)
    filePath = property(lambda self: self._file_path)


__all__ = [
    'Queue', 'Delay', 'Source', 'Sink', 'Fork', 'Join',
    'Router', 'ClassSwitch', 'Cache', 'ReplacementStrategy',
    'Place', 'Transition', 'Mode', 'TimingStrategy', 'Logger'
]
