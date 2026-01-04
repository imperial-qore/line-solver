"""
Base classes for LINE native Python implementation.

This module provides the foundation for native Python queueing network models,
including Element, NetworkElement, Node, StatefulNode, and Station classes.
Ported from MATLAB implementation in matlab/src/lang/@*
"""

from abc import ABC, abstractmethod
from enum import IntEnum
from typing import Optional, Dict, List, Any, Tuple
import numpy as np


class ElementType(IntEnum):
    """Enumeration of element types in LINE networks."""
    MODEL = 0
    NODE = 1
    CLASS = 2
    REWARD = 3
    ITEM = 4


class NodeType(IntEnum):
    """Enumeration of node types in queueing networks."""
    SOURCE = 0
    SINK = 1
    QUEUE = 2
    DELAY = 3
    FORK = 4
    JOIN = 5
    CACHE = 6
    ROUTER = 7
    CLASSSWITCH = 8
    PLACE = 9
    TRANSITION = 10


class SchedStrategy(IntEnum):
    """Enumeration of scheduling strategies.

    Values aligned with api/sn/network_struct.py for compatibility.
    """
    FCFS = 0      # First-Come First-Served
    LCFS = 1      # Last-Come First-Served
    LCFSPR = 2    # LCFS Preemptive Resume
    LCFSPI = 3    # LCFS Preemptive Identical
    PS = 4        # Processor Sharing
    DPS = 5       # Discriminatory PS
    GPS = 6       # Generalized PS
    INF = 7       # Infinite Server (Delay)
    RAND = 8      # Random
    HOL = 9       # Head of Line
    SEPT = 10     # Shortest Expected Processing Time
    LEPT = 11     # Longest Expected Processing Time
    SIRO = 12     # Service in Random Order
    SJF = 13      # Shortest Job First
    LJF = 14      # Longest Job First
    POLLING = 15  # Polling
    EXT = 16      # External
    LPS = 17      # Least Progress Scheduling
    SETF = 18     # Shortest Elapsed Time First
    DPSPRIO = 19  # DPS with Priority
    GPSPRIO = 20  # GPS with Priority
    PSPRIO = 21   # PS with Priority
    FCFSPR = 22   # FCFS with Preemption (for compatibility)
    EDF = 23      # Earliest Deadline First
    FORK = 24     # Fork node
    JOIN = 25     # Join node


class RoutingStrategy(IntEnum):
    """Enumeration of routing strategies.

    Values aligned with api/sn/network_struct.py for JMT compatibility.
    """
    DISABLED = 0
    RAND = 1      # Random routing
    PROB = 2      # Probabilistic routing
    RROBIN = 3    # Round-robin
    WRROBIN = 4   # Weighted round-robin
    KCHOICES = 7  # K-choices policy


class JobClassType(IntEnum):
    """Enumeration of job class types."""
    OPEN = 0
    CLOSED = 1
    SIGNAL = 2


class SchedStrategyType(IntEnum):
    """Enumeration of scheduling policy types."""
    NP = 0  # Non-preemptive
    PR = 1  # Preemptive


class JoinStrategy(IntEnum):
    """Enumeration of join strategies."""
    STD = 0       # Standard (AND-join)
    QUORUM = 1    # Quorum
    CANDJOIN = 2  # Cache AND-join


class DropStrategy(IntEnum):
    """Enumeration of drop strategies for finite capacity."""
    DROP = 0      # Drop arriving job
    BAS = 1       # Block After Service
    WaitingQueue = 2  # Block in queue
    WAITQ = 2     # MATLAB-compatible alias for WaitingQueue


class ReplacementStrategy(IntEnum):
    """
    Cache replacement strategies.

    Determines which item to evict when the cache is full and a new item
    needs to be stored.

    Attributes:
        RR: Random Replacement - evict a random item
        FIFO: First-In-First-Out - evict oldest item
        SFIFO: Segmented FIFO - FIFO with segmented structure
        LRU: Least Recently Used - evict least recently accessed item
    """
    RR = 0      # Random Replacement
    FIFO = 1    # First-In-First-Out
    SFIFO = 2   # Segmented FIFO
    LRU = 3     # Least Recently Used


class HeteroSchedPolicy(IntEnum):
    """
    Scheduling policies for heterogeneous multiserver queues.

    These policies determine how jobs are assigned to server types in
    heterogeneous multiserver queues when a job's class is compatible
    with multiple server types.

    Attributes:
        ORDER: Assign to first available compatible server type (in definition order)
        ALIS: Assign Longest Idle Server (round-robin with busy servers at back)
        ALFS: Assign Longest Free Server (fairness sorting by coverage)
        FAIRNESS: Fair distribution across compatible server types
        FSF: Fastest Server First (based on expected service time)
        RAIS: Random Available Idle Server
    """
    ORDER = 0      # First available compatible server type
    ALIS = 1       # Assign Longest Idle Server
    ALFS = 2       # Assign Longest Free Server
    FAIRNESS = 3   # Fair distribution
    FSF = 4        # Fastest Server First
    RAIS = 5       # Random Available Idle Server

    @classmethod
    def from_text(cls, text: str) -> 'HeteroSchedPolicy':
        """Convert text to HeteroSchedPolicy constant."""
        mapping = {
            'ORDER': cls.ORDER,
            'ALIS': cls.ALIS,
            'ALFS': cls.ALFS,
            'FAIRNESS': cls.FAIRNESS,
            'FSF': cls.FSF,
            'RAIS': cls.RAIS,
        }
        upper_text = text.upper()
        if upper_text in mapping:
            return mapping[upper_text]
        raise ValueError(f"Unknown HeteroSchedPolicy: {text}")

    def to_text(self) -> str:
        """Convert HeteroSchedPolicy constant to text."""
        return self.name

    # Aliases for from_text
    from_string = from_text
    fromString = from_text


class Element(ABC):
    """Abstract base class for all LINE network elements."""

    def __init__(self, element_type: ElementType, name: str = ""):
        """
        Initialize an element.

        Args:
            element_type: Type of element (from ElementType enum)
            name: Name of the element
        """
        self._element_type = element_type
        self._name = name
        self._index = -1  # 0-based index, -1 if not assigned
        self._java_obj = None  # Cached object

    @property
    def name(self) -> str:
        """Get the element name."""
        return self._name

    @name.setter
    def name(self, value: str) -> None:
        """Set the element name."""
        self._name = value
        self._invalidate_java()

    def get_name(self) -> str:
        """Get the element name ."""
        return self._name

    def getName(self) -> str:
        """Get the element name ."""
        return self._name

    @property
    def element_type(self) -> ElementType:
        """Get the element type."""
        return self._element_type

    def get_index(self) -> int:
        """Get 1-based index (MATLAB compatibility)."""
        return self._index + 1 if self._index >= 0 else -1

    def get_index0(self) -> int:
        """Get 0-based index (Python native)."""
        return self._index

    def _set_index(self, idx: int) -> None:
        """Set 0-based index (internal use only)."""
        self._index = idx

    def _invalidate_java(self) -> None:
        """Invalidate cached object."""
        self._java_obj = None

    def to_java(self):
        """Convert to object for interop (future use)"""
        raise NotImplementedError("Subclass must implement to_java()")

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}('{self._name}', index={self.get_index()})"


class NetworkElement(Element):
    """Base class for elements that belong to a network."""

    def __init__(self, element_type: ElementType, name: str = ""):
        """
        Initialize a network element.

        Args:
            element_type: Type of element
            name: Name of the element
        """
        super().__init__(element_type, name)


class JobClass(NetworkElement):
    """Abstract base class for job classes."""

    def __init__(self, model_or_type, name_or_jobclass_type: str = None, jobclass_type_or_prio = None, prio_or_deadline = None):
        """
        Initialize a job class.

        Supports two signatures for compatibility:
        1. (model, name, JobClassType, priority) - from classes.py
        2. (jobclass_type, name, priority, deadline) - native implementation

        Args:
            model_or_type: Model instance or JobClassType
            name_or_jobclass_type: Class name or JobClassType
            jobclass_type_or_prio: JobClassType or priority
            prio_or_deadline: Priority or deadline
        """
        # Determine which signature is being used
        if isinstance(model_or_type, JobClassType):
            # Signature 2: Native (jobclass_type, name, priority, deadline)
            jobclass_type = model_or_type
            name = name_or_jobclass_type
            priority = jobclass_type_or_prio if jobclass_type_or_prio is not None else 0
            deadline = prio_or_deadline if prio_or_deadline is not None else np.inf
            model = None
        else:
            # Signature 1: From classes.py (model, name, jobclass_type, priority)
            model = model_or_type
            name = name_or_jobclass_type
            jobclass_type = jobclass_type_or_prio
            priority = prio_or_deadline if prio_or_deadline is not None else 0
            deadline = np.inf

        super().__init__(ElementType.CLASS, name)
        self._model = model
        self._jobclass_type = jobclass_type
        self._priority = priority
        self._deadline = deadline
        self._refstat = None  # Reference station (Queue for closed, Source for open)
        self._patience_distribution = None
        self._patience_type = None
        self._reply_signal_class = None
        self._completes = True  # Whether this class completes a visit (for cycle time calculation)

    @property
    def jobclass_type(self) -> JobClassType:
        """Get the job class type."""
        return self._jobclass_type

    def __index__(self) -> int:
        """Return zero-based index for Python array indexing."""
        return self._index

    @property
    def priority(self) -> int:
        """Get scheduling priority."""
        return self._priority

    @priority.setter
    def priority(self, value: int) -> None:
        """Set scheduling priority."""
        self._priority = value
        self._invalidate_java()

    @property
    def deadline(self) -> float:
        """Get relative deadline."""
        return self._deadline

    @deadline.setter
    def deadline(self, value: float) -> None:
        """Set relative deadline."""
        self._deadline = value
        self._invalidate_java()

    @property
    def completes(self) -> bool:
        """
        Get whether this class completes a visit.

        When completes=True (default), visiting this class contributes to
        the system response time. When completes=False, the class represents
        an intermediate step that doesn't count towards completion.
        """
        return self._completes

    @completes.setter
    def completes(self, value: bool) -> None:
        """Set whether this class completes a visit."""
        self._completes = value
        self._invalidate_java()

    def set_reference_station(self, station) -> None:
        """
        Set reference station for this class.

        Args:
            station: Reference Station node
        """
        self._refstat = station
        self._invalidate_java()

    def get_reference_station(self):
        """Get reference station."""
        return self._refstat

    def is_reference_station(self, node) -> bool:
        """Check if a node is the reference station."""
        return self._refstat is node

    def set_patience(self, patience_type: str, distribution) -> None:
        """
        Set patience distribution for reneging/balking.

        Args:
            patience_type: 'reneging' or 'balking'
            distribution: Patience distribution
        """
        self._patience_type = patience_type
        self._patience_distribution = distribution
        self._invalidate_java()

    def get_patience(self):
        """Get patience distribution (if any)."""
        return self._patience_distribution

    def has_patience(self) -> bool:
        """Check if class has patience."""
        return self._patience_distribution is not None

    def set_reply_signal_class(self, reply_class) -> None:
        """Set reply signal class for synchronous calls."""
        self._reply_signal_class = reply_class
        self._invalidate_java()

    def get_reply_signal_class(self):
        """Get reply signal class."""
        return self._reply_signal_class

    def getPriority(self) -> int:
        """Get scheduling priority ."""
        return self._priority

    def setPriority(self, value: int) -> None:
        """Set scheduling priority ."""
        self._priority = value
        self._invalidate_java()

    def getDeadline(self) -> float:
        """Get relative deadline ."""
        return self._deadline

    def setDeadline(self, value: float) -> None:
        """Set relative deadline ."""
        self._deadline = value
        self._invalidate_java()

    def getCompletes(self) -> bool:
        """Get whether this class completes a visit ."""
        return self._completes

    def setCompletes(self, value: bool) -> None:
        """Set whether this class completes a visit ."""
        self._completes = value
        self._invalidate_java()

    def getIndex(self) -> int:
        """Alias for get_index (CamelCase)."""
        return self.get_index()

    def setReferenceStation(self, station) -> None:
        """Alias for set_reference_station (CamelCase)."""
        return self.set_reference_station(station)

    def getReferenceStation(self):
        """Alias for get_reference_station (CamelCase)."""
        return self.get_reference_station()

    def isReferenceStation(self, station) -> bool:
        """Alias for is_reference_station (CamelCase)."""
        return self.is_reference_station(station)

    def setPatience(self, dist) -> None:
        """Alias for set_patience (CamelCase)."""
        return self.set_patience(dist)

    def getPatience(self):
        """Alias for get_patience (CamelCase)."""
        return self.get_patience()

    def hasPatience(self) -> bool:
        """Alias for has_patience (CamelCase)."""
        return self.has_patience()

    def setReplySignalClass(self, class_obj) -> None:
        """Alias for set_reply_signal_class (CamelCase)."""
        return self.set_reply_signal_class(class_obj)

    def getReplySignalClass(self):
        """Alias for get_reply_signal_class (CamelCase)."""
        return self.get_reply_signal_class()

    # Snake_case aliases
    get_priority = getPriority
    set_priority = setPriority
    get_deadline = getDeadline
    set_deadline = setDeadline
    get_completes = getCompletes
    set_completes = setCompletes


class Node(NetworkElement):
    """Abstract base class for network nodes."""

    def __init__(self, node_type: NodeType, name: str = ""):
        """
        Initialize a node.

        Args:
            node_type: Type of node (from NodeType enum)
            name: Node name
        """
        super().__init__(ElementType.NODE, name)
        self._node_type = node_type
        self._model = None
        self._routing_strategies = {}  # Dict[JobClass, RoutingStrategy]
        self._routing_params = {}  # Dict[JobClass, routing_parameters]

    def __index__(self) -> int:
        """Return node/station index for Python array indexing.

        Returns station index if available (for Station subclasses),
        otherwise returns node index.
        """
        if hasattr(self, '_station_index') and self._station_index >= 0:
            return self._station_index
        return self._index

    @property
    def node_type(self) -> NodeType:
        """Get the node type."""
        return self._node_type

    def set_model(self, model) -> None:
        """
        Link this node to a network model.

        Args:
            model: Network instance
        """
        self._model = model
        self._invalidate_java()

    def get_model(self):
        """Get the parent network model."""
        return self._model

    def is_stateful(self) -> bool:
        """Check if node is stateful (can have jobs)."""
        return False

    def is_station(self) -> bool:
        """Check if node is a station (can serve jobs)."""
        return False

    def link(self, node_to) -> None:
        """
        Create a link from this node to another node.

        Args:
            node_to: Destination node
        """
        if self._model is None:
            raise ValueError(f"[{self.name}] Node not linked to network")
        self._model.add_link(self, node_to)

    def set_routing(self, jobclass: JobClass, strategy: RoutingStrategy, *params) -> None:
        """
        Set routing strategy for a job class.

        Args:
            jobclass: Job class
            strategy: Routing strategy (RAND, PROB, etc.)
            params: Additional parameters (depends on strategy)
        """
        self._routing_strategies[jobclass] = strategy
        if params:
            self._routing_params[jobclass] = params
        self._invalidate_java()

    def get_routing(self, jobclass: JobClass) -> Tuple[RoutingStrategy, tuple]:
        """
        Get routing strategy for a job class.

        Args:
            jobclass: Job class

        Returns:
            Tuple of (strategy, parameters)
        """
        strategy = self._routing_strategies.get(jobclass, RoutingStrategy.RAND)
        params = self._routing_params.get(jobclass, ())
        return strategy, params

    # CamelCase aliases
    setModel = set_model
    getModel = get_model
    isStateful = is_stateful
    isStation = is_station
    setRouting = set_routing
    getRouting = get_routing


class StatefulNode(Node):
    """Abstract base class for nodes that maintain state."""

    def __init__(self, node_type: NodeType, name: str = ""):
        """
        Initialize a stateful node.

        Args:
            node_type: Type of node
            name: Node name
        """
        super().__init__(node_type, name)
        self._state = None  # Current state vector
        self._state_space = None  # Possible states

    @property
    def state(self) -> Optional[np.ndarray]:
        """Get current state."""
        return self._state

    @state.setter
    def state(self, value: np.ndarray) -> None:
        """Set current state."""
        self._state = np.array(value)
        self._invalidate_java()

    def get_state_space(self) -> Optional[List[np.ndarray]]:
        """Get state space."""
        return self._state_space

    def set_state_space(self, space: List[np.ndarray]) -> None:
        """
        Set state space for this node.

        Args:
            space: List of possible states
        """
        self._state_space = space
        self._invalidate_java()

    def is_stateful(self) -> bool:
        """Check if node is stateful."""
        return True

    def set_prob_routing(self, jobclass: JobClass, destination, prob: float) -> None:
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

    def get_prob_routing(self, jobclass: JobClass):
        """Get probabilistic routing for a job class."""
        if not hasattr(self, '_prob_routing'):
            return {}
        return self._prob_routing.get(jobclass, {})

    def setState(self, value: np.ndarray) -> None:
        """Set current state ."""
        self._state = np.array(value)
        self._invalidate_java()

    def getState(self) -> Optional[np.ndarray]:
        """Get current state ."""
        return self._state

    def getStateSpace(self) -> Optional[List[np.ndarray]]:
        """Get state space ."""
        return self._state_space

    def setStateSpace(self, space: List[np.ndarray]) -> None:
        """Set state space ."""
        self._state_space = space
        self._invalidate_java()

    # CamelCase aliases
    setProbRouting = set_prob_routing
    getProbRouting = get_prob_routing

    # Snake_case aliases
    set_state = setState
    get_state = getState
    set_state_space = setStateSpace


class Station(StatefulNode):
    """Abstract base class for queueing stations."""

    def __init__(self, node_type: NodeType, name: str = ""):
        """
        Initialize a station.

        Args:
            node_type: Type of station
            name: Station name
        """
        super().__init__(node_type, name)
        self._number_of_servers = 1
        self._capacity = np.inf
        self._class_capacity = {}  # Dict[JobClass, capacity]
        self._drop_rule = DropStrategy.DROP
        self._load_depend_scaling = None
        self._class_depend_scaling = None
        self._joint_depend_scaling = None
        self._station_index = -1  # Index among stations

    @property
    def number_of_servers(self) -> int:
        """Get number of servers."""
        return self._number_of_servers

    @number_of_servers.setter
    def number_of_servers(self, value: int) -> None:
        """
        Set number of servers.

        Args:
            value: Number of servers (1, >1, or np.inf for infinite)
        """
        if value < 1 and not np.isinf(value):
            raise ValueError(f"[{self.name}] Number of servers must be >= 1, got {value}")
        self._number_of_servers = value
        self._invalidate_java()

    @property
    def capacity(self) -> float:
        """Get total capacity."""
        return self._capacity

    @capacity.setter
    def capacity(self, value: float) -> None:
        """
        Set total capacity.

        Args:
            value: Total capacity (>0 or infinity)
        """
        if value <= 0 and not np.isinf(value):
            raise ValueError(f"[{self.name}] Capacity must be > 0, got {value}")
        self._capacity = value
        self._invalidate_java()

    def set_capacity(self, capacity: float) -> None:
        """Set total capacity (MATLAB compatibility)."""
        self.capacity = capacity

    def set_number_of_servers(self, value: int) -> None:
        """
        Set number of servers.

        Args:
            value: Number of servers (1, >1, or np.inf for infinite)
        """
        self.number_of_servers = value

    # CamelCase alias
    setNumberOfServers = set_number_of_servers
    set_num_servers = set_number_of_servers  # Short alias

    def get_capacity(self) -> float:
        """Get total capacity (MATLAB compatibility)."""
        return self._capacity

    # CamelCase aliases for capacity
    setCapacity = set_capacity
    getCapacity = get_capacity
    set_cap = set_capacity
    setCap = set_capacity

    def set_class_capacity(self, jobclass: JobClass, capacity: float) -> None:
        """
        Set per-class capacity limit.

        Args:
            jobclass: Job class
            capacity: Capacity for this class
        """
        if capacity <= 0 and not np.isinf(capacity):
            raise ValueError(f"[{self.name}] Class capacity must be > 0, got {capacity}")
        self._class_capacity[jobclass] = capacity
        self._invalidate_java()

    def get_class_capacity(self, jobclass: JobClass) -> float:
        """Get per-class capacity limit."""
        return self._class_capacity.get(jobclass, self._capacity)

    def set_load_dependence(self, alpha: np.ndarray) -> None:
        """
        Set load-dependent service rates.

        Args:
            alpha: Scaling factors indexed by population
        """
        self._load_depend_scaling = np.array(alpha)
        self._invalidate_java()

    def get_load_dependence(self) -> Optional[np.ndarray]:
        """Get load-dependent scaling."""
        return self._load_depend_scaling

    def set_class_dependence(self, beta) -> None:
        """
        Set class-dependent service rates.

        Args:
            beta: Scaling function or array indexed by class
        """
        self._class_depend_scaling = beta
        self._invalidate_java()

    def get_class_dependence(self):
        """Get class-dependent scaling."""
        return self._class_depend_scaling

    def set_joint_dependence(self, scaling_table: np.ndarray, cutoffs: List[int]) -> None:
        """
        Set joint load-class-dependent rates.

        Args:
            scaling_table: Joint scaling factors
            cutoffs: Population cutoff points
        """
        self._joint_depend_scaling = np.array(scaling_table)
        self._joint_cutoffs = cutoffs
        self._invalidate_java()

    def get_joint_dependence(self) -> Tuple[Optional[np.ndarray], Optional[List[int]]]:
        """Get joint scaling and cutoffs."""
        return self._joint_depend_scaling, getattr(self, '_joint_cutoffs', None)

    def set_drop_rule(self, rule: DropStrategy) -> None:
        """
        Set drop strategy for finite capacity.

        Args:
            rule: Drop strategy (DROP, BAS, or WaitingQueue)
        """
        self._drop_rule = rule
        self._invalidate_java()

    def get_drop_rule(self) -> DropStrategy:
        """Get drop strategy."""
        return self._drop_rule

    def is_station(self) -> bool:
        """Check if node is a station."""
        return True

    def get_station_index(self) -> int:
        """Get 1-based station index."""
        return self._station_index + 1 if self._station_index >= 0 else -1

    def get_station_index0(self) -> int:
        """Get 0-based station index."""
        return self._station_index

    def _set_station_index(self, idx: int) -> None:
        """Set 0-based station index (internal)."""
        self._station_index = idx
        self._invalidate_java()

    def getNumberOfServers(self) -> int:
        """Get number of servers ."""
        return self._number_of_servers

    # CamelCase aliases
    setClassCapacity = set_class_capacity
    getClassCapacity = get_class_capacity
    setLoadDependence = set_load_dependence
    getLoadDependence = get_load_dependence
    setClassDependence = set_class_dependence
    getClassDependence = get_class_dependence
    setJointDependence = set_joint_dependence
    getJointDependence = get_joint_dependence
    setDropRule = set_drop_rule
    getDropRule = get_drop_rule
    getStationIndex = get_station_index

    # Snake_case aliases
    get_number_of_servers = getNumberOfServers


# Forward declarations for type hints
class Network:
    """Placeholder for Network class (defined in network.py)."""
    pass


__all__ = [
    'ElementType', 'NodeType', 'SchedStrategy', 'RoutingStrategy',
    'JobClassType', 'SchedStrategyType', 'JoinStrategy', 'DropStrategy',
    'Element', 'NetworkElement', 'JobClass', 'Node', 'StatefulNode', 'Station'
]
