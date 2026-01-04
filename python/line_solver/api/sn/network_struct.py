"""
NetworkStruct - Native Python implementation.

This dataclass summarizes the characteristics of a Network object,
providing all parameters needed for queueing network analysis.

Ported from MATLAB implementation.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple, Union
from enum import IntEnum


class MatrixArray(np.ndarray):
    """
    Numpy array subclass with .get() and .set() methods for API compatibility.

    This class provides compatibility with the wrapper mode that uses JLine's
    Matrix class which has get(i, j) and set(i, j, value) methods.
    """

    def __new__(cls, input_array):
        """Create MatrixArray from existing array."""
        obj = np.asarray(input_array).view(cls)
        return obj

    def __array_finalize__(self, obj):
        """Handle view casting and new-from-template."""
        pass

    def __getitem__(self, key):
        """
        Override indexing to handle 2D indexing on 1D arrays.

        This provides compatibility with MATLAB-style row/column vectors
        where a 1D array can be indexed as (0, j) or (i, 0).
        """
        # Handle 2D tuple indexing on 1D arrays
        if isinstance(key, tuple) and len(key) == 2 and self.ndim == 1:
            i, j = key
            # For row vector style (0, j) -> return element j
            if i == 0:
                return super().__getitem__(j)
            # For column vector style (i, 0) -> return element i
            elif j == 0:
                return super().__getitem__(i)
        return super().__getitem__(key)

    def __setitem__(self, key, value):
        """
        Override item setting to handle 2D indexing on 1D arrays.
        """
        # Handle 2D tuple indexing on 1D arrays
        if isinstance(key, tuple) and len(key) == 2 and self.ndim == 1:
            i, j = key
            # For row vector style (0, j) -> set element j
            if i == 0:
                super().__setitem__(j, value)
                return
            # For column vector style (i, 0) -> set element i
            elif j == 0:
                super().__setitem__(i, value)
                return
        super().__setitem__(key, value)

    def get(self, i, j=None):
        """
        Get element at index (i, j) or just i if 1D.

        Args:
            i: Row index (or element index for 1D)
            j: Column index (optional, for 2D arrays)

        Returns:
            Element value at the specified index
        """
        if j is None:
            return self[i]
        # For 1D arrays, handle MATLAB-style row/column vector indexing
        if self.ndim == 1:
            # Row vector style: (0, j) -> element j
            if i == 0:
                return self[j]
            # Column vector style: (i, 0) -> element i
            elif j == 0:
                return self[i]
            # Otherwise just return element at first index
            return self[i]
        return self[i, j]

    def set(self, i, j, value=None):
        """
        Set element at index (i, j) or just i if 1D.

        Args:
            i: Row index (or element index for 1D)
            j: Column index or value (for 1D arrays)
            value: Value to set (optional, for 2D arrays)
        """
        if value is None:
            # Called as set(i, value) for 1D
            self[i] = j
        else:
            # Called as set(i, j, value) for 2D
            self[i, j] = value


class NodeType(IntEnum):
    """Node types in a queueing network.

    NOTE: Values must match lang/base.py NodeType enum.
    """
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
    LOGGER = 11
    FINITE_CAPACITY_REGION = 12

    @staticmethod
    def toText(node_type: 'NodeType') -> str:
        """Convert node type to text representation."""
        names = {
            NodeType.SOURCE: 'Source',
            NodeType.SINK: 'Sink',
            NodeType.QUEUE: 'Queue',
            NodeType.DELAY: 'Delay',
            NodeType.FORK: 'Fork',
            NodeType.JOIN: 'Join',
            NodeType.CACHE: 'Cache',
            NodeType.ROUTER: 'Router',
            NodeType.CLASSSWITCH: 'ClassSwitch',
            NodeType.PLACE: 'Place',
            NodeType.TRANSITION: 'Transition',
            NodeType.LOGGER: 'Logger',
            NodeType.FINITE_CAPACITY_REGION: 'Region',
        }
        return names.get(node_type, f'Unknown({node_type})')


class SchedStrategy(IntEnum):
    """Scheduling strategies."""
    FCFS = 0    # First-Come First-Served
    LCFS = 1    # Last-Come First-Served
    LCFSPR = 2  # LCFS Preemptive Resume
    LCFSPI = 3  # LCFS Preemptive Identical
    PS = 4      # Processor Sharing
    DPS = 5     # Discriminatory PS
    GPS = 6     # Generalized PS
    INF = 7     # Infinite Server (Delay)
    RAND = 8    # Random
    HOL = 9     # Head of Line
    SEPT = 10   # Shortest Expected Processing Time
    LEPT = 11   # Longest Expected Processing Time
    SIRO = 12   # Service in Random Order
    SJF = 13    # Shortest Job First
    LJF = 14    # Longest Job First
    POLLING = 15
    EXT = 16    # External
    LPS = 17    # Least Progress Scheduling
    SETF = 18   # Shortest Elapsed Time First
    DPSPRIO = 19  # DPS with Priority
    GPSPRIO = 20  # GPS with Priority
    PSPRIO = 21   # PS with Priority


class RoutingStrategy(IntEnum):
    """Routing strategies."""
    DISABLED = 0
    RAND = 1    # Random
    PROB = 2    # Probabilistic
    RROBIN = 3  # Round Robin
    WRROBIN = 4 # Weighted Round Robin
    JSQ = 5     # Join Shortest Queue
    RL = 6      # Reinforcement Learning
    KCHOICES = 7 # K-Choices
    CLASS_SWITCH = 8


class DropStrategy(IntEnum):
    """Drop strategies for finite capacity."""
    WAITQ = 0   # Wait in queue
    DROP = 1    # Drop on full
    BAS = 2     # Blocking after service


@dataclass
class NetworkStruct:
    """
    Data structure summarizing network characteristics.

    This class is the Python equivalent in native Python.
    It contains all parameters needed by solvers to analyze a queueing network.

    Attributes:
        nstations: Number of stations (queues, delays, sources, joins, places)
        nstateful: Number of stateful nodes
        nnodes: Total number of nodes
        nclasses: Number of job classes
        nchains: Number of chains (routing chains)
        nclosedjobs: Total number of jobs in closed classes

        njobs: (1, K) Population per class (inf for open classes)
        nservers: (M, 1) Number of servers per station
        rates: (M, K) Service rates
        scv: (M, K) Squared coefficient of variation

        visits: Dict[int, ndarray] - Chain ID -> (M, K) visit ratios
        inchain: Dict[int, ndarray] - Chain ID -> class indices in chain
        chains: (K, 1) Chain membership per class
        refstat: (K, 1) Reference station per class
        refclass: (1, C) Reference class per chain

        sched: Dict[int, SchedStrategy] - Station ID -> scheduling strategy
        routing: (N, K) Routing strategy matrix
        rt: Routing probability matrix

        nodetype: List[NodeType] - Node types
        isstation: (N, 1) Boolean mask for stations
        isstateful: (N, 1) Boolean mask for stateful nodes
        nodeToStation: (N, 1) Node index -> station index mapping
        nodeToStateful: (N, 1) Node index -> stateful index mapping
        stationToNode: (M, 1) Station index -> node index mapping
        stationToStateful: (M, 1) Station index -> stateful index mapping
        statefulToNode: (S, 1) Stateful index -> node index mapping
        statefulToStation: (S, 1) Stateful index -> station index mapping

        state: Dict State per stateful node
        lldscaling: (M, Nmax) Load-dependent scaling matrix
        cdscaling: Class-dependent scaling functions
        cap: (M, 1) Station capacities
        classcap: (M, K) Per-class capacities

        connmatrix: (N, N) Connection matrix
        nodenames: List[str] - Node names
        classnames: List[str] - Class names
    """

    # Dimensions
    nstations: int = 0
    nstateful: int = 0
    nnodes: int = 0
    nclasses: int = 0
    nchains: int = 0
    nclosedjobs: int = 0

    # Population and capacity
    njobs: np.ndarray = field(default_factory=lambda: np.array([]))
    nservers: np.ndarray = field(default_factory=lambda: np.array([]))
    cap: Optional[np.ndarray] = None
    classcap: Optional[np.ndarray] = None

    # Service parameters
    rates: np.ndarray = field(default_factory=lambda: np.array([]))
    scv: np.ndarray = field(default_factory=lambda: np.array([]))
    phases: Optional[np.ndarray] = None
    phasessz: Optional[np.ndarray] = None
    phaseshift: Optional[np.ndarray] = None

    # Network structure - chains
    visits: Dict[int, np.ndarray] = field(default_factory=dict)
    nodevisits: Dict[int, np.ndarray] = field(default_factory=dict)
    inchain: Dict[int, np.ndarray] = field(default_factory=dict)
    chains: np.ndarray = field(default_factory=lambda: np.array([]))
    refstat: np.ndarray = field(default_factory=lambda: np.array([]))
    refclass: np.ndarray = field(default_factory=lambda: np.array([]))

    # Scheduling and routing
    sched: Dict[int, int] = field(default_factory=dict)  # station_id -> SchedStrategy
    schedparam: Optional[np.ndarray] = None
    routing: np.ndarray = field(default_factory=lambda: np.array([]))
    rt: Optional[np.ndarray] = None
    rtnodes: Optional[np.ndarray] = None

    # Node classification
    nodetype: List[int] = field(default_factory=list)
    isstation: np.ndarray = field(default_factory=lambda: np.array([]))
    isstateful: np.ndarray = field(default_factory=lambda: np.array([]))
    isstatedep: Optional[np.ndarray] = None  # (N, 3) - buffer, srv, routing

    # Node mappings (0-indexed)
    nodeToStation: np.ndarray = field(default_factory=lambda: np.array([]))
    nodeToStateful: np.ndarray = field(default_factory=lambda: np.array([]))
    stationToNode: np.ndarray = field(default_factory=lambda: np.array([]))
    stationToStateful: np.ndarray = field(default_factory=lambda: np.array([]))
    statefulToNode: np.ndarray = field(default_factory=lambda: np.array([]))
    statefulToStation: np.ndarray = field(default_factory=lambda: np.array([]))

    # State information
    state: Dict[int, np.ndarray] = field(default_factory=dict)
    stateprior: Dict[int, np.ndarray] = field(default_factory=dict)
    space: Dict[int, np.ndarray] = field(default_factory=dict)

    # Load-dependent scaling
    lldscaling: Optional[np.ndarray] = None  # (M, Nmax) - limited load dependence
    cdscaling: Optional[Dict] = None  # class-dependent scaling functions

    # Class properties
    classprio: Optional[np.ndarray] = None  # (1, K) - class priorities
    classdeadline: Optional[np.ndarray] = None  # (1, K) - class deadlines
    isslc: Optional[np.ndarray] = None  # (K, 1) - is self-looping class
    issignal: Optional[np.ndarray] = None  # (K, 1) - is signal class
    signaltype: Optional[List] = None  # signal types
    syncreply: Optional[np.ndarray] = None  # (K, 1) - sync reply class mapping

    # Connectivity
    connmatrix: Optional[np.ndarray] = None  # (N, N) connection matrix

    # Naming
    nodenames: List[str] = field(default_factory=list)
    classnames: List[str] = field(default_factory=list)

    # Process information
    mu: Optional[Dict] = None  # station -> class -> service rate matrix
    phi: Optional[Dict] = None  # station -> class -> service phase probs
    proc: Optional[Dict] = None  # station -> class -> process cell
    pie: Optional[Dict] = None  # station -> class -> initial phase probs
    procid: Optional[Dict] = None  # station -> class -> process type
    lst: Optional[Dict] = None  # station -> class -> Laplace-Stieltjes transform

    # Fork-join
    fj: Optional[np.ndarray] = None  # fork-join topology matrix

    # Drop rules
    droprule: Optional[Dict] = None  # station -> class -> drop strategy

    # Finite capacity regions
    nregions: int = 0
    region: Optional[List] = None
    regionrule: Optional[np.ndarray] = None
    regionweight: Optional[np.ndarray] = None  # Matrix(F, K) - class weights per region
    regionsz: Optional[np.ndarray] = None  # Matrix(F, K) - class sizes per region

    # Synchronization
    sync: Optional[Dict] = None
    gsync: Optional[Dict] = None

    # Node parameters
    nodeparam: Optional[Dict] = None

    # Reward functions
    reward: Optional[Dict] = None

    # Original routing (for caching)
    rtorig: Optional[Dict] = None

    # Additional matrices
    csmask: Optional[np.ndarray] = None  # class-switching mask
    nvars: Optional[np.ndarray] = None  # number of variables

    # Fields that should be wrapped with MatrixArray for API compatibility
    _MATRIX_ARRAY_FIELDS = frozenset([
        'njobs', 'nservers', 'rates', 'scv', 'chains', 'refstat', 'refclass',
        'routing', 'isstation', 'isstateful', 'nodeToStation', 'nodeToStateful',
        'stationToNode', 'stationToStateful', 'statefulToNode', 'statefulToStation',
        'classprio', 'rt', 'rtnodes'
    ])

    def __setattr__(self, name: str, value):
        """Override to convert numpy arrays to MatrixArray for API compatibility."""
        # Convert numpy arrays to MatrixArray for specified fields
        if name in NetworkStruct._MATRIX_ARRAY_FIELDS:
            if value is not None and not isinstance(value, MatrixArray):
                if isinstance(value, np.ndarray) or isinstance(value, (list, tuple)):
                    value = MatrixArray(value)
        object.__setattr__(self, name, value)

    def __post_init__(self):
        """Ensure arrays are MatrixArray (numpy arrays with .get()/.set() methods)."""
        for fname in NetworkStruct._MATRIX_ARRAY_FIELDS:
            val = getattr(self, fname)
            if val is not None:
                # Convert to MatrixArray for API compatibility (.get()/.set() methods)
                if not isinstance(val, MatrixArray):
                    object.__setattr__(self, fname, MatrixArray(val))

    def validate(self) -> None:
        """
        Validate structural consistency.

        Raises:
            ValueError: If structural consistency is violated
        """
        # Check basic counts
        if self.nstations < 0 or self.nstateful < 0 or self.nnodes < 0:
            raise ValueError("Node counts must be non-negative")

        if self.nstations > self.nstateful:
            raise ValueError("Number of stations cannot exceed number of stateful nodes")

        if self.nstateful > self.nnodes:
            raise ValueError("Number of stateful nodes cannot exceed total number of nodes")

        # Validate matrix dimensions
        if self.isstation is not None and len(self.isstation) > 0:
            if len(self.isstation) != self.nnodes:
                raise ValueError(f"isstation length {len(self.isstation)} != nnodes {self.nnodes}")
            if np.sum(self.isstation) != self.nstations:
                raise ValueError("nstations must equal sum of isstation")

        if self.isstateful is not None and len(self.isstateful) > 0:
            if len(self.isstateful) != self.nnodes:
                raise ValueError(f"isstateful length {len(self.isstateful)} != nnodes {self.nnodes}")
            if np.sum(self.isstateful) != self.nstateful:
                raise ValueError("nstateful must equal sum of isstateful")

        # Validate hierarchy: all stations must be stateful
        if (self.isstation is not None and self.isstateful is not None and
                len(self.isstation) > 0 and len(self.isstateful) > 0):
            for i in range(self.nnodes):
                if self.isstation[i] > 0 and self.isstateful[i] == 0:
                    raise ValueError(f"All stations must be stateful nodes (violation at node {i})")

        # Validate rates matrix dimensions
        if self.rates is not None and len(self.rates) > 0:
            if self.rates.ndim == 2:
                if self.rates.shape[0] != self.nstations:
                    raise ValueError(f"rates rows {self.rates.shape[0]} != nstations {self.nstations}")
                if self.rates.shape[1] != self.nclasses:
                    raise ValueError(f"rates cols {self.rates.shape[1]} != nclasses {self.nclasses}")

    def is_valid(self) -> bool:
        """
        Check if structure is valid.

        Returns:
            True if structure passes validation, False otherwise
        """
        try:
            self.validate()
            return True
        except ValueError:
            return False

    def get_chain_population(self, chain_id: int) -> float:
        """
        Get total population in a chain.

        Args:
            chain_id: Chain index (0-based)

        Returns:
            Total number of jobs in the chain
        """
        if chain_id not in self.inchain:
            return 0.0
        class_indices = self.inchain[chain_id]
        return np.sum(self.njobs.flat[class_indices.astype(int)])

    def is_closed_chain(self, chain_id: int) -> bool:
        """
        Check if a chain is closed (finite population).

        Args:
            chain_id: Chain index (0-based)

        Returns:
            True if chain is closed, False if open
        """
        if chain_id not in self.inchain:
            return False
        class_indices = self.inchain[chain_id]
        # Chain is closed if all classes have finite population
        return np.all(np.isfinite(self.njobs.flat[class_indices.astype(int)]))

    def is_open_chain(self, chain_id: int) -> bool:
        """
        Check if a chain is open (infinite population).

        Args:
            chain_id: Chain index (0-based)

        Returns:
            True if chain is open, False if closed
        """
        if chain_id not in self.inchain:
            return False
        class_indices = self.inchain[chain_id]
        # Chain is open if any class has infinite population
        return np.any(np.isinf(self.njobs.flat[class_indices.astype(int)]))

    def get_station_indices(self) -> np.ndarray:
        """
        Get indices of station nodes.

        Returns:
            Array of node indices that are stations
        """
        return np.where(self.isstation > 0)[0]

    def get_stateful_indices(self) -> np.ndarray:
        """
        Get indices of stateful nodes.

        Returns:
            Array of node indices that are stateful
        """
        return np.where(self.isstateful > 0)[0]

    def get_scheduling_at_station(self, station_id: int) -> int:
        """
        Get scheduling strategy at a station.

        Args:
            station_id: Station index (0-based)

        Returns:
            SchedStrategy value
        """
        return self.sched.get(station_id, SchedStrategy.FCFS)

    def has_multi_server(self) -> bool:
        """Check if any station has multiple servers."""
        if self.nservers is None or len(self.nservers) == 0:
            return False
        return np.any(self.nservers > 1)

    def has_load_dependence(self) -> bool:
        """Check if model has load-dependent service rates."""
        return self.lldscaling is not None and np.any(self.lldscaling != 1.0)

    def has_class_dependence(self) -> bool:
        """Check if model has class-dependent scaling."""
        return self.cdscaling is not None and len(self.cdscaling) > 0

    def has_open_classes(self) -> bool:
        """Check if model has open (infinite population) classes."""
        return np.any(np.isinf(self.njobs))

    def has_closed_classes(self) -> bool:
        """Check if model has closed (finite population) classes."""
        return np.any(np.isfinite(self.njobs) & (self.njobs > 0))

    def get_total_population(self) -> float:
        """Get total population across all closed classes."""
        return np.sum(self.njobs[np.isfinite(self.njobs)])

    def get_open_class_indices(self) -> np.ndarray:
        """Get indices of open classes."""
        return np.where(np.isinf(self.njobs.flatten()))[0]

    def get_closed_class_indices(self) -> np.ndarray:
        """Get indices of closed classes."""
        return np.where(np.isfinite(self.njobs.flatten()))[0]

    def copy(self) -> 'NetworkStruct':
        """Create a deep copy of this NetworkStruct."""
        import copy
        return copy.deepcopy(self)

    def __repr__(self) -> str:
        """String representation."""
        return (f"NetworkStruct(nstations={self.nstations}, nstateful={self.nstateful}, "
                f"nnodes={self.nnodes}, nclasses={self.nclasses}, nchains={self.nchains}, "
                f"nclosedjobs={self.nclosedjobs})")

    @property
    def obj(self):
        """Return self for compatibility with wrapper code that accesses .obj"""
        return self
