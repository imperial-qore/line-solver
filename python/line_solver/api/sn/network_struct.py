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
    FCFSPR = 22   # FCFS with Preemption (for compatibility)
    EDF = 23      # Earliest Deadline First
    FORK = 24     # Fork node
    JOIN = 25     # Join node
    REF = 26      # Reference task
    EDD = 27      # Earliest Due Date
    SRPT = 28     # Shortest Remaining Processing Time
    SRPTPRIO = 29 # SRPT with Priorities
    LCFSPRIO = 30 # LCFS with Priorities
    LCFSPRPRIO = 31  # LCFSPR with Priorities
    LCFSPIPRIO = 32  # LCFSPI with Priorities
    FCFSPRPRIO = 33  # FCFSPR with Priorities
    FCFSPIPRIO = 34  # FCFSPI with Priorities
    FCFSPRIO = 35 # FCFS with Priorities
    FSP = 36      # Fair Sojourn Protocol (virtual-PS finish time ranking)
    PAS = 37      # Pass-and-swap (order-independent queue with class swap graph)
    OI = 38       # Order-independent (pass-and-swap specialization with empty/zero swap graph)


class RoutingStrategy(IntEnum):
    """Routing strategies.

    Values must match MATLAB's RoutingStrategy constants for JMT compatibility.
    """
    RAND = 0      # Random (uniform among destinations)
    PROB = 1      # Probabilistic (explicit probabilities)
    RROBIN = 2    # Round Robin
    WRROBIN = 3   # Weighted Round Robin
    JSQ = 4       # Join Shortest Queue
    FIRING = 5    # Firing (for Petri nets)
    SQ = 6  # K-Choices
    RL = 7        # Reinforcement Learning
    DISABLED = -1 # Disabled routing


class DropStrategy(IntEnum):
    """Drop strategies for finite capacity.

    Values match the MATLAB DropStrategy constants and the JAR
    jline.lang.constant.DropStrategy ids, which are the interchange encoding of
    sn.droprule and sn.regionrule. Keep the three Python definitions of this
    enum (here, lang/base.py, constants.py) numerically identical: they are
    written and read by different modules over the same sn fields.
    """
    WAITQ = -1  # Wait in queue; also the default when no rule is set
    DROP = 1    # Drop on full
    BAS = 2     # Blocking after service
    BBS = 3     # Blocking before service
    RSRD = 4    # Re-service on rejection
    RETRIAL = 5             # Job moves to orbit and retries after delay
    RETRIAL_WITH_LIMIT = 6  # Job retries up to max attempts, then drops


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
        isfunction: (M, 1) Boolean mask, STATION-indexed: queue stations that
            carry setup/delay-off times (function stations)
        nodeToStation: (N, 1) Node index -> station index mapping
        nodeToStateful: (N, 1) Node index -> stateful index mapping
        stationToNode: (M, 1) Station index -> node index mapping
        stationToStateful: (M, 1) Station index -> stateful index mapping
        statefulToNode: (S, 1) Stateful index -> node index mapping
        statefulToStation: (S, 1) Stateful index -> station index mapping

        state: Dict State per stateful node
        lldscaling: (M, Nmax) Load-dependent scaling matrix
        cdscaling: Class-dependent (product-form) scaling functions beta_{i,r}
        jdscaling: Joint-dependent (non-product-form) scaling functions eta_i
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
    isfunction: np.ndarray = field(default_factory=lambda: np.array([]))

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
    cdscaling: Optional[Dict] = None  # class-dependent (product-form) scaling functions beta_{i,r}
    cdscalingpeak: Optional[np.ndarray] = None  # (M, R) declared peak class-dependent rate scaling, for Util=T*S/peak
    jdscaling: Optional[Dict] = None  # joint-dependent (non-product-form) scaling functions eta_i
    jdscalingpeak: Optional[np.ndarray] = None  # (M, R) declared peak joint-dependent rate scaling, for Util=T*S/peak

    # Class properties
    classprio: Optional[np.ndarray] = None  # (1, K) - class priorities
    classdeadline: Optional[np.ndarray] = None  # (1, K) - class deadlines
    isslc: Optional[np.ndarray] = None  # (K, 1) - is self-looping class
    issignal: Optional[np.ndarray] = None  # (K, 1) - is signal class
    signaltarget: Optional[np.ndarray] = None  # (K,) 0-based target class index a signal removes; -1 if none
    signaltype: Optional[List] = None  # signal types
    syncreply: Optional[np.ndarray] = None  # (K, 1) - sync reply class mapping
    immfeed: Optional[np.ndarray] = None  # (M, K) - immediate feedback enabled per station/class
    signalremdist: Optional[List] = None  # removal distribution per signal class
    signalrempolicy: Optional[List] = None  # removal policy per signal class
    iscatastrophe: Optional[np.ndarray] = None  # (K, 1) - is catastrophe signal class

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
    fjsync: Optional[List] = None  # fork firing synchronizations (FJ-augmented structs only, see ModelAdapter.fjtag)
    fjclassmap: Optional[np.ndarray] = None  # (nclasses,) original class of each FJ auxiliary class, 0 for originals
    isfjaugmented: bool = False  # True on FJ tag-augmented structs (Join/Fork carry count-vector states)

    # Drop rules
    droprule: Optional[Dict] = None  # station -> class -> drop strategy

    # Finite capacity regions
    nregions: int = 0
    region: Optional[List] = None
    regionrule: Optional[np.ndarray] = None
    regionweight: Optional[np.ndarray] = None  # Matrix(F, K) - class weights per region
    regionsz: Optional[np.ndarray] = None  # Matrix(F, K) - class sizes per region
    regionmaxmem: Optional[List] = None  # List of F Matrix(M,1): region global memory budget replicated on member station rows, -1 = unbounded
    # see _kb/04-networkstruct.md for rationale (regionmembers field semantics)
    regionmembers: Optional[List] = None

    # Synchronization
    sync: Optional[Dict] = None
    gsync: Optional[Dict] = None

    # Node parameters
    nodeparam: Optional[Dict] = None

    # Routing weights for WRROBIN strategy
    # Dict mapping (node_idx, class_idx) -> {dest_node_idx: weight}
    routingweights: Optional[Dict] = None

    # Reward functions
    reward: Optional[Dict] = None

    # Original routing (for caching)
    rtorig: Optional[Dict] = None

    # Additional matrices
    csmask: Optional[np.ndarray] = None  # class-switching mask
    nvars: Optional[np.ndarray] = None  # number of variables
    isbasblocking: Optional[np.ndarray] = None  # (nnodes,) 1 iff node is the blocking side of a true-BAS relation (BUG-83)
    isbasdestination: Optional[np.ndarray] = None  # (nstations,nclasses) True iff refusing an arrival here must block an upstream BAS station

    # Balking and retrial (populated from node objects during refresh_struct)
    impatienceType: Optional[np.ndarray] = None      # (M, K) ProcessType id of patience dist (reneging); 0 = none
    impatienceMu: Optional[np.ndarray] = None         # (M, K) patience rate (1/mean) for reneging
    impatienceClass: Optional[np.ndarray] = None      # (M, K) ImpatienceType id (RENEGING/BALKING); 0 = none
    impatiencePhi: Optional[np.ndarray] = None        # (M, K) patience SCV
    impatiencePhases: Optional[np.ndarray] = None     # (M, K) phases of the patience distribution
    impatienceProc: Optional[List[List[Any]]] = None  # M x K list-of-lists; (D0, D1) tuple or None
    impatiencePie: Optional[List[List[Any]]] = None   # M x K list-of-lists; entry probability vector or None
    impatienceDist: Optional[List[List[Any]]] = None  # M x K list-of-lists; native Distribution objects (JMT export)
    balkingStrategy: Optional[np.ndarray] = None     # (M, K) BalkingStrategy id, 0 = none
    balkingThresholds: Optional[List[List[Any]]] = None  # M x K list-of-lists; entry is list of (minJobs,maxJobs,prob) tuples or None
    retrialType: Optional[np.ndarray] = None         # (M, K) ProcessType id, 0 = none (use retrialProc to disambiguate)
    retrialMu: Optional[np.ndarray] = None           # (M, K) retrial delay rate (1/mean)
    retrialPhi: Optional[np.ndarray] = None          # (M, K) retrial delay SCV
    retrialProc: Optional[List[List[Any]]] = None    # M x K list-of-lists; entry is (D0, D1) tuple or None
    retrialMaxAttempts: Optional[np.ndarray] = None  # (M, K) max retrial attempts (-1 = unlimited)
    retrialPolicy: Optional[np.ndarray] = None       # (M, K) RetrialPolicy: LINEAR = per-customer rate n*nu, CONSTANT = orbit-wide rate nu
    orbitMaxJobs: Optional[np.ndarray] = None        # (M, K) orbit capacity (-1 = unbounded); a job finding the orbit full is lost
    orbitImpatience: Optional[List[List[Any]]] = None  # M x K list-of-lists; orbit abandonment (D0, D1) tuple or None

    # Server breakdown / repair (populated from node objects during refresh_struct)
    hasbreakdown: Optional[np.ndarray] = None        # (nnodes,) 1 iff the node's server is subject to breakdowns
    breakdownMu: Optional[np.ndarray] = None         # (M,) failure rate of an up server (0 = never fails)
    repairMu: Optional[np.ndarray] = None            # (M,) repair rate of a down server (0 = never repaired)
    breakdownProc: Optional[List[Any]] = None        # (M,) list; failure-time (D0, D1) tuple or None
    repairProc: Optional[List[Any]] = None           # (M,) list; repair-time (D0, D1) tuple or None
    downServiceRates: Optional[np.ndarray] = None    # (M, K) service rate while the server is down (0 = no service)

    # Cache item state tracking (mirrors JAR varsparam)
    varsparam: Optional[np.ndarray] = None  # (nnodes, 1) item indices, -1 = none

    # Marked (MMAP) source arrivals (mirrors MATLAB/JAR sn.markidx)
    markidx: Optional[np.ndarray] = None  # (M, K) mark index (1-based) of class r at source station i, -1 = not marked

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

    def has_joint_dependence(self) -> bool:
        """Check if model has joint-dependent (non-product-form) scaling."""
        return self.jdscaling is not None and len(self.jdscaling) > 0

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
