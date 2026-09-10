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
    LOGGER = 11


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
    # The size-based disciplines were missing HERE while constants.SchedStrategy
    # carried them, so _normalize_sched_strategy -- which maps BY NAME -- refused
    # a PSJF/FB/LRPT model outright. Appended at free ids rather than at the
    # constants values 36/37/38/39, which are already FSP/PAS/OI here: the two
    # enums are reconciled by NAME, never by value, and nothing indexes this one
    # positionally.
    PSJF = 39     # Preemptive Shortest Job First
    FB = 40       # Foreground-Background
    LAS = 41      # Least Attained Service (the FB alias)
    LRPT = 42     # Longest Remaining Processing Time
    # appended, not slotted next to FCFSPR: 43 is the first id free in all
    # three Python SchedStrategy enums, and constants.toID is positional
    FCFSPI = 43   # FCFS Preemptive Identical

    @staticmethod
    def to_feature(strategy: 'SchedStrategy') -> str:
        """Convert a scheduling strategy to its feature name for solver support checking."""
        # see _kb/01-model-classes.md (Enumerations) for rationale
        return _SCHED_FEATURE_BY_NAME.get(_strategy_name(strategy), '')


def _strategy_name(strategy) -> str:
    """Name of an enum-like strategy value, for cross-class feature lookup."""
    name = getattr(strategy, 'name', None)
    if name is None:
        name = str(strategy).split('.')[-1]
    return name


# The preempt family, whose station buffer holds [class, phase] pairs rather
# than bare class tags. It factors as {FCFS,LCFS} x {PR,PI} x {plain,PRIO}, and
# the subsets below are those factors. Any site that enumerates the family must
# cover all eight: listing only some is how the PI disciplines were left without
# a departure arm (MATLAB State.toMarginal:193, fromMarginal:245 list all).
SCHED_PREEMPT_LCFS = (SchedStrategy.LCFSPR, SchedStrategy.LCFSPRPRIO,
                      SchedStrategy.LCFSPI, SchedStrategy.LCFSPIPRIO)
SCHED_PREEMPT_FCFS = (SchedStrategy.FCFSPR, SchedStrategy.FCFSPRPRIO,
                      SchedStrategy.FCFSPI, SchedStrategy.FCFSPIPRIO)
SCHED_PREEMPT_PRIO = (SchedStrategy.LCFSPRPRIO, SchedStrategy.LCFSPIPRIO,
                      SchedStrategy.FCFSPRPRIO, SchedStrategy.FCFSPIPRIO)
# PR resumes at the stored phase, PI discards it and restarts from pie
SCHED_PREEMPT_RESUME = (SchedStrategy.LCFSPR, SchedStrategy.LCFSPRPRIO,
                        SchedStrategy.FCFSPR, SchedStrategy.FCFSPRPRIO)
SCHED_PREEMPT = SCHED_PREEMPT_LCFS + SCHED_PREEMPT_FCFS

# Buffer ENCODING groups, the counterpart of the C++ state_detail::buffer_is_*
# predicates. Every encode and decode of a station's local state must branch on
# these rather than on an inline list: five functions in api/state/marginal.py
# each repeated their own list, and all five had drifted to the same omissions,
# so LCFSPRIO, SRPT and SRPTPRIO were counted server-only.
# POLLING is deliberately ABSENT. It carries a per-class count in most callers
# but has its own branch in fromMarginal and none at all in
# fromMarginalAndRunning, so folding it in here would silently change what that
# function does. It stays named explicitly at the sites that already handle it.
SCHED_BUFFER_CLASS_TAG = (SchedStrategy.FCFS, SchedStrategy.FCFSPRIO,
                          SchedStrategy.HOL, SchedStrategy.LCFS,
                          SchedStrategy.LCFSPRIO)
SCHED_BUFFER_PER_CLASS_COUNT = (SchedStrategy.SIRO, SchedStrategy.SEPT,
                                SchedStrategy.LEPT, SchedStrategy.SRPT,
                                SchedStrategy.SRPTPRIO)

# see _kb/01-model-classes.md (Enumerations) for rationale
_SCHED_FEATURE_BY_NAME = {
    'INF': 'SchedStrategy_INF',
    'FCFS': 'SchedStrategy_FCFS',
    'FCFSPR': 'SchedStrategy_FCFSPR',
    'FCFSPI': 'SchedStrategy_FCFSPI',
    'LCFS': 'SchedStrategy_LCFS',
    'LCFSPR': 'SchedStrategy_LCFSPR',
    'LCFSPI': 'SchedStrategy_LCFSPI',
    'POLLING': 'SchedStrategy_POLLING',
    'RAND': 'SchedStrategy_SIRO',
    'SIRO': 'SchedStrategy_SIRO',
    'SJF': 'SchedStrategy_SJF',
    'LJF': 'SchedStrategy_LJF',
    'PS': 'SchedStrategy_PS',
    'DPS': 'SchedStrategy_DPS',
    'GPS': 'SchedStrategy_GPS',
    'SEPT': 'SchedStrategy_SEPT',
    'LEPT': 'SchedStrategy_LEPT',
    'SRPT': 'SchedStrategy_SRPT',
    'SRPTPRIO': 'SchedStrategy_SRPTPRIO',
    'HOL': 'SchedStrategy_HOL',
    # see _kb/01-model-classes.md (Enumerations) for rationale
    'FCFSPRIO': 'SchedStrategy_HOL',
    'EXT': 'SchedStrategy_EXT',
    'REF': 'SchedStrategy_REF',
    'PSPRIO': 'SchedStrategy_PSPRIO',
    'DPSPRIO': 'SchedStrategy_DPSPRIO',
    'GPSPRIO': 'SchedStrategy_GPSPRIO',
    'LCFSPRIO': 'SchedStrategy_LCFSPRIO',
    'LCFSPRPRIO': 'SchedStrategy_LCFSPRPRIO',
    'LCFSPIPRIO': 'SchedStrategy_LCFSPIPRIO',
    'FCFSPRPRIO': 'SchedStrategy_FCFSPRPRIO',
    'FCFSPIPRIO': 'SchedStrategy_FCFSPIPRIO',
    'EDD': 'SchedStrategy_EDD',
    'EDF': 'SchedStrategy_EDF',
    'LPS': 'SchedStrategy_LPS',
    'PSJF': 'SchedStrategy_PSJF',
    'FB': 'SchedStrategy_FB',
    # Least attained service is feedback scheduling; MATLAB defines LAS as an
    # alias of FB, native python keeps a separate member.
    'LAS': 'SchedStrategy_FB',
    'LRPT': 'SchedStrategy_LRPT',
    'SETF': 'SchedStrategy_SETF',
    'SET': 'SchedStrategy_SETF',
    'FSP': 'SchedStrategy_FSP',
    'PAS': 'SchedStrategy_PAS',
    'OI': 'SchedStrategy_OI',
}


class RoutingStrategy(IntEnum):
    """Enumeration of routing strategies.

    Values must match MATLAB's RoutingStrategy constants for JMT compatibility.
    """
    RAND = 0      # Random routing (uniform among destinations)
    PROB = 1      # Probabilistic routing (explicit probabilities)
    RROBIN = 2    # Round-robin
    WRROBIN = 3   # Weighted round-robin
    JSQ = 4       # Join-Shortest-Queue
    FIRING = 5    # Firing (for Petri nets)
    SQ = 6  # K-choices policy
    SDR = 7       # Krzesinski (1987) product-form state-dependent routing
    DISABLED = -1 # Disabled routing

    @staticmethod
    def to_feature(strategy: 'RoutingStrategy') -> str:
        """Convert a routing strategy to its feature name for solver support checking."""
        # see _kb/01-model-classes.md (Enumerations) for rationale
        return _ROUTING_FEATURE_BY_NAME.get(_strategy_name(strategy), '')


# see _kb/01-model-classes.md (Enumerations) for rationale
_ROUTING_FEATURE_BY_NAME = {
    'RAND': 'RoutingStrategy_RAND',
    'PROB': 'RoutingStrategy_PROB',
    'RROBIN': 'RoutingStrategy_RROBIN',
    'WRROBIN': 'RoutingStrategy_WRROBIN',
    'JSQ': 'RoutingStrategy_JSQ',
    'SQ': 'RoutingStrategy_SQ',
    'SDR': 'RoutingStrategy_SDR',
}


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
    """Enumeration of join strategies.

    PARTIAL is the canonical name of the k-of-n join, as in the MATLAB
    JoinStrategy class, so that `.name` serializes straight to the interchange
    spelling; QUORUM is kept as an alias for it because that is how the JAR
    (jline.lang.constant.JoinStrategy.Quorum) spells the same strategy. Only
    the names cross codebases, in the JSON `joinStrategy` field and the .lqnx
    path; the numeric values are compared solely against this enum, so they
    need not match the MATLAB or JAR ids.
    """
    STD = 0        # Standard (AND-join)
    PARTIAL = 1    # k-of-n join (MATLAB spelling, canonical)
    QUORUM = 1     # Alias of PARTIAL (JAR spelling)
    CANDJOIN = 2   # Cache AND-join


class DropStrategy(IntEnum):
    """Enumeration of drop strategies for finite capacity.

    Values match the MATLAB DropStrategy constants and the JAR
    jline.lang.constant.DropStrategy ids, which are the interchange encoding of
    sn.droprule and sn.regionrule. Keep the three Python definitions of this
    enum (here, api/sn/network_struct.py, constants.py) numerically identical:
    they are written and read by different modules over the same sn fields.
    """
    WAITQ = -1    # Wait in queue; also the default when no rule is set
    WaitingQueue = -1  # Alias for WAITQ
    Queue = -1    # Alias for WAITQ
    DROP = 1      # Drop arriving job
    BAS = 2       # Block after service
    BBS = 3       # Block before service
    RSRD = 4      # Re-service on rejection
    RETRIAL = 5           # Unlimited retries
    RETRIAL_WITH_LIMIT = 6  # Retries up to max attempts


class BalkingStrategy(IntEnum):
    """Enumeration of balking strategies."""
    QUEUE_LENGTH = 1   # Balk based on queue length ranges
    EXPECTED_WAIT = 2  # Balk based on expected waiting time
    COMBINED = 3       # Both conditions (OR logic)

    @staticmethod
    def to_text(strategy):
        if strategy == BalkingStrategy.QUEUE_LENGTH:
            return "QUEUE_LENGTH"
        elif strategy == BalkingStrategy.EXPECTED_WAIT:
            return "EXPECTED_WAIT"
        elif strategy == BalkingStrategy.COMBINED:
            return "COMBINED"
        return strategy.name

    @classmethod
    def from_text(cls, text):
        text = text.upper().replace(" ", "_")
        return cls[text]


class RetrialPolicy(IntEnum):
    """Enumeration of retrial policies for an orbiting population.

    The policy fixes how the aggregate rate at which the orbit attempts to
    re-enter the station depends on the orbit size n:

      LINEAR   - each orbiting customer retries independently at rate nu, so
                 the aggregate retrial rate is n*nu. This is the classical
                 retrial queue of Falin and Templeton.
      CONSTANT - the orbit as a whole retries at rate nu whenever it is
                 non-empty, independently of n. This models a single retrial
                 controller shared by the orbit rather than per-customer
                 timers.

    Values match MATLAB RetrialPolicy.
    """
    LINEAR = 1
    CONSTANT = 2

    @staticmethod
    def to_text(policy):
        if policy == RetrialPolicy.LINEAR:
            return 'linear'
        elif policy == RetrialPolicy.CONSTANT:
            return 'constant'
        raise ValueError('Unrecognized retrial policy type.')

    @classmethod
    def from_text(cls, text):
        key = str(text).lower()
        if key in ('linear', 'per-customer', 'falin'):
            return cls.LINEAR
        if key in ('constant', 'fixed', 'controller'):
            return cls.CONSTANT
        raise ValueError("Unrecognized retrial policy '%s'." % text)

    toText = to_text
    fromText = from_text


class ImpatienceType(IntEnum):
    """Enumeration of customer impatience types."""
    RENEGING = 1  # Customer abandons after joining (timer-based)
    BALKING = 2   # Customer refuses to join based on queue state
    RETRIAL = 3   # Customer moves to orbit and retries after delay

    @staticmethod
    def to_text(imp_type):
        if imp_type == ImpatienceType.RENEGING:
            return "reneging"
        elif imp_type == ImpatienceType.BALKING:
            return "balking"
        elif imp_type == ImpatienceType.RETRIAL:
            return "retrial"
        return imp_type.name.lower()

    @classmethod
    def from_text(cls, text):
        text = text.upper()
        return cls[text]


class ReplacementStrategy(IntEnum):
    """
    Cache replacement strategies.

    Determines which item to evict when the cache is full and a new item
    needs to be stored.

    Attributes:
        RR: Random Replacement - evict a random item
        FIFO: First-In-First-Out - evict oldest item
        SFIFO: Strict FIFO - FIFO eviction, no reinsertion/promotion on hit
        LRU: Least Recently Used - evict least recently accessed item
    """
    RR = 0      # Random Replacement
    FIFO = 1    # First-In-First-Out
    SFIFO = 2   # Strict FIFO
    LRU = 3     # Least Recently Used
    HLRU = 4    # hierarchical/k-LRU: h lists, LRU discipline, promote i->i+1 on hit
    CLIMB = 5   # move-up-one-position on hit (transposition rule)
    QLRU = 6    # q-LRU: LRU discipline with probabilistic admission q on a miss

    @staticmethod
    def to_string(strategy: 'ReplacementStrategy') -> str:
        """Convert a replacement strategy to its string representation."""
        if strategy == ReplacementStrategy.RR:
            return 'rr'
        elif strategy == ReplacementStrategy.FIFO:
            return 'fifo'
        elif strategy == ReplacementStrategy.SFIFO:
            return 'sfifo'
        elif strategy == ReplacementStrategy.LRU:
            return 'lru'
        elif strategy == ReplacementStrategy.HLRU:
            return 'hlru'
        elif strategy == ReplacementStrategy.CLIMB:
            return 'climb'
        elif strategy == ReplacementStrategy.QLRU:
            return 'qlru'
        else:
            return str(strategy.name).lower()

    @staticmethod
    def to_feature(strategy: 'ReplacementStrategy') -> str:
        """Convert a replacement strategy to its feature name for solver support checking."""
        if strategy == ReplacementStrategy.RR:
            return 'ReplacementStrategy_RR'
        elif strategy == ReplacementStrategy.FIFO:
            return 'ReplacementStrategy_FIFO'
        elif strategy == ReplacementStrategy.SFIFO:
            return 'ReplacementStrategy_SFIFO'
        elif strategy == ReplacementStrategy.LRU:
            return 'ReplacementStrategy_LRU'
        else:
            return f'ReplacementStrategy_{strategy.name}'


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

    def to_jmt_text(self) -> str:
        """Convert to JMT's descriptive identifier.

        JMT's CommonConstants.STATION_SCHEDULING_POLICY_* expects the long
        human-readable form, which the bare names do not match, so JSIM XML must
        be written with this method (mirrors MATLAB HeteroSchedPolicy.toJMTText).
        """
        mapping = {
            HeteroSchedPolicy.ORDER: 'Order (Assign according to order below)',
            HeteroSchedPolicy.ALIS: 'ALIS (Assign Longest Idle Server)',
            HeteroSchedPolicy.ALFS: 'ALFS (Assign Least Flexible Server)',
            HeteroSchedPolicy.FAIRNESS: 'Fairness (Move back server type when used)',
            HeteroSchedPolicy.FSF: 'FSF (Fastest Servers First)',
            HeteroSchedPolicy.RAIS: 'RAIS (Random Assignment to Idle Servers)',
        }
        return mapping.get(self, 'Order (Assign according to order below)')

    # Aliases for from_text
    from_string = from_text
    fromString = from_text
    toText = to_text
    toJMTText = to_jmt_text


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

    def _invalidate_struct(self) -> None:
        """Discard the model's cached NetworkStruct after a rate-scaling change.

        Without this a setter called AFTER the first getStruct() is silently
        dropped: the solver reads the stale sn and solves the model unscaled,
        with no error. Twin of MATLAB Station.invalidateStruct.
        """
        model = getattr(self, '_model', None)
        if model is not None and hasattr(model, 'reset_struct'):
            model.reset_struct()

    def _invalidate_java(self) -> None:
        """Invalidate cached object."""
        self._java_obj = None

    def to_java(self):
        """Convert to Java object for JVM interoperability.

        Not available in the native Python implementation; the native Python
        solver operates independently of the JVM.
        """
        raise NotImplementedError(
            "to_java() is not available in the native Python implementation. "
            "The native Python solver operates independently of the JVM."
        )

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}('{self._name}', index={self.get_index()})"

    # API-parity fallback: resolve documented alias spellings
    # (snake/camel/accessor-prefix) to the implemented methods (see _aliasing.py)
    def __getattr__(self, name):
        from .._aliasing import alias_getattr
        return alias_getattr(self, name)


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
        self._spawn_class = None  # Class injected at the same station on each completion (LQN phase-2)
        self._completes = True  # Whether this class completes a visit (for cycle time calculation)
        self._is_reference_class = False  # Whether this class is the reference class for its chain (for WN computation)
        self._immediate_feedback = False  # Whether immediate feedback is enabled globally for this class

    @property
    def jobclass_type(self) -> JobClassType:
        """Get the job class type."""
        return self._jobclass_type

    def summary(self) -> None:
        """Print the one-line class summary, twin of MATLAB JobClass.summary."""
        print(f"Class ({self._jobclass_type.name}): {self.name}")

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

    @property
    def is_reference_class(self) -> bool:
        """
        Get whether this class is the reference class for its chain.

        When is_reference_class=True, this class's visits at the reference station
        are used as the denominator for WN (residence time) computation.
        Typically set to True for the TASK class in LN layer models.
        """
        return self._is_reference_class

    @is_reference_class.setter
    def is_reference_class(self, value: bool) -> None:
        """Set whether this class is the reference class for its chain."""
        self._is_reference_class = value
        self._invalidate_java()

    def setReferenceClass(self, value: bool) -> None:
        """Set whether this class is the reference class for its chain (MATLAB-compatible name)."""
        self.is_reference_class = value

    def set_reference_class(self, value: bool) -> None:
        """Set whether this class is the reference class (snake_case alias)."""
        self.setReferenceClass(value)

    def isReferenceClass(self) -> bool:
        """Get whether this class is the reference class (MATLAB-compatible name)."""
        return self._is_reference_class

    @property
    def immediateFeedback(self) -> bool:
        """Get whether immediate feedback is enabled globally for this class."""
        return self._immediate_feedback

    @immediateFeedback.setter
    def immediateFeedback(self, value: bool) -> None:
        """Set whether immediate feedback is enabled globally for this class."""
        self._immediate_feedback = bool(value)
        self._invalidate_java()

    def setImmediateFeedback(self, value: bool) -> None:
        """Set immediate feedback for this class (MATLAB-compatible name)."""
        self.immediateFeedback = value

    def hasImmediateFeedback(self) -> bool:
        """Check if immediate feedback is enabled for this class."""
        return self._immediate_feedback

    # Snake_case aliases
    set_immediate_feedback = setImmediateFeedback
    has_immediate_feedback = hasImmediateFeedback

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

    def set_patience(self, *args) -> None:
        """
        Set the class-level patience (impatience) distribution.

        Mirrors MATLAB ``JobClass.setPatience`` and Java ``JobClass.setPatience``:

        - ``set_patience(distribution)``: defaults the impatience type to
          RENEGING.
        - ``set_patience(impatience_type, distribution)``: sets both.

        A node-scoped patience set with ``Queue.set_patience(jobclass, ...)``
        takes precedence over the class-level one.

        Args:
            impatience_type: ImpatienceType (RENEGING or BALKING). The legacy
                strings 'reneging'/'balking' are also accepted.
            distribution: Patience distribution.
        """
        from .base import ImpatienceType  # local import: defined in this module
        if len(args) == 1:
            impatience_type, distribution = ImpatienceType.RENEGING, args[0]
        elif len(args) == 2:
            impatience_type, distribution = args
        else:
            raise TypeError(
                "set_patience expects (distribution) or "
                "(impatience_type, distribution), got %d arguments" % len(args))
        if isinstance(impatience_type, str):
            impatience_type = ImpatienceType.from_text(impatience_type)
        if not isinstance(impatience_type, ImpatienceType):
            raise TypeError(
                "impatience_type must be an ImpatienceType, got %r"
                % (impatience_type,))
        self._patience_type = impatience_type
        self._patience_distribution = distribution
        self._invalidate_java()

    def get_patience(self):
        """Get patience distribution (if any)."""
        return self._patience_distribution

    def get_impatience_type(self):
        """Get the class-level ImpatienceType (or None if no patience is set)."""
        return self._patience_type

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

    def set_spawn_class(self, spawn_class) -> None:
        """Set the class injected at the same station on each service
        completion of this class (LQN phase-2 continuation)."""
        self._spawn_class = spawn_class
        self._invalidate_java()

    def get_spawn_class(self):
        """Get the spawn-on-completion class."""
        return self._spawn_class

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

    def getNumberOfJobs(self):
        """Number of jobs in this class: population for closed classes
        (overridden there), 0 for open classes (matching the JAR)."""
        return 0

    def setNumberOfJobs(self, njobs):
        """Set the class population; only meaningful for closed classes."""
        raise NotImplementedError(
            "setNumberOfJobs is only available for closed classes")

    get_number_of_jobs = getNumberOfJobs
    set_number_of_jobs = setNumberOfJobs

    def setReferenceStation(self, station) -> None:
        """Alias for set_reference_station (CamelCase)."""
        return self.set_reference_station(station)

    def getReferenceStation(self):
        """Alias for get_reference_station (CamelCase)."""
        return self.get_reference_station()

    def isReferenceStation(self, station) -> bool:
        """Alias for is_reference_station (CamelCase)."""
        return self.is_reference_station(station)

    def setPatience(self, *args) -> None:
        """Alias for set_patience (CamelCase).

        Accepts ``setPatience(distribution)`` or
        ``setPatience(impatience_type, distribution)``. It previously forwarded
        a single argument into the two-argument-only set_patience, so every
        call raised TypeError.
        """
        return self.set_patience(*args)

    def getPatience(self):
        """Alias for get_patience (CamelCase)."""
        return self.get_patience()

    def getImpatienceType(self):
        """Alias for get_impatience_type (CamelCase)."""
        return self.get_impatience_type()

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

    @property
    def _node_index(self) -> int:
        """Get the 0-based node index (alias for _index for routing compatibility)."""
        return self._index

    @_node_index.setter
    def _node_index(self, value: int) -> None:
        """Set the 0-based node index."""
        self._index = value

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

    def remove_job_class(self, jobclass) -> None:
        """
        Remove all per-class configuration referencing a job class from this node.

        Every per-class setting on a node (service and arrival processes,
        scheduling parameters, class capacities, routing strategies, balking,
        retrial, patience, ...) is stored in a dictionary keyed by the JobClass
        object, or by a tuple containing it. Dropping those entries therefore
        removes the class from this node exhaustively, and stays correct when a
        new per-class dictionary is added to a node type. Subclasses override
        this to handle per-class state that is not dictionary-based, such as the
        class-switching matrix of a ClassSwitch.

        Called by Network.remove_class; mirrors Node.removeJobClass in the JAR.

        Args:
            jobclass: the job class being removed from the model
        """
        def _refers(key) -> bool:
            if key is jobclass:
                return True
            if isinstance(key, tuple):
                return any(element is jobclass for element in key)
            return False

        for _attr, value in list(vars(self).items()):
            if not isinstance(value, dict):
                continue
            for key in [k for k in list(value.keys()) if _refers(k)]:
                del value[key]
            # see _kb/01-model-classes.md (Model-to-model transformations) for rationale
            for key in [k for k, v in list(value.items()) if v is jobclass]:
                del value[key]

    removeJobClass = remove_job_class

    def remap_job_classes(self, new_classes) -> None:
        """
        Re-key every per-class setting of this node onto a new list of classes.

        A deep copy of a node carries orphan JobClass copies as dictionary keys,
        which no lookup with the copied model's own class objects can ever hit,
        so the setting silently reverts to its default. This remaps keys (and
        values that are themselves classes, such as the hit and miss classes of
        a cache) by class index, covering every per-class dictionary of every
        node type. Called by Network.copy.

        Args:
            new_classes: the class list the node must be re-keyed onto
        """
        def _remap(obj):
            if isinstance(obj, JobClass):
                idx = getattr(obj, '_index', None)
                if idx is not None and 0 <= idx < len(new_classes):
                    return new_classes[idx]
                return obj
            if isinstance(obj, tuple):
                return tuple(_remap(element) for element in obj)
            if isinstance(obj, dict):
                return {_remap(k): _remap(v) for k, v in obj.items()}
            if isinstance(obj, list):
                return [_remap(element) for element in obj]
            return obj

        def _holds_class(value) -> bool:
            for k, v in value.items():
                if isinstance(k, JobClass) or isinstance(v, JobClass):
                    return True
                if isinstance(k, tuple) and any(isinstance(e, JobClass) for e in k):
                    return True
                if isinstance(v, dict) and _holds_class(v):
                    return True
            return False

        def _holds_class_seq(value) -> bool:
            return any(isinstance(element, JobClass) for element in value)

        state = vars(self)
        for attr, value in list(state.items()):
            if isinstance(value, dict) and value and _holds_class(value):
                state[attr] = _remap(value)
            elif isinstance(value, (list, tuple)) and value and _holds_class_seq(value):
                # Per-class SEQUENCES are orphaned by a deep copy exactly like
                # per-class dictionaries: Source._marked_classes is a list, and
                # leaving it keyed on the original classes makes the copy's
                # markidx refresh fail to find them.
                state[attr] = _remap(value)

    remapJobClasses = remap_job_classes

    def is_stateful(self) -> bool:
        """Check if node is stateful (can have jobs)."""
        return False

    def is_station(self) -> bool:
        """Check if node is a station (can serve jobs)."""
        return False

    def has_class_switching(self) -> bool:
        """Check if this node switches the class of the jobs traversing it.

        Twin of the MATLAB Node.hasClassSwitching, which tests whether the
        service section is a ClassSwitcher; ClassSwitch and Cache nodes are the
        two that carry one.
        """
        return False

    def summary(self) -> None:
        """Print the one-line node summary, twin of MATLAB Node.summary."""
        print(f"\nNode: {self.name}")

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
                    For WRROBIN: (destination_node, weight) - can be called multiple times
                    For SQ: (d,) - number of sampled destinations
        """
        self._routing_strategies[jobclass] = strategy
        if params:
            # For WRROBIN, store weights per destination instead of overwriting
            # Use value comparison to handle different RoutingStrategy enum definitions
            strategy_value = strategy.value if hasattr(strategy, 'value') else int(strategy)
            if strategy_value == RoutingStrategy.SQ.value:
                # see _kb/01-model-classes.md (Enumerations) for rationale
                d_par = int(params[0]) if len(params) >= 1 and params[0] is not None else 2
                if d_par < 1:
                    raise ValueError("SQ parameter d must be a positive integer")
                if len(params) >= 2 and params[1] is not None and int(params[1]) != 0:
                    raise ValueError(
                        "SQ with dispatcher memory is not supported. Only SQ(d) is available.")
                self._routing_params[jobclass] = (d_par,)
            elif strategy_value == RoutingStrategy.WRROBIN.value and len(params) >= 2:
                destination, weight = params[0], params[1]
                # see _kb/01-model-classes.md (Enumerations) for rationale
                if not hasattr(self, '_routing_weights') or self._routing_weights is None:
                    self._routing_weights = {}
                if jobclass not in self._routing_weights:
                    self._routing_weights[jobclass] = {}
                self._routing_weights[jobclass][destination] = weight
            else:
                self._routing_params[jobclass] = params
        self._invalidate_java()

    def set_state_dep_routing(self, jobclass, departure, branches, level, C, d) -> None:
        """Declare this node the entry center of a Krzesinski SDR subnetwork.

        Krzesinski, A. E., "Multiclass Queueing Networks with State-Dependent
        Routing", Performance Evaluation 7(2):125-143, 1987.

        departure is the departure center d of Q(V,V), which may be this node
        itself in a central server model. branches follows the paper's own
        indexing: branches[0] must be empty because branch index 1 denotes the
        complement M-V, and branches[b] lists the nodes of branch b with its
        entry center first and its departure center last. level[b] is the index
        t of the subnetwork with B_b in V_t - V_{t+1}. C is the length-T list of
        coefficients C_t and d the T x B matrix of coefficients d_tb, read for
        1 <= t <= level[b] and b >= 1 in 0-based terms.

        Negative C_t and positive d_tb make the routing prefer the least
        congested branches and impose the population bounds
        m_b <= d_tb/(-C_t) and v_t <= D_tt/(-C_t).
        """
        import numpy as _np
        if not branches or branches[0]:
            raise ValueError('branches[0] must be empty: branch index 1 denotes '
                             'the complement M-V')
        B = len(branches)
        T = len(C)
        if len(level) != B:
            raise ValueError('level must have one entry per branch index, '
                             'including the unused index 0')
        dm = _np.atleast_2d(_np.asarray(d, dtype=float))
        if dm.shape[0] < T or dm.shape[1] < B:
            raise ValueError('d must be at least %dx%d' % (T, B))
        branch, entry_of, departure_of = [None] * B, [None] * B, [None] * B
        for b in range(1, B):
            bn = branches[b]
            if not isinstance(bn, (list, tuple)):
                bn = [bn]
            if not bn:
                raise ValueError('branch %d is empty' % (b + 1))
            branch[b] = list(bn)
            entry_of[b] = bn[0]
            departure_of[b] = bn[-1]
        self._routing_strategies[jobclass] = RoutingStrategy.SDR
        self._routing_params[jobclass] = ({
            'departure': departure,
            'branch': branch,
            'entryOf': entry_of,
            'departureOf': departure_of,
            'level': list(level),
            'C': list(C),
            'd': dm,
        },)
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

    def get_routing_weight(self, jobclass: JobClass, destination) -> float:
        """Get the WRROBIN routing weight for a class and destination (default 1.0)."""
        w = getattr(self, '_routing_weights', None)
        if not w or jobclass not in w:
            return 1.0
        return w[jobclass].get(destination, 1.0)

    getRoutingWeight = get_routing_weight

    def get_prob_routing(self, jobclass: JobClass):
        """Get probabilistic routing for a job class as {destination: prob}."""
        if not hasattr(self, '_prob_routing'):
            return {}
        return dict(self._prob_routing.get(jobclass, {}))

    getProbRouting = get_prob_routing

    def set_prob_routing(self, jobclass: JobClass, destination, prob: float) -> None:
        """
        Set probabilistic routing to a destination node.

        Args:
            jobclass: Job class
            destination: Destination node
            prob: Routing probability (0 to 1)
        """
        self.set_routing(jobclass, RoutingStrategy.PROB, destination, prob)

    # CamelCase aliases
    setModel = set_model
    getModel = get_model
    isStateful = is_stateful
    isStation = is_station
    hasClassSwitching = has_class_switching
    setRouting = set_routing
    getRouting = get_routing
    setProbRouting = set_prob_routing


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
        self._state_prior = None  # Prior probability over initial states

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

    def getStatePrior(self):
        """Get the prior probability over initial states."""
        return self._state_prior

    def setStatePrior(self, prior):
        """Set the prior probability over initial states."""
        if prior is not None:
            self._state_prior = np.array(prior).flatten()
        else:
            self._state_prior = None

    get_state_prior = getStatePrior
    set_state_prior = setStatePrior

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
        self._state_explicitly_set = True
        self._invalidate_java()
        # Invalidate the model struct so _refresh_state() runs again
        model = self.get_model() if hasattr(self, 'get_model') else None
        if model is not None and hasattr(model, '_has_struct'):
            model._has_struct = False
            model._sn = None

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
    """Abstract base class for stations (nodes that buffer jobs)."""

    def get_balking(self, jobclass):
        """Get balking config for a job class. Returns (strategy, thresholds) or (None, [])."""
        if not hasattr(self, '_balking_strategies'):
            return None, []
        return (self._balking_strategies.get(jobclass),
                self._balking_thresholds.get(jobclass, []))

    def has_balking(self, jobclass):
        """Check if a job class has balking configured."""
        strategy, _ = self.get_balking(jobclass)
        return strategy is not None

    def get_retrial(self, jobclass):
        """Get retrial config for a job class. Returns (delay_dist, max_attempts) or (None, -1)."""
        if not hasattr(self, '_retrial_delays'):
            return None, -1
        return (self._retrial_delays.get(jobclass),
                self._retrial_max_attempts.get(jobclass, -1))

    def has_retrial(self, jobclass):
        """Check if a job class has retrial configured."""
        delay_dist, _ = self.get_retrial(jobclass)
        return delay_dist is not None

    def get_limit(self):
        """LPS concurrency limit, or the number of servers when unset."""
        lim = getattr(self, '_lps_limit', None)
        if lim is None:
            return self.number_of_servers
        return lim

    def get_polling_type(self):
        """Polling type for this station (None if unset)."""
        return getattr(self, '_polling_type', None)

    def get_total_num_of_servers(self):
        """Total number of servers across all server types."""
        if not getattr(self, '_server_types', None):
            return self.number_of_servers
        return sum(st.get_num_of_servers() for st in self._server_types)

    getBalking = get_balking
    hasBalking = has_balking
    getRetrial = get_retrial
    hasRetrial = has_retrial
    getLimit = get_limit
    getPollingType = get_polling_type
    getTotalNumOfServers = get_total_num_of_servers

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
        # UNSET, not DROP. MATLAB leaves Station.dropRule empty and the JAR
        # returns null from getDropRule; both derive the rule from the final
        # capacity and class type in refreshCapacity. Presetting DROP here made
        # the matching derivation in Network._refresh_capacity dead code for
        # every station nobody configured, so sn.droprule read DROP where the
        # other two codebases read WAITQ.
        self._drop_rule = None
        self._load_depend_scaling = None
        self._class_depend_scaling = None
        self._class_depend_scaling_peak = None
        self._joint_depend_scaling = None
        self._joint_depend_scaling_peak = None
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
        self._invalidate_struct()

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
        # sn.cap/sn.classcap are DERIVED (refresh_capacity folds this value
        # together with the per-class caps and the chain population), so a
        # cached struct does not see the new buffer. Without this a
        # set_capacity called after the first getStruct() was silently dropped
        # and every sn-reading solver answered the UNBOUNDED model -- SolverCTMC
        # returned 1.1475 jobs for a buffer of 1. Twin of MATLAB
        # Station.setCapacity.
        self._invalidate_struct()

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
        self._invalidate_struct()

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
        self._invalidate_struct()

    def get_load_dependence(self) -> Optional[np.ndarray]:
        """Get load-dependent scaling."""
        return self._load_depend_scaling

    def set_class_dependence(self, beta, peak_rate_per_class=None) -> None:
        """
        Set class-dependent service rates beta_{i,r}(n).

        Args:
            beta: Callable of the per-class population vector n at this station.
                It may return either a scalar, i.e. a chain-independent scaling
                beta_i(n) shared by every class, or an array of length R giving
                the per-class scalings [beta_{i,1}(n), ..., beta_{i,R}(n)]. The
                per-class form expresses Sauer's chain-dependent service rates
                mu_{r,i}(n) (Sauer 1983, eq. (40)), and subsumes the joint
                dependence (LJD/LJCD) tables of earlier releases: a tabulated
                dependence is obtained by having the callable index the table
                (see api.pfqn.ljd.fes_beta_handle).
            peak_rate_per_class: REQUIRED peak (maximum) rate scaling per class
                (e.g. the effective number of servers). Used to normalize
                utilization as Util = T*S/peak, matching the T*S/c convention of
                ordinary multiserver stations. A scalar is broadcast to every
                class.
        """
        if peak_rate_per_class is None:
            raise ValueError(
                "Class dependence requires an explicit peak rate: "
                "set_class_dependence(beta, peak_rate_per_class). Pass a scalar "
                "(identical peak for every class) or a per-class vector.")
        peak = np.atleast_1d(np.asarray(peak_rate_per_class, dtype=float)).ravel()
        if peak.size == 0 or np.any(peak <= 0):
            raise ValueError("peak_rate_per_class must be a positive scalar or per-class vector.")
        self._class_depend_scaling = beta
        self._class_depend_scaling_peak = peak
        self._invalidate_java()
        self._invalidate_struct()

    def get_class_dependence(self):
        """Get class-dependent scaling."""
        return self._class_depend_scaling

    def get_class_dependence_peak(self):
        """Get the declared peak class-dependent rate scaling (per class)."""
        return getattr(self, '_class_depend_scaling_peak', None)

    def set_joint_dependence(self, eta, peak_rate_per_class=None) -> None:
        """
        Set joint-dependent service rates eta_i(n) (NON-product-form).

        Unlike set_class_dependence (product-form beta_{i,r} which must depend
        on the own-class marginal n_{i,r} and preserves the BCMP product form),
        eta may read the joint per-class population vector arbitrarily (e.g.
        min(ni[0], c), which reads a foreign class marginal). It is therefore
        NON-product-form: solvers treat it as an approximation with no
        exactness/uniqueness guarantee. Evaluated numerically identically to
        the class-dependence handle; the distinction is semantic and is carried
        by the separate sn.jdscaling field.

        Args:
            eta: Callable of the joint per-class population vector n at this
                station. It may return either a scalar (shared by every class)
                or an array of length R giving per-class scalings [eta_{i,1}(n),
                ..., eta_{i,R}(n)] (Sauer's chain-dependent rate mu_{r,i}(n)).
            peak_rate_per_class: REQUIRED peak (maximum) rate scaling per class.
                Used to normalize utilization as Util = T*S/peak. A scalar is
                broadcast to every class.
        """
        if peak_rate_per_class is None:
            raise ValueError(
                "Joint dependence requires an explicit peak rate: "
                "set_joint_dependence(eta, peak_rate_per_class). Pass a scalar "
                "(identical peak for every class) or a per-class vector.")
        peak = np.atleast_1d(np.asarray(peak_rate_per_class, dtype=float)).ravel()
        if peak.size == 0 or np.any(peak <= 0):
            raise ValueError("peak_rate_per_class must be a positive scalar or per-class vector.")
        self._joint_depend_scaling = eta
        self._joint_depend_scaling_peak = peak
        self._invalidate_java()
        self._invalidate_struct()

    def get_joint_dependence(self):
        """Get joint-dependent scaling."""
        return getattr(self, '_joint_depend_scaling', None)

    def get_joint_dependence_peak(self):
        """Get the declared peak joint-dependent rate scaling (per class)."""
        return getattr(self, '_joint_depend_scaling_peak', None)

    def set_drop_rule(self, *args) -> None:
        """
        Set drop strategy for finite capacity.

        Two forms are supported, mirroring the MATLAB/Java API:

        - ``set_drop_rule(rule)``: set the same drop strategy for all classes.
        - ``set_drop_rule(jobclass, rule)``: set the drop strategy for a single
          class. Per-class rules are stored in a dict keyed by the job class.

        Args:
            rule: Drop strategy (DROP, BAS, or WaitingQueue).
            jobclass: (per-class form) the job class the rule applies to.
        """
        if len(args) == 1:
            rule = args[0]
            self._drop_rule = rule
        elif len(args) == 2:
            jobclass, rule = args
            if not isinstance(self._drop_rule, dict):
                # a station-level rule already set is the DEFAULT for the other
                # classes, so seed the dict with it rather than discarding it
                previous = self._drop_rule
                self._drop_rule = {}
                if previous is not None and self._model is not None:
                    for k in getattr(self._model, '_classes', []):
                        self._drop_rule[k] = previous
            self._drop_rule[jobclass] = rule
        else:
            raise TypeError(
                "set_drop_rule expects (rule) or (jobclass, rule), "
                "got %d arguments" % len(args))
        self._invalidate_java()
        self._invalidate_struct()

    def get_drop_rule(self, jobclass=None) -> Optional[DropStrategy]:
        """Get drop strategy, or None when the station has none set.

        With no argument returns the station-level drop rule (or the dict of
        per-class rules). With a ``jobclass`` argument returns the rule for
        that class, falling back to the station-level rule.

        None is the unset state, matching MATLAB's empty ``dropRule`` and the
        JAR's null: the effective rule is then derived from the capacity and
        the class type when the struct is refreshed, so a caller that needs the
        effective one must read ``sn.droprule`` rather than defaulting here.
        """
        if jobclass is not None and isinstance(self._drop_rule, dict):
            return self._drop_rule.get(jobclass, None)
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
