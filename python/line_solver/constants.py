"""
Constants and enumerations for LINE queueing network models.

This module defines the various constants, enumerations, and strategies
used throughout LINE for specifying model behavior, including:

- Scheduling strategies (FCFS, LCFS, PS, etc.)
- Routing strategies (PROB, RAND, etc.)
- Node types (SOURCE, QUEUE, SINK, etc.)
- Job class types (OPEN, CLOSED)
- Solver types and options
- Activity precedence types for layered networks
- Call types and drop strategies

These constants ensure type safety and consistency across the API.
"""

from enum import Enum, auto


class ActivityPrecedenceType(Enum):
    """
    Types of activity precedence relationships in layered networks.

    These specify how activities are ordered and synchronized:
    - PRE_SEQ: Sequential prerequisite (must complete before)
    - PRE_AND: AND prerequisite (all must complete before)
    - PRE_OR: OR prerequisite (any must complete before)
    - POST_SEQ: Sequential post-condition
    - POST_AND: AND post-condition
    - POST_OR: OR post-condition
    - POST_LOOP: Loop post-condition
    - POST_CACHE: Cache post-condition
    """
    PRE_SEQ = auto()
    PRE_AND = auto()
    PRE_OR = auto()
    POST_SEQ = auto()
    POST_AND = auto()
    POST_OR = auto()
    POST_LOOP = auto()
    POST_CACHE = auto()


class CallType(Enum):
    """
    Types of calls between tasks in layered networks.

    - SYNC: Synchronous call (caller waits for response)
    - ASYNC: Asynchronous call (caller continues immediately)
    - FWD: Forward call (caller terminates, response goes to caller's caller)
    """
    SYNC = auto()
    ASYNC = auto()
    FWD = auto()


class DropStrategy(Enum):
    """
    Strategies for handling queue overflow and capacity limits.

    - WaitingQueue: Jobs wait in a waiting queue when capacity is exceeded
    - Queue: Alias for WaitingQueue
    - Drop: Jobs are dropped (lost) when capacity is exceeded
    - BlockingAfterService: Jobs are blocked after service completion
    - BlockingBeforeService: Jobs are blocked before service starts
    - ReServiceOnRejection: Rejected jobs are re-served at the upstream station

    Values match the MATLAB DropStrategy constants and the JAR
    jline.lang.constant.DropStrategy ids, which are the interchange encoding of
    sn.droprule and sn.regionrule. Keep the three Python definitions of this
    enum (here, lang/base.py, api/sn/network_struct.py) numerically identical:
    they are written and read by different modules over the same sn fields.
    """
    WaitingQueue = -1
    Drop = 1
    BlockingAfterService = 2
    BlockingBeforeService = 3
    ReServiceOnRejection = 4
    Queue = -1  # alias for WaitingQueue


class DepartureDiscipline(Enum):
    """
    Departure disciplines for the depository of a queueing place (QPN semantics).

    A queueing place serves tokens in its embedded queue and, on service
    completion, moves them to a depository from which they become available to
    the output transitions. The departure discipline governs the order in which
    depository tokens become available.

    - NORMAL: tokens available immediately upon service completion (standard QPN)
    - FIFO: tokens available in their order of arrival to the depository
    """
    NORMAL = 0
    FIFO = 1


class SignalType(Enum):
    """
    Types of signals for signal classes in G-networks and related models.

    This is the single canonical definition: lang/classes.py re-exports it
    rather than defining a second enum. Two coexisting definitions used to be
    disambiguated only by the import order in line_solver/__init__.py, and the
    losing definition carried auto() ordinals (1-based) that would have
    mis-decoded against the 0-based MATLAB/Java enums on the JSON wire.

    The member values are the lowercase names used on the JSON wire.

    Attributes:
        NEGATIVE: Removes a job from the destination queue (G-network negative customer)
        REPLY: Triggers a reply action
        CATASTROPHE: Removes ALL jobs from the destination queue
    """
    NEGATIVE = 'negative'
    REPLY = 'reply'
    CATASTROPHE = 'catastrophe'


class RemovalPolicy(Enum):
    """
    Removal policies for negative signals in G-networks.

    Single canonical definition; see the note on SignalType above. The member
    values are the lowercase names used on the JSON wire.

    Attributes:
        RANDOM: Select job uniformly at random from all jobs at the station
        FCFS: Remove the oldest job (first arrived)
        LCFS: Remove the newest job (last arrived)
    """
    RANDOM = 'random'
    FCFS = 'fcfs'
    LCFS = 'lcfs'


class EventType(Enum):
    """
    Types of events in discrete-event simulation.

    - INIT: Initialization event
    - LOCAL: Local processing event
    - ARV: Job arrival event
    - DEP: Job departure event
    - PHASE: Phase transition event in multi-phase processes
    - READ: Cache read event
    - STAGE: Staging area event
    """
    INIT = auto()
    LOCAL = auto()
    ARV = auto()
    DEP = auto()
    PHASE = auto()
    READ = auto()
    STAGE = auto()
    ENABLE = auto()
    FIRE = auto()
    PRE = auto()
    POST = auto()
    RENEGE = auto()  # a waiting job abandons the queue (impatience)
    RETRY = auto()   # an orbiting job retries entry into a retrial station
    SWITCH = auto()  # the server of a polling station advances its switchover timer
    FAILURE = auto()  # the server of a station breaks down (goes from up to down)
    REPAIR = auto()   # the server of a station is repaired (goes from down to up)


class JobClassType(Enum):
    """
    Types of job classes in queueing networks.

    - OPEN: Open class (jobs arrive from outside the system)
    - CLOSED: Closed class (fixed population circulating in the system)
    - DISABLED: Disabled class (not currently active)
    """
    OPEN = auto()
    CLOSED = auto()
    DISABLED = auto()


def __getattr__(name):
    """Forward JoinStrategy to its single definition in lang.base.

    It used to be redefined here as a second, independent enum. Nothing stored
    it: a Join node holds lang.base.JoinStrategy, so the copy compared equal to
    no strategy any model carried. The forward is lazy because importing
    lang.base at module level would close an import cycle back through
    lang.network.
    """
    if name == 'JoinStrategy':
        from .lang.base import JoinStrategy
        return JoinStrategy
    raise AttributeError("module %r has no attribute %r" % (__name__, name))


class MetricType(Enum):
    """
    Types of performance metrics that can be computed.
    """
    ResidT = auto()
    RespT = auto()
    DropRate = auto()
    QLen = auto()
    QueueT = auto()
    FCRWeight = auto()
    FCRMemOcc = auto()
    FJQLen = auto()
    FJRespT = auto()
    RespTSink = auto()
    SysDropR = auto()
    SysQLen = auto()
    SysPower = auto()
    SysRespT = auto()
    SysTput = auto()
    Tput = auto()
    ArvR = auto()
    TputSink = auto()
    Util = auto()
    TranQLen = auto()
    TranUtil = auto()
    TranTput = auto()
    TranRespT = auto()
    Tard = auto()
    SysTard = auto()


class Metric:
    """An output metric of a Solver, such as a performance index."""

    def __init__(self, metric_type, job_class, station=None):
        self.type = metric_type
        self.job_class = job_class
        self.station = station
        self.disabled = False
        self.transient = False


class TranResult:
    """Container for transient result time series with attribute access."""

    def __init__(self, t, metric):
        self.t = t
        self.metric = metric


class NodeType(Enum):
    """
    Types of nodes in queueing network models.
    """
    Transition = auto()
    Place = auto()
    Fork = auto()
    Router = auto()
    Cache = auto()
    Logger = auto()
    ClassSwitch = auto()
    Delay = auto()
    Source = auto()
    Sink = auto()
    Join = auto()
    Queue = auto()

    @staticmethod
    def from_line(obj):
        obj_str = str(obj)
        return getattr(NodeType, obj_str, None)


class ProcessType(Enum):
    """
    Types of stochastic processes for arrivals and service times.
    """
    EXP = auto()
    ERLANG = auto()
    DISABLED = auto()
    IMMEDIATE = auto()
    HYPEREXP = auto()
    APH = auto()
    COXIAN = auto()
    PH = auto()
    MAP = auto()
    UNIFORM = auto()
    DET = auto()
    GAMMA = auto()
    PARETO = auto()
    WEIBULL = auto()
    LOGNORMAL = auto()
    MMPP2 = auto()
    REPLAYER = auto()
    TRACE = auto()
    COX2 = auto()
    BINOMIAL = auto()
    POISSON = auto()
    BMAP = auto()
    MMAP = auto()
    # see _kb/07-cross-language-parity.md (ProcessType enum numeric mismatch) for rationale
    DUNIFORM = auto()
    BERNOULLI = auto()
    PRIOR = auto()
    GEOMETRIC = auto()
    ME = auto()
    RAP = auto()
    DISCRETESAMPLER = auto()
    ZIPF = auto()
    DMAP = auto()
    EMPIRICALCDF = auto()
    NHPP = auto()

    @staticmethod
    def isMarkovian(t):
        """True when ``sn.proc`` carries an exact matrix representation of a
        process of type ``t``: a genuine (D0, D1) pair, or its
        matrix-exponential analogue for ME/RAP.

        The distinction is about ``sn.proc``, NOT about what ``getProcess``
        returns. ``getProcess`` hands back raw distribution PARAMETERS for
        several non-Markovian families -- Gamma, Weibull, Lognormal, Pareto and
        Uniform return two scalars (Pareto ``{alpha, k}``, Uniform
        ``{min, max}``) -- and the network refresh replaces those with
        ``map_erlang(mean, n)`` before storing them in ``sn.proc``, where
        ``n = ceil(1/SCV)`` capped at 100 (``n = 20`` when ``SCV < CoarseTol``).
        That fit matches the mean, and matches the SCV only when ``SCV <= 1``:
        Pareto with SCV 64 gives ``n = 1``, a single exponential of SCV 1. So
        for these types ``sn.proc`` is an approximation, not the law that was
        requested, and nothing on the cell says so -- the only signal is
        ``sn.procid``.

        Solvers that read ``sn.proc`` as if it were the exact law must gate on
        this predicate. It is the procid-level counterpart of the JAR's
        ``Distribution.isMarkovian()``, i.e. of the Markovian class hierarchy,
        so the two lists must stay in step.
        """
        return t in (ProcessType.EXP, ProcessType.ERLANG, ProcessType.HYPEREXP,
                     ProcessType.PH, ProcessType.APH, ProcessType.MAP,
                     ProcessType.COXIAN, ProcessType.COX2, ProcessType.MMPP2,
                     ProcessType.ME, ProcessType.RAP, ProcessType.DMAP,
                     ProcessType.BMAP, ProcessType.MMAP)

    @staticmethod
    def fromString(obj):
        mapping = {
            "Exp": ProcessType.EXP,
            "Erlang": ProcessType.ERLANG,
            "HyperExp": ProcessType.HYPEREXP,
            "PH": ProcessType.PH,
            "APH": ProcessType.APH,
            "MAP": ProcessType.MAP,
            "BMAP": ProcessType.BMAP,
            "Uniform": ProcessType.UNIFORM,
            "Det": ProcessType.DET,
            "Coxian": ProcessType.COXIAN,
            "Gamma": ProcessType.GAMMA,
            "Pareto": ProcessType.PARETO,
            "MMPP2": ProcessType.MMPP2,
            "Replayer": ProcessType.REPLAYER,
            "Trace": ProcessType.TRACE,
            "Immediate": ProcessType.IMMEDIATE,
            "Disabled": ProcessType.DISABLED,
            "Cox2": ProcessType.COX2,
            "Weibull": ProcessType.WEIBULL,
            "Lognormal": ProcessType.LOGNORMAL,
            "Poisson": ProcessType.POISSON,
            "Binomial": ProcessType.BINOMIAL,
            "NHPP": ProcessType.NHPP,
            "MarkedMAP": ProcessType.MMAP,
            "MarkedMMPP": ProcessType.MMAP,
            "MMAP": ProcessType.MMAP,
            # see _kb/07-cross-language-parity.md (EmpiricalCdf naming) for rationale
            "DiscreteUniform": ProcessType.DUNIFORM,
            "Bernoulli": ProcessType.BERNOULLI,
            "Prior": ProcessType.PRIOR,
            "Geometric": ProcessType.GEOMETRIC,
            "ME": ProcessType.ME,
            # see _kb/07-cross-language-parity.md (ProcessType.fromString CME entry) for rationale
            "CME": ProcessType.ME,
            "RAP": ProcessType.RAP,
            "DiscreteSampler": ProcessType.DISCRETESAMPLER,
            "Zipf": ProcessType.ZIPF,
            "DMAP": ProcessType.DMAP,
            "EmpiricalCdf": ProcessType.EMPIRICALCDF,
            "EmpiricalCDF": ProcessType.EMPIRICALCDF,
        }
        return mapping.get(str(obj))


# ReplacementStrategy is defined in lang/base.py to avoid circular imports
# Import it from there: from line_solver.lang.base import ReplacementStrategy


class RoutingStrategy(Enum):
    """
    Strategies for routing jobs between network nodes.
    Values must match MATLAB's RoutingStrategy constants for JMT compatibility.
    """
    RAND = 0
    PROB = 1
    RROBIN = 2
    WRROBIN = 3
    JSQ = 4
    FIRING = 5
    SQ = 6  # KCHOICES is now SQ: shortest queue of d, SQ(d)
    RL = 7
    DISABLED = -1


class SchedStrategy(Enum):
    """
    Scheduling strategies for service stations.

    Values match lang/base.py SchedStrategy for consistency.
    """
    FCFS = 0       # First-Come First-Served
    LCFS = 1       # Last-Come First-Served
    LCFSPR = 2     # LCFS with Preemptive Resume
    LCFSPI = 3     # LCFS with Preemptive Interrupt
    PS = 4         # Processor Sharing
    DPS = 5        # Discriminatory Processor Sharing
    GPS = 6        # Generalized Processor Sharing
    INF = 7        # Infinite Server (Delay)
    RAND = 8       # Random
    HOL = 9        # Head of Line
    SEPT = 10      # Shortest Expected Processing Time
    LEPT = 11      # Longest Expected Processing Time
    SIRO = 12      # Service In Random Order
    SJF = 13       # Shortest Job First
    LJF = 14       # Longest Job First
    POLLING = 15   # Polling
    EXT = 16       # External arrival stream
    LPS = 17       # Limited Processor Sharing
    SETF = 18      # Shortest Elapsed Time First
    DPSPRIO = 19   # DPS with Priorities
    GPSPRIO = 20   # GPS with Priorities
    PSPRIO = 21    # PS with Priorities
    FCFSPR = 22    # FCFS with Preemptive Resume
    EDF = 23       # Earliest Deadline First
    FORK = 24      # Fork node
    JOIN = 25      # Join node
    REF = 26       # Reference task
    EDD = 27       # Earliest Due Date
    SRPT = 28      # Shortest Remaining Processing Time
    SRPTPRIO = 29  # SRPT with Priorities
    LCFSPRIO = 30  # LCFS with Priorities
    LCFSPRPRIO = 31  # LCFSPR with Priorities
    LCFSPIPRIO = 32  # LCFSPI with Priorities
    FCFSPRPRIO = 33  # FCFSPR with Priorities
    FCFSPIPRIO = 34  # FCFSPI with Priorities
    FCFSPRIO = 35  # FCFS with Priorities
    PSJF = 36      # Preemptive Shortest Job First
    FB = 37        # Foreground-Background
    LAS = 38       # Least Attained Service
    LRPT = 39      # Longest Remaining Processing Time
    FSP = 40       # Fair Sojourn Protocol (virtual-PS finish time ranking)
    PAS = 41       # Pass-and-swap (order-independent queue with class swap graph)
    OI = 42        # Order-independent (pass-and-swap specialization with empty/zero swap graph)

    @staticmethod
    def fromString(obj):
        obj_str = str(obj)
        return getattr(SchedStrategy, obj_str, None)

    @staticmethod
    def fromLINEString(sched: str):
        return SchedStrategy.fromString(sched.upper())

    @staticmethod
    def toID(sched):
        return list(SchedStrategy).index(sched)


class SchedStrategyType(Enum):
    """
    Categories of scheduling strategies by preemption behavior.
    """
    PR = auto()
    PNR = auto()
    NP = auto()
    NPPrio = auto()


class ServiceStrategy(Enum):
    """
    Service strategies defining service time dependence.
    """
    LI = auto()
    LD = auto()
    CD = auto()
    SD = auto()


class SolverType(Enum):
    """
    Types of solvers available in LINE.
    """
    AUTO = auto()
    BA = auto()
    CTMC = auto()
    LDES = auto()
    ENV = auto()
    FLUID = auto()
    JMT = auto()
    LN = auto()
    LQNS = auto()
    MAM = auto()
    MVA = auto()
    NC = auto()
    QNS = auto()
    SSA = auto()


class TimingStrategy(Enum):
    """
    Timing strategies for transitions in Petri nets.
    """
    TIMED = auto()
    IMMEDIATE = auto()


class VerboseLevel(Enum):
    """
    Verbosity levels for LINE solver output.
    """
    SILENT = 0
    STD = 1
    DEBUG = 2


class PollingType(Enum):
    """
    Polling strategies for polling systems.
    """
    GATED = auto()
    EXHAUSTIVE = auto()
    KLIMITED = auto()
    DECREMENTING = auto()

    @staticmethod
    def fromString(obj):
        obj_str = str(obj).upper()
        return getattr(PollingType, obj_str, None)


class HeteroSchedPolicy(Enum):
    """
    Scheduling policies for heterogeneous multiserver queues.
    """
    ORDER = auto()
    ALIS = auto()
    ALFS = auto()
    FAIRNESS = auto()
    FSF = auto()
    RAIS = auto()

    @staticmethod
    def fromString(obj):
        obj_str = str(obj).upper()
        return getattr(HeteroSchedPolicy, obj_str, None)


class GlobalConstants:
    """
    Global constants and configuration for the LINE solver.
    """
    Zero = 1e-14
    CoarseTol = 1e-3  # Match MATLAB/JAR default (1.0e-03)
    FineTol = 1e-8  # Match MATLAB's default
    Immediate = 1e8  # 1/FineTol - large but finite rate for immediate service (matches MATLAB)
    MaxInt = 2**31 - 1
    Version = "3.0.7"
    DummyMode = False
    # Latch for the once-per-session library attribution printed by the solvers
    # (mirrors MATLAB GlobalConstants and jline.lang.GlobalConstants).
    _libraryAttributionShown = False

    @classmethod
    def isLibraryAttributionShown(cls) -> bool:
        """True if the library attribution has already been printed."""
        return cls._libraryAttributionShown

    @classmethod
    def setLibraryAttributionShown(cls, value: bool = True) -> None:
        """Record that the library attribution has been printed."""
        cls._libraryAttributionShown = bool(value)

    _instance = None
    _verbose = VerboseLevel.STD

    def __repr__(self):
        return f"GlobalConstants(Version={self.Version}, Verbose={self.getVerbose()})"

    @classmethod
    def getInstance(cls):
        """Get the singleton instance of GlobalConstants."""
        if cls._instance is None:
            cls._instance = cls()
        return cls._instance

    get_instance = getInstance

    @classmethod
    def getVerbose(cls):
        """Get the current verbosity level."""
        return cls._verbose

    get_verbose = getVerbose

    @classmethod
    def setVerbose(cls, verbosity):
        """Set the verbosity level for solver output."""
        if isinstance(verbosity, VerboseLevel):
            cls._verbose = verbosity
        else:
            raise ValueError(f"Invalid verbosity level: {verbosity}")

    set_verbose = setVerbose

    @classmethod
    def getConstants(cls):
        """Get a dictionary of all global constants."""
        return {
            'Zero': cls.Zero,
            'CoarseTol': cls.CoarseTol,
            'FineTol': cls.FineTol,
            'Immediate': cls.Immediate,
            'MaxInt': cls.MaxInt,
            'Version': cls.Version,
            'DummyMode': cls.DummyMode,
            'Verbose': cls.getVerbose()
        }

    get_constants = getConstants


def default_verbose() -> bool:
    """Default solver verbosity, inherited from GlobalConstants.

    True unless the global verbosity level is SILENT, so solver banners
    print by default as in MATLAB and Java.
    """
    return GlobalConstants.getVerbose() != VerboseLevel.SILENT
