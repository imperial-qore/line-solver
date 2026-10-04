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

from enum import Enum, IntEnum, auto


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


class CallType(IntEnum):
    """
    Types of calls between tasks in layered networks.

    THE single definition; layered.py and the LN solver import it. Values are
    MATLAB's (`matlab/src/lang/constant/CallType.m`) and are what `lsn.calltype`
    carries. layered.py used to declare a third copy with STRING values, and
    being the last one the package imported it was the one `line_solver.CallType`
    resolved to -- so a caller comparing it against a member from either of the
    other two got False whichever way round they wrote it.

    - SYNC: Synchronous call (caller waits for response)
    - ASYNC: Asynchronous call (caller continues immediately)
    - FWD: Forward call (caller terminates, response goes to caller's caller)
    """
    SYNC = 1
    ASYNC = 2
    FWD = 3


class DropStrategy(IntEnum):
    """Strategies for handling queue overflow and capacity limits.

    THE single definition; lang/base.py and api/sn/network_struct.py import it.
    All three used to declare it, each carrying a docstring asking the next
    maintainer to keep them numerically identical -- and two of the three were
    plain `Enum` against `IntEnum` here, which never compare equal whatever the
    numbers say. Values are MATLAB's DropStrategy constants and the JAR
    jline.lang.constant.DropStrategy ids, which are the interchange encoding of
    sn.droprule and sn.regionrule.

    The short names are MATLAB's spelling and the CamelCase ones the JAR's, and
    both are kept as MATLAB keeps both: they are the same member, so a name
    from either codebase reads the same model.
    """
    WAITQ = -1              # Wait in queue; also the default when no rule is set
    WaitingQueue = -1       # Alias for WAITQ (the JAR's spelling)
    Queue = -1              # Alias for WAITQ (MATLAB's second spelling)
    DROP = 1                # Drop arriving job
    Drop = 1                # Alias for DROP
    BAS = 2                 # Block after service
    BlockingAfterService = 2    # Alias for BAS
    BBS = 3                 # Block before service
    BlockingBeforeService = 3   # Alias for BBS
    RSRD = 4                # Re-service on rejection
    ReServiceOnRejection = 4    # Alias for RSRD
    RETRIAL = 5             # Job moves to orbit and retries after delay
    Retrial = 5             # Alias for RETRIAL
    RETRIAL_WITH_LIMIT = 6  # Job retries up to max attempts, then drops
    RetrialWithLimit = 6    # Alias for RETRIAL_WITH_LIMIT


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
    # START and PREEMPT are instantaneous tags on the arc of the ARV or DEP
    # transition that causes them, never the active half of an sn.sync entry:
    # they carry no clock, add no state and change no numerical result.
    # PREEMPT is spelled in full because PRE already names the Petri-net
    # pre-arc. REPAIR emits no START: the supported breakdown model resumes the
    # held job rather than restarting it.
    START = auto()    # a job begins or resumes holding a server
    PREEMPT = auto()  # a job holding a server is pushed back into the buffer


class LayeredNetworkElement(IntEnum):
    """Element types in layered queueing networks.

    THE single definition; layered.py, the LN solver and api/lsn import it.
    Values are the JAR's `jline.lang.layered.LayeredNetworkElement` constants,
    which is what `lsn.type` carries (see LayeredNetwork.getStruct). MATLAB has
    no such enumeration: `LayeredNetworkElement.m` is the BASE CLASS the layered
    elements derive from, so the JAR is ground truth here.

    api/lsn/max_multiplicity.py used to declare a third copy numbering
    PROCESSOR 4 and inventing a HOST 5, against the `lsn.type` it reads; it only
    ever compared TASK and ENTRY, which happen to agree, so the wrong numbering
    never showed.
    """
    PROCESSOR = 0
    HOST = 0      # Alias for PROCESSOR (the JAR's second spelling)
    TASK = 1
    ENTRY = 2
    ACTIVITY = 3
    CALL = 4


class JobClassType(IntEnum):
    """Types of job classes in queueing networks.

    THE single definition; lang/base.py imports it. Values are MATLAB's
    (`matlab/src/lang/constant/JobClassType.m`), which is ground truth. This
    module used to declare a plain `Enum` with `auto()` values (OPEN=1) and
    lang/base.py an `IntEnum` with OPEN=0, CLOSED=1 and a SIGNAL member no
    codebase carries -- and no DISABLED, so `DisabledClass` raised
    AttributeError on every construction while this copy, which had DISABLED,
    was the one nothing imported. SIGNAL is gone with the copy: a Signal class
    is stored as CLOSED (see Network._refresh_signals).

    - OPEN: jobs arrive from outside the system and depart it
    - CLOSED: a fixed population circulating in the system
    - DISABLED: an inactive class the solvers do not process
    """
    DISABLED = -1
    OPEN = 0
    CLOSED = 1


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


class NodeType(IntEnum):
    """Node types in a queueing network.

    THE single definition; lang/base.py and api/sn/network_struct.py import it,
    and it is the numbering `sn.nodetype` carries. All three used to declare it:
    the two above agreed, with a note asking the next maintainer to keep them
    that way, while this module declared a CamelCase plain `Enum` with `auto()`
    values that matched neither -- and a plain `Enum` member never compares equal
    to an `IntEnum` one whatever the numbers say, so the five sites that imported
    it silently answered False. They are fixed with this collapse: the console's
    model kind never reported an SPN or a cache, SolverCTMC's WAITQ-plus-SPN
    refusal never fired, SolverFLD's hide_immediate treated every SPN as a plain
    queueing network (int() on a plain Enum raises, and the raise was swallowed),
    and SolverAG's is_source was always False.

    The NUMBERING is native python's own and deliberately not MATLAB's, which
    runs Queue=0, Source=1, ..., Sink=-1, Join=-2. Nothing crosses a codebase
    boundary as an integer here: the JSON wire carries the node type as the
    STRING `toText` returns.
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
    MAPT = auto()
    PHT = auto()
    # The MARKED families. A mark is a label carried by an event, and at a Source
    # it selects the class of the arriving job (Source.set_marked_arrival,
    # sn.markidx). MPH is the renewal special case of MMAP, obtained from a PH
    # with K marked exits by D0 = S and D1k = s_k*alpha; MMAPT is MMAP with a
    # piecewise-constant matrix schedule, and MPHT is the same lowering applied
    # segment by segment. MPH keeps an id of its own rather than aliasing MMAP,
    # so a solver that cannot honour a marked renewal process refuses it
    # explicitly instead of inheriting MMAP's support.
    MPH = auto()
    MMAPT = auto()
    MPHT = auto()
    # BMMAPT crosses the BATCH axis with the two above: a block is indexed by
    # segment, mark and batch size, and one epoch releases a batch of jobs that
    # all carry the same mark. It reduces to MMAPT when every batch size is 1 and
    # to BMAP when the schedule is flat, and like MPH it keeps an id of its own so
    # a solver that cannot release batches refuses it by name instead of
    # inheriting MMAPT's support.
    BMMAPT = auto()

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
                     ProcessType.BMAP, ProcessType.MMAP, ProcessType.MPH,
                     ProcessType.MMAPT, ProcessType.MPHT, ProcessType.BMMAPT)

    @staticmethod
    def isMarked(t):
        """True when the type carries PER-MARK arrival blocks, i.e. an event of
        this process is labelled and the label is meaningful to the model. At a
        Source the label selects the class of the arriving job
        (``Source.set_marked_arrival``, ``sn.markidx``).
        """
        return ProcessType.isMarkedStationary(t) or ProcessType.isMarkedSchedule(t)

    @staticmethod
    def isMarkedStationary(t):
        """True when ``sn.proc`` holds the STATIONARY marked cell, the M3A layout
        ``{D0, D1agg, D11, ..., D1K}``.

        MPH is the renewal special case of MMAP and lowers to exactly that cell
        (``D0 = S``, ``D1k = s_k*alpha``), so every consumer that reads the M3A
        layout serves both and must gate on this predicate rather than on
        equality with MMAP.
        """
        return t in (ProcessType.MMAP, ProcessType.MPH)

    @staticmethod
    def isMarkedSchedule(t):
        """True when ``sn.proc`` holds the MARKED SCHEDULE slot
        ``{breakpoints, D0segs, D1aggsegs, cyclic, markSegs}``.

        MPHt is stored lowered to MMAPt form segment by segment, so one walk
        serves both, exactly as one MAPt walk serves MAPt and PHt.

        BMMAPT IS INCLUDED, and its slot is that one with the batch blocks
        appended, so a consumer gated on this predicate reads a BMMAPt as the
        MMAPt it aggregates down to. That is right for anything time-blind or
        batch-blind and WRONG for anything that releases jobs: an arrival or
        service sampler must branch on :meth:`isBatch` as well, or it silently
        delivers one job per epoch.
        """
        return t in (ProcessType.MMAPT, ProcessType.MPHT, ProcessType.BMMAPT)

    @staticmethod
    def isBatch(t):
        """True when an EVENT of this process releases (or, as a service process,
        completes) a BATCH of jobs whose size the process itself carries in its
        blocks.

        This is the batch twin of :meth:`isMarkedStationary` and
        :meth:`isMarkedSchedule`, and the same rule applies: a procid test that
        should serve every batch family is a MEMBERSHIP test, never equality with
        BMAP.

        It is disjoint from ``sn.arrivalbatch``, which is a SEPARATE batch-size
        law bolted onto a renewal stream by ``Source.set_arrival_batch``. A
        process that is isBatch already carries its own sizes, so the two are
        mutually exclusive by construction.
        """
        return t in (ProcessType.BMAP, ProcessType.BMMAPT)

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
            "MAPt": ProcessType.MAPT,
            "PHt": ProcessType.PHT,
            "MarkedMAP": ProcessType.MMAP,
            "MarkedMMPP": ProcessType.MMAP,
            "MMAP": ProcessType.MMAP,
            "MPH": ProcessType.MPH,
            "MarkedPH": ProcessType.MPH,
            "MMAPt": ProcessType.MMAPT,
            "MPHt": ProcessType.MPHT,
            "BMMAPt": ProcessType.BMMAPT,
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


class RoutingStrategy(IntEnum):
    """
    Strategies for routing jobs between network nodes.

    THE single definition. lang/base.py and api/sn/network_struct.py import
    this class; until 2026-09-17 each declared a twin of its own. All three
    agreed on every member AND every value, so the two `IntEnum` twins did
    compare equal to each other -- but this one was a plain `Enum`, which
    compares equal to neither, and it is the copy the package exports
    (`from line_solver import RoutingStrategy`). So a member handed in by a
    caller tested False against every `RoutingStrategy.X` inside lang/ and
    api/sn/, silently.

    Values must match MATLAB's RoutingStrategy constants for JMT compatibility.
    """
    RAND = 0       # Random routing (uniform among destinations)
    PROB = 1       # Probabilistic routing (explicit probabilities)
    RROBIN = 2     # Round-robin
    WRROBIN = 3    # Weighted round-robin
    JSQ = 4        # Join-Shortest-Queue
    FIRING = 5     # Firing (for Petri nets)
    SQ = 6         # KCHOICES is now SQ: shortest queue of d, SQ(d)
    SDR = 7        # Krzesinski (1987) product-form state-dependent routing
    DISABLED = -1  # Disabled routing

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


class SchedStrategy(IntEnum):
    """
    Scheduling strategies for service stations.

    THE single definition. lang/base.py and api/sn/network_struct.py import
    this class; until 2026-09-16 each declared a twin of its own and seven
    members carried different integers here than there (PSJF, FB, LAS, LRPT,
    FSP, PAS, OI), so an `Enum` member from this module never compared equal
    to the `IntEnum` member of the same name. The numbering kept below is the
    one lang/base.py and network_struct.py already agreed on, because that is
    the one every FSP/PAS/OI caller was reading.

    Declared in ascending value order on purpose: `toID` is positional.
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
    FSP = 36       # Fair Sojourn Protocol (virtual-PS finish time ranking)
    PAS = 37       # Pass-and-swap (order-independent queue with class swap graph)
    OI = 38        # Order-independent (pass-and-swap specialization with empty/zero swap graph)
    PSJF = 39      # Preemptive Shortest Job First
    FB = 40        # Foreground-Background
    LAS = 41       # Least Attained Service (the FB alias)
    LRPT = 42      # Longest Remaining Processing Time
    # appended, not slotted next to FCFSPR: 43 was the first id free in every
    # Python SchedStrategy when FCFSPI was added, and toID below is positional
    FCFSPI = 43    # FCFS Preemptive Identical

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


class SchedStrategyType(IntEnum):
    """
    Categories of scheduling strategies by preemption behavior.

    THE single definition; lang/base.py imports it. Both the values and the
    member set are MATLAB's (`matlab/src/lang/constant/SchedStrategyType.m`),
    which is ground truth: this module used to declare `auto()` values in the
    JAR's declaration order (PR=1, PNR=2, NP=3, NPPrio=4) and lang/base.py a
    two-member `IntEnum` (NP=0, PR=1), so the two disagreed on every value they
    shared. `lang/base.py`'s numbering was the MATLAB one and is kept.
    """
    NP = 0        # Non-preemptive
    PR = 1        # Preemptive resume
    PNR = 2       # Preemptive non-resume
    NPPrio = 3    # Non-preemptive priority
    PRPrio = 4    # Preemptive resume priority
    PNRPrio = 5   # Preemptive non-resume priority

    @staticmethod
    def get_type_id(strategy) -> 'SchedStrategyType':
        """Classify a scheduling discipline by preemption and priority.

        Port of MATLAB `SchedStrategyType.getTypeId`, and the classification
        the Queue constructor stores in `_sched_policy`. The three codebases
        keep ONE table each and all three call it: MATLAB's `Queue.m`, the
        JAR's `Queue.java` and `Queue` below all ask this question rather than
        answering it, because until 2026-09-17 they answered it separately and
        disagreed on eleven disciplines.

        `RAND` and `LAS` are classified here but not in MATLAB, which has no
        RAND member and spells LAS as the same id as FB; native python declares
        all four, so both aliases need naming.

        Raises ValueError on a discipline no Queue accepts (the node markers
        EXT, FORK, JOIN and REF), as MATLAB raises there.
        """
        name = getattr(strategy, 'name', None)
        if name is None:
            try:
                name = SchedStrategy(int(strategy)).name
            except (ValueError, TypeError):
                name = str(strategy).split('.')[-1].upper()
        policy = _SCHED_POLICY_BY_NAME.get(name)
        if policy is None:
            raise ValueError("Unrecognized scheduling strategy type: %s" % (name,))
        return policy

    @staticmethod
    def to_text(type_) -> str:
        """Port of MATLAB `SchedStrategyType.toText`."""
        text = _SCHED_POLICY_TEXT.get(getattr(type_, 'name', None))
        if text is None:
            raise ValueError("Unrecognized scheduling strategy type: %s" % (type_,))
        return text



# The classification stored in schedPolicy, keyed by NAME as everything that
# crosses a codebase boundary is. MATLAB's SchedStrategyType.getTypeId and the
# JAR's are the same table, and all three Queue constructors call theirs rather
# than repeating it. Three separate tables is what let them disagree on eleven
# disciplines: MATLAB answered PR for the preemptive-identical family and NP
# for PSJF/FB/LRPT, and neither MATLAB nor the JAR used PRPrio or PNRPrio at
# all. Nothing in any codebase READS schedPolicy, so the disagreement was never
# a wrong number, only a wrong label, which is why it survived so long.
#
# The preemption axis follows each discipline's own definition in SchedStrategy:
# PI is "preemptive independent", i.e. the preempted job RESTARTS, so it is
# non-resume and not resume; SRPT, PSJF, FB, LRPT, FSP and EDF preempt, against
# SJF, LJF, SEPT, LEPT, EDD and SETF, which rank the same jobs without
# preempting. A discipline carrying a priority order gets its Prio type.
_SCHED_POLICY_BY_NAME = {}

for _name in ('INF', 'FCFS', 'LCFS', 'SIRO', 'RAND', 'SJF', 'LJF', 'SEPT',
              'LEPT', 'EDD', 'SETF', 'PAS', 'OI', 'POLLING'):
    _SCHED_POLICY_BY_NAME[_name] = SchedStrategyType.NP

for _name in ('PS', 'DPS', 'GPS', 'LPS', 'LCFSPR', 'FCFSPR', 'EDF', 'SRPT',
              'PSJF', 'FB', 'LAS', 'LRPT', 'FSP'):
    _SCHED_POLICY_BY_NAME[_name] = SchedStrategyType.PR

for _name in ('LCFSPI', 'FCFSPI'):
    _SCHED_POLICY_BY_NAME[_name] = SchedStrategyType.PNR

for _name in ('HOL', 'FCFSPRIO', 'LCFSPRIO'):
    _SCHED_POLICY_BY_NAME[_name] = SchedStrategyType.NPPrio

for _name in ('PSPRIO', 'DPSPRIO', 'GPSPRIO', 'LCFSPRPRIO', 'FCFSPRPRIO',
              'SRPTPRIO'):
    _SCHED_POLICY_BY_NAME[_name] = SchedStrategyType.PRPrio

for _name in ('LCFSPIPRIO', 'FCFSPIPRIO'):
    _SCHED_POLICY_BY_NAME[_name] = SchedStrategyType.PNRPrio

del _name

_SCHED_POLICY_TEXT = {
    'NP': 'NonPreemptive',
    'PR': 'PreemptiveResume',
    'PNR': 'PreemptiveNonResume',
    'NPPrio': 'NonPreemptivePriority',
    'PRPrio': 'PreemptiveResumePriority',
    'PNRPrio': 'PreemptiveNonResumePriority',
}


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
    SSA = auto()


# TimingStrategy is NOT defined here. It lives in line_solver.lang.nodes, whose
# values (TIMED = 0, IMMEDIATE = 1) are the ones MATLAB and the JAR use and the
# ones every mode comparison in the model layer is written against. A second enum
# here declared TIMED = 1, IMMEDIATE = 2 through auto(), so a transition set from
# it never compared equal to an immediate mode and was silently served as timed.
# Import it from line_solver.lang.nodes (or from the package root).


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


class HeteroSchedPolicy(IntEnum):
    """
    Scheduling policies for heterogeneous multiserver queues.

    THE single definition; lang/base.py imports it. Values are MATLAB's
    (`matlab/src/lang/constant/HeteroSchedPolicy.m`), which is ground truth.
    This module used to declare a plain `Enum` twin with `auto()` values
    (ORDER=1..RAIS=6) that agreed with nothing; it was shadowed at package
    level by the copy below and so was never the one anybody got.

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


class GlobalConstants:
    """
    Global constants and configuration for the LINE solver.
    """
    Zero = 1e-14
    # Magnitude above which an off-diagonal generator entry counts as an arc.
    # Sign is NOT a criterion: an ME generator embeds genuinely negative
    # off-diagonal entries -- see _kb/11-conventions-and-gotchas.md
    ArcTol = 1e-12
    CoarseTol = 1e-3  # Match MATLAB/JAR default (1.0e-03)
    FineTol = 1e-8  # Match MATLAB's default
    Immediate = 1e8  # 1/FineTol - large but finite rate for immediate service (matches MATLAB)
    MaxInt = 2**31 - 1
    Version = "3.0.8"
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
            'ArcTol': cls.ArcTol,
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
