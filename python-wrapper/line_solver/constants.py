
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

from enum_tools import Enum
import jpype
import jpype.imports
from line_solver.distributions import *

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
    
    def __repr__(self):
        return str(self.value)
    PRE_SEQ = jpype.JPackage('jline').lang.constant.ActivityPrecedenceType.PRE_SEQ
    PRE_AND = jpype.JPackage('jline').lang.constant.ActivityPrecedenceType.PRE_AND
    PRE_OR = jpype.JPackage('jline').lang.constant.ActivityPrecedenceType.PRE_OR
    POST_SEQ = jpype.JPackage('jline').lang.constant.ActivityPrecedenceType.POST_SEQ
    POST_AND = jpype.JPackage('jline').lang.constant.ActivityPrecedenceType.POST_AND
    POST_OR = jpype.JPackage('jline').lang.constant.ActivityPrecedenceType.POST_OR
    POST_LOOP = jpype.JPackage('jline').lang.constant.ActivityPrecedenceType.POST_LOOP
    POST_CACHE = jpype.JPackage('jline').lang.constant.ActivityPrecedenceType.POST_CACHE

class CallType(Enum):
    """
    Types of calls between tasks in layered networks.
    
    - SYNC: Synchronous call (caller waits for response)
    - ASYNC: Asynchronous call (caller continues immediately)
    - FWD: Forward call (caller terminates, response goes to caller's caller)
    """
    
    def __repr__(self):
        return str(self.value)
    SYNC = jpype.JPackage('jline').lang.constant.CallType.SYNC
    ASYNC = jpype.JPackage('jline').lang.constant.CallType.ASYNC
    FWD = jpype.JPackage('jline').lang.constant.CallType.FWD

class DropStrategy(Enum):
    """
    Strategies for handling queue overflow and capacity limits.

    - WaitingQueue: Jobs wait in a waiting queue when capacity is exceeded
    - Drop: Jobs are dropped (lost) when capacity is exceeded
    - BlockingAfterService: Jobs are blocked after service completion
    """
    def __repr__(self):
        return str(self.value)
    WaitingQueue = jpype.JPackage('jline').lang.constant.DropStrategy.WaitingQueue
    Drop = jpype.JPackage('jline').lang.constant.DropStrategy.Drop
    BlockingAfterService = jpype.JPackage('jline').lang.constant.DropStrategy.BlockingAfterService

class SignalType(Enum):
    """
    Types of signals for signal classes in G-networks and related models.

    Signal classes are specialized open classes that can have special effects
    on queues they visit, such as removing jobs (negative signals).

    Attributes:
        NEGATIVE: Removes a job from the destination queue (G-network negative customer)
        REPLY: Triggers a reply action
        CATASTROPHE: Removes ALL jobs from the destination queue

    Reference: Gelenbe, E. (1991). "Product-form queueing networks with
               negative and positive customers", Journal of Applied Probability
    """
    def __repr__(self):
        return str(self.value)
    NEGATIVE = jpype.JPackage('jline').lang.constant.SignalType.NEGATIVE
    REPLY = jpype.JPackage('jline').lang.constant.SignalType.REPLY
    CATASTROPHE = jpype.JPackage('jline').lang.constant.SignalType.CATASTROPHE


class RemovalPolicy(Enum):
    """
    Removal policies for negative signals in G-networks.

    When a negative signal removes jobs from a queue, the removal policy
    determines which job(s) are selected for removal.

    Attributes:
        RANDOM: Select job uniformly at random from all jobs at the station
        FCFS: Remove the oldest job (first arrived)
        LCFS: Remove the newest job (last arrived)

    Example:
        >>> from line_solver import Signal, SignalType, RemovalPolicy, Geometric
        >>> neg_class = Signal(model, "Negative", SignalType.NEGATIVE)
        >>> neg_class.removal_policy = RemovalPolicy.LCFS
        >>> neg_class.removal_distribution = Geometric(0.5)  # Mean 2 removals
    """
    def __repr__(self):
        return str(self.value)
    RANDOM = jpype.JPackage('jline').lang.constant.RemovalPolicy.RANDOM
    FCFS = jpype.JPackage('jline').lang.constant.RemovalPolicy.FCFS
    LCFS = jpype.JPackage('jline').lang.constant.RemovalPolicy.LCFS

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
    def __repr__(self):
        return str(self.value)
    INIT = jpype.JPackage('jline').lang.constant.EventType.INIT
    LOCAL = jpype.JPackage('jline').lang.constant.EventType.LOCAL
    ARV = jpype.JPackage('jline').lang.constant.EventType.ARV
    DEP = jpype.JPackage('jline').lang.constant.EventType.DEP
    PHASE = jpype.JPackage('jline').lang.constant.EventType.PHASE
    READ = jpype.JPackage('jline').lang.constant.EventType.READ
    STAGE = jpype.JPackage('jline').lang.constant.EventType.STAGE

class JobClassType(Enum):
    """
    Types of job classes in queueing networks.
    
    - OPEN: Open class (jobs arrive from outside the system)
    - CLOSED: Closed class (fixed population circulating in the system)
    - DISABLED: Disabled class (not currently active)
    """
    
    def __repr__(self):
        return str(self.value)
    OPEN = jpype.JPackage('jline').lang.constant.JobClassType.OPEN
    CLOSED = jpype.JPackage('jline').lang.constant.JobClassType.CLOSED
    DISABLED = jpype.JPackage('jline').lang.constant.JobClassType.DISABLED

class JoinStrategy(Enum):
    """
    Strategies for join node synchronization in fork-join networks.
    
    - STD: Standard join (wait for all parallel branches)
    - PARTIAL: Partial join (proceed when some branches complete)
    - Quorum: Quorum-based join (wait for minimum number of branches)
    - Guard: Guard condition join (custom completion criteria)
    """
    def __repr__(self):
        return str(self.value)
    STD = jpype.JPackage('jline').lang.constant.JoinStrategy.STD
    PARTIAL = jpype.JPackage('jline').lang.constant.JoinStrategy.PARTIAL
    Quorum = jpype.JPackage('jline').lang.constant.JoinStrategy.Quorum
    Guard = jpype.JPackage('jline').lang.constant.JoinStrategy.Guard

class MetricType(Enum):
    """
    Types of performance metrics that can be computed.

    Includes response times, queue lengths, throughput, utilization,
    drop rates, and system-level metrics for comprehensive performance analysis.

    Attributes:
        ResidT: Residence time (total time in system)
        RespT: Response time (time from arrival to departure)
        DropRate: Rate of dropped/blocked jobs
        QLen: Queue length (number of jobs waiting)
        QueueT: Queueing time (time spent waiting)
        FCRWeight: FCR (Fair Chance Routing) weight metric
        FCRMemOcc: FCR memory occupancy metric
        FJQLen: Fork-Join number of customers
        FJRespT: Fork-Join response time
        RespTSink: Response time measured at sink
        SysQLen: System-wide queue length
        SysRespT: System-wide response time
        SysTput: System-wide throughput
        Tput: Throughput (jobs per unit time)
        ArvR: Arrival rate
        TputSink: Throughput measured at sink
        Util: Utilization (fraction of time busy)
        TranQLen: Transient queue length
        TranUtil: Transient utilization
        TranTput: Transient throughput
        TranRespT: Transient response time
        Tard: Tardiness (time beyond soft deadline)
        SysTard: System-wide tardiness
    """
    def __repr__(self):
        return str(self.value)
    ResidT = jpype.JPackage('jline').lang.constant.MetricType.ResidT
    RespT = jpype.JPackage('jline').lang.constant.MetricType.RespT
    DropRate = jpype.JPackage('jline').lang.constant.MetricType.DropRate
    QLen = jpype.JPackage('jline').lang.constant.MetricType.QLen
    QueueT = jpype.JPackage('jline').lang.constant.MetricType.QueueT
    FCRWeight = jpype.JPackage('jline').lang.constant.MetricType.FCRWeight
    FCRMemOcc = jpype.JPackage('jline').lang.constant.MetricType.FCRMemOcc
    FJQLen = jpype.JPackage('jline').lang.constant.MetricType.FJQLen
    FJRespT = jpype.JPackage('jline').lang.constant.MetricType.FJRespT
    RespTSink = jpype.JPackage('jline').lang.constant.MetricType.RespTSink
    SysQLen = jpype.JPackage('jline').lang.constant.MetricType.SysQLen
    SysRespT = jpype.JPackage('jline').lang.constant.MetricType.SysRespT
    SysTput = jpype.JPackage('jline').lang.constant.MetricType.SysTput
    Tput = jpype.JPackage('jline').lang.constant.MetricType.Tput
    ArvR = jpype.JPackage('jline').lang.constant.MetricType.ArvR
    TputSink = jpype.JPackage('jline').lang.constant.MetricType.TputSink
    Util = jpype.JPackage('jline').lang.constant.MetricType.Util
    TranQLen = jpype.JPackage('jline').lang.constant.MetricType.TranQLen
    TranUtil = jpype.JPackage('jline').lang.constant.MetricType.TranUtil
    TranTput = jpype.JPackage('jline').lang.constant.MetricType.TranTput
    TranRespT = jpype.JPackage('jline').lang.constant.MetricType.TranRespT
    Tard = jpype.JPackage('jline').lang.constant.MetricType.Tard
    SysTard = jpype.JPackage('jline').lang.constant.MetricType.SysTard

class NodeType(Enum):
    """
    Types of nodes in queueing network models.
    
    Defines the various types of nodes that can be used to construct
    queueing networks, from basic service stations to advanced modeling
    elements like caches and fork-join synchronization.
    
    Attributes:
        Source: External arrival source for generating jobs
        Queue: Service station with waiting area and server
        Sink: External departure destination for completed jobs
        Delay: Infinite server (no waiting, immediate service)
        Fork: Splits jobs into parallel tasks for concurrent processing
        Join: Synchronizes parallel tasks back into single job
        Router: Probabilistic routing decisions between paths
        Cache: Cache model with hit/miss behavior and replacement policies
        Logger: Event logging and monitoring for performance analysis
        ClassSwitch: Job class switching for multi-class networks
        Transition: Petri net transition element
        Place: Petri net place element
    """
    def __repr__(self):
        return str(self.value)
    Transition = jpype.JPackage('jline').lang.constant.NodeType.Transition
    Place = jpype.JPackage('jline').lang.constant.NodeType.Place
    Fork = jpype.JPackage('jline').lang.constant.NodeType.Fork
    Router = jpype.JPackage('jline').lang.constant.NodeType.Router
    Cache = jpype.JPackage('jline').lang.constant.NodeType.Cache
    Logger = jpype.JPackage('jline').lang.constant.NodeType.Logger
    ClassSwitch = jpype.JPackage('jline').lang.constant.NodeType.ClassSwitch
    Delay = jpype.JPackage('jline').lang.constant.NodeType.Delay
    Source = jpype.JPackage('jline').lang.constant.NodeType.Source
    Sink = jpype.JPackage('jline').lang.constant.NodeType.Sink
    Join = jpype.JPackage('jline').lang.constant.NodeType.Join
    Queue = jpype.JPackage('jline').lang.constant.NodeType.Queue

    @staticmethod
    def fromJLine(obj):
        obj_str = str(obj)
        if obj_str == 'Transition':
            return NodeType.Transition
        elif obj_str == 'Place':
            return NodeType.Place
        elif obj_str == 'Fork':
            return NodeType.Fork
        elif obj_str == 'Router':
            return NodeType.Router
        elif obj_str == 'Cache':
            return NodeType.Cache
        elif obj_str == 'Logger':
            return NodeType.Logger
        elif obj_str == 'ClassSwitch':
            return NodeType.ClassSwitch
        elif obj_str == 'Delay':
            return NodeType.Delay
        elif obj_str == 'Source':
            return NodeType.Source
        elif obj_str == 'Sink':
            return NodeType.Sink
        elif obj_str == 'Join':
            return NodeType.Join
        elif obj_str == 'Queue':
            return NodeType.Queue

class ProcessType(Enum):
    """
    Types of stochastic processes for arrivals and service times.
    
    Comprehensive collection of probability distributions and processes
    for modeling arrivals, service times, and other random phenomena in
    queueing networks.
    
    Attributes:
        EXP: Exponential distribution (memoryless, rate parameter)
        ERLANG: Erlang distribution (sum of exponentials)
        HYPEREXP: Hyperexponential distribution (mixture of exponentials)
        GAMMA: Gamma distribution (shape and scale parameters)
        UNIFORM: Uniform distribution (constant probability over interval)
        DET: Deterministic distribution (constant value, no randomness)
        PARETO: Pareto distribution (heavy-tailed, power law)
        WEIBULL: Weibull distribution (reliability and failure modeling)
        LOGNORMAL: Log-normal distribution (log of values is normal)
        PH: Phase-type distribution (Markovian representation)
        APH: Acyclic phase-type distribution (acyclic PH)
        COXIAN: Coxian distribution (series of exponential phases)
        COX2: Two-phase Coxian distribution
        MAP: Markovian Arrival Process (correlated arrivals)
        MMPP2: Two-state Markov Modulated Poisson Process
        BINOMIAL: Binomial distribution (discrete, n trials)
        POISSON: Poisson distribution (discrete, rate parameter)
        TRACE: Empirical trace distribution (replay recorded data)
        REPLAYER: Trace replayer for empirical data
        IMMEDIATE: Immediate transition (zero delay)
        DISABLED: Disabled process (no arrivals/service)
    """
    def __repr__(self):
        return str(self.value)

    EXP = jpype.JPackage('jline').lang.constant.ProcessType.EXP
    ERLANG = jpype.JPackage('jline').lang.constant.ProcessType.ERLANG
    DISABLED = jpype.JPackage('jline').lang.constant.ProcessType.DISABLED
    IMMEDIATE = jpype.JPackage('jline').lang.constant.ProcessType.IMMEDIATE
    HYPEREXP = jpype.JPackage('jline').lang.constant.ProcessType.HYPEREXP
    APH = jpype.JPackage('jline').lang.constant.ProcessType.APH
    COXIAN = jpype.JPackage('jline').lang.constant.ProcessType.COXIAN
    PH = jpype.JPackage('jline').lang.constant.ProcessType.PH
    MAP = jpype.JPackage('jline').lang.constant.ProcessType.MAP
    UNIFORM = jpype.JPackage('jline').lang.constant.ProcessType.UNIFORM
    DET = jpype.JPackage('jline').lang.constant.ProcessType.DET
    GAMMA = jpype.JPackage('jline').lang.constant.ProcessType.GAMMA
    PARETO = jpype.JPackage('jline').lang.constant.ProcessType.PARETO
    WEIBULL = jpype.JPackage('jline').lang.constant.ProcessType.WEIBULL
    LOGNORMAL = jpype.JPackage('jline').lang.constant.ProcessType.LOGNORMAL
    MMPP2 = jpype.JPackage('jline').lang.constant.ProcessType.MMPP2
    REPLAYER = jpype.JPackage('jline').lang.constant.ProcessType.REPLAYER
    TRACE = jpype.JPackage('jline').lang.constant.ProcessType.TRACE
    COX2 = jpype.JPackage('jline').lang.constant.ProcessType.COX2
    BINOMIAL = jpype.JPackage('jline').lang.constant.ProcessType.BINOMIAL
    POISSON = jpype.JPackage('jline').lang.constant.ProcessType.POISSON

    @staticmethod
    def fromString(obj):
        obj_str = str(obj)
        if obj_str == "Exp":
            return ProcessType.EXP
        elif obj_str == "Erlang":
            return ProcessType.ERLANG
        elif obj_str == "HyperExp":
            return ProcessType.HYPEREXP
        elif obj_str == "PH":
            return ProcessType.PH
        elif obj_str == "APH":
            return ProcessType.APH
        elif obj_str == "MAP":
            return ProcessType.MAP
        elif obj_str == "Uniform":
            return ProcessType.UNIFORM
        elif obj_str == "Det":
            return ProcessType.DET
        elif obj_str == "Coxian":
            return ProcessType.COXIAN
        elif obj_str == "Gamma":
            return ProcessType.GAMMA
        elif obj_str == "Pareto":
            return ProcessType.PARETO
        elif obj_str == "MMPP2":
            return ProcessType.MMPP2
        elif obj_str == "Replayer":
            return ProcessType.REPLAYER
        elif obj_str == "Trace":
            return ProcessType.TRACE
        elif obj_str == "Immediate":
            return ProcessType.IMMEDIATE
        elif obj_str == "Disabled":
            return ProcessType.DISABLED
        elif obj_str == "Cox2":
            return ProcessType.COX2
        elif obj_str == "Weibull":
            return ProcessType.WEIBULL
        elif obj_str == "Lognormal":
            return ProcessType.LOGNORMAL
        elif obj_str == "Poisson":
            return ProcessType.POISSON
        elif obj_str == "Binomial":
            return ProcessType.BINOMIAL
        else:
            raise ValueError(f"Unsupported ProcessType: {obj}")

    def toDistribution(process_type, *args):
        if process_type == ProcessType.EXP:
            return Exp
        elif process_type == ProcessType.DET:
            return Det
        elif process_type == ProcessType.ERLANG:
            return Erlang
        elif process_type == ProcessType.HYPEREXP:
            return HyperExp
        elif process_type == ProcessType.UNIFORM:
            return Uniform
        elif process_type == ProcessType.IMMEDIATE:
            return Immediate
        elif process_type == ProcessType.DISABLED:
            return Disabled
        else:
            raise ValueError(f"Unsupported ProcessType: {process_type}")

class ReplacementStrategy(Enum):
    """
    Cache replacement strategies for evicting items.
    
    Defines policies for selecting which cached items to evict when
    the cache reaches capacity and new items need to be stored.
    
    Attributes:
        RR: Round Robin replacement (cyclical eviction)
        FIFO: First In, First Out replacement (oldest item evicted)
        SFIFO: Segmented FIFO replacement (partitioned FIFO policy)
        LRU: Least Recently Used replacement (evict least recently accessed)
    """
    def __repr__(self):
        return str(self.value)
    RR = jpype.JPackage('jline').lang.constant.ReplacementStrategy.RR
    FIFO = jpype.JPackage('jline').lang.constant.ReplacementStrategy.FIFO
    SFIFO = jpype.JPackage('jline').lang.constant.ReplacementStrategy.SFIFO
    LRU = jpype.JPackage('jline').lang.constant.ReplacementStrategy.LRU

    def ordinal(self):
        return self.value.ordinal()

class RoutingStrategy(Enum):
    """
    Strategies for routing jobs between network nodes.
    
    Defines how jobs are routed from one node to another in queueing
    networks, supporting load balancing and performance optimization.
    
    Attributes:
        RAND: Random routing (uniform selection among destinations)
        PROB: Probabilistic routing based on routing matrix
        RROBIN: Round robin routing (cyclical destination selection)
        WRROBIN: Weighted round robin routing (weighted cyclical selection)
        JSQ: Join Shortest Queue routing (route to least congested)
        KCHOICES: Power-of-k-choices routing (select best of k random)
        DISABLED: Disabled routing (no routing, for Petri nets)
        FIRING: Firing-based routing (for Petri net transitions)
    """
    def __repr__(self):
        return str(self.value)
    RAND = jpype.JPackage('jline').lang.constant.RoutingStrategy.RAND
    PROB = jpype.JPackage('jline').lang.constant.RoutingStrategy.PROB
    RROBIN = jpype.JPackage('jline').lang.constant.RoutingStrategy.RROBIN
    WRROBIN = jpype.JPackage('jline').lang.constant.RoutingStrategy.WRROBIN
    JSQ = jpype.JPackage('jline').lang.constant.RoutingStrategy.JSQ
    DISABLED = jpype.JPackage('jline').lang.constant.RoutingStrategy.DISABLED
    FIRING = jpype.JPackage('jline').lang.constant.RoutingStrategy.FIRING
    KCHOICES = jpype.JPackage('jline').lang.constant.RoutingStrategy.KCHOICES

class SchedStrategy(Enum):
    """
    Scheduling strategies for service stations.
    
    Comprehensive collection of queueing disciplines and scheduling policies
    for determining the order in which jobs receive service at stations.
    
    Attributes:
        INF: Infinite servers (delay station, no waiting)
        FCFS: First Come, First Served (arrival order service)
        LCFS: Last Come, First Served (stack-like service)
        LCFSPR: Last Come, First Served with Preemptive Resume
        SIRO: Service in Random Order
        SJF: Shortest Job First (service time based)
        LJF: Longest Job First (service time based)
        EDD: Earliest Due Date (deadline-based, non-preemptive)
        EDF: Earliest Deadline First (deadline-based, preemptive)
        PS: Processor Sharing (equal time sharing among jobs)
        DPS: Discriminatory Processor Sharing (weighted sharing)
        GPS: Generalized Processor Sharing (generalized weighted sharing)
        SEPT: Shortest Expected Processing Time
        LEPT: Longest Expected Processing Time
        SRPT: Shortest Remaining Processing Time
        SRPTPRIO: Shortest Remaining Processing Time with Priority
        HOL: Head of Line priority scheduling
        FORK: Fork node scheduling (job splitting)
        EXT: External scheduling (user-defined)
        REF: Reference node scheduling
        POLLING: Polling-based scheduling
        PSPRIO: Processor Sharing with Priority
        DPSPRIO: Discriminatory Processor Sharing with Priority
        GPSPRIO: Generalized Processor Sharing with Priority
    """

    def __repr__(self):
        return str(self.value)

    INF = jpype.JPackage('jline').lang.constant.SchedStrategy.INF
    FCFS = jpype.JPackage('jline').lang.constant.SchedStrategy.FCFS
    LCFS = jpype.JPackage('jline').lang.constant.SchedStrategy.LCFS
    LCFSPR = jpype.JPackage('jline').lang.constant.SchedStrategy.LCFSPR
    SIRO = jpype.JPackage('jline').lang.constant.SchedStrategy.SIRO
    SJF = jpype.JPackage('jline').lang.constant.SchedStrategy.SJF
    LJF = jpype.JPackage('jline').lang.constant.SchedStrategy.LJF
    EDD = jpype.JPackage('jline').lang.constant.SchedStrategy.EDD
    EDF = jpype.JPackage('jline').lang.constant.SchedStrategy.EDF
    PS = jpype.JPackage('jline').lang.constant.SchedStrategy.PS
    DPS = jpype.JPackage('jline').lang.constant.SchedStrategy.DPS
    GPS = jpype.JPackage('jline').lang.constant.SchedStrategy.GPS
    SEPT = jpype.JPackage('jline').lang.constant.SchedStrategy.SEPT
    LEPT = jpype.JPackage('jline').lang.constant.SchedStrategy.LEPT
    SRPT = jpype.JPackage('jline').lang.constant.SchedStrategy.SRPT
    SRPTPRIO = jpype.JPackage('jline').lang.constant.SchedStrategy.SRPTPRIO
    HOL = jpype.JPackage('jline').lang.constant.SchedStrategy.HOL
    FORK = jpype.JPackage('jline').lang.constant.SchedStrategy.FORK
    EXT = jpype.JPackage('jline').lang.constant.SchedStrategy.EXT
    REF = jpype.JPackage('jline').lang.constant.SchedStrategy.REF
    POLLING = jpype.JPackage('jline').lang.constant.SchedStrategy.POLLING
    PSPRIO = jpype.JPackage('jline').lang.constant.SchedStrategy.PSPRIO
    DPSPRIO = jpype.JPackage('jline').lang.constant.SchedStrategy.DPSPRIO
    GPSPRIO = jpype.JPackage('jline').lang.constant.SchedStrategy.GPSPRIO
    LCFSPRIO = jpype.JPackage('jline').lang.constant.SchedStrategy.LCFSPRIO
    LCFSPRPRIO = jpype.JPackage('jline').lang.constant.SchedStrategy.LCFSPRPRIO
    LCFSPIPRIO = jpype.JPackage('jline').lang.constant.SchedStrategy.LCFSPIPRIO
    FCFSPRPRIO = jpype.JPackage('jline').lang.constant.SchedStrategy.FCFSPRPRIO
    FCFSPIPRIO = jpype.JPackage('jline').lang.constant.SchedStrategy.FCFSPIPRIO
    FCFSPRIO = jpype.JPackage('jline').lang.constant.SchedStrategy.FCFSPRIO
    PSJF = jpype.JPackage('jline').lang.constant.SchedStrategy.PSJF
    FB = jpype.JPackage('jline').lang.constant.SchedStrategy.FB
    LAS = jpype.JPackage('jline').lang.constant.SchedStrategy.FB
    LRPT = jpype.JPackage('jline').lang.constant.SchedStrategy.LRPT
    SETF = jpype.JPackage('jline').lang.constant.SchedStrategy.SETF

    @staticmethod
    def fromString(obj):
        obj_str = str(obj)
        if obj_str == 'INF':
            return SchedStrategy.INF
        elif obj_str == 'FCFS':
            return SchedStrategy.FCFS
        elif obj_str == 'LCFS':
            return SchedStrategy.LCFS
        elif obj_str == 'LCFSPR':
            return SchedStrategy.LCFSPR
        elif obj_str == 'SIRO':
            return SchedStrategy.SIRO
        elif obj_str == 'SJF':
            return SchedStrategy.SJF
        elif obj_str == 'LJF':
            return SchedStrategy.LJF
        elif obj_str == 'EDD':
            return SchedStrategy.EDD
        elif obj_str == 'EDF':
            return SchedStrategy.EDF
        elif obj_str == 'PS':
            return SchedStrategy.PS
        elif obj_str == 'DPS':
            return SchedStrategy.DPS
        elif obj_str == 'GPS':
            return SchedStrategy.GPS
        elif obj_str == 'SEPT':
            return SchedStrategy.SEPT
        elif obj_str == 'LEPT':
            return SchedStrategy.LEPT
        elif obj_str == 'SRPT':
            return SchedStrategy.SRPT
        elif obj_str == 'SRPTPRIO':
            return SchedStrategy.SRPTPRIO
        elif obj_str == 'HOL':
            return SchedStrategy.HOL
        elif obj_str == 'FORK':
            return SchedStrategy.FORK
        elif obj_str == 'EXT':
            return SchedStrategy.EXT
        elif obj_str == 'REF':
            return SchedStrategy.REF
        elif obj_str == 'POLLING':
            return SchedStrategy.POLLING
        elif obj_str == 'PSPRIO':
            return SchedStrategy.PSPRIO
        elif obj_str == 'DPSPRIO':
            return SchedStrategy.DPSPRIO
        elif obj_str == 'GPSPRIO':
            return SchedStrategy.GPSPRIO
        elif obj_str == 'LCFSPRIO':
            return SchedStrategy.LCFSPRIO
        elif obj_str == 'LCFSPRPRIO':
            return SchedStrategy.LCFSPRPRIO
        elif obj_str == 'LCFSPIPRIO':
            return SchedStrategy.LCFSPIPRIO
        elif obj_str == 'FCFSPRPRIO':
            return SchedStrategy.FCFSPRPRIO
        elif obj_str == 'FCFSPIPRIO':
            return SchedStrategy.FCFSPIPRIO
        elif obj_str == 'FCFSPRIO':
            return SchedStrategy.FCFSPRIO
        elif obj_str == 'PSJF':
            return SchedStrategy.PSJF
        elif obj_str == 'FB':
            return SchedStrategy.FB
        elif obj_str == 'LAS':
            return SchedStrategy.LAS
        elif obj_str == 'LRPT':
            return SchedStrategy.LRPT
        elif obj_str == 'SETF':
            return SchedStrategy.SETF
        else:
            raise ValueError(f"Unsupported SchedStrategy: {obj}")

    @staticmethod
    def fromLINEString(sched: str):
        return jpype.JPackage('jline').lang.constant.SchedStrategy.fromLINEString(sched.lower())

    @staticmethod
    def toID(sched):
        return int(jpype.JPackage('jline').lang.constant.SchedStrategy.toID(sched))

class SchedStrategyType(Enum):
    """
    Categories of scheduling strategies by preemption behavior.
    
    Classifies scheduling strategies based on how they handle job
    interruption and resumption for priority management.
    
    Attributes:
        PR: Preemptive Resume (higher priority jobs preempt, resume later)
        PNR: Preemptive Non-Resume (higher priority jobs preempt, restart)
        NP: Non-Preemptive (jobs cannot be interrupted once started)
        NPPrio: Non-Preemptive with Priority (priority order, no preemption)
    """
    def __repr__(self):
        return str(self.value)
    PR = jpype.JPackage('jline').lang.constant.SchedStrategyType.PR
    PNR = jpype.JPackage('jline').lang.constant.SchedStrategyType.PNR
    NP = jpype.JPackage('jline').lang.constant.SchedStrategyType.NP
    NPPrio = jpype.JPackage('jline').lang.constant.SchedStrategyType.NPPrio

class ServiceStrategy(Enum):
    """
    Service strategies defining service time dependence.
    
    Determines how service times vary based on system conditions,
    load, job characteristics, or other state variables.
    
    Attributes:
        LI: Load Independent (service time independent of queue length)
        LD: Load Dependent (service time depends on queue length)
        CD: Class Dependent (service time depends on job class)
        SD: State Dependent (service time depends on system state)
    """
    def __repr__(self):
        return str(self.value)
    LI = jpype.JPackage('jline').lang.constant.ServiceStrategy.LI
    LD = jpype.JPackage('jline').lang.constant.ServiceStrategy.LD
    CD = jpype.JPackage('jline').lang.constant.ServiceStrategy.CD
    SD = jpype.JPackage('jline').lang.constant.ServiceStrategy.SD

class SolverType(Enum):
    """
    Types of solvers available in LINE.
    
    Comprehensive collection of analysis methods and algorithms
    for solving queueing network models with different techniques.
    
    Attributes:
        AUTO: Automatic solver selection based on model characteristics
        CTMC: Continuous Time Markov Chain solver for exact analysis
        ENV: Ensemble environment solver for random environments and ensembles
        FLUID: FLD - Fluid/Mean-Field Approximation solver for large-scale systems
        JMT: Java Modelling Tools integration for simulation
        LN: Layered Network solver for layered queueing networks
        LQNS: Layered Queueing Network Solver (external tool integration)
        MAM: Matrix Analytic Methods for MAP/PH models
        MVA: Mean Value Analysis for product-form networks
        NC: Normalizing Constant method for closed networks
        QNS: Queueing Network Solver for general networks
        SSA: Stochastic Simulation Algorithms for discrete-event simulation
    """
    def __repr__(self):
        return str(self.value)
    AUTO = jpype.JPackage('jline').lang.constant.SolverType.AUTO
    CTMC = jpype.JPackage('jline').lang.constant.SolverType.CTMC
    ENV = jpype.JPackage('jline').lang.constant.SolverType.ENV
    FLUID = jpype.JPackage('jline').lang.constant.SolverType.FLUID
    JMT = jpype.JPackage('jline').lang.constant.SolverType.JMT
    LN = jpype.JPackage('jline').lang.constant.SolverType.LN
    LQNS = jpype.JPackage('jline').lang.constant.SolverType.LQNS
    MAM = jpype.JPackage('jline').lang.constant.SolverType.MAM
    MVA = jpype.JPackage('jline').lang.constant.SolverType.MVA
    NC = jpype.JPackage('jline').lang.constant.SolverType.NC
    QNS = jpype.JPackage('jline').lang.constant.SolverType.QNS
    SSA = jpype.JPackage('jline').lang.constant.SolverType.SSA

class TimingStrategy(Enum):
    """
    Timing strategies for transitions in Petri nets.
    
    Defines the temporal behavior of Petri net transitions,
    controlling when enabled transitions fire.
    
    Attributes:
        TIMED: Transition has exponentially distributed firing time
        IMMEDIATE: Transition fires immediately when enabled (zero delay)
    """
    def __repr__(self):
        return str(self.value)
    TIMED = jpype.JPackage('jline').lang.constant.TimingStrategy.TIMED
    IMMEDIATE = jpype.JPackage('jline').lang.constant.TimingStrategy.IMMEDIATE

class VerboseLevel(Enum):
    """
    Verbosity levels for LINE solver output.
    
    Controls the amount of information displayed during solver
    execution for monitoring progress and debugging.
    
    Attributes:
        SILENT: No output except errors (minimal logging)
        STD: Standard output level with basic information
        DEBUG: Detailed debug output for troubleshooting and analysis
    """
    def __repr__(self):
        return str(self.value)
    SILENT = jpype.JPackage('jline').VerboseLevel.SILENT
    STD = jpype.JPackage('jline').VerboseLevel.STD
    DEBUG = jpype.JPackage('jline').VerboseLevel.DEBUG

class PollingType(Enum):
    """
    Polling strategies for polling systems.

    Defines how a single server visits and serves multiple queues
    in polling systems with switchover times between queues.

    Attributes:
        GATED: Server polls each queue once per cycle (jobs arriving during service wait)
        EXHAUSTIVE: Server empties each queue before moving to next (serves all present jobs)
        KLIMITED: Server serves at most k jobs from each queue per visit
    """
    def __repr__(self):
        return str(self.value)
    GATED = jpype.JPackage('jline').lang.constant.PollingType.GATED
    EXHAUSTIVE = jpype.JPackage('jline').lang.constant.PollingType.EXHAUSTIVE
    KLIMITED = jpype.JPackage('jline').lang.constant.PollingType.KLIMITED

    @staticmethod
    def fromString(obj):
        obj_str = str(obj).upper()
        if obj_str == 'GATED':
            return PollingType.GATED
        elif obj_str == 'EXHAUSTIVE':
            return PollingType.EXHAUSTIVE
        elif obj_str == 'KLIMITED' or obj_str == 'K-LIMITED':
            return PollingType.KLIMITED
        else:
            raise ValueError(f"Unsupported PollingType: {obj}")

class HeteroSchedPolicy(Enum):
    """
    Scheduling policies for heterogeneous multiserver queues.

    Defines how jobs are assigned to server types in heterogeneous
    multiserver queues where different server types may have different
    speeds and class compatibilities.

    Attributes:
        ORDER: Order-based assignment (servers assigned by type order)
        ALIS: Assign to Longest Idle Server (favors servers idle longest)
        ALFS: Assign to Least Frequently Selected (load balancing)
        FAIRNESS: Fairness-based assignment (balanced utilization)
        FSF: Fastest Server First (prefer faster server types)
        RAIS: Random Assignment with Idle Selection
    """
    def __repr__(self):
        return str(self.value)
    ORDER = jpype.JPackage('jline').lang.constant.HeteroSchedPolicy.ORDER
    ALIS = jpype.JPackage('jline').lang.constant.HeteroSchedPolicy.ALIS
    ALFS = jpype.JPackage('jline').lang.constant.HeteroSchedPolicy.ALFS
    FAIRNESS = jpype.JPackage('jline').lang.constant.HeteroSchedPolicy.FAIRNESS
    FSF = jpype.JPackage('jline').lang.constant.HeteroSchedPolicy.FSF
    RAIS = jpype.JPackage('jline').lang.constant.HeteroSchedPolicy.RAIS

    @staticmethod
    def fromString(obj):
        """Convert a string to HeteroSchedPolicy enum value."""
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
            raise ValueError(f"Unsupported HeteroSchedPolicy: {obj}")

class GlobalConstants:
    """
    Global constants and configuration for the LINE solver.
    
    Provides access to numerical tolerances, system limits, version information,
    and global settings like verbosity level. This class follows the singleton
    pattern to ensure consistent configuration across the solver.
    
    Attributes:
        Zero: Numerical zero threshold
        CoarseTol: Coarse numerical tolerance
        FineTol: Fine numerical tolerance  
        Immediate: Value representing immediate events
        MaxInt: Maximum integer value
        Version: LINE solver version string
        DummyMode: Dummy mode flag
    """

    def __repr__(self):
        return f"GlobalConstants(Version={self.Version}, Verbose={self.getVerbose()})"

    Zero = jpype.JPackage('jline').GlobalConstants.Zero
    CoarseTol = jpype.JPackage('jline').GlobalConstants.CoarseTol
    FineTol = jpype.JPackage('jline').GlobalConstants.FineTol
    Immediate = jpype.JPackage('jline').GlobalConstants.Immediate
    MaxInt = jpype.JPackage('jline').GlobalConstants.MaxInt
    Version = jpype.JPackage('jline').GlobalConstants.Version
    DummyMode = jpype.JPackage('jline').GlobalConstants.DummyMode

    _instance = None

    @classmethod
    def getInstance(cls):
        """
        Get the singleton instance of GlobalConstants.
        
        Returns:
            GlobalConstants: The singleton instance.
        """
        if cls._instance is None:
            cls._instance = cls()
        return cls._instance

    get_instance = getInstance

    @staticmethod
    def getVerbose():
        """
        Get the current verbosity level.
        
        Returns:
            VerboseLevel: Current verbosity setting (SILENT, STD, or DEBUG).
        """
        java_gc = jpype.JPackage('jline').GlobalConstants.getInstance()
        java_verbose = java_gc.getVerbose()

        if java_verbose == jpype.JPackage('jline').VerboseLevel.STD:
            return VerboseLevel.STD
        elif java_verbose == jpype.JPackage('jline').VerboseLevel.DEBUG:
            return VerboseLevel.DEBUG
        elif java_verbose == jpype.JPackage('jline').VerboseLevel.SILENT:
            return VerboseLevel.SILENT
        else:
            return VerboseLevel.STD

    get_verbose = getVerbose

    @staticmethod
    def setVerbose(verbosity):
        """
        Set the verbosity level for solver output.
        
        Args:
            verbosity (VerboseLevel): Desired verbosity level.
            
        Raises:
            ValueError: If verbosity level is invalid.
        """
        java_gc = jpype.JPackage('jline').GlobalConstants.getInstance()

        if verbosity == VerboseLevel.STD:
            java_gc.setVerbose(jpype.JPackage('jline').VerboseLevel.STD)
        elif verbosity == VerboseLevel.DEBUG:
            java_gc.setVerbose(jpype.JPackage('jline').VerboseLevel.DEBUG)
        elif verbosity == VerboseLevel.SILENT:
            java_gc.setVerbose(jpype.JPackage('jline').VerboseLevel.SILENT)
        else:
            raise ValueError(f"Invalid verbosity level: {verbosity}. Must be one of VerboseLevel.SILENT, VerboseLevel.STD, or VerboseLevel.DEBUG")

    set_verbose = setVerbose

    @classmethod
    def getConstants(cls):
        """
        Get a dictionary of all global constants.
        
        Returns:
            dict: Dictionary containing all constant values and current settings.
        """
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
