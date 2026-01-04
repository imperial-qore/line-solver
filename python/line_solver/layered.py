"""
Native Python implementation of Layered Queueing Network (LQN) models.

This module provides pure Python classes for defining and analyzing
layered queueing networks.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any, Union
from enum import Enum

from .constants import SchedStrategy


class CallType(Enum):
    """Types of calls between tasks."""
    SYNC = 'SYNC'       # Synchronous (blocking) call
    ASYNC = 'ASYNC'     # Asynchronous (non-blocking) call
    FWD = 'FWD'         # Forwarding call


class PrecedenceType(Enum):
    """Types of activity precedence patterns."""
    SERIAL = 'SERIAL'   # Sequential execution
    PARALLEL = 'PARALLEL'  # Parallel execution (AND-fork/join)
    CHOICE = 'CHOICE'      # Probabilistic choice (OR-fork/join)
    LOOP = 'LOOP'          # Repeated execution
    CACHE_ACCESS = 'CACHE_ACCESS'  # Cache access pattern (hit/miss)


class ReplacementStrategy(Enum):
    """Cache replacement strategies."""
    FIFO = 'FIFO'  # First-In First-Out
    LRU = 'LRU'    # Least Recently Used
    RR = 'RR'      # Random Replacement
    SFIFO = 'SFIFO'  # Segmented FIFO


def _get_dist_mean(dist) -> float:
    """Get mean from distribution (handles dataclass and native distributions)."""
    if dist is None:
        return 0.0
    if hasattr(dist, 'mean') and not callable(dist.mean):
        # Dataclass Distribution with .mean field
        return dist.mean
    elif hasattr(dist, 'getMean'):
        # Native distribution with getMean() method
        return dist.getMean()
    elif hasattr(dist, 'get_mean'):
        return dist.get_mean()
    else:
        return 0.0


def _get_dist_scv(dist) -> float:
    """Get SCV from distribution (handles dataclass and native distributions)."""
    if dist is None:
        return 1.0
    if hasattr(dist, 'scv') and not callable(dist.scv):
        # Dataclass Distribution with .scv field
        return dist.scv
    elif hasattr(dist, 'getSCV'):
        # Native distribution with getSCV() method
        return dist.getSCV()
    elif hasattr(dist, 'get_scv'):
        return dist.get_scv()
    else:
        return 1.0


@dataclass
class ActivityPrecedence:
    """
    Activity precedence constraint for layered queueing networks.

    Precedence constraints define the execution order of activities within a task,
    supporting serial, parallel, and conditional execution patterns.
    """
    prec_type: PrecedenceType
    activities: List['Activity'] = field(default_factory=list)
    pre_activities: List['Activity'] = field(default_factory=list)
    post_activities: List['Activity'] = field(default_factory=list)
    probabilities: List[float] = field(default_factory=list)
    count: float = 1.0  # For loops

    @staticmethod
    def Serial(*args) -> 'ActivityPrecedence':
        """
        Create a serial (sequential) precedence for a list of activities.

        Activities execute one after another in the given order.

        Supports two calling conventions:
            - Serial(a1, a2, a3): Multiple activity arguments (MATLAB-style)
            - Serial([a1, a2, a3]): Single list argument

        Args:
            *args: Either multiple Activity objects or a single list of activities

        Returns:
            ActivityPrecedence object representing serial composition

        Example:
            >>> task.add_precedence(ActivityPrecedence.Serial(a1, a2, a3))
            >>> task.add_precedence(ActivityPrecedence.Serial([a1, a2, a3]))
        """
        # Handle both forms: Serial(a1, a2, ...) and Serial([a1, a2, ...])
        if len(args) == 1 and isinstance(args[0], (list, tuple)):
            activities = list(args[0])
        else:
            activities = list(args)
        return ActivityPrecedence(
            prec_type=PrecedenceType.SERIAL,
            activities=activities
        )

    # Python snake_case aliases for compatibility
    @staticmethod
    def serial(*args) -> 'ActivityPrecedence':
        """Python snake_case alias for Serial()."""
        return ActivityPrecedence.Serial(*args)

    @staticmethod
    def AndFork(pre_act: 'Activity', post_acts: List['Activity']) -> 'ActivityPrecedence':
        """
        Create an AND-fork precedence (parallel split).

        All post-activities start executing when pre_act completes.
        Used together with AndJoin() to model parallel execution.

        Args:
            pre_act: Activity that triggers the fork
            post_acts: List of Activity objects to execute in parallel

        Returns:
            ActivityPrecedence object

        Example:
            >>> task.add_precedence(ActivityPrecedence.AndFork(start, [branch1, branch2]))
        """
        return ActivityPrecedence(
            prec_type=PrecedenceType.PARALLEL,
            pre_activities=[pre_act],
            post_activities=list(post_acts)
        )

    # Python snake_case alias
    @staticmethod
    def and_fork(pre_act: 'Activity', post_acts: List['Activity']) -> 'ActivityPrecedence':
        """Python snake_case alias for AndFork()."""
        return ActivityPrecedence.AndFork(pre_act, post_acts)

    @staticmethod
    def AndJoin(pre_acts: List['Activity'], post_act: 'Activity') -> 'ActivityPrecedence':
        """
        Create an AND-join precedence (synchronization).

        Post-activity starts when ALL pre-activities have completed.
        Used together with AndFork() to model parallel execution.

        Args:
            pre_acts: List of Activity objects to synchronize on
            post_act: Activity that executes after synchronization

        Returns:
            ActivityPrecedence object

        Example:
            >>> task.add_precedence(ActivityPrecedence.AndJoin([branch1, branch2], end))
        """
        return ActivityPrecedence(
            prec_type=PrecedenceType.PARALLEL,
            pre_activities=list(pre_acts),
            post_activities=[post_act]
        )

    # Python snake_case alias
    @staticmethod
    def and_join(pre_acts: List['Activity'], post_act: 'Activity') -> 'ActivityPrecedence':
        """Python snake_case alias for AndJoin()."""
        return ActivityPrecedence.AndJoin(pre_acts, post_act)

    @staticmethod
    def OrFork(pre_act: 'Activity', post_acts: List['Activity'],
                probs: List[float]) -> 'ActivityPrecedence':
        """
        Create an OR-fork precedence (probabilistic branching).

        Exactly one post-activity is selected based on probabilities when
        pre_act completes.

        Args:
            pre_act: Activity that triggers the fork
            post_acts: List of Activity objects as branch options
            probs: List of probabilities for each branch (must sum to 1.0)

        Returns:
            ActivityPrecedence object

        Example:
            >>> task.add_precedence(ActivityPrecedence.OrFork(start, [fast, slow], [0.7, 0.3]))
        """
        return ActivityPrecedence(
            prec_type=PrecedenceType.CHOICE,
            pre_activities=[pre_act],
            post_activities=list(post_acts),
            probabilities=list(probs)
        )

    # Python snake_case alias
    @staticmethod
    def or_fork(pre_act: 'Activity', post_acts: List['Activity'],
                probs: List[float]) -> 'ActivityPrecedence':
        """Python snake_case alias for OrFork()."""
        return ActivityPrecedence.OrFork(pre_act, post_acts, probs)

    @staticmethod
    def OrJoin(pre_acts: List['Activity'], post_act: 'Activity') -> 'ActivityPrecedence':
        """
        Create an OR-join precedence (merge).

        Post-activity starts when ANY of the pre-activities complete.
        Used together with OrFork() to model probabilistic branching.

        Args:
            pre_acts: List of Activity objects to merge
            post_act: Activity that executes after merge

        Returns:
            ActivityPrecedence object

        Example:
            >>> task.add_precedence(ActivityPrecedence.OrJoin([fast, slow], end))
        """
        return ActivityPrecedence(
            prec_type=PrecedenceType.CHOICE,
            pre_activities=list(pre_acts),
            post_activities=[post_act]
        )

    # Python snake_case alias
    @staticmethod
    def or_join(pre_acts: List['Activity'], post_act: 'Activity') -> 'ActivityPrecedence':
        """Python snake_case alias for OrJoin()."""
        return ActivityPrecedence.OrJoin(pre_acts, post_act)

    @staticmethod
    def Loop(pre_act: 'Activity', loop_acts: List['Activity'],
             count: float) -> 'ActivityPrecedence':
        """
        Create a loop precedence for repeated execution.

        Loop activities execute a specified number of times before continuing.

        Args:
            pre_act: Activity that triggers the loop
            loop_acts: List of Activity objects in the loop body
            count: Number of loop iterations (can be fractional for geometric mean)

        Returns:
            ActivityPrecedence object

        Example:
            >>> task.add_precedence(ActivityPrecedence.Loop(init, [compute], 5))
        """
        return ActivityPrecedence(
            prec_type=PrecedenceType.LOOP,
            pre_activities=[pre_act],
            activities=list(loop_acts),
            count=float(count)
        )

    # Python snake_case alias
    @staticmethod
    def loop(pre_act: 'Activity', loop_acts: List['Activity'],
             count: float) -> 'ActivityPrecedence':
        """Python snake_case alias for Loop()."""
        return ActivityPrecedence.Loop(pre_act, loop_acts, count)

    @staticmethod
    def CacheAccess(access_act: 'Activity',
                    outcome_acts: List['Activity']) -> 'ActivityPrecedence':
        """
        Create a cache access precedence pattern.

        Models cache hit/miss behavior where access_act performs the cache lookup
        and outcome_acts contains [hit_activity, miss_activity].

        Args:
            access_act: Activity that performs cache access
            outcome_acts: List of [hit_activity, miss_activity]

        Returns:
            ActivityPrecedence object

        Example:
            >>> task.add_precedence(ActivityPrecedence.CacheAccess(lookup, [hit, miss]))
        """
        return ActivityPrecedence(
            prec_type=PrecedenceType.CACHE_ACCESS,
            pre_activities=[access_act],
            post_activities=list(outcome_acts)
        )

    # Python snake_case alias
    @staticmethod
    def cache_access(access_act: 'Activity',
                     outcome_acts: List['Activity']) -> 'ActivityPrecedence':
        """Python snake_case alias for CacheAccess()."""
        return ActivityPrecedence.CacheAccess(access_act, outcome_acts)


@dataclass
class Distribution:
    """Simple distribution representation for service times."""
    mean: float
    scv: float = 1.0  # Squared coefficient of variation (1.0 = exponential)

    @classmethod
    def exponential(cls, mean: float) -> 'Distribution':
        """Create exponential distribution with given mean."""
        return cls(mean=mean, scv=1.0)

    @classmethod
    def deterministic(cls, value: float) -> 'Distribution':
        """Create deterministic (constant) distribution."""
        return cls(mean=value, scv=0.0)


class Activity:
    """
    Activity in a layered queueing network.

    An activity represents a unit of work performed by a task.
    Activities have service time distributions and can make calls
    to other entries.

    Supports two calling conventions:
    1. Activity(name, host_demand)
    2. Activity(model, name, host_demand)
    """

    def __init__(self, model_or_name, name_or_demand=None, demand=None):
        """Initialize an Activity with flexible arguments."""
        if demand is not None:
            # 3-arg call: (model, name, host_demand)
            self._model = model_or_name
            self.name = name_or_demand
            self.host_demand = self._convert_demand(demand)
            # Register with model
            if hasattr(self._model, 'add_activity'):
                self._model.add_activity(self)
        else:
            # 2-arg call: (name, host_demand)
            self._model = None
            self.name = model_or_name
            self.host_demand = self._convert_demand(name_or_demand)

        self.task = None
        self.bound_entry = None
        self.reply_entry = None
        self.calls = []  # List of (entry, mean_calls, call_type)

    def _convert_demand(self, demand):
        """Convert external distribution to internal Distribution if needed."""
        if demand is None:
            return None
        if isinstance(demand, Distribution):
            return demand
        # If it's an Exp or other external distribution, wrap it
        if hasattr(demand, 'getMean'):
            mean = demand.getMean()
            return Distribution.exponential(mean)
        return demand

    def __hash__(self):
        return id(self)

    def __eq__(self, other):
        return self is other

    @property
    def obj(self):
        """Return self for compatibility with wrapper code that accesses .obj"""
        return self

    def on(self, task: 'Task') -> 'Activity':
        """Assign this activity to a task."""
        self.task = task
        task.activities.append(self)
        return self

    def bound_to(self, entry: 'Entry') -> 'Activity':
        """Bind this activity to an entry (first activity of entry)."""
        self.bound_entry = entry
        entry.bound_activity = self
        return self

    def synch_call(self, entry: 'Entry', mean_calls: float = 1.0) -> 'Activity':
        """Add a synchronous call to another entry."""
        self.calls.append((entry, mean_calls, CallType.SYNC))
        return self

    def asynch_call(self, entry: 'Entry', mean_calls: float = 1.0) -> 'Activity':
        """Add an asynchronous call to another entry."""
        self.calls.append((entry, mean_calls, CallType.ASYNC))
        return self

    def replies_to(self, entry: 'Entry') -> 'Activity':
        """Mark this activity as replying to an entry."""
        self.reply_entry = entry
        return self

    def getHostDemand(self):
        """Get host demand distribution."""
        return self.host_demand

    def getHostDemandMean(self) -> float:
        """Get host demand mean."""
        if self.host_demand is None:
            return 0.0
        return getattr(self.host_demand, 'mean', 0.0)

    def getHostDemandSCV(self) -> float:
        """Get host demand SCV."""
        if self.host_demand is None:
            return 1.0
        return getattr(self.host_demand, 'scv', 1.0)

    def getCallOrder(self) -> list:
        """Get call order."""
        return self.calls

    def getBoundToEntry(self):
        """Get bound entry."""
        return self.bound_entry

    def getParent(self):
        """Get parent task."""
        return self.task

    def getSyncCallDests(self) -> list:
        """Get synchronous call destinations."""
        return [entry for entry, _, call_type in self.calls if call_type == CallType.SYNC]

    def getSyncCallMeans(self) -> list:
        """Get synchronous call mean counts."""
        return [mean for _, mean, call_type in self.calls if call_type == CallType.SYNC]

    def getAsyncCallDests(self) -> list:
        """Get asynchronous call destinations."""
        return [entry for entry, _, call_type in self.calls if call_type == CallType.ASYNC]

    def getAsyncCallMeans(self) -> list:
        """Get asynchronous call mean counts."""
        return [mean for _, mean, call_type in self.calls if call_type == CallType.ASYNC]

    def getThinkTimeMean(self) -> float:
        """Get think time mean."""
        return 0.0  # Activities don't have think time

    # Snake_case aliases
    get_host_demand = getHostDemand
    get_host_demand_mean = getHostDemandMean
    get_host_demand_scv = getHostDemandSCV
    get_call_order = getCallOrder
    get_bound_to_entry = getBoundToEntry
    get_parent = getParent
    get_sync_call_dests = getSyncCallDests
    get_sync_call_means = getSyncCallMeans
    get_async_call_dests = getAsyncCallDests
    get_async_call_means = getAsyncCallMeans
    get_think_time_mean = getThinkTimeMean


class Entry:
    """
    Entry in a layered queueing network.

    An entry is a service interface provided by a task. Entries
    are called by other tasks and define the work performed
    through their bound activities.

    Supports two calling conventions:
    1. Entry(name)
    2. Entry(model, name)
    """

    def __init__(self, model_or_name, name=None):
        """Initialize an Entry with flexible arguments."""
        if name is not None:
            # 2-arg call: (model, name)
            self._model = model_or_name
            self.name = name
            # Register with model
            if hasattr(self._model, 'add_entry'):
                self._model.add_entry(self)
        else:
            # 1-arg call: (name,)
            self._model = None
            self.name = model_or_name

        self.task = None
        self.bound_activity = None
        self._forwarding_dests = []  # List of target entry names
        self._forwarding_probs = []  # List of forwarding probabilities
        self._arrival = None  # Arrival distribution for open arrivals

    def __hash__(self):
        return id(self)

    def __eq__(self, other):
        return self is other

    @property
    def obj(self):
        """Return self for compatibility with wrapper code that accesses .obj"""
        return self

    def on(self, task: 'Task') -> 'Entry':
        """Assign this entry to a task."""
        self.task = task
        task.entries.append(self)
        return self

    def getBoundToActivity(self):
        """Get the bound activity."""
        return self.bound_activity

    def getReplyActivity(self):
        """Get the reply activity (same as bound activity for simple entries)."""
        return self.bound_activity

    def getParent(self):
        """Get the parent task."""
        return self.task

    def getForwardingDests(self) -> list:
        """Get forwarding destinations."""
        return self._forwarding_dests

    def getForwardingProbs(self) -> list:
        """Get forwarding probabilities."""
        return self._forwarding_probs

    def getArrival(self):
        """Get arrival distribution for open arrivals."""
        return self._arrival

    def setArrival(self, arrival) -> 'Entry':
        """Set arrival distribution for open arrivals."""
        self._arrival = arrival
        return self

    def addForwarding(self, target_entry: 'Entry', prob: float = 1.0) -> 'Entry':
        """Add forwarding to another entry.

        Args:
            target_entry: Entry to forward requests to
            prob: Forwarding probability (default 1.0)

        Returns:
            self for method chaining
        """
        self._forwarding_dests.append(target_entry.name if hasattr(target_entry, 'name') else target_entry)
        self._forwarding_probs.append(prob)
        return self

    # Snake_case aliases
    get_bound_to_activity = getBoundToActivity
    get_reply_activity = getReplyActivity
    get_parent = getParent
    get_forwarding_dests = getForwardingDests
    get_forwarding_probs = getForwardingProbs
    get_arrival = getArrival
    set_arrival = setArrival
    add_forwarding = addForwarding


class Task:
    """
    Task in a layered queueing network.

    A task represents a software process or thread that provides
    services through entries. Tasks are deployed on processors
    and have a multiplicity (number of instances/threads).

    Supports two calling conventions:
    1. Task(name, multiplicity, sched_strategy)
    2. Task(model, name, multiplicity, sched_strategy)
    """

    def __init__(self, model_or_name, name_or_mult=None, mult_or_sched=None, sched=None):
        """Initialize a Task with flexible arguments."""
        # Detect calling convention
        if sched is not None:
            # 4-arg call: (model, name, multiplicity, sched_strategy)
            self._model = model_or_name
            self.name = name_or_mult
            self.multiplicity = mult_or_sched
            self.sched_strategy = self._convert_sched_strategy(sched)
            # Register with model
            if hasattr(self._model, 'add_task'):
                self._model.add_task(self)
        else:
            # 3-arg call: (name, multiplicity, sched_strategy)
            self._model = None
            self.name = model_or_name
            self.multiplicity = name_or_mult
            self.sched_strategy = self._convert_sched_strategy(mult_or_sched)

        self.processor = None
        self.think_time = None
        self.setup_time = None
        self.delay_off_time = None
        self.entries = []
        self.activities = []
        self.precedences = []

    def _convert_sched_strategy(self, sched):
        """Convert SchedStrategy to SchedStrategy if needed."""
        if hasattr(sched, 'name') and hasattr(SchedStrategy, sched.name):
            return getattr(SchedStrategy, sched.name)
        return sched

    def __hash__(self):
        return id(self)

    def __eq__(self, other):
        return self is other

    @property
    def obj(self):
        """Return self for compatibility with wrapper code that accesses .obj"""
        return self

    def on(self, processor: 'Processor') -> 'Task':
        """Deploy this task on a processor."""
        self.processor = processor
        processor.tasks.append(self)
        return self

    def set_think_time(self, think_time: Union[Distribution, float]) -> 'Task':
        """Set the think time for this task."""
        if isinstance(think_time, (int, float)):
            self.think_time = Distribution.exponential(float(think_time))
        else:
            self.think_time = think_time
        return self

    def setSetupTime(self, setup_time) -> 'Task':
        """Set the setup time (cold start delay) for this task."""
        self.setup_time = setup_time
        return self

    def set_setup_time(self, setup_time) -> 'Task':
        """Set the setup time (cold start delay) for this task."""
        self.setup_time = setup_time
        return self

    def setDelayOffTime(self, delay_off_time) -> 'Task':
        """Set the delay-off time (teardown delay) for this task."""
        self.delay_off_time = delay_off_time
        return self

    def set_delay_off_time(self, delay_off_time) -> 'Task':
        """Set the delay-off time (teardown delay) for this task."""
        self.delay_off_time = delay_off_time
        return self

    def add_precedence(self, precedence: Union['ActivityPrecedence', List['ActivityPrecedence']]) -> 'Task':
        """
        Add activity precedence constraint(s) to this task.

        Precedence constraints define the execution order of activities within
        the task, supporting serial, parallel, and conditional patterns.

        Args:
            precedence: An ActivityPrecedence object or list of ActivityPrecedence objects,
                       created using serial(), and_fork(), and_join(), or_fork(), or_join(), or loop()

        Returns:
            self (for method chaining)

        Example:
            >>> task.add_precedence(ActivityPrecedence.serial([a1, a2, a3]))
            >>> task.add_precedence([
            ...     ActivityPrecedence.serial([a1, a2]),
            ...     ActivityPrecedence.and_fork(a1, [a3, a4])
            ... ])
        """
        if isinstance(precedence, (list, tuple)):
            self.precedences.extend(precedence)
        else:
            self.precedences.append(precedence)
        return self

    def getMultiplicity(self) -> int:
        """Get multiplicity."""
        return self.multiplicity

    def getReplication(self) -> int:
        """Get replication level."""
        return getattr(self, '_replication', 1)

    def setReplication(self, replication: int) -> 'Task':
        """Set replication level."""
        self._replication = replication
        return self

    def getScheduling(self):
        """Get scheduling strategy."""
        return self.sched_strategy

    def getThinkTimeMean(self) -> float:
        """Get think time mean."""
        if self.think_time is None:
            return 0.0
        return getattr(self.think_time, 'mean', 0.0)

    def getThinkTimeSCV(self) -> float:
        """Get think time SCV."""
        if self.think_time is None:
            return 1.0
        return getattr(self.think_time, 'scv', 1.0)

    def getParent(self):
        """Get parent processor."""
        return self.processor

    def getPrecedences(self) -> list:
        """Get activity precedences."""
        return self.precedences

    def getSetupTimeMean(self) -> float:
        """Get setup time mean."""
        if self.setup_time is None:
            return 0.0
        return getattr(self.setup_time, 'mean', 0.0)

    def getDelayOffTimeMean(self) -> float:
        """Get delay-off time mean."""
        if self.delay_off_time is None:
            return 0.0
        return getattr(self.delay_off_time, 'mean', 0.0)

    def setThinkTime(self, think_time) -> 'Task':
        """Set think time."""
        return self.set_think_time(think_time)

    # Snake_case aliases
    get_multiplicity = getMultiplicity
    get_replication = getReplication
    set_replication = setReplication
    get_scheduling = getScheduling
    get_think_time_mean = getThinkTimeMean
    get_think_time_scv = getThinkTimeSCV
    get_parent = getParent
    get_precedences = getPrecedences
    get_setup_time_mean = getSetupTimeMean
    get_delay_off_time_mean = getDelayOffTimeMean


class Processor:
    """
    Processor in a layered queueing network.

    A processor represents a computing resource (CPU, server, etc.)
    that hosts tasks. Processors have a multiplicity (number of
    identical resources) and a scheduling strategy.

    Supports two calling conventions:
    1. Processor(name, multiplicity, sched_strategy)
    2. Processor(model, name, multiplicity, sched_strategy)
    """

    def __init__(self, model_or_name, name_or_mult=None, mult_or_sched=None, sched=None):
        """Initialize a Processor with flexible arguments."""
        # Detect calling convention
        if sched is not None:
            # 4-arg call: (model, name, multiplicity, sched_strategy)
            self._model = model_or_name
            self.name = name_or_mult
            self.multiplicity = mult_or_sched
            self.sched_strategy = self._convert_sched_strategy(sched)
            # Register with model
            if hasattr(self._model, 'add_processor'):
                self._model.add_processor(self)
        else:
            # 3-arg call: (name, multiplicity, sched_strategy)
            self._model = None
            self.name = model_or_name
            self.multiplicity = name_or_mult
            self.sched_strategy = self._convert_sched_strategy(mult_or_sched)

        self.tasks = []
        self._quantum = 0.0  # Time quantum for PS scheduling
        self._speed_factor = 1.0  # Processor speed factor

    def _convert_sched_strategy(self, sched):
        """Convert SchedStrategy to SchedStrategy if needed."""
        if hasattr(sched, 'name') and hasattr(SchedStrategy, sched.name):
            return getattr(SchedStrategy, sched.name)
        return sched

    def __hash__(self):
        return id(self)

    def __eq__(self, other):
        return self is other

    @property
    def obj(self):
        """Return self for compatibility with wrapper code that accesses .obj"""
        return self

    def getMultiplicity(self) -> int:
        """Get multiplicity."""
        return self.multiplicity

    def getReplication(self) -> int:
        """Get replication level."""
        return getattr(self, '_replication', 1)

    def setReplication(self, replication: int) -> 'Processor':
        """Set replication level."""
        self._replication = replication
        return self

    def getScheduling(self):
        """Get scheduling strategy."""
        return self.sched_strategy

    def getQuantum(self) -> float:
        """Get quantum for PS scheduling."""
        return self._quantum

    def setQuantum(self, quantum: float) -> 'Processor':
        """Set quantum for PS scheduling."""
        self._quantum = quantum
        return self

    def getSpeedFactor(self) -> float:
        """Get speed factor."""
        return self._speed_factor

    def setSpeedFactor(self, speed_factor: float) -> 'Processor':
        """Set speed factor."""
        self._speed_factor = speed_factor
        return self

    # Snake_case aliases
    get_multiplicity = getMultiplicity
    get_replication = getReplication
    set_replication = setReplication
    get_scheduling = getScheduling
    get_quantum = getQuantum
    set_quantum = setQuantum
    get_speed_factor = getSpeedFactor
    set_speed_factor = setSpeedFactor


class CacheTask(Task):
    """
    Cache task in a layered queueing network.

    A CacheTask models a caching service that stores items in a limited
    capacity cache. It tracks cache hits and misses based on a replacement
    strategy.
    """

    def __init__(self, model, name: str, total_items: int, cache_capacity: int,
                 replacement_strategy: ReplacementStrategy, multiplicity: float = 1):
        """
        Create a cache task.

        Args:
            model: Parent LayeredNetwork
            name: Name of the cache task
            total_items: Total number of distinct items that can be requested
            cache_capacity: Maximum number of items the cache can hold
            replacement_strategy: Cache replacement policy (FIFO, LRU, RR, etc.)
            multiplicity: Number of task instances
        """
        # Initialize base Task fields
        self._model = model
        self.name = name
        self.multiplicity = multiplicity
        self.sched_strategy = SchedStrategy.FCFS
        self.processor = None
        self.think_time = None
        self.entries = []
        self.activities = []
        self.precedences = []
        # CacheTask specific fields
        self.total_items = total_items
        self.cache_capacity = cache_capacity
        self.replacement_strategy = replacement_strategy
        # Register with model
        if model is not None and hasattr(model, 'tasks'):
            model.tasks.append(self)

    def __hash__(self):
        return id(self)

    def __eq__(self, other):
        return self is other

    def on(self, processor: 'Processor') -> 'CacheTask':
        """Deploy this cache task on a processor."""
        self.processor = processor
        processor.tasks.append(self)
        return self


class ItemEntry(Entry):
    """
    Item entry for a cache task.

    An ItemEntry represents the interface to request items from a cache.
    It specifies the total number of items and their access probabilities.
    """

    def __init__(self, model, name: str, total_items: int, access_prob):
        """
        Create an item entry.

        Args:
            model: Parent LayeredNetwork
            name: Name of the entry
            total_items: Total number of distinct items
            access_prob: Access probability distribution (DiscreteSampler or list)
        """
        self._model = model
        self.name = name
        self.task = None
        self.bound_activity = None
        self.total_items = total_items
        self.access_prob = access_prob
        # Register with model
        if model is not None and hasattr(model, 'entries'):
            model.entries.append(self)

    def __hash__(self):
        return id(self)

    def __eq__(self, other):
        return self is other

    def on(self, task: 'CacheTask') -> 'ItemEntry':
        """Assign this item entry to a cache task."""
        self.task = task
        task.entries.append(self)
        return self


# Convenience aliases
CacheTask = CacheTask
ItemEntry = ItemEntry


@dataclass
class LayeredNetworkStruct:
    """
    Internal structure representation of a layered queueing network.

    This structure contains all the numerical data needed for analysis,
    extracted from the high-level model objects.
    """
    # Counts
    nhosts: int = 0
    ntasks: int = 0
    nentries: int = 0
    nacts: int = 0
    ncalls: int = 0
    nidx: int = 0

    # Index shifts (for mapping indices)
    hshift: int = 0
    tshift: int = 0
    eshift: int = 0
    ashift: int = 0

    # Matrices
    mult: np.ndarray = None          # Multiplicities
    graph: np.ndarray = None         # Adjacency graph
    iscaller: np.ndarray = None      # Caller matrix
    issynccaller: np.ndarray = None  # Sync caller matrix
    isasynccaller: np.ndarray = None # Async caller matrix
    isref: np.ndarray = None         # Reference task flags
    callpair: np.ndarray = None      # Call pairs (caller_act, callee_entry, mean_calls)
    parent: np.ndarray = None        # Parent relationships
    replygraph: np.ndarray = None    # Reply graph (nacts x nentries)
    actphase: np.ndarray = None      # Activity phase (1 or 2) for each activity

    # Names
    names: np.ndarray = None
    hashnames: np.ndarray = None

    # Mappings
    tasksof: Dict[int, List[int]] = field(default_factory=dict)
    entriesof: Dict[int, List[int]] = field(default_factory=dict)
    actsof: Dict[int, List[int]] = field(default_factory=dict)
    callsof: Dict[int, List[int]] = field(default_factory=dict)

    # Service demands and think times
    hostdem: Dict[int, float] = field(default_factory=dict)
    think: Dict[int, float] = field(default_factory=dict)

    # Scheduling strategies
    sched: Dict[int, str] = field(default_factory=dict)


class LayeredNetwork:
    """
    Native Python implementation of a Layered Queueing Network.

    This class provides a pure Python way to define and analyze
    layered queueing networks.

    Example:
        >>> model = LayeredNetwork('ClientServer')
        >>> P1 = model.add_processor('ClientProc', 1, SchedStrategy.PS)
        >>> P2 = model.add_processor('ServerProc', 1, SchedStrategy.PS)
        >>> T1 = model.add_task('Client', 5, SchedStrategy.REF, P1)
        >>> T1.set_think_time(2.0)
        >>> T2 = model.add_task('Server', float('inf'), SchedStrategy.INF, P2)
        >>> E1 = model.add_entry('ClientEntry', T1)
        >>> E2 = model.add_entry('ServerEntry', T2)
        >>> A1 = model.add_activity('ClientAct', 0.5, T1)
        >>> A1.bound_to(E1).synch_call(E2, 1.0)
        >>> A2 = model.add_activity('ServerAct', 1.0, T2)
        >>> A2.bound_to(E2).replies_to(E2)
    """

    def __init__(self, name: str):
        """Initialize a new layered queueing network."""
        self.name = name
        self.processors: List[Processor] = []
        self.tasks: List[Task] = []
        self.entries: List[Entry] = []
        self.activities: List[Activity] = []

        # Index mappings (built when getStruct is called)
        self._proc_idx: Dict[Processor, int] = {}
        self._task_idx: Dict[Task, int] = {}
        self._entry_idx: Dict[Entry, int] = {}
        self._act_idx: Dict[Activity, int] = {}

    @property
    def _wrapper_nodes(self) -> set:
        """Return all nodes in the network for compatibility with wrapper code."""
        return set(self.processors + self.tasks + self.entries + self.activities)

    def add_processor(self, name_or_proc, multiplicity: float = None,
                     sched_strategy: SchedStrategy = None) -> Processor:
        """
        Add a processor to the network.

        Supports two calling conventions:
        1. add_processor(Processor_instance)
        2. add_processor(name, multiplicity, sched_strategy)
        """
        if isinstance(name_or_proc, Processor):
            proc = name_or_proc
        else:
            proc = Processor(name_or_proc, multiplicity, sched_strategy)
        if proc not in self.processors:
            self.processors.append(proc)
        return proc

    def add_task(self, name_or_task, multiplicity: float = None,
                sched_strategy: SchedStrategy = None,
                processor: Processor = None) -> Task:
        """
        Add a task to the network.

        Supports two calling conventions:
        1. add_task(Task_instance)
        2. add_task(name, multiplicity, sched_strategy, processor)
        """
        if isinstance(name_or_task, Task):
            task = name_or_task
        else:
            task = Task(name_or_task, multiplicity, sched_strategy)
            if processor is not None:
                task.on(processor)
        if task not in self.tasks:
            self.tasks.append(task)
        return task

    def add_entry(self, name_or_entry, task: Task = None) -> Entry:
        """
        Add an entry to the network.

        Supports two calling conventions:
        1. add_entry(Entry_instance)
        2. add_entry(name, task)
        """
        if isinstance(name_or_entry, Entry):
            entry = name_or_entry
        else:
            entry = Entry(name_or_entry)
            if task is not None:
                entry.on(task)
        if entry not in self.entries:
            self.entries.append(entry)
        return entry

    def add_activity(self, name_or_activity, host_demand: Union[Distribution, float] = None,
                    task: Task = None) -> Activity:
        """
        Add an activity to the network.

        Supports two calling conventions:
        1. add_activity(Activity_instance)
        2. add_activity(name, host_demand, task)
        """
        if isinstance(name_or_activity, Activity):
            activity = name_or_activity
        else:
            if isinstance(host_demand, (int, float)):
                demand = Distribution.exponential(float(host_demand))
            else:
                demand = host_demand
            activity = Activity(name_or_activity, demand)
            if task is not None:
                activity.on(task)
        if activity not in self.activities:
            self.activities.append(activity)
        return activity

    def _build_indices(self):
        """Build index mappings for all elements."""
        # Index layout: [hosts | tasks | entries | activities]
        idx = 1  # 1-based indexing

        # Processors (hosts)
        for proc in self.processors:
            self._proc_idx[proc] = idx
            idx += 1

        tshift = idx - 1

        # Tasks
        for task in self.tasks:
            self._task_idx[task] = idx
            idx += 1

        eshift = idx - 1

        # Entries
        for entry in self.entries:
            self._entry_idx[entry] = idx
            idx += 1

        ashift = idx - 1

        # Activities
        for act in self.activities:
            self._act_idx[act] = idx
            idx += 1

        return tshift, eshift, ashift, idx - 1

    def getStruct(self) -> LayeredNetworkStruct:
        """
        Build and return the internal structure representation.

        This method converts the high-level model objects into
        numerical arrays and matrices suitable for analysis.
        """
        tshift, eshift, ashift, nidx = self._build_indices()

        lqn = LayeredNetworkStruct()
        lqn.nhosts = len(self.processors)
        lqn.ntasks = len(self.tasks)
        lqn.nentries = len(self.entries)
        lqn.nacts = len(self.activities)
        lqn.nidx = nidx

        lqn.hshift = 0
        lqn.tshift = tshift
        lqn.eshift = eshift
        lqn.ashift = ashift

        # Build names arrays
        lqn.names = np.empty(nidx + 1, dtype=object)
        lqn.hashnames = np.empty(nidx + 1, dtype=object)
        lqn.names[0] = ''
        lqn.hashnames[0] = ''

        for proc in self.processors:
            idx = self._proc_idx[proc]
            lqn.names[idx] = proc.name
            lqn.hashnames[idx] = proc.name

        for task in self.tasks:
            idx = self._task_idx[task]
            lqn.names[idx] = task.name
            lqn.hashnames[idx] = task.name

        for entry in self.entries:
            idx = self._entry_idx[entry]
            lqn.names[idx] = entry.name
            lqn.hashnames[idx] = entry.name

        for act in self.activities:
            idx = self._act_idx[act]
            lqn.names[idx] = act.name
            lqn.hashnames[idx] = act.name

        # Build multiplicity matrix
        lqn.mult = np.zeros((1, nidx + 1))
        for proc in self.processors:
            lqn.mult[0, self._proc_idx[proc]] = proc.multiplicity
        for task in self.tasks:
            mult = task.multiplicity
            if mult == float('inf'):
                mult = 1e9  # Large number for infinite server
            lqn.mult[0, self._task_idx[task]] = mult

        # Build reference task flags
        lqn.isref = np.zeros((nidx + 1, 1))
        for task in self.tasks:
            if task.sched_strategy == SchedStrategy.REF:
                lqn.isref[self._task_idx[task], 0] = 1

        # Build caller matrices
        lqn.iscaller = np.zeros((nidx + 1, nidx + 1))
        lqn.issynccaller = np.zeros((nidx + 1, nidx + 1))
        lqn.isasynccaller = np.zeros((nidx + 1, nidx + 1))

        # Build call pairs and caller relationships
        calls = []
        for act in self.activities:
            # Skip activities without assigned tasks
            if act.task is None:
                continue

            act_idx = self._act_idx[act]
            caller_task = act.task
            caller_task_idx = self._task_idx[caller_task]

            for target_entry, mean_calls, call_type in act.calls:
                # Skip calls to entries without tasks
                if target_entry.task is None:
                    continue

                target_entry_idx = self._entry_idx[target_entry]
                target_task = target_entry.task
                target_task_idx = self._task_idx[target_task]

                # Mark caller relationship
                lqn.iscaller[caller_task_idx, target_task_idx] = 1

                if call_type == CallType.SYNC:
                    lqn.issynccaller[caller_task_idx, target_task_idx] = 1
                else:
                    lqn.isasynccaller[caller_task_idx, target_task_idx] = 1

                calls.append((act_idx, target_entry_idx, mean_calls))

        lqn.ncalls = len(calls)
        if calls:
            lqn.callpair = np.zeros((lqn.ncalls + 1, 4))
            for i, (act_idx, entry_idx, mean_calls) in enumerate(calls, 1):
                lqn.callpair[i, 1] = act_idx
                lqn.callpair[i, 2] = entry_idx
                lqn.callpair[i, 3] = mean_calls
        else:
            lqn.callpair = np.zeros((1, 4))

        # Build parent relationships
        lqn.parent = np.zeros((nidx + 1, 1))
        for task in self.tasks:
            task_idx = self._task_idx[task]
            if task.processor:
                lqn.parent[task_idx, 0] = self._proc_idx[task.processor]
        for entry in self.entries:
            entry_idx = self._entry_idx[entry]
            if entry.task:
                lqn.parent[entry_idx, 0] = self._task_idx[entry.task]
        for act in self.activities:
            act_idx = self._act_idx[act]
            if act.task:
                lqn.parent[act_idx, 0] = self._task_idx[act.task]

        # Build tasksof mapping (tasks on each host)
        lqn.tasksof = {}
        for hidx in range(1, lqn.nhosts + 1):
            lqn.tasksof[hidx] = []
        for task in self.tasks:
            task_idx = self._task_idx[task]
            if task.processor:
                proc_idx = self._proc_idx[task.processor]
                lqn.tasksof[proc_idx].append(task_idx)

        # Build entriesof mapping (entries of each task)
        lqn.entriesof = {}
        for t in range(1, lqn.ntasks + 1):
            tidx = lqn.tshift + t
            lqn.entriesof[tidx] = []
        for entry in self.entries:
            entry_idx = self._entry_idx[entry]
            if entry.task:
                task_idx = self._task_idx[entry.task]
                if task_idx not in lqn.entriesof:
                    lqn.entriesof[task_idx] = []
                lqn.entriesof[task_idx].append(entry_idx)

        # Build actsof mapping (activities of each entry/task)
        lqn.actsof = {}
        for entry in self.entries:
            entry_idx = self._entry_idx[entry]
            lqn.actsof[entry_idx] = []
        for task in self.tasks:
            task_idx = self._task_idx[task]
            lqn.actsof[task_idx] = []

        for act in self.activities:
            act_idx = self._act_idx[act]
            if act.bound_entry:
                entry_idx = self._entry_idx[act.bound_entry]
                if entry_idx not in lqn.actsof:
                    lqn.actsof[entry_idx] = []
                lqn.actsof[entry_idx].append(act_idx)
            if act.task:
                task_idx = self._task_idx[act.task]
                if task_idx not in lqn.actsof:
                    lqn.actsof[task_idx] = []
                if act_idx not in lqn.actsof[task_idx]:
                    lqn.actsof[task_idx].append(act_idx)

        # Build host demands
        lqn.hostdem = {}
        for act in self.activities:
            act_idx = self._act_idx[act]
            lqn.hostdem[act_idx] = _get_dist_mean(act.host_demand)

        # Build think times
        lqn.think = {}
        for task in self.tasks:
            task_idx = self._task_idx[task]
            if task.think_time:
                lqn.think[task_idx] = _get_dist_mean(task.think_time)
            else:
                lqn.think[task_idx] = 0.0

        # Build scheduling strategies
        lqn.sched = {}
        for proc in self.processors:
            lqn.sched[self._proc_idx[proc]] = proc.sched_strategy.value
        for task in self.tasks:
            lqn.sched[self._task_idx[task]] = task.sched_strategy.value

        # Build graph (adjacency matrix)
        lqn.graph = np.zeros((nidx + 1, nidx + 1))
        # Add edges from processors to tasks
        for task in self.tasks:
            task_idx = self._task_idx[task]
            if task.processor:
                proc_idx = self._proc_idx[task.processor]
                lqn.graph[proc_idx, task_idx] = 1
        # Add edges from tasks to entries
        for entry in self.entries:
            entry_idx = self._entry_idx[entry]
            if entry.task:
                task_idx = self._task_idx[entry.task]
                lqn.graph[task_idx, entry_idx] = 1
        # Add edges from entries to activities
        for act in self.activities:
            act_idx = self._act_idx[act]
            if act.bound_entry:
                entry_idx = self._entry_idx[act.bound_entry]
                lqn.graph[entry_idx, act_idx] = 1
        # Add edges for calls
        for act in self.activities:
            act_idx = self._act_idx[act]
            for target_entry, _, _ in act.calls:
                entry_idx = self._entry_idx[target_entry]
                lqn.graph[act_idx, entry_idx] = 1

        # Build replygraph (nacts x nentries) - which activities reply to which entries
        lqn.replygraph = np.zeros((lqn.nacts, lqn.nentries))
        for act in self.activities:
            act_idx = self._act_idx[act] - ashift  # Convert to 0-based activity index
            if hasattr(act, 'reply_entries') and act.reply_entries:
                for reply_entry in act.reply_entries:
                    entry_idx = self._entry_idx[reply_entry] - eshift  # Convert to 0-based entry index
                    if 0 <= act_idx < lqn.nacts and 0 <= entry_idx < lqn.nentries:
                        lqn.replygraph[act_idx, entry_idx] = 1

        # Build actphase - phase number (1 or 2) for each activity
        lqn.actphase = np.ones(lqn.nacts)  # Default phase is 1
        for act in self.activities:
            act_idx = self._act_idx[act] - ashift  # Convert to 0-based activity index
            if hasattr(act, 'phase') and act.phase is not None:
                lqn.actphase[act_idx] = act.phase

        # Validation: Check for non-terminal reply activities
        # An activity that replies to an entry should not have Phase 1 successor activities
        # Phase 2 successors are allowed (post-reply processing)
        for a in range(lqn.nacts):
            if np.any(lqn.replygraph[a, :] > 0):  # activity 'a' replies to some entry
                aidx = ashift + a  # global activity index
                for succ in range(nidx + 1):
                    if lqn.graph[aidx, succ] != 0:  # successor exists
                        if succ > eshift + lqn.nentries:  # successor is an activity
                            succ_act_idx = succ - ashift  # convert to activity array index
                            if 0 <= succ_act_idx < lqn.nacts and lqn.actphase[succ_act_idx] == 1:
                                raise ValueError(
                                    f"Unsupported replyTo in non-terminal activity: "
                                    f"activity at index {a} has Phase 1 successor"
                                )

        return lqn

    def writeXML(self, filename: str, use_abstract_names: bool = False) -> None:
        """
        Write the layered network to LQNX XML format.

        This method generates an LQNX file compatible with the lqns/lqsim
        command-line tools, matching the MATLAB writeXML implementation.

        Args:
            filename: Path to write the LQNX XML file
            use_abstract_names: If True, use abstract names (P1, T1, E1, A1...)
                              instead of actual element names

        Example:
            >>> model = LayeredNetwork('ClientServer')
            >>> # ... build model ...
            >>> model.writeXML('model.lqnx')
        """
        import xml.etree.ElementTree as ET
        from xml.dom import minidom

        # Build name mapping
        name_map = {}
        if use_abstract_names:
            pctr, tctr, ectr, actr = 1, 1, 1, 1
            for proc in self.processors:
                name_map[proc.name] = f'P{pctr}'
                pctr += 1
                for task in proc.tasks:
                    name_map[task.name] = f'T{tctr}'
                    tctr += 1
                    for entry in task.entries:
                        name_map[entry.name] = f'E{ectr}'
                        ectr += 1
                    for act in task.activities:
                        name_map[act.name] = f'A{actr}'
                        actr += 1
        else:
            for proc in self.processors:
                name_map[proc.name] = proc.name
                for task in proc.tasks:
                    name_map[task.name] = task.name
                    for entry in task.entries:
                        name_map[entry.name] = entry.name
                    for act in task.activities:
                        name_map[act.name] = act.name

        def sched_to_text(sched):
            """Convert sched_strategy to text (handles both enum and int)."""
            if isinstance(sched, SchedStrategy):
                return sched.name.lower()
            elif hasattr(sched, 'value'):
                return str(sched.value).lower()
            else:
                # It's an integer, need to map to name
                for s in SchedStrategy:
                    if s.value == sched:
                        return s.name.lower()
                return 'fcfs'  # Default

        def is_inf_sched(sched):
            """Check if scheduling is INF."""
            if isinstance(sched, SchedStrategy):
                return sched == SchedStrategy.INF
            elif hasattr(sched, 'value'):
                return sched.value == SchedStrategy.INF.value
            else:
                return sched == SchedStrategy.INF.value

        def is_ref_sched(sched):
            """Check if scheduling is REF."""
            if isinstance(sched, SchedStrategy):
                return sched == SchedStrategy.REF
            elif hasattr(sched, 'value'):
                return sched.value == SchedStrategy.REF.value
            else:
                return sched == SchedStrategy.REF.value

        def get_dist_mean(dist):
            """Get mean from distribution (handles dataclass and native distributions)."""
            if dist is None:
                return 0.0
            if hasattr(dist, 'mean') and not callable(dist.mean):
                # Dataclass Distribution with .mean field
                return dist.mean
            elif hasattr(dist, 'getMean'):
                # Native distribution with getMean() method
                return dist.getMean()
            elif hasattr(dist, 'get_mean'):
                return dist.get_mean()
            else:
                return 0.0

        def get_dist_scv(dist):
            """Get SCV from distribution (handles dataclass and native distributions)."""
            if dist is None:
                return 1.0
            if hasattr(dist, 'scv') and not callable(dist.scv):
                # Dataclass Distribution with .scv field
                return dist.scv
            elif hasattr(dist, 'getSCV'):
                # Native distribution with getSCV() method
                return dist.getSCV()
            elif hasattr(dist, 'get_scv'):
                return dist.get_scv()
            else:
                return 1.0

        # Create root element
        root = ET.Element('lqn-model')
        root.set('xmlns:xsi', 'http://www.w3.org/2001/XMLSchema-instance')
        root.set('xsi:noNamespaceSchemaLocation', 'lqn.xsd')
        root.set('name', self.name)

        # Write processors
        for proc in self.processors:
            proc_elem = ET.SubElement(root, 'processor')
            proc_elem.set('name', name_map[proc.name])
            proc_elem.set('scheduling', sched_to_text(proc.sched_strategy))

            if not is_inf_sched(proc.sched_strategy):
                mult = proc.multiplicity
                if np.isinf(mult):
                    mult = 1
                proc_elem.set('multiplicity', str(int(mult)))

            proc_elem.set('speed-factor', '1')

            # Write tasks for this processor
            for task in proc.tasks:
                task_elem = ET.SubElement(proc_elem, 'task')
                task_elem.set('name', name_map[task.name])
                task_elem.set('scheduling', sched_to_text(task.sched_strategy))

                if not is_inf_sched(task.sched_strategy):
                    mult = task.multiplicity
                    if np.isinf(mult):
                        mult = 1000000  # Large number for infinite
                    task_elem.set('multiplicity', str(int(mult)))

                if is_ref_sched(task.sched_strategy) and task.think_time:
                    task_elem.set('think-time', str(get_dist_mean(task.think_time)))

                # Write entries for this task
                for entry in task.entries:
                    entry_elem = ET.SubElement(task_elem, 'entry')
                    entry_elem.set('name', name_map[entry.name])
                    entry_elem.set('type', 'NONE')

                # Write task-activities
                task_acts_elem = ET.SubElement(task_elem, 'task-activities')

                for act in task.activities:
                    act_elem = ET.SubElement(task_acts_elem, 'activity')
                    act_elem.set('name', name_map[act.name])
                    act_elem.set('host-demand-mean', str(get_dist_mean(act.host_demand)))
                    act_elem.set('host-demand-cvsq', str(get_dist_scv(act.host_demand)))
                    act_elem.set('call-order', 'STOCHASTIC')

                    if act.bound_entry:
                        act_elem.set('bound-to-entry', name_map[act.bound_entry.name])

                    # Write calls
                    for target_entry, mean_calls, call_type in act.calls:
                        if call_type == CallType.SYNC:
                            call_elem = ET.SubElement(act_elem, 'synch-call')
                        elif call_type == CallType.ASYNC:
                            call_elem = ET.SubElement(act_elem, 'asynch-call')
                        else:  # FWD
                            call_elem = ET.SubElement(act_elem, 'asynch-call')

                        call_elem.set('dest', name_map[target_entry.name])
                        call_elem.set('calls-mean', str(mean_calls))

                # Write precedences
                for prec in task.precedences:
                    # Handle serial precedences specially - they need N-1 precedence elements for N activities
                    if prec.prec_type == PrecedenceType.SERIAL and prec.activities and len(prec.activities) >= 2:
                        # Serial chain: A → B → C needs two precedence elements: A→B and B→C
                        for i in range(len(prec.activities) - 1):
                            prec_elem = ET.SubElement(task_acts_elem, 'precedence')
                            pre_elem = ET.SubElement(prec_elem, 'pre')
                            act_ref = ET.SubElement(pre_elem, 'activity')
                            act_ref.set('name', name_map[prec.activities[i].name])
                            post_elem = ET.SubElement(prec_elem, 'post')
                            act_ref = ET.SubElement(post_elem, 'activity')
                            act_ref.set('name', name_map[prec.activities[i + 1].name])
                        continue

                    prec_elem = ET.SubElement(task_acts_elem, 'precedence')

                    # Write pre-activities
                    if prec.pre_activities:
                        if len(prec.pre_activities) == 1:
                            pre_elem = ET.SubElement(prec_elem, 'pre')
                        else:
                            if prec.prec_type == PrecedenceType.PARALLEL:
                                pre_elem = ET.SubElement(prec_elem, 'pre-AND')
                            else:
                                pre_elem = ET.SubElement(prec_elem, 'pre-OR')

                        for pre_act in prec.pre_activities:
                            act_ref = ET.SubElement(pre_elem, 'activity')
                            act_ref.set('name', name_map[pre_act.name])

                    # Write post-activities
                    # For Loop precedence, use activities list as post_activities if post_activities is empty
                    post_acts = prec.post_activities
                    if not post_acts and prec.prec_type == PrecedenceType.LOOP and prec.activities:
                        post_acts = prec.activities

                    if post_acts:
                        if len(post_acts) == 1 and prec.prec_type != PrecedenceType.LOOP:
                            post_elem = ET.SubElement(prec_elem, 'post')
                        else:
                            if prec.prec_type == PrecedenceType.PARALLEL:
                                post_elem = ET.SubElement(prec_elem, 'post-AND')
                            elif prec.prec_type == PrecedenceType.CHOICE:
                                post_elem = ET.SubElement(prec_elem, 'post-OR')
                            elif prec.prec_type == PrecedenceType.LOOP:
                                post_elem = ET.SubElement(prec_elem, 'post-LOOP')
                                # Set end attribute to the last activity name
                                post_elem.set('end', name_map[post_acts[-1].name])
                            else:
                                post_elem = ET.SubElement(prec_elem, 'post')

                        # For LOOP, iterate up to but not including the last activity
                        loop_acts_range = post_acts[:-1] if prec.prec_type == PrecedenceType.LOOP and len(post_acts) > 1 else post_acts
                        for i, post_act in enumerate(loop_acts_range):
                            act_ref = ET.SubElement(post_elem, 'activity')
                            act_ref.set('name', name_map[post_act.name])
                            if prec.prec_type == PrecedenceType.CHOICE and i < len(prec.probabilities):
                                act_ref.set('prob', str(prec.probabilities[i]))
                            elif prec.prec_type == PrecedenceType.LOOP:
                                act_ref.set('count', str(prec.count))

                # Write reply-entry elements (for non-reference tasks)
                if task.sched_strategy != SchedStrategy.REF:
                    for entry in task.entries:
                        # Find activities that reply to this entry
                        reply_activities = [act for act in task.activities if act.reply_entry == entry]
                        if reply_activities:
                            reply_entry_elem = ET.SubElement(task_acts_elem, 'reply-entry')
                            reply_entry_elem.set('name', name_map[entry.name])
                            for reply_act in reply_activities:
                                reply_act_elem = ET.SubElement(reply_entry_elem, 'reply-activity')
                                reply_act_elem.set('name', name_map[reply_act.name])

        # Write to file with pretty formatting
        xml_str = ET.tostring(root, encoding='unicode')
        dom = minidom.parseString(xml_str)
        pretty_xml = dom.toprettyxml(indent='  ')

        # Remove extra blank lines
        lines = [line for line in pretty_xml.split('\n') if line.strip()]
        pretty_xml = '\n'.join(lines)

        with open(filename, 'w', encoding='utf-8') as f:
            f.write(pretty_xml)

    def summary(self) -> str:
        """Get a text summary of the layered network structure."""
        lines = [f"Layered Network: {self.name}"]
        lines.append(f"  Processors: {len(self.processors)}")
        lines.append(f"  Tasks: {len(self.tasks)}")
        lines.append(f"  Entries: {len(self.entries)}")
        lines.append(f"  Activities: {len(self.activities)}")

        lines.append("\nProcessors:")
        for proc in self.processors:
            lines.append(f"  - {proc.name} (mult={proc.multiplicity}, sched={proc.sched_strategy.value})")
            for task in proc.tasks:
                lines.append(f"      Task: {task.name}")

        lines.append("\nTasks:")
        for task in self.tasks:
            lines.append(f"  - {task.name} (mult={task.multiplicity}, sched={task.sched_strategy.value})")
            if task.think_time:
                lines.append(f"      Think time: {_get_dist_mean(task.think_time)}")
            for entry in task.entries:
                lines.append(f"      Entry: {entry.name}")

        lines.append("\nActivities and Calls:")
        for act in self.activities:
            lines.append(f"  - {act.name} (demand={_get_dist_mean(act.host_demand)})")
            for entry, mean_calls, call_type in act.calls:
                lines.append(f"      -> {entry.name} ({call_type.value}, calls={mean_calls})")

        return "\n".join(lines)

    def getNodeCount(self) -> int:
        """Get total number of nodes (processors + tasks + entries + activities)."""
        return len(self.processors) + len(self.tasks) + len(self.entries) + len(self.activities)

    get_node_count = getNodeCount

    def getNodes(self) -> list:
        """Get all nodes in the network."""
        return list(self.processors) + list(self.tasks) + list(self.entries) + list(self.activities)

    get_nodes = getNodes

    def getNodeByName(self, name: str):
        """
        Get a node by its name.

        Args:
            name: Name of the node to find

        Returns:
            The node with the given name, or None if not found
        """
        for proc in self.processors:
            if proc.name == name:
                return proc
        for task in self.tasks:
            if task.name == name:
                return task
        for entry in self.entries:
            if entry.name == name:
                return entry
        for act in self.activities:
            if act.name == name:
                return act
        return None

    get_node_by_name = getNodeByName

    def registerNode(self, node):
        """
        Register a node with the network (compatibility method).

        In native mode, nodes are auto-registered when created with the model.
        This method is provided for API compatibility.

        Args:
            node: The node to register (Processor, Task, Entry, or Activity)
        """
        if isinstance(node, Processor) and node not in self.processors:
            self.processors.append(node)
        elif isinstance(node, Task) and node not in self.tasks:
            self.tasks.append(node)
        elif isinstance(node, Entry) and node not in self.entries:
            self.entries.append(node)
        elif isinstance(node, Activity) and node not in self.activities:
            self.activities.append(node)

    register_node = registerNode

    def getHosts(self) -> list:
        """Get all processors (hosts)."""
        return list(self.processors)

    def getTasks(self) -> list:
        """Get all tasks."""
        return list(self.tasks)

    def getEntries(self) -> list:
        """Get all entries."""
        return list(self.entries)

    def getActivities(self) -> list:
        """Get all activities."""
        return list(self.activities)

    def getNodeIndex(self, node) -> int:
        """Get the index of a node."""
        all_nodes = self.getNodes()
        try:
            return all_nodes.index(node)
        except ValueError:
            return -1

    def getNodeNames(self) -> list:
        """Get names of all nodes."""
        return [n.name for n in self.getNodes()]

    def getEnsemble(self):
        """Get the ensemble (returns self for single-model LQNs)."""
        return self

    def getLayers(self) -> list:
        """Get layers (returns list of tasks grouped by layer)."""
        # Build layers based on task call graph depth
        if not self.tasks:
            return []
        # For now, return a simple layering: reference tasks in layer 0, others by depth
        layers = []
        ref_tasks = [t for t in self.tasks if t.sched_strategy == SchedStrategy.REF]
        if ref_tasks:
            layers.append(ref_tasks)
        non_ref_tasks = [t for t in self.tasks if t.sched_strategy != SchedStrategy.REF]
        if non_ref_tasks:
            layers.append(non_ref_tasks)
        return layers

    def getNumberOfLayers(self) -> int:
        """Get number of layers."""
        return len(self.getLayers())

    def getNumberOfModels(self) -> int:
        """Get number of models (returns 1 for single LQN)."""
        return 1

    # Note: getStruct() is defined earlier in the class with proper implementation

    # Snake_case aliases for getters
    get_hosts = getHosts
    get_tasks = getTasks
    get_entries = getEntries
    get_activities = getActivities
    get_node_index = getNodeIndex
    get_node_names = getNodeNames
    get_ensemble = getEnsemble
    get_layers = getLayers
    get_number_of_layers = getNumberOfLayers
    get_number_of_models = getNumberOfModels

    def copy(self) -> 'LayeredNetwork':
        """
        Create a deep copy of this layered network.

        Returns:
            A new LayeredNetwork instance with the same structure
        """
        import copy as copy_module

        # Create new network
        new_model = LayeredNetwork(self.name)

        # Map old objects to new objects
        proc_map = {}
        task_map = {}
        entry_map = {}
        act_map = {}

        # Copy processors
        for proc in self.processors:
            new_proc = new_model.add_processor(proc.name, proc.multiplicity, proc.sched_strategy)
            proc_map[proc] = new_proc

        # Copy tasks (preserve CacheTask type)
        for task in self.tasks:
            if task.processor:
                new_proc = proc_map.get(task.processor)
            else:
                new_proc = None

            # Check if this is a CacheTask
            if isinstance(task, CacheTask):
                new_task = CacheTask(
                    new_model, task.name, task.total_items,
                    task.cache_capacity, task.replacement_strategy,
                    task.multiplicity
                )
                # CacheTask already adds itself to model.tasks
                new_model.tasks.remove(new_task)  # Remove auto-added
            else:
                # Regular Task uses positional args: (name, multiplicity, sched_strategy)
                new_task = Task(task.name, task.multiplicity, task.sched_strategy)

            if new_proc:
                new_task.on(new_proc)
            if task.think_time:
                new_task.set_think_time(_get_dist_mean(task.think_time))
            new_model.tasks.append(new_task)
            task_map[task] = new_task

        # Copy entries (preserve ItemEntry type)
        for entry in self.entries:
            if entry.task:
                new_task = task_map.get(entry.task)
            else:
                new_task = None

            # Check if this is an ItemEntry
            if isinstance(entry, ItemEntry):
                new_entry = ItemEntry(
                    new_task if new_task else None,
                    entry.name,
                    entry.total_items,
                    copy_module.deepcopy(entry.access_prob) if entry.access_prob else None
                )
                # ItemEntry may already be registered
                if new_entry in new_model.entries:
                    new_model.entries.remove(new_entry)
            else:
                # Entry uses positional arg: (name)
                new_entry = Entry(entry.name)
                if new_task:
                    new_entry.on(new_task)

            new_model.entries.append(new_entry)
            entry_map[entry] = new_entry

        # Copy activities
        for act in self.activities:
            if act.task:
                new_task = task_map.get(act.task)
            else:
                new_task = None
            # Activity uses positional args: (name, host_demand)
            new_act = Activity(act.name, copy_module.deepcopy(act.host_demand))
            if new_task:
                new_act.on(new_task)
            new_model.activities.append(new_act)
            act_map[act] = new_act

            # Copy bound entry
            if act.bound_entry:
                new_act.bound_entry = entry_map.get(act.bound_entry)

            # Copy reply entry
            if act.reply_entry:
                new_act.reply_entry = entry_map.get(act.reply_entry)

        # Copy calls (after all activities are created)
        for act in self.activities:
            new_act = act_map[act]
            for target_entry, mean_calls, call_type in act.calls:
                new_target = entry_map.get(target_entry)
                if new_target:
                    new_act.calls.append((new_target, mean_calls, call_type))

        # Copy precedences for each task
        for task in self.tasks:
            new_task = task_map[task]
            for prec in task.precedences:
                new_activities = [act_map.get(a) for a in prec.activities if a in act_map]
                new_activities = [a for a in new_activities if a is not None]
                if new_activities:
                    new_prec = ActivityPrecedence(
                        prec_type=prec.prec_type,
                        activities=new_activities,
                        probabilities=list(prec.probabilities) if prec.probabilities else [],
                        count=prec.count
                    )
                    new_task.precedences.append(new_prec)

        return new_model

    # Provide obj property for compatibility with wrapper-based code
    @property
    def obj(self):
        """Return self for compatibility with wrapper code that accesses .obj"""
        return self

    @classmethod
    def parse_xml(cls, filename: str, verbose: bool = False) -> 'LayeredNetwork':
        """
        Parse an LQNX XML file and create a LayeredNetwork model.

        This method parses layered queueing network XML files in LQNX format
        and constructs the corresponding Python model.

        Args:
            filename: Path to the LQNX XML file
            verbose: If True, print parsing progress

        Returns:
            LayeredNetwork model

        Example:
            >>> model = LayeredNetwork.parse_xml('model.lqnx')
        """
        import xml.etree.ElementTree as ET
        import os

        # Handle relative paths
        if not os.path.isabs(filename):
            if not os.path.exists(filename):
                for base in ['.', os.getcwd()]:
                    full_path = os.path.join(base, filename)
                    if os.path.exists(full_path):
                        filename = full_path
                        break

        if not os.path.exists(filename):
            raise FileNotFoundError(f"File cannot be found: {filename}")

        tree = ET.parse(filename)
        root = tree.getroot()

        if verbose:
            print(f"Parsing LQN file: {filename}")

        model_name = root.get('name', filename.replace('_', r'\_'))
        model = cls(model_name)

        entry_map = {}
        activity_map = {}

        for proc_elem in root.findall('.//processor'):
            proc_name = proc_elem.get('name', '')
            scheduling = proc_elem.get('scheduling', 'fcfs').upper()

            mult_str = proc_elem.get('multiplicity', '1')
            if mult_str.lower() == 'inf':
                multiplicity = float('inf')
            else:
                try:
                    multiplicity = float(mult_str)
                except ValueError:
                    multiplicity = 1.0

            if scheduling.upper() == 'INF':
                multiplicity = float('inf')
                sched = SchedStrategy.INF
            else:
                sched = cls._parse_sched_strategy(scheduling)

            processor = model.add_processor(proc_name, multiplicity, sched)

            for task_elem in proc_elem.findall('./task'):
                task_name = task_elem.get('name', '')
                task_sched = task_elem.get('scheduling', 'fcfs').upper()

                task_mult_str = task_elem.get('multiplicity', '1')
                if task_mult_str.lower() == 'inf':
                    task_mult = float('inf')
                else:
                    try:
                        task_mult = float(task_mult_str)
                    except ValueError:
                        task_mult = 1.0

                if task_sched.upper() == 'INF':
                    task_mult = float('inf')
                    task_sched_enum = SchedStrategy.INF
                else:
                    task_sched_enum = cls._parse_sched_strategy(task_sched)

                task = model.add_task(task_name, task_mult, task_sched_enum, processor)

                think_time_str = task_elem.get('think-time', '0')
                try:
                    think_time = float(think_time_str)
                    if think_time > 0:
                        task.set_think_time(think_time)
                except ValueError:
                    pass

                for entry_elem in task_elem.findall('./entry'):
                    entry_name = entry_elem.get('name', '')
                    entry = model.add_entry(entry_name, task)
                    entry_map[entry_name] = entry

                for ta_elem in task_elem.findall('./task-activities'):
                    for act_elem in ta_elem.findall('./activity'):
                        act_name = act_elem.get('name', '')
                        demand_mean = float(act_elem.get('host-demand-mean', '0'))
                        demand_scv = float(act_elem.get('host-demand-cvsq', '1.0'))

                        if demand_mean <= 0:
                            demand = Distribution(mean=0.0, scv=1.0)
                        else:
                            demand = Distribution(mean=demand_mean, scv=demand_scv)

                        activity = model.add_activity(act_name, demand, task)
                        activity_map[act_name] = activity

                        bound_entry_name = act_elem.get('bound-to-entry', '')
                        if bound_entry_name and bound_entry_name in entry_map:
                            activity.bound_to(entry_map[bound_entry_name])

                        for call_elem in act_elem.findall('./synch-call'):
                            dest = call_elem.get('dest', '')
                            mean_calls = float(call_elem.get('calls-mean', '1'))
                            if not hasattr(activity, '_pending_calls'):
                                activity._pending_calls = []
                            activity._pending_calls.append((dest, mean_calls, CallType.SYNC))

                        for call_elem in act_elem.findall('./asynch-call'):
                            dest = call_elem.get('dest', '')
                            mean_calls = float(call_elem.get('calls-mean', '1'))
                            if not hasattr(activity, '_pending_calls'):
                                activity._pending_calls = []
                            activity._pending_calls.append((dest, mean_calls, CallType.ASYNC))

                    for prec_elem in ta_elem.findall('./precedence'):
                        pre_acts = []
                        post_acts = []
                        pre_type = 'pre'
                        post_type = 'post'

                        for pre_tag in ['pre', 'pre-AND', 'pre-OR']:
                            pre_elem = prec_elem.find(f'./{pre_tag}')
                            if pre_elem is not None:
                                pre_type = pre_tag
                                for act_ref in pre_elem.findall('./activity'):
                                    act_name = act_ref.get('name', '')
                                    if act_name in activity_map:
                                        pre_acts.append(activity_map[act_name])
                                break

                        for post_tag in ['post', 'post-AND', 'post-OR', 'post-LOOP']:
                            post_elem = prec_elem.find(f'./{post_tag}')
                            if post_elem is not None:
                                post_type = post_tag
                                for act_ref in post_elem.findall('./activity'):
                                    act_name = act_ref.get('name', '')
                                    if act_name in activity_map:
                                        post_acts.append(activity_map[act_name])
                                break

                        if pre_acts and post_acts:
                            if post_type == 'post-AND':
                                prec = ActivityPrecedence.AndFork(pre_acts[0], post_acts)
                            elif post_type == 'post-OR':
                                probs = []
                                for act_ref in post_elem.findall('./activity'):
                                    prob_str = act_ref.get('prob', '')
                                    if prob_str:
                                        probs.append(float(prob_str))
                                    else:
                                        probs.append(1.0 / len(post_acts))
                                prec = ActivityPrecedence.OrFork(pre_acts[0], post_acts, probs)
                            elif pre_type == 'pre-AND':
                                prec = ActivityPrecedence.AndJoin(pre_acts, post_acts[0])
                            elif pre_type == 'pre-OR':
                                prec = ActivityPrecedence.OrJoin(pre_acts, post_acts[0])
                            else:
                                prec = ActivityPrecedence.Serial(pre_acts + post_acts)
                            task.add_precedence(prec)

                    for reply_elem in ta_elem.findall('./reply-entry'):
                        reply_entry_name = reply_elem.get('name', '')
                        if reply_entry_name in entry_map:
                            reply_entry = entry_map[reply_entry_name]
                            for reply_act_elem in reply_elem.findall('./reply-activity'):
                                reply_act_name = reply_act_elem.get('name', '')
                                if reply_act_name in activity_map:
                                    activity_map[reply_act_name].replies_to(reply_entry)

        for activity in model.activities:
            if hasattr(activity, '_pending_calls'):
                for dest, mean_calls, call_type in activity._pending_calls:
                    if dest in entry_map:
                        if call_type == CallType.SYNC:
                            activity.synch_call(entry_map[dest], mean_calls)
                        else:
                            activity.asynch_call(entry_map[dest], mean_calls)
                delattr(activity, '_pending_calls')

        return model

    parseXML = parse_xml
    readXML = parse_xml
    load = parse_xml

    @staticmethod
    def _parse_sched_strategy(sched_str: str) -> SchedStrategy:
        """Parse scheduling strategy string to enum."""
        sched_upper = sched_str.upper()
        if sched_upper in ('FCFS', 'FIFO'):
            return SchedStrategy.FCFS
        elif sched_upper == 'PS':
            return SchedStrategy.PS
        elif sched_upper == 'INF':
            return SchedStrategy.INF
        elif sched_upper == 'REF':
            return SchedStrategy.REF
        elif sched_upper == 'HOL':
            return SchedStrategy.HOL
        elif sched_upper in ('LCFS', 'LIFO'):
            return SchedStrategy.LCFS
        else:
            return SchedStrategy.FCFS


# Convenience aliases for compatibility
Processor = Processor
Task = Task
FunctionTask = Task  # Alias for backward compatibility (FaaS/serverless modeling)
Entry = Entry
Activity = Activity
LayeredNetwork = LayeredNetwork


__all__ = [
    'LayeredNetwork',
    'Processor',
    'Task',
    'Entry',
    'Activity',
    'CacheTask',
    'ItemEntry',
    'LayeredNetworkStruct',
    'ActivityPrecedence',
    'PrecedenceType',
    'SchedStrategy',
    'CallType',
    'Distribution',
    'ReplacementStrategy',
    # Convenience aliases
    'LayeredNetwork',
    'Processor',
    'Task',
    'FunctionTask',
    'Entry',
    'Activity',
    'CacheTask',
    'ItemEntry',
]
