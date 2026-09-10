"""
Native Python implementation of Layered Queueing Network (LQN) models.

This module provides pure Python classes for defining and analyzing
layered queueing networks.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any, Union
from enum import Enum, IntEnum

from .constants import SchedStrategy, RoutingStrategy


class LayeredNetworkElement(IntEnum):
    """Element types in layered queueing networks (matches MATLAB enum values)."""
    PROCESSOR = 0
    TASK = 1
    ENTRY = 2
    ACTIVITY = 3
    CALL = 4
from .lang.base import ReplacementStrategy
from .distributions import Immediate, Exp, Bernoulli, Geometric


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


def _get_dist_mean(dist) -> float:
    """Get mean from distribution (handles dataclass and native distributions)."""
    if dist is None:
        return 0.0
    # Handle numeric types directly
    if isinstance(dist, (int, float)):
        return float(dist)
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


def _call_count_dist(mean_calls):
    """Distribution of the number of calls issued per invocation.

    A mean below 1 is a call that either happens or does not, hence Bernoulli.
    Geometric(1/m) is undefined there: its parameter would exceed 1 and its SCV
    (1-p) would come out negative.
    """
    m = float(mean_calls)
    if not np.isfinite(m) or m <= 0:
        return Immediate()
    if m < 1.0:
        return Bernoulli(m)
    return Geometric(1.0 / m)


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
    # Quorum count k of an AND-join, as a one-element array. None means the join waits
    # for all its predecessors, which is the LQN default.
    pre_params: Optional[Any] = None

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
    def AndJoin(pre_acts: List['Activity'], post_act: 'Activity',
                quorum: Optional[int] = None) -> 'ActivityPrecedence':
        """
        Create an AND-join precedence (synchronization).

        Post-activity starts once quorum of the pre-activities have completed. With no
        quorum the join waits for ALL of them, which is the LQN default.
        Used together with AndFork() to model parallel execution.

        Args:
            pre_acts: List of Activity objects to synchronize on
            post_act: Activity that executes after synchronization
            quorum: Number k of pre-activities required to fire the join. Defaults to
                None, meaning all of them. Values outside [1, len(pre_acts)] are ignored.

        Returns:
            ActivityPrecedence object

        Example:
            >>> task.add_precedence(ActivityPrecedence.AndJoin([branch1, branch2], end))
            >>> task.add_precedence(ActivityPrecedence.AndJoin([b1, b2, b3], end, quorum=2))
        """
        return ActivityPrecedence(
            prec_type=PrecedenceType.PARALLEL,
            pre_activities=list(pre_acts),
            post_activities=[post_act],
            pre_params=None if quorum is None else np.array([quorum])
        )

    # Python snake_case alias
    @staticmethod
    def and_join(pre_acts: List['Activity'], post_act: 'Activity',
                 quorum: Optional[int] = None) -> 'ActivityPrecedence':
        """Python snake_case alias for AndJoin()."""
        return ActivityPrecedence.AndJoin(pre_acts, post_act, quorum)

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


def _dist_from_mean_scv(mean, scv):
    """Materialize a native distribution from (mean, scv), as MATLAB parseXML
    does; the internal Distribution placeholder is not a recognized process
    type downstream (SolverLN layer construction resolves by class name)."""
    from .distributions import Exp, Erlang, HyperExp, Immediate
    if mean <= 0:
        return Immediate()
    if abs(scv - 1.0) < 1e-12:
        return Exp(1.0 / mean)
    if scv < 1.0:
        return Erlang.fit_mean_and_scv(mean, scv)
    return HyperExp.fit_mean_and_scv(mean, scv)


# --- LINE .lqnx dialect: cache, item entry, setup / delay-off -----------------
#
# See scratchpad LQNX_CACHE_SPEC.md. The stock LQN schema has no element for any
# of these, so a model carrying them used to export as a DIFFERENT model. The
# field set is the JSON interchange's (linemodel_save.m), so the two transports
# carry the same information and a model round-trips through either.


def _callgroup_to_lqnx(strategy):
    """RoutingStrategy -> the wire enum name, spelled as the JSON interchange
    spells it. Only the two strategies a call group can be built with are named:
    WRROBIN would need per-target weights the group API does not take, and the
    remaining strategies are not dispatch policies at all, so an unnamed one is
    an error rather than a silent PROB."""
    if strategy == RoutingStrategy.RROBIN:
        return 'RROBIN'
    if strategy == RoutingStrategy.JSQ:
        return 'JSQ'
    raise ValueError('Call groups carry RROBIN or JSQ; routing strategy %s cannot '
                     'be written to .lqnx' % strategy)


def _callgroup_from_lqnx(name, act_name):
    """Wire enum name -> RoutingStrategy, the inverse of _callgroup_to_lqnx."""
    key = (name or '').strip().upper()
    if key == 'RROBIN':
        return RoutingStrategy.RROBIN
    if key == 'JSQ':
        return RoutingStrategy.JSQ
    raise ValueError('Activity "%s" declares a call group with an unrecognized '
                     'strategy "%s"; the dialect spells them RROBIN and JSQ'
                     % (act_name, name))


def _parse_call_groups(act_elem, activity):
    """Read the LINE dialect <call-group> children of an activity element.

    The member calls are ordinary synch-call elements and have already been read
    into _pending_calls, so only the grouping is recorded here; issuing them again
    would double the call rate. Target names are resolved to Entry objects with
    the calls, once every entry exists.
    """
    for grp_elem in act_elem.findall('./call-group'):
        strategy = _callgroup_from_lqnx(grp_elem.get('strategy', ''), activity.name)
        dests = [d.get('name', '') for d in grp_elem.findall('./dest')]
        if not hasattr(activity, '_pending_call_groups'):
            activity._pending_call_groups = []
        activity._pending_call_groups.append((strategy, dests))


def _cap_list(cap):
    """Per-list capacities as a list. `itemLevelCap` is an ARRAY in MATLAB and a
    multi-list cache is the normal case, so a scalar is the one-list special
    case rather than the other way round."""
    if cap is None:
        return []
    if isinstance(cap, (int, float)):
        return [int(cap)]
    try:
        return [int(c) for c in cap]
    except TypeError:
        return [int(cap)]


def _write_cache_elem(task_elem, task):
    """<cache items= replacement= retrieval=> with one <level capacity=> per list."""
    if type(task).__name__ != 'CacheTask':
        return
    import xml.etree.ElementTree as ET
    cache = ET.SubElement(task_elem, 'cache')
    cache.set('items', str(int(task.total_items)))
    rs = task.replacement_strategy
    # The wire spelling is the enum NAME, the same mapping the JSON interchange
    # uses (linemodel_io: `rs.name`); no second spelling table.
    cache.set('replacement', rs.name if hasattr(rs, 'name') else str(rs))
    if getattr(task, 'retrieval', False):
        cache.set('retrieval', 'true')
    for cap in _cap_list(task.cache_capacity):
        ET.SubElement(cache, 'level').set('capacity', str(int(cap)))


def _write_setup_elems(task_elem, task):
    """<setup mean= scv=/> and <delay-off mean= scv=/>, emitted only when set."""
    import xml.etree.ElementTree as ET
    for tag, dist in (('setup', getattr(task, 'setup_time', None)),
                      ('delay-off', getattr(task, 'delay_off_time', None))):
        if dist is None:
            continue
        mean = _get_dist_mean(dist)
        if mean <= 0:
            continue
        elem = ET.SubElement(task_elem, tag)
        elem.set('mean', repr(float(mean)))
        elem.set('scv', repr(float(_get_dist_scv(dist))))


def _write_item_entry_elem(entry_elem, entry):
    """<item-entry cardinality=> with an <access-popularity> of the class name and
    its constructor parameters, in order, one <parameter value=> each.

    A VECTOR constructor parameter is flattened, which is unambiguous here because
    the reader knows the cardinality: a DiscreteSampler written with n parameters
    is p alone, with 2n it is p then x. Flagged to the coordinator as the one
    place the spec's flat <parameter> list needs a stated rule.
    """
    if type(entry).__name__ != 'ItemEntry':
        return
    import xml.etree.ElementTree as ET
    ie = ET.SubElement(entry_elem, 'item-entry')
    ie.set('cardinality', str(int(entry.total_items)))
    dist = getattr(entry, 'access_prob', None)
    if dist is None:
        return
    params = _popularity_params(dist)
    ap = ET.SubElement(ie, 'access-popularity')
    ap.set('name', type(dist).__name__)
    for value in params:
        ET.SubElement(ap, 'parameter').set('value', repr(float(value)))


def _popularity_params(dist):
    """Constructor parameters of an access-popularity distribution, in order.

    A class this encoding cannot express is REFUSED BY NAME. That is not the
    warning-and-drop this whole change exists to remove: the three constructs the
    dialect must carry (cache, item entry, setup/delay-off) are carried
    unconditionally, and writing a popularity element without the parameters that
    define it would reproduce exactly the bug being fixed -- a file that loads as
    a different access law.
    """
    import numpy as _np
    name = type(dist).__name__
    if name == 'DiscreteSampler':
        probs = list(_np.asarray(dist._probs).flatten())
        values = getattr(dist, '_values', None)
        # The default support 1..n is reconstructible, so it is omitted and the
        # common uniform sampler writes exactly the n probabilities of the spec
        # example; a non-default support is appended, giving 2n. The reader
        # splits on the cardinality already declared on <item-entry>.
        if values is not None:
            xs = list(_np.asarray(values).flatten())
            if len(xs) == len(probs) and not _np.allclose(xs, _np.arange(1, len(xs) + 1)):
                return probs + xs
        return probs
    if name == 'Zipf':
        # Constructor order is (s, n). NOTE for the other codebases: MATLAB's
        # Zipf stores params 1=p, 2=x, 3=s, 4=n, so "constructor order" there
        # means params 3 and 4 -- writing p and x would be lossy, since s cannot
        # be recovered from them. Python's Zipf holds `_s`/`_n` directly, so the
        # trap does not arise here, but the WIRE form is the same two values.
        return [float(getattr(dist, 's', getattr(dist, '_s', 1.0))),
                float(getattr(dist, 'n', getattr(dist, '_n', 1)))]
    raise RuntimeError(
        "the .lqnx access-popularity encoding carries a DiscreteSampler or a Zipf, "
        "and this ItemEntry uses a '%s'. Writing it as anything else would load as a "
        "different access law, so the export is refused by name. Use model.json "
        "(line_solver.io.save_model), whose distribution encoding is general."
        % name)


def _popularity_from_params(name, params, cardinality):
    """Rebuild an access-popularity distribution from its written parameters.

    THE `name` ATTRIBUTE SELECTS THE RULE, so there is no ambiguity between a
    Zipf and a 2-item DiscreteSampler. A parameter count the named class cannot
    take is REFUSED BY NAME rather than guessed: silently substituting a uniform
    law would load as a different access law, which is the defect this dialect
    exists to remove.
    """
    import numpy as _np
    from .distributions import DiscreteSampler, Zipf
    n = int(cardinality) if cardinality else 0
    if name == 'DiscreteSampler':
        # The cardinality split: n parameters is p over the default support
        # 1..n, 2n is p followed by an explicit support x.
        if n > 0 and len(params) == n:
            return DiscreteSampler(_np.array(params))
        if n > 0 and len(params) == 2 * n:
            return DiscreteSampler(_np.array(params[:n]), _np.array(params[n:]))
        if n <= 0 and params:
            return DiscreteSampler(_np.array(params))
        raise RuntimeError(
            "access-popularity 'DiscreteSampler' carries %d parameters, which is neither "
            "the cardinality %d (probabilities over the default support) nor twice it "
            "(probabilities followed by an explicit support); the file is malformed."
            % (len(params), n))
    if name == 'Zipf':
        # Two parameters, s then n. The cardinality split does NOT apply here.
        if len(params) != 2:
            raise RuntimeError(
                "access-popularity 'Zipf' takes exactly two parameters, the exponent s "
                "and the support size n, and this one carries %d; the file is malformed."
                % len(params))
        return Zipf(float(params[0]), int(params[1]))
    raise RuntimeError(
        "access-popularity names the class '%s', which the .lqnx dialect does not carry "
        "(it carries DiscreteSampler and Zipf). Reading it as anything else would load a "
        "different access law." % name)


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
        # Call groups dispatched by a routing strategy instead of independently.
        # Each element is (RoutingStrategy, [entry, ...]); the calls themselves
        # stay in self.calls so every consumer that ignores dispatch order still
        # sees the same aggregate call means.
        self.call_groups = []
        self.think_time = 0.0  # Activity-level think time (LQNX think-time attribute)
        self.phase = 1  # Phase number (1 or 2), default=1

    def _convert_demand(self, demand):
        """Normalize the host-demand argument.

        Native LINE distributions (Exp, Erlang, HyperExp, APH, Coxian, MAP, ...)
        are preserved verbatim so their SCV / phase-type structure survives into
        getStruct (lqn.hostdem_proc) and reaches the layer models — required for
        moment2/moment3 and any SCV-aware layer solver. Only bare Distribution
        placeholders and scalars are passed through unchanged; nothing is
        collapsed to an exponential of matching mean here.
        """
        if demand is None:
            return None
        # Preserve any object exposing a mean (native distributions, internal
        # Distribution dataclass) so higher moments are not discarded.
        if isinstance(demand, Distribution) or hasattr(demand, 'getMean') or hasattr(demand, 'get_mean'):
            return demand
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

    def synch_call_rrobin(self, entries, mean_calls: float = 1.0) -> 'Activity':
        """Dispatch synchronous calls round-robin over a set of target entries.

        MEAN_CALLS is the total mean number of calls the activity issues per
        invocation; successive calls go to the targets in cyclic order, so each
        target receives MEAN_CALLS/len(ENTRIES) of them. The probabilistic model
        with the same per-target means is the ungrouped equivalent: what the
        group adds is the deterministic interleaving, not a different call rate.

        Only the squashed ('flat') layering can represent this, because under
        'srvn' the targets never share a submodel. See SolverLN._assert_call_groups.
        """
        return self._add_call_group(RoutingStrategy.RROBIN, entries, mean_calls,
                                    'synch_call_rrobin')

    def synch_call_jsq(self, entries, mean_calls: float = 1.0) -> 'Activity':
        """Dispatch synchronous calls to the least loaded of a set of target entries.

        Same contract as SYNCH_CALL_RROBIN, with the cyclic pointer replaced by
        join-the-shortest-queue: each call goes to the target task whose station
        holds the fewest jobs at dispatch time, ties split uniformly. The
        probabilistic twin with MEAN_CALLS/len(ENTRIES) per target is again the
        ungrouped equivalent.

        Only the squashed ('flat') layering can represent this, and only a layer
        solver with state-dependent routing can honour it; see
        SolverLN._assert_call_groups.
        """
        return self._add_call_group(RoutingStrategy.JSQ, entries, mean_calls,
                                    'synch_call_jsq')

    def _add_call_group(self, strategy, entries, mean_calls, caller: str) -> 'Activity':
        """Record a routed call group and its per-target call means."""
        if entries is None or len(entries) < 2:
            raise ValueError('%s needs at least two target entries' % caller)
        share = float(mean_calls) / len(entries)
        for entry in entries:
            self.calls.append((entry, share, CallType.SYNC))
        return self.record_call_group(strategy, entries)

    def record_call_group(self, strategy, entries) -> 'Activity':
        """Record the grouping of synchronous calls this activity ALREADY declares.

        _ADD_CALL_GROUP issues the member calls and then records them; the .lqnx
        reader has read them back as ordinary synch-call elements, so it records
        the grouping alone and must not issue them a second time.
        """
        if entries is None or len(entries) < 2:
            raise ValueError('A call group needs at least two target entries')
        self.call_groups.append((strategy, list(entries)))
        return self

    def replies_to(self, entry: 'Entry') -> 'Activity':
        """Mark this activity as replying to an entry."""
        self.reply_entry = entry
        return self

    def setPhase(self, phase_num: int) -> 'Activity':
        """Set the phase number for this activity (1, 2 or 3).

        Phase 1 is the default; phases 2 and 3 mark post-reply activities whose
        demand is incurred after the entry has replied. The range is 1..3
        because lqn-core.xsd bounds the phase attribute there, and every
        consumer of lqn.actphase tests phase > 1, so 3 is served as 2 is.
        """
        if not isinstance(phase_num, (int, float)) or phase_num < 1 or phase_num > 3:
            raise ValueError('Phase number must be 1, 2 or 3')
        self.phase = int(phase_num)
        return self

    def set_phase(self, phase_num: int) -> 'Activity':
        """Set the phase number (snake_case alias for setPhase)."""
        return self.setPhase(phase_num)

    def getPhase(self) -> int:
        """Get the phase number (1..3)."""
        return self.phase

    def get_phase(self) -> int:
        """Get the phase number (snake_case alias for getPhase)."""
        return self.getPhase()

    def setThinkTime(self, think_time) -> 'Activity':
        """Set an activity-level think time, separate from the host demand and
        from the task-level think time. Accepts a numeric mean (converted to an
        exponential) or any Distribution."""
        if isinstance(think_time, (int, float)):
            self.think_time = _dist_from_mean_scv(float(think_time), 1.0)
        else:
            self.think_time = think_time
        return self

    def set_think_time(self, think_time) -> 'Activity':
        """snake_case alias for setThinkTime."""
        return self.setThinkTime(think_time)

    def setHostDemand(self, value) -> 'Activity':
        """Set the mean host demand. A numeric mean is converted to an
        exponential distribution (SCV=1); a Distribution is stored verbatim.
        Used by the LQN parameter identification routines (see infer_lqn)."""
        if isinstance(value, (int, float)):
            self.host_demand = _dist_from_mean_scv(float(value), 1.0)
        else:
            self.host_demand = self._convert_demand(value)
        return self

    def set_host_demand(self, value) -> 'Activity':
        """snake_case alias for setHostDemand."""
        return self.setHostDemand(value)

    def getHostDemand(self):
        """Get host demand distribution."""
        return self.host_demand

    def getHostDemandMean(self) -> float:
        """Get host demand mean (handles native distributions and the internal
        Distribution dataclass)."""
        if self.host_demand is None:
            return 0.0
        return _get_dist_mean(self.host_demand)

    def getHostDemandSCV(self) -> float:
        """Get host demand SCV (handles native distributions and the internal
        Distribution dataclass)."""
        if self.host_demand is None:
            return 1.0
        return _get_dist_scv(self.host_demand)

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
        if self.think_time is None:
            return 0.0
        return _get_dist_mean(self.think_time)

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

    # camelCase aliases, under the names MATLAB and the JAR use, so an LQN
    # script transliterates unchanged. The setters already had them; the
    # activity-graph builders did not, which is what a caller reaches for first.
    boundTo = bound_to
    synchCall = synch_call
    asynchCall = asynch_call
    synchCallRRobin = synch_call_rrobin
    synchCallJSQ = synch_call_jsq
    repliesTo = replies_to
    recordCallGroup = record_call_group


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

    def set_arrival(self, arrival) -> 'Entry':
        """Set arrival distribution for open arrivals (snake_case alias for setArrival)."""
        return self.setArrival(arrival)

    def setArrival(self, arrival) -> 'Entry':
        """Set arrival distribution for open arrivals."""
        self._arrival = arrival
        return self

    def add_forwarding(self, target_entry: 'Entry', prob: float = 1.0) -> 'Entry':
        """Add forwarding to another entry (snake_case alias for addForwarding)."""
        return self.addForwarding(target_entry, prob)

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
    forward = addForwarding  # MATLAB/JAR API name


class AdmissionConstrained:
    """
    Admission constraints on the layer station of a Task or a Processor.

    A server declares ``A*n <= b`` on the station that represents it in its
    layer, where ``n`` counts the jobs in service or queueing there. Only a Task
    or a Processor becomes a server station, so only those carry constraints.
    """

    @property
    def lincon_a(self):
        """Positional constraint matrix, or None."""
        return getattr(self, '_lincon_a', None)

    @property
    def lincon_b(self):
        """Positional capacity vector, or None."""
        return getattr(self, '_lincon_b', None)

    @property
    def lincon_rows(self):
        """Rows declared by operand name, as (names, coeffs, cap) tuples."""
        if not hasattr(self, '_lincon_rows'):
            self._lincon_rows = []
        return self._lincon_rows

    def addConstraint(self, operands, coeffs=None, cap=None):
        """
        Append one admission constraint row naming its operands, so the meaning
        does not depend on declaration order::

            t2.addConstraint([e2, e3], [1, 1], 2)  # n(E2) + n(E3) <= 2
            t2.addConstraint(e3, 1, 1)             # n(E3) <= 1

        Operands are the entries of a Task, or the tasks of a Processor, given as
        objects, names, or a mix. Names are resolved against the model in
        LayeredNetwork.getStruct, where an operand that does not belong to this
        server is an error rather than a silent mis-mapping.
        """
        names = _operand_names(operands)
        n = len(names)
        if coeffs is None:
            coeffs = np.ones(n)
        coeffs = np.atleast_1d(np.asarray(coeffs, dtype=float)).ravel()
        if coeffs.size == 1 and n > 1:
            coeffs = np.full(n, coeffs[0])
        if coeffs.size != n:
            raise ValueError(f"Admission constraint has {n} operands but {coeffs.size} coefficients.")
        if not np.all(np.isfinite(coeffs)) or np.any(coeffs < 0):
            raise ValueError("Admission constraint coefficients must be finite and non-negative.")
        if np.all(coeffs == 0):
            raise ValueError("Admission constraint has all-zero coefficients, which constrains nothing.")
        if len(set(names)) != n:
            raise ValueError("Admission constraint names the same operand more than once; give it a single combined coefficient instead.")
        if cap is None or not np.isscalar(cap) or not np.isfinite(cap) or cap < 1:
            raise ValueError("Admission constraint capacity must be a finite scalar of at least 1.")
        self.lincon_rows.append((names, coeffs, float(cap)))
        return self

    def setConstraint(self, A, b):
        """
        Raw form of addConstraint, for programmatic construction. Columns of A
        are indexed positionally by the entries of a Task, or by the tasks of a
        Processor, in declaration order, so the mapping shifts if an entry is
        added later; prefer addConstraint, which names its operands. Only the
        column count is checked, in LayeredNetwork.getStruct, since entries may
        be added after this call. Rows from both forms are concatenated.
        """
        A = np.atleast_2d(np.asarray(A, dtype=float))
        b = np.asarray(b, dtype=float).ravel()
        if A.size == 0 or b.size == 0:
            raise ValueError("Constraint matrix A and capacity vector b must be non-empty.")
        if A.ndim > 2:
            raise ValueError("Constraint matrix A must be two-dimensional.")
        if A.shape[0] != b.size:
            raise ValueError("A and b must have matching number of rows.")
        if not np.all(np.isfinite(A)) or np.any(A < 0):
            raise ValueError("Constraint matrix A must be finite and non-negative.")
        if not np.all(np.isfinite(b)) or np.any(b < 1):
            raise ValueError("Capacity vector b must be finite and at least 1.")
        if np.any(np.all(A == 0, axis=1)):
            raise ValueError("Constraint matrix A has an all-zero row, which constrains nothing.")
        self._lincon_a = A
        self._lincon_b = b
        return self

    def getLinearConstraints(self):
        """Positional constraint pair declared on this element, before name resolution."""
        return self.lincon_a, self.lincon_b

    def hasLinearConstraints(self):
        """Whether this element declares any admission constraint, in either form."""
        return (self.lincon_a is not None and self.lincon_b is not None) or bool(self.lincon_rows)

    add_constraint = addConstraint
    set_constraint = setConstraint
    get_linear_constraints = getLinearConstraints
    has_linear_constraints = hasLinearConstraints


class RateDependent:
    """
    Service-rate dependences on the layer station of a Task or a Processor.

    A server declares how the rate of the station that represents it in its
    layer scales with the jobs held there: with the total population
    (``setLoadDependence``), with the per-operand population in product form
    (``setClassDependence``), or jointly and outside product form
    (``setJointDependence``). An operand is task j of a Processor or entry j of
    a Task, in the same order as the columns of ``setConstraint``.
    """

    @property
    def lld_scaling(self):
        """Load-dependence vector alpha(n), or None."""
        return getattr(self, '_lld_scaling', None)

    @property
    def lcd_scaling(self):
        """Class-dependence handle beta(n) over this server's operands, or None."""
        return getattr(self, '_lcd_scaling', None)

    @property
    def lcd_scaling_peak(self):
        """Per-operand peak rate scaling of the class dependence, or None."""
        return getattr(self, '_lcd_scaling_peak', None)

    @property
    def ljd_scaling(self):
        """Joint-dependence handle eta(n) over this server's operands, or None."""
        return getattr(self, '_ljd_scaling', None)

    @property
    def ljd_scaling_peak(self):
        """Per-operand peak rate scaling of the joint dependence, or None."""
        return getattr(self, '_ljd_scaling_peak', None)

    @property
    def server_pools(self):
        """Declared compatibility pools, a list of dicts, empty when absent."""
        return getattr(self, '_server_pools', [])

    def addServerType(self, server_type):
        """
        Declare one pool of ``server_type.get_num_of_servers()`` identical
        servers, each running at ``server_type.get_rate()``, eligible only for
        the operands in ``server_type.get_compatible_classes()``::

            P1.addServerType(ServerType('Fast', 2, [T2]))
            P1.addServerType(ServerType('Shared', 1, [T2, T3]))

        The operands are the tasks of a Processor, or the entries of a Task,
        given as objects or names. They are resolved against the model in
        getStruct, where an operand that does not belong to this server is an
        error rather than a silent mis-mapping, exactly as for addConstraint.

        SolverLN lowers the whole declaration to the activated-server rate of
        sn_compat_rate, carried onto the layer station as a joint dependence, so
        the pools are an APPROXIMATION in a layer for the same reason
        setJointDependence is.
        """
        self._assert_rate_dependent('Compatibility')
        if self.ljd_scaling is not None:
            raise ValueError(f"{self.name} already declares a joint dependence, so it cannot "
                             f"also declare server pools, which are a rate law of their own.")
        if server_type.get_num_of_servers() < 1:
            raise ValueError(f"Server pool '{server_type.name}' must hold at least one server")
        compat = server_type.get_compatible_classes()
        if not compat:
            raise ValueError(f"Server pool '{server_type.name}' is compatible with no operand, "
                             f"so it can never serve")
        pools = getattr(self, '_server_pools', None)
        if pools is None:
            pools = []
            self._server_pools = pools
        for p in pools:
            if p['name'] == server_type.name:
                raise ValueError(f"Server pool '{server_type.name}' is already declared on "
                                 f"{self.name}")
        # Resolved to names here, as addConstraint does: getStruct reads strings.
        compat_names = [c if isinstance(c, str) else c.name for c in compat]
        if len(set(compat_names)) != len(compat_names):
            raise ValueError(f"Server pool '{server_type.name}' names the same operand more "
                             f"than once")
        server_type.set_id(len(pools))
        pools.append({'name': server_type.name,
                      'count': float(server_type.get_num_of_servers()),
                      'rate': float(server_type.get_rate()),
                      'compatible': compat_names})
        return self

    def getServerTypes(self):
        """Declared compatibility pools."""
        return self.server_pools

    def hasServerPools(self):
        """Whether this element declares compatibility pools."""
        return bool(self.server_pools)

    add_server_type = addServerType
    get_server_types = getServerTypes
    has_server_pools = hasServerPools

    def setLoadDependence(self, alpha):
        """
        alpha[n] is the service-rate scaling of the station that represents this
        server in its layer when that station holds n jobs in total, as in
        Queue.setLoadDependence. The scaling multiplies the station rate on top
        of its multiplicity, so a multi-server host applies min(n,m)*alpha[n].
        """
        self._assert_rate_dependent('Load')
        alpha = np.atleast_1d(np.asarray(alpha, dtype=float)).ravel()
        if alpha.size == 0 or not np.all(np.isfinite(alpha)) or np.any(alpha <= 0):
            raise ValueError("Load-dependence scalings must be finite and positive.")
        self._lld_scaling = alpha
        return self

    def setClassDependence(self, beta, peak_rate_per_operand=None):
        """
        beta(n) takes the per-operand population vector of this server: n[j]
        counts the jobs held on behalf of operand j, which is task j of a
        Processor or entry j of a Task, in the same order as the columns of
        setConstraint. It returns a scalar shared by every operand, or a
        per-operand vector. peak_rate_per_operand is REQUIRED (scalar or
        per-operand vector) and normalizes Util = T*S/peak. Product form holds
        only where an operand occupies the layer station through a single job
        class; otherwise SolverLN emits the equivalent joint dependence, which is
        numerically identical but carries no exactness guarantee.
        """
        self._assert_rate_dependent('Class')
        _assert_dependence_handle(beta, peak_rate_per_operand, 'Class')
        self._lcd_scaling = beta
        self._lcd_scaling_peak = np.atleast_1d(np.asarray(peak_rate_per_operand, dtype=float)).ravel()
        return self

    def setJointDependence(self, eta, peak_rate_per_operand=None):
        """
        eta(n) reads the per-operand population vector of this server arbitrarily
        (e.g. min(n[0],c)) and is therefore non-product-form: solvers treat it as
        an approximation. Operand order and the required peak_rate_per_operand
        are as in setClassDependence.
        """
        self._assert_rate_dependent('Joint')
        _assert_dependence_handle(eta, peak_rate_per_operand, 'Joint')
        if self.server_pools:
            raise ValueError(f"{self.name} already declares server pools, which are themselves a "
                             f"rate law, so it cannot also take a joint dependence.")
        self._ljd_scaling = eta
        self._ljd_scaling_peak = np.atleast_1d(np.asarray(peak_rate_per_operand, dtype=float)).ravel()
        return self

    def hasRateDependence(self):
        """Whether this element declares any service-rate dependence."""
        return (self.lld_scaling is not None or self.lcd_scaling is not None
                or self.ljd_scaling is not None or bool(self.server_pools))

    def _assert_rate_dependent(self, what):
        """Reject servers whose layer station does not admit a rate scaling."""
        sched = getattr(self, 'sched_strategy', None)
        if sched not in (SchedStrategy.PS, SchedStrategy.FCFS):
            raise ValueError(f"{what}-dependence supported only for processor sharing (PS) and "
                             f"first-come first-serve (FCFS) servers, but {self.name} is scheduled {sched}.")

    set_load_dependence = setLoadDependence
    set_class_dependence = setClassDependence
    set_joint_dependence = setJointDependence
    has_rate_dependence = hasRateDependence


def _assert_dependence_handle(f, peak, what):
    """Common validation of a class- or joint-dependence declaration."""
    if not callable(f):
        raise ValueError(f"{what} dependence must be specified through a function handle.")
    if peak is None:
        raise ValueError(f"{what} dependence requires an explicit peak rate: pass a scalar "
                         f"(identical peak for every operand) or a per-operand vector.")
    peak = np.atleast_1d(np.asarray(peak, dtype=float)).ravel()
    if peak.size == 0 or not np.all(np.isfinite(peak)) or np.any(peak <= 0):
        raise ValueError("peak_rate_per_operand must be a finite positive scalar or per-operand vector.")


def _expand_peak(peak, ncols, elemname, colwhat, what):
    """Per-operand peak rate scaling of a class- or joint-dependence declaration."""
    peak = np.atleast_1d(np.asarray(peak, dtype=float)).ravel()
    if peak.size == 1:
        return np.full(ncols, peak[0])
    if peak.size != ncols:
        raise ValueError(f"{what}-dependence peak rate on {elemname} has {peak.size} entries "
                         f"but there are {ncols} {colwhat}.")
    return peak


def _operand_names(operands):
    """Operand names from objects, names, or a sequence mixing the two."""
    if operands is None:
        raise ValueError("Admission constraint requires at least one operand.")
    if isinstance(operands, str):
        return [operands]
    if not isinstance(operands, (list, tuple, np.ndarray)):
        operands = [operands]
    if len(operands) == 0:
        raise ValueError("Admission constraint requires at least one operand.")
    names = []
    for k, op in enumerate(operands):
        if isinstance(op, str):
            names.append(op)
        elif hasattr(op, 'name'):
            names.append(op.name)
        else:
            raise ValueError(f"Admission constraint operand {k + 1} is neither a name nor a LayeredNetwork element.")
    return names


class Task(AdmissionConstrained, RateDependent):
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
            # MATLAB sets multiplicity to Inf when scheduling is INF with finite multiplicity
            # This is necessary for correct njobs computation in SolverLN
            if self.sched_strategy == SchedStrategy.INF and np.isfinite(self.multiplicity):
                self.multiplicity = np.inf
            # Register with model
            if hasattr(self._model, 'add_task'):
                self._model.add_task(self)
        else:
            # 3-arg call: (name, multiplicity, sched_strategy)
            self._model = None
            self.name = model_or_name
            self.multiplicity = name_or_mult
            self.sched_strategy = self._convert_sched_strategy(mult_or_sched)
            # MATLAB sets multiplicity to Inf when scheduling is INF with finite multiplicity
            if self.sched_strategy == SchedStrategy.INF and np.isfinite(self.multiplicity):
                self.multiplicity = np.inf

        self.processor = None
        self.think_time = None
        self.setup_time = None
        self.delay_off_time = None
        # Scheduling priority; lower is served first, the convention sn.classprio
        # uses in a Network. Read by the HOL disciplines of a host or a task.
        self.priority = 0
        self.entries = []
        self.activities = []
        self.precedences = []
        self._fan_in = {}   # {source_task_name: value}
        self._fan_out = {}  # {dest_task_name: value}

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

    def set_priority(self, priority: int) -> 'Task':
        """Set the scheduling priority of this task (lower is served first)."""
        self.priority = int(priority)
        return self

    def get_priority(self) -> int:
        """Scheduling priority of this task."""
        return self.priority

    setPriority = set_priority
    getPriority = get_priority

    def on(self, processor: 'Processor') -> 'Task':
        """Deploy this task on a processor."""
        self.processor = processor
        processor.tasks.append(self)
        return self

    def set_think_time(self, think_time: Union[Distribution, float]) -> 'Task':
        """Set the think time for this task."""
        if isinstance(think_time, (int, float)):
            self.think_time = _dist_from_mean_scv(float(think_time), 1.0)
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

    def set_replication(self, replication: int) -> 'Task':
        """Set replication level (snake_case alias for setReplication)."""
        return self.setReplication(replication)

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
        return _get_dist_mean(self.think_time)

    def getThinkTimeSCV(self) -> float:
        """Get think time SCV."""
        if self.think_time is None:
            return 1.0
        return _get_dist_scv(self.think_time)

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
        return _get_dist_mean(self.setup_time)

    def getDelayOffTimeMean(self) -> float:
        """Get delay-off time mean."""
        if self.delay_off_time is None:
            return 0.0
        return _get_dist_mean(self.delay_off_time)

    def setThinkTime(self, think_time) -> 'Task':
        """Set think time."""
        return self.set_think_time(think_time)

    def set_fan_in(self, source: str, value: int) -> 'Task':
        """Set fan-in from a source task (snake_case alias for setFanIn)."""
        return self.setFanIn(source, value)

    def setFanIn(self, source: str, value: int) -> 'Task':
        """Set fan-in from a source task (for replication load distribution)."""
        self._fan_in[source] = value
        return self

    def getFanIn(self) -> dict:
        """Get fan-in mapping {source_task_name: value}."""
        return self._fan_in

    def set_fan_out(self, dest: str, value: int) -> 'Task':
        """Set fan-out to a destination task (snake_case alias for setFanOut)."""
        return self.setFanOut(dest, value)

    def setFanOut(self, dest: str, value: int) -> 'Task':
        """Set fan-out to a destination task (for replication load distribution)."""
        self._fan_out[dest] = value
        return self

    def getFanOut(self) -> dict:
        """Get fan-out mapping {dest_task_name: value}."""
        return self._fan_out

    # Snake_case aliases
    get_multiplicity = getMultiplicity
    get_replication = getReplication
    set_replication = setReplication
    get_scheduling = getScheduling
    get_think_time_mean = getThinkTimeMean
    get_think_time_scv = getThinkTimeSCV
    get_parent = getParent
    get_precedences = getPrecedences
    # camelCase alias, under the name MATLAB and the JAR use: an LQN script
    # reaches for addPrecedence right after building the activity graph.
    addPrecedence = add_precedence
    get_setup_time_mean = getSetupTimeMean
    get_delay_off_time_mean = getDelayOffTimeMean
    set_fan_in = setFanIn
    get_fan_in = getFanIn
    set_fan_out = setFanOut
    get_fan_out = getFanOut

    def has_setup_delayoff(self) -> bool:
        """Return False for regular Task. SetupTask overrides this."""
        return False


class Processor(AdmissionConstrained, RateDependent):
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

    def set_replication(self, replication: int) -> 'Processor':
        """Set replication level (snake_case alias for setReplication)."""
        return self.setReplication(replication)

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
        super().__init__(model, name, multiplicity, SchedStrategy.FCFS)
        # CacheTask specific fields
        self.total_items = total_items
        self.cache_capacity = cache_capacity
        self.replacement_strategy = replacement_strategy
        self.retrieval = False

    def __hash__(self):
        return id(self)

    def __eq__(self, other):
        return self is other

    def on(self, processor: 'Processor') -> 'CacheTask':
        """Deploy this cache task on a processor."""
        self.processor = processor
        processor.tasks.append(self)
        return self

    def set_retrieval(self, retrieval: bool = True) -> 'CacheTask':
        """Enable/disable a delayed-hit retrieval system on the cache miss path.

        When set, concurrent misses for the same item arriving while a fetch (the
        miss-branch activity and its backend calls) is in flight are parked and
        released together as delayed hits when the fetch completes, instead of each
        triggering an independent fetch. Mirrors MATLAB CacheTask.setRetrieval.
        """
        self.retrieval = bool(retrieval)
        return self

    def has_retrieval(self) -> bool:
        return bool(getattr(self, 'retrieval', False))


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
        self._forwarding_dests = []
        self._forwarding_probs = []
        self._arrival = None
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
    repl: np.ndarray = None          # Replication factors
    prio: np.ndarray = None          # Task scheduling priority (0 elsewhere); lower is served first
    maxmult: np.ndarray = None       # Max multiplicities (for servers)
    graph: np.ndarray = None         # Adjacency graph
    iscaller: np.ndarray = None      # Caller matrix
    issynccaller: np.ndarray = None  # Sync caller matrix
    isasynccaller: np.ndarray = None # Async caller matrix
    isref: np.ndarray = None         # Reference task flags
    iscache: np.ndarray = None       # Cache task flags (True if task is a CacheTask)
    hasretrieval: np.ndarray = None  # 1 on a CacheTask with a delayed-hit retrieval miss path
    nitems: np.ndarray = None        # Number of items for CacheTask/ItemEntry
    itemcap: Dict[int, np.ndarray] = field(default_factory=dict)  # Cache capacity per level
    replacestrat: np.ndarray = None  # Cache replacement strategy
    itemproc: Dict[int, object] = field(default_factory=dict)     # Item popularity distribution
    callpair: np.ndarray = None      # Call pairs (caller_act, callee_entry, mean_calls)
    # Calls dispatched as a group by a routing strategy, as
    # (caller_activity_idx, RoutingStrategy, [target_entry_idx, ...]).
    # Requires the squashed layering; see SolverLN._assert_call_groups.
    callgroups: List = field(default_factory=list)
    parent: np.ndarray = None        # Parent relationships
    replygraph: np.ndarray = None    # Reply graph (nacts x nentries)
    actphase: np.ndarray = None      # Activity phase (1..3) for each activity
    actposttype: np.ndarray = None   # Activity post-precedence type (POST_AND, POST_OR, POST_SEQ)
    actpretype: np.ndarray = None    # Activity pre-precedence type (PRE_AND, PRE_OR, PRE_SEQ)
    # Quorum count of an AND-join, indexed by the join target activity. Equals the number
    # of predecessors when the join waits for all of them.
    actquorum: np.ndarray = None

    # Names
    # lincon[i] is the (A, b) pair of the admission constraint A*n <= b on the layer
    # station of host or task i. Columns are that host's tasks, or that task's
    # entries, in tasksof/entriesof order. Absent where unconstrained -- see
    # _kb/04-networkstruct.md
    lincon: dict = None

    # Service-rate dependences on the layer station of host or task i, keyed by absolute index i.
    # lldscaling[i] is the vector alpha(n) applied at total population n; cdscaling and jdscaling
    # are handles over that server's operands, in the tasksof/entriesof order used by lincon
    # columns, with cdscalingpeak/jdscalingpeak their per-operand peak rate scaling. Absent where
    # undeclared -- see _kb/04-networkstruct.md
    lldscaling: dict = None
    cdscaling: dict = None
    cdscalingpeak: dict = None
    jdscaling: dict = None
    jdscalingpeak: dict = None
    # pools[i] is dict(names, counts, rates, compat) of the compatibility pools declared on
    # that server; compat[t, j] nonzero = pool t may serve operand j. Absent where not declared
    pools: dict = None

    names: np.ndarray = None
    hashnames: np.ndarray = None

    # Mappings
    tasksof: Dict[int, List[int]] = field(default_factory=dict)
    entriesof: Dict[int, List[int]] = field(default_factory=dict)
    actsof: Dict[int, List[int]] = field(default_factory=dict)
    callsof: Dict[int, List[int]] = field(default_factory=dict)

    # Service demands and think times
    hostdem: Dict[int, float] = field(default_factory=dict)
    # think and actthink hold the Distribution when one is set (0.0 otherwise),
    # as in MATLAB lsn.think{} / lsn.actthink{}, so the SCV survives.
    think: Dict[int, object] = field(default_factory=dict)
    think_scv: Dict[int, float] = field(default_factory=dict)
    actthink: Dict[int, object] = field(default_factory=dict)

    # Full host-demand Distribution objects keyed by absolute index (mirrors
    # MATLAB lqn.hostdem{} / hostdem_proc). lqn.hostdem stays a scalar mean for
    # the numeric consumers; hostdem_proc preserves SCV/phase-type so SolverLN
    # can set the real service distribution on the layer models (moment2/moment3).
    hostdem_proc: Dict[int, object] = field(default_factory=dict)
    hostdem_scv: Dict[int, float] = field(default_factory=dict)

    # Entry-level open arrival distributions (absolute entry index -> Distribution)
    arrival: Dict[int, object] = field(default_factory=dict)

    # Scheduling strategies
    sched: Dict[int, str] = field(default_factory=dict)

    # Fan-out matrix: fanout[source_task_idx, dest_task_idx] = fan-out value
    fanout: np.ndarray = None


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


    def findSolver(self, metric: str = '', showAll: bool = False):
        """Which solvers and solver methods can analyze THIS model.

            model.findSolver()                # every (solver, method) pair that runs
            model.findSolver('cdf')           # ... that returns a passage-time law
            model.findSolver('getCdfRespT')   # the same question, asked by accessor
            model.findSolver('', True)        # also the pairs that are refused, and why

        The returned DataFrame has one row per pair, with columns Solver,
        Method, Runnable, Class ('exact', 'approx', 'bound' or 'simulation'),
        Metrics and Reason. Method is the method name to pass as a solver method, so
        a row can be acted on directly::

            T = model.findSolver('cdf')
            solver = LINE(model, T.Method[0])

        findMethod and help are aliases of this method.

        Args:
            metric: measure group ('cdf') or accessor ('getCdfRespT') to narrow
                the report to; '' or 'any' keeps every pair.
            showAll: also list the refused pairs, with the reason each was
                refused.

        Returns:
            pandas.DataFrame with the six columns above.
        """
        # The gate lives in SolverAUTO, which is the class that already knows
        # every family, how to build one and what each refuses. Asking it here
        # rather than reimplementing the walk is what keeps the model's answer
        # and AUTO's own dispatch from being two opinions.
        from .solvers.solver_auto.solver_auto import SolverAUTO
        # silenced() wraps the CONSTRUCTION too: it probes every candidate
        # with supports(model), which warns on a model one of them refuses.
        with SolverAUTO.silenced():
            auto = SolverAUTO(self, verbose=False)
        return auto.findSolver(metric, showAll)

    def findMethod(self, metric: str = '', showAll: bool = False):
        """Alias of findSolver: which solvers and solver methods can analyze
        this model.

        The two names exist because the question is asked both ways round --
        "which solver do I use" and "which method do I pass" -- and the answer
        is the same table, whose Method column carries the method name either caller
        needs.
        """
        return self.findSolver(metric, showAll)

    def help(self, metric: str = '', showAll: bool = False):
        """Alias of findSolver: what can this model be solved with?"""
        return self.findSolver(metric, showAll)

    find_solver = findSolver
    find_method = findMethod

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
                demand = _dist_from_mean_scv(float(host_demand), 1.0)
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
        # Elements are 0-based: local h,t,e,a run 0..n-1 and the absolute index is
        # shift+local, so 0..nidx-1. MATLAB keeps the 1-based numbering.
        idx = 0

        # Processors (hosts)
        for proc in self.processors:
            self._proc_idx[proc] = idx
            idx += 1

        tshift = idx

        # Tasks
        for task in self.tasks:
            self._task_idx[task] = idx
            idx += 1

        eshift = idx

        # Entries
        for entry in self.entries:
            self._entry_idx[entry] = idx
            idx += 1

        ashift = idx

        # Activities
        for act in self.activities:
            self._act_idx[act] = idx
            idx += 1

        return tshift, eshift, ashift, idx

    def _lsn_max_multiplicity(self, lqn: 'LayeredNetworkStruct', nidx: int) -> np.ndarray:
        """
        Compute maximum multiplicity (throughput capacity) for each node.

        This implements the MATLAB lsn_max_multiplicity algorithm using
        Kahn's algorithm for topological sorting to propagate flow constraints
        through the network.

        Args:
            lqn: The LayeredNetworkStruct being built
            nidx: Total number of indices

        Returns:
            Array of maxmult values indexed by node index
        """
        # Use DAG for flow propagation (needs to be built or use graph)
        # The DAG excludes loop-back edges
        n = nidx

        # Get the adjacency graph (use dag if available, else graph)
        if hasattr(lqn, 'dag') and lqn.dag is not None:
            ag = (lqn.dag > 0).astype(int)
        elif hasattr(lqn, 'graph') and lqn.graph is not None:
            ag = (lqn.graph > 0).astype(int)
        else:
            # If no graph yet, return mult as fallback
            maxmult = np.zeros((1, n))
            for proc in self.processors:
                idx = self._proc_idx[proc]
                mult_val = lqn.mult[0, idx] if idx < lqn.mult.shape[1] else 1.0
                if np.isinf(mult_val):
                    maxmult[0, idx] = 0
                else:
                    maxmult[0, idx] = mult_val
            for task in self.tasks:
                idx = self._task_idx[task]
                maxmult[0, idx] = lqn.mult[0, idx] if idx < lqn.mult.shape[1] else 1.0
            return maxmult

        # Ensure ag is the right size
        if ag.shape[0] < n:
            new_ag = np.zeros((n, n))
            new_ag[:ag.shape[0], :ag.shape[1]] = ag
            ag = new_ag

        # Get mult and type arrays
        # MATLAB (lsn_max_multiplicity.m lines 59-61):
        # if length(mult) < n
        #     mult(end+1:n) = Inf;
        # end
        # This pads mult with Inf for entries and activities, allowing flow to pass through
        mult = lqn.mult.flatten() if lqn.mult is not None else np.ones(n)
        if len(mult) < n:
            new_mult = np.ones(n) * np.inf
            new_mult[:len(mult)] = mult
            mult = new_mult
        # Replace zeros with Inf for entries and activities (they don't limit flow)
        # MATLAB only defines mult for hosts+tasks; entries/activities get Inf by padding
        for i in range(n):
            if mult[i] == 0:
                mult[i] = np.inf

        # Build isref array (need to set it before this function is called)
        isref = np.zeros(n)
        for task in self.tasks:
            if task.sched_strategy == SchedStrategy.REF:
                isref[self._task_idx[task]] = 1

        # Build type array
        type_arr = np.zeros(n)
        for proc in self.processors:
            type_arr[self._proc_idx[proc]] = LayeredNetworkElement.PROCESSOR
        for task in self.tasks:
            type_arr[self._task_idx[task]] = LayeredNetworkElement.TASK
        for entry in self.entries:
            type_arr[self._entry_idx[entry]] = LayeredNetworkElement.ENTRY
        for act in self.activities:
            type_arr[self._act_idx[act]] = LayeredNetworkElement.ACTIVITY

        # Topological sort using Kahn's algorithm
        order = self._kahn_topological_sort(ag)

        # Initially load ref task multiplicity into inflow
        inflow = np.zeros(n)
        for i in range(n):
            if type_arr[i] == LayeredNetworkElement.TASK and isref[i]:
                inflow[i] = mult[i]
            # Also account for entries with open arrivals
            elif type_arr[i] == LayeredNetworkElement.ENTRY:
                if hasattr(lqn, 'arrival') and lqn.arrival is not None:
                    if isinstance(lqn.arrival, dict) and i in lqn.arrival:
                        inflow[i] = 1
                    elif isinstance(lqn.arrival, list) and i < len(lqn.arrival) and lqn.arrival[i] is not None:
                        inflow[i] = 1

        outflow = np.zeros(n)

        # see _kb/04-networkstruct.md (Python native layered.py port notes) for rationale
        has_setup = np.zeros(n)
        if getattr(lqn, 'hassetup', None) is not None:
            isf = np.asarray(lqn.hassetup).ravel()
            has_setup[:min(n, isf.size)] = isf[:min(n, isf.size)]

        # Propagate flow through DAG in topological order
        for k in range(len(order)):
            i = order[k]
            if has_setup[i] and inflow[i] > 0:
                outflow[i] = mult[i]
            else:
                outflow[i] = min(inflow[i], mult[i])
            for j in range(n):
                if j != i and ag[i, j]:
                    inflow[j] = inflow[j] + outflow[i]

        # For non-ref tasks with INF mult, set outflow = Inf (MATLAB lines 73-77)
        for i in range(n):
            if type_arr[i] == LayeredNetworkElement.TASK and np.isinf(mult[i]) and not isref[i]:
                outflow[i] = np.inf

        # see _kb/04-networkstruct.md (Python native layered.py port notes) for rationale
        for i in range(n):
            if type_arr[i] == LayeredNetworkElement.TASK and not isref[i]:
                if outflow[i] == 0 and not np.isinf(mult[i]) and mult[i] > 0:
                    outflow[i] = mult[i]

        # see _kb/04-networkstruct.md (Python native layered.py port notes) for rationale
        for proc in self.processors:
            proc_idx = self._proc_idx[proc]
            if outflow[proc_idx] == 0 and not np.isinf(mult[proc_idx]) and mult[proc_idx] > 0:
                outflow[proc_idx] = mult[proc_idx]

        return outflow.reshape(1, -1)

    def _kahn_topological_sort(self, ag: np.ndarray) -> List[int]:
        """
        Perform topological sort using Kahn's algorithm.

        Args:
            ag: Adjacency matrix (ag[i,j] = 1 means edge from i to j)

        Returns:
            List of node indices in topological order
        """
        n = ag.shape[0]
        in_degree = np.sum(ag, axis=0).astype(int)  # Count incoming edges for each node
        queue = [i for i in range(n) if in_degree[i] == 0]
        order = []

        while queue:
            node = queue.pop(0)
            order.append(node)
            for j in range(n):
                if ag[node, j]:
                    in_degree[j] -= 1
                    if in_degree[j] == 0:
                        queue.append(j)

        # If not all nodes are in order, there's a cycle - return partial order
        if len(order) < n:
            # Add remaining nodes
            for i in range(n):
                if i not in order:
                    order.append(i)

        return order

    def _get_prec_type(self, prec) -> PrecedenceType:
        """
        Get the PrecedenceType from either layered.py or workflow.py ActivityPrecedence.

        This handles the class name collision where both modules define ActivityPrecedence
        with different attribute names (prec_type vs pre_type/post_type).
        """
        # If it has prec_type attribute, use it directly (layered.py ActivityPrecedence)
        if hasattr(prec, 'prec_type') and prec.prec_type is not None:
            return prec.prec_type

        # workflow.py ActivityPrecedence uses pre_type and post_type
        # Map post_type to our PrecedenceType enum
        if hasattr(prec, 'post_type') and prec.post_type is not None:
            post_type = prec.post_type
            # Handle both string and enum values
            post_type_str = post_type if isinstance(post_type, str) else str(post_type)
            post_type_upper = post_type_str.upper().replace('-', '_')

            if 'POST_AND' in post_type_upper or 'AND' in post_type_upper:
                return PrecedenceType.PARALLEL
            elif 'POST_OR' in post_type_upper or '_OR' in post_type_upper:
                return PrecedenceType.CHOICE
            elif 'POST_LOOP' in post_type_upper or 'LOOP' in post_type_upper:
                return PrecedenceType.LOOP
            else:
                return PrecedenceType.SERIAL

        # Default to SERIAL
        return PrecedenceType.SERIAL

    def _get_prec_activities(self, prec):
        """
        Get the activities list from either layered.py or workflow.py ActivityPrecedence.

        For workflow.py ActivityPrecedence, this returns post_acts which contains
        [loop_body_activities..., end_activity] for loops, or just post_activities
        for other types.
        """
        # layered.py uses 'activities' for serial and loop body
        if hasattr(prec, 'activities') and prec.activities:
            return prec.activities
        # workflow.py uses 'post_acts' for the activities list
        if hasattr(prec, 'post_acts') and prec.post_acts:
            return prec.post_acts
        # Also try post_activities as fallback
        if hasattr(prec, 'post_activities') and prec.post_activities:
            return prec.post_activities
        return []

    def _get_prec_pre_activities(self, prec):
        """Get pre_activities from either ActivityPrecedence type."""
        if hasattr(prec, 'pre_activities') and prec.pre_activities:
            return prec.pre_activities
        # workflow.py uses 'pre_acts'
        if hasattr(prec, 'pre_acts') and prec.pre_acts:
            return prec.pre_acts
        return []

    def _get_prec_post_activities(self, prec):
        """Get post_activities from either ActivityPrecedence type."""
        if hasattr(prec, 'post_activities') and prec.post_activities:
            return prec.post_activities
        # workflow.py uses 'post_acts'
        if hasattr(prec, 'post_acts') and prec.post_acts:
            return prec.post_acts
        return []

    def _get_prec_count(self, prec) -> float:
        """Get loop count from either ActivityPrecedence type."""
        if hasattr(prec, 'count') and prec.count is not None:
            return float(prec.count)
        if hasattr(prec, 'post_params') and prec.post_params is not None:
            params = prec.post_params
            if isinstance(params, np.ndarray) and params.size > 0:
                return float(params.flat[0])
            elif isinstance(params, (list, tuple)) and len(params) > 0:
                return float(params[0])
        return 1.0

    def _get_prec_probabilities(self, prec):
        """Get probabilities from either ActivityPrecedence type."""
        if hasattr(prec, 'probabilities') and prec.probabilities:
            return prec.probabilities
        if hasattr(prec, 'post_params') and prec.post_params is not None:
            params = prec.post_params
            if isinstance(params, np.ndarray):
                return params.flatten().tolist()
            elif isinstance(params, (list, tuple)):
                return list(params)
        return []

    def _get_act_idx(self, act):
        """
        Get the index for an activity from _act_idx.

        Handles both Activity objects and string names.
        """
        # Direct lookup if it's an Activity object
        if act in self._act_idx:
            return self._act_idx[act]

        # If it's a string, look up by name
        if isinstance(act, str):
            for activity, idx in self._act_idx.items():
                if hasattr(activity, 'name') and activity.name == act:
                    return idx

        return None

    def _act_in_idx(self, act) -> bool:
        """Check if an activity (object or name) is in _act_idx."""
        return self._get_act_idx(act) is not None

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
        lqn.names = np.empty(nidx, dtype=object)
        lqn.hashnames = np.empty(nidx, dtype=object)

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

        # Build type array (matches MATLAB LayeredNetworkElement enum)
        # 0=PROCESSOR, 1=TASK, 2=ENTRY, 3=ACTIVITY
        lqn.type = np.zeros(nidx, dtype=int)
        for proc in self.processors:
            lqn.type[self._proc_idx[proc]] = 0  # PROCESSOR
        for task in self.tasks:
            lqn.type[self._task_idx[task]] = 1  # TASK
        for entry in self.entries:
            lqn.type[self._entry_idx[entry]] = 2  # ENTRY
        for act in self.activities:
            lqn.type[self._act_idx[act]] = 3  # ACTIVITY

        # Build multiplicity matrix
        lqn.mult = np.zeros((1, nidx + 1))
        for proc in self.processors:
            lqn.mult[0, self._proc_idx[proc]] = proc.multiplicity
        for task in self.tasks:
            mult = task.multiplicity
            # Keep float('inf') as np.inf - solver_ln.py handles infinite servers
            if mult == float('inf'):
                mult = np.inf
            lqn.mult[0, self._task_idx[task]] = mult

        # Build replication factors (defaults to 1)
        # MATLAB: lsn.repl is indexed by host/task index (1:tshift+ntasks)
        # We need to cover all indices up to nidx for consistent indexing
        lqn.repl = np.ones((1, nidx + 1))

        # Populate replication factors from processor/task attributes
        for proc in self.processors:
            if proc in self._proc_idx:
                repl = proc.getReplication()
                if repl > 1:
                    lqn.repl[0, self._proc_idx[proc]] = repl
        for task in self.tasks:
            if task in self._task_idx:
                repl = task.getReplication()
                if repl > 1:
                    lqn.repl[0, self._task_idx[task]] = repl

        # Task scheduling priority, read by the HOL disciplines; 0 elsewhere.
        lqn.prio = np.zeros((1, nidx + 1))
        for task in self.tasks:
            if task in self._task_idx:
                lqn.prio[0, self._task_idx[task]] = task.get_priority()

        # Note: maxmult is computed later after the graph is built

        # Build reference task flags
        lqn.isref = np.zeros((nidx, 1))
        for task in self.tasks:
            if task.sched_strategy == SchedStrategy.REF:
                lqn.isref[self._task_idx[task], 0] = 1

        # Build cache task flags (matches MATLAB: lsn.iscache = lsn.nitems > 0)
        lqn.iscache = np.zeros((nidx, 1))
        # Delayed-hit retrieval flag: 1 on a CacheTask whose miss path is a
        # retrieval system (CacheTask.set_retrieval). Mirrors JAR lsn.hasretrieval.
        lqn.hasretrieval = np.zeros((nidx, 1))
        for task in self.tasks:
            if isinstance(task, CacheTask):
                lqn.iscache[self._task_idx[task], 0] = 1
                if task.has_retrieval():
                    lqn.hasretrieval[self._task_idx[task], 0] = 1

        # Build cache-related arrays (matches MATLAB getStruct.m lines 66-68, 135-137, 183-184)
        lqn.nitems = np.zeros((nidx, 1))
        lqn.itemcap = {}
        lqn.replacestrat = np.zeros((nidx, 1), dtype=int)
        lqn.itemproc = {}

        for task in self.tasks:
            idx = self._task_idx[task]
            if isinstance(task, CacheTask):
                lqn.nitems[idx, 0] = task.total_items
                # Convert cache_capacity to array if needed
                cap = task.cache_capacity
                if isinstance(cap, (list, tuple)):
                    lqn.itemcap[idx] = np.array(cap, dtype=int)
                else:
                    lqn.itemcap[idx] = np.array([int(cap)])
                # Store replacement strategy as integer
                lqn.replacestrat[idx, 0] = task.replacement_strategy.value if hasattr(task.replacement_strategy, 'value') else int(task.replacement_strategy)

        for entry in self.entries:
            idx = self._entry_idx[entry]
            if isinstance(entry, ItemEntry):
                lqn.nitems[idx, 0] = entry.total_items
                lqn.itemproc[idx] = entry.access_prob

        # Build caller matrices
        lqn.iscaller = np.zeros((nidx, nidx))
        lqn.issynccaller = np.zeros((nidx, nidx))
        lqn.isasynccaller = np.zeros((nidx, nidx))

        # Build call pairs and caller relationships
        calls = []
        # Call groups dispatched by a routing strategy, as
        # (caller_activity_idx, RoutingStrategy, [target_entry_idx, ...]).
        # Empty for every model that does not use Activity.synch_call_rrobin.
        lqn.callgroups = []
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

                # Mark caller relationships (matches MATLAB getStruct.m lines 294-297)
                # MATLAB sets: iscaller(tidx, target_tidx), iscaller(aidx, target_tidx),
                #              iscaller(tidx, target_eidx), iscaller(aidx, target_eidx)
                lqn.iscaller[caller_task_idx, target_task_idx] = 1
                lqn.iscaller[act_idx, target_task_idx] = 1
                lqn.iscaller[caller_task_idx, target_entry_idx] = 1
                lqn.iscaller[act_idx, target_entry_idx] = 1

                if call_type == CallType.SYNC:
                    lqn.issynccaller[caller_task_idx, target_task_idx] = 1
                    lqn.issynccaller[act_idx, target_task_idx] = 1
                    lqn.issynccaller[caller_task_idx, target_entry_idx] = 1
                    lqn.issynccaller[act_idx, target_entry_idx] = 1
                else:
                    lqn.isasynccaller[caller_task_idx, target_task_idx] = 1
                    lqn.isasynccaller[act_idx, target_task_idx] = 1
                    lqn.isasynccaller[caller_task_idx, target_entry_idx] = 1
                    lqn.isasynccaller[act_idx, target_entry_idx] = 1

                calls.append((act_idx, target_entry_idx, mean_calls, call_type))

            for strategy, group_entries in getattr(act, 'call_groups', []):
                gidx = [self._entry_idx[e] for e in group_entries
                        if e.task is not None and e in self._entry_idx]
                if len(gidx) >= 2:
                    lqn.callgroups.append((act_idx, strategy, gidx))

        # see _kb/04-networkstruct.md (Python native layered.py port notes) for rationale

        # Collect forwarding calls from entries
        fwd_calls = []
        for entry in self.entries:
            if entry._forwarding_dests:
                source_eidx = self._entry_idx[entry]
                for i, dest in enumerate(entry._forwarding_dests):
                    prob = entry._forwarding_probs[i]
                    if isinstance(dest, str):
                        # Find target entry by name
                        target_entry = None
                        for e in self.entries:
                            if e.name == dest:
                                target_entry = e
                                break
                        if target_entry is None:
                            continue
                    else:
                        target_entry = dest
                    target_eidx = self._entry_idx[target_entry]
                    fwd_calls.append((source_eidx, target_eidx, prob))

        lqn.ncalls = len(calls) + len(fwd_calls)
        if lqn.ncalls > 0:
            lqn.callpair = np.zeros((lqn.ncalls, 3))
            lqn.calltype = np.zeros(lqn.ncalls, dtype=int)
            lqn.callproc = [None] * lqn.ncalls
            # Add sync/async calls
            for i, (act_idx, entry_idx, mean_calls, call_type) in enumerate(calls):
                lqn.callpair[i, 0] = act_idx
                lqn.callpair[i, 1] = entry_idx
                lqn.callpair[i, 2] = mean_calls
                # calltype: 1=SYNC, 2=ASYNC (matches MATLAB CallType enum)
                lqn.calltype[i] = 1 if call_type == CallType.SYNC else 2
                # callproc: distribution for call multiplicity
                lqn.callproc[i] = _call_count_dist(mean_calls)
            # Add forwarding calls (calltype=3=FWD, matching JAR CallType.FWD)
            # Forwarding callpair: col0=source_entry, col1=target_entry, col2=fwd_prob
            # Note: NOT added to issynccaller/isasynccaller (forwarding is not a call
            # dependency); the graph entry->entry edge is added after graph construction
            # below (matches MATLAB getStruct.m lsn.graph(eidx, target_eidx) = 1)
            fwd_start = len(calls)
            for j, (source_eidx, target_eidx, fwd_prob) in enumerate(fwd_calls):
                cidx = fwd_start + j
                lqn.callpair[cidx, 0] = source_eidx
                lqn.callpair[cidx, 1] = target_eidx
                lqn.callpair[cidx, 2] = fwd_prob
                lqn.calltype[cidx] = 3  # CallType.FWD
                lqn.callproc[cidx] = _call_count_dist(fwd_prob)
        else:
            lqn.callpair = np.zeros((0, 3))
            lqn.calltype = np.zeros(0, dtype=int)
            lqn.callproc = []

        # Build callsof mapping (activity -> list of call indices)
        # This maps each source activity to the calls it makes
        lqn.callsof = {}
        for cidx in range(lqn.ncalls):
            if cidx < lqn.callpair.shape[0]:
                src_aidx = int(lqn.callpair[cidx, 0])
                if src_aidx >= 0:
                    if src_aidx not in lqn.callsof:
                        lqn.callsof[src_aidx] = []
                    lqn.callsof[src_aidx].append(cidx)

        # Build parent relationships
        lqn.parent = np.full((nidx, 1), -1.0)
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

        # Adjust task replication to account for host processor replication.
        # In LQN, task repl >= host repl. If task repl is 1 (default), inherit host repl.
        for task in self.tasks:
            task_idx = self._task_idx[task]
            host_idx = int(lqn.parent[task_idx, 0])
            if host_idx > 0:
                lqn.repl[0, task_idx] = max(lqn.repl[0, task_idx], lqn.repl[0, host_idx])

        # Build tasksof mapping (tasks on each host)
        lqn.tasksof = {}
        for hidx in range(lqn.nhosts):
            lqn.tasksof[hidx] = []
        for task in self.tasks:
            task_idx = self._task_idx[task]
            if task.processor:
                proc_idx = self._proc_idx[task.processor]
                lqn.tasksof[proc_idx].append(task_idx)

        # Build entriesof mapping (entries of each task)
        lqn.entriesof = {}
        for t in range(lqn.ntasks):
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

        # Admission constraint columns are only resolvable once tasksof/entriesof exist
        # -- see _kb/04-networkstruct.md
        lqn.lincon = {}
        lqn.lldscaling = {}
        lqn.cdscaling = {}
        lqn.cdscalingpeak = {}
        lqn.jdscaling = {}
        lqn.jdscalingpeak = {}
        lqn.pools = {}
        for cidx in range(lqn.nhosts + lqn.ntasks):
            if cidx < lqn.tshift:
                elem = self.processors[cidx - lqn.hshift]
                col_idx = lqn.tasksof.get(cidx, [])
                colwhat = 'tasks on this host'
            else:
                elem = self.tasks[cidx - lqn.tshift]
                col_idx = lqn.entriesof.get(cidx, [])
                colwhat = 'entries of this task'
            ncols = len(col_idx)
            # Rate dependences share the operand order of the constraint columns
            if elem.lld_scaling is not None:
                lqn.lldscaling[cidx] = elem.lld_scaling
            if elem.lcd_scaling is not None:
                lqn.cdscaling[cidx] = elem.lcd_scaling
                lqn.cdscalingpeak[cidx] = _expand_peak(elem.lcd_scaling_peak, ncols,
                                                       lqn.names[cidx], colwhat, 'Class')
            if elem.ljd_scaling is not None:
                lqn.jdscaling[cidx] = elem.ljd_scaling
                lqn.jdscalingpeak[cidx] = _expand_peak(elem.ljd_scaling_peak, ncols,
                                                       lqn.names[cidx], colwhat, 'Joint')
            # Compatibility pools name the operands they may serve, so the names
            # become columns only here, on the same operand order as the
            # constraints below.
            if elem.server_pools:
                col_names = [lqn.names[c] for c in col_idx]
                npools = len(elem.server_pools)
                compat = np.zeros((npools, ncols))
                counts = np.zeros(npools)
                rates = np.zeros(npools)
                pool_names = []
                for t, pool in enumerate(elem.server_pools):
                    pool_names.append(pool['name'])
                    counts[t] = pool['count']
                    rates[t] = pool['rate']
                    for cname in pool['compatible']:
                        if cname not in col_names:
                            raise ValueError(f"Server pool '{pool['name']}' on "
                                             f"{lqn.names[cidx]} names {cname}, which is not one "
                                             f"of the {colwhat}.")
                        compat[t, col_names.index(cname)] = 1
                # An operand no pool can serve would be served at rate zero and
                # never complete, so it is a declaration error, not an empty column.
                unserved = np.flatnonzero(~np.any(compat != 0, axis=0))
                if unserved.size:
                    raise ValueError(f"{col_names[int(unserved[0])]} on {lqn.names[cidx]} is "
                                     f"compatible with no server pool, so it can never be served.")
                # The pools describe HOW the declared servers are shared, not how
                # many there are, so the two statements have to agree. Letting
                # them diverge would leave the layer station sized by the
                # multiplicity and scaled by a peak taken over a different number
                # of servers, reporting a utilization against a denominator the
                # model never declared.
                mult_c = float(lqn.mult[0, cidx])
                if np.isfinite(mult_c) and float(np.sum(counts)) != mult_c:
                    raise ValueError(f"Server pools on {lqn.names[cidx]} hold "
                                     f"{float(np.sum(counts)):g} servers but its multiplicity is "
                                     f"{mult_c:g}; the pools partition the declared servers, so "
                                     f"the two must agree.")
                lqn.pools[cidx] = {'names': pool_names, 'counts': counts,
                                   'rates': rates, 'compat': compat}
            if not elem.hasLinearConstraints():
                continue
            A_pos, b_pos = elem.getLinearConstraints()
            if A_pos is not None and A_pos.shape[1] != ncols:
                raise ValueError(f"Admission constraint on {lqn.names[cidx]} has {A_pos.shape[1]} "
                                 f"columns but there are {ncols} {colwhat}.")
            blocks_a = [] if A_pos is None else [A_pos]
            blocks_b = [] if b_pos is None else [np.asarray(b_pos, dtype=float).ravel()]
            if elem.lincon_rows:
                # resolve rows declared by operand name against this server's columns
                col_names = [lqn.names[j] for j in col_idx]
                a_named = np.zeros((len(elem.lincon_rows), ncols))
                b_named = np.zeros(len(elem.lincon_rows))
                for r, (row_names, row_coeffs, row_cap) in enumerate(elem.lincon_rows):
                    for k, nm in enumerate(row_names):
                        if nm not in col_names:
                            raise ValueError(f"Admission constraint on {lqn.names[cidx]} names {nm}, "
                                             f"which is not one of the {colwhat}.")
                        a_named[r, col_names.index(nm)] = row_coeffs[k]
                    b_named[r] = row_cap
                blocks_a.append(a_named)
                blocks_b.append(b_named)
            lqn.lincon[cidx] = (np.vstack(blocks_a), np.concatenate(blocks_b))

        # Build entry-level open arrival distributions
        # Mirrors JAR LayeredNetwork.java:1158-1167 — entries that had setArrival(dist)
        # called record the Distribution keyed by entry absolute index.
        lqn.arrival = {}
        for entry in self.entries:
            arrival_dist = entry.getArrival()
            if arrival_dist is not None:
                entry_idx = self._entry_idx[entry]
                lqn.arrival[entry_idx] = arrival_dist

        # Build host demands
        # MATLAB sets hostdem for tasks (mean=0), entries (Immediate/mean=0), and activities.
        # lqn.hostdem carries the scalar mean (consumed numerically throughout
        # solver_ln/io); lqn.hostdem_proc carries the full Distribution object so
        # the real service SCV/phase-type reaches the layer models (moment2/moment3).
        lqn.hostdem = {}
        lqn.hostdem_proc = {}
        lqn.hostdem_scv = {}
        # Tasks have hostdem = 0 (they don't consume CPU directly)
        for task in self.tasks:
            task_idx = self._task_idx[task]
            lqn.hostdem[task_idx] = Immediate()
            lqn.hostdem_proc[task_idx] = Immediate()
            lqn.hostdem_scv[task_idx] = 0.0
        # Entries have hostdem = Immediate (mean=0)
        for entry in self.entries:
            entry_idx = self._entry_idx[entry]
            lqn.hostdem[entry_idx] = Immediate()
            lqn.hostdem_proc[entry_idx] = Immediate()
            lqn.hostdem_scv[entry_idx] = 0.0
        # Activities have actual host demands
        for act in self.activities:
            act_idx = self._act_idx[act]
            hd = act.host_demand
            lqn.hostdem[act_idx] = _get_dist_mean(hd)
            lqn.hostdem_proc[act_idx] = hd if hd is not None else Immediate()
            lqn.hostdem_scv[act_idx] = _get_dist_scv(hd)

        # Build activity think times (LQNX think-time on activities)
        # setThinkTime always stores a Distribution, including for a scalar
        # argument, so the duration must be read through its mean: comparing the
        # object itself against 0 raises TypeError and takes down every model
        # carrying an activity think time. (Restored after being clobbered by a
        # wholesale block rewrite in 6cc8ad3d2.)
        lqn.actthink = {}
        for act in self.activities:
            act_idx = self._act_idx[act]
            zt = act.think_time
            if zt is not None and hasattr(zt, 'getMean'):
                lqn.actthink[act_idx] = zt if zt.getMean() > 1e-8 else 0.0
            elif zt:
                lqn.actthink[act_idx] = float(zt) if float(zt) > 1e-8 else 0.0
            else:
                lqn.actthink[act_idx] = 0.0

        # Build think times
        lqn.think = {}
        for task in self.tasks:
            task_idx = self._task_idx[task]
            if task.think_time:
                # The Distribution itself, as in MATLAB lsn.think{} and the JAR:
                # a mean would drop the SCV of a non-exponential think time.
                lqn.think[task_idx] = task.think_time
                lqn.think_scv[task_idx] = _get_dist_scv(task.think_time)
            else:
                lqn.think[task_idx] = 0.0

        # Build setup times and delay-off times (setup/delay-off tasks)
        lqn.setuptime = {}
        lqn.delayofftime = {}
        lqn.hassetup = np.zeros((1, nidx + 1))
        for task in self.tasks:
            task_idx = self._task_idx[task]
            has_setup = task.setup_time is not None and _get_dist_mean(task.setup_time) > 0
            has_delayoff = task.delay_off_time is not None and _get_dist_mean(task.delay_off_time) > 0
            if has_setup or has_delayoff or task.has_setup_delayoff():
                lqn.hassetup[0, task_idx] = 1
                if task.setup_time:
                    lqn.setuptime[task_idx] = task.setup_time
                if task.delay_off_time:
                    lqn.delayofftime[task_idx] = task.delay_off_time

        # Build fan-out matrix from Task objects' fan-out maps
        lqn.fanout = np.zeros((nidx, nidx))
        task_name_to_idx = {}
        for task in self.tasks:
            task_name_to_idx[task.name] = self._task_idx[task]
        for task in self.tasks:
            tidx = self._task_idx[task]
            for dest_name, fo_val in task.getFanOut().items():
                if dest_name in task_name_to_idx:
                    dest_idx = task_name_to_idx[dest_name]
                    lqn.fanout[tidx, dest_idx] = fo_val

        # Build scheduling strategies
        lqn.sched = {}
        for proc in self.processors:
            lqn.sched[self._proc_idx[proc]] = proc.sched_strategy.value
        for task in self.tasks:
            lqn.sched[self._task_idx[task]] = task.sched_strategy.value

        # Build activity precedence type arrays. As in MATLAB getStruct they are
        # indexed by GLOBAL element index over 1..nidx, not by local activity
        # number: the graph, replygraph and every consumer are global too, so a
        # local space here would need an ashift correction at each read.
        # Values are the ActivityPrecedenceType ids shared with MATLAB and the JAR:
        # PRE_SEQ=1, PRE_AND=2, PRE_OR=3, POST_SEQ=11, POST_AND=12, POST_OR=13.
        lqn.actposttype = np.ones(nidx) * 11  # Default: POST_SEQ
        lqn.actpretype = np.ones(nidx) * 1    # Default: PRE_SEQ
        lqn.actquorum = np.zeros(nidx)

        # Populate actposttype and actpretype from precedence constraints
        for task in self.tasks:
            for prec in task.precedences:
                prec_type = self._get_prec_type(prec)
                if prec_type == PrecedenceType.PARALLEL:
                    # AND-fork: post_activities are POST_AND (only if multiple targets = fork, not join)
                    post_activities = self._get_prec_post_activities(prec)
                    if post_activities and len(post_activities) > 1:
                        for post_act in post_activities:
                            if post_act in self._act_idx:
                                act_idx = self._act_idx[post_act]
                                if ashift <= act_idx < ashift + lqn.nacts:
                                    lqn.actposttype[act_idx] = 12  # ID_POST_AND
                    # AND-join: pre_activities are PRE_AND
                    pre_activities = self._get_prec_pre_activities(prec)
                    if pre_activities and len(pre_activities) > 1:
                        for pre_act in pre_activities:
                            if pre_act in self._act_idx:
                                act_idx = self._act_idx[pre_act]
                                if ashift <= act_idx < ashift + lqn.nacts:
                                    lqn.actpretype[act_idx] = 2  # ID_PRE_AND
                        # Record the quorum on the join target. A missing or out-of-range
                        # value means the join waits for all its predecessors.
                        nbranches = len(pre_activities)
                        quorum = nbranches
                        pre_params = getattr(prec, 'pre_params', None)
                        if pre_params is not None and np.size(pre_params) == 1:
                            q = int(round(float(np.ravel(pre_params)[0])))
                            if 1 <= q <= nbranches:
                                quorum = q
                        for post_act in self._get_prec_post_activities(prec):
                            if post_act in self._act_idx:
                                act_idx = self._act_idx[post_act]
                                if ashift <= act_idx < ashift + lqn.nacts:
                                    lqn.actquorum[act_idx] = quorum
                elif prec_type == PrecedenceType.CHOICE:
                    # OR-fork: post_activities are POST_OR (only if multiple targets = fork, not join)
                    post_activities = self._get_prec_post_activities(prec)
                    if post_activities and len(post_activities) > 1:
                        for post_act in post_activities:
                            if post_act in self._act_idx:
                                act_idx = self._act_idx[post_act]
                                if ashift <= act_idx < ashift + lqn.nacts:
                                    lqn.actposttype[act_idx] = 13  # ID_POST_OR
                    # OR-join: pre_activities are PRE_OR
                    pre_activities = self._get_prec_pre_activities(prec)
                    if pre_activities and len(pre_activities) > 1:
                        for pre_act in pre_activities:
                            if pre_act in self._act_idx:
                                act_idx = self._act_idx[pre_act]
                                if ashift <= act_idx < ashift + lqn.nacts:
                                    lqn.actpretype[act_idx] = 3  # ID_PRE_OR

        # Build graph (adjacency matrix)
        # MATLAB convention: graph[child, parent] = 1 (element points to its parent/owner)
        lqn.graph = np.zeros((nidx, nidx))
        # Track loop-back edges for DAG construction (matches MATLAB loop_back_edges)
        loop_back_edges = np.zeros((nidx, nidx), dtype=bool)
        # Add edges from tasks to processors (task points to its parent processor)
        for task in self.tasks:
            task_idx = self._task_idx[task]
            if task.processor:
                proc_idx = self._proc_idx[task.processor]
                lqn.graph[task_idx, proc_idx] = 1  # task -> processor (parent)
        # Add edges from tasks to entries (task points to its entries)
        for entry in self.entries:
            entry_idx = self._entry_idx[entry]
            if entry.task:
                task_idx = self._task_idx[entry.task]
                lqn.graph[task_idx, entry_idx] = 1  # task -> entry (ownership)
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

        # Note: We DO NOT add entry->task edges because combined with task->entry edges,
        # they would create cycles that break the topological sort in lsn_max_multiplicity.
        # Instead, the task maxmult is handled specially in lsn_max_multiplicity.

        # Add edges for precedence constraints (serial, and-fork, and-join, etc.)
        for task in self.tasks:
            for prec in task.precedences:
                prec_type = self._get_prec_type(prec)
                activities = self._get_prec_activities(prec)
                pre_activities = self._get_prec_pre_activities(prec)
                post_activities = self._get_prec_post_activities(prec)
                probabilities = self._get_prec_probabilities(prec)

                if prec_type == PrecedenceType.SERIAL and activities:
                    # Serial chain: A -> B -> C creates edges A->B and B->C
                    for i in range(len(activities) - 1):
                        pre_act = activities[i]
                        post_act = activities[i + 1]
                        if self._act_in_idx(pre_act) and self._act_in_idx(post_act):
                            pre_aidx = self._get_act_idx(pre_act)
                            post_aidx = self._get_act_idx(post_act)
                            lqn.graph[pre_aidx, post_aidx] = 1.0
                elif prec_type == PrecedenceType.PARALLEL:
                    # AND-fork/join: pre_activities -> post_activities
                    if pre_activities and post_activities:
                        for pre_act in pre_activities:
                            for post_act in post_activities:
                                if self._act_in_idx(pre_act) and self._act_in_idx(post_act):
                                    pre_aidx = self._get_act_idx(pre_act)
                                    post_aidx = self._get_act_idx(post_act)
                                    lqn.graph[pre_aidx, post_aidx] = 1.0
                elif prec_type == PrecedenceType.CHOICE:
                    # OR-fork/join: pre_activities -> post_activities with probabilities
                    if pre_activities and post_activities:
                        probs = probabilities if probabilities else [1.0 / len(post_activities)] * len(post_activities)
                        for pre_act in pre_activities:
                            for i, post_act in enumerate(post_activities):
                                if self._act_in_idx(pre_act) and self._act_in_idx(post_act):
                                    pre_aidx = self._get_act_idx(pre_act)
                                    post_aidx = self._get_act_idx(post_act)
                                    prob = probs[i] if i < len(probs) else 1.0 / len(post_activities)
                                    lqn.graph[pre_aidx, post_aidx] = prob
                elif prec_type == PrecedenceType.CACHE_ACCESS:
                    # Cache access: access_act -> [hit_act, miss_act] with probabilities
                    # Compute expected hit/miss probabilities from CacheTask configuration
                    hit_prob = 0.5  # Default if cache config not available
                    miss_prob = 0.5
                    if isinstance(task, CacheTask):
                        # see _kb/04-networkstruct.md (Python native layered.py port notes) for rationale
                        cache_cap = task.cache_capacity
                        if isinstance(cache_cap, (list, tuple)):
                            cache_cap = sum(cache_cap) if cache_cap else 0
                        total_items = task.total_items
                        if isinstance(total_items, (list, tuple)):
                            total_items = total_items[0] if total_items else 1
                        if total_items > 0:
                            hit_prob = min(float(cache_cap) / float(total_items), 1.0)
                            miss_prob = 1.0 - hit_prob
                    if pre_activities and post_activities:
                        for pre_act in pre_activities:
                            # post_activities is [hit_activity, miss_activity]
                            for i, post_act in enumerate(post_activities):
                                if self._act_in_idx(pre_act) and self._act_in_idx(post_act):
                                    pre_aidx = self._get_act_idx(pre_act)
                                    post_aidx = self._get_act_idx(post_act)
                                    # First post_activity is hit, second is miss
                                    prob = hit_prob if i == 0 else miss_prob
                                    lqn.graph[pre_aidx, post_aidx] = prob
                elif prec_type == PrecedenceType.LOOP:
                    # see _kb/04-networkstruct.md (Python native layered.py port notes) for rationale
                    counts = self._get_prec_count(prec)
                    loop_pre_activities = pre_activities if pre_activities else self._get_prec_pre_activities(prec)
                    loop_activities = activities if activities else self._get_prec_activities(prec)
                    if loop_pre_activities and loop_activities and len(loop_activities) >= 1:
                        loop_entry_act = loop_pre_activities[0]
                        # activities = [loop_body_activities..., end_activity]
                        # For Loop(A1, [A2, A3], 3): loop body = [A2], end = A3
                        loop_body_acts = loop_activities[:-1] if len(loop_activities) > 1 else []
                        loop_end_act = loop_activities[-1]

                        if self._act_in_idx(loop_entry_act):
                            loop_entry_aidx = self._get_act_idx(loop_entry_act)

                            if counts < 1:
                                # When expected iterations < 1, we may skip loop entirely
                                # E[iterations] = counts means P(enter loop) = counts
                                if loop_body_acts and self._act_in_idx(loop_body_acts[0]):
                                    first_loop_aidx = self._get_act_idx(loop_body_acts[0])
                                    lqn.graph[loop_entry_aidx, first_loop_aidx] = counts
                                if self._act_in_idx(loop_end_act):
                                    loop_end_aidx = self._get_act_idx(loop_end_act)
                                    lqn.graph[loop_entry_aidx, loop_end_aidx] = 1.0 - counts
                                # Connect loop body in serial
                                cur_aidx = first_loop_aidx if loop_body_acts else loop_entry_aidx
                                for i in range(1, len(loop_body_acts)):
                                    if self._act_in_idx(loop_body_acts[i]):
                                        next_aidx = self._get_act_idx(loop_body_acts[i])
                                        lqn.graph[cur_aidx, next_aidx] = 1.0
                                        cur_aidx = next_aidx
                                # After loop body, exit to end (no looping back for counts < 1)
                                if loop_body_acts and self._act_in_idx(loop_end_act):
                                    lqn.graph[cur_aidx, self._get_act_idx(loop_end_act)] = 1.0
                            else:
                                # see _kb/04-networkstruct.md (Python native layered.py port notes) for rationale
                                cur_aidx = loop_entry_aidx
                                for loop_act in loop_body_acts:
                                    if self._act_in_idx(loop_act):
                                        next_aidx = self._get_act_idx(loop_act)
                                        lqn.graph[cur_aidx, next_aidx] = 1.0
                                        cur_aidx = next_aidx

                                # If no loop body, entry connects directly to end
                                if not loop_body_acts:
                                    if self._act_in_idx(loop_end_act):
                                        lqn.graph[loop_entry_aidx, self._get_act_idx(loop_end_act)] = 1.0
                                else:
                                    # Loop back edge: last_body -> first_body with prob (1 - 1/counts)
                                    first_loop_aidx = self._get_act_idx(loop_body_acts[0])
                                    loop_back_prob = 1.0 - 1.0 / counts
                                    exit_prob = 1.0 / counts
                                    lqn.graph[cur_aidx, first_loop_aidx] = loop_back_prob
                                    # Mark this as a loop-back edge for DAG construction
                                    loop_back_edges[cur_aidx, first_loop_aidx] = True
                                    # Exit edge: last_body -> end with prob 1/counts
                                    if self._act_in_idx(loop_end_act):
                                        loop_end_aidx = self._get_act_idx(loop_end_act)
                                        lqn.graph[cur_aidx, loop_end_aidx] = exit_prob

        # Build DAG (graph without loop-back edges) - matches MATLAB getStruct.m line 631
        # Add forwarding edges to the element graph (MATLAB getStruct.m:
        # lsn.graph(eidx, target_eidx) = 1). The actsof BFS below filters by
        # parent task, so the cross-task FWD edge does not leak activities.
        for cidx in range(lqn.ncalls):
            if cidx < len(lqn.calltype) and int(lqn.calltype[cidx]) == 3:  # CallType.FWD
                fwd_src_eidx = int(lqn.callpair[cidx, 0])
                fwd_tgt_eidx = int(lqn.callpair[cidx, 1])
                if fwd_src_eidx > 0 and fwd_tgt_eidx > 0:
                    lqn.graph[fwd_src_eidx, fwd_tgt_eidx] = 1

        # Refine actsof for entries into the entry-to-activity reachability
        # closure, now that graph carries the precedences. The bound activity
        # alone is not the entry's activity set: a Serial chain puts the rest of
        # the work on successors that are still executed by the entry, so LQNS
        # and SolverLN both account them to it. Restricting the closure to the
        # entry's own task drops the callee entries reached through call edges.
        # Matches MATLAB getStruct.m and the JAR/C++ readers.
        for entry in self.entries:
            eidx = self._entry_idx[entry]
            tidx = int(lqn.parent[eidx, 0])
            visited = np.zeros(nidx, dtype=bool)
            visited[eidx] = True
            stack = [eidx]
            while stack:
                v = stack.pop()
                for w in np.flatnonzero(lqn.graph[v, :]):
                    if not visited[w]:
                        visited[w] = True
                        stack.append(int(w))
            lqn.actsof[eidx] = [int(i) for i in range(nidx)
                                if visited[i]
                                and lqn.type[i] == LayeredNetworkElement.ACTIVITY
                                and int(lqn.parent[i, 0]) == tidx]

        lqn.dag = lqn.graph.copy()
        lqn.dag[loop_back_edges] = 0

        # Reverse edges from TASK to ENTRY for non-reference tasks
        # This enables proper flow propagation in _lsn_max_multiplicity
        # Matches MATLAB getStruct.m lines 598-606
        for i in range(nidx):
            if lqn.type[i] == LayeredNetworkElement.TASK and lqn.isref[i, 0] == 0:
                for j in range(nidx):
                    if lqn.type[j] == LayeredNetworkElement.ENTRY and lqn.dag[i, j] != 0:
                        lqn.dag[i, j] = 0
                        lqn.dag[j, i] = 1

        # Build replygraph (nacts x nentries) - which activities reply to which entries
        lqn.replygraph = np.zeros((lqn.nacts, lqn.nentries))
        for act in self.activities:
            act_idx = self._act_idx[act] - ashift  # activity-local index
            # Check for reply_entry (singular) - set by replies_to() method
            if hasattr(act, 'reply_entry') and act.reply_entry is not None:
                reply_entry = act.reply_entry
                if reply_entry in self._entry_idx:
                    entry_idx = self._entry_idx[reply_entry] - eshift  # entry-local index
                    if 0 <= act_idx < lqn.nacts and 0 <= entry_idx < lqn.nentries:
                        lqn.replygraph[act_idx, entry_idx] = 1
            # Also check reply_entries (plural) for backwards compatibility
            if hasattr(act, 'reply_entries') and act.reply_entries:
                for reply_entry in act.reply_entries:
                    if reply_entry in self._entry_idx:
                        entry_idx = self._entry_idx[reply_entry] - eshift  # entry-local index
                        if 0 <= act_idx < lqn.nacts and 0 <= entry_idx < lqn.nentries:
                            lqn.replygraph[act_idx, entry_idx] = 1

        # Build actphase - phase number (1..3) for each activity
        lqn.actphase = np.ones(lqn.nacts)  # Default phase is 1
        for act in self.activities:
            act_idx = self._act_idx[act] - ashift  # activity-local index
            if hasattr(act, 'phase') and act.phase is not None:
                if 0 <= act_idx < lqn.nacts:
                    lqn.actphase[act_idx] = act.phase

        # Rebuild actsof for entries using BFS through graph
        # Matches MATLAB getStruct.m lines 640-658
        # The initial actsof only includes bound activities, but entries may have
        # additional reachable activities (e.g., phase 2 post-reply activities)
        for eoff in range(lqn.nentries):
            eidx = eshift + eoff
            tidx = int(lqn.parent[eidx, 0])
            visited = set()
            stack = [eidx]
            visited.add(eidx)
            while stack:
                v = stack.pop()
                for nbr in range(nidx):
                    if nbr not in visited and lqn.graph[v, nbr] != 0:
                        visited.add(nbr)
                        stack.append(nbr)
            acts = [idx for idx in visited
                    if lqn.type[idx] == LayeredNetworkElement.ACTIVITY
                    and int(lqn.parent[idx, 0]) == tidx]
            lqn.actsof[eidx] = sorted(acts)

        # Correct multiplicity for INF-scheduled tasks (MATLAB getStruct.m lines 576-589)
        # For INF tasks, set mult = sum of caller task multiplicities
        # This must happen BEFORE maxmult computation
        # Process in index order so callers are corrected before callees (matches MATLAB)
        inf_task_indices = []
        for task in self.tasks:
            tidx = self._task_idx[task]
            if task.sched_strategy == SchedStrategy.INF and lqn.type[tidx] == LayeredNetworkElement.TASK:
                    inf_task_indices.append(tidx)
        inf_task_indices.sort()  # Process in increasing index order
        for tidx in inf_task_indices:
            # Find caller tasks using iscaller matrix
            caller_tasks = []
            for other_task in self.tasks:
                other_tidx = self._task_idx[other_task]
                if other_tidx != tidx and lqn.iscaller[other_tidx, tidx] > 0:
                    caller_tasks.append(other_tidx)
            # Forwarding sources also feed jobs into the target task (MATLAB
            # counts them via the taskgraph FWD edge; python has no taskgraph)
            for cidx in range(lqn.ncalls):
                if cidx < len(lqn.calltype) and int(lqn.calltype[cidx]) == 3:  # CallType.FWD
                    fwd_tgt_tidx = int(lqn.parent[int(lqn.callpair[cidx, 1]), 0])
                    if fwd_tgt_tidx == tidx:
                        fwd_src_tidx = int(lqn.parent[int(lqn.callpair[cidx, 0]), 0])
                        if fwd_src_tidx != tidx and fwd_src_tidx not in caller_tasks:
                            caller_tasks.append(fwd_src_tidx)
            if len(caller_tasks) > 0:
                caller_mults = [lqn.mult[0, ct] for ct in caller_tasks]
                has_inf_caller = any(np.isinf(m) for m in caller_mults)
                if has_inf_caller:
                    # Heuristic: sum finite + inf_count * max(all mult)
                    finite_sum = sum(m for m in caller_mults if not np.isinf(m))
                    inf_count = sum(1 for m in caller_mults if np.isinf(m))
                    max_mult = float(np.max(lqn.mult[0, :]))
                    lqn.mult[0, tidx] = finite_sum + inf_count * max_mult
                else:
                    lqn.mult[0, tidx] = sum(caller_mults)

        # Build maxmult using lsn_max_multiplicity algorithm (matches MATLAB)
        # This computes the maximum throughput capacity for each node
        # Must be called AFTER graph is built since it uses the DAG for flow propagation
        lqn.maxmult = self._lsn_max_multiplicity(lqn, nidx)

        # Validation: every entry must have a boundTo activity.
        # getStruct.m's guard, absent here until 2026-08-15. An entry with an
        # empty <entry-phase-activities> reaches no activity, so it has no
        # service and no reply; it used to solve and report a row of NaN
        # instead of being named. Checked before the reply guard below, the
        # order getStruct.m refuses them in.
        for e in range(lqn.nentries):
            eidx = eshift + e
            bound = False
            for succ in range(ashift, nidx):
                if lqn.graph[eidx, succ] != 0:
                    bound = True
                    break
            if not bound:
                raise ValueError('An entry does not have any boundTo activity.')

        # Validation: Check for non-terminal reply activities
        # An activity that replies to an entry should not have Phase 1 successor activities
        # Phase 2 successors are allowed (post-reply processing)
        for a in range(lqn.nacts):
            if np.any(lqn.replygraph[a, :] > 0):  # activity 'a' replies to some entry
                aidx = ashift + a  # global activity index
                for succ in range(nidx):
                    if lqn.graph[aidx, succ] != 0:  # successor exists
                        if succ >= ashift:  # successor is an activity
                            succ_act_idx = succ - ashift  # activity-local index
                            if 0 <= succ_act_idx < lqn.nacts and lqn.actphase[succ_act_idx] == 1:
                                raise ValueError(
                                    f"Unsupported replyTo in non-terminal activity: "
                                    f"activity '{lqn.names[aidx]}' has Phase 1 successor"
                                )

        return lqn

    def write_xml(self, filename: str, use_abstract_names: bool = False) -> None:
        """Write the layered network to LQNX XML (snake_case alias for writeXML)."""
        return self.writeXML(filename, use_abstract_names)

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
            # FCFSPRPRIO is LINE's reading of the LQN "pri" discipline; its own
            # name is not valid LQN, so write back the spelling lqns accepts
            sched_val = sched.value if hasattr(sched, 'value') else sched
            if sched_val == SchedStrategy.FCFSPRPRIO.value:
                return 'pri'
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

        def is_ps_sched(sched):
            """Check if scheduling is PS or PSPRIO (the quantum-bearing hosts)."""
            ps_vals = (SchedStrategy.PS.value, SchedStrategy.PSPRIO.value)
            if isinstance(sched, SchedStrategy):
                return sched.value in ps_vals
            elif hasattr(sched, 'value'):
                return sched.value in ps_vals
            else:
                return sched in ps_vals

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

            # replication is read back by parseXML and by the C++ lqnx reader,
            # so dropping it here silently solves an unreplicated model
            if proc.getReplication() > 1:
                proc_elem.set('replication', str(int(proc.getReplication())))

            if not is_inf_sched(proc.sched_strategy):
                mult = proc.multiplicity
                if np.isinf(mult):
                    mult = 1
                proc_elem.set('multiplicity', str(int(mult)))

            # only when set: parseXML falls back to 0.001 on an absent quantum,
            # and writing python's 0.0 default would overwrite that fallback
            if is_ps_sched(proc.sched_strategy) and proc.getQuantum() > 0:
                proc_elem.set('quantum', str(proc.getQuantum()))

            proc_elem.set('speed-factor', '1')

            # Write tasks for this processor
            for task in proc.tasks:
                task_elem = ET.SubElement(proc_elem, 'task')
                task_elem.set('name', name_map[task.name])
                task_elem.set('scheduling', sched_to_text(task.sched_strategy))
                if task.get_priority() != 0:
                    task_elem.set('priority', str(int(task.get_priority())))

                if task.getReplication() > 1:
                    task_elem.set('replication', str(int(task.getReplication())))

                if not is_inf_sched(task.sched_strategy):
                    mult = task.multiplicity
                    if np.isinf(mult):
                        mult = 1000000  # Large number for infinite
                    task_elem.set('multiplicity', str(int(mult)))

                # think-time may only be written for a reference task: lqns
                # rejects the file outright ('Task "X" is not a reference task;
                # it cannot have think time'). LINE does give a non-reference
                # task's think time to its callers, so that value cannot survive
                # this format and the loss is reported rather than left silent.
                if task.think_time:
                    if is_ref_sched(task.sched_strategy):
                        task_elem.set('think-time', str(get_dist_mean(task.think_time)))
                    elif get_dist_mean(task.think_time) > 0:
                        from .api.io.logging import line_warning
                        line_warning('writeXML',
                                     'Task %s is not a reference task, so its think time (%g) '
                                     'is not written: the LQN XML schema accepts think-time '
                                     'on reference tasks only.',
                                     task.name, get_dist_mean(task.think_time))

                # LINE .lqnx dialect: a cache, and a setup / delay-off time, have
                # no element in the stock LQN schema, so the file used to describe
                # a DIFFERENT model -- a CacheTask degenerated to a plain task and
                # its hit/miss split to an unweighted 50/50 <post>. The three
                # elements below carry the same field set the JSON interchange
                # carries (linemodel_save.m: totalItems, cacheCapacity,
                # replacementStrategy, entryType=ItemEntry, accessProb), so the
                # two transports are informationally equal. lqns/lqsim ignore
                # unknown children, and so does every LINE reader.
                _write_cache_elem(task_elem, task)
                _write_setup_elems(task_elem, task)

                # fan-out/fan-in are parsed by every reader (parseXML, the JAR,
                # cpp/lqn_reader.h) and were written by none, so a replicated
                # model lost its call multiplicities on every round trip. The
                # schema (lqn-core.xsd, TaskType) places them before the entries.
                for dest_name, fo_val in task.getFanOut().items():
                    fo_elem = ET.SubElement(task_elem, 'fan-out')
                    fo_elem.set('dest', name_map.get(dest_name, dest_name))
                    fo_elem.set('value', str(int(fo_val)))
                for src_name, fi_val in task.getFanIn().items():
                    fi_elem = ET.SubElement(task_elem, 'fan-in')
                    fi_elem.set('source', name_map.get(src_name, src_name))
                    fi_elem.set('value', str(int(fi_val)))

                # Write entries for this task
                for entry in task.entries:
                    entry_elem = ET.SubElement(task_elem, 'entry')
                    entry_elem.set('name', name_map[entry.name])
                    entry_elem.set('type', 'NONE')
                    # LINE dialect: an ItemEntry's cardinality and access
                    # popularity, without which the entry loads as a plain Entry.
                    _write_item_entry_elem(entry_elem, entry)
                    # Emit open-arrival-rate if entry has an arrival distribution set
                    arrival_dist = entry.getArrival()
                    if arrival_dist is not None:
                        mean = get_dist_mean(arrival_dist)
                        if mean and mean > 0:
                            entry_elem.set('open-arrival-rate', str(1.0 / mean))
                    # Emit forwarding elements (matches MATLAB writeXML; without
                    # these lqns solves a model with the forwarding dropped)
                    for fwd_i, fwd_dest in enumerate(entry._forwarding_dests):
                        fwd_elem = ET.SubElement(entry_elem, 'forwarding')
                        fwd_elem.set('dest', name_map.get(fwd_dest, fwd_dest))
                        fwd_elem.set('prob', str(entry._forwarding_probs[fwd_i]))

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

                    # LINE dialect: a routed call group names which of the
                    # synch-calls above one dispatcher issues, and under which
                    # strategy. The member calls stay ordinary synch-calls, so a
                    # reader that ignores this element still sees the same
                    # aggregate call means -- which is what lqns and lqsim,
                    # having no dispatcher, should see.
                    for strategy, group_entries in getattr(act, 'call_groups', []):
                        grp_elem = ET.SubElement(act_elem, 'call-group')
                        grp_elem.set('strategy', _callgroup_to_lqnx(strategy))
                        for target_entry in group_entries:
                            dest_elem = ET.SubElement(grp_elem, 'dest')
                            dest_elem.set('name', name_map[target_entry.name])

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
                                # Emit the quorum only for a genuine quorum k < n; a join that
                                # waits for all predecessors is the LQN default.
                                pre_params = getattr(prec, 'pre_params', None)
                                if pre_params is not None and np.size(pre_params) == 1:
                                    q = int(round(float(np.ravel(pre_params)[0])))
                                    if 1 <= q < len(prec.pre_activities):
                                        pre_elem.set('quorum', str(q))
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
                            elif prec.prec_type == PrecedenceType.CACHE_ACCESS:
                                # A cache access wrote a bare <post>, which every
                                # reader takes as an unweighted split -- that is
                                # where a spurious "exact 0.5" hit ratio came from.
                                post_elem = ET.SubElement(prec_elem, 'post-CACHE')
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
                            elif prec.prec_type == PrecedenceType.CACHE_ACCESS:
                                # explicit, not positional: hit first, miss second
                                act_ref.set('cache-result', 'hit' if i == 0 else 'miss')

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

    def getName(self) -> str:
        """Get the model name, as Model.getName does in MATLAB and the JAR."""
        return self.name

    get_name = getName

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

    # Aggregate flat interface (Network-compatible) used by SolverENV.
    # see _kb/04-networkstruct.md (Python native layered.py port notes) for rationale
    def _layer_ensemble(self):
        """Return the layer submodels, building them via SolverLN if needed."""
        ens = getattr(self, 'ensemble', None)
        if ens is None or len(ens) == 0:
            from .solvers.solver_ln.solver_ln import SolverLN
            SolverLN(self)  # _build_layers sets self.ensemble
            ens = getattr(self, 'ensemble', None)
        return ens if ens is not None else []

    def _layer_blocks(self):
        """(Roff, Coff, Msz, Ksz): 0-based row/col offset and station/class
        count of each layer's block in the aggregate layout."""
        ens = self._layer_ensemble()
        Msz = [int(e.get_number_of_stations()) for e in ens]
        Ksz = [int(e.get_number_of_classes()) for e in ens]
        Roff = [0] * len(ens)
        Coff = [0] * len(ens)
        for e in range(1, len(ens)):
            Roff[e] = Roff[e - 1] + Msz[e - 1]
            Coff[e] = Coff[e - 1] + Ksz[e - 1]
        return Roff, Coff, Msz, Ksz

    def get_number_of_stations(self) -> int:
        """Aggregate (sum over layers) station count."""
        _, _, Msz, _ = self._layer_blocks()
        return int(sum(Msz))

    def get_number_of_classes(self) -> int:
        """Aggregate (sum over layers) class count."""
        _, _, _, Ksz = self._layer_blocks()
        return int(sum(Ksz))

    def get_number_of_nodes(self) -> int:
        """Aggregate (sum over layers) node count."""
        return int(sum(e.get_number_of_nodes() for e in self._layer_ensemble()))

    def get_number_of_stateful_nodes(self) -> int:
        """Aggregate (sum over layers) stateful-node count."""
        total = 0
        for e in self._layer_ensemble():
            if hasattr(e, 'get_number_of_stateful_nodes'):
                total += int(e.get_number_of_stateful_nodes())
            elif hasattr(e, 'getNumberOfStatefulNodes'):
                total += int(e.getNumberOfStatefulNodes())
        return total

    def get_tran_handles(self):
        """Block-diagonal aggregate transient handles over the layer networks.
        Off-block cells are None; SolverLN.getTranAvg reproduces the same
        layout so the handles only fix the aggregate M x K shape."""
        ens = self._layer_ensemble()
        Roff, Coff, Msz, Ksz = self._layer_blocks()
        M = int(sum(Msz))
        K = int(sum(Ksz))
        Qt = [[None] * K for _ in range(M)]
        Ut = [[None] * K for _ in range(M)]
        Tt = [[None] * K for _ in range(M)]
        for e, layer in enumerate(ens):
            Qe, Ue, Te = layer.get_tran_handles()
            if Qe is None:
                continue
            for i in range(Msz[e]):
                for r in range(Ksz[e]):
                    Qt[Roff[e] + i][Coff[e] + r] = Qe[i][r]
                    if Ue is not None:
                        Ut[Roff[e] + i][Coff[e] + r] = Ue[i][r]
                    if Te is not None:
                        Tt[Roff[e] + i][Coff[e] + r] = Te[i][r]
        return Qt, Ut, Tt

    def init_from_marginal(self, n, options=None) -> None:
        """Split the aggregate (M x K) marginal queue-length matrix into
        per-layer blocks and warm-start each layer network. Mirrors MATLAB
        LayeredNetwork.initFromMarginal; the cross-switch state-continuity
        mechanism for SolverENV."""
        import numpy as _np
        n = _np.atleast_2d(_np.asarray(n, dtype=float))
        ens = self._layer_ensemble()
        Roff, Coff, Msz, Ksz = self._layer_blocks()
        blocks = []
        for e, layer in enumerate(ens):
            block = n[Roff[e]:Roff[e] + Msz[e], Coff[e]:Coff[e] + Ksz[e]]
            if options is not None:
                layer.init_from_marginal(block, options)
            else:
                layer.init_from_marginal(block)
            blocks.append(block.copy())
        self._init_marginal_blocks = blocks

    def get_init_marginal_blocks(self):
        """Per-layer blocks of the warm start supplied by the last
        init_from_marginal call, or None if none was supplied.

        The layered fixed point hard-resets every layer when it detects
        convergence, which discards the warm start, so SolverLN replays these
        blocks just before running the layer transients; that is what lets a
        stage of SolverENV resume from the marginal handed over at the
        environment switch."""
        return getattr(self, '_init_marginal_blocks', None)

    # MATLAB-style camelCase aliases
    getNumberOfStations = get_number_of_stations
    getNumberOfClasses = get_number_of_classes
    getNumberOfNodes = get_number_of_nodes
    getNumberOfStatefulNodes = get_number_of_stateful_nodes
    getTranHandles = get_tran_handles
    initFromMarginal = init_from_marginal
    getInitMarginalBlocks = get_init_marginal_blocks

    def plotGraph(self, method: str = 'nodes', ax=None, show: bool = True):
        """Plot the layered-network call graph (native equivalent of MATLAB
        LayeredNetwork.plotGraph).

        Nodes are laid out in type layers and colored by element type: hosts
        (black), tasks (magenta; reference tasks gold), entries (red),
        activities (blue).

        Args:
            method: label source -- 'nodes'/'names' use hashnames/names, 'ids' uses indices.
            ax: optional matplotlib Axes to draw into.
            show: call plt.show() when True.

        Returns:
            The matplotlib Axes containing the plot.
        """
        import numpy as np
        import matplotlib.pyplot as plt
        import networkx as nx

        lqn = self.getStruct()
        T = np.asarray(lqn.graph)
        n = int(lqn.nidx)

        def _label(idx):
            if method == 'ids':
                return str(idx)
            arr = lqn.names if method == 'names' else lqn.hashnames
            if arr is not None and idx < len(arr) and arr[idx] is not None:
                return str(arr[idx])
            return str(idx)

        def _kind(idx):
            if lqn.hshift <= idx < lqn.hshift + lqn.nhosts:
                return 'host'
            if lqn.tshift <= idx < lqn.tshift + lqn.ntasks:
                return 'task'
            if lqn.eshift <= idx < lqn.eshift + lqn.nentries:
                return 'entry'
            if lqn.ashift <= idx < lqn.ashift + lqn.nacts:
                return 'act'
            return 'other'

        layer = {'host': 0, 'task': 1, 'entry': 2, 'act': 3, 'other': 4}
        color = {'host': 'black', 'task': 'magenta', 'entry': 'red', 'act': 'tab:blue', 'other': 'gray'}

        G = nx.DiGraph()
        for idx in range(n):
            k = _kind(idx)
            G.add_node(idx, layer=layer[k], kind=k)
        rows, cols = np.nonzero(T)
        for i, j in zip(rows, cols):
            if 0 <= i < n and 0 <= j < n:
                G.add_edge(int(i), int(j))

        try:
            pos = nx.multipartite_layout(G, subset_key='layer')
        except Exception:
            pos = nx.spring_layout(G, seed=0)

        node_colors = []
        for idx in G.nodes():
            k = G.nodes[idx]['kind']
            if k == 'task' and lqn.isref is not None and idx < len(lqn.isref) and lqn.isref[idx]:
                node_colors.append('gold')
            else:
                node_colors.append(color[k])

        if ax is None:
            _, ax = plt.subplots(figsize=(9, 6))
        nx.draw_networkx_edges(G, pos, ax=ax, arrows=True, edge_color='0.5')
        nx.draw_networkx_nodes(G, pos, ax=ax, node_color=node_colors, node_size=500)
        nx.draw_networkx_labels(G, pos, ax=ax, labels={idx: _label(idx) for idx in G.nodes()}, font_size=8)
        ax.set_title('Model: %s' % self.name)
        ax.axis('off')
        if show:
            plt.show()
        return ax

    def plot(self, show_task_graph: bool = False, show: bool = True):
        """Plot the layered-network graph (native equivalent of MATLAB plot)."""
        return self.plotGraph(show=show)

    def view(self, show: bool = True):
        """Visualize the layered network (native equivalent of MATLAB view)."""
        return self.plotGraph(show=show)

    def get_layers(self) -> list:
        """Get layers (snake_case alias for getLayers)."""
        return self.getLayers()

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
        cls._validate_input_model(root)

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
                # MATLAB overrides multiplicity to Inf for INF scheduling processors
                # (shows a warning but still uses Inf multiplicity)
                multiplicity = float('inf')
                sched = SchedStrategy.INF
            else:
                sched = cls._parse_sched_strategy(scheduling)

            processor = model.add_processor(proc_name, multiplicity, sched)

            # Parse quantum (default 0.001 per Java)
            quantum_str = proc_elem.get('quantum', '')
            if quantum_str:
                try:
                    processor.setQuantum(float(quantum_str))
                except ValueError:
                    pass

            # Parse speed-factor (default 1.0)
            speed_str = proc_elem.get('speed-factor', '')
            if speed_str:
                try:
                    processor.setSpeedFactor(float(speed_str))
                except ValueError:
                    pass

            # Parse processor replication (default 1)
            proc_repl_str = proc_elem.get('replication', '1')
            try:
                proc_replication = int(float(proc_repl_str))
            except ValueError:
                proc_replication = 1
            if proc_replication > 1:
                processor.setReplication(proc_replication)

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
                    # MATLAB overrides task_mult to Inf for INF scheduling tasks
                    # (shows a warning but still uses Inf multiplicity)
                    # This is important for njobs calculation which uses mult for population
                    task_mult = float('inf')
                    task_sched_enum = SchedStrategy.INF
                else:
                    task_sched_enum = cls._parse_sched_strategy(task_sched)

                # LINE .lqnx dialect: the presence of <cache> is what makes the
                # task a CacheTask, as taskType=CacheTask does on the JSON wire.
                cache_elem = task_elem.find('./cache')
                if cache_elem is not None:
                    caps = [int(float(lv.get('capacity', '1')))
                            for lv in cache_elem.findall('./level')] or [1]
                    task = CacheTask(model, task_name,
                                     int(float(cache_elem.get('items', '1'))),
                                     caps if len(caps) > 1 else caps[0],
                                     cls._parse_replacement(cache_elem.get('replacement', 'FIFO')),
                                     task_mult)
                    task.retrieval = cache_elem.get('retrieval', 'false').lower() == 'true'
                    task.on(processor)
                else:
                    task = model.add_task(task_name, task_mult, task_sched_enum, processor)

                # LINE .lqnx dialect: setup / delay-off times, rebuilt through the
                # same setters the model API uses, so a mean plus an SCV yields the
                # distribution family the setter would have produced.
                setup_elem = task_elem.find('./setup')
                if setup_elem is not None:
                    task.setSetupTime(_dist_from_mean_scv(
                        float(setup_elem.get('mean', '0')),
                        float(setup_elem.get('scv', '1'))))
                off_elem = task_elem.find('./delay-off')
                if off_elem is not None:
                    task.setDelayOffTime(_dist_from_mean_scv(
                        float(off_elem.get('mean', '0')),
                        float(off_elem.get('scv', '1'))))

                # Task scheduling priority (lower is served first)
                task_prio_str = task_elem.get('priority', '')
                if task_prio_str:
                    try:
                        task.set_priority(int(float(task_prio_str)))
                    except ValueError:
                        pass

                # Parse task replication (Java uses processor replication for task)
                task_repl_str = task_elem.get('replication', '')
                if task_repl_str:
                    try:
                        task_repl = int(float(task_repl_str))
                    except ValueError:
                        task_repl = proc_replication
                else:
                    task_repl = proc_replication
                if task_repl > 1:
                    task.setReplication(task_repl)

                # Parse fan-in elements
                for fan_in_elem in task_elem.findall('./fan-in'):
                    source = fan_in_elem.get('source', '')
                    value_str = fan_in_elem.get('value', '')
                    if source and value_str:
                        try:
                            task.setFanIn(source, int(value_str))
                        except ValueError:
                            pass

                # Parse fan-out elements
                for fan_out_elem in task_elem.findall('./fan-out'):
                    dest = fan_out_elem.get('dest', '')
                    value_str = fan_out_elem.get('value', '')
                    if dest and value_str:
                        try:
                            task.setFanOut(dest, int(value_str))
                        except ValueError:
                            pass

                think_time_str = task_elem.get('think-time', '0')
                try:
                    think_time = float(think_time_str)
                    if think_time > 0:
                        task.set_think_time(think_time)
                except ValueError:
                    pass

                for entry_elem in task_elem.findall('./entry'):
                    entry_name = entry_elem.get('name', '')
                    # LINE .lqnx dialect: the presence of <item-entry> is what
                    # makes the entry an ItemEntry, as entryType=ItemEntry does on
                    # the JSON wire. Without it a cache's request interface loads
                    # as a plain Entry and the cache has nothing to serve.
                    ie_elem = entry_elem.find('./item-entry')
                    if ie_elem is not None:
                        cardinality = int(float(ie_elem.get('cardinality', '1')))
                        pop_elem = ie_elem.find('./access-popularity')
                        if pop_elem is not None:
                            params = [float(p.get('value', '0'))
                                      for p in pop_elem.findall('./parameter')]
                            popularity = _popularity_from_params(
                                pop_elem.get('name', 'DiscreteSampler'), params, cardinality)
                        else:
                            # linemodel_load.m defaults an absent popularity to uniform
                            from .distributions import DiscreteSampler as _DS
                            popularity = _DS(np.ones(cardinality) / cardinality)
                        entry = ItemEntry(model, entry_name, cardinality, popularity)
                        entry.on(task)
                    else:
                        entry = model.add_entry(entry_name, task)
                    entry_map[entry_name] = entry

                    # Parse open-arrival-rate (mirrors JAR LayeredNetwork.java:357-360)
                    open_arrival_str = entry_elem.get('open-arrival-rate', '')
                    if open_arrival_str:
                        try:
                            open_arrival_rate = float(open_arrival_str)
                            if open_arrival_rate > 0:
                                entry.setArrival(Exp.fitMean(1.0 / open_arrival_rate))
                        except ValueError:
                            pass

                    # Parse forwarding calls (must be deferred until all entries are created)
                    for fwd_elem in entry_elem.findall('./forwarding'):
                        dest_name = fwd_elem.get('dest', '')
                        prob_str = fwd_elem.get('prob', '1.0')
                        try:
                            prob = float(prob_str)
                        except ValueError:
                            prob = 1.0
                        if dest_name:
                            if not hasattr(entry, '_pending_forwards'):
                                entry._pending_forwards = []
                            entry._pending_forwards.append((dest_name, prob))

                    # Handle entry-phase-activities (phase-based entries)
                    # Convert to activity-graph format for consistency
                    epa_elem = entry_elem.find('./entry-phase-activities')
                    if epa_elem is not None:
                        phase_activities = []
                        for act_elem in epa_elem.findall('./activity'):
                            act_name = act_elem.get('name', '')
                            demand_mean = float(act_elem.get('host-demand-mean', '0'))
                            phase = int(act_elem.get('phase', '1'))

                            demand = _dist_from_mean_scv(demand_mean, 1.0)

                            activity = model.add_activity(act_name, demand, task)
                            activity.phase = phase  # Store phase number (1 or 2)
                            activity_map[act_name] = activity
                            phase_activities.append((phase, activity))

                            # Parse activity think-time (LQNX think-time attribute)
                            act_think_time = act_elem.get('think-time', '')
                            if act_think_time:
                                activity.think_time = float(act_think_time)

                            # Parse calls within phase activity
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

                            _parse_call_groups(act_elem, activity)

                        # Sort by phase and set up binding/reply
                        phase_activities.sort(key=lambda x: x[0])
                        if phase_activities:
                            # Phase 1 activity is bound to entry (matches MATLAB parseXML line 223-224)
                            first_activity = phase_activities[0][1]
                            first_activity.bound_to(entry)

                            # Last PHASE 1 activity replies to entry (not last overall)
                            # MATLAB parseXML line 276: newEntry.replyActivity{end+1} = name{1}
                            # Phase 2 activities are post-reply processing
                            if task_sched_enum != SchedStrategy.REF:
                                # Find last phase 1 activity
                                last_ph1 = None
                                for ph, act in phase_activities:
                                    if ph == 1:
                                        last_ph1 = act
                                if last_ph1 is not None:
                                    last_ph1.replies_to(entry)

                            # Create serial precedence if multiple phases
                            if len(phase_activities) > 1:
                                acts = [pa[1] for pa in phase_activities]
                                task.add_precedence(ActivityPrecedence.Serial(acts))

                for ta_elem in task_elem.findall('./task-activities'):
                    for act_elem in ta_elem.findall('./activity'):
                        act_name = act_elem.get('name', '')
                        demand_mean = float(act_elem.get('host-demand-mean', '0'))
                        demand_scv = float(act_elem.get('host-demand-cvsq', '1.0'))

                        demand = _dist_from_mean_scv(demand_mean, demand_scv)

                        activity = model.add_activity(act_name, demand, task)
                        activity_map[act_name] = activity

                        # Parse activity think-time (LQNX think-time attribute)
                        act_think_time = act_elem.get('think-time', '')
                        if act_think_time:
                            activity.think_time = float(act_think_time)

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

                        _parse_call_groups(act_elem, activity)

                    for prec_elem in ta_elem.findall('./precedence'):
                        pre_acts = []
                        post_acts = []
                        pre_type = 'pre'
                        post_type = 'post'

                        pre_quorum = None
                        for pre_tag in ['pre', 'pre-AND', 'pre-OR']:
                            pre_elem = prec_elem.find(f'./{pre_tag}')
                            if pre_elem is not None:
                                pre_type = pre_tag
                                # An AND-join may declare a quorum k: it fires once k of its
                                # predecessors complete, rather than waiting for all of them.
                                quorum_str = pre_elem.get('quorum', '')
                                if quorum_str:
                                    try:
                                        pre_quorum = int(float(quorum_str))
                                    except ValueError:
                                        pre_quorum = None
                                for act_ref in pre_elem.findall('./activity'):
                                    act_name = act_ref.get('name', '')
                                    if act_name in activity_map:
                                        pre_acts.append(activity_map[act_name])
                                break

                        # The post side is minOccurs="0" in lqn-core.xsd: a precedence
                        # carrying only a pre element declares a TERMINAL activity and no
                        # successor, so it contributes no edge.
                        post_elem = None
                        for post_tag in ['post', 'post-AND', 'post-OR', 'post-LOOP',
                                         'post-CACHE']:
                            post_elem = prec_elem.find(f'./{post_tag}')
                            if post_elem is not None:
                                post_type = post_tag
                                for act_ref in post_elem.findall('./activity'):
                                    act_name = act_ref.get('name', '')
                                    if act_name in activity_map:
                                        post_acts.append(activity_map[act_name])
                                break
                        if post_elem is None:
                            continue

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
                            elif post_type == 'post-CACHE':
                                # cache-result is explicit; a file written before
                                # the attribute existed falls back to document
                                # order, which is hit first and miss second.
                                ordered = {}
                                for i, act_ref in enumerate(post_elem.findall('./activity')):
                                    act_name = act_ref.get('name', '')
                                    if act_name not in activity_map:
                                        continue
                                    result = (act_ref.get('cache-result')
                                              or ('hit' if i == 0 else 'miss')).lower()
                                    ordered.setdefault(result, activity_map[act_name])
                                outcome = [a for a in (ordered.get('hit'), ordered.get('miss'))
                                           if a is not None] or post_acts
                                prec = ActivityPrecedence.CacheAccess(pre_acts[0], outcome)
                            elif post_type == 'post-LOOP':
                                # Parse LOOP: post_acts = loop body activities
                                # end attribute = name of end activity
                                # count attribute on <activity> = mean loop iterations
                                loop_count = 1.0
                                for act_ref in post_elem.findall('./activity'):
                                    count_str = act_ref.get('count', '')
                                    if count_str:
                                        loop_count = float(count_str)
                                        break
                                end_act_name = post_elem.get('end', '')
                                if end_act_name and end_act_name in activity_map:
                                    # Loop(pre_act, [body_acts..., end_act], count)
                                    loop_acts = list(post_acts) + [activity_map[end_act_name]]
                                    prec = ActivityPrecedence.Loop(pre_acts[0], loop_acts, loop_count)
                                else:
                                    prec = ActivityPrecedence.Serial(pre_acts + post_acts)
                            elif pre_type == 'pre-AND':
                                prec = ActivityPrecedence.AndJoin(pre_acts, post_acts[0], pre_quorum)
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
            if hasattr(activity, '_pending_call_groups'):
                for strategy, dests in activity._pending_call_groups:
                    entries = [entry_map[d] for d in dests if d in entry_map]
                    if len(entries) >= 2:
                        activity.record_call_group(strategy, entries)
                delattr(activity, '_pending_call_groups')

        # Resolve pending forwarding calls (must be after all entries created)
        for entry in model.entries:
            if hasattr(entry, '_pending_forwards'):
                for dest_name, prob in entry._pending_forwards:
                    if dest_name in entry_map:
                        entry.addForwarding(entry_map[dest_name], prob)
                delattr(entry, '_pending_forwards')

        return model

    @staticmethod
    def _validate_input_model(root) -> None:
        """
        Reject a structurally inconsistent LQN document.

        Run on the parsed document before any object is built, so that a defective
        input is named at its source instead of surfacing as a downstream failure.
        The same checks, in the same order and with the same messages, are applied
        by the MATLAB, JAR and C++ readers.

        Args:
            root: root element of the parsed LQN document

        Raises:
            RuntimeError: on the first inconsistency found
        """
        from .api.io.logging import line_error

        def num(text, dflt):
            if text is None or text == '':
                return dflt
            try:
                return float(text)
            except ValueError:
                return float('nan')

        tol = 1e-6
        proc_names = []
        task_names = []
        entry_names = []
        entry_owner = []  # task owning entry_names[k]
        is_ref_entry = []
        call_dests = []
        reply_entries = []
        has_ref_task = False
        has_open_arrival = False

        for proc_elem in root.iter('processor'):
            proc_name = proc_elem.get('name', '')
            if proc_name in proc_names:
                line_error('parse_xml', f'Duplicate processor name "{proc_name}".')
            proc_names.append(proc_name)

            for task_elem in proc_elem.iter('task'):
                task_name = task_elem.get('name', '')
                if task_name in task_names:
                    line_error('parse_xml', f'Duplicate task name "{task_name}".')
                task_names.append(task_name)
                is_ref = (task_elem.get('scheduling', '') or '').lower() == 'ref'
                has_ref_task = has_ref_task or is_ref

                entry_elems = list(task_elem.iter('entry'))
                if not entry_elems:
                    line_error('parse_xml', f'Task "{task_name}" has no entries.')
                for entry_elem in entry_elems:
                    entry_name = entry_elem.get('name', '')
                    if entry_name in entry_names:
                        line_error('parse_xml', f'Duplicate entry name "{entry_name}".')
                    entry_names.append(entry_name)
                    entry_owner.append(task_name)
                    is_ref_entry.append(is_ref)

                    open_arrival_rate = num(entry_elem.get('open-arrival-rate'), float('nan'))
                    if open_arrival_rate > 0:
                        has_open_arrival = True
                        if is_ref:
                            line_error('parse_xml', f'Entry "{entry_name}" belongs to reference task "{task_name}" and cannot have open arrivals.')

                    fwd_elems = list(entry_elem.iter('forwarding'))
                    if is_ref and fwd_elems:
                        line_error('parse_xml', f'Entry "{entry_name}" belongs to reference task "{task_name}" and cannot forward requests.')
                    fwd_total = 0.0
                    for fwd_elem in fwd_elems:
                        prob = num(fwd_elem.get('prob'), 1.0)
                        if prob != prob or prob < 0.0 or prob > 1.0:
                            line_error('parse_xml', f'Forwarding from entry "{entry_name}" to entry "{fwd_elem.get("dest", "")}" has an invalid probability of {prob:g}.')
                        fwd_total += prob
                    if fwd_total > 1.0 + tol:
                        line_error('parse_xml', f'Entry "{entry_name}" has a total forwarding probability of {fwd_total:g}.')

                # activity names are unique within their task; a name under a pre or post list is a reference, not a declaration
                act_names = []
                parent_of = {}
                for parent in task_elem.iter():
                    for child in parent:
                        parent_of[child] = parent
                for act_elem in task_elem.iter('activity'):
                    parent = parent_of.get(act_elem)
                    if parent is None or parent.tag not in ('task-activities', 'entry-phase-activities'):
                        continue
                    act_name = act_elem.get('name', '')
                    if act_name in act_names:
                        line_error('parse_xml', f'Duplicate activity name "{act_name}" in task "{task_name}".')
                    act_names.append(act_name)

                for call_elem in task_elem.iter('synch-call'):
                    call_dests.append(call_elem.get('dest', ''))
                for call_elem in task_elem.iter('asynch-call'):
                    call_dests.append(call_elem.get('dest', ''))
                for fwd_elem in task_elem.iter('forwarding'):
                    call_dests.append(fwd_elem.get('dest', ''))

                for or_elem in task_elem.iter('post-OR'):
                    branch_total = 0.0
                    for branch_elem in or_elem.iter('activity'):
                        prob = num(branch_elem.get('prob'), 1.0)
                        if prob != prob or prob < 0.0 or prob > 1.0:
                            line_error('parse_xml', f'Activity "{branch_elem.get("name", "")}" in task "{task_name}" has an invalid branch probability of {prob:g}.')
                        branch_total += prob
                    if abs(branch_total - 1.0) > tol:
                        line_error('parse_xml', f'Branch probabilities of an OR-fork in task "{task_name}" sum to {branch_total:g} instead of 1.')

                for reply_elem in task_elem.iter('reply-entry'):
                    reply_entries.append(reply_elem.get('name', ''))

        for dest in call_dests:
            if dest in entry_names:
                idx = entry_names.index(dest)
                if is_ref_entry[idx]:
                    line_error('parse_xml', f'Entry "{entry_names[idx]}" belongs to reference task "{entry_owner[idx]}" and cannot receive requests.')

        for reply_name in reply_entries:
            if reply_name in entry_names:
                idx = entry_names.index(reply_name)
                if is_ref_entry[idx]:
                    line_error('parse_xml', f'Entry "{entry_names[idx]}" belongs to reference task "{entry_owner[idx]}" and cannot be replied to.')

        if not has_ref_task and not has_open_arrival:
            line_error('parse_xml', 'The model has no reference task and no open arrivals.')

    parseXML = parse_xml
    readXML = parse_xml
    load = parse_xml

    @staticmethod
    def _parse_replacement(name: str) -> ReplacementStrategy:
        """Wire spelling -> ReplacementStrategy, the same names the JSON
        interchange uses (linemodel_save.m `repl_to_str`). An unknown spelling
        falls back to FIFO, as `linemodel_load.m` does, rather than raising."""
        try:
            return getattr(ReplacementStrategy, str(name).upper())
        except AttributeError:
            return ReplacementStrategy.FIFO

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
        elif sched_upper in ('PRI', 'PP'):
            # LQNS preemptive priority resume (SCHEDULE_PPR). lqns spells it
            # 'pri' (LQIO::SCHEDULE::PPR); 'pp' is the stale lqn-core.xsd
            # spelling, absent from the lqns 6.2.31 sources
            return SchedStrategy.FCFSPRPRIO
        elif sched_upper in ('LCFS', 'LIFO'):
            return SchedStrategy.LCFS
        else:
            return SchedStrategy.FCFS


class SetupTask(Task):
    """
    Task whose servers are switched off while idle.

    A server resuming from the off state pays a setup (activation) time before
    serving the request that woke it up, and stays available for a delay-off
    (idle) period after emptying its queue before switching off again:
    - setup_time: activation time paid on resuming from the off state
    - delay_off_time: idle period before a server powers off

    These are the setup and close-down times of a server with vacations, e.g.
    on-demand virtual machines and containers, power-managed servers under a
    timeout policy, warm-up delays, serverless cold start / keep-alive. Both
    times are declared on the base Task, so this subclass is a naming
    convenience.
    """

    def __init__(self, model_or_name, name_or_mult=None, mult_or_sched=None, sched=None):
        """Initialize a SetupTask with flexible arguments."""
        super().__init__(model_or_name, name_or_mult, mult_or_sched, sched)
        self._is_setup_task = True

    def has_setup_delayoff(self) -> bool:
        """Return True to indicate this is a SetupTask."""
        return True


class FunctionTask(SetupTask):
    """
    Former name of SetupTask, kept for backward compatibility.

    Setup and delay-off times are not specific to serverless
    (function-as-a-service) platforms, so the class carrying them is now named
    after the modelling primitive rather than after that application domain.
    """
    pass


# Convenience aliases for compatibility
Processor = Processor
Task = Task
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
    'SetupTask',
    'FunctionTask',
    'Entry',
    'Activity',
    'CacheTask',
    'ItemEntry',
]
