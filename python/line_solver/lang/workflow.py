"""
Computational Workflow for Phase-Type Distribution Conversion.

This module provides classes for modeling computational workflows that can
be converted to phase-type (APH) distributions for queueing analysis.

Key Classes:
    Workflow: A computational workflow with activities and precedences
    WorkflowActivity: A computational activity in a workflow
    ActivityPrecedence: Precedence relationship between activities

The Workflow class supports:
- Serial, parallel, loop, and branching structures
- Conversion to phase-type distributions via toPH()
- Loading from WfCommons JSON format

A workflow whose precedence graph is series-parallel is reduced exactly, by
recursive composition of the series-parallel tree, which handles arbitrary
nesting (a fork inside a loop, a branch that is itself a fork-join). Graphs
that are not series-parallel fall back to the block-based composition, which is
a heuristic. A loop repeats its body a geometric number of times of mean COUNT,
the semantics of the POST_LOOP precedence of an activity graph.

References:
    Original Java: jar/src/main/java/jline/lang/workflow/Workflow.java
"""

import numpy as np
from typing import List, Dict, Any, Optional, Tuple, Union

from ..constants import GlobalConstants
from dataclasses import dataclass, field
from enum import Enum
import json


class ActivityPrecedenceType:
    """Constants for activity precedence types."""
    PRE_SEQ = "pre"
    PRE_AND = "pre-AND"
    PRE_OR = "pre-OR"
    POST_SEQ = "post"
    POST_AND = "post-AND"
    POST_OR = "post-OR"
    POST_LOOP = "post-LOOP"
    POST_CACHE = "post-CACHE"

    # Numeric IDs
    ID_PRE_SEQ = 1
    ID_PRE_AND = 2
    ID_PRE_OR = 3
    ID_POST_SEQ = 11
    ID_POST_AND = 12
    ID_POST_OR = 13
    ID_POST_LOOP = 14
    ID_POST_CACHE = 15


@dataclass
class ActivityPrecedence:
    """
    A precedence relationship between workflow activities.

    Supports various precedence types:
    - Serial: simple sequence (pre -> post)
    - AndFork: parallel split (pre -> all post in parallel)
    - AndJoin: synchronization (all pre -> post)
    - OrFork: probabilistic branch (pre -> one of post with probability)
    - OrJoin: merge (first of pre -> post)
    - Loop: iteration (pre -> loop activities -> end)
    """
    pre_acts: List[str]
    post_acts: List[str]
    pre_type: str = ActivityPrecedenceType.PRE_SEQ
    post_type: str = ActivityPrecedenceType.POST_SEQ
    pre_params: Optional[np.ndarray] = None
    post_params: Optional[np.ndarray] = None

    @staticmethod
    def Serial(pre_act: Union[str, 'WorkflowActivity'],
               post_act: Union[str, 'WorkflowActivity']) -> 'ActivityPrecedence':
        """Create a serial precedence: pre_act -> post_act."""
        pre_name = pre_act if isinstance(pre_act, str) else pre_act.name
        post_name = post_act if isinstance(post_act, str) else post_act.name
        return ActivityPrecedence(
            pre_acts=[pre_name],
            post_acts=[post_name],
            pre_type=ActivityPrecedenceType.PRE_SEQ,
            post_type=ActivityPrecedenceType.POST_SEQ
        )

    @staticmethod
    def SerialSequence(*activities) -> List['ActivityPrecedence']:
        """Create serial precedences for a sequence of activities."""
        precedences = []
        for i in range(len(activities) - 1):
            precedences.append(ActivityPrecedence.Serial(activities[i], activities[i + 1]))
        return precedences

    @staticmethod
    def AndFork(pre_act: Union[str, 'WorkflowActivity'],
                post_acts: List[Union[str, 'WorkflowActivity']],
                fanout: Optional[np.ndarray] = None) -> 'ActivityPrecedence':
        """Create an AND-fork: pre_act -> all post_acts in parallel."""
        pre_name = pre_act if isinstance(pre_act, str) else pre_act.name
        post_names = [a if isinstance(a, str) else a.name for a in post_acts]
        if fanout is None:
            fanout = np.ones(len(post_acts))
        return ActivityPrecedence(
            pre_acts=[pre_name],
            post_acts=post_names,
            pre_type=ActivityPrecedenceType.PRE_SEQ,
            post_type=ActivityPrecedenceType.POST_AND,
            post_params=fanout
        )

    @staticmethod
    def AndJoin(pre_acts: List[Union[str, 'WorkflowActivity']],
                post_act: Union[str, 'WorkflowActivity'],
                quorum: Optional[np.ndarray] = None) -> 'ActivityPrecedence':
        """Create an AND-join: all pre_acts -> post_act (synchronization)."""
        pre_names = [a if isinstance(a, str) else a.name for a in pre_acts]
        post_name = post_act if isinstance(post_act, str) else post_act.name
        # An absent quorum means a full join, as in layered.ActivityPrecedence
        return ActivityPrecedence(
            pre_acts=pre_names,
            post_acts=[post_name],
            pre_type=ActivityPrecedenceType.PRE_AND,
            post_type=ActivityPrecedenceType.POST_SEQ,
            pre_params=None if quorum is None else np.atleast_1d(np.asarray(quorum, dtype=float))
        )

    @staticmethod
    def OrFork(pre_act: Union[str, 'WorkflowActivity'],
               post_acts: List[Union[str, 'WorkflowActivity']],
               probs: np.ndarray) -> 'ActivityPrecedence':
        """Create an OR-fork: pre_act -> one of post_acts with probability."""
        pre_name = pre_act if isinstance(pre_act, str) else pre_act.name
        post_names = [a if isinstance(a, str) else a.name for a in post_acts]
        return ActivityPrecedence(
            pre_acts=[pre_name],
            post_acts=post_names,
            pre_type=ActivityPrecedenceType.PRE_SEQ,
            post_type=ActivityPrecedenceType.POST_OR,
            post_params=np.asarray(probs)
        )

    @staticmethod
    def OrJoin(pre_acts: List[Union[str, 'WorkflowActivity']],
               post_act: Union[str, 'WorkflowActivity']) -> 'ActivityPrecedence':
        """Create an OR-join: first of pre_acts -> post_act."""
        pre_names = [a if isinstance(a, str) else a.name for a in pre_acts]
        post_name = post_act if isinstance(post_act, str) else post_act.name
        return ActivityPrecedence(
            pre_acts=pre_names,
            post_acts=[post_name],
            pre_type=ActivityPrecedenceType.PRE_OR,
            post_type=ActivityPrecedenceType.POST_SEQ
        )

    @staticmethod
    def Loop(pre_act: Union[str, 'WorkflowActivity'],
             loop_acts: List[Union[str, 'WorkflowActivity']],
             end_act: Optional[Union[str, 'WorkflowActivity', float]] = None,
             count: float = 1.0) -> 'ActivityPrecedence':
        """
        Create a loop: pre_act -> loop_acts (repeated count times) -> end_act.

        The body repeats a geometric number of times of mean COUNT, the
        POST_LOOP semantics of an activity graph. Both call conventions are
        accepted: Loop(pre, [body, end], count), as in MATLAB and the JAR and
        in layered.ActivityPrecedence, and Loop(pre, [body], end, count).
        """
        if isinstance(end_act, (int, float)) and not isinstance(end_act, bool):
            # Loop(pre, [body..., end], count): the third positional is COUNT
            count = float(end_act)
            end_act = None
        pre_name = pre_act if isinstance(pre_act, str) else pre_act.name
        post_names = [a if isinstance(a, str) else a.name for a in loop_acts]
        if end_act is not None:
            end_name = end_act if isinstance(end_act, str) else end_act.name
            post_names.append(end_name)
        return ActivityPrecedence(
            pre_acts=[pre_name],
            post_acts=post_names,
            pre_type=ActivityPrecedenceType.PRE_SEQ,
            post_type=ActivityPrecedenceType.POST_LOOP,
            post_params=np.array([count])
        )



def _scale_distribution_rate(distrib, factor: float):
    """
    Time-scale a distribution by FACTOR, or None when it cannot be scaled.

    The import is deferred because line_solver.distributions imports this
    module's package, so a module-level import would close a cycle.
    """
    try:
        from ..distributions import dist_scale_rate
    except ImportError:
        return None
    try:
        return dist_scale_rate(distrib, factor)
    except (ValueError, TypeError, NotImplementedError):
        return None


class WorkflowActivity:
    """
    A computational activity in a Workflow.

    Represents a single activity with a host demand (service time distribution).
    Activities can be composed into workflows using precedence relationships.

    Unlike Activity in LayeredNetwork it carries no call list: an external call
    is represented as an activity whose host demand is the law of the call
    response time, so that a synchronous call and a local computation compose in
    the same way. An asynchronous call blocks the caller for no time and is
    simply left out of the workflow.

    Args:
        workflow: Parent Workflow object
        name: Activity name
        host_demand: Mean service time or Distribution object

    Example:
        >>> wf = Workflow("MyWorkflow")
        >>> a = wf.addActivity("A", 1.0)  # Exp(1.0) service
        >>> b = wf.addActivity("B", 2.0)
    """

    def __init__(self, workflow: 'Workflow', name: str,
                 host_demand: Union[float, Any]):
        self._workflow = workflow
        self._name = name
        self._index = -1
        self._metadata: Dict[str, Any] = {}
        self._distribution = None
        self._host_demand_mean = 1.0
        self._host_demand_scv = 1.0
        self.setHostDemand(host_demand)

    @property
    def name(self) -> str:
        return self._name

    @property
    def index(self) -> int:
        return self._index

    @index.setter
    def index(self, value: int):
        self._index = value

    @property
    def host_demand_mean(self) -> float:
        return self._host_demand_mean

    @property
    def host_demand_scv(self) -> float:
        return self._host_demand_scv

    @property
    def metadata(self) -> Dict[str, Any]:
        return self._metadata

    @metadata.setter
    def metadata(self, value: Dict[str, Any]):
        self._metadata = value

    def setHostDemand(self, value: Union[float, Any]) -> None:
        """Set the host demand (service time)."""
        if isinstance(value, (int, float)):
            if float(value) <= GlobalConstants.FineTol:
                # Zero-time activity, as in the MATLAB and JAR twins
                self._host_demand_mean = GlobalConstants.FineTol
                self._host_demand_scv = GlobalConstants.FineTol
                self._distribution = None
            else:
                self._host_demand_mean = float(value)
                self._host_demand_scv = 1.0
                self._distribution = None
        else:
            self._distribution = value
            self._host_demand_mean = value.getMean() if hasattr(value, 'getMean') else 1.0
            self._host_demand_scv = value.getSCV() if hasattr(value, 'getSCV') else 1.0
        self._invalidateParent()

    def setHostDemandMean(self, mean_value: float) -> None:
        """
        Change the mean, preserving the shape.

        Scales the current law in time rather than refitting it, so the SCV,
        the skewness and the order are preserved and the cached series-parallel
        tree keeps its shape.
        """
        if not (mean_value > 0) or np.isinf(mean_value):
            raise ValueError("The activity mean must be a positive finite scalar.")

        old_mean = self._host_demand_mean
        if self._distribution is None or not (old_mean > 0) or np.isinf(old_mean) \
                or old_mean <= GlobalConstants.FineTol:
            self.setHostDemand(float(mean_value))
            return

        factor = old_mean / mean_value
        scaled = _scale_distribution_rate(self._distribution, factor)
        if scaled is None:
            self.setHostDemand(float(mean_value))
            return
        self._distribution = scaled
        self._host_demand_mean = float(mean_value)
        # The SCV is invariant under a time scaling
        if self._workflow is not None and self._index >= 0:
            self._workflow.rescaleActivityLeaf(self._index, factor)

    def _invalidateParent(self) -> None:
        # The parent caches the composed law, so the leaf must be marked dirty
        # here as well as on a topology change
        if self._workflow is not None and self._index >= 0:
            self._workflow.invalidateActivity(self._index)

    def getPHRepresentation(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Get phase-type representation of this activity.

        Returns:
            Tuple of (alpha, T) where:
                alpha: Initial probability vector (1 x n)
                T: Sub-generator matrix (n x n)
        """
        # Handle immediate (zero service time), on the same threshold as the
        # MATLAB and JAR twins
        if self._host_demand_mean <= GlobalConstants.FineTol:
            return np.array([[1.0]]), np.array([[-GlobalConstants.Immediate]])

        # Handle Markovian distribution
        if self._distribution is not None and hasattr(self._distribution, 'getInitProb'):
            alpha = self._distribution.getInitProb()
            # Try different methods for getting the sub-generator matrix
            T = None
            if hasattr(self._distribution, 'getD0'):
                T = self._distribution.getD0()
            elif hasattr(self._distribution, 'D'):
                T = self._distribution.D(0)
            if T is not None:
                if alpha.ndim == 1:
                    alpha = alpha.reshape(1, -1)
                elif alpha.shape[0] > 1 and alpha.shape[1] == 1:
                    alpha = alpha.T
                return alpha, T

        # Default: exponential distribution
        rate = 1.0 / self._host_demand_mean
        alpha = np.array([[1.0]])
        T = np.array([[-rate]])
        return alpha, T

    def getNumberOfPhases(self) -> int:
        """Get number of phases in the PH representation."""
        _, T = self.getPHRepresentation()
        return T.shape[0]

    def getName(self) -> str:
        return self._name

    def getHostDemandMean(self) -> float:
        return self._host_demand_mean

    def getHostDemandSCV(self) -> float:
        return self._host_demand_scv

    def __repr__(self) -> str:
        return f"WorkflowActivity('{self._name}', mean={self._host_demand_mean:.4f})"


class Workflow:
    """
    A computational workflow that can be converted to a phase-type distribution.

    Workflow models a directed acyclic graph of activities with precedence
    relationships. The workflow can be converted to an APH (Acyclic Phase-Type)
    distribution for use in queueing network analysis.

    Supports:
    - Serial composition (sequence of activities)
    - Parallel composition (AND-fork/join)
    - Probabilistic branching (OR-fork/join)
    - Loops whose body repeats a geometric number of times

    A precedence graph that is series-parallel is reduced exactly, by recursive
    composition of the series-parallel tree, which handles arbitrary nesting.
    Graphs that are not series-parallel fall back to the block composition,
    which is a heuristic.

    Example:
        >>> wf = Workflow("ServiceWorkflow")
        >>> a = wf.addActivity("A", 1.0)
        >>> b = wf.addActivity("B", 2.0)
        >>> c = wf.addActivity("C", 0.5)
        >>> wf.addPrecedence(ActivityPrecedence.Serial(a, b))
        >>> wf.addPrecedence(ActivityPrecedence.Serial(b, c))
        >>> alpha, T = wf.toPH()

    References:
        Original Java: jar/src/main/java/jline/lang/workflow/Workflow.java
    """

    # Class-level aliases for precedence constructors
    Serial = staticmethod(ActivityPrecedence.Serial)
    SerialSequence = staticmethod(ActivityPrecedence.SerialSequence)
    AndFork = staticmethod(ActivityPrecedence.AndFork)
    AndJoin = staticmethod(ActivityPrecedence.AndJoin)
    OrFork = staticmethod(ActivityPrecedence.OrFork)
    OrJoin = staticmethod(ActivityPrecedence.OrJoin)
    Loop = staticmethod(ActivityPrecedence.Loop)

    def __init__(self, name: str):
        self._name = name
        self._activities: List[WorkflowActivity] = []
        self._activity_map: Dict[str, int] = {}
        self._precedences: List[ActivityPrecedence] = []
        self._cached_ph: Optional[Tuple[np.ndarray, np.ndarray]] = None
        self._sp_tree: Optional[Dict[str, Any]] = None
        self._sp_failed: bool = False

    @property
    def name(self) -> str:
        return self._name

    def addActivity(self, name: str,
                    host_demand: Union[float, Any] = 1.0) -> WorkflowActivity:
        """
        Add an activity to the workflow.

        Args:
            name: Activity name
            host_demand: Mean service time or Distribution object

        Returns:
            The created WorkflowActivity
        """
        act = WorkflowActivity(self, name, host_demand)
        self._activities.append(act)
        act.index = len(self._activities) - 1
        self._activity_map[name] = act.index
        self.invalidateTopology()
        return act

    # Snake_case alias
    add_activity = addActivity

    def addPrecedence(self, prec: Union[ActivityPrecedence, List[ActivityPrecedence]]) -> None:
        """
        Add precedence relationship(s) to the workflow.

        Args:
            prec: ActivityPrecedence or list of ActivityPrecedence objects
        """
        if isinstance(prec, list):
            for p in prec:
                self._precedences.append(p)
        else:
            self._precedences.append(prec)
        self.invalidateTopology()

    # Snake_case alias
    add_precedence = addPrecedence

    def getActivity(self, name: str) -> Optional[WorkflowActivity]:
        """Get an activity by name."""
        idx = self._activity_map.get(name)
        if idx is None:
            return None
        return self._activities[idx]

    def getActivityIndex(self, name: str) -> int:
        """Get the zero-based index of an activity, or -1 if absent."""
        idx = self._activity_map.get(name)
        return -1 if idx is None else idx

    def getActivities(self) -> List[WorkflowActivity]:
        """Get all activities."""
        return list(self._activities)

    def getPrecedences(self) -> List[ActivityPrecedence]:
        """Get all precedences."""
        return list(self._precedences)

    def validate(self) -> Tuple[bool, str]:
        """
        Validate the workflow structure.

        Returns:
            Tuple of (is_valid, error_message)
        """
        if not self._activities:
            return False, "Workflow must have at least one activity."

        # Check all referenced activities exist
        for prec in self._precedences:
            for act_name in prec.pre_acts:
                if act_name not in self._activity_map:
                    return False, f"Activity '{act_name}' referenced in precedence not found."
            for act_name in prec.post_acts:
                if act_name not in self._activity_map:
                    return False, f"Activity '{act_name}' referenced in precedence not found."

        # Check OR-fork probabilities
        for prec in self._precedences:
            if prec.post_type == ActivityPrecedenceType.POST_OR:
                if prec.post_params is None:
                    return False, "OR-fork must have probabilities."
                prob_sum = np.sum(prec.post_params)
                if abs(prob_sum - 1.0) > GlobalConstants.FineTol:
                    return False, ("OR-fork probabilities must sum to 1 "
                                   "(got %.4f)." % prob_sum)

        # Check loop counts are a single positive number
        for prec in self._precedences:
            if prec.post_type == ActivityPrecedenceType.POST_LOOP:
                if prec.post_params is None or np.size(prec.post_params) != 1:
                    return False, "Loop count must be a single positive number."
                if float(np.ravel(prec.post_params)[0]) <= 0:
                    return False, "Loop count must be a positive number."

        # Partial (quorum) AND-joins are refused rather than silently served as
        # full joins, which would be a different law
        for prec in self._precedences:
            if prec.pre_type == ActivityPrecedenceType.PRE_AND and prec.pre_params is not None:
                params = np.ravel(np.asarray(prec.pre_params, dtype=float))
                if params.size == 0:
                    continue
                quorum = params[0]
                if 0 < quorum < len(prec.pre_acts):
                    return False, (
                        "AND-join with quorum %d of %d is not supported by Workflow: "
                        "a partial join is not the maximum of the branches. Use a full "
                        "join, or SolverLN with method='default', which routes the join "
                        "explicitly." % (int(quorum), len(prec.pre_acts)))

        return True, ""

    def toPH(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Convert workflow to a phase-type representation.

        A series-parallel precedence graph is reduced exactly; any other graph
        falls back to the block composition. The generator is acyclic (an APH)
        unless a geometric loop closes a cycle over a multi-phase body, which
        isAcyclicGenerator reports.

        Returns:
            Tuple of (alpha, T) where:
                alpha: Initial probability vector (1 x n)
                T: Sub-generator matrix (n x n)

        Raises:
            ValueError: If workflow is invalid
        """
        if self._cached_ph is not None:
            return self._cached_ph

        is_valid, error = self.validate()
        if not is_valid:
            raise ValueError(error)

        result = self._composeSeriesParallel()
        if result is None:
            # Not series-parallel: fall back to the block composition
            result = self._buildCTMC()

        self._cached_ph = result
        return result

    def refreshPH(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Recompose the workflow law after a demand change.

        Recomposes only the series-parallel nodes on the path from a dirty leaf
        to the root; nodes whose subtree is unchanged keep their cached
        (alpha, T). The topology is not rebuilt. When the workflow is not
        series-parallel this degrades to a full recomposition through the block
        path.
        """
        return self.toPH()

    def setActivityDemand(self, name: str, host_demand: Union[float, Any]) -> 'Workflow':
        """
        Change the host demand of one activity.

        Marks only that leaf dirty, so that the next toPH/refreshPH recomposes
        the path from the leaf to the root and reuses every other cached block.
        This is the entry point used by iterative solvers that update
        call-response laws at each iteration.
        """
        act = self.getActivity(name)
        if act is None:
            raise ValueError("Activity '%s' not found in workflow." % name)
        act.setHostDemand(host_demand)
        return self

    def setActivityDemandMean(self, name: str, mean_value: float) -> 'Workflow':
        """
        Change only the mean of one activity.

        Rescales the activity law in time, so its SCV and its whole shape are
        preserved. The leaf keeps its order and its initial probability vector,
        and the cached series-parallel tree keeps its shape.
        """
        act = self.getActivity(name)
        if act is None:
            raise ValueError("Activity '%s' not found in workflow." % name)
        act.setHostDemandMean(mean_value)
        return self

    def invalidateTopology(self) -> None:
        """
        Discard the cached law and decomposition.

        Called when an activity or a precedence is added, which can change the
        shape of the series-parallel tree.
        """
        self._cached_ph = None
        self._sp_tree = None
        self._sp_failed = False

    def invalidateActivity(self, act_idx: int) -> None:
        """
        Mark one activity law as dirty.

        Invalidates the leaf of ACT_IDX and its ancestors in the cached
        series-parallel tree, keeping every other cached block. The topology is
        untouched.
        """
        self._cached_ph = None
        if self._sp_tree is None:
            return
        leaf_of = self._sp_tree['leaf_of']
        if act_idx < 0 or act_idx >= len(leaf_of) or leaf_of[act_idx] < 0:
            # Activity outside the decomposition: rebuild it entirely
            self._sp_tree = None
            self._sp_failed = False
            return
        Workflow._invalidateBranch(self._sp_tree, leaf_of[act_idx])

    def rescaleActivityLeaf(self, act_idx: int, factor: float) -> None:
        """
        Time-scale a cached leaf in place.

        The leaf law is scaled as T -> T*FACTOR with ALPHA fixed, which divides
        its mean by FACTOR and leaves its SCV and its order untouched.
        Ancestors still recompose, because they mix phases of several leaves,
        but the tree keeps its shape and no APH is refitted.
        """
        self._cached_ph = None
        if self._sp_tree is None:
            return
        leaf_of = self._sp_tree['leaf_of']
        if act_idx < 0 or act_idx >= len(leaf_of) or leaf_of[act_idx] < 0:
            self._sp_tree = None
            self._sp_failed = False
            return
        k = leaf_of[act_idx]
        Workflow._invalidateBranch(self._sp_tree, k)
        if self._sp_tree['T'][k] is not None:
            self._sp_tree['T'][k] = self._sp_tree['T'][k] * factor
            self._sp_tree['valid'][k] = True

    def getSPTree(self) -> Optional[Dict[str, Any]]:
        """
        Return the cached series-parallel decomposition.

        Returns the flat series-parallel tree, or None if the precedence graph
        is not series-parallel. Key 'execs' carries the expected number of
        executions of each node per workflow execution, which is what an LQN
        metric reconstruction splits the layer results by.
        """
        if self._sp_tree is None and not self._sp_failed:
            self._buildSPTree()
        return self._sp_tree

    # ------------------------------------------------------------------
    # Series-parallel decomposition
    # ------------------------------------------------------------------

    def _composeSeriesParallel(self) -> Optional[Tuple[np.ndarray, np.ndarray]]:
        """
        Exact reduction of a series-parallel graph.

        Decomposes the precedence graph into a series-parallel tree and
        composes it bottom-up, reusing every cached node whose subtree is
        unchanged. Returns None when the graph is not series-parallel, which is
        the signal to fall back to the block composition.
        """
        if self._sp_tree is None:
            if self._sp_failed:
                return None
            self._buildSPTree()
            if self._sp_tree is None:
                return None
        return self._composeNode(self._sp_tree['root'])

    def _buildSPTree(self) -> bool:
        """
        Decompose the precedence graph into a series-parallel tree.

        On success self._sp_tree holds a flat tree whose nodes are 'leaf',
        'serial', 'par', 'or' or 'loop'. On failure self._sp_tree stays None
        and self._sp_failed is set, so the decomposition is attempted once per
        topology.
        """
        self._sp_tree = None
        self._sp_failed = True

        n = len(self._activities)
        if n == 0:
            return False

        S = {
            'out_p': [-1] * n,
            'in_p': [-1] * n,
            'consumed': [False] * n,
            'type': [], 'act': [], 'kids': [], 'parent': [],
            'probs': [], 'count': [], 'alpha': [], 'T': [], 'valid': [],
        }

        # An activity may head at most one precedence and be reached by at most
        # one precedence; otherwise the graph is not series-parallel
        for p, prec in enumerate(self._precedences):
            for name in prec.pre_acts:
                i = self.getActivityIndex(name)
                if i < 0 or S['out_p'][i] >= 0:
                    return False
                S['out_p'][i] = p
            for name in prec.post_acts:
                j = self.getActivityIndex(name)
                if j < 0 or S['in_p'][j] >= 0:
                    return False
                S['in_p'][j] = p

        starts = [i for i in range(n) if S['in_p'][i] < 0]
        if len(starts) != 1:
            return False

        kids, status, _, ok = self._spParseSeq(S, starts[0], [])
        if not ok or status != 'end' or not all(S['consumed']):
            return False

        root = self._spSerialNode(S, kids)
        if root < 0:
            return False

        tree = {
            'type': S['type'], 'act': S['act'], 'kids': S['kids'],
            'parent': S['parent'], 'probs': S['probs'], 'count': S['count'],
            'alpha': S['alpha'], 'T': S['T'], 'valid': S['valid'],
            'root': root,
        }
        leaf_of = [-1] * n
        for k, t in enumerate(tree['type']):
            if t == 'leaf':
                leaf_of[tree['act'][k]] = k
        tree['leaf_of'] = leaf_of
        tree['execs'] = Workflow._spExecutionCounts(tree)

        self._sp_tree = tree
        self._sp_failed = False
        return True

    def _spParseSeq(self, S: Dict[str, Any], cur: int,
                    stop_set: List[int]) -> Tuple[List[int], str, int, bool]:
        """
        Parse a maximal sequence of blocks starting at CUR.

        Returns the nodes of the sequence and how it terminated:
            'end'  - no successor
            'stop' - reached an activity owned by the caller (stop_at)
            'join' - reached a join precedence (stop_at is its index)
        """
        kids: List[int] = []

        while True:
            if cur < 0:
                return kids, 'end', -1, True
            if stop_set and cur in stop_set:
                return kids, 'stop', cur, True
            if S['consumed'][cur]:
                return kids, 'end', -1, False
            S['consumed'][cur] = True
            kids.append(self._spAddNode(S, 'leaf', cur, [], None, 0.0))

            p = S['out_p'][cur]
            if p < 0:
                return kids, 'end', -1, True
            prec = self._precedences[p]
            if len(prec.pre_acts) > 1:
                # CUR is the tail of a branch: the caller composes the join
                return kids, 'join', p, True

            post_inds = self._spIndicesOf(prec.post_acts)
            if any(i < 0 for i in post_inds):
                return kids, 'end', -1, False

            if prec.post_type == ActivityPrecedenceType.POST_AND:
                knode, cur, ok = self._spParseFork(S, post_inds, None, stop_set, True)
                if not ok:
                    return kids, 'end', -1, False
                kids.append(knode)
            elif prec.post_type == ActivityPrecedenceType.POST_OR:
                probs = None if prec.post_params is None else np.ravel(
                    np.asarray(prec.post_params, dtype=float))
                if probs is None or probs.size != len(post_inds):
                    return kids, 'end', -1, False
                knode, cur, ok = self._spParseFork(S, post_inds, probs, stop_set, False)
                if not ok:
                    return kids, 'end', -1, False
                kids.append(knode)
            elif prec.post_type == ActivityPrecedenceType.POST_LOOP:
                knode, cur, ok = self._spParseLoop(S, post_inds, prec.post_params, stop_set)
                if not ok:
                    return kids, 'end', -1, False
                kids.append(knode)
            elif prec.post_type == ActivityPrecedenceType.POST_SEQ:
                if len(post_inds) != 1:
                    return kids, 'end', -1, False
                cur = post_inds[0]
            else:
                # POST_CACHE and any other pattern is not a workflow
                # composition rule
                return kids, 'end', -1, False

    def _spParseFork(self, S: Dict[str, Any], branch_heads: List[int],
                     probs: Optional[np.ndarray], stop_set: List[int],
                     is_and: bool) -> Tuple[int, int, bool]:
        """Parse the branches of a fork and their join."""
        nb = len(branch_heads)
        branch_nodes = [0] * nb
        bstatus = [''] * nb
        bstop = [-1] * nb

        for b in range(nb):
            bkids, st, sa, okb = self._spParseSeq(S, branch_heads[b], stop_set)
            if not okb:
                return -1, -1, False
            bn = self._spSerialNode(S, bkids)
            if bn < 0:
                return -1, -1, False
            branch_nodes[b] = bn
            bstatus[b] = st
            bstop[b] = sa

        if all(st == 'join' for st in bstatus):
            if any(x != bstop[0] for x in bstop):
                return -1, -1, False
            join_prec = self._precedences[bstop[0]]
            if len(join_prec.pre_acts) != nb:
                return -1, -1, False
            if is_and:
                if join_prec.pre_type != ActivityPrecedenceType.PRE_AND:
                    return -1, -1, False
            else:
                if join_prec.pre_type != ActivityPrecedenceType.PRE_OR:
                    return -1, -1, False
            post_inds = self._spIndicesOf(join_prec.post_acts)
            if len(post_inds) != 1 or post_inds[0] < 0:
                return -1, -1, False
            next_act = post_inds[0]
        elif all(st == 'end' for st in bstatus):
            # Branches terminate the workflow. An AND-fork with no join still
            # synchronises at the end of the workflow
            next_act = -1
        elif (not is_and) and all(st == 'stop' for st in bstatus) \
                and all(x == bstop[0] for x in bstop):
            next_act = bstop[0]
        else:
            return -1, -1, False

        if is_and:
            knode = self._spAddNode(S, 'par', -1, branch_nodes, None, 0.0)
        else:
            knode = self._spAddNode(S, 'or', -1, branch_nodes, probs, 0.0)
        return knode, next_act, True

    def _spParseLoop(self, S: Dict[str, Any], post_inds: List[int],
                     counts: Any, stop_set: List[int]) -> Tuple[int, int, bool]:
        """
        Parse a loop block.

        The last post activity is the continuation after the loop; the others
        form the loop body, which repeats geometrically.
        """
        if counts is None or np.size(counts) != 1:
            return -1, -1, False
        count = float(np.ravel(np.asarray(counts, dtype=float))[0])

        if len(post_inds) >= 2:
            body_acts = list(post_inds[:-1])
            end_act = post_inds[-1]
        else:
            body_acts = [post_inds[0]]
            end_act = -1

        loop_stop = list(stop_set) + list(body_acts)
        if end_act >= 0:
            loop_stop.append(end_act)

        body_kids: List[int] = []
        j = 0
        while j < len(body_acts):
            a = body_acts[j]
            if S['consumed'][a]:
                j += 1
                continue
            this_stop = [x for x in loop_stop if x != a]
            kk, st, sa, okb = self._spParseSeq(S, a, this_stop)
            if not okb:
                return -1, -1, False
            body_kids.extend(kk)
            if st == 'end':
                j += 1
            elif st == 'stop':
                if sa in body_acts:
                    j = body_acts.index(sa)
                elif end_act >= 0 and sa == end_act:
                    j = len(body_acts)
                else:
                    return -1, -1, False
            else:
                # A join reached from inside the body crosses the loop
                # boundary, so the graph is not series-parallel
                return -1, -1, False

        body_node = self._spSerialNode(S, body_kids)
        if body_node < 0:
            return -1, -1, False

        knode = self._spAddNode(S, 'loop', -1, [body_node], None, count)
        return knode, end_act, True

    def _spSerialNode(self, S: Dict[str, Any], kids: List[int]) -> int:
        """
        Wrap a list of nodes in a serial node.

        A single node is returned as is, so the tree carries no trivial
        one-child serial nodes.
        """
        if not kids:
            return -1
        if len(kids) == 1:
            return kids[0]
        return self._spAddNode(S, 'serial', -1, kids, None, 0.0)

    @staticmethod
    def _spAddNode(S: Dict[str, Any], node_type: str, act: int,
                   kids: List[int], probs: Optional[np.ndarray],
                   count: float) -> int:
        """Append a node to the series-parallel tree."""
        k = len(S['type'])
        S['type'].append(node_type)
        S['act'].append(act)
        S['kids'].append(list(kids))
        S['probs'].append(probs)
        S['count'].append(count)
        S['alpha'].append(None)
        S['T'].append(None)
        S['valid'].append(False)
        S['parent'].append(-1)
        for c in kids:
            S['parent'][c] = k
        return k

    def _spIndicesOf(self, names: List[str]) -> List[int]:
        """Map a list of activity names to indices."""
        return [self.getActivityIndex(nm) for nm in names]

    def _composeNode(self, k: int) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compose the law of one series-parallel node.

        Cached nodes are returned untouched, so a demand change only recomposes
        the path from the dirty leaf to the root.
        """
        tree = self._sp_tree
        if tree['valid'][k]:
            return tree['alpha'][k], tree['T'][k]

        kids = tree['kids'][k]
        ntype = tree['type'][k]
        if ntype == 'leaf':
            alpha, T = self._activities[tree['act'][k]].getPHRepresentation()
        elif ntype == 'serial':
            alpha, T = self._composeNode(kids[0])
            for i in range(1, len(kids)):
                a2, T2 = self._composeNode(kids[i])
                alpha, T = Workflow._composeSerial(alpha, T, a2, T2)
        elif ntype == 'par':
            alpha, T = self._composeNode(kids[0])
            for i in range(1, len(kids)):
                a2, T2 = self._composeNode(kids[i])
                alpha, T = Workflow._composeParallel(alpha, T, a2, T2)
        elif ntype == 'or':
            alphas = []
            Ts = []
            for c in kids:
                ai, Ti = self._composeNode(c)
                alphas.append(ai)
                Ts.append(Ti)
            alpha, T = Workflow._composeMixture(alphas, Ts, tree['probs'][k])
        elif ntype == 'loop':
            a1, T1 = self._composeNode(kids[0])
            alpha, T = Workflow._composeLoopGeometric(a1, T1, tree['count'][k])
        else:
            raise ValueError("Unknown series-parallel node type '%s'." % ntype)

        tree['alpha'][k] = alpha
        tree['T'][k] = T
        tree['valid'][k] = True
        return alpha, T

    @staticmethod
    def _invalidateBranch(tree: Dict[str, Any], k: int) -> None:
        """Invalidate a node and all its ancestors."""
        while k >= 0:
            tree['valid'][k] = False
            k = tree['parent'][k]

    @staticmethod
    def _spExecutionCounts(tree: Dict[str, Any]) -> List[float]:
        """
        Expected executions of each node per workflow run.

        Serial and parallel children inherit the count of their parent, an OR
        branch is weighted by its probability, and a loop body is weighted by
        the loop count. This is the weight by which a layer result is split
        back over entries, activities and calls.
        """
        n_nodes = len(tree['type'])
        execs = [0.0] * n_nodes
        execs[tree['root']] = 1.0
        stack = [tree['root']]
        while stack:
            k = stack.pop()
            for i, c in enumerate(tree['kids'][k]):
                if tree['type'][k] == 'or':
                    execs[c] = execs[k] * float(tree['probs'][k][i])
                elif tree['type'][k] == 'loop':
                    execs[c] = execs[k] * float(tree['count'][k])
                else:
                    execs[c] = execs[k]
                stack.append(c)
        return execs

    def _buildCTMC(self) -> Tuple[np.ndarray, np.ndarray]:
        """Build the CTMC representation of the workflow."""
        n = len(self._activities)

        if n == 1:
            return self._activities[0].getPHRepresentation()

        # Analyze workflow structure
        structure = self._analyzeStructure()

        # Check if simple serial workflow
        if not structure['forks'] and not structure['joins'] and not structure['loops']:
            return self._composeSerialWorkflow(structure['adj_list'])

        return self._composeComplexWorkflow(structure)

    def _analyzeStructure(self) -> Dict[str, Any]:
        """Analyze the workflow structure to identify patterns."""
        n = len(self._activities)

        adj_list = [[] for _ in range(n)]
        in_degree = [0] * n
        out_degree = [0] * n
        forks = []
        joins = []
        loops = []

        for prec in self._precedences:
            pre_indices = [self._activity_map[name] for name in prec.pre_acts]
            post_indices = [self._activity_map[name] for name in prec.post_acts]

            # Build adjacency list
            for pre_idx in pre_indices:
                for post_idx in post_indices:
                    adj_list[pre_idx].append(post_idx)
                    out_degree[pre_idx] += 1
                    in_degree[post_idx] += 1

            # Identify patterns
            if prec.post_type == ActivityPrecedenceType.POST_AND:
                forks.append({
                    'type': 'and',
                    'pre_act': pre_indices[0],
                    'post_acts': post_indices
                })
            elif prec.post_type == ActivityPrecedenceType.POST_OR:
                forks.append({
                    'type': 'or',
                    'pre_act': pre_indices[0],
                    'post_acts': post_indices,
                    'probs': prec.post_params
                })
            elif prec.post_type == ActivityPrecedenceType.POST_LOOP:
                loop_acts = post_indices[:-1] if len(post_indices) > 1 else post_indices
                end_act = post_indices[-1] if len(post_indices) > 1 else -1
                count = prec.post_params[0] if prec.post_params is not None else 1.0
                loops.append({
                    'pre_act': pre_indices[0],
                    'loop_acts': loop_acts,
                    'end_act': end_act,
                    'count': count
                })

            if prec.pre_type == ActivityPrecedenceType.PRE_AND:
                joins.append({
                    'type': 'and',
                    'pre_acts': pre_indices,
                    'post_act': post_indices[0]
                })
            elif prec.pre_type == ActivityPrecedenceType.PRE_OR:
                joins.append({
                    'type': 'or',
                    'pre_acts': pre_indices,
                    'post_act': post_indices[0]
                })

        return {
            'adj_list': adj_list,
            'in_degree': in_degree,
            'out_degree': out_degree,
            'forks': forks,
            'joins': joins,
            'loops': loops
        }

    def _topologicalSort(self, adj_list: List[List[int]]) -> List[int]:
        """Perform topological sort on the activity graph."""
        n = len(self._activities)
        in_deg = [0] * n

        for neighbors in adj_list:
            for j in neighbors:
                in_deg[j] += 1

        queue = [i for i in range(n) if in_deg[i] == 0]
        order = []

        while queue:
            curr = queue.pop(0)
            order.append(curr)

            for next_node in adj_list[curr]:
                in_deg[next_node] -= 1
                if in_deg[next_node] == 0:
                    queue.append(next_node)

        # Add any remaining nodes (for cycles or disconnected)
        order_set = set(order)
        for i in range(n):
            if i not in order_set:
                order.append(i)

        return order

    def _composeSerialWorkflow(self, adj_list: List[List[int]]) -> Tuple[np.ndarray, np.ndarray]:
        """Compose a simple serial workflow."""
        order = self._topologicalSort(adj_list)

        alpha, T = self._activities[order[0]].getPHRepresentation()

        for i in range(1, len(order)):
            next_alpha, next_T = self._activities[order[i]].getPHRepresentation()
            alpha, T = self._composeSerial(alpha, T, next_alpha, next_T)

        return alpha, T

    def _composeComplexWorkflow(self, structure: Dict[str, Any]) -> Tuple[np.ndarray, np.ndarray]:
        """Compose a complex workflow with forks, joins, and loops."""
        n = len(self._activities)

        # Initialize block representations
        block_alpha = [None] * n
        block_T = [None] * n
        # is_absorbed[i] = True means activity i was absorbed into another block and should be skipped
        is_absorbed = [False] * n
        # block_updated[i] = True means block_alpha[i] contains a composed block (not just the activity)
        block_updated = [False] * n

        for i in range(n):
            block_alpha[i], block_T[i] = self._activities[i].getPHRepresentation()

        # Process loops
        for loop in structure['loops']:
            pre_idx = loop['pre_act']
            loop_acts = loop['loop_acts']
            end_act = loop['end_act']
            count = float(loop['count'])

            # Compose loop activities
            if len(loop_acts) == 1:
                alpha_loop, T_loop = self._activities[loop_acts[0]].getPHRepresentation()
            else:
                alpha_loop, T_loop = self._activities[loop_acts[0]].getPHRepresentation()
                for j in range(1, len(loop_acts)):
                    next_alpha, next_T = self._activities[loop_acts[j]].getPHRepresentation()
                    alpha_loop, T_loop = self._composeSerial(alpha_loop, T_loop, next_alpha, next_T)

            # POST_LOOP repeats the body a geometric number of times
            conv_alpha, conv_T = self._composeLoopGeometric(alpha_loop, T_loop, count)

            # Compose with pre-activity
            result_alpha, result_T = self._composeSerial(
                block_alpha[pre_idx], block_T[pre_idx],
                conv_alpha, conv_T
            )

            # Compose with end activity
            if end_act >= 0:
                end_alpha, end_T = self._activities[end_act].getPHRepresentation()
                result_alpha, result_T = self._composeSerial(result_alpha, result_T, end_alpha, end_T)
                is_absorbed[end_act] = True

            block_alpha[pre_idx] = result_alpha
            block_T[pre_idx] = result_T
            block_updated[pre_idx] = True
            for idx in loop_acts:
                is_absorbed[idx] = True

        # Process AND-forks with matching joins
        for fork in structure['forks']:
            if fork['type'] == 'and':
                matching_join = self._findMatchingJoin(fork['post_acts'], structure['joins'], 'and')

                if matching_join is not None:
                    pre_idx = fork['pre_act']
                    post_idx = matching_join['post_act']

                    # Compose parallel block
                    par_alpha, par_T = self._composeAndForkBlock(
                        fork['post_acts'], block_alpha, block_T
                    )

                    if not block_updated[pre_idx]:
                        result_alpha, result_T = self._composeSerial(
                            block_alpha[pre_idx], block_T[pre_idx],
                            par_alpha, par_T
                        )
                    else:
                        result_alpha, result_T = self._composeSerial(
                            block_alpha[pre_idx], block_T[pre_idx],
                            par_alpha, par_T
                        )

                    if not is_absorbed[post_idx]:
                        result_alpha, result_T = self._composeSerial(
                            result_alpha, result_T,
                            block_alpha[post_idx], block_T[post_idx]
                        )

                    block_alpha[pre_idx] = result_alpha
                    block_T[pre_idx] = result_T
                    block_updated[pre_idx] = True
                    for idx in fork['post_acts']:
                        is_absorbed[idx] = True
                    is_absorbed[post_idx] = True

        # Process OR-forks
        for fork in structure['forks']:
            if fork['type'] == 'or':
                matching_join = self._findMatchingJoin(fork['post_acts'], structure['joins'], 'or')

                pre_idx = fork['pre_act']

                # Compose OR-fork block
                or_alpha, or_T = self._composeOrForkBlock(
                    fork['post_acts'], fork['probs'], block_alpha, block_T
                )

                if not block_updated[pre_idx]:
                    result_alpha, result_T = self._composeSerial(
                        block_alpha[pre_idx], block_T[pre_idx],
                        or_alpha, or_T
                    )
                else:
                    result_alpha, result_T = self._composeSerial(
                        block_alpha[pre_idx], block_T[pre_idx],
                        or_alpha, or_T
                    )

                if matching_join is not None:
                    post_idx = matching_join['post_act']
                    if not is_absorbed[post_idx]:
                        result_alpha, result_T = self._composeSerial(
                            result_alpha, result_T,
                            block_alpha[post_idx], block_T[post_idx]
                        )
                        is_absorbed[post_idx] = True

                block_alpha[pre_idx] = result_alpha
                block_T[pre_idx] = result_T
                block_updated[pre_idx] = True
                for idx in fork['post_acts']:
                    is_absorbed[idx] = True

        # Compose remaining activities in topological order
        # Skip absorbed activities, include all others (whether their block was updated or not)
        order = self._topologicalSort(structure['adj_list'])
        alpha = None
        T = None

        for idx in order:
            if not is_absorbed[idx]:
                if alpha is None:
                    alpha = block_alpha[idx]
                    T = block_T[idx]
                else:
                    alpha, T = self._composeSerial(alpha, T, block_alpha[idx], block_T[idx])

        if alpha is None:
            alpha, T = self._activities[0].getPHRepresentation()

        return alpha, T

    def _findMatchingJoin(self, post_acts: List[int],
                          joins: List[Dict], join_type: str) -> Optional[Dict]:
        """Find a matching join for a fork's post activities."""
        post_set = set(post_acts)

        for join in joins:
            if join['type'] == join_type:
                pre_set = set(join['pre_acts'])
                if post_set == pre_set:
                    return join
        return None

    def _composeAndForkBlock(self, parallel_inds: List[int],
                             block_alpha: List[np.ndarray],
                             block_T: List[np.ndarray]) -> Tuple[np.ndarray, np.ndarray]:
        """Compose parallel activities (AND-fork)."""
        alpha = block_alpha[parallel_inds[0]]
        T = block_T[parallel_inds[0]]

        for i in range(1, len(parallel_inds)):
            alpha, T = self._composeParallel(
                alpha, T,
                block_alpha[parallel_inds[i]], block_T[parallel_inds[i]]
            )

        return alpha, T

    def _composeOrForkBlock(self, branch_inds: List[int],
                            probs: np.ndarray,
                            block_alpha: List[np.ndarray],
                            block_T: List[np.ndarray]) -> Tuple[np.ndarray, np.ndarray]:
        """Compose branching activities (OR-fork)."""
        alphas = [block_alpha[idx] for idx in branch_inds]
        Ts = [block_T[idx] for idx in branch_inds]
        return Workflow._composeMixture(alphas, Ts, np.ravel(np.asarray(probs, dtype=float)))

    @staticmethod
    def _composeSerial(alpha1: np.ndarray, T1: np.ndarray,
                       alpha2: np.ndarray, T2: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Compose two PH distributions in series."""
        ZERO = 1e-14
        n1 = T1.shape[0]
        n2 = T2.shape[0]

        # Absorption rates from first distribution
        e1 = np.ones((n1, 1))
        abs_rate1 = -T1 @ e1

        # Build combined T matrix
        T_out = np.zeros((n1 + n2, n1 + n2))
        T_out[:n1, :n1] = T1
        T_out[n1:, n1:] = T2

        # Transition from T1 to T2 via absorption
        alpha2_flat = alpha2.flatten()
        for r in range(n1):
            for c in range(n2):
                val = abs_rate1[r, 0] * alpha2_flat[c]
                if abs(val) > ZERO:
                    T_out[r, n1 + c] = val

        # A defective alpha1 carries an atom at zero, which starts the second
        # law immediately; this is aph_simplify pattern 1
        alpha1_flat = alpha1.flatten()
        defect1 = 1.0 - float(np.sum(alpha1_flat))
        alpha_out = np.zeros((1, n1 + n2))
        alpha_out[0, :n1] = alpha1_flat
        alpha_out[0, n1:] = defect1 * alpha2_flat

        return alpha_out, T_out

    @staticmethod
    def _composeParallel(alpha1: np.ndarray, T1: np.ndarray,
                         alpha2: np.ndarray, T2: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Compose two PH distributions in parallel (synchronization)."""
        ZERO = 1e-14
        n1 = T1.shape[0]
        n2 = T2.shape[0]

        e1 = np.ones((n1, 1))
        e2 = np.ones((n2, 1))
        abs_rate1 = -T1 @ e1
        abs_rate2 = -T2 @ e2

        n_both = n1 * n2
        n_only1 = n1
        n_only2 = n2
        n_total = n_both + n_only1 + n_only2

        T_out = np.zeros((n_total, n_total))

        # Kronecker sum for both running
        T_both = np.kron(T1, np.eye(n2)) + np.kron(np.eye(n1), T2)
        T_out[:n_both, :n_both] = T_both

        # Transitions when one completes
        for i in range(n1):
            for j in range(n2):
                both_idx = i * n2 + j
                only1_idx = n_both + i
                only2_idx = n_both + n_only1 + j

                # T2 completes, T1 continues
                T_out[both_idx, only1_idx] += abs_rate2[j, 0]
                # T1 completes, T2 continues
                T_out[both_idx, only2_idx] += abs_rate1[i, 0]

        # Sub-matrices for when only one is running
        T_out[n_both:n_both + n_only1, n_both:n_both + n_only1] = T1
        T_out[n_both + n_only1:, n_both + n_only1:] = T2

        # Initial distribution (both start together)
        alpha_out = np.zeros((1, n_total))
        alpha1_flat = alpha1.flatten()
        alpha2_flat = alpha2.flatten()
        for i in range(n1):
            for j in range(n2):
                both_idx = i * n2 + j
                alpha_out[0, both_idx] = alpha1_flat[i] * alpha2_flat[j]

        return alpha_out, T_out

    @staticmethod
    def _composeMixture(alphas: List[np.ndarray], Ts: List[np.ndarray],
                        probs: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
        Probabilistic mixture of several PH laws.

        Block-diagonal generator whose initial vector picks branch I with
        probability PROBS[I]. This is aph_simplify pattern 3 generalised to any
        number of branches.
        """
        sizes = [Ti.shape[0] for Ti in Ts]
        total = int(sum(sizes))

        T_out = np.zeros((total, total))
        alpha_out = np.zeros((1, total))

        offset = 0
        for i, Ti in enumerate(Ts):
            ni = sizes[i]
            T_out[offset:offset + ni, offset:offset + ni] = Ti
            alpha_out[0, offset:offset + ni] = float(probs[i]) * np.ravel(alphas[i])
            offset += ni

        return alpha_out, T_out

    @staticmethod
    def _composeLoopGeometric(alpha: np.ndarray, T: np.ndarray,
                              count: float) -> Tuple[np.ndarray, np.ndarray]:
        """
        Geometric repetition of a PH law.

        Implements the LayeredNetwork POST_LOOP semantics, in which the number
        of executions of the loop body is geometric of mean COUNT. For COUNT>=1
        the body runs at least once and repeats on absorption with probability
        P = 1-1/COUNT, so

            T_OUT = T + P/D * (-T*e)*ALPHA,   ALPHA_OUT = ALPHA/D

        with D = 1 - P*(1-ALPHA*e) the correction for an atom at zero in ALPHA.
        The order of the law is that of the body, unlike the COUNT-fold
        convolution _composeRepeat, and the mean is COUNT times the mean of the
        body in both cases.

        For COUNT<1 the body is executed at most once, with probability COUNT,
        which is how a fractional loop count is read when the LQN activity
        graph is built; the skipped branch is an immediate phase, the
        representation of a zero-time activity used throughout this class.
        """
        alpha = np.asarray(alpha, dtype=float).reshape(1, -1)
        T = np.asarray(T, dtype=float)
        n = T.shape[0]
        e = np.ones((n, 1))

        if count <= 0:
            # zero-time branch
            return np.array([[1.0]]), np.array([[-GlobalConstants.Immediate]])

        if abs(count - 1.0) <= GlobalConstants.FineTol:
            return alpha.copy(), T.copy()

        if count < 1:
            # Executed with probability COUNT, skipped otherwise
            alpha_out = np.zeros((1, n + 1))
            alpha_out[0, :n] = count * alpha[0, :]
            alpha_out[0, n] = 1.0 - count
            T_out = np.zeros((n + 1, n + 1))
            T_out[:n, :n] = T
            T_out[n, n] = -GlobalConstants.Immediate  # zero-time skip branch
            return alpha_out, T_out

        p = 1.0 - 1.0 / count
        defect = 1.0 - float(np.sum(alpha))
        denom = 1.0 - p * defect
        alpha_out = alpha / denom
        T_out = T + (p / denom) * ((-T @ e) @ alpha)
        return alpha_out, T_out

    @staticmethod
    def isAcyclicGenerator(T: np.ndarray) -> bool:
        """
        True if the phase graph of T has no cycle.

        A geometric loop over a body of two or more phases closes a cycle, so
        the composed law is a PH and not an APH.
        """
        T = np.asarray(T, dtype=float)
        n = T.shape[0]
        A = np.abs(T) > GlobalConstants.ArcTol
        np.fill_diagonal(A, False)
        in_deg = A.sum(axis=0).astype(int)
        queue = [i for i in range(n) if in_deg[i] == 0]
        visited = 0
        while queue:
            curr = queue.pop(0)
            visited += 1
            for sidx in np.flatnonzero(A[curr, :]):
                in_deg[sidx] -= 1
                if in_deg[sidx] == 0:
                    queue.append(int(sidx))
        return visited == n

    @staticmethod
    def _composeRepeat(alpha: np.ndarray, T: np.ndarray,
                       count: int) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compose a PH distribution convolved COUNT times.

        This is the deterministic fold, kept for a caller that genuinely needs
        an exact repetition count. POST_LOOP is geometric and uses
        _composeLoopGeometric instead.
        """
        count = int(count)
        if count <= 0:
            return np.array([[1.0]]), np.array([[-GlobalConstants.Immediate]])

        if count == 1:
            return alpha.copy(), T.copy()

        alpha_out = alpha.copy()
        T_out = T.copy()

        for _ in range(1, count):
            alpha_out, T_out = Workflow._composeSerial(alpha_out, T_out, alpha, T)

        return alpha_out, T_out

    def getMean(self) -> float:
        """Get the mean execution time of the workflow."""
        alpha, T = self.toPH()
        # Mean = -alpha * T^(-1) * e
        e = np.ones((T.shape[0], 1))
        try:
            T_inv = np.linalg.inv(T)
            mean = -alpha @ T_inv @ e
            return float(mean[0, 0])
        except np.linalg.LinAlgError:
            # Fallback: sum of activity means
            return sum(a.host_demand_mean for a in self._activities)

    def getSCV(self) -> float:
        """Get the squared coefficient of variation."""
        alpha, T = self.toPH()
        e = np.ones((T.shape[0], 1))
        try:
            T_inv = np.linalg.inv(T)
            mean = -alpha @ T_inv @ e
            mean2 = 2 * alpha @ T_inv @ T_inv @ e
            variance = float(mean2[0, 0]) - float(mean[0, 0]) ** 2
            return variance / (float(mean[0, 0]) ** 2) if float(mean[0, 0]) > 0 else 1.0
        except np.linalg.LinAlgError:
            return 1.0

    def sample(self, n: int = 1) -> np.ndarray:
        """
        Generate random samples from the workflow's phase-type distribution.

        Uses the PH representation (alpha, T) to simulate absorption times.
        For a PH distribution with sub-generator T:
        - T[i,i] < 0: exit rate from state i is -T[i,i]
        - T[i,j] >= 0 for i != j: transition rate from i to j
        - Absorption rate from i: -sum(T[i,:])

        Args:
            n: Number of samples to generate

        Returns:
            numpy array of n samples (absorption times)
        """
        alpha, T = self.toPH()
        n_phases = T.shape[0]
        samples = np.zeros(n)

        alpha_flat = alpha.flatten()

        # Precompute rates
        # Exit rate from each state = -T[i,i]
        exit_rates = -np.diag(T)
        # Absorption rate from each state = -sum(T[i,:])
        abs_rates = -np.sum(T, axis=1)

        for i in range(n):
            time = 0.0
            # Choose initial state according to alpha
            state = np.random.choice(n_phases, p=alpha_flat)

            while True:
                # Time in current state (exponential with rate = exit_rate)
                rate = exit_rates[state]
                time += np.random.exponential(1.0 / rate)

                # Determine next event: transition or absorption
                # Probability of absorption = abs_rates[state] / exit_rates[state]
                abs_prob = abs_rates[state] / rate
                if np.random.random() < abs_prob:
                    # Absorbed
                    break
                else:
                    # Transition to another state
                    # Get transition rates (off-diagonal elements of row state)
                    trans_rates = T[state, :].copy()
                    trans_rates[state] = 0  # No self-transitions
                    trans_sum = np.sum(trans_rates)
                    if trans_sum > 0:
                        trans_probs = trans_rates / trans_sum
                        state = np.random.choice(n_phases, p=trans_probs)
                    else:
                        # No transitions possible, absorb
                        break

            samples[i] = time

        return samples

    @staticmethod
    def fromWfCommons(json_file: str) -> 'Workflow':
        """
        Load a workflow from a WfCommons JSON file.

        WfCommons (https://github.com/wfcommons/workflow-schema)
        is a standard format for representing scientific workflow traces.

        Args:
            json_file: Path to the WfCommons JSON file

        Returns:
            Workflow object
        """
        with open(json_file, 'r') as f:
            data = json.load(f)

        name = data.get('name', 'WfCommons_Workflow')
        wf = Workflow(name)

        # Get jobs from workflow
        jobs = data.get('workflow', {}).get('jobs', [])
        if not jobs:
            jobs = data.get('jobs', [])

        # Create activities
        job_map = {}
        for job in jobs:
            job_name = job.get('name', job.get('id'))
            runtime = job.get('runtime', 1.0)
            act = wf.addActivity(job_name, runtime)
            job_map[job_name] = act

            # Store metadata
            act.metadata = {
                'files': job.get('files', []),
                'machine': job.get('machine'),
                'args': job.get('args', [])
            }

        # Create precedences from dependencies
        for job in jobs:
            job_name = job.get('name', job.get('id'))
            parents = job.get('parents', [])
            for parent_name in parents:
                if parent_name in job_map and job_name in job_map:
                    wf.addPrecedence(ActivityPrecedence.Serial(parent_name, job_name))

        return wf

    def __repr__(self) -> str:
        return f"Workflow('{self._name}', activities={len(self._activities)}, precedences={len(self._precedences)})"


# Convenience aliases
Serial = ActivityPrecedence.Serial
SerialSequence = ActivityPrecedence.SerialSequence
AndFork = ActivityPrecedence.AndFork
AndJoin = ActivityPrecedence.AndJoin
OrFork = ActivityPrecedence.OrFork
OrJoin = ActivityPrecedence.OrJoin
Loop = ActivityPrecedence.Loop
