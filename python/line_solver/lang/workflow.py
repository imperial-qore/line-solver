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
- Conversion to phase-type (APH) distributions via toPH()
- Loading from WfCommons JSON format

References:
    Original Java: jar/src/main/kotlin/jline/lang/workflow/Workflow.java
"""

import numpy as np
from typing import List, Dict, Any, Optional, Tuple, Union
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
        if quorum is None:
            quorum = np.ones(len(pre_acts))
        return ActivityPrecedence(
            pre_acts=pre_names,
            post_acts=[post_name],
            pre_type=ActivityPrecedenceType.PRE_AND,
            post_type=ActivityPrecedenceType.POST_SEQ,
            pre_params=quorum
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
             end_act: Optional[Union[str, 'WorkflowActivity']] = None,
             count: float = 1.0) -> 'ActivityPrecedence':
        """Create a loop: pre_act -> loop_acts (repeated count times) -> end_act."""
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


class WorkflowActivity:
    """
    A computational activity in a Workflow.

    Represents a single activity with a host demand (service time distribution).
    Activities can be composed into workflows using precedence relationships.

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

        # Set host demand
        if isinstance(host_demand, (int, float)):
            self._host_demand_mean = float(host_demand)
            self._host_demand_scv = 1.0  # Exponential
            self._distribution = None
        else:
            # Assume it's a Distribution object
            self._distribution = host_demand
            self._host_demand_mean = host_demand.getMean() if hasattr(host_demand, 'getMean') else 1.0
            self._host_demand_scv = host_demand.getSCV() if hasattr(host_demand, 'getSCV') else 1.0

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
            self._host_demand_mean = float(value)
            self._host_demand_scv = 1.0
            self._distribution = None
        else:
            self._distribution = value
            self._host_demand_mean = value.getMean() if hasattr(value, 'getMean') else 1.0
            self._host_demand_scv = value.getSCV() if hasattr(value, 'getSCV') else 1.0

    def getPHRepresentation(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Get phase-type representation of this activity.

        Returns:
            Tuple of (alpha, T) where:
                alpha: Initial probability vector (1 x n)
                T: Sub-generator matrix (n x n)
        """
        FINE_TOL = 1e-14

        # Handle immediate (zero service time)
        if self._host_demand_mean <= FINE_TOL:
            return np.array([[1.0]]), np.array([[-1e10]])

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
    - Loops with fixed iteration counts

    Example:
        >>> wf = Workflow("ServiceWorkflow")
        >>> a = wf.addActivity("A", 1.0)
        >>> b = wf.addActivity("B", 2.0)
        >>> c = wf.addActivity("C", 0.5)
        >>> wf.addPrecedence(ActivityPrecedence.Serial(a, b))
        >>> wf.addPrecedence(ActivityPrecedence.Serial(b, c))
        >>> alpha, T = wf.toPH()

    References:
        Original Java: jar/src/main/kotlin/jline/lang/workflow/Workflow.java
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
        self._cached_ph = None
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
        self._cached_ph = None

    # Snake_case alias
    add_precedence = addPrecedence

    def getActivity(self, name: str) -> Optional[WorkflowActivity]:
        """Get an activity by name."""
        idx = self._activity_map.get(name)
        if idx is None:
            return None
        return self._activities[idx]

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
        FINE_TOL = 1e-14
        for prec in self._precedences:
            if prec.post_type == ActivityPrecedenceType.POST_OR:
                if prec.post_params is None:
                    return False, "OR-fork must have probabilities."
                prob_sum = np.sum(prec.post_params)
                if abs(prob_sum - 1.0) > FINE_TOL:
                    return False, "OR-fork probabilities must sum to 1."

        return True, ""

    def toPH(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Convert workflow to phase-type (APH) representation.

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

        alpha, T = self._buildCTMC()
        self._cached_ph = (alpha, T)
        return alpha, T

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
            count = int(loop['count'])

            # Compose loop activities
            if len(loop_acts) == 1:
                alpha_loop, T_loop = self._activities[loop_acts[0]].getPHRepresentation()
            else:
                alpha_loop, T_loop = self._activities[loop_acts[0]].getPHRepresentation()
                for j in range(1, len(loop_acts)):
                    next_alpha, next_T = self._activities[loop_acts[j]].getPHRepresentation()
                    alpha_loop, T_loop = self._composeSerial(alpha_loop, T_loop, next_alpha, next_T)

            # Repeat for count iterations
            conv_alpha, conv_T = self._composeRepeat(alpha_loop, T_loop, count)

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
        ZERO = 1e-14

        total_phases = sum(block_T[idx].shape[0] for idx in branch_inds)

        T = np.zeros((total_phases, total_phases))
        alpha = np.zeros((1, total_phases))

        offset = 0
        for i, idx in enumerate(branch_inds):
            alpha_i = block_alpha[idx]
            T_i = block_T[idx]
            n_i = T_i.shape[0]

            # Copy T block
            T[offset:offset + n_i, offset:offset + n_i] = T_i

            # Set initial probabilities
            for j in range(alpha_i.size):
                alpha[0, offset + j] = probs[i] * alpha_i.flat[j]

            offset += n_i

        return alpha, T

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

        # Combined initial distribution
        alpha_out = np.zeros((1, n1 + n2))
        alpha1_flat = alpha1.flatten()
        for i in range(n1):
            alpha_out[0, i] = alpha1_flat[i]

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
    def _composeRepeat(alpha: np.ndarray, T: np.ndarray,
                       count: int) -> Tuple[np.ndarray, np.ndarray]:
        """Compose a PH distribution repeated count times."""
        if count <= 0:
            return np.array([[1.0]]), np.array([[-1e10]])

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
