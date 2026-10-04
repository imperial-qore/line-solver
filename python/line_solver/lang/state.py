"""
State classes for LINE queueing network models (pure Python).

This module provides state representation for network analysis.
"""

from typing import Any, Dict, List, Optional, Union, TYPE_CHECKING
import numpy as np

if TYPE_CHECKING:
    from .network import Network


class State:
    """
    State representation for stochastic network models.

    Represents the system state including job populations at each node,
    phase information for multi-phase processes, and other state variables.
    """

    def __init__(self, network: Optional['Network'] = None):
        """
        Initialize a state for the network.

        Args:
            network: The network this state belongs to.
        """
        self._network = network
        self._state: Dict[int, np.ndarray] = {}  # stateful_idx -> state vector

    def get(self, stateful_idx: int) -> Optional[np.ndarray]:
        """Get state for a stateful node."""
        return self._state.get(stateful_idx)

    def set(self, stateful_idx: int, state: np.ndarray):
        """Set state for a stateful node."""
        self._state[stateful_idx] = np.array(state)

    def toArray(self) -> np.ndarray:
        """Convert to a flat array representation."""
        if not self._state:
            return np.array([])
        # Concatenate all state vectors
        return np.concatenate([self._state[i] for i in sorted(self._state.keys())])

    @staticmethod
    def fromMarginal(model: 'Network', node_idx: int,
                     n: Union[List[int], np.ndarray]) -> np.ndarray:
        """
        Generate state space with specific marginal job counts at a node.

        Creates all possible network states where the specified node has
        exactly n[r] jobs of class r.

        Args:
            model: Network model
            node_idx: Node index (0-based)
            n: Vector of job counts per class

        Returns:
            State space matrix where each row is a valid state.
        """
        from ..api.state.marginal import fromMarginal as _fromMarginal

        # Get network struct from model
        if hasattr(model, 'get_struct'):
            sn = model.get_struct()
        elif hasattr(model, 'getStruct'):
            sn = model.getStruct()
        elif hasattr(model, '_sn'):
            sn = model._sn
        else:
            raise ValueError("Cannot get network struct from model")

        return _fromMarginal(sn, State._node_index(model, node_idx), n)

    @staticmethod
    def fromMarginalAndRunning(model: 'Network', node_idx: int,
                                n: Union[List[int], np.ndarray],
                                s: Union[List[int], np.ndarray]) -> np.ndarray:
        """
        Generate state space with specific marginal and running job counts.

        Creates states where node has n[r] jobs of class r total,
        with s[r] jobs of class r currently in service (running).

        Args:
            model: Network model
            node_idx: Node index (0-based)
            n: Vector of total job counts per class
            s: Vector of running job counts per class

        Returns:
            State space matrix where each row is a valid state.
        """
        from ..api.state.marginal import fromMarginalAndRunning as _fromMarginalAndRunning

        # Get network struct from model
        if hasattr(model, 'get_struct'):
            sn = model.get_struct()
        elif hasattr(model, 'getStruct'):
            sn = model.getStruct()
        elif hasattr(model, '_sn'):
            sn = model._sn
        else:
            raise ValueError("Cannot get network struct from model")

        return _fromMarginalAndRunning(sn, State._node_index(model, node_idx), n, s)

    @staticmethod
    def fromMarginalAndStarted(model: 'Network', node_idx: int,
                                n: Union[List[int], np.ndarray],
                                s: Union[List[int], np.ndarray]) -> np.ndarray:
        """
        Generate state space with specific marginal and started job counts.

        Creates states where node has n[r] jobs of class r total,
        with s[r] jobs of class r that have started service.

        Args:
            model: Network model
            node_idx: Node index (0-based)
            n: Vector of total job counts per class
            s: Vector of started job counts per class

        Returns:
            State space matrix where each row is a valid state.
        """
        from ..api.state.marginal import fromMarginalAndStarted as _fromMarginalAndStarted

        # Get network struct from model
        if hasattr(model, 'get_struct'):
            sn = model.get_struct()
        elif hasattr(model, 'getStruct'):
            sn = model.getStruct()
        elif hasattr(model, '_sn'):
            sn = model._sn
        else:
            raise ValueError("Cannot get network struct from model")

        return _fromMarginalAndStarted(sn, State._node_index(model, node_idx), n, s)

    @staticmethod
    def fromMarg(model: 'Network', node_idx: int, ntot: int) -> np.ndarray:
        """
        Generate the state space with a given TOTAL queue length at a node.

        Class-summed counterpart of fromMarginal: it fixes only how many jobs
        the node holds ALTOGETHER, and returns the union of fromMarginal over
        every class split of ntot the node can hold. Classes disabled at the
        station are excluded from the split enumeration through classcap.

        Args:
            model: Network model
            node_idx: Node index (0-based)
            ntot: Total number of jobs at the node, all classes summed

        Returns:
            State space matrix with the requested total.

        References:
            MATLAB: matlab/src/lang/+State/fromMarg.m
        """
        from ..api.state.marginal import fromMarg as _fromMarg

        return _fromMarg(State._struct(model), State._node_index(model, node_idx), ntot)

    @staticmethod
    def fromMargAndStarted(model: 'Network', node_idx: int, ntot: int, stot: int) -> np.ndarray:
        """
        Generate the states with a given TOTAL queue length and a given TOTAL
        number of started jobs.

        Returns the union of fromMarginalAndStarted over every (n,s) pair with
        sum(n)=ntot, sum(s)=stot and s <= n elementwise.

        Args:
            model: Network model
            node_idx: Node index (0-based)
            ntot: Total number of jobs at the node, all classes summed
            stot: Total number of jobs that have started service

        Returns:
            State space matrix with the requested totals.

        References:
            MATLAB: matlab/src/lang/+State/fromMargAndStarted.m
        """
        from ..api.state.marginal import fromMargAndStarted as _fromMargAndStarted

        return _fromMargAndStarted(State._struct(model),
                                   State._node_index(model, node_idx), ntot, stot)

    @staticmethod
    def _struct(model):
        """The NetworkStruct of a model, however the model exposes it."""
        if hasattr(model, 'get_struct'):
            return model.get_struct()
        if hasattr(model, 'getStruct'):
            return model.getStruct()
        if hasattr(model, '_sn'):
            return model._sn
        raise ValueError("Cannot get network struct from model")

    @staticmethod
    def _node_index(model, node):
        """A node argument as the 0-based node index the api layer expects.

        MATLAB gets this for free: Node defines `subsindex` (lang/nodes/Node.m),
        so `sn.nodetype(node{2})` IS `sn.nodetype(2)` and every +State function
        takes a node object without knowing it. Python has no such protocol, so
        the identical example call (`State.fromMarginalAndStarted(model,
        node{2}, ...)` -> `State.from_marginal_and_started(model, node[1], ...)`)
        arrives holding a Queue. Normalizing here rather than in api/ keeps that
        layer's `ind: int` contract, which is also what the MATLAB +State
        functions document.

        `get_node_index` is 1-based, as its docstring says; these entry points
        are 0-based.
        """
        if isinstance(node, (int, np.integer)):
            return int(node)
        getter = getattr(model, 'get_node_index', None) or \
            getattr(model, 'getNodeIndex', None)
        if getter is None:
            return node  # a bare struct cannot resolve an object; leave it be
        return int(getter(node)) - 1

    @staticmethod
    def _struct_of(model):
        """Resolve a NetworkStruct from a model/struct argument."""
        if hasattr(model, 'get_struct'):
            return model.get_struct()
        if hasattr(model, 'getStruct'):
            return model.getStruct()
        if hasattr(model, '_sn'):
            return model._sn
        # Already a struct
        return model

    @staticmethod
    def toMarginal(model: 'Network', node_idx: int,
                   state_i: np.ndarray):
        """
        Extract marginal job statistics from a state for a specific node.

        Args:
            model: Network model (or NetworkStruct).
            node_idx: Node index (0-based).
            state_i: State vector/matrix for the node.

        Returns:
            Tuple ``(ni, nir, sir, kir)``:
            - ni: total jobs in the node,
            - nir: total jobs per class,
            - sir: jobs in service per class,
            - kir: jobs in service per class and phase.

        References:
            MATLAB: matlab/src/lang/+State/toMarginal.m
        """
        from ..api.state.marginal import toMarginal as _toMarginal
        sn = State._struct_of(model)
        return _toMarginal(sn, State._node_index(model, node_idx), state_i)

    @staticmethod
    def fromMarginalBounds(model: 'Network', node_idx: int,
                           lb, ub, cap=None) -> np.ndarray:
        """
        Generate all valid states whose per-class marginal lies between ``lb``
        and ``ub`` (inclusive), subject to the node capacity ``cap``.

        Args:
            model: Network model (or NetworkStruct).
            node_idx: Node index (0-based).
            lb: Lower bound on the number of resident jobs (scalar total or
                per-class vector). ``None`` is treated as all-zero.
            ub: Upper bound on the number of resident jobs (scalar total or
                per-class vector).
            cap: Total capacity at the node (defaults to the station's total
                capacity, or infinity if unbounded).

        Returns:
            State space matrix where each row is a valid state.

        References:
            MATLAB: matlab/src/lang/+State/fromMarginalBounds.m
        """
        from ..api.state.ctmc_ssg import _from_marginal_bounds
        from ..api.state.marginal import toMarginal as _toMarginal
        sn = State._struct_of(model)
        node_idx = State._node_index(model, node_idx)
        R = sn.nclasses

        ub_vec = np.atleast_1d(ub).astype(int)
        if ub_vec.size == 1:
            ub_vec = np.full(R, int(ub_vec[0]))

        # Determine capacity if not supplied.
        if cap is None:
            ist = int(sn.nodeToStation[node_idx]) if sn.isstation[node_idx] else -1
            if ist >= 0 and sn.classcap is not None and np.all(np.isfinite(sn.classcap[ist])):
                cap = int(np.sum(sn.classcap[ist]))
            else:
                cap = float('inf')

        # The building block enumerates 0 <= n <= ub and filters by capacity.
        space = _from_marginal_bounds(sn, node_idx, ub_vec, cap)
        if not isinstance(space, np.ndarray) or space.size == 0:
            return space

        # Apply the lower bound by filtering on the per-class marginal.
        if lb is None:
            return space
        lb_vec = np.atleast_1d(lb).astype(int)
        if lb_vec.size == 1:
            if lb_vec[0] <= 0:
                return space
            lb_vec = np.full(R, int(lb_vec[0]))
        if np.all(lb_vec <= 0):
            return space

        keep = []
        for row in range(space.shape[0]):
            _, nir, _, _ = _toMarginal(sn, node_idx, space[row, :])
            nir = np.atleast_1d(np.asarray(nir).ravel())
            if nir.size >= R and np.all(nir[:R] >= lb_vec):
                keep.append(row)
        return space[keep, :] if keep else np.zeros((0, space.shape[1]))

    @staticmethod
    def isValid(model: 'Network', n: np.ndarray, s: np.ndarray = None,
                options=None) -> bool:
        """
        Validate a network state against capacity and scheduling constraints.

        Args:
            model: Network model (or NetworkStruct).
            n: (nstations x nclasses) matrix of resident jobs per class.
            s: (nstations x nclasses) matrix of running jobs per class
               (optional).
            options: Unused; accepted for API compatibility.

        Returns:
            True if the state satisfies all capacity, server, disabled-process
            and chain-population constraints; False otherwise.

        References:
            MATLAB: matlab/src/lang/+State/isValid.m
        """
        from ..api.state.marginal import toMarginal as _toMarginal
        from .base import SchedStrategy, NodeType
        from ..constants import GlobalConstants
        sn = State._struct_of(model)

        if n is None and s is not None:
            return False

        # If n is given as a list of per-stateful state vectors, reduce it to
        # the (nstations x nclasses) job/running matrices via toMarginal.
        if isinstance(n, (list, tuple)):
            ncell = n
            n = np.zeros((sn.nstations, sn.nclasses))
            s = np.zeros((sn.nstations, sn.nclasses))
            for isf in range(len(ncell)):
                ist = int(sn.statefulToStation[isf])
                ind = int(sn.statefulToNode[isf])
                _, nir, sir, _ = _toMarginal(sn, ind, ncell[isf])
                if ist >= 0:
                    n[ist, :] = np.atleast_1d(np.asarray(nir).ravel())[:sn.nclasses]
                    s[ist, :] = np.atleast_1d(np.asarray(sir).ravel())[:sn.nclasses]

        n = np.atleast_2d(np.asarray(n, dtype=float))
        R = sn.nclasses

        # Disabled-process and per-class capacity checks. A class whose service
        # process is disabled has a NaN service rate (mirrors the NaN in the
        # MATLAB (D0,D1) representation); a state with jobs of such a class at a
        # non-Place station is invalid.
        for ist in range(sn.nstations):
            node_idx = int(sn.stationToNode[ist])
            is_place = sn.nodetype[node_idx] == NodeType.PLACE
            for r in range(R):
                if not is_place and sn.rates is not None \
                        and np.isnan(sn.rates[ist, r]) and n[ist, r] > 0:
                    return False
            if sn.classcap is not None and np.any(n[ist, :] > sn.classcap[ist, :]):
                return False

        # Server/running-job constraints.
        if s is not None and np.size(s) > 0:
            s = np.atleast_2d(np.asarray(s, dtype=float))
            nonpreemptive = {SchedStrategy.FCFS, SchedStrategy.SIRO,
                             SchedStrategy.LCFS, SchedStrategy.HOL,
                             SchedStrategy.POLLING}
            for ist in range(sn.nstations):
                if sn.nservers[ist] > 0:
                    if np.sum(s[ist, :]) > sn.nservers[ist]:
                        sched_ist = sn.sched.get(ist) if isinstance(sn.sched, dict) else sn.sched[ist]
                        if sched_ist in nonpreemptive:
                            return False
                    if np.any(n < s):
                        return False

        # Closed-chain population conservation.
        for nc in range(sn.nchains):
            members = np.where(np.asarray(sn.chains[nc, :]).ravel() != 0)[0]
            njobs_chain = np.sum(np.asarray(sn.njobs).ravel()[members])
            if np.isfinite(njobs_chain):
                statejobs_chain = np.sum(n[:, members])
                if statejobs_chain == 0:
                    if njobs_chain != 0:
                        return False
                elif abs(1 - njobs_chain / statejobs_chain) > GlobalConstants.CoarseTol:
                    return False

        return True

    # snake_case aliases for consistency
    from_marginal = fromMarginal
    from_marginal_and_running = fromMarginalAndRunning
    from_marginal_and_started = fromMarginalAndStarted
    from_marg = fromMarg
    from_marg_and_started = fromMargAndStarted
    to_marginal = toMarginal
    from_marginal_bounds = fromMarginalBounds
    is_valid = isValid
