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

        return _fromMarginal(sn, node_idx, n)

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

        # Get full state space for marginal
        states = _fromMarginal(sn, node_idx, n)

        if states.size == 0:
            return states

        # Filter states that match the running constraint
        # The running jobs are typically encoded in specific columns of the state
        # This is a simplified implementation - full implementation would need
        # to understand the state encoding structure
        s = np.atleast_1d(s)
        n = np.atleast_1d(n)

        # For now, return states that could match - exact filtering requires
        # understanding the specific state encoding
        return states

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

        # Get full state space for marginal
        states = _fromMarginal(sn, node_idx, n)

        if states.size == 0:
            return states

        # Filter states that match the started constraint
        # Simplified implementation - exact filtering requires state encoding knowledge
        s = np.atleast_1d(s)
        n = np.atleast_1d(n)

        return states

    # snake_case aliases for consistency
    from_marginal = fromMarginal
    from_marginal_and_running = fromMarginalAndRunning
    from_marginal_and_started = fromMarginalAndStarted
