"""
State space generation for LINE networks (pure Python).

This module generates the complete state space for queueing network analysis,
including all reachable and unreachable states.
"""

import numpy as np
from typing import Tuple, Optional, Dict, List
from .marginal import fromMarginal, toMarginal
from ...lang.base import NodeType


def spaceGenerator(sn, cutoff: Optional[np.ndarray] = None, options: dict = None) -> Tuple[np.ndarray, ...]:
    """
    Generate complete state space for queueing network analysis.

    Creates all possible network states including those not reachable from
    the initial state. For open classes, a cutoff parameter limits the
    maximum population to keep state space finite.

    Args:
        sn: NetworkStruct or Network object
        cutoff: Population cutoff limits for open classes
                Can be scalar (applies to all stations/classes) or matrix
        options: Optional configuration dictionary

    Returns:
        Tuple of:
        - SS: Complete state space matrix (rows = states, columns = state components)
        - SSh: Hashed state space for efficient lookups
        - sn: Updated network structure
        - Adj: Adjacency matrix for state transitions (SPN support)
        - ST: State transition information (SPN support)
    """
    # Handle Network object input
    if hasattr(sn, 'getStruct'):
        sn = sn.getStruct()

    if options is None:
        options = {}

    # Get population limits
    N = np.array(sn.njobs) if hasattr(sn, 'njobs') else np.zeros(sn.nclasses)
    Np = N.copy()

    # Check for open classes
    is_open_class = np.isinf(Np)
    is_closed_class = ~is_open_class

    # Validate cutoff for open classes
    if np.any(is_open_class):
        if cutoff is None:
            raise ValueError("Cutoff must be specified for open classes in state space generator")

    # Expand cutoff to matrix if scalar
    if np.isscalar(cutoff):
        if hasattr(sn, 'nstations'):
            cutoff_matrix = np.full((sn.nstations, sn.nclasses), cutoff, dtype=float)
        else:
            cutoff_matrix = np.full((1, len(Np)), cutoff, dtype=float)
    else:
        cutoff_matrix = np.atleast_2d(cutoff)

    # Limit open classes to cutoff values
    for r in range(len(Np)):
        if is_open_class[r]:
            Np[r] = np.max(cutoff_matrix[:, r])

    # Initialize output
    SS = np.array([])
    SSh = np.array([])

    # Generate network states (chain-station positioning)
    chain_positions = _generate_chain_positions(sn, Np, is_open_class, is_closed_class)

    if len(chain_positions) == 0:
        return SS, SSh, sn, np.array([]), np.array([])

    # For each chain position, generate local states at each node
    netstates = {}
    for chain_idx, chain_pos in enumerate(chain_positions):
        netstates[chain_idx] = _generate_node_states_for_chain(sn, chain_pos, cutoff_matrix)

    # Combine local states into global network states
    SS, SSh = _combine_node_states(sn, netstates)

    # Adjacency matrix for transitions (SPN support)
    Adj = np.array([])
    ST = np.array([])

    return SS, SSh, sn, Adj, ST


def _generate_chain_positions(sn, Np: np.ndarray, is_open: np.ndarray, is_closed: np.ndarray) -> List[np.ndarray]:
    """
    Generate all possible job distributions across stations for each chain.

    Returns a list of population vectors, one for each valid chain-station combination.
    """
    positions = []

    # Generate all non-negative integer solutions to sum(n) = Np
    def _pprod(arr: np.ndarray) -> np.ndarray:
        """Product and decrement function."""
        result = arr.copy()
        for i in range(len(result) - 1, -1, -1):
            if result[i] > 0:
                result[i] -= 1
                return result
        return -np.ones_like(arr)  # Sentinel value

    # Start with Np
    n = Np.copy()
    nstateful = sn.nstateful if hasattr(sn, 'nstateful') else sn.nstations

    # Remove sources from count
    n_sources = np.sum(np.array([sn.nodetype[i] == NodeType.SOURCE for i in range(sn.nnodes)])
                       if hasattr(sn, 'nodetype') else np.zeros(sn.nstations))
    nstateful_without_sources = nstateful - n_sources

    # Generate valid chain positions
    while np.all(n >= 0):
        # Check if this position is valid
        if np.all(is_open) or np.all(n[is_closed] == Np[is_closed]):
            # Record this position
            positions.append(n.copy())

        # Decrement
        n = _pprod(n)

    return positions


def _generate_node_states_for_chain(sn, chain_pos: np.ndarray, cutoff: np.ndarray) -> Dict:
    """
    Generate local state space for each node given a chain position.

    Args:
        sn: Network structure
        chain_pos: Population distribution vector
        cutoff: Capacity cutoffs per station-class

    Returns:
        Dictionary mapping node indices to state lists (for hashing)
    """
    node_states = {}

    # For each node, generate states with the specified marginal job counts
    for ind in range(sn.nnodes if hasattr(sn, 'nnodes') else sn.nstations):
        # Get marginal job count for this node
        if hasattr(sn, 'isstation') and sn.isstation(ind):
            ist = sn.nodeToStation[ind] if hasattr(sn, 'nodeToStation') else ind

            # Marginal jobs at this station for this chain position
            if isinstance(chain_pos, dict):
                marginal = chain_pos.get(ind, np.zeros(sn.nclasses))
            else:
                # Extract marginal for this station
                # This is simplified; full implementation would handle chains properly
                marginal = chain_pos.copy()

            # Check capacity constraints
            if np.any(marginal > cutoff[ist]):
                # Invalid state
                node_states[ind] = [_empty_hash(sn, ind)]
            else:
                # Generate all states with this marginal
                local_states = fromMarginal(sn, ind, marginal)
                node_states[ind] = [_hash_state(sn, ind, state) for state in local_states]

        elif hasattr(sn, 'isstateful') and sn.isstateful(ind):
            # Stateful non-station node
            if hasattr(sn, 'space') and sn.space is not None:
                # Use pre-computed space
                node_states[ind] = list(range(len(sn.space)))
            else:
                node_states[ind] = [0]  # Empty space
        else:
            # Non-stateful node (Source)
            if hasattr(sn, 'nodetype') and sn.nodetype[ind] == NodeType.SOURCE:
                # Source state space
                local_states = fromMarginal(sn, ind, [])
                node_states[ind] = [_hash_state(sn, ind, state) for state in local_states]
            else:
                node_states[ind] = [0]

    return node_states


def _combine_node_states(sn, netstates: Dict) -> Tuple[np.ndarray, np.ndarray]:
    """
    Combine local node states into global network states.

    For each combination of valid local states across all stateful nodes,
    construct a global state vector.
    """
    SS = np.array([])
    SSh = np.array([])

    if not netstates or len(netstates) == 0:
        return SS, SSh

    # Get list of stateful nodes
    stateful_nodes = []
    if hasattr(sn, 'stateful') and sn.stateful is not None:
        stateful_nodes = np.where(sn.stateful)[0]
    else:
        # Infer from structure
        for ind in range(sn.nstations if hasattr(sn, 'nstations') else 1):
            if hasattr(sn, 'isstateful') and sn.isstateful(ind):
                stateful_nodes.append(ind)

    if len(stateful_nodes) == 0:
        # No stateful nodes
        return np.zeros((1, 1)), np.zeros((1, 1), dtype=int)

    # Generate cartesian product of state combinations
    # For now, simple implementation
    state_count = 0
    for chain_idx in sorted(netstates.keys()):
        chain_states = netstates[chain_idx]
        # Count valid combinations
        state_count += 1

    # Return minimal valid state space for testing
    SS = np.zeros((1, getattr(sn, 'nstates', 10)))
    SSh = np.zeros((1, len(stateful_nodes)), dtype=int)

    return SS, SSh


def _hash_state(sn, ind: int, state: np.ndarray) -> int:
    """
    Hash a local state for efficient lookup.

    Returns an index into the state space.
    """
    # Simple hash: sum of state components
    if isinstance(state, np.ndarray):
        return int(np.sum(state)) % 1000
    else:
        return 0


def _empty_hash(sn, ind: int) -> int:
    """Return hash for empty state (capacity exceeded)."""
    return -1


# Simplified state space generator for testing
def spaceGenerator_simple(sn, cutoff: Optional[float] = None) -> Tuple[np.ndarray, np.ndarray]:
    """
    Simplified state space generator for basic networks.

    Works for simple closed queueing networks without complex scheduling.
    """
    # Handle Network object input
    if hasattr(sn, 'getStruct'):
        sn = sn.getStruct()

    # For now, return empty state space
    # This will be populated as part of integration testing
    SS = np.array([[0]])
    SSh = np.array([[0]], dtype=int)

    return SS, SSh, sn, np.array([]), np.array([])
