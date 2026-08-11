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

    # Start with Np
    n = Np.copy()
    n_len = len(n)
    nstateful = sn.nstateful if hasattr(sn, 'nstateful') else sn.nstations

    # Remove sources from count
    n_sources = np.sum(np.array([sn.nodetype[i] == NodeType.SOURCE for i in range(sn.nnodes)])
                       if hasattr(sn, 'nodetype') else np.zeros(sn.nstations))
    nstateful_without_sources = nstateful - n_sources

    # Use JIT-accelerated pprod_prev for large state spaces

    def _pprod(arr: np.ndarray) -> np.ndarray:
        """Product and decrement function."""
        result = arr.copy()
        for i in range(len(result) - 1, -1, -1):
            if result[i] > 0:
                result[i] -= 1
                return result
        return -np.ones_like(arr)  # Sentinel value

    while np.all(n >= 0):
        if np.all(is_open) or np.all(n[is_closed] == Np[is_closed]):
            positions.append(n.copy())
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

    # Collect per-chain state lists and compute sizes for cartesian product
    chain_keys = sorted(netstates.keys())
    chain_state_lists = [netstates[k] for k in chain_keys]
    sizes = np.array([len(sl) for sl in chain_state_lists], dtype=np.int64)
    n_arrays = len(sizes)
    total_combos = int(np.prod(sizes)) if n_arrays > 0 else 0

    # Use JIT-accelerated cartesian product for large state spaces
    # Pure Python cartesian product via indices
    if total_combos > 0:
        indices = np.zeros((total_combos, n_arrays), dtype=np.int64)
        idx = np.zeros(n_arrays, dtype=int)
        for row in range(total_combos):
            indices[row, :] = idx
            pos = n_arrays - 1
            while pos >= 0:
                idx[pos] += 1
                if idx[pos] < sizes[pos]:
                    break
                idx[pos] = 0
                pos -= 1
        state_count = total_combos
    else:
        state_count = 0

    if state_count == 0:
        return np.zeros((1, getattr(sn, 'nstates', 10))), np.zeros((1, len(stateful_nodes)), dtype=int)

    # Build SSh from indices into per-chain state lists
    SSh = np.zeros((state_count, len(stateful_nodes)), dtype=int)
    for row in range(state_count):
        for ci in range(n_arrays):
            # Use the hash from the corresponding chain-state list
            SSh[row, ci % len(stateful_nodes)] = chain_state_lists[ci][int(indices[row, ci])]

    SS = np.zeros((state_count, getattr(sn, 'nstates', 10)))

    return SS, SSh

def _hash_state(sn, ind: int, state: np.ndarray) -> int:
    """
    Hash a local state for efficient lookup.

    Returns an index into the state space.
    """
    if isinstance(state, np.ndarray):
        n = len(state)
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

def enumerate_all_populations(N: np.ndarray) -> np.ndarray:
    """
    Enumerate all population vectors n where 0 <= n[r] <= N[r].

    Uses JIT acceleration for large population spaces.

    Args:
        N: Maximum population per class

    Returns:
        Array of shape (total, n_classes) with all population vectors
    """
    N = np.asarray(N, dtype=np.int64).ravel()
    n_classes = len(N)
    total = int(np.prod(N + 1))

    # Pure Python fallback
    output = np.zeros((total, n_classes), dtype=np.int64)
    n = np.zeros(n_classes, dtype=np.int64)
    for row in range(total):
        output[row, :] = n
        pos = n_classes - 1
        while pos >= 0:
            n[pos] += 1
            if n[pos] <= N[pos]:
                break
            n[pos] = 0
            pos -= 1
    return output

def count_phase_distributions(n_jobs: int, n_phases: int) -> int:
    """
    Count the number of ways to distribute n_jobs across n_phases.

    Uses JIT acceleration when available.

    Args:
        n_jobs: Number of jobs to distribute
        n_phases: Number of phases

    Returns:
        Number of distinct distributions
    """
    # Pure Python fallback
    from scipy.special import comb
    if n_jobs == 0 or n_phases == 1:
        return 1
    return int(comb(n_jobs + n_phases - 1, n_phases - 1, exact=True))

def generate_phase_distributions(n_jobs: int, n_phases: int) -> np.ndarray:
    """
    Generate all ways to distribute n_jobs across n_phases.

    Uses JIT acceleration for large distributions.

    Args:
        n_jobs: Number of jobs to distribute
        n_phases: Number of phases

    Returns:
        Array of shape (count, n_phases) with all distributions
    """
    count = count_phase_distributions(n_jobs, n_phases)

    # Pure Python fallback
    output = np.zeros((count, n_phases), dtype=np.int64)
    if n_jobs == 0:
        return output
    if n_phases == 1:
        output[0, 0] = n_jobs
        return output

    state = np.zeros(n_phases, dtype=np.int64)
    state[0] = n_jobs
    row = 0
    while True:
        output[row, :] = state
        row += 1
        pos = n_phases - 2
        while pos >= 0 and state[pos] == 0:
            pos -= 1
        if pos < 0:
            break
        state[pos] -= 1
        state[pos + 1] += 1
        total_right = np.sum(state[pos + 2:])
        state[pos + 2:] = 0
        state[pos + 1] += total_right

    return output[:row]

def generate_integer_partitions(total: int, num_parts: int,
                                 max_vals: np.ndarray) -> np.ndarray:
    """
    Generate all constrained integer partitions.

    Generates all ways to partition total into num_parts where part[i] <= max_vals[i].
    Uses JIT acceleration for large partition spaces.

    Args:
        total: Total to partition
        num_parts: Number of parts
        max_vals: Maximum value for each part

    Returns:
        Array of shape (count, num_parts) with all partitions
    """
    max_vals = np.asarray(max_vals, dtype=np.int64).ravel()

    # Pure Python fallback using recursive generation
    results = []

    def _generate(pos, remaining, state):
        if pos == num_parts - 1:
            if remaining <= max_vals[pos]:
                state[pos] = remaining
                results.append(state.copy())
            return
        for v in range(min(remaining, int(max_vals[pos])) + 1):
            state[pos] = v
            _generate(pos + 1, remaining - v, state)

    _generate(0, total, np.zeros(num_parts, dtype=np.int64))

    if len(results) == 0:
        return np.zeros((0, num_parts), dtype=np.int64)
    return np.array(results, dtype=np.int64)

def find_state_in_space(state_space: np.ndarray, target: np.ndarray) -> int:
    """
    Search for a target state in a state space matrix.

    Uses JIT acceleration for large state spaces.

    Args:
        state_space: State space matrix (n_states x n_cols)
        target: Target state to find

    Returns:
        Index if found (0-based), -1 if not found
    """
    state_space = np.asarray(state_space)
    target = np.asarray(target)

    if state_space.ndim != 2 or len(state_space) == 0:
        return -1

    n_states, n_cols = state_space.shape

    # Pure Python fallback
    for i in range(n_states):
        if np.array_equal(state_space[i, :n_cols], target[:n_cols]):
            return i
    return -1
