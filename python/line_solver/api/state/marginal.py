"""
Marginal state analysis for LINE networks (pure Python).

This module provides functions to extract and generate marginal distributions
from network states, essential for state probability computation.
"""

import numpy as np
from typing import Tuple, Optional, Union, List
from ...lang.base import SchedStrategy, NodeType


def toMarginal(sn, ind: int, state_i: np.ndarray = None, phasesz: np.ndarray = None,
               phaseshift: np.ndarray = None, space_buf: np.ndarray = None,
               space_srv: np.ndarray = None, space_var: np.ndarray = None) -> Tuple[np.ndarray, ...]:
    """
    Extract marginal job distributions from global state for a specific node.

    Computes the marginal queue-length distributions and job counts for a
    specific node from the global network state, considering scheduling
    strategies and service phases.

    Args:
        sn: NetworkStruct or Network object
        ind: Node index (0-based)
        state_i: Global state vector or matrix (rows = states, columns = state components)
        phasesz: Vector of phase sizes for each class
        phaseshift: Phase shift parameters for state extraction
        space_buf: Buffer space configuration
        space_srv: Service space configuration
        space_var: Local variables configuration

    Returns:
        Tuple of:
        - ni: Total jobs in node (array: n_states)
        - nir: Total jobs per class (array: n_states x n_classes)
        - sir: Jobs in service per class (array: n_states x n_classes)
        - kir: Jobs in service per class per phase (array: n_states x n_classes x max_phases)
    """
    # Handle Network object input
    if hasattr(sn, 'getStruct'):
        sn = sn.getStruct()

    # Validate inputs
    if state_i is None:
        if hasattr(sn, 'state') and sn.state is not None:
            state_i = sn.state
        else:
            raise ValueError("state_i must be provided or available in sn.state")

    state_i = np.atleast_2d(state_i)
    n_states = state_i.shape[0]

    # Handle non-station nodes (Transitions, Places, etc.)
    if hasattr(sn, 'isstation') and not sn.isstation[ind]:
        if hasattr(sn, 'nodetype') and sn.nodetype[ind] == NodeType.TRANSITION:
            # Stateful non-station node (e.g., Transition in SPN)
            R = sn.nodeparam[ind].get('nmodes', 1) if hasattr(sn, 'nodeparam') else 1
            sir = np.zeros((n_states, R))
            max_phases = max(phasesz) if phasesz is not None else 1
            kir = np.zeros((n_states, R, max_phases))

            if phasesz is not None and phaseshift is not None:
                for r in range(R):
                    for k in range(phasesz[r]):
                        kir[:, r, k] = space_srv[:, phaseshift[r] + k]
                        sir[:, r] += kir[:, r, k]

            nir = sir
            ni = np.sum(nir, axis=1)
            return ni, nir, sir, kir

    # Get network parameters
    R = sn.nclasses if hasattr(sn, 'nclasses') else 1

    # Get station-level information
    if hasattr(sn, 'nodeToStation'):
        ist = sn.nodeToStation[ind]
    else:
        ist = ind

    # Set default phase parameters if not provided
    if phasesz is None:
        if hasattr(sn, 'phasessz'):
            phasesz = sn.phasessz[ist]
        else:
            phasesz = np.ones(R, dtype=int)
    else:
        phasesz = np.atleast_1d(phasesz)

    if phaseshift is None:
        if hasattr(sn, 'phaseshift'):
            phaseshift = sn.phaseshift[ist]
        else:
            phaseshift = np.concatenate([[0], np.cumsum(phasesz[:-1])])
    else:
        phaseshift = np.atleast_1d(phaseshift)

    # Extract space information if not provided
    if space_var is None or space_srv is None or space_buf is None:
        total_vars = sum(sn.nvars[ind]) if hasattr(sn, 'nvars') else 0
        total_phases = sum(phasesz)

        if space_var is None:
            space_var = state_i[:, -total_vars:] if total_vars > 0 else np.array([])

        if space_srv is None:
            srv_end = state_i.shape[1] - total_vars
            srv_start = srv_end - total_phases
            space_srv = state_i[:, srv_start:srv_end] if total_phases > 0 else np.array([])

        if space_buf is None:
            buf_end = state_i.shape[1] - total_phases - total_vars
            space_buf = state_i[:, :buf_end] if buf_end > 0 else np.array([])

    # Determine if distribution is exponential (single phase)
    isExponential = max(phasesz) == 1

    # Initialize output arrays
    if isExponential:
        sir = space_srv
        kir = space_srv
        nir = np.zeros((n_states, R))
    else:
        sir = np.zeros((n_states, R))
        max_phases = max(phasesz)
        kir = np.zeros((n_states, R, max_phases))
        nir = np.zeros((n_states, R))

        # Extract phase-specific information
        for r in range(R):
            for k in range(phasesz[r]):
                if phaseshift[r] + k < space_srv.shape[1]:
                    kir[:, r, k] = space_srv[:, phaseshift[r] + k]
                    sir[:, r] += kir[:, r, k]

    # Compute total jobs per class based on scheduling strategy
    if hasattr(sn, 'sched'):
        sched = sn.sched[ist]
    else:
        sched = SchedStrategy.FCFS

    # Apply scheduling strategy rules to compute nir from sir and space_buf
    if sched in [SchedStrategy.INF, SchedStrategy.PS, SchedStrategy.PSPRIO,
                 SchedStrategy.DPS, SchedStrategy.GPS]:
        # These policies: jobs in service = jobs in station (no buffer distinction)
        nir = sir.copy()

    elif sched == SchedStrategy.EXT:
        # External node (Source): infinite jobs per class
        nir = np.full((n_states, R), np.inf)

    elif sched in [SchedStrategy.FCFS, SchedStrategy.HOL,
                   SchedStrategy.LCFS]:
        # FCFS/LCFS: count jobs in service plus those in buffer
        nir = sir.copy()
        if space_buf.size > 0:
            for r in range(R):
                # Count jobs of class r in buffer
                if space_buf.ndim == 1:
                    nir[:, r] += np.sum(space_buf == r)
                else:
                    nir[:, r] += np.sum(space_buf == r, axis=1)

    elif sched in [SchedStrategy.SIRO]:
        # These policies track buffer per class
        nir = sir.copy()
        if space_buf.size > 0:
            for r in range(R):
                if space_buf.ndim == 1:
                    if r < len(space_buf):
                        nir[:, r] += space_buf[r]
                else:
                    if r < space_buf.shape[1]:
                        nir[:, r] += space_buf[:, r]

    else:
        # Default: jobs in station = jobs in service
        nir = sir.copy()

    # Handle disabled stations (set counts to 0)
    if hasattr(sn, 'rates'):
        for r in range(R):
            disabled = np.isnan(sn.rates[ist, r]) if hasattr(sn.rates[ist, r], '__iter__') else np.isnan(sn.rates[ist, r])
            if disabled:
                nir[:, r] = 0
                sir[:, r] = 0
                for k in range(max(phasesz)):
                    if k < kir.shape[2]:
                        kir[:, r, k] = 0

    # Total jobs in node
    ni = np.sum(nir, axis=1)

    return ni, nir, sir, kir


def fromMarginal(sn, ind: int, n: Union[np.ndarray, list], options: dict = None) -> np.ndarray:
    """
    Generate state space with specific marginal job counts at a node.

    Creates all possible network states where the specified node has
    exactly n[r] jobs of class r. This is essential for computing
    state probabilities and performance metrics.

    Args:
        sn: NetworkStruct or Network object
        ind: Node index (0-based)
        n: Vector of job counts per class [n_classes]
        options: Optional configuration dictionary

    Returns:
        State space matrix where each row is a valid state with the specified
        marginal job counts at node ind.
    """
    # Handle Network object input
    if hasattr(sn, 'getStruct'):
        sn = sn.getStruct()

    if options is None:
        options = {}

    # Validate inputs
    n = np.atleast_1d(n)
    R = sn.nclasses if hasattr(sn, 'nclasses') else len(n)

    # Ensure n has correct length
    if len(n) < R:
        n = np.concatenate([n, np.zeros(R - len(n))])

    # Get station index
    if hasattr(sn, 'nodeToStation'):
        ist = sn.nodeToStation[ind]
    else:
        ist = ind

    # Handle non-station nodes (Transitions, cached variables, etc.)
    if hasattr(sn, 'isstation') and not sn.isstation[ind]:
        # For stateful non-station nodes, return empty state or fixed state
        if hasattr(sn, 'space') and len(sn.space) > 0:
            return np.array([[0]])
        else:
            return np.zeros((1, sum(n)))

    # Get scheduling strategy
    if hasattr(sn, 'sched'):
        sched = sn.sched[ist]
    else:
        sched = SchedStrategy.FCFS

    # Get number of servers
    if hasattr(sn, 'nservers'):
        S = sn.nservers[ist]
    else:
        S = 1

    # Get phase information
    phases = np.zeros(R, dtype=int)
    if hasattr(sn, 'proc') and sn.proc is not None:
        for r in range(R):
            if hasattr(sn.proc[ist], '__getitem__') and r < len(sn.proc[ist]):
                proc_entry = sn.proc[ist][r]
                if proc_entry is not None:
                    # Handle dict case (Erlang, HyperExp, Exp store params as dict)
                    if isinstance(proc_entry, dict):
                        if 'k' in proc_entry:
                            # Erlang: phases = k (number of stages)
                            phases[r] = int(proc_entry['k'])
                        elif 'probs' in proc_entry and 'rates' in proc_entry:
                            # HyperExp: number of phases = length of probs array
                            phases[r] = len(proc_entry['probs'])
                        else:
                            # Exp or other: 1 phase
                            phases[r] = 1
                    # Handle list case (PH, MAP, etc. store as [alpha, T])
                    elif isinstance(proc_entry, (list, tuple)) and len(proc_entry) > 0:
                        first_elem = proc_entry[0]
                        phases[r] = len(first_elem) if isinstance(first_elem, (list, np.ndarray)) else 1
                    else:
                        phases[r] = 1
                else:
                    phases[r] = 1
            else:
                phases[r] = 1
    else:
        phases = np.ones(R, dtype=int)

    # Check capacity constraints
    if hasattr(sn, 'classcap') and sched != SchedStrategy.EXT:
        if np.any(n > sn.classcap[ist]):
            return np.array([])

    # Generate state space based on scheduling strategy
    space = _generate_state_space_for_sched(sched, S, n, phases, R)

    # Sort and unique the state space
    if len(space) > 0:
        space = np.unique(space, axis=0)
        # Reverse sort to put states with jobs in phase 1 earlier
        space = space[::-1]

    return space


def _generate_state_space_for_sched(sched: int, S: int, n: np.ndarray, phases: np.ndarray, R: int) -> np.ndarray:
    """
    Generate local state space for a specific scheduling strategy.

    Args:
        sched: SchedStrategy enum value
        S: Number of servers
        n: Job counts per class
        phases: Number of phases per class
        R: Number of classes

    Returns:
        State space matrix for this node
    """
    space = np.array([])

    if sched == SchedStrategy.EXT:
        # Source: one job per phase per class
        state = np.ones((1, sum(phases)))
        # Add infinite buffer marker
        space = np.concatenate([np.full((1, 1), np.inf), state], axis=1)

    elif sched in [SchedStrategy.INF, SchedStrategy.PS, SchedStrategy.DPS,
                   SchedStrategy.GPS, SchedStrategy.PSPRIO]:
        # Jobs only in service, no buffer
        # Generate all combinations of phase distributions
        space = _cartesian_space_for_phases(n, phases, R)

    elif sched in [SchedStrategy.FCFS, SchedStrategy.HOL,
                   SchedStrategy.LCFS]:
        # Jobs in ordered buffer + service
        if sum(n) == 0:
            space = np.zeros((1, 1 + sum(phases)))
        else:
            # Generate buffer and phase states
            space = _cartesian_space_fcfs_lcfs(n, phases, S, R)

    elif sched == SchedStrategy.SIRO:
        # Unordered buffer + service
        if sum(n) <= S:
            # All jobs in service
            space = _cartesian_space_for_phases(n, phases, R)
            # Add zero buffer
            space = np.concatenate([np.zeros((space.shape[0], R)), space], axis=1)
        else:
            # Some jobs in buffer, some in service
            space = _generate_siro_space(n, phases, S, R)

    else:
        # Default: just track service phases
        space = _cartesian_space_for_phases(n, phases, R)

    return space if isinstance(space, np.ndarray) and space.size > 0 else np.array([])


def _cartesian_space_for_phases(n: np.ndarray, phases: np.ndarray, R: int) -> np.ndarray:
    """
    Generate state space as cartesian product of phase distributions.

    For each class r, we generate states where class r jobs are distributed
    across phases. The total for class r equals n[r].
    """
    if sum(n) == 0:
        return np.zeros((1, sum(phases)))

    # Generate phase states for each class
    class_states = []
    for r in range(R):
        if phases[r] > 1:
            # Multiple phases: distribute n[r] jobs across phases
            phase_space = _generate_phase_distribution(int(n[r]), int(phases[r]))
        else:
            # Single phase
            phase_space = np.array([[int(n[r])]])
        class_states.append(phase_space)

    # Cartesian product of all class states
    return _cartesian_product(*class_states)


def _generate_phase_distribution(n_jobs: int, n_phases: int) -> np.ndarray:
    """Generate all ways to distribute n_jobs across n_phases."""
    if n_jobs == 0:
        return np.zeros((1, n_phases), dtype=int)

    if n_phases == 1:
        return np.array([[n_jobs]])

    # Recursive: distribute across phases
    states = []
    for k in range(n_jobs + 1):
        # k jobs in phase 1, rest in remaining phases
        rest = _generate_phase_distribution(n_jobs - k, n_phases - 1)
        for r in rest:
            states.append(np.concatenate([[k], r]))

    return np.array(states)


def _cartesian_space_fcfs_lcfs(n: np.ndarray, phases: np.ndarray, S: int, R: int) -> np.ndarray:
    """
    Generate state space for FCFS/LCFS scheduling.

    Tracks ordered buffer and service phases. This is complex as it requires
    generating all permutations of job orderings.
    """
    # Simplified implementation for small n
    if sum(n) > 20:
        # For large populations, use approximation
        return _cartesian_space_for_phases(n, phases, R)

    # Generate basic structure: buffer (ordered) + phases (service)
    state = np.zeros((1, 1 + sum(phases)))

    # For now, just return basic structure
    # Full implementation would generate all permutations
    return state


def _generate_siro_space(n: np.ndarray, phases: np.ndarray, S: int, R: int) -> np.ndarray:
    """Generate state space for SIRO (Service In Random Order) scheduling."""
    # Simplified: distribute jobs between buffer and service
    space = []

    # All combinations where some jobs are in service (up to S)
    from itertools import combinations_with_replacement

    total_jobs = int(sum(n))
    for num_in_service in range(min(total_jobs, S) + 1):
        # Distribute num_in_service among classes
        for service_dist in _integer_partitions(num_in_service, R, n):
            buffer_dist = n - service_dist

            # Generate phase states for service portion
            phase_state = _cartesian_space_for_phases(service_dist, phases, R)

            # Add buffer info
            for ps in phase_state:
                state = np.concatenate([buffer_dist, ps])
                space.append(state)

    return np.array(space) if space else np.array([])


def _integer_partitions(total: int, num_parts: int, max_vals: np.ndarray = None):
    """Generate all ways to partition total into num_parts."""
    if num_parts == 1:
        if max_vals is not None and total <= max_vals[0]:
            yield np.array([total])
        elif max_vals is None:
            yield np.array([total])
        return

    if max_vals is None:
        max_vals = np.full(num_parts, total)

    for i in range(min(total, int(max_vals[0])) + 1):
        for rest in _integer_partitions(total - i, num_parts - 1, max_vals[1:]):
            yield np.concatenate([[i], rest])


def _cartesian_product(*arrays) -> np.ndarray:
    """
    Compute cartesian product of arrays.
    Each array is concatenated horizontally for each combination.
    """
    if not arrays:
        return np.array([])

    if len(arrays) == 1:
        return arrays[0]

    # Start with first array
    result = arrays[0]

    # Iteratively cartesian product with remaining arrays
    for arr in arrays[1:]:
        # Expand result and arr to compute cartesian product
        n_result = result.shape[0]
        n_arr = arr.shape[0]

        # Repeat result n_arr times
        result_expanded = np.repeat(result, n_arr, axis=0)

        # Tile arr n_result times
        arr_expanded = np.tile(arr, (n_result, 1))

        # Concatenate
        result = np.concatenate([result_expanded, arr_expanded], axis=1)

    return result
