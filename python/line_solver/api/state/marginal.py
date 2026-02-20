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
        if hasattr(sn, 'phasessz') and sn.phasessz is not None:
            phasesz = sn.phasessz[ist]
        else:
            phasesz = np.ones(R, dtype=int)
    else:
        phasesz = np.atleast_1d(phasesz)

    if phaseshift is None:
        if hasattr(sn, 'phaseshift') and sn.phaseshift is not None:
            phaseshift = sn.phaseshift[ist]
        else:
            phaseshift = np.concatenate([[0], np.cumsum(phasesz[:-1])])
    else:
        phaseshift = np.atleast_1d(phaseshift)

    # Extract space information if not provided
    if space_var is None or space_srv is None or space_buf is None:
        total_vars = sum(sn.nvars[ind]) if hasattr(sn, 'nvars') and sn.nvars is not None else 0
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
        # Buffer entries use 1-based class indices (MATLAB convention):
        # 0 = empty, 1 = class 0, 2 = class 1, etc.
        nir = sir.copy()
        if space_buf.size > 0:
            for r in range(R):
                # Count jobs of class r in buffer (1-based: class r stored as r+1)
                if space_buf.ndim == 1:
                    nir[:, r] += np.sum(space_buf == (r + 1))
                else:
                    nir[:, r] += np.sum(space_buf == (r + 1), axis=1)

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


# Maximum number of unique permutations to generate before falling back to
# a single representative state. Prevents combinatorial explosion for large queues.
_MAX_UNIQUE_PERMS = 5000


def _uniqueperms(vec):
    """
    Generate all unique permutations of a vector with possibly replicate elements.

    Efficient multiset permutation algorithm matching MATLAB's uniqueperms.
    Avoids the O(n!) overhead of set(permutations(vec)) by directly generating
    only unique permutations.

    For single-value vectors (e.g., [1,1,1,...]), returns immediately with one row.
    For large multisets where the number of unique permutations exceeds _MAX_UNIQUE_PERMS,
    returns a single representative (sorted) permutation to prevent combinatorial explosion.

    Args:
        vec: List or 1D array of elements (possibly with repeats)

    Returns:
        List of tuples, each a unique permutation
    """
    from math import factorial

    if len(vec) == 0:
        return []

    # Count multiplicities
    counts = {}
    for v in vec:
        counts[v] = counts.get(v, 0) + 1

    unique_vals = sorted(counts.keys())

    # Single unique value: only one permutation
    if len(unique_vals) == 1:
        return [tuple(vec)]

    # Estimate number of unique permutations: n! / (a1! * a2! * ... * ak!)
    n = len(vec)
    try:
        num_perms = factorial(n)
        for c in counts.values():
            num_perms //= factorial(c)
    except (OverflowError, ValueError):
        num_perms = _MAX_UNIQUE_PERMS + 1

    if num_perms > _MAX_UNIQUE_PERMS:
        # Too many permutations — return single sorted representative
        return [tuple(sorted(vec))]

    # Generate unique permutations recursively (matching MATLAB uniqueperms algorithm)
    return _uniqueperms_recursive(sorted(vec))


def _uniqueperms_recursive(vec):
    """Recursive multiset permutation generator."""
    if len(vec) <= 1:
        return [tuple(vec)]

    result = []
    seen = set()
    for i, val in enumerate(vec):
        if val in seen:
            continue
        seen.add(val)
        rest = vec[:i] + vec[i+1:]
        for perm in _uniqueperms_recursive(rest):
            result.append((val,) + perm)

    return result


def _cartesian_space_fcfs_lcfs(n: np.ndarray, phases: np.ndarray, S: int, R: int) -> np.ndarray:
    """
    Generate state space for FCFS/LCFS scheduling.

    Tracks ordered buffer and service phases. The state format matches MATLAB:
    - Buffer positions: class ID of job in each buffer position (0=empty)
    - Server phases: phase distribution per class for jobs in service

    Args:
        n: Array of job counts per class
        phases: Array of number of phases per class
        S: Number of servers
        R: Number of classes

    Returns:
        State space matrix where each row is [buf_positions..., phase_counts...]
    """
    total_jobs = int(sum(n))

    if total_jobs == 0:
        # Empty state: 1 buffer column (0) + phase columns
        return np.zeros((1, 1 + int(sum(phases))))

    # Build list of job classes with repetition
    # e.g., n=[2,1] -> vi=[1,1,2] (using 1-based class IDs)
    vi = []
    for r in range(R):
        vi.extend([r + 1] * int(n[r]))  # 1-based class IDs

    # Generate all unique permutations of job ordering
    # Uses efficient multiset permutation algorithm (matching MATLAB's uniqueperms)
    # instead of set(permutations(vi)) which is O(n!) even with duplicates
    mi = _uniqueperms(vi)
    mi = np.array(mi) if len(mi) > 0 else np.array([vi])

    if len(mi) == 0:
        return np.zeros((1, 1 + int(sum(phases))))

    # Build states for each permutation
    all_states = []

    for perm_idx in range(mi.shape[0]):
        perm = mi[perm_idx]

        # mi_buf: class of job in buffer position i (jobs not yet in service)
        # mi_srv: class of job in server (jobs being served)
        num_in_service = int(min(total_jobs, S))
        num_in_buffer = int(total_jobs - num_in_service)

        # Buffer: first (total_jobs - S) positions
        # Server: last S positions
        if num_in_buffer > 0:
            mi_buf = perm[:num_in_buffer]
        else:
            mi_buf = np.array([0])  # Empty buffer marker

        mi_srv = perm[max(0, num_in_buffer):]

        # Count jobs of each class in service: si[r] = number of class (r+1) jobs in server
        si = np.zeros(R, dtype=int)
        for class_id in mi_srv:
            si[int(class_id) - 1] += 1  # Convert to 0-based index

        # Generate phase distributions for jobs in service
        # kstate = Cartesian product of phase distributions per class
        kstate = _cartesian_space_for_phases(si, phases, R)

        # Build full states: [mi_buf, kstate]
        for ks in kstate:
            state = np.concatenate([mi_buf, ks])
            all_states.append(state)

    if not all_states:
        return np.zeros((1, 1 + int(sum(phases))))

    # Stack and ensure consistent column count
    space = np.array(all_states)

    # Sort and remove duplicates
    space = np.unique(space, axis=0)

    return space


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


def fromMarginalAndRunning(sn, ind: int, n: Union[np.ndarray, list],
                           s: Union[np.ndarray, list], options: dict = None) -> np.ndarray:
    """
    Generate state space with specific marginal and running job counts.

    Creates states where node has n[r] jobs of class r total,
    with s[r] jobs of class r currently in service (running).

    Args:
        sn: NetworkStruct or Network object
        ind: Node index (0-based)
        n: Vector of total job counts per class
        s: Vector of running job counts per class
        options: Optional configuration dictionary

    Returns:
        State space matrix where each row is a valid state.
    """
    # Handle Network object input
    if hasattr(sn, 'getStruct'):
        sn = sn.getStruct()

    if options is None:
        options = {'force': False}

    n = np.atleast_1d(n).astype(int)
    s = np.atleast_1d(s).astype(int)

    R = sn.nclasses if hasattr(sn, 'nclasses') else len(n)

    # Ensure arrays have correct length
    if len(n) < R:
        n = np.concatenate([n, np.zeros(R - len(n), dtype=int)])
    if len(s) < R:
        s = np.concatenate([s, np.zeros(R - len(s), dtype=int)])

    # Get station index
    if hasattr(sn, 'nodeToStation'):
        ist = sn.nodeToStation[ind]
    else:
        ist = ind

    # Get number of servers
    if hasattr(sn, 'nservers'):
        S = int(sn.nservers[ist])
    else:
        S = 1

    # Get phase information (same logic as fromMarginal)
    phases = np.zeros(R, dtype=int)
    if hasattr(sn, 'proc') and sn.proc is not None:
        for r in range(R):
            if hasattr(sn.proc[ist], '__getitem__') and r < len(sn.proc[ist]):
                proc_entry = sn.proc[ist][r]
                if proc_entry is not None:
                    if isinstance(proc_entry, dict):
                        if 'k' in proc_entry:
                            phases[r] = int(proc_entry['k'])
                        elif 'probs' in proc_entry and 'rates' in proc_entry:
                            phases[r] = len(proc_entry['probs'])
                        else:
                            phases[r] = 1
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
    if hasattr(sn, 'classcap'):
        if np.any(n > sn.classcap[ist]):
            return np.array([])

    # Check running constraint: running jobs cannot exceed servers
    if S > 0 and sum(s) > S:
        return np.array([])

    # Get scheduling strategy
    if hasattr(sn, 'sched'):
        sched = sn.sched[ist]
    else:
        sched = SchedStrategy.FCFS

    # Generate state space based on scheduling strategy
    space = _generate_state_space_with_running(sched, S, n, s, phases, R)

    # Sort and unique the state space
    if len(space) > 0:
        space = np.unique(space, axis=0)
        # Reverse sort to put states with jobs in phase 1 earlier
        space = space[::-1]

    return space


def _generate_state_space_with_running(sched: int, S: int, n: np.ndarray, s: np.ndarray,
                                       phases: np.ndarray, R: int) -> np.ndarray:
    """
    Generate local state space with specific running job counts.

    Args:
        sched: SchedStrategy enum value
        S: Number of servers
        n: Total job counts per class
        s: Running job counts per class
        phases: Number of phases per class
        R: Number of classes

    Returns:
        State space matrix for this node
    """
    from itertools import permutations

    if sched in [SchedStrategy.INF, SchedStrategy.PS, SchedStrategy.DPS,
                 SchedStrategy.GPS, SchedStrategy.PSPRIO]:
        # For these strategies, all jobs are "in service"
        # Running = total for these schedulers
        return _cartesian_space_for_phases(s, phases, R)

    elif sched in [SchedStrategy.FCFS, SchedStrategy.HOL, SchedStrategy.LCFS]:
        # FCFS/LCFS: ordered buffer + service
        total_jobs = int(sum(n))
        running_jobs = int(sum(s))

        if total_jobs == 0:
            return np.zeros((1, 1 + int(sum(phases))))

        # Jobs in buffer = n - s per class
        buffer_jobs = n - s

        # Build list of job classes in buffer (with repetition)
        # Using 1-based class IDs to match MATLAB convention
        inbuf = []
        for r in range(R):
            if buffer_jobs[r] > 0:
                inbuf.extend([r + 1] * int(buffer_jobs[r]))  # 1-based class IDs

        # Generate all unique permutations of buffer ordering
        if len(inbuf) > 0:
            mi_buf_list = list(set(permutations(inbuf)))
            mi_buf = np.array(mi_buf_list)
        else:
            mi_buf = np.array([[0]])  # Empty buffer marker

        # si is exactly s (the running constraint)
        si = s

        all_states = []
        for b in range(mi_buf.shape[0]):
            # Generate phase distributions for running jobs
            # kstate = cartesian product of phase distributions per class
            kstate = _cartesian_space_for_phases(si, phases, R)

            # Build full states: [mi_buf[b], kstate]
            for ks in kstate:
                state = np.concatenate([mi_buf[b], ks])
                all_states.append(state)

        if not all_states:
            return np.zeros((1, 1 + int(sum(phases))))

        return np.array(all_states)

    elif sched == SchedStrategy.SIRO:
        # SIRO: unordered buffer + service
        total_jobs = int(sum(n))

        if total_jobs <= S:
            # All jobs in service
            space = _cartesian_space_for_phases(s, phases, R)
            return np.concatenate([np.zeros((space.shape[0], R)), space], axis=1)
        else:
            # Buffer jobs = n - s
            mi_buf = n - s

            all_states = []
            # Generate phase distributions for running jobs
            kstate = _cartesian_space_for_phases(s, phases, R)

            for ks in kstate:
                state = np.concatenate([mi_buf, ks])
                all_states.append(state)

            return np.array(all_states) if all_states else np.array([])

    else:
        # Default: just track service phases
        return _cartesian_space_for_phases(s, phases, R)


def fromMarginalAndStarted(sn, ind: int, n: Union[np.ndarray, list],
                           s: Union[np.ndarray, list], options: dict = None) -> np.ndarray:
    """
    Generate state space with specific marginal and started job counts.

    Creates states where node has n[r] jobs of class r total,
    with s[r] jobs of class r that have started service.
    Started jobs are placed in phase 1 (the initial phase).

    Args:
        sn: NetworkStruct or Network object
        ind: Node index (0-based)
        n: Vector of total job counts per class
        s: Vector of started job counts per class
        options: Optional configuration dictionary

    Returns:
        State space matrix where each row is a valid state.
    """
    # Handle Network object input
    if hasattr(sn, 'getStruct'):
        sn = sn.getStruct()

    if options is None:
        options = {'force': True}

    n = np.atleast_1d(n).astype(int)
    s = np.atleast_1d(s).astype(int)

    R = sn.nclasses if hasattr(sn, 'nclasses') else len(n)

    # Ensure arrays have correct length
    if len(n) < R:
        n = np.concatenate([n, np.zeros(R - len(n), dtype=int)])
    if len(s) < R:
        s = np.concatenate([s, np.zeros(R - len(s), dtype=int)])

    # Get station index
    if hasattr(sn, 'nodeToStation'):
        ist = sn.nodeToStation[ind]
    else:
        ist = ind

    # Get number of servers
    if hasattr(sn, 'nservers'):
        S = int(sn.nservers[ist])
    else:
        S = 1

    # Get phase information
    phases = np.zeros(R, dtype=int)
    if hasattr(sn, 'proc') and sn.proc is not None:
        for r in range(R):
            if hasattr(sn.proc[ist], '__getitem__') and r < len(sn.proc[ist]):
                proc_entry = sn.proc[ist][r]
                if proc_entry is not None:
                    if isinstance(proc_entry, dict):
                        if 'k' in proc_entry:
                            phases[r] = int(proc_entry['k'])
                        elif 'probs' in proc_entry and 'rates' in proc_entry:
                            phases[r] = len(proc_entry['probs'])
                        else:
                            phases[r] = 1
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
    if hasattr(sn, 'classcap'):
        if np.any(n > sn.classcap[ist]):
            return np.array([])

    # Check started constraint: started jobs cannot exceed servers
    if S > 0 and sum(s) > S:
        return np.array([])

    # Get scheduling strategy
    if hasattr(sn, 'sched'):
        sched = sn.sched[ist]
    else:
        sched = SchedStrategy.FCFS

    # Generate state space based on scheduling strategy
    space = _generate_state_space_with_started(sched, S, n, s, phases, R)

    # Sort and unique the state space
    if len(space) > 0:
        space = np.unique(space, axis=0)
        # Reverse sort to put states with jobs in phase 1 earlier
        space = space[::-1]

    return space


def _generate_state_space_with_started(sched: int, S: int, n: np.ndarray, s: np.ndarray,
                                       phases: np.ndarray, R: int) -> np.ndarray:
    """
    Generate local state space with specific started job counts.

    Started jobs are placed in phase 1 only (the initial phase).

    Args:
        sched: SchedStrategy enum value
        S: Number of servers
        n: Total job counts per class
        s: Started job counts per class
        phases: Number of phases per class
        R: Number of classes

    Returns:
        State space matrix for this node
    """
    from itertools import permutations

    if sched in [SchedStrategy.INF, SchedStrategy.PS, SchedStrategy.DPS,
                 SchedStrategy.GPS, SchedStrategy.PSPRIO]:
        # For these strategies, started jobs go in phase 1
        # Build single state with all started jobs in phase 1
        kstate = np.zeros((1, int(sum(phases))))
        col = 0
        for r in range(R):
            # Put all s[r] jobs in phase 1 (column col)
            kstate[0, col] = s[r]
            col += int(phases[r])
        return kstate

    elif sched in [SchedStrategy.FCFS, SchedStrategy.HOL, SchedStrategy.LCFS]:
        # FCFS/LCFS: ordered buffer + service
        total_jobs = int(sum(n))

        if total_jobs == 0:
            return np.zeros((1, 1 + int(sum(phases))))

        # Jobs in buffer = n - s per class
        buffer_jobs = n - s

        # Build list of job classes in buffer (with repetition)
        inbuf = []
        for r in range(R):
            if buffer_jobs[r] > 0:
                inbuf.extend([r + 1] * int(buffer_jobs[r]))  # 1-based class IDs

        # Generate all unique permutations of buffer ordering
        if len(inbuf) > 0:
            mi_buf_list = list(set(permutations(inbuf)))
            mi_buf = np.array(mi_buf_list)
        else:
            mi_buf = np.array([[0]])  # Empty buffer marker

        # For started, all jobs are in phase 1 (no phase enumeration)
        kstate = np.zeros((1, int(sum(phases))))
        col = 0
        for r in range(R):
            kstate[0, col] = s[r]  # All started jobs in phase 1
            col += int(phases[r])

        all_states = []
        for b in range(mi_buf.shape[0]):
            state = np.concatenate([mi_buf[b], kstate[0]])
            all_states.append(state)

        if not all_states:
            return np.zeros((1, 1 + int(sum(phases))))

        return np.array(all_states)

    elif sched == SchedStrategy.SIRO:
        # SIRO: unordered buffer + service
        total_jobs = int(sum(n))

        # Buffer jobs = n - s
        mi_buf = n - s

        # Started jobs in phase 1
        kstate = np.zeros((1, int(sum(phases))))
        col = 0
        for r in range(R):
            kstate[0, col] = s[r]
            col += int(phases[r])

        state = np.concatenate([mi_buf, kstate[0]])
        return np.array([state])

    else:
        # Default: started jobs in phase 1
        kstate = np.zeros((1, int(sum(phases))))
        col = 0
        for r in range(R):
            kstate[0, col] = s[r]
            col += int(phases[r])
        return kstate
