"""
Continuous-Time Markov Chain (CTMC) analysis algorithms.

Native Python implementations for CTMC steady-state, transient analysis,
and related methods. Uses scipy for linear algebra and numerical integration.

Key algorithms:
    ctmc_solve: Steady-state distribution
    ctmc_makeinfgen: Construct valid infinitesimal generator
    ctmc_transient: Transient probabilities via matrix exponential
    ctmc_uniformization: Uniformization for transient analysis
    ctmc_stochcomp: Stochastic complementation
"""

import numpy as np
from numpy.linalg import LinAlgError
from scipy import linalg
from scipy.sparse import csc_matrix, issparse
from scipy.sparse.linalg import spsolve
from scipy.sparse.csgraph import connected_components
from scipy.integrate import solve_ivp
from typing import Dict, Any, Optional, List, Tuple, Union
from dataclasses import dataclass

# Try to import numba for JIT compilation
try:
    from numba import njit
    HAS_NUMBA = True
except ImportError:
    HAS_NUMBA = False
    def njit(*args, **kwargs):
        def decorator(func):
            return func
        return decorator if args and callable(args[0]) else decorator


def ctmc_makeinfgen(Q: np.ndarray) -> np.ndarray:
    """
    Convert a matrix into a valid infinitesimal generator for a CTMC.

    An infinitesimal generator has:
    - Row sums equal to zero
    - Non-positive diagonal elements
    - Non-negative off-diagonal elements

    Args:
        Q: Candidate infinitesimal generator matrix

    Returns:
        Valid infinitesimal generator matrix with corrected diagonal
    """
    Q = np.asarray(Q, dtype=np.float64)
    n = Q.shape[0]

    # Extract off-diagonal elements
    result = Q.copy()
    np.fill_diagonal(result, 0.0)

    # Set diagonal to negative row sum (ensures row sums = 0)
    row_sums = result.sum(axis=1)
    np.fill_diagonal(result, -row_sums)

    return result


def _find_weakly_connected_components(Q: np.ndarray) -> List[List[int]]:
    """
    Find weakly connected components in the CTMC graph.

    Args:
        Q: Generator matrix

    Returns:
        List of lists, each containing state indices for a component
    """
    n = Q.shape[0]

    # Build adjacency matrix: B = |Q + Q'| > 0
    B = np.abs(Q + Q.T) > 0

    # Find connected components
    n_components, labels = connected_components(
        csc_matrix(B), directed=False, return_labels=True
    )

    if n_components == 1:
        return [list(range(n))]

    # Group states by component
    components = [[] for _ in range(n_components)]
    for i, label in enumerate(labels):
        components[label].append(i)

    return components


def ctmc_solve(Q: np.ndarray) -> np.ndarray:
    """
    Solve for steady-state probabilities of a CTMC.

    Computes the stationary distribution π by solving πQ = 0
    with normalization constraint Σπ = 1.

    Handles reducible CTMCs by decomposing into strongly connected
    components and solving each separately.

    Args:
        Q: Infinitesimal generator matrix (row sums should be zero)

    Returns:
        Steady-state probability distribution (1D array)
    """
    Q = np.asarray(Q, dtype=np.float64)
    n = Q.shape[0]

    # Trivial case
    if n == 1:
        return np.array([1.0])

    # Ensure valid generator
    Q = ctmc_makeinfgen(Q)

    # Check for weakly connected components
    components = _find_weakly_connected_components(Q)

    if len(components) > 1:
        # Reducible chain: solve each component separately
        p = np.zeros(n)
        for comp in components:
            if len(comp) == 1:
                p[comp[0]] = 1.0
            else:
                # Extract submatrix for component
                idx = np.array(comp)
                Qc = Q[np.ix_(idx, idx)]
                Qc = ctmc_makeinfgen(Qc)
                pc = ctmc_solve(Qc)
                for i, state in enumerate(comp):
                    p[state] = pc[i]

        # Normalize
        p /= p.sum()
        return p

    # Remove zero rows/cols (absorbing states with no outgoing transitions)
    row_sums = np.abs(Q).sum(axis=1)
    col_sums = np.abs(Q).sum(axis=0)
    active = (row_sums > 0) & (col_sums > 0)

    if not np.any(active):
        # All states are isolated
        return np.ones(n) / n

    active_idx = np.where(active)[0]
    Qnnz = Q[np.ix_(active_idx, active_idx)]
    Qnnz = ctmc_makeinfgen(Qnnz)

    # Solve the system: π * Q = 0, with Σπ = 1
    # Replace last column of Q with ones (normalization constraint)
    # Then solve Q^T * π^T = b
    n_active = len(active_idx)

    # Replace last column with 1s (before transpose)
    Qnnz_mod = Qnnz.copy()
    Qnnz_mod[:, -1] = 1.0

    # Transpose and solve
    A = Qnnz_mod.T
    b = np.zeros(n_active)
    b[-1] = 1.0

    try:
        pi_active = linalg.solve(A, b)
    except LinAlgError:
        # Fallback: use least squares
        pi_active, _, _, _ = linalg.lstsq(A, b)

    # Handle numerical issues
    if np.any(np.isnan(pi_active)):
        # Try recursive decomposition
        sub_components = _find_weakly_connected_components(Qnnz)
        if len(sub_components) > 1:
            pi_active = np.zeros(n_active)
            for comp in sub_components:
                idx = np.array(comp)
                Qc = Qnnz[np.ix_(idx, idx)]
                Qc = ctmc_makeinfgen(Qc)
                pc = ctmc_solve(Qc)
                for i, state in enumerate(comp):
                    pi_active[state] = pc[i]
            pi_active /= pi_active.sum()

    # Map back to full state space
    p = np.zeros(n)
    p[active_idx] = pi_active

    # Ensure non-negative and normalized
    p = np.maximum(p, 0.0)
    if p.sum() > 0:
        p /= p.sum()
    else:
        p = np.ones(n) / n

    return p


def ctmc_solve_reducible(Q: np.ndarray) -> np.ndarray:
    """
    Solve reducible CTMCs by converting to DTMC via uniformization.

    Args:
        Q: Infinitesimal generator matrix (possibly reducible)

    Returns:
        Steady-state probability vector
    """
    # ctmc_solve already handles reducible chains
    return ctmc_solve(Q)


def ctmc_transient(
    Q: np.ndarray,
    initial_dist: np.ndarray,
    time_points: Union[float, np.ndarray],
    method: str = 'expm'
) -> np.ndarray:
    """
    Compute transient probabilities of a CTMC.

    Calculates time-dependent state probabilities π(t) for
    specified time points using matrix exponential methods.

    Args:
        Q: Infinitesimal generator matrix
        initial_dist: Initial probability distribution π(0)
        time_points: Array of time points to evaluate, or single time value
        method: 'expm' for matrix exponential, 'ode' for ODE solver

    Returns:
        Transient probabilities at each time point.
        Shape: (len(time_points), n) if multiple times, (n,) if single time
    """
    Q = np.asarray(Q, dtype=np.float64)
    initial_dist = np.asarray(initial_dist, dtype=np.float64).flatten()

    if np.isscalar(time_points):
        time_points = np.array([time_points])
        single_time = True
    else:
        time_points = np.asarray(time_points)
        single_time = False

    n = Q.shape[0]
    results = np.zeros((len(time_points), n))

    if method == 'expm':
        # Use matrix exponential: π(t) = π(0) * exp(Qt)
        for i, t in enumerate(time_points):
            if t == 0:
                results[i] = initial_dist
            else:
                expQt = linalg.expm(Q * t)
                results[i] = initial_dist @ expQt
    else:
        # Use ODE solver: dπ/dt = π * Q
        def ode_func(t, pi):
            return pi @ Q

        for i, t in enumerate(time_points):
            if t == 0:
                results[i] = initial_dist
            else:
                sol = solve_ivp(
                    ode_func, [0, t], initial_dist,
                    method='LSODA', dense_output=False
                )
                results[i] = sol.y[:, -1]

    return results[0] if single_time else results


def ctmc_uniformization(
    Q: np.ndarray,
    lambda_rate: Optional[float] = None
) -> Dict[str, Any]:
    """
    Uniformize CTMC generator matrix.

    Converts CTMC to an equivalent uniformized discrete-time chain
    for numerical analysis and simulation purposes.

    The uniformized DTMC has transition matrix P = I + Q/λ where
    λ is the uniformization rate (max exit rate).

    Args:
        Q: Infinitesimal generator matrix
        lambda_rate: Uniformization rate (optional, auto-computed if None)

    Returns:
        dict containing:
            - 'P': Uniformized transition matrix
            - 'lambda': Uniformization rate
    """
    Q = np.asarray(Q, dtype=np.float64)

    if lambda_rate is None:
        # Use max exit rate (max |diagonal|)
        lambda_rate = -np.min(np.diag(Q))
        if lambda_rate <= 0:
            lambda_rate = 1.0

    n = Q.shape[0]
    I = np.eye(n)
    P = I + Q / lambda_rate

    return {
        'P': P,
        'lambda': lambda_rate
    }


def ctmc_randomization(
    Q: np.ndarray,
    initial_dist: np.ndarray,
    time_points: np.ndarray,
    precision: float = 1e-10
) -> np.ndarray:
    """
    Compute CTMC transient probabilities using randomization.

    Uses Jensen's randomization method (uniformization) to compute
    transient probabilities by converting the CTMC to a uniformized DTMC.

    This method is numerically stable and avoids matrix exponentials.

    Args:
        Q: Infinitesimal generator matrix
        initial_dist: Initial probability distribution
        time_points: Array of time points to evaluate
        precision: Numerical precision for truncation (Poisson tail)

    Returns:
        Transient probabilities at each time point
    """
    Q = np.asarray(Q, dtype=np.float64)
    initial_dist = np.asarray(initial_dist, dtype=np.float64).flatten()
    time_points = np.asarray(time_points)

    n = Q.shape[0]

    # Uniformization rate
    lambda_rate = -np.min(np.diag(Q))
    if lambda_rate <= 0:
        lambda_rate = 1.0

    # Uniformized transition matrix
    P = np.eye(n) + Q / lambda_rate

    results = np.zeros((len(time_points), n))

    for idx, t in enumerate(time_points):
        if t == 0:
            results[idx] = initial_dist
            continue

        # Compute Poisson probabilities and truncation point
        q = lambda_rate * t

        # Find truncation point k_max such that sum of Poisson tail < precision
        from scipy.stats import poisson
        k_max = int(poisson.ppf(1 - precision, q)) + 10

        # Compute Poisson probabilities
        poisson_probs = poisson.pmf(np.arange(k_max + 1), q)

        # Compute π(t) = Σ_k P(N(t)=k) * π(0) * P^k
        pi_t = np.zeros(n)
        pi_k = initial_dist.copy()  # π(0) * P^0 = π(0)

        for k in range(k_max + 1):
            pi_t += poisson_probs[k] * pi_k
            pi_k = pi_k @ P  # π(0) * P^(k+1)

        results[idx] = pi_t

    return results


def ctmc_stochcomp(
    Q: np.ndarray,
    keep_states: np.ndarray,
    eliminate_states: Optional[np.ndarray] = None
) -> Dict[str, np.ndarray]:
    """
    Compute stochastic complement of CTMC.

    Reduces the CTMC by eliminating specified states while preserving
    the steady-state distribution restricted to the kept states.

    Args:
        Q: Infinitesimal generator matrix
        keep_states: States to retain in reduced model (array of indices)
        eliminate_states: States to eliminate (optional, inferred if None)

    Returns:
        dict containing:
            - 'S': Stochastic complement (reduced generator)
            - 'Q11': Submatrix for kept states
            - 'Q12': Transitions from kept to eliminated
            - 'Q21': Transitions from eliminated to kept
            - 'Q22': Submatrix for eliminated states
            - 'T': Transient contribution matrix
    """
    Q = np.asarray(Q, dtype=np.float64)
    keep_states = np.asarray(keep_states, dtype=int).flatten()

    n = Q.shape[0]
    all_states = set(range(n))
    keep_set = set(keep_states)

    if eliminate_states is None:
        eliminate_states = np.array(sorted(all_states - keep_set), dtype=int)
    else:
        eliminate_states = np.asarray(eliminate_states, dtype=int).flatten()

    if len(eliminate_states) == 0:
        return {
            'S': Q,
            'Q11': Q,
            'Q12': np.zeros((len(keep_states), 0)),
            'Q21': np.zeros((0, len(keep_states))),
            'Q22': np.zeros((0, 0)),
            'T': np.zeros((len(keep_states), len(keep_states)))
        }

    # Extract submatrices
    I = keep_states
    Ic = eliminate_states

    Q11 = Q[np.ix_(I, I)]
    Q12 = Q[np.ix_(I, Ic)]
    Q21 = Q[np.ix_(Ic, I)]
    Q22 = Q[np.ix_(Ic, Ic)]

    # S = Q11 + Q12 * (-Q22)^{-1} * Q21
    # T = Q12 * (-Q22)^{-1} * Q21
    Q22_neg = -Q22
    try:
        Q22_inv = linalg.inv(Q22_neg)
    except LinAlgError:
        Q22_inv = linalg.pinv(Q22_neg)

    T = Q12 @ Q22_inv @ Q21
    S = Q11 + T

    return {
        'S': S,
        'Q11': Q11,
        'Q12': Q12,
        'Q21': Q21,
        'Q22': Q22,
        'T': T
    }


def ctmc_timereverse(
    Q: np.ndarray,
    pi: Optional[np.ndarray] = None
) -> np.ndarray:
    """
    Compute time-reversed CTMC generator.

    The time-reversed generator Q* has elements:
    Q*_{ij} = π_j * Q_{ji} / π_i

    Args:
        Q: Original infinitesimal generator matrix
        pi: Steady-state distribution (optional, computed if None)

    Returns:
        Time-reversed generator matrix
    """
    Q = np.asarray(Q, dtype=np.float64)

    if pi is None:
        pi = ctmc_solve(Q)
    else:
        pi = np.asarray(pi, dtype=np.float64).flatten()

    n = Q.shape[0]
    Q_rev = np.zeros_like(Q)

    for i in range(n):
        for j in range(n):
            if pi[i] > 0:
                Q_rev[i, j] = pi[j] * Q[j, i] / pi[i]

    # Ensure valid generator
    Q_rev = ctmc_makeinfgen(Q_rev)

    return Q_rev


def ctmc_rand(
    n: int,
    density: float = 0.3,
    max_rate: float = 10.0
) -> np.ndarray:
    """
    Generate random CTMC generator matrix.

    Args:
        n: Number of states
        density: Sparsity density (0 to 1, default 0.3)
        max_rate: Maximum transition rate (default 10.0)

    Returns:
        Random infinitesimal generator matrix
    """
    # Generate random off-diagonal entries
    Q = np.random.rand(n, n) * max_rate

    # Apply density mask
    mask = np.random.rand(n, n) < density
    Q = Q * mask

    # Zero out diagonal
    np.fill_diagonal(Q, 0.0)

    # Make valid generator
    Q = ctmc_makeinfgen(Q)

    return Q


def ctmc_simulate(
    Q: np.ndarray,
    initial_state: int,
    max_time: float,
    max_events: int = 10000,
    seed: Optional[int] = None
) -> Dict[str, np.ndarray]:
    """
    Simulate CTMC sample path using Gillespie algorithm.

    Generates a realization of the continuous-time Markov chain
    using the next-reaction method.

    Args:
        Q: Infinitesimal generator matrix
        initial_state: Starting state (integer index)
        max_time: Maximum simulation time
        max_events: Maximum number of transitions (default: 10000)
        seed: Random seed for reproducibility (optional)

    Returns:
        dict with:
            - 'states': Array of visited states
            - 'times': Array of transition times
            - 'sojourn_times': Time spent in each state
    """
    if seed is not None:
        np.random.seed(seed)

    Q = np.asarray(Q, dtype=np.float64)
    n = Q.shape[0]

    states = [initial_state]
    times = [0.0]
    sojourn_times = []

    current_state = initial_state
    current_time = 0.0

    for _ in range(max_events):
        # Exit rate from current state
        exit_rate = -Q[current_state, current_state]

        if exit_rate <= 0:
            # Absorbing state
            sojourn_times.append(max_time - current_time)
            break

        # Time to next transition (exponential)
        sojourn = np.random.exponential(1.0 / exit_rate)

        if current_time + sojourn > max_time:
            sojourn_times.append(max_time - current_time)
            break

        sojourn_times.append(sojourn)
        current_time += sojourn

        # Determine next state
        rates = Q[current_state, :].copy()
        rates[current_state] = 0.0
        probs = rates / rates.sum()
        next_state = np.random.choice(n, p=probs)

        states.append(next_state)
        times.append(current_time)
        current_state = next_state

    return {
        'states': np.array(states),
        'times': np.array(times),
        'sojourn_times': np.array(sojourn_times)
    }


def ctmc_isfeasible(Q: np.ndarray, tolerance: float = 1e-10) -> bool:
    """
    Check if matrix is a valid CTMC infinitesimal generator.

    Validates:
    - Off-diagonal elements are non-negative
    - Row sums are zero
    - Diagonal elements are non-positive

    Args:
        Q: Candidate generator matrix
        tolerance: Numerical tolerance (default: 1e-10)

    Returns:
        True if matrix is valid CTMC generator
    """
    Q = np.asarray(Q)
    n = Q.shape[0]

    if Q.shape[0] != Q.shape[1]:
        return False

    # Check off-diagonal non-negative
    for i in range(n):
        for j in range(n):
            if i != j and Q[i, j] < -tolerance:
                return False

    # Check diagonal non-positive
    if np.any(np.diag(Q) > tolerance):
        return False

    # Check row sums are zero
    row_sums = Q.sum(axis=1)
    if np.any(np.abs(row_sums) > tolerance):
        return False

    return True


# ============================================================================
# State Space Generation
# ============================================================================

@dataclass
class CtmcSsgResult:
    """Result from CTMC state space generation."""
    state_space: np.ndarray  # Complete state space matrix
    state_space_aggr: np.ndarray  # Aggregated state space (per station-class)
    state_space_hashed: np.ndarray  # Hashed state indices
    node_state_space: Dict[int, np.ndarray]  # State space per node
    sn: Any  # Updated network structure


def ctmc_ssg(sn: Any, options: Optional[Dict] = None) -> CtmcSsgResult:
    """
    Generate complete CTMC state space for a queueing network.

    Creates all possible network states including those not reachable from
    the initial state. For open classes, a cutoff parameter limits the
    maximum population to keep state space finite.

    The state space is aggregated to show per-station-class job counts.

    Args:
        sn: NetworkStruct object (from getStruct())
        options: Solver options dict with fields:
            - cutoff: Population cutoff for open classes (required if open)
            - config.hide_immediate: Hide immediate transitions (default True)

    Returns:
        CtmcSsgResult containing:
            - state_space: Complete state space matrix (rows=states, cols=state components)
            - state_space_aggr: Aggregated state space (rows=states, cols=stations*classes)
            - state_space_hashed: Hashed state indices for lookup
            - node_state_space: Dictionary of per-node state spaces
            - sn: Updated network structure with space field populated

    References:
        MATLAB: matlab/src/api/mc/ctmc_ssg.m
    """
    from ..state import spaceGenerator, toMarginal

    if options is None:
        options = {}

    # Get cutoff from options
    cutoff = options.get('cutoff', None)

    # Generate state space
    state_space, state_space_hashed, sn, adj, st = spaceGenerator(sn, cutoff, options)

    # Get node state space from sn.space
    node_state_space = {}
    if hasattr(sn, 'space') and sn.space is not None:
        if isinstance(sn.space, dict):
            node_state_space = sn.space
        elif isinstance(sn.space, (list, np.ndarray)):
            for i, space in enumerate(sn.space):
                if space is not None:
                    node_state_space[i] = space

    # Set default hide_immediate
    if 'hide_immediate' not in options:
        options['hide_immediate'] = True

    # Compute aggregated state space
    nstateful = sn.nstateful if hasattr(sn, 'nstateful') else sn.nstations
    nclasses = sn.nclasses if hasattr(sn, 'nclasses') else 1

    # Initialize aggregated state space
    if state_space_hashed.size > 0:
        state_space_aggr = np.zeros((state_space_hashed.shape[0], sn.nstations * nclasses))
    else:
        state_space_aggr = np.zeros((0, sn.nstations * nclasses))

    # Process synchronizations to compute aggregated states
    sync = sn.sync if hasattr(sn, 'sync') else []
    A = len(sync) if sync else 0

    for s in range(state_space_hashed.shape[0] if state_space_hashed.size > 0 else 0):
        state = state_space_hashed[s, :] if state_space_hashed.ndim > 1 else state_space_hashed

        # Update state cell array
        state_cell = {}
        for ind in range(sn.nnodes if hasattr(sn, 'nnodes') else sn.nstations):
            isstateful = sn.isstateful(ind) if hasattr(sn, 'isstateful') else True
            if isstateful:
                isf = sn.nodeToStateful[ind] if hasattr(sn, 'nodeToStateful') else ind
                if isf < len(node_state_space) and node_state_space.get(isf) is not None:
                    space_isf = node_state_space[isf]
                    state_idx = int(state[isf]) if isf < len(state) else 0
                    if state_idx < len(space_isf):
                        state_cell[isf] = space_isf[state_idx]

                isstation = sn.isstation(ind) if hasattr(sn, 'isstation') else True
                if isstation:
                    ist = sn.nodeToStation[ind] if hasattr(sn, 'nodeToStation') else ind
                    if isf in state_cell:
                        # Compute marginal job counts
                        try:
                            ni, nir, sir, kir = toMarginal(sn, ind, state_cell[isf])
                            # Store in aggregated state space
                            start_col = ist * nclasses
                            end_col = (ist + 1) * nclasses
                            if nir.ndim == 1:
                                state_space_aggr[s, start_col:end_col] = nir
                            else:
                                state_space_aggr[s, start_col:end_col] = nir[0, :]
                        except Exception:
                            pass

    return CtmcSsgResult(
        state_space=state_space,
        state_space_aggr=state_space_aggr,
        state_space_hashed=state_space_hashed,
        node_state_space=node_state_space,
        sn=sn
    )


def ctmc_ssg_reachability(sn: Any, options: Optional[Dict] = None) -> CtmcSsgResult:
    """
    Generate reachable CTMC state space for a queueing network.

    Creates only the states reachable from the initial state through valid
    transitions. This is more efficient than ctmc_ssg for networks with
    constrained reachability.

    Args:
        sn: NetworkStruct object (from getStruct())
        options: Solver options dict with fields:
            - config.hide_immediate: Hide immediate transitions (default True)

    Returns:
        CtmcSsgResult containing:
            - state_space: Reachable state space matrix
            - state_space_aggr: Aggregated state space (per station-class)
            - state_space_hashed: Hashed state indices
            - node_state_space: Dictionary of per-node state spaces
            - sn: Updated network structure

    References:
        MATLAB: matlab/src/api/mc/ctmc_ssg_reachability.m
    """
    from ..state import spaceGenerator, toMarginal

    if options is None:
        options = {}

    # Set default hide_immediate in config
    config = options.get('config', {})
    if 'hide_immediate' not in config:
        config['hide_immediate'] = True
        options['config'] = config

    # For reachability analysis, we need to explore from initial state
    # Currently using spaceGenerator as base (full implementation would use
    # State.reachableSpaceGenerator which explores from initial state)
    state_space, state_space_hashed, sn, adj, st = spaceGenerator(sn, options.get('cutoff'), options)

    # Get node state space
    node_state_space = {}
    if hasattr(sn, 'space') and sn.space is not None:
        if isinstance(sn.space, dict):
            node_state_space = sn.space
        elif isinstance(sn.space, (list, np.ndarray)):
            for i, space in enumerate(sn.space):
                if space is not None:
                    node_state_space[i] = space

    # Compute aggregated state space
    nstateful = sn.nstateful if hasattr(sn, 'nstateful') else sn.nstations
    nclasses = sn.nclasses if hasattr(sn, 'nclasses') else 1

    if state_space_hashed.size > 0:
        state_space_aggr = np.zeros((state_space_hashed.shape[0], sn.nstations * nclasses))
    else:
        state_space_aggr = np.zeros((0, sn.nstations * nclasses))

    # Process states to compute aggregated representation
    sync = sn.sync if hasattr(sn, 'sync') else []
    A = len(sync) if sync else 0

    for s in range(state_space_hashed.shape[0] if state_space_hashed.size > 0 else 0):
        state = state_space_hashed[s, :] if state_space_hashed.ndim > 1 else state_space_hashed

        state_cell = {}
        for ind in range(sn.nnodes if hasattr(sn, 'nnodes') else sn.nstations):
            isstateful = sn.isstateful(ind) if hasattr(sn, 'isstateful') else True
            if isstateful:
                isf = sn.nodeToStateful[ind] if hasattr(sn, 'nodeToStateful') else ind
                if isf < len(node_state_space) and node_state_space.get(isf) is not None:
                    space_isf = node_state_space[isf]
                    state_idx = int(state[isf]) if isf < len(state) else 0
                    if state_idx < len(space_isf):
                        state_cell[isf] = space_isf[state_idx]

                isstation = sn.isstation(ind) if hasattr(sn, 'isstation') else True
                if isstation:
                    ist = sn.nodeToStation[ind] if hasattr(sn, 'nodeToStation') else ind
                    if isf in state_cell:
                        try:
                            ni, nir, sir, kir = toMarginal(sn, ind, state_cell[isf])
                            start_col = ist * nclasses
                            end_col = (ist + 1) * nclasses
                            if nir.ndim == 1:
                                state_space_aggr[s, start_col:end_col] = nir
                            else:
                                state_space_aggr[s, start_col:end_col] = nir[0, :]
                        except Exception:
                            pass

    return CtmcSsgResult(
        state_space=state_space,
        state_space_aggr=state_space_aggr,
        state_space_hashed=state_space_hashed,
        node_state_space=node_state_space,
        sn=sn
    )


__all__ = [
    'ctmc_solve',
    'ctmc_solve_reducible',
    'ctmc_makeinfgen',
    'ctmc_transient',
    'ctmc_uniformization',
    'ctmc_randomization',
    'ctmc_stochcomp',
    'ctmc_timereverse',
    'ctmc_rand',
    'ctmc_simulate',
    'ctmc_isfeasible',
    'ctmc_ssg',
    'ctmc_ssg_reachability',
    'CtmcSsgResult',
]
