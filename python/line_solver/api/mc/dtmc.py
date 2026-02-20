"""
Discrete-Time Markov Chain (DTMC) analysis algorithms.

Native Python implementations for DTMC steady-state and related methods.
Leverages the CTMC solver by converting (P - I) to a generator-like matrix.

Key algorithms:
    dtmc_solve: Steady-state distribution
    dtmc_makestochastic: Convert matrix to stochastic
    dtmc_isfeasible: Validate transition matrix
    dtmc_simulate: Sample path simulation
"""

import numpy as np
from numpy.linalg import LinAlgError
from scipy import linalg
from typing import Dict, Any, Optional, List

from .ctmc import ctmc_solve, _find_weakly_connected_components


def dtmc_solve(P: np.ndarray) -> np.ndarray:
    """
    Solve for steady-state probabilities of a DTMC.

    Computes the stationary distribution π by solving π(P - I) = 0
    with normalization constraint Σπ = 1.

    This leverages the CTMC solver by treating (P - I) as an
    infinitesimal generator.

    Args:
        P: Transition probability matrix (row stochastic)

    Returns:
        Steady-state probability distribution (1D array)
    """
    P = np.asarray(P, dtype=np.float64)

    # Convert to "generator" form: Q = P - I
    # This has the property that πQ = π(P-I) = πP - π = 0 when πP = π
    n = P.shape[0]
    Q = P.copy()
    for i in range(n):
        Q[i, i] -= 1.0

    return ctmc_solve(Q)


def dtmc_solve_reducible(P: np.ndarray, pin: np.ndarray = None) -> np.ndarray:
    """
    Solve reducible DTMCs with transient states.

    Handles DTMCs with multiple recurrent classes and transient states by:
    1. Decomposing into strongly connected components (SCCs)
    2. Identifying recurrent vs transient SCCs
    3. Computing limiting distribution considering absorption from transient states

    For a reducible DTMC with a single transient SCC, this computes the limiting
    distribution when starting from the transient states (e.g., class switching
    networks where jobs start in a transient class).

    Args:
        P: Transition probability matrix (possibly reducible)
        pin: Initial probability vector (optional)

    Returns:
        Steady-state probability vector
    """
    from scipy.sparse.csgraph import connected_components
    from scipy.sparse import csc_matrix

    P = np.asarray(P, dtype=np.float64)
    n = P.shape[0]

    if n == 1:
        return np.array([1.0])

    # Find strongly connected components
    # Use 'strong' connection to properly identify recurrent vs transient
    n_components, scc_labels = connected_components(
        csc_matrix(P > 1e-10), directed=True, connection='strong', return_labels=True
    )

    if n_components == 1:
        # Irreducible chain - use standard solve
        return dtmc_solve(P)

    # Identify recurrent vs transient SCCs
    # A SCC is recurrent if all outgoing transitions stay within the SCC
    num_scc = n_components
    scc_idx = [np.where(scc_labels == i)[0] for i in range(num_scc)]

    is_rec = np.zeros(num_scc, dtype=bool)
    for i in range(num_scc):
        states_i = scc_idx[i]
        # Check if all outgoing transitions from SCC i stay within SCC i
        outgoing_prob = 0.0
        for s in states_i:
            for j in range(num_scc):
                if j != i:
                    states_j = scc_idx[j]
                    outgoing_prob += np.sum(P[s, states_j])
        # If no outgoing transitions to other SCCs, this is recurrent
        is_rec[i] = outgoing_prob < 1e-10

    # Build lumped transition matrix between SCCs
    Pl = np.zeros((num_scc, num_scc))
    for i in range(num_scc):
        states_i = scc_idx[i]
        for j in range(num_scc):
            if i != j:
                states_j = scc_idx[j]
                Pl[i, j] = np.sum(P[np.ix_(states_i, states_j)])

    Pl = dtmc_makestochastic(Pl)

    # Recurrent SCCs have self-loops with probability 1
    for i in range(num_scc):
        if is_rec[i]:
            Pl[i, :] = 0.0
            Pl[i, i] = 1.0

    # Compute initial distribution over SCCs
    # This matches MATLAB's dtmc_solve_reducible lines 72-89
    if pin is None:
        # Start with uniform over all SCCs
        pinl = np.ones(num_scc)

        # Find states with near-zero column sums (absorbing/unreachable states)
        # These states have no incoming transitions from other states
        col_sums = np.sum(P, axis=0)
        z_cols = np.where(col_sums < 1e-12)[0]

        # Zero out SCCs that contain only absorbing states
        for j in z_cols:
            pinl[scc_labels[j]] = 0

        # Normalize
        if np.sum(pinl) > 0:
            pinl = pinl / np.sum(pinl)
        else:
            # Fallback to uniform
            pinl = np.ones(num_scc) / num_scc
    else:
        # Lump initial vector by SCC
        pinl = np.zeros(num_scc)
        for i in range(num_scc):
            pinl[i] = np.sum(pin[scc_idx[i]])

    # Compute limiting distribution of lumped chain via power iteration
    PI = _compute_limiting_matrix_power(Pl)

    # Compute stationary distribution for each starting SCC
    pis = np.zeros((num_scc, n))
    for i in range(num_scc):
        if pinl[i] > 1e-10:
            # Compute limiting distribution starting from SCC i
            pi0_l = np.zeros(num_scc)
            pi0_l[i] = 1.0
            pil_i = pi0_l @ PI

            # For each SCC, solve internal stationary distribution
            for j in range(num_scc):
                if pil_i[j] > 1e-10:
                    states_j = scc_idx[j]
                    if len(states_j) == 1:
                        pis[i, states_j[0]] = pil_i[j]
                    else:
                        # Solve internal DTMC
                        Pj = P[np.ix_(states_j, states_j)]
                        pij = dtmc_solve(Pj)
                        pis[i, states_j] = pil_i[j] * pij

    # Combine distributions weighted by initial SCC probabilities
    pi = np.zeros(n)
    for i in range(num_scc):
        if pinl[i] > 1e-10:
            pi += pis[i, :] * pinl[i]

    # Special case: single transient SCC with no explicit initial distribution
    transient_sccs = np.where(~is_rec)[0]
    if len(transient_sccs) == 1 and pin is None:
        pi = pis[transient_sccs[0], :]

    return pi


def _compute_limiting_matrix_power(P: np.ndarray, max_iter: int = 1000, tol: float = 1e-10) -> np.ndarray:
    """
    Compute limiting matrix using power iteration.

    For DTMCs, P^n converges to a matrix where each row gives the limiting
    distribution starting from that state.

    Args:
        P: Stochastic transition matrix
        max_iter: Maximum iterations (default: 1000)
        tol: Convergence tolerance (default: 1e-10)

    Returns:
        Limiting matrix (each row is limiting distribution from that start state)
    """
    Pk = P.copy()
    for _ in range(max_iter):
        Pk1 = Pk @ P
        max_diff = np.max(np.abs(Pk1 - Pk))
        Pk = Pk1
        if max_diff < tol:
            break
    return Pk


def dtmc_makestochastic(A: np.ndarray) -> np.ndarray:
    """
    Convert matrix to row-stochastic transition matrix.

    Normalizes each row to sum to 1. Rows with zero sum are
    replaced with uniform distribution.

    Args:
        A: Input matrix to normalize

    Returns:
        Row-stochastic matrix
    """
    P = np.asarray(A, dtype=np.float64).copy()
    n = P.shape[0]

    for i in range(n):
        row_sum = np.sum(P[i, :])
        if row_sum > 0:
            P[i, :] = P[i, :] / row_sum
            P[i, i] = min(max(0.0, 1.0 - (np.sum(P[i, :]) - P[i, i])), 1.0)
        else:
            P[i, :] = 0.0
            P[i, i] = 1.0

    return P


def dtmc_isfeasible(P: np.ndarray, tolerance: float = 1e-10) -> bool:
    """
    Check if matrix is a valid DTMC transition matrix.

    Validates:
    - All elements are non-negative
    - All row sums equal 1

    Args:
        P: Candidate transition matrix
        tolerance: Numerical tolerance (default: 1e-10)

    Returns:
        True if matrix is valid DTMC transition matrix
    """
    P = np.asarray(P)

    if P.shape[0] != P.shape[1]:
        return False

    # Check non-negative
    if np.any(P < -tolerance):
        return False

    # Check row sums = 1
    row_sums = P.sum(axis=1)
    if np.any(np.abs(row_sums - 1.0) > tolerance):
        return False

    return True


def dtmc_simulate(
    P: np.ndarray,
    initial_state: int,
    num_steps: int,
    seed: Optional[int] = None
) -> np.ndarray:
    """
    Simulate DTMC sample path.

    Generates a realization of the discrete-time Markov chain
    for a specified number of steps.

    Args:
        P: Transition probability matrix
        initial_state: Starting state index
        num_steps: Number of simulation steps

    Returns:
        Array of visited states (length num_steps + 1)
    """
    if seed is not None:
        np.random.seed(seed)

    P = np.asarray(P, dtype=np.float64)
    n = P.shape[0]

    states = np.zeros(num_steps + 1, dtype=int)
    states[0] = initial_state
    current_state = initial_state

    for step in range(num_steps):
        probs = P[current_state, :]
        next_state = np.random.choice(n, p=probs)
        states[step + 1] = next_state
        current_state = next_state

    return states


def dtmc_rand(
    n: int,
    density: float = 0.5
) -> np.ndarray:
    """
    Generate random DTMC transition matrix.

    Args:
        n: Number of states
        density: Sparsity density (0 to 1, default 0.5)

    Returns:
        Random transition probability matrix
    """
    # Generate random entries
    P = np.random.rand(n, n)

    # Apply density mask
    mask = np.random.rand(n, n) < density
    P = P * mask

    # Ensure at least one non-zero per row
    for i in range(n):
        if P[i, :].sum() == 0:
            j = np.random.randint(n)
            P[i, j] = np.random.rand()

    # Normalize to stochastic
    return dtmc_makestochastic(P)


def dtmc_timereverse(
    P: np.ndarray,
    pi: Optional[np.ndarray] = None
) -> np.ndarray:
    """
    Compute time-reversed DTMC transition matrix.

    The time-reversed chain has transition probabilities:
    P*_{ij} = π_j * P_{ji} / π_i

    Args:
        P: Original transition matrix
        pi: Steady-state distribution (optional, computed if None)

    Returns:
        Time-reversed transition matrix
    """
    P = np.asarray(P, dtype=np.float64)

    if pi is None:
        pi = dtmc_solve(P)
    else:
        pi = np.asarray(pi, dtype=np.float64).flatten()

    n = P.shape[0]
    P_rev = np.zeros_like(P)

    for i in range(n):
        for j in range(n):
            if pi[i] > 0:
                P_rev[i, j] = pi[j] * P[j, i] / pi[i]

    # Ensure stochastic
    P_rev = dtmc_makestochastic(P_rev)

    return P_rev


def dtmc_stochcomp(
    P: np.ndarray,
    keep_states: np.ndarray,
    eliminate_states: Optional[np.ndarray] = None
) -> np.ndarray:
    """
    Compute stochastic complement of DTMC.

    Reduces the DTMC by eliminating specified states while preserving
    the steady-state distribution restricted to the kept states.

    Args:
        P: Transition probability matrix
        keep_states: States to retain in reduced model
        eliminate_states: States to eliminate (optional, inferred if None)

    Returns:
        Reduced transition matrix (stochastic complement)
    """
    P = np.asarray(P, dtype=np.float64)
    keep_states = np.asarray(keep_states, dtype=int).flatten()

    n = P.shape[0]
    all_states = set(range(n))
    keep_set = set(keep_states)

    if eliminate_states is None:
        eliminate_states = np.array(sorted(all_states - keep_set), dtype=int)
    else:
        eliminate_states = np.asarray(eliminate_states, dtype=int).flatten()

    if len(eliminate_states) == 0:
        return P

    # Extract submatrices
    I = keep_states
    Ic = eliminate_states

    P11 = P[np.ix_(I, I)]
    P12 = P[np.ix_(I, Ic)]
    P21 = P[np.ix_(Ic, I)]
    P22 = P[np.ix_(Ic, Ic)]

    # Stochastic complement: S = P11 + P12 * (I - P22)^{-1} * P21
    n_elim = len(eliminate_states)
    I_mat = np.eye(n_elim)

    try:
        inv_term = linalg.inv(I_mat - P22)
    except LinAlgError:
        inv_term = linalg.pinv(I_mat - P22)

    S = P11 + P12 @ inv_term @ P21

    return S


def dtmc_transient(
    P: np.ndarray,
    initial_dist: np.ndarray,
    steps: int
) -> np.ndarray:
    """
    Compute transient probabilities of a DTMC.

    Calculates π(n) = π(0) * P^n for each step from 0 to steps.

    Args:
        P: Transition probability matrix
        initial_dist: Initial probability distribution π(0)
        steps: Number of time steps

    Returns:
        Array of shape (steps+1, n) with transient probabilities
    """
    P = np.asarray(P, dtype=np.float64)
    initial_dist = np.asarray(initial_dist, dtype=np.float64).flatten()

    n = P.shape[0]
    results = np.zeros((steps + 1, n))
    results[0] = initial_dist

    pi = initial_dist.copy()
    for k in range(1, steps + 1):
        pi = pi @ P
        results[k] = pi

    return results


def dtmc_hitting_time(
    P: np.ndarray,
    target_states: np.ndarray
) -> np.ndarray:
    """
    Compute mean hitting times to target states.

    Calculates the expected number of steps to reach any target state
    from each starting state.

    Args:
        P: Transition probability matrix
        target_states: Array of target state indices

    Returns:
        Array of mean hitting times from each state
    """
    P = np.asarray(P, dtype=np.float64)
    target_states = np.asarray(target_states, dtype=int).flatten()

    n = P.shape[0]
    target_set = set(target_states)
    non_target = [i for i in range(n) if i not in target_set]

    if len(non_target) == 0:
        return np.zeros(n)

    # For non-target states: h = 1 + P_NT * h_NT
    # (I - P_NT) * h_NT = 1
    P_NT = P[np.ix_(non_target, non_target)]
    A = np.eye(len(non_target)) - P_NT
    b = np.ones(len(non_target))

    try:
        h_NT = linalg.solve(A, b)
    except LinAlgError:
        h_NT, _, _, _ = linalg.lstsq(A, b)

    h = np.zeros(n)
    h[non_target] = h_NT

    return h


__all__ = [
    'dtmc_solve',
    'dtmc_solve_reducible',
    'dtmc_makestochastic',
    'dtmc_isfeasible',
    'dtmc_simulate',
    'dtmc_rand',
    'dtmc_timereverse',
    'dtmc_stochcomp',
    'dtmc_transient',
    'dtmc_hitting_time',
]
