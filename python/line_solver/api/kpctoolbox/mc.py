"""
Markov Chain analysis functions for KPC-Toolbox.

Native Python implementations of Continuous-Time (CTMC) and Discrete-Time (DTMC)
Markov Chain analysis functions.
"""

import numpy as np
from numpy.linalg import solve, inv
from scipy.linalg import lu_factor, lu_solve
from typing import Tuple, NamedTuple, Union
from collections import deque


class ConnectedComponents(NamedTuple):
    """Result of connected component analysis."""
    num_components: int
    component_assignment: np.ndarray


class CTMCSolveResult(NamedTuple):
    """Result of CTMC solving."""
    equilibrium_distribution: np.ndarray
    generator: np.ndarray
    num_components: int
    component_assignment: np.ndarray


# ============================================================================
# CTMC Functions
# ============================================================================

def ctmc_makeinfgen(Q: np.ndarray) -> np.ndarray:
    """
    Normalize a matrix to be a valid infinitesimal generator.

    Sets diagonal elements such that row sums are zero.

    Args:
        Q: Input matrix

    Returns:
        Valid infinitesimal generator matrix where rows sum to zero
    """
    Q = np.asarray(Q, dtype=np.float64)
    n = Q.shape[0]
    result = Q.copy()

    for i in range(n):
        # Zero out diagonal first
        result[i, i] = 0.0
        # Set diagonal to make row sum zero
        result[i, i] = -np.sum(result[i, :])

    return result


def ctmc_rand(n: int, seed: int = None) -> np.ndarray:
    """
    Generate a random infinitesimal generator matrix.

    Args:
        n: Size of the matrix (n x n)
        seed: Random seed (optional)

    Returns:
        Random infinitesimal generator
    """
    if seed is not None:
        np.random.seed(seed)

    Q = np.random.rand(n, n)
    return ctmc_makeinfgen(Q)


def weaklyconncomp(G: np.ndarray) -> ConnectedComponents:
    """
    Find weakly connected components in a directed graph.

    Args:
        G: Adjacency matrix of the graph

    Returns:
        ConnectedComponents with number of components and component assignment
    """
    G = np.asarray(G, dtype=np.float64)
    n = G.shape[0]

    # Make symmetric (undirected) for weakly connected components
    adj = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if (not np.isnan(G[i, j]) and G[i, j] != 0.0) or \
               (not np.isnan(G[j, i]) and G[j, i] != 0.0):
                adj[i, j] = 1.0
                adj[j, i] = 1.0

    # BFS-based connected component finding
    visited = np.zeros(n, dtype=bool)
    component = np.zeros(n, dtype=int)
    num_components = 0

    for start in range(n):
        if not visited[start]:
            num_components += 1
            queue = deque([start])
            visited[start] = True
            component[start] = num_components

            while queue:
                current = queue.popleft()
                for neighbor in range(n):
                    if not visited[neighbor] and adj[current, neighbor] != 0.0:
                        visited[neighbor] = True
                        component[neighbor] = num_components
                        queue.append(neighbor)

    return ConnectedComponents(num_components, component)


def ctmc_solve(Q: np.ndarray) -> np.ndarray:
    """
    Compute the equilibrium distribution of a continuous-time Markov chain.

    Args:
        Q: Infinitesimal generator matrix

    Returns:
        Equilibrium distribution vector
    """
    return ctmc_solve_full(Q).equilibrium_distribution


def ctmc_solve_full(Q: np.ndarray) -> CTMCSolveResult:
    """
    Compute the equilibrium distribution with full details.

    Args:
        Q: Infinitesimal generator matrix

    Returns:
        CTMCSolveResult with equilibrium distribution and component info
    """
    Q = np.asarray(Q, dtype=np.float64)
    n = Q.shape[0]

    # Handle trivial case
    if n == 1:
        return CTMCSolveResult(
            np.array([1.0]),
            Q,
            1,
            np.array([1])
        )

    # Normalize to valid infinitesimal generator
    normalized_Q = ctmc_makeinfgen(Q)

    # Check for connected components
    sym_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if np.abs(normalized_Q[i, j]) + np.abs(normalized_Q[j, i]) > 0:
                sym_matrix[i, j] = 1.0

    cc = weaklyconncomp(sym_matrix)

    if cc.num_components > 1:
        # Reducible generator - solve each component recursively
        p = np.zeros(n)
        for c in range(1, cc.num_components + 1):
            indices = np.where(cc.component_assignment == c)[0]
            m = len(indices)
            Qc = np.zeros((m, m))
            for ii, i in enumerate(indices):
                for jj, j in enumerate(indices):
                    Qc[ii, jj] = normalized_Q[i, j]
            pc = ctmc_solve(ctmc_makeinfgen(Qc))
            for ii, i in enumerate(indices):
                p[i] = pc[ii]

        # Normalize
        total = np.sum(p)
        if total > 0:
            p /= total

        return CTMCSolveResult(p, normalized_Q, cc.num_components, cc.component_assignment)

    # Check if all zero
    if np.all(normalized_Q == 0.0):
        p = np.ones(n) / n
        return CTMCSolveResult(p, normalized_Q, 1, np.ones(n, dtype=int))

    # Solve the system: p * Q = 0, sum(p) = 1
    # Modified to: Q' * p = b where b = [0,0,...,1]
    # and last ROW of Q' is replaced with ones (normalization constraint)

    Q_mod = normalized_Q.T.copy()

    # Replace last row with ones for normalization constraint
    Q_mod[n - 1, :] = 1.0

    # Build RHS vector
    b = np.zeros(n)
    b[n - 1] = 1.0

    # Solve using LU decomposition
    try:
        lu, piv = lu_factor(Q_mod)
        p = lu_solve((lu, piv), b)
    except Exception:
        # If singular, return uniform
        p = np.ones(n) / n

    # Ensure non-negative
    p = np.maximum(p, 0.0)

    # Renormalize
    total = np.sum(p)
    if total > 0:
        p /= total

    return CTMCSolveResult(p, normalized_Q, 1, np.ones(n, dtype=int))


def ctmc_timereverse(Q: np.ndarray) -> np.ndarray:
    """
    Compute the time-reversed generator of a CTMC.

    Args:
        Q: Infinitesimal generator matrix

    Returns:
        Infinitesimal generator matrix of the time-reversed process
    """
    Q = np.asarray(Q, dtype=np.float64)
    n = Q.shape[0]
    pie = ctmc_solve(Q)

    Q_rev = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i != j and pie[j] > 1e-14:
                Q_rev[j, i] = Q[i, j] * pie[i] / pie[j]

    # Normalize to valid infinitesimal generator (row sums = 0)
    return ctmc_makeinfgen(Q_rev)


def ctmc_randomization(Q: np.ndarray, q: float = None,
                       seed: int = None) -> Tuple[np.ndarray, float]:
    """
    Apply uniformization (randomization) to transform a CTMC into a DTMC.

    Args:
        Q: Infinitesimal generator matrix
        q: Uniformization rate (if None, uses max|diag(Q)| + random)
        seed: Random seed (optional)

    Returns:
        Tuple of (uniformized stochastic matrix, uniformization rate)
    """
    Q = np.asarray(Q, dtype=np.float64)
    n = Q.shape[0]

    if seed is not None:
        np.random.seed(seed)

    # Compute uniformization rate
    if q is None:
        q = np.max(np.abs(np.diag(Q))) + np.random.rand()

    # P = Q/q + I
    P = Q / q + np.eye(n)

    return dtmc_makestochastic(P), q


def ctmc_uniformization(pi0: np.ndarray, Q: np.ndarray, t: float,
                        tol: float = 1e-12, maxiter: int = 1000
                        ) -> Tuple[np.ndarray, int]:
    """
    Compute transient probabilities using uniformization method.

    Args:
        pi0: Initial probability distribution vector
        Q: Infinitesimal generator matrix
        t: Time point for transient analysis
        tol: Error tolerance (default: 1e-12)
        maxiter: Maximum iterations (default: 1000)

    Returns:
        Tuple of (probability distribution at time t, number of iterations)
    """
    Q = np.asarray(Q, dtype=np.float64)
    pi0 = np.asarray(pi0, dtype=np.float64)
    n = Q.shape[0]

    # For very large t, return equilibrium distribution
    if t > 1e6:
        return ctmc_solve(Q), 0

    # Uniformization rate
    max_diag = np.max(np.abs(np.diag(Q)))
    q = 1.1 * max_diag

    if q == 0 or t == 0:
        return pi0.copy(), 0

    # Uniformized matrix: Qs = I + Q/q
    Qs = np.eye(n) + Q / q

    # For large q*t, use more iterations and different approach
    qt = q * t

    # Compute using Fox-Glynn algorithm or direct summation
    # Find truncation point using Poisson tail bound
    if qt > 700:  # Avoid overflow in exp
        # For very large qt, just return equilibrium
        return ctmc_solve(Q), 0

    # Find number of terms needed using Poisson CDF
    # P(N > k) < tol where N ~ Poisson(qt)
    from math import ceil
    import scipy.stats as stats

    # Use inverse Poisson CDF to find truncation point
    kmax = int(ceil(stats.poisson.ppf(1 - tol, qt))) + 10
    kmax = min(kmax, maxiter)
    kmax = max(kmax, int(qt + 6 * np.sqrt(qt)))  # Ensure enough terms

    # Compute transient probability using matrix powers
    # pi(t) = sum_{k=0}^{inf} exp(-qt) * (qt)^k / k! * pi0 * Qs^k

    # Initialize
    pi = np.zeros(n)
    P = pi0.copy()  # Current power: pi0 * Qs^k

    # k=0 term
    poisson_weight = np.exp(-qt)  # exp(-qt)
    pi += poisson_weight * P

    for k in range(1, kmax + 1):
        P = P @ Qs
        poisson_weight *= (qt / k)
        pi += poisson_weight * P

        # Check convergence
        if poisson_weight < tol:
            break

    # Normalize to ensure valid probability distribution
    total = np.sum(pi)
    if total > 0:
        pi = pi / total

    return pi, min(k, kmax)


def ctmc_relsolve(Q: np.ndarray, refstate: int = 0) -> np.ndarray:
    """
    Compute the equilibrium distribution relative to a reference state.

    Solves global balance equations with p(refstate) = 1 normalization.

    Args:
        Q: Infinitesimal generator matrix
        refstate: Index of reference state (0-based, default: 0)

    Returns:
        Relative equilibrium distribution vector
    """
    Q = np.asarray(Q, dtype=np.float64)
    n = Q.shape[0]

    if n == 1:
        return np.array([1.0])

    normalized_Q = ctmc_makeinfgen(Q)

    # Modify system: set ROW refstate to have 1 at refstate (normalization constraint)
    Q_mod = normalized_Q.T.copy()

    # Replace reference row (not column!)
    Q_mod[refstate, :] = 0.0
    Q_mod[refstate, refstate] = 1.0

    # Build RHS vector
    b = np.zeros(n)
    b[refstate] = 1.0

    # Solve
    try:
        lu, piv = lu_factor(Q_mod)
        p = lu_solve((lu, piv), b)
    except Exception:
        p = np.ones(n) / n

    return p


# ============================================================================
# DTMC Functions
# ============================================================================

def dtmc_makestochastic(P: np.ndarray) -> np.ndarray:
    """
    Normalize a matrix to be a valid stochastic matrix.

    Rescales rows to sum to 1. Rows with zero sum get 1 at diagonal.

    Args:
        P: Input matrix

    Returns:
        Valid stochastic transition matrix
    """
    P = np.asarray(P, dtype=np.float64)
    n = P.shape[0]
    result = np.zeros((n, n))

    for i in range(n):
        row_sum = np.sum(P[i, :])

        if row_sum > 0:
            # Normalize row
            result[i, :] = P[i, :] / row_sum
            # Ensure diagonal adjustment for numerical stability
            off_diag_sum = np.sum(result[i, :]) - result[i, i]
            result[i, i] = max(0.0, min(1.0, 1.0 - off_diag_sum))
        else:
            # Set absorbing state
            result[i, i] = 1.0

    return result


def dtmc_isfeasible(P: np.ndarray) -> int:
    """
    Check feasibility of a stochastic matrix.

    Verifies if row sums are close to 1 and elements are non-negative.

    Args:
        P: Matrix to check

    Returns:
        Precision level (1-15) if feasible, or 0 if not feasible
    """
    P = np.asarray(P, dtype=np.float64)
    n = P.shape[0]

    # Compute row sums and min element
    row_sums = np.sum(P, axis=1)
    min_element = np.min(P)
    min_row_sum = np.min(row_sums)
    max_row_sum = np.max(row_sums)

    result = 0
    for tol_exp in range(1, 16):
        tolerance = 10.0 ** (-tol_exp)
        if (min_row_sum > 1 - tolerance and
            max_row_sum < 1 + tolerance and
            min_element > -tolerance):
            result = tol_exp

    return result


def dtmc_solve(P: np.ndarray) -> np.ndarray:
    """
    Compute the equilibrium distribution of a discrete-time Markov chain.

    Args:
        P: Stochastic transition matrix

    Returns:
        Equilibrium distribution vector
    """
    P = np.asarray(P, dtype=np.float64)
    n = P.shape[0]

    # Convert DTMC to CTMC: Q = P - I
    Q = P - np.eye(n)

    return ctmc_solve(Q)


def dtmc_rand(n: int, seed: int = None) -> np.ndarray:
    """
    Generate a random stochastic transition matrix.

    Args:
        n: Size of the matrix (n x n)
        seed: Random seed (optional)

    Returns:
        Random stochastic transition matrix
    """
    P, _ = ctmc_randomization(ctmc_rand(n, seed), seed=seed)
    return P


def dtmc_simulate(P: np.ndarray, pi0: np.ndarray, n_steps: int,
                  seed: int = None) -> np.ndarray:
    """
    Simulate a trajectory of a discrete-time Markov chain.

    Args:
        P: Stochastic transition matrix
        pi0: Initial probability distribution vector
        n_steps: Number of steps to simulate
        seed: Random seed (optional)

    Returns:
        Vector of state indices visited in the simulation (0-based)
    """
    P = np.asarray(P, dtype=np.float64)
    pi0 = np.asarray(pi0, dtype=np.float64)
    n = P.shape[0]

    if seed is not None:
        np.random.seed(seed)

    states = np.zeros(n_steps, dtype=int)

    # Sample initial state from pi0
    cum_sum = np.cumsum(pi0)
    rnd = np.random.rand()
    current_state = np.searchsorted(cum_sum, rnd)
    current_state = min(current_state, n - 1)

    # Precompute cumulative transition probabilities
    cum_P = np.cumsum(P, axis=1)

    # Simulate trajectory
    for step in range(n_steps):
        states[step] = current_state

        # Check for absorbing state
        if P[current_state, current_state] == 1.0:
            states[step:] = current_state
            break

        # Sample next state
        rnd = np.random.rand()
        current_state = np.searchsorted(cum_P[current_state], rnd)
        current_state = min(current_state, n - 1)

    return states


def dtmc_stochcomp(P: np.ndarray, I: np.ndarray = None) -> np.ndarray:
    """
    Compute the stochastic complement of a DTMC partition.

    Args:
        P: Stochastic transition matrix
        I: Indices of states to retain (0-based, default: first half)

    Returns:
        Stochastic complement matrix for the subset I
    """
    P = np.asarray(P, dtype=np.float64)
    n = P.shape[0]

    # Default: first half of states
    if I is None:
        I = np.arange((n + 1) // 2)
    else:
        I = np.asarray(I)

    # Complement of I
    Ic = np.array([i for i in range(n) if i not in I])

    m1 = len(I)
    m2 = len(Ic)

    if m2 == 0:
        return P[np.ix_(I, I)]

    # Extract submatrices
    P11 = P[np.ix_(I, I)]
    P12 = P[np.ix_(I, Ic)]
    P21 = P[np.ix_(Ic, I)]
    P22 = P[np.ix_(Ic, Ic)]

    # S = P11 + P12 * (I - P22)^{-1} * P21
    try:
        I_minus_P22 = np.eye(m2) - P22
        inv_I_minus_P22 = inv(I_minus_P22)
        S = P11 + P12 @ inv_I_minus_P22 @ P21
    except Exception:
        return P11

    return S


def dtmc_timereverse(P: np.ndarray) -> np.ndarray:
    """
    Compute the time-reversed transition matrix of a DTMC.

    Args:
        P: Stochastic transition matrix of the original process

    Returns:
        Stochastic transition matrix of the time-reversed process
    """
    P = np.asarray(P, dtype=np.float64)
    n = P.shape[0]
    pie = dtmc_solve(P)

    P_rev = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if pie[j] > 1e-14:
                P_rev[j, i] = P[i, j] * pie[i] / pie[j]

    # Normalize to valid stochastic matrix (row sums = 1)
    return dtmc_makestochastic(P_rev)


def dtmc_uniformization(pi0: np.ndarray, P: np.ndarray, t: float = 1e4,
                        tol: float = 1e-12, maxiter: int = 100
                        ) -> Tuple[np.ndarray, int]:
    """
    Compute transient probabilities for a DTMC using uniformization.

    Args:
        pi0: Initial probability distribution vector
        P: Transition probability matrix
        t: Time point for transient analysis (default: 1e4)
        tol: Error tolerance (default: 1e-12)
        maxiter: Maximum iterations (default: 100)

    Returns:
        Tuple of (probability distribution at time t, number of iterations)
    """
    P = np.asarray(P, dtype=np.float64)
    n = P.shape[0]

    # Convert to CTMC generator: Q = P - I
    Q = P - np.eye(n)

    return ctmc_uniformization(pi0, ctmc_makeinfgen(Q), t, tol, maxiter)


__all__ = [
    # CTMC functions
    'ctmc_makeinfgen',
    'ctmc_solve',
    'ctmc_solve_full',
    'ctmc_rand',
    'ctmc_timereverse',
    'ctmc_randomization',
    'ctmc_uniformization',
    'ctmc_relsolve',
    'weaklyconncomp',
    # DTMC functions
    'dtmc_makestochastic',
    'dtmc_isfeasible',
    'dtmc_solve',
    'dtmc_rand',
    'dtmc_simulate',
    'dtmc_stochcomp',
    'dtmc_timereverse',
    'dtmc_uniformization',
    # Data classes
    'ConnectedComponents',
    'CTMCSolveResult',
]
