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
            if (np.isnan(G[i, j]) or G[i, j] != 0.0) or \
               (np.isnan(G[j, i]) or G[j, i] != 0.0):
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

    # Check if system is singular before solving (matches MATLAB backslash behavior)
    # MATLAB's \ returns NaN for singular systems; numpy solve may return wrong answers
    if np.linalg.matrix_rank(Q_mod) < n:
        p = np.full(n, np.nan)
    else:
        try:
            lu, piv = lu_factor(Q_mod)
            p = lu_solve((lu, piv), b)
        except Exception:
            p = np.full(n, np.nan)

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

    # Modify system: replace last row of Q^T (= last column of Q) with normalization constraint
    # MATLAB: Qnnz(:,end) = 0; Qnnz(refstate,end) = 1; bnnz(end) = 1; then solves Qnnz' \ bnnz
    Q_mod = normalized_Q.T.copy()

    # Replace last row of Q^T
    last_row = n - 1
    Q_mod[last_row, :] = 0.0
    Q_mod[last_row, refstate] = 1.0

    # Build RHS vector
    b = np.zeros(n)
    b[last_row] = 1.0

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


def ctmc_transient(Q, pi0=None, t0=None, t1=None, useStiff=False, reltol=1e-3, timestep=None):
    """
    Transient analysis of a continuous-time Markov chain using ODE solvers.

    Computes the transient state probabilities by solving the ODE system
    dpi/dt = pi * Q using scipy ODE solvers.

    Matches MATLAB: matlab/lib/kpctoolbox/mc/ctmc_transient.m

    Calling conventions (matching MATLAB):
        ctmc_transient(Q, t1)             - uniform pi0, t0=0
        ctmc_transient(Q, pi0, t1)        - given pi0, t0=0
        ctmc_transient(Q, pi0, t0, t1)    - given pi0 and time range

    Parameters
    ----------
    Q : array_like
        Infinitesimal generator matrix.
    pi0 : array_like or float, optional
        Initial probability distribution vector. If a scalar is given and
        t0/t1 are not provided, it is interpreted as t1 (MATLAB nargin==2).
    t0 : float, optional
        Start time. Default is 0.
    t1 : float, optional
        End time.
    useStiff : bool, optional
        If True, use stiff ODE solver (BDF, analogous to MATLAB ode15s).
        Default is False (uses RK23, analogous to MATLAB ode23).
    reltol : float, optional
        Relative tolerance for the ODE solver. Default is 1e-3.
    timestep : float or None, optional
        Fixed time step for the output. If None, uses adaptive stepping.

    Returns
    -------
    pi : numpy.ndarray
        State probability vectors at each time point, shape (len(t), n).
    t : numpy.ndarray
        Time points corresponding to the probabilities, shape (len(t),).
    """
    from scipy.integrate import solve_ivp

    Q = np.asarray(Q, dtype=np.float64)
    n = Q.shape[0]

    # Handle MATLAB-style calling conventions
    if t1 is None and t0 is None:
        # nargin==2: ctmc_transient(Q, t1) => pi0 is actually t1
        t1_val = float(pi0)
        t0_val = 0.0
        pi0_vec = np.ones(n) / n
    elif t1 is None:
        # nargin==3: ctmc_transient(Q, pi0, t1) => t0 is actually t1
        pi0_vec = np.asarray(pi0, dtype=np.float64).ravel()
        t1_val = float(t0)
        t0_val = 0.0
    else:
        # nargin>=4: ctmc_transient(Q, pi0, t0, t1)
        pi0_vec = np.asarray(pi0, dtype=np.float64).ravel()
        t0_val = float(t0)
        t1_val = float(t1)

    # ODE: dpi/dt = pi * Q
    def ode_func(t, pi):
        return (pi.reshape(1, -1) @ Q).ravel()

    method = 'BDF' if useStiff else 'RK23'

    if timestep is not None and timestep > 0:
        # Fixed time steps
        t_eval = np.arange(t0_val, t1_val, timestep)
        if len(t_eval) == 0 or t_eval[-1] != t1_val:
            t_eval = np.append(t_eval, t1_val)

        sol = solve_ivp(
            ode_func, [t0_val, t1_val], pi0_vec,
            method=method, t_eval=t_eval, rtol=reltol
        )
        t_out = sol.t
        pi_out = sol.y.T
    else:
        # Adaptive stepping
        if np.isinf(t1_val):
            nonZeroRates = np.abs(Q[Q != 0])
            nonZeroRates = nonZeroRates[nonZeroRates > 0]
            T_end = abs(100 / np.min(nonZeroRates))
            sol = solve_ivp(
                ode_func, [t0_val, T_end], pi0_vec,
                method=method, rtol=reltol
            )
        else:
            sol = solve_ivp(
                ode_func, [t0_val, t1_val], pi0_vec,
                method=method, rtol=reltol
            )
        t_out = sol.t
        pi_out = sol.y.T

    return pi_out, t_out


def stronglyconncomp(P: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Find strongly connected components and classify as recurrent/transient.

    Uses Tarjan's algorithm (via scipy) to find SCCs of a directed graph
    defined by a transition matrix or generator matrix.

    Args:
        P: Transition probability matrix (DTMC) or adjacency matrix.
           For generators, pass (Q with diagonal zeroed and positive entries).

    Returns:
        Tuple of:
            scc: 1-based SCC label for each state (1-indexed for MATLAB parity)
            isrec: Boolean array where isrec[i] is True if SCC i is recurrent
    """
    from scipy.sparse import csc_matrix
    from scipy.sparse.csgraph import connected_components

    P = np.asarray(P, dtype=np.float64)
    n = P.shape[0]

    # Build adjacency: entry (i,j) is 1 if transition i->j is possible
    adj = (P > 0).astype(np.float64)
    np.fill_diagonal(adj, 0.0)

    n_comp, labels = connected_components(
        csc_matrix(adj), directed=True, connection='strong',
        return_labels=True
    )

    # Convert to 1-based indexing (MATLAB convention)
    scc = labels + 1

    # Classify: an SCC is recurrent if it has no outgoing edges to other SCCs
    isrec = np.zeros(n_comp, dtype=bool)
    for c in range(n_comp):
        states_c = np.where(labels == c)[0]
        outgoing = 0.0
        for s in states_c:
            for c2 in range(n_comp):
                if c2 != c:
                    states_c2 = np.where(labels == c2)[0]
                    outgoing += np.sum(adj[s, states_c2])
        isrec[c] = outgoing < 1e-10

    return scc, isrec


def ctmc_simulate(Q: np.ndarray, pi0: np.ndarray, n_events: int,
                  seed: int = None) -> Tuple[np.ndarray, np.ndarray]:
    """
    Simulate a trajectory of a continuous-time Markov chain.

    Generates a sample path by repeatedly sampling sojourn times
    (exponential with rate -Q[st,st]) and next states.

    Args:
        Q: Infinitesimal generator matrix
        pi0: Initial probability distribution vector (if empty, use uniform)
        n_events: Number of transition events to simulate
        seed: Random seed (optional)

    Returns:
        Tuple of:
            sojourn_times: Sojourn time in each state visited
            states: Sequence of state indices visited (0-based)

    MATLAB: matlab/src/api/mc/ctmc_simulate.m
    """
    Q = np.asarray(Q, dtype=np.float64)
    n = Q.shape[0]

    if seed is not None:
        np.random.seed(seed)

    if pi0 is None or len(pi0) == 0:
        pi0 = np.random.rand(n)
        pi0 /= pi0.sum()
    else:
        pi0 = np.asarray(pi0, dtype=np.float64)

    # Sample initial state from pi0
    cum = np.cumsum(pi0)
    st = np.searchsorted(cum, np.random.rand())
    st = min(st, n - 1)

    # Build cumulative transition probabilities from off-diagonal rates
    Q_off = Q.copy()
    np.fill_diagonal(Q_off, 0.0)
    F = np.cumsum(Q_off, axis=1)
    # Normalize rows
    for i in range(n):
        row_total = F[i, -1]
        if row_total > 0:
            F[i, :] /= row_total

    sojourn_times = np.zeros(n_events)
    states = np.zeros(n_events, dtype=int)

    for i in range(n_events):
        states[i] = st
        # Sojourn time: exponential with rate -Q[st,st]
        rate = -Q[st, st]
        if rate > 0:
            sojourn_times[i] = np.random.exponential(1.0 / rate)
        else:
            sojourn_times[i] = np.inf

        # Next state: sample from cumulative transition probs
        r = np.random.rand()
        next_st = np.searchsorted(F[st, :], r)
        next_st = min(next_st, n - 1)
        st = next_st

    return sojourn_times, states


def ctmc_stochcomp(Q: np.ndarray, I: np.ndarray = None) -> np.ndarray:
    """
    Compute the stochastic complement of a CTMC generator.

    For a generator partitioned as:
        Q = [[Q11, Q12],
             [Q21, Q22]]
    the stochastic complement for subset I is:
        S = Q11 + Q12 * (-Q22)^{-1} * Q21

    Args:
        Q: Infinitesimal generator matrix
        I: Indices of states to retain (0-based). Default: first half.

    Returns:
        Stochastic complement generator for the subset I.

    MATLAB: matlab/src/api/mc/ctmc_stochcomp.m
    """
    Q = np.asarray(Q, dtype=np.float64)
    n = Q.shape[0]

    if I is None:
        I = np.arange((n + 1) // 2)
    else:
        I = np.asarray(I)

    # Complement of I
    Ic = np.array([i for i in range(n) if i not in I])

    if len(Ic) == 0:
        return Q[np.ix_(I, I)]

    Q11 = Q[np.ix_(I, I)]
    Q12 = Q[np.ix_(I, Ic)]
    Q21 = Q[np.ix_(Ic, I)]
    Q22 = Q[np.ix_(Ic, Ic)]

    try:
        inv_neg_Q22 = inv(-Q22)
        S = Q11 + Q12 @ inv_neg_Q22 @ Q21
    except Exception:
        S = Q11

    return S


def ctmc_courtois(Q: np.ndarray, MS: list,
                  q: float = None) -> Tuple:
    """
    Courtois decomposition for nearly completely decomposable CTMCs.

    Solves a CTMC by decomposing it into macrostates, solving each
    macrostate independently (micro-probabilities), then computing
    macro-probabilities and combining them.

    Args:
        Q: Infinitesimal generator matrix
        MS: List of lists, where MS[i] contains 0-based state indices
            in macrostate i
        q: Randomization coefficient (optional; auto-computed if None)

    Returns:
        Tuple of (p, Qperm, Qdec, eps, epsMAX, P, B, C, q) where:
            p: Approximate steady-state probability vector
            Qperm: Q reordered by macrostates
            Qdec: Block-diagonal decomposed generator
            eps: NCD (nearly complete decomposability) index
            epsMAX: Maximum acceptable eps
            P: Randomized probability matrix
            B: Coupling part of P
            C: Difference matrix (Qperm - Qdec, set to 0)
            q: Randomization coefficient used

    MATLAB: matlab/src/api/mc/ctmc_courtois.m
    """
    Q = np.asarray(Q, dtype=np.float64)
    n = Q.shape[0]
    nMacroStates = len(MS)

    # Rearrange generator according to macrostates
    v = []
    for ms in MS:
        v.extend(ms)
    v = np.array(v, dtype=int)

    Qperm = Q[np.ix_(v, v)]
    Qdec = Qperm.copy()

    # Zero out off-block-diagonal entries
    procRows = 0
    for i in range(nMacroStates):
        ms_len = len(MS[i])
        if procRows > 0:
            Qdec[procRows:procRows + ms_len, :procRows] = 0.0
        Qdec[procRows:procRows + ms_len, procRows + ms_len:] = 0.0
        procRows += ms_len

    # Make each diagonal block a valid generator
    Qdec = ctmc_makeinfgen(Qdec)

    # NCD error
    C_mat = 0
    eps_val = 1.0

    # Randomization
    if q is None:
        q = 1.05 * np.max(np.abs(Qperm))

    P, _ = ctmc_randomization(Qperm, q)
    A = P.copy()

    procRows = 0
    for i in range(nMacroStates):
        ms_len = len(MS[i])
        if procRows > 0:
            A[procRows:procRows + ms_len, :procRows] = 0.0
        A[procRows:procRows + ms_len, procRows + ms_len:] = 0.0
        procRows += ms_len

    B_mat = P - A
    eps_val = np.max(np.sum(B_mat, axis=0))

    # Compute epsMAX from second-largest eigenvalues of diagonal blocks
    A_stoch = A.copy()
    procRows = 0
    for i in range(nMacroStates):
        ms_len = len(MS[i])
        for j in range(ms_len):
            row = procRows + j
            others = list(range(procRows, procRows + ms_len))
            others.remove(row)
            A_stoch[row, row] = 1.0 - np.sum(A_stoch[row, others])
        procRows += ms_len

    eigMS = np.zeros(nMacroStates)
    procRows = 0
    for i in range(nMacroStates):
        ms_len = len(MS[i])
        block = A_stoch[procRows:procRows + ms_len, procRows:procRows + ms_len]
        eigs = np.sort(np.abs(np.linalg.eigvals(block)))
        if len(eigs) > 1:
            eigMS[i] = eigs[-2]
        procRows += ms_len

    epsMAX = (1.0 - np.max(eigMS)) / 2.0

    # Compute micro-probabilities (solve each block)
    pmicro = np.zeros(n)
    procRows = 0
    for i in range(nMacroStates):
        ms_len = len(MS[i])
        Qblock = Qdec[procRows:procRows + ms_len, procRows:procRows + ms_len]
        try:
            pi_block = ctmc_solve(Qblock)
        except Exception:
            pi_block = np.ones(ms_len) / ms_len
        pmicro[procRows:procRows + ms_len] = pi_block
        procRows += ms_len

    # Compute macro-probabilities
    G = np.zeros((nMacroStates, nMacroStates))
    procRows = 0
    for i in range(nMacroStates):
        ms_len_i = len(MS[i])
        procCols = 0
        for j in range(nMacroStates):
            ms_len_j = len(MS[j])
            if i != j:
                for iState in range(ms_len_i):
                    G[i, j] += pmicro[procRows + iState] * np.sum(
                        P[procRows + iState, procCols:procCols + ms_len_j]
                    )
            procCols += ms_len_j
        procRows += ms_len_i

    for i in range(nMacroStates):
        G[i, i] = 1.0 - np.sum(G[i, :])

    pMacro = dtmc_solve(G)

    # Combine micro and macro
    p_perm = np.zeros(n)
    procRows = 0
    for i in range(nMacroStates):
        ms_len = len(MS[i])
        p_perm[procRows:procRows + ms_len] = pMacro[i] * pmicro[procRows:procRows + ms_len]
        procRows += ms_len

    # Re-order back to original state ordering
    p = np.zeros(n)
    for i, vi in enumerate(v):
        p[vi] = p_perm[i]

    return p, Qperm, Qdec, eps_val, epsMAX, P, B_mat, C_mat, q


def ctmc_kms(Q: np.ndarray, MS: list,
             numSteps: int) -> Tuple:
    """
    Koury-McAllister-Stewart aggregation-disaggregation method.

    Iterative method that starts from the Courtois decomposition
    solution and refines it through aggregation-disaggregation steps.

    Args:
        Q: Infinitesimal generator matrix
        MS: List of lists, where MS[i] contains 0-based state indices
            in macrostate i
        numSteps: Number of iterative steps

    Returns:
        Tuple of (p, p_1, Qperm, eps, epsMAX, pcourt) where:
            p: Estimated steady-state probability vector
            p_1: Previous iteration result
            Qperm: Permuted Q from courtois
            eps: NCD index
            epsMAX: Maximum acceptable eps
            pcourt: Courtois initial solution

    MATLAB: matlab/src/api/mc/ctmc_kms.m
    """
    Q = np.asarray(Q, dtype=np.float64)
    nMacroStates = len(MS)
    nStates = Q.shape[0]

    # Start from Courtois solution
    pcourt, Qperm, Qdec, eps_val, epsMAX, P, B_mat, C_mat, q_val = ctmc_courtois(Q, MS)

    pn = pcourt.copy()
    p_1 = pn.copy()

    # Reorder MS indices for permuted state space
    v = []
    for ms in MS:
        v.extend(ms)

    for step in range(numSteps):
        pn_1 = pn.copy()
        p_1 = pn_1.copy()

        # Aggregation step: compute conditional probabilities
        pcondn_1 = pn_1.copy()
        procRows = 0
        for I in range(nMacroStates):
            ms_len = len(MS[I])
            block_sum = np.sum(pn_1[procRows:procRows + ms_len])
            if block_sum > 0:
                pcondn_1[procRows:procRows + ms_len] = (
                    pn_1[procRows:procRows + ms_len] / block_sum
                )
            procRows += ms_len

        G = np.zeros((nMacroStates, nMacroStates))
        procCols = 0
        for I in range(nMacroStates):
            ms_len_I = len(MS[I])
            procRows = 0
            for J in range(nMacroStates):
                ms_len_J = len(MS[J])
                # G(I,J) = pcondn_1_J^T * P(J, I) * ones(I)
                pblock = pcondn_1[procRows:procRows + ms_len_J]
                Pblock = P[procRows:procRows + ms_len_J, procCols:procCols + ms_len_I]
                G[I, J] = pblock @ Pblock @ np.ones(ms_len_I)
                procRows += ms_len_J
            procCols += ms_len_I

        w = dtmc_solve(G.T)

        # Disaggregation step
        z = pcondn_1.copy()
        L_mat = np.zeros((nStates, nStates))
        D_mat = np.zeros((nStates, nStates))
        U_mat = np.zeros((nStates, nStates))

        zn = np.zeros(nStates)
        procRows = 0
        for I in range(nMacroStates):
            ms_len_I = len(MS[I])
            zn[procRows:procRows + ms_len_I] = w[I] * z[procRows:procRows + ms_len_I]

            procCols = 0
            for J in range(nMacroStates):
                ms_len_J = len(MS[J])
                rows = slice(procRows, procRows + ms_len_I)
                cols = slice(procCols, procCols + ms_len_J)
                if I > J:
                    L_mat[rows, cols] = P[rows, cols]
                elif I == J:
                    D_mat[rows, cols] = np.eye(ms_len_I) - P[rows, cols]
                else:
                    U_mat[rows, cols] = P[rows, cols]
                procCols += ms_len_J
            procRows += ms_len_I

        M = D_mat - U_mat
        MPL = (nStates - 2) // 2

        A_block = M[:MPL + 1, :MPL + 1]
        B_block = M[:MPL + 1, MPL + 1:]
        C_block = M[MPL + 1:, MPL + 1:]

        try:
            invA = inv(A_block)
            invC = inv(C_block)
            top = np.hstack([invA, -invA @ B_block @ invC])
            bot_rows = M.shape[0] - (MPL + 1)
            bottom = np.hstack([np.zeros((bot_rows, MPL + 1)), invC])
            combined = np.vstack([top, bottom])
            pn = zn @ L_mat @ combined
        except Exception:
            pn = zn

    # Re-order back to original state ordering
    p_out = np.zeros(nStates)
    v_arr = np.array(v, dtype=int)
    for i in range(len(v_arr)):
        p_out[v_arr[i]] = pn[i]

    return p_out, p_1, Qperm, eps_val, epsMAX, pcourt


def ctmc_takahashi(Q: np.ndarray, MS: list,
                   numSteps: int) -> Tuple:
    """
    Takahashi aggregation-disaggregation method.

    Iterative refinement starting from Courtois decomposition. Uses
    Gauss-Seidel style disaggregation.

    Args:
        Q: Infinitesimal generator matrix
        MS: List of lists, where MS[i] contains 0-based state indices
            in macrostate i
        numSteps: Number of iterative steps

    Returns:
        Tuple of (p, p_1, pcourt, Qperm, eps, epsMAX) where:
            p: Estimated steady-state probability vector
            p_1: Previous iteration result
            pcourt: Courtois initial solution
            Qperm: Permuted Q from courtois
            eps: NCD index
            epsMAX: Maximum acceptable eps

    MATLAB: matlab/src/api/mc/ctmc_takahashi.m
    """
    Q = np.asarray(Q, dtype=np.float64)
    nMacroStates = len(MS)
    nStates = Q.shape[0]

    # Start from Courtois solution
    pcourt, Qperm, Qdec, eps_val, epsMAX, _, _, _, q_val = ctmc_courtois(Q, MS)
    P, _ = ctmc_randomization(Q)

    pn = pcourt.copy()
    p_1 = pn.copy()

    # Reorder MS indices
    v = []
    for ms in MS:
        v.extend(ms)

    for step in range(numSteps):
        pn_1 = pn.copy()
        p_1 = pn_1.copy()

        # Aggregation step
        G = np.zeros((nMacroStates, nMacroStates))
        procRows = 0
        for I in range(nMacroStates):
            ms_len_I = len(MS[I])
            S = np.sum(pn_1[procRows:procRows + ms_len_I])
            procCols = 0
            for J in range(nMacroStates):
                ms_len_J = len(MS[J])
                if I != J:
                    for ii in range(ms_len_I):
                        for jj in range(ms_len_J):
                            if S > 1e-14:
                                G[I, J] += P[procRows + ii, procCols + jj] * pn_1[procRows + ii] / S
                procCols += ms_len_J
            procRows += ms_len_I

        for i in range(nMacroStates):
            G[i, i] = 1.0 - np.sum(G[i, :])

        gamma = dtmc_solve(G)

        # Disaggregation step
        procRows = 0
        GI = np.zeros((nMacroStates, nStates))
        for I in range(nMacroStates):
            ms_len_I = len(MS[I])
            S = np.sum(pn_1[procRows:procRows + ms_len_I])
            if S > 1e-14:
                for j in range(nStates):
                    GI[I, j] += np.sum(
                        P[procRows:procRows + ms_len_I, j] *
                        pn_1[procRows:procRows + ms_len_I]
                    ) / S
            procRows += ms_len_I

        procRows = 0
        for I in range(nMacroStates):
            ms_len_I = len(MS[I])
            A_block = np.eye(ms_len_I)
            b_vec = np.zeros(ms_len_I)

            for ii in range(ms_len_I):
                for jj in range(ms_len_I):
                    A_block[ii, jj] -= P[procRows + jj, procRows + ii]
                for K in range(nMacroStates):
                    if K != I:
                        b_vec[ii] += gamma[K] * GI[K, procRows + ii]

            try:
                pn[procRows:procRows + ms_len_I] = np.linalg.solve(A_block, b_vec)
            except Exception:
                pn[procRows:procRows + ms_len_I] = np.linalg.lstsq(A_block, b_vec, rcond=None)[0]
            procRows += ms_len_I

        # Normalize
        total = np.sum(pn)
        if total > 0:
            pn /= total

    return pn, p_1, pcourt, Qperm, eps_val, epsMAX


def ctmc_multi(Q: np.ndarray, MS: list, MSS: list) -> Tuple:
    """
    Multigrid aggregation-disaggregation method.

    Two-level multigrid: decomposes into macrostates (MS), then solves
    the macrostate chain using another Courtois decomposition (MSS).

    Args:
        Q: Infinitesimal generator matrix
        MS: List of lists, MS[i] = state indices in macrostate i
        MSS: List of lists, MSS[i] = macrostate indices in macro-macrostate i
             (indices into the G matrix, 0-based)

    Returns:
        Tuple of (p, p_1, pcourt, Qperm, eps, epsMAX) where:
            p: Approximate steady-state probability vector
            p_1: (unused, same as p)
            pcourt: Courtois-level solution
            Qperm: Permuted Q
            eps: NCD index
            epsMAX: Maximum acceptable eps

    MATLAB: matlab/src/api/mc/ctmc_multi.m
    """
    Q = np.asarray(Q, dtype=np.float64)
    n = Q.shape[0]
    nMacroStates = len(MS)

    # Rearrange generator according to macrostates
    v = []
    for ms in MS:
        v.extend(ms)
    v_arr = np.array(v, dtype=int)

    Qperm = Q[np.ix_(v_arr, v_arr)]
    Qdec = Qperm.copy()

    procRows = 0
    for i in range(nMacroStates):
        ms_len = len(MS[i])
        if procRows > 0:
            Qdec[procRows:procRows + ms_len, :procRows] = 0.0
        Qdec[procRows:procRows + ms_len, procRows + ms_len:] = 0.0
        procRows += ms_len

    Qdec = ctmc_makeinfgen(Qdec)

    # Randomization
    q_val = 1.05 * np.max(np.abs(Qperm))
    P, _ = ctmc_randomization(Qperm, q_val)
    A = P.copy()

    procRows = 0
    for i in range(nMacroStates):
        ms_len = len(MS[i])
        if procRows > 0:
            A[procRows:procRows + ms_len, :procRows] = 0.0
        A[procRows:procRows + ms_len, procRows + ms_len:] = 0.0
        procRows += ms_len

    B_mat = P - A
    eps_val = np.max(np.sum(B_mat, axis=0))

    # epsMAX
    A_stoch = A.copy()
    procRows = 0
    for i in range(nMacroStates):
        ms_len = len(MS[i])
        for j in range(ms_len):
            row = procRows + j
            others = list(range(procRows, procRows + ms_len))
            others.remove(row)
            A_stoch[row, row] = 1.0 - np.sum(A_stoch[row, others])
        procRows += ms_len

    eigMS = np.zeros(nMacroStates)
    procRows = 0
    for i in range(nMacroStates):
        ms_len = len(MS[i])
        block = A_stoch[procRows:procRows + ms_len, procRows:procRows + ms_len]
        eigs = np.sort(np.abs(np.linalg.eigvals(block)))
        if len(eigs) > 1:
            eigMS[i] = eigs[-2]
        procRows += ms_len

    epsMAX = (1.0 - np.max(eigMS)) / 2.0

    # Micro-probabilities
    pmicro = np.zeros(n)
    procRows = 0
    for i in range(nMacroStates):
        ms_len = len(MS[i])
        Qblock = Qdec[procRows:procRows + ms_len, procRows:procRows + ms_len]
        try:
            pi_block = ctmc_solve(Qblock)
        except Exception:
            pi_block = np.ones(ms_len) / ms_len
        pmicro[procRows:procRows + ms_len] = pi_block
        procRows += ms_len

    # Macro-probabilities
    G = np.zeros((nMacroStates, nMacroStates))
    procRows = 0
    for i in range(nMacroStates):
        ms_len_i = len(MS[i])
        procCols = 0
        for j in range(nMacroStates):
            ms_len_j = len(MS[j])
            if i != j:
                for iState in range(ms_len_i):
                    G[i, j] += pmicro[procRows + iState] * np.sum(
                        P[procRows + iState, procCols:procCols + ms_len_j]
                    )
            procCols += ms_len_j
        procRows += ms_len_i

    for i in range(nMacroStates):
        G[i, i] = 1.0 - np.sum(G[i, :])

    # Solve macrostate chain using another Courtois decomposition
    pMacro, _, _, _, _, _, _, _, _ = ctmc_courtois(
        ctmc_makeinfgen(G - np.eye(nMacroStates)), MSS
    )

    p_perm = np.zeros(n)
    procRows = 0
    for i in range(nMacroStates):
        ms_len = len(MS[i])
        p_perm[procRows:procRows + ms_len] = pMacro[i] * pmicro[procRows:procRows + ms_len]
        procRows += ms_len

    # Re-order back
    p = np.zeros(n)
    for i, vi in enumerate(v_arr):
        p[vi] = p_perm[i]

    pcourt = p  # single-level Courtois is the same in this 2-level

    return p, p, pcourt, Qperm, eps_val, epsMAX


def ctmc_solve_reducible(Q: np.ndarray, pi0: np.ndarray = None) -> np.ndarray:
    """
    Solve a reducible CTMC by converting to DTMC via randomization.

    Converts the CTMC to a DTMC using randomization, then calls
    dtmc_solve_reducible.

    Args:
        Q: Infinitesimal generator matrix (possibly reducible)
        pi0: Initial distribution vector (optional)

    Returns:
        Steady-state probability vector

    MATLAB: matlab/src/api/mc/ctmc_solve_reducible.m
    """
    P, _ = ctmc_randomization(Q)
    pi, _, _, _, _ = dtmc_solve_reducible(P, pi0)
    return pi


def ctmc_solve_reducible_blkdecomp(Q: np.ndarray, pi0: np.ndarray = None,
                                    tol: float = 1e-12) -> Tuple:
    """
    Solve a reducible CTMC using direct block decomposition.

    Algorithm:
      1. Decompose states into transient and recurrent classes via SCC
      2. For transient states: solve n * Q_tt = -p0_t for expected sojourn
      3. Compute hitting probabilities: h = n * Q_ta + p0_r
      4. For each recurrent class: solve pi_c * Q_cc = 0, scale by hitting prob

    Args:
        Q: Infinitesimal generator matrix
        pi0: Initial distribution (optional)
        tol: Numerical tolerance

    Returns:
        Tuple of (pi, pis, pi0_out, scc, isrec) where:
            pi: Limiting distribution
            pis: Per-SCC limiting distributions (numSCC x N)
            pi0_out: Starting distribution for each SCC
            scc: SCC label per state (1-based)
            isrec: Boolean array of recurrent SCCs

    MATLAB: matlab/src/api/mc/ctmc_solve_reducible_blkdecomp.m
    """
    Q = np.asarray(Q, dtype=np.float64)
    N = Q.shape[0]

    # Ensure valid generator
    Q = ctmc_makeinfgen(Q)

    # Build adjacency
    Adj = Q.copy()
    np.fill_diagonal(Adj, 0.0)

    scc, isrec = stronglyconncomp(Adj > 0)
    numSCC = int(np.max(scc))

    # Irreducible case
    if numSCC == 1:
        pi = ctmc_solve(Q)
        return pi, pi.reshape(1, -1), np.array([]), scc, isrec

    # Build SCC index sets (convert 1-based scc to 0-based lists)
    scc_idx = []
    for i in range(1, numSCC + 1):
        scc_idx.append(np.where(scc == i)[0])

    # Classify SCCs
    trans_scc_ids = [i for i in range(numSCC) if not isrec[i]]
    rec_scc_ids = [i for i in range(numSCC) if isrec[i]]

    trans_states = np.sort(np.concatenate([scc_idx[i] for i in trans_scc_ids])) if trans_scc_ids else np.array([], dtype=int)
    rec_states = np.sort(np.concatenate([scc_idx[i] for i in rec_scc_ids])) if rec_scc_ids else np.array([], dtype=int)

    nt = len(trans_states)
    nr = len(rec_states)

    # Extract Q sub-blocks
    Q_tt = Q[np.ix_(trans_states, trans_states)] if nt > 0 else None
    Q_ta = Q[np.ix_(trans_states, rec_states)] if nt > 0 and nr > 0 else None

    # Per-SCC limiting distributions
    pis = np.zeros((numSCC, N))
    pi0_out = np.zeros((numSCC, N))

    for s in range(numSCC):
        p0 = np.zeros(N)
        p0[scc_idx[s]] = 1.0 / len(scc_idx[s])
        pi0_out[s, :] = p0

        hit = np.zeros(nr)
        if nt > 0 and Q_tt is not None and Q_ta is not None:
            p0_t = p0[trans_states]
            if np.any(np.abs(p0_t) > 0):
                try:
                    sojourn = np.linalg.solve(Q_tt.T, -p0_t)
                except np.linalg.LinAlgError:
                    sojourn = np.linalg.lstsq(Q_tt.T, -p0_t, rcond=None)[0]
                hit = sojourn @ Q_ta

        hit = hit + p0[rec_states]

        for c in rec_scc_ids:
            idx_c = scc_idx[c]
            # Find positions of idx_c in rec_states
            loc = np.searchsorted(rec_states, idx_c)
            reachprob = np.sum(hit[loc])
            if reachprob < 1e-15:
                continue

            if len(idx_c) == 1:
                pis[s, idx_c] = reachprob
            else:
                pi_c = ctmc_solve(Q[np.ix_(idx_c, idx_c)])
                pis[s, idx_c] = pi_c * reachprob

    # Compute weighted average
    if pi0 is None:
        pinl = np.ones(numSCC)
        col_sums = np.sum(np.abs(Q), axis=0)
        for j in np.where(col_sums < 1e-12)[0]:
            pinl[int(scc[j]) - 1] = 0
        total = np.sum(pinl)
        if total > 0:
            pinl /= total
        else:
            pinl = np.ones(numSCC) / numSCC
    else:
        pinl = np.zeros(numSCC)
        for i in range(numSCC):
            pinl[i] = np.sum(pi0[scc_idx[i]])

    pi = np.zeros(N)
    for i in range(numSCC):
        if pinl[i] > 0:
            pi += pis[i, :] * pinl[i]

    # Special case: single transient SCC without explicit initial
    if len(trans_scc_ids) == 1 and pi0 is None:
        pi = pis[trans_scc_ids[0], :]

    total = np.sum(pi)
    if total > 0:
        pi /= total

    return pi, pis, pi0_out, scc, isrec


def dtmc_solve_reducible(P: np.ndarray, pin: np.ndarray = None,
                         tol: float = 1e-12) -> Tuple:
    """
    Solve a reducible DTMC using SCC decomposition.

    Decomposes the chain into strongly connected components, computes
    the lumped transition matrix, finds the limiting matrix via spectral
    decomposition, and combines per-SCC solutions.

    Args:
        P: Stochastic transition matrix (possibly reducible)
        pin: Initial distribution vector (optional)
        tol: Numerical tolerance

    Returns:
        Tuple of (pi, pis, pi0, scc, isrec) where:
            pi: Limiting distribution
            pis: Per-SCC limiting distributions
            pi0: Starting distribution per SCC
            scc: SCC label per state (1-based)
            isrec: Boolean array of recurrent SCCs

    MATLAB: matlab/src/api/mc/dtmc_solve_reducible.m
    """
    P = np.asarray(P, dtype=np.float64)
    n = P.shape[0]

    scc, isrec = stronglyconncomp(P)
    numSCC = int(np.max(scc))

    if numSCC == 1:
        pi = dtmc_solve(P)
        return pi, pi.reshape(1, -1), np.array([]), scc, isrec

    # Build SCC index sets
    scc_idx = []
    for i in range(1, numSCC + 1):
        scc_idx.append(np.where(scc == i)[0])

    # Build lumped transition matrix
    Pl = np.zeros((numSCC, numSCC))
    for i in range(numSCC):
        for j in range(numSCC):
            if i != j:
                Pl[i, j] = np.sum(P[np.ix_(scc_idx[i], scc_idx[j])])

    Pl = dtmc_makestochastic(Pl)

    # Ensure recurrent SCCs have self-loops
    for i in range(numSCC):
        if isrec[i]:
            Pl[i, :] = 0.0
            Pl[i, i] = 1.0

    # Compute pinl
    if pin is None:
        pinl = np.ones(numSCC)
        col_sums = np.sum(P, axis=0)
        for j in np.where(col_sums < 1e-12)[0]:
            pinl[int(scc[j]) - 1] = 0
        total = np.sum(pinl)
        if total > 0:
            pinl /= total
        else:
            pinl = np.ones(numSCC) / numSCC
    else:
        pinl = np.zeros(numSCC)
        for i in range(numSCC):
            pinl[i] = np.sum(pin[scc_idx[i]])

    # Compute limiting matrix via power iteration
    PI = _compute_limiting_matrix_power(Pl)

    pi0 = np.zeros((numSCC, numSCC))
    pil = np.zeros((numSCC, numSCC))
    pis = np.zeros((numSCC, n))
    pi = np.zeros(n)

    for i in range(numSCC):
        if pinl[i] > 0:
            pi0[i, :] = 0.0
            pi0[i, i] = 1.0
            pil[i, :] = pi0[i, :] @ PI

            for j in range(numSCC):
                states_j = scc_idx[j]
                if len(states_j) > 0:
                    pi_j = dtmc_solve(P[np.ix_(states_j, states_j)])
                    pis[i, states_j] = pil[i, j] * pi_j

            pi += pis[i, :] * pinl[i]

    # Special case: single transient SCC
    transient_sccs = [i for i in range(numSCC) if not isrec[i]]
    if len(transient_sccs) == 1 and pin is None:
        pi = pis[transient_sccs[0], :]

    return pi, pis, pi0, scc, isrec


def _compute_limiting_matrix_power(P: np.ndarray, max_iter: int = 1000,
                                    tol: float = 1e-10) -> np.ndarray:
    """
    Compute limiting matrix P^inf using power iteration.

    For an ergodic DTMC, P^n converges to a matrix where each row
    is the stationary distribution.

    Args:
        P: Stochastic transition matrix
        max_iter: Maximum number of iterations
        tol: Convergence tolerance

    Returns:
        Limiting matrix
    """
    Pk = P.copy()
    for _ in range(max_iter):
        Pk1 = Pk @ P
        max_diff = np.max(np.abs(Pk1 - Pk))
        Pk = Pk1
        if max_diff < tol:
            break
    return Pk


__all__ = [
    # CTMC functions
    'ctmc_makeinfgen',
    'ctmc_solve',
    'ctmc_solve_full',
    'ctmc_rand',
    'ctmc_timereverse',
    'ctmc_randomization',
    'ctmc_uniformization',
    'ctmc_transient',
    'ctmc_relsolve',
    'ctmc_simulate',
    'ctmc_stochcomp',
    'ctmc_courtois',
    'ctmc_kms',
    'ctmc_takahashi',
    'ctmc_multi',
    'ctmc_solve_reducible',
    'ctmc_solve_reducible_blkdecomp',
    'weaklyconncomp',
    'stronglyconncomp',
    # DTMC functions
    'dtmc_makestochastic',
    'dtmc_isfeasible',
    'dtmc_solve',
    'dtmc_rand',
    'dtmc_simulate',
    'dtmc_stochcomp',
    'dtmc_timereverse',
    'dtmc_uniformization',
    'dtmc_solve_reducible',
    # Data classes
    'ConnectedComponents',
    'CTMCSolveResult',
]
