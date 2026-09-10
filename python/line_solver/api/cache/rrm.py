"""
Random Replacement Model (RRM) Mean-Field Cache Analysis.

Native Python implementations of mean-field methods for analyzing cache
systems with random replacement policies.

Key functions:
    cache_rrm_meanfield_ode: ODE function for RRM dynamics
    cache_rrm_meanfield: Solve RRM mean-field steady state

References:
    Original MATLAB: matlab/src/api/cache/cache_rrm_meanfield*.m
"""

import numpy as np
from scipy.integrate import solve_ivp
from typing import Tuple, Optional


def cache_rrm_meanfield_ode(t: float, x: np.ndarray,
                            lambd: np.ndarray, m: np.ndarray,
                            n: int, h: int) -> np.ndarray:
    """
    ODE function for RRM mean-field cache dynamics.

    Defines the differential equations for the Random Replacement Model
    mean-field cache dynamics.

    Args:
        t: Time variable (unused, for ODE solver compatibility)
        x: State vector of length n*(h+1), representing probabilities
        lambd: Arrival rates per item (n,)
        m: Cache capacity vector (h,)
        n: Number of items
        h: Number of cache levels

    Returns:
        Time derivative of state vector

    References:
        Original MATLAB: matlab/src/api/cache/cache_rrm_meanfield_ode.m
    """
    lambd = np.asarray(lambd, dtype=np.float64).ravel()
    m = np.asarray(m, dtype=np.float64).ravel()

    # Reshape state vector to matrix form: x[k, s]
    x = x.reshape((n, 1 + h))

    dxdt = np.zeros((n, 1 + h))

    for k in range(n):
        for s in range(1, h + 1):  # s = 1, ..., h (1-indexed in formulation)
            # First term: promotion from list s-1
            sum1 = 0.0
            for k1 in range(n):
                sum1 += lambd[k1] / m[s-1] * x[k1, s-1] * x[k, s]

            # Second term: demotion from list s+1
            if s < h:
                sum2 = 0.0
                for k1 in range(n):
                    sum2 += lambd[k1] / m[s] * x[k1, s] * x[k, s+1]
                sum2 = sum2 - lambd[k] * x[k, s]
            else:
                sum2 = 0.0

            # Drift component
            dxdt[k, s] = lambd[k] * x[k, s-1] - sum1 + sum2

        # Case s=0: conservation of probability
        dxdt[k, 0] = -np.sum(dxdt[k, 1:h+1])

    return dxdt.flatten()


def cache_rrm_meanfield(lambd: np.ndarray, m: np.ndarray,
                        t_end: float = 10000.0,
                        seed: int = 23000
                        ) -> Tuple[np.ndarray, float, float]:
    """
    Solve RRM mean-field steady state using ODE integration.

    Computes the steady-state probability distribution for a cache with
    random replacement policy using mean-field ODE dynamics.

    Args:
        lambd: Arrival rates per item (n,)
        m: Cache capacity vector (h,)
        t_end: End time for ODE integration (default: 10000.0)
        seed: Random seed for initial conditions (default: 23000)

    Returns:
        Tuple of (prob, missrate, missratio) where:
            - prob: Steady-state probability matrix (n x h+1)
            - missrate: Global miss rate (lambda * miss_prob)
            - missratio: Miss ratio (missrate / sum(lambda))

    References:
        Original MATLAB: matlab/src/api/cache/cache_rrm_meanfield.m
    """
    np.random.seed(seed)

    lambd = np.asarray(lambd, dtype=np.float64).ravel()
    m = np.asarray(m, dtype=np.float64).ravel()

    n = len(lambd)
    h = len(m)

    # Initial condition: all items start in state 0 (miss)
    x0 = np.zeros((n, 1 + h))
    x0[:, 0] = 1.0

    # Solve ODE
    def ode_func(t, x):
        return cache_rrm_meanfield_ode(t, x, lambd, m, n, h)

    sol = solve_ivp(ode_func, [0, t_end], x0.flatten(),
                    method='BDF',  # Stiff solver like ode23s
                    rtol=1e-6, atol=1e-9)

    # Extract final state
    x_final = sol.y[:, -1].reshape((n, 1 + h))

    # Compute miss metrics
    missrate = np.dot(lambd, x_final[:, 0])
    lambda_sum = np.sum(lambd)
    missratio = missrate / lambda_sum if lambda_sum > 0 else 0.0

    return x_final, missrate, missratio


def cache_gamma_lp(lambd: np.ndarray, R: list) -> Tuple[np.ndarray, int, int, int]:
    """
    Compute gamma parameters for cache models using linear programming approach.

    Computes item popularity probabilities at each cache level based on
    arrival rates and routing probabilities.

    Args:
        lambd: Arrival rates per user per item per list (u x n x h+1)
        R: Routing probability structure (list of lists, R[v][i] is matrix for user v, item i)

    Returns:
        Tuple of (gamma, u, n, h, parent) where:
            - gamma: Item popularity probabilities at each level (n x h)
            - u: Number of users
            - n: Number of items
            - h: Number of cache levels
            - parent: Parent list of each list, 0-based with -1 for the lists
              rooted in the miss list (h,)

    References:
        Original MATLAB: matlab/src/api/cache/cache_gamma_lp.m
    """
    lambd = np.asarray(lambd, dtype=np.float64)

    u = lambd.shape[0]  # number of users
    n = lambd.shape[1]  # number of items
    h = lambd.shape[2] - 1  # number of lists

    gamma = np.zeros((n, h))

    def find_parent(Rvi, j):
        """Find parent of node j in routing matrix."""
        if j == 0:
            return None
        parents = np.where(Rvi[:j, j] > 0)[0]
        if len(parents) == 0:
            return None
        if len(parents) > 1:
            raise ValueError("Cache has a list with more than one parent, but structure must be a tree.")
        return parents[0]

    for i in range(n):
        for j in range(h):
            # Compute gamma(i, j)
            # Sum routing matrices across users
            Rvi = np.zeros_like(R[0][i])
            for v in range(u):
                Rvi = Rvi + R[v][i]

            # Build path from root to level j+1 (0-indexed level is j, but +1 for 1-indexed list)
            target = j + 1  # 1-indexed list
            Pij = [target]

            # Trace back to root
            pr_j = find_parent(Rvi, target)
            while pr_j is not None:
                Pij.insert(0, pr_j)
                pr_j = find_parent(Rvi, pr_j)

            if len(Pij) < 2:
                gamma[i, j] = 0.0
            else:
                gamma[i, j] = 1.0
                for li in range(1, len(Pij)):
                    l_1 = Pij[li - 1]
                    l = Pij[li]
                    y = 0.0
                    for v in range(u):
                        # In Python, levels are 0-indexed: 0=miss, 1..h=cache levels
                        # MATLAB uses 1-indexed: 1=miss, 2..h+1=cache levels
                        # For a path segment from l_1 to l:
                        # - l_1 is the source level (0-indexed in Python)
                        # - We need to sum over all arrival rates up to and including l_1
                        # This means range(l_1 + 1) which gives [0, 1, ..., l_1]
                        for t in range(l_1 + 1):
                            y += lambd[v, i, t] * R[v][i][t, l]
                    gamma[i, j] *= y

    # tree structure of the lists, read off the routing matrix of item 0
    # aggregated over users -- the same matrix the gamma loop walks, so a
    # per-item tree can never be contaminated by another item's matrix
    Rtot = np.zeros_like(R[0][0])
    for v in range(u):
        Rtot = Rtot + R[v][0]
    parent = np.full(h, -1, dtype=int)
    for j in range(h):
        pj = find_parent(Rtot, j + 1)
        parent[j] = -1 if pj is None else int(pj) - 1

    return gamma, u, n, h, parent


def _bfs_path(A: np.ndarray, src: int, dst: int):
    """Breadth-first shortest path over the strictly positive entries of A."""
    nn = A.shape[0]
    if src >= nn or dst >= nn or src < 0 or dst < 0:
        return []
    visited = [False] * nn
    parent = [-1] * nn
    queue = [src]
    visited[src] = True
    while queue:
        current = queue.pop(0)
        if current == dst:
            path = []
            node = dst
            while node != -1:
                path.insert(0, node)
                node = parent[node]
            return path
        for nxt in range(nn):
            if not visited[nxt] and A[current, nxt] > 0:
                visited[nxt] = True
                parent[nxt] = current
                queue.append(nxt)
    return []


def cache_gamma(lambd: np.ndarray, R: list) -> Tuple[np.ndarray, int, int, int]:
    """Access factors of a multi-list cache whose lists form a general graph.

    Companion of `cache_gamma_lp`, which requires the access structure to be a
    tree and walks the unique parent relation. Here the structure is only
    required to be reachable: the path to list j is the BREADTH-FIRST shortest
    path in the access graph of item i, so a list with several parents is
    admissible and the first shortest path found in node order is the one
    taken. Along that path,

        gamma[i,j] = (sum_v lambd[v,i,0]) prod_edges (a,b) sum_v lambd[v,i,a] R[v][i][a,b]

    THREE DIVERGENCES FROM cache_gamma_lp, all of which change the number, so
    the two are not substitutes: the destination of column j is node j and not
    node j+1, so column 0 carries no edge factor at all; the leading factor is
    the aggregate miss-node request rate rather than one; and each edge factor
    reads lambd at the SOURCE node a alone rather than summing over every
    t <= a. Use cache_gamma_lp for the access factors a cache solver consumes.

    The adjacency is read from user 0 only, so a model whose users route an
    item differently is analysed on the first user's graph. An unreachable list
    gives gamma[i,j] = 0.

    Args:
        lambd: (u, n, h+1) request rate of user v for item i while at node t.
        R: (u, n) list of (h+1, h+1) routing matrices.

    Returns:
        Tuple of (gamma, u, n, h) with gamma of shape (n, h).

    References:
        jar/src/main/java/jline/api/cache/Cache_gamma.java
    """
    lambd = np.asarray(lambd, dtype=np.float64)
    u = lambd.shape[0]
    n = lambd.shape[1]
    h = lambd.shape[2] - 1

    gamma = np.zeros((n, h))
    for i in range(n):
        graph = np.asarray(R[0][i], dtype=np.float64)
        for j in range(h):
            path = _bfs_path(graph, 0, j)
            if not path:
                continue
            g = 0.0
            for v in range(u):
                g += lambd[v, i, 0]
            for li in range(1, len(path)):
                a = path[li - 1]
                b = path[li]
                y = 0.0
                for v in range(u):
                    y += lambd[v, i, a] * R[v][i][a, b]
                g *= y
            gamma[i, j] = g
    return gamma, u, n, h


__all__ = [
    'cache_rrm_meanfield_ode',
    'cache_rrm_meanfield',
    'cache_gamma_lp',
    'cache_gamma',
]
