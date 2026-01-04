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
        Tuple of (gamma, u, n, h) where:
            - gamma: Item popularity probabilities at each level (n x h)
            - u: Number of users
            - n: Number of items
            - h: Number of cache levels

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
                        for t in range(l_1 + 1):  # t from 0 to l_1 (1-indexed: 1 to l_1)
                            y += lambd[v, i, t] * R[v][i][t, l]
                    gamma[i, j] *= y

    return gamma, u, n, h


__all__ = [
    'cache_rrm_meanfield_ode',
    'cache_rrm_meanfield',
    'cache_gamma_lp',
]
