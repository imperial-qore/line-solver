"""
Load-dependent MVA for Product-Form Queueing Networks.

Implements exact MVA for closed queueing networks with load-dependent
service rates, where the service rate at each station may depend on
the number of jobs present.
"""

import numpy as np
from typing import Tuple, Optional


def _pprod_init(N: np.ndarray) -> np.ndarray:
    """Initialize population product iterator.

    Args:
        N: Maximum population vector.

    Returns:
        Initial state (zeros).
    """
    return np.zeros_like(N, dtype=int)


def _pprod_next(n: np.ndarray, N: np.ndarray) -> Optional[np.ndarray]:
    """Get next population vector in lexicographic order.

    Iterates through all vectors n where 0 <= n <= N.

    Args:
        n: Current population vector.
        N: Maximum population vector.

    Returns:
        Next population vector, or None if done.
    """
    n = n.copy()
    R = len(N)

    # Check if at maximum
    if np.all(n == N):
        return None

    # Find rightmost index that can be incremented
    s = R - 1
    while s >= 0 and n[s] == N[s]:
        n[s] = 0
        s -= 1

    if s < 0:
        return None

    n[s] += 1
    return n


def _hashpop(n: np.ndarray, N: np.ndarray) -> int:
    """Hash population vector to linear index.

    Converts a multi-dimensional population vector to a unique
    linear index for array storage.

    Args:
        n: Population vector.
        N: Maximum population vector.

    Returns:
        Linear index (0-based).
    """
    index = 0
    for r in range(len(N)):
        prod = 1
        for j in range(r):
            prod *= N[j] + 1
        index += prod * n[r]
    return index


def _oner(n: np.ndarray, s: int) -> np.ndarray:
    """Create population vector with one less job in class s.

    Args:
        n: Population vector.
        s: Class index to decrement.

    Returns:
        Modified population vector.
    """
    n_copy = n.copy()
    if s >= 0 and s < len(n_copy):
        n_copy[s] = max(0, n_copy[s] - 1)
    return n_copy


def pfqn_mvald(
    L: np.ndarray,
    N: np.ndarray,
    Z: np.ndarray,
    mu: np.ndarray,
    stabilize: bool = True
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, bool, np.ndarray]:
    """
    Exact MVA for load-dependent closed queueing networks.

    This algorithm extends standard MVA to handle stations where
    the service rate depends on the number of jobs present.

    Args:
        L: Service demand matrix (M x R) - rows are stations, columns are classes.
        N: Population vector (R,).
        Z: Think time vector (R,).
        mu: Load-dependent rate matrix (M x Ntot) where mu[i,k] is the
            service rate at station i when k jobs are present.
        stabilize: Force non-negative probabilities (default: True).

    Returns:
        Tuple of (XN, QN, UN, CN, lGN, isNumStable, pi):
            XN: Class throughputs (1 x R).
            QN: Mean queue lengths (M x R).
            UN: Utilizations (M x R).
            CN: Cycle times (1 x R).
            lGN: Log normalizing constant evolution.
            isNumStable: True if numerically stable.
            pi: Marginal queue-length probabilities (M x Ntot+1).
    """
    L = np.atleast_2d(np.asarray(L, dtype=float))
    N = np.asarray(N, dtype=int).ravel()
    Z = np.asarray(Z, dtype=float).ravel()
    mu = np.atleast_2d(np.asarray(mu, dtype=float))

    M, R = L.shape
    Ntot = int(np.sum(N))

    # Initialize outputs
    XN = np.zeros((1, R))
    QN = np.zeros((M, R))
    UN = np.zeros((M, R))
    CN = np.zeros((1, R))
    WN = np.zeros((M, R))

    # Handle negative population
    if np.any(N < 0):
        lGN = np.array([-np.inf])
        pi = np.zeros((M, Ntot + 1))
        return XN, QN, UN, CN, lGN, True, pi

    # Number of population states
    num_states = int(np.prod(N + 1))

    # Throughputs for each population
    Xs = np.zeros((R, num_states))

    # Marginal queue-length probabilities: pi[i,k,pop] = P(queue i has k jobs | pop)
    pi_all = np.ones((M, Ntot + 1, num_states))

    is_num_stable = True
    warn = True
    lGN_list = [0.0]

    # Iterate through all population vectors
    n = _pprod_init(N)
    while n is not None:
        pop_idx = _hashpop(n, N)
        WN_local = np.zeros((M, R))

        # Compute waiting times for each class
        for s in range(R):
            if n[s] > 0:
                n_minus_s = _oner(n, s)
                n_minus_s_idx = _hashpop(n_minus_s, N)

                for ist in range(M):
                    for k in range(1, int(np.sum(n)) + 1):
                        if k <= mu.shape[1]:
                            mu_val = mu[ist, k - 1]  # mu is 0-indexed, k starts at 1
                            if mu_val > 0:
                                WN_local[ist, s] += (L[ist, s] / mu_val) * k * pi_all[ist, k - 1, n_minus_s_idx]

                # Compute throughput
                denom = Z[s] + np.sum(WN_local[:, s])
                if denom > 0:
                    Xs[s, pop_idx] = n[s] / denom

        # Compute pi(k|n) for k >= 1
        for k in range(1, int(np.sum(n)) + 1):
            for ist in range(M):
                pi_all[ist, k, pop_idx] = 0

            for s in range(R):
                if n[s] > 0:
                    n_minus_s = _oner(n, s)
                    n_minus_s_idx = _hashpop(n_minus_s, N)

                    for ist in range(M):
                        if k <= mu.shape[1]:
                            mu_val = mu[ist, k - 1]
                            if mu_val > 0:
                                pi_all[ist, k, pop_idx] += (L[ist, s] / mu_val) * Xs[s, pop_idx] * pi_all[ist, k - 1, n_minus_s_idx]

        # Compute pi(0|n)
        for ist in range(M):
            sum_pi = np.sum(pi_all[ist, 1:int(np.sum(n)) + 1, pop_idx])
            p0 = 1.0 - sum_pi

            if p0 < 0:
                if warn:
                    warn = False
                    is_num_stable = False
                if stabilize:
                    pi_all[ist, 0, pop_idx] = np.finfo(float).eps
                else:
                    pi_all[ist, 0, pop_idx] = p0
            else:
                pi_all[ist, 0, pop_idx] = p0

        # Track log normalizing constant evolution
        last_nnz = -1
        for r in range(R - 1, -1, -1):
            if n[r] > 0:
                last_nnz = r
                break

        if last_nnz >= 0:
            # Check condition: sum(n[0:last_nnz-1]) == sum(N[0:last_nnz-1]) and sum(n[last_nnz+1:]) == 0
            cond1 = np.sum(n[:last_nnz]) == np.sum(N[:last_nnz]) if last_nnz > 0 else True
            cond2 = np.sum(n[last_nnz + 1:]) == 0 if last_nnz < R - 1 else True
            if cond1 and cond2:
                X_val = Xs[last_nnz, pop_idx]
                if X_val > 0:
                    lGN_list.append(lGN_list[-1] - np.log(X_val))

        # Next population
        n = _pprod_next(n, N)

    # Get final results for full population N
    N_idx = _hashpop(N, N)
    XN = Xs[:, N_idx].reshape(1, -1)
    pi = pi_all[:, :, N_idx]

    # Compute waiting times at full population
    WN = np.zeros((M, R))
    for s in range(R):
        if N[s] > 0:
            n_minus_s = _oner(N, s)
            n_minus_s_idx = _hashpop(n_minus_s, N)

            for ist in range(M):
                for k in range(1, Ntot + 1):
                    if k <= mu.shape[1]:
                        mu_val = mu[ist, k - 1]
                        if mu_val > 0:
                            WN[ist, s] += (L[ist, s] / mu_val) * k * pi_all[ist, k - 1, n_minus_s_idx]

    # Compute metrics
    QN = WN * np.tile(XN, (M, 1))
    UN = 1 - pi[:, 0].reshape(-1, 1)  # Utilization from pi(0)
    UN = np.tile(UN, (1, R))  # Broadcast to all classes
    # Zero out utilization for classes with zero demand (Disabled services)
    UN[L == 0] = 0.0

    # Compute cycle times
    CN = np.zeros((1, R))
    for r in range(R):
        if XN[0, r] > 0:
            CN[0, r] = N[r] / XN[0, r] - Z[r]

    lGN = np.array(lGN_list)

    return XN, QN, UN, CN, lGN, is_num_stable, pi


def pfqn_mvams(
    lambda_arr: np.ndarray,
    L: np.ndarray,
    N: np.ndarray,
    Z: np.ndarray,
    mi: Optional[np.ndarray] = None,
    S: Optional[np.ndarray] = None
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float]:
    """
    General-purpose MVA for mixed networks with multiserver nodes.

    This function handles networks with open/closed classes and multi-server
    stations, routing to the appropriate specialized algorithm.

    Args:
        lambda_arr: Arrival rate vector (R,). Use 0 for closed classes.
        L: Service demand matrix (M x R).
        N: Population vector (R,). Use np.inf for open classes.
        Z: Think time vector (R,).
        mi: Queue replication factors (M,) (default: ones).
        S: Number of servers per station (M,) (default: ones).

    Returns:
        Tuple of (XN, QN, UN, CN, lG):
            XN: Class throughputs (1 x R).
            QN: Mean queue lengths (M x R).
            UN: Utilizations (M x R).
            CN: Cycle times (1 x R).
            lG: Log normalizing constant.

    References:
        Original MATLAB: matlab/src/api/pfqn/pfqn_mvams.m
    """
    from .mva import pfqn_mva
    from .mixed import pfqn_mvamx
    from .mvaldmx import pfqn_mvaldms

    L = np.atleast_2d(np.asarray(L, dtype=float))
    N = np.asarray(N, dtype=float).ravel()  # Keep as float to allow inf
    Z = np.asarray(Z, dtype=float).ravel()
    lambda_arr = np.asarray(lambda_arr, dtype=float).ravel()

    M, R = L.shape

    # Compute Ntot only from finite (closed class) populations
    Ntot = int(np.sum(N[np.isfinite(N)]))

    if S is None:
        S = np.ones(M, dtype=int)
    else:
        S = np.asarray(S, dtype=int).ravel()

    if mi is None:
        mi = np.ones(M)
    else:
        mi = np.asarray(mi, dtype=float).ravel()

    # If Z is empty, use zeros
    if len(Z) == 0:
        Z = np.zeros(R)

    # Build load-dependent rate matrix for multi-server stations
    if Ntot > 0:
        mu = np.ones((M, Ntot))
        for ist in range(M):
            for k in range(1, Ntot + 1):
                mu[ist, k - 1] = min(k, S[ist])
    else:
        mu = np.ones((M, 1))

    # Get max servers (only from finite S values)
    S_finite = S[np.isfinite(S)]
    max_S = int(np.max(S_finite)) if len(S_finite) > 0 else 1

    if max_S == 1:  # No multi-server nodes
        if np.any(np.isinf(N)):  # Open or mixed model
            XN, QN, UN, CN, lG = pfqn_mvamx(lambda_arr, L, N, Z, mi)
        else:  # Closed model
            N_int = N.astype(int)
            XN, CN_out, QN, UN, RN, TN, AN = pfqn_mva(L, N_int, Z, mi)
            # Compute log normalizing constant
            from .nc import pfqn_ca
            _, lG = pfqn_ca(L, N_int, Z)
            CN = CN_out  # Use cycle times from MVA
    else:  # Model has multi-server nodes
        if np.any(np.isinf(N)):  # Open or mixed model
            if np.max(mi) == 1:
                lG = np.nan  # NC not available in this case
                XN, QN, UN, CN, _ = pfqn_mvaldms(lambda_arr, L, N, Z, S)
            else:
                raise ValueError("Queue replicas not available in exact MVA for mixed models.")
        else:  # Closed model with multi-server
            N_int = N.astype(int)
            XN, QN, UN, CN, lGN, _, _ = pfqn_mvald(L, N_int, Z, mu)
            lG = lGN[-1] if len(lGN) > 0 else np.nan

    return XN, QN, UN, CN, lG


__all__ = [
    'pfqn_mvald',
    'pfqn_mvams',
]
