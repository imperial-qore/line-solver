"""
Loss Network Analysis via Erlang Fixed Point.

Native Python implementations of fixed-point algorithms for analyzing
loss networks using Erlang formulas. Provides efficient computation
of blocking probabilities for circuit-switched networks.

References:
    Kelly, F.P., "Blocking probabilities in large circuit-switched networks."
    Advances in Applied Probability, 1986.
"""

import numpy as np
from scipy import special
from typing import Tuple


def _erlang_b(nu: float, C: int) -> float:
    """
    Calculate the Erlang B formula for blocking probability.

    The Erlang B formula computes the blocking probability in a loss system
    where calls arrive according to a Poisson process.

    Args:
        nu: Traffic intensity (offered load = arrival rate * service time).
        C: Number of servers (capacity).

    Returns:
        Blocking probability.
    """
    if C < 0:
        return 1.0
    if nu <= 0:
        return 0.0

    # Use log-space computation for numerical stability
    # log(nu^C / C!) - log(sum_{i=0}^C nu^i / i!)
    log_terms = np.array([i * np.log(nu) - special.gammaln(i + 1) for i in range(C + 1)])
    log_max = np.max(log_terms)
    log_sum = log_max + np.log(np.sum(np.exp(log_terms - log_max)))

    log_block = C * np.log(nu) - special.gammaln(C + 1) - log_sum
    return np.exp(log_block)


def lossn_erlangfp(nu: np.ndarray, A: np.ndarray, c: np.ndarray,
                   tol: float = 1e-8, max_iter: int = 1000
                   ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """
    Erlang fixed point approximation for loss networks.

    Calls (jobs) on route (class) r arrive according to Poisson rate nu_r.
    Call service times on route r have unit mean.

    The link capacity requirements are:
        sum_r A[j,r] * n[j,r] < c[j]
    for all links j, where n[j,r] counts calls on route r on link j.

    Args:
        nu: Arrival rates vector (R,) for each route.
        A: Capacity requirement matrix (J, R) - A[j,r] is capacity
           required on link j by route r.
        c: Capacity vector (J,) - c[j] is capacity of link j.
        tol: Convergence tolerance.
        max_iter: Maximum iterations.

    Returns:
        Tuple of (qlen, loss, eblock, niter) where:
            - qlen: Mean queue-length for each route (R,)
            - loss: Loss probability for each route (R,)
            - eblock: Blocking probability for each link (J,)
            - niter: Number of iterations

    Example:
        >>> nu = np.array([0.3, 0.1])
        >>> A = np.array([[1, 1], [1, 4]])  # 2 links, 2 routes
        >>> c = np.array([1, 3])
        >>> qlen, loss, eblock, niter = lossn_erlangfp(nu, A, c)
    """
    nu = np.asarray(nu, dtype=np.float64).ravel()
    A = np.asarray(A, dtype=np.float64)
    c = np.asarray(c, dtype=np.float64).ravel()

    R = len(nu)  # Number of routes
    J = len(c)   # Number of links

    # Ensure A has correct shape
    if A.shape != (J, R):
        raise ValueError(f"A must have shape ({J}, {R}), got {A.shape}")

    # Initialize blocking probabilities
    E = np.full(J, 0.5)
    E_old = np.zeros(J)
    niter = 0

    # Fixed point iteration
    while np.linalg.norm(E - E_old) > tol and niter < max_iter:
        niter += 1
        E_old = E.copy()

        for j in range(J):
            rho_j = 0.0

            for r in range(R):
                if A[j, r] > 0:
                    term = nu[r] * A[j, r]

                    # Multiply by (1 - E_i)^A[i,r] for all links i
                    for i in range(J):
                        if A[i, r] > 1e-10:
                            term *= (1 - E_old[i]) ** A[i, r]

                    rho_j += term

            # Adjust for self-blocking
            if E_old[j] < 1:
                rho_j /= (1 - E_old[j])

            E[j] = _erlang_b(rho_j, int(c[j]))

    # Compute queue lengths: effective arrival rate * (1 - blocking)
    qlen = nu.copy()
    for r in range(R):
        for j in range(J):
            qlen[r] *= (1 - E[j]) ** A[j, r]

    qlen = np.maximum(qlen, 0.0)

    # Compute loss probabilities
    loss = np.zeros(R)
    for r in range(R):
        if nu[r] > 0:
            loss[r] = 1 - qlen[r] / nu[r]
        else:
            loss[r] = 0.0

    return qlen, loss, E, niter


def erlang_b(offered_load: float, servers: int) -> float:
    """
    Compute Erlang B blocking probability.

    The Erlang B formula gives the probability that an arriving call
    is blocked in an M/M/c/c loss system.

    Args:
        offered_load: Traffic intensity (arrival rate * service time).
        servers: Number of servers.

    Returns:
        Blocking probability.

    Example:
        >>> erlang_b(10.0, 12)  # Offered load 10 Erlang, 12 channels
        0.1054...
    """
    return _erlang_b(offered_load, servers)


def erlang_c(offered_load: float, servers: int) -> float:
    """
    Compute Erlang C delay probability.

    The Erlang C formula gives the probability that an arriving call
    must wait in an M/M/c queue.

    Args:
        offered_load: Traffic intensity (arrival rate * service time).
        servers: Number of servers.

    Returns:
        Delay probability (probability of waiting).

    Example:
        >>> erlang_c(10.0, 12)  # Offered load 10 Erlang, 12 agents
    """
    if servers < 1:
        return 1.0
    if offered_load <= 0:
        return 0.0

    rho = offered_load / servers
    if rho >= 1:
        return 1.0

    # Erlang C = c * B(c, A) / (c - A * (1 - B(c, A)))
    B = _erlang_b(offered_load, servers)
    C = servers * B / (servers - offered_load * (1 - B))

    return min(1.0, max(0.0, C))


__all__ = [
    'lossn_erlangfp',
    'erlang_b',
    'erlang_c',
]
