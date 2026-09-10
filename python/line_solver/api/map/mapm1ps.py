"""
MAP/M/1-PS Sojourn Time Distribution.

Computes the complementary distribution function of sojourn time in a
MAP/M/1 processor-sharing queue.

References:
    Masuyama, H., & Takine, T. (2003). Sojourn time distribution in a
    MAP/M/1 processor-sharing queue. Operations Research Letters, 31(6), 406-412.
"""

import numpy as np
from scipy import linalg
from typing import Union

ArrayLike = Union[np.ndarray, list]


def _compute_stationary_distribution(C: np.ndarray, D: np.ndarray) -> np.ndarray:
    """
    Compute stationary distribution of MAP.

    Solves pi * Q = 0 with pi * e = 1.

    Args:
        C: MAP C matrix (transitions without arrivals).
        D: MAP D matrix (transitions with arrivals).

    Returns:
        Stationary distribution as row vector.
    """
    M = C.shape[0]
    Q = C + D

    # Solve pi * Q = 0 with pi * e = 1
    # Convert to Q^T * pi^T = 0, replace last equation with sum = 1
    A = Q.T.copy()

    # Replace last row with normalization constraint
    A[-1, :] = 1.0

    b = np.zeros(M)
    b[-1] = 1.0

    # Solve system
    pi_T = linalg.solve(A, b)
    return pi_T.reshape(1, -1)


def _compute_r_matrix(C: np.ndarray, D: np.ndarray, mu: float) -> np.ndarray:
    """
    Compute R matrix (minimal nonnegative solution).

    Solves: D + R(C - mu*I) + mu*R^2 = 0

    Args:
        C: MAP C matrix.
        D: MAP D matrix.
        mu: Service rate.

    Returns:
        R matrix.
    """
    M = C.shape[0]

    if M == 1:
        # Scalar case: use quadratic formula
        # mu*R^2 + (C - mu)*R + D = 0
        a = mu
        b = C[0, 0] - mu
        c = D[0, 0]

        discriminant = b * b - 4 * a * c
        if discriminant < 0:
            raise ValueError("No real solution for R matrix")

        R1 = (-b - np.sqrt(discriminant)) / (2 * a)
        R2 = (-b + np.sqrt(discriminant)) / (2 * a)

        # Choose minimal nonnegative solution
        if R1 >= 0 and R1 < 1:
            Rval = R1
        elif R2 >= 0 and R2 < 1:
            Rval = R2
        else:
            raise ValueError("No valid R in [0, 1)")

        return np.array([[Rval]])

    # Matrix case: Newton iteration
    I = np.eye(M)
    CminusMuI = C - mu * I

    # Initial guess: R = -D * inv(C - mu*I)
    R = -D @ linalg.inv(CminusMuI)

    max_iter = 1000
    tol = 1e-12

    for _ in range(max_iter):
        Rold = R.copy()

        # F(R) = D + R(C - mu*I) + mu*R^2
        F = D + R @ CminusMuI + mu * R @ R

        # F'(R) = (C - mu*I) + 2*mu*R
        Fprime = CminusMuI + 2 * mu * R

        # Newton step: R_new = R_old - F(R) * inv(F'(R))
        R = R - F @ linalg.inv(Fprime)

        # Check convergence
        if np.max(np.abs(R - Rold)) < tol:
            break

    # Ensure non-negativity
    R = np.maximum(R, 0.0)

    return R


def _determine_n_epsilon(
    pi0: np.ndarray,
    R: np.ndarray,
    D: np.ndarray,
    lam: float,
    epsilon: float
) -> int:
    """
    Determine N(epsilon) - minimum n such that truncated sum > 1 - epsilon.
    """
    M = R.shape[0]
    e = np.ones((M, 1))
    cumsum_prob = 0.0
    Rpower = np.eye(M)

    for n in range(1001):
        cumsum_prob += (1.0 / lam) * (pi0 @ Rpower @ D @ e).item()
        if cumsum_prob > 1 - epsilon:
            return n
        Rpower = Rpower @ R

    return 100  # Default fallback


def _poisson_pmf(k: int, lam: float) -> float:
    """
    Compute Poisson PMF using log-space computation to avoid overflow.
    """
    if lam == 0.0:
        return 1.0 if k == 0 else 0.0

    # P(k; lambda) = exp(k*log(lambda) - lambda - log(k!))
    log_prob = k * np.log(lam) - lam

    # Subtract log(k!)
    for i in range(1, k + 1):
        log_prob -= np.log(float(i))

    return np.exp(log_prob)


def _compute_h_recursive(
    C: np.ndarray,
    D: np.ndarray,
    mu: float,
    N: int,
    K: int,
    theta: float
) -> list:
    """
    Compute h_{n,k} recursively using Algorithm from Theorem 1.

    h_{n,0} = e for all n
    h_{n,k+1} = 1/(theta+mu) * [n*mu/(n+1) h_{n-1,k} + (theta*I + C)h_{n,k} + D h_{n+1,k}]

    Returns:
        2D list h[n][k] of (M,1) arrays.
    """
    M = C.shape[0]
    I = np.eye(M)
    e = np.ones((M, 1))

    theta_plus_mu = theta + mu
    thetaI_plus_C = theta * I + C

    # Initialize h array: h[n][k] where n in [0,N], k in [0,K]
    h = [[np.zeros((M, 1)) for _ in range(K + 1)] for _ in range(N + 1)]

    # Base case: h_{n,0} = e for all n
    for n in range(N + 1):
        h[n][0] = e.copy()

    # Recursive computation: fill column by column (increasing k)
    for k in range(K):
        for n in range(N + 1):
            term1 = np.zeros((M, 1))
            term2 = thetaI_plus_C @ h[n][k]
            term3 = np.zeros((M, 1))

            # First term: n*mu/(n+1) * h_{n-1,k}
            if n > 0:
                term1 = (n * mu / (n + 1)) * h[n - 1][k]

            # Third term: D * h_{n+1,k}
            if n < N:
                term3 = D @ h[n + 1][k]

            h[n][k + 1] = (term1 + term2 + term3) / theta_plus_mu

    return h


def map_m1ps_cdf_respt(
    C: ArrayLike,
    D: ArrayLike,
    mu: float,
    x: ArrayLike,
    epsilon: float = 1e-11,
    epsilon_prime: float = 1e-10
) -> np.ndarray:
    """
    Compute complementary sojourn time CDF for MAP/M/1-PS queue.

    Based on:
        Masuyama, H., & Takine, T. (2003). Sojourn time distribution in a
        MAP/M/1 processor-sharing queue.

    Args:
        C: MAP C matrix (transitions without arrivals).
        D: MAP D matrix (transitions with arrivals).
        mu: Service rate (must be > 0).
        x: Time points for CDF evaluation.
        epsilon: Queue length truncation parameter.
        epsilon_prime: Uniformization truncation parameter.

    Returns:
        Complementary CDF values W_bar(x) = Pr[W > x] at each point in x.

    Example:
        >>> # Exponential arrivals (Poisson process with rate lambda=0.5)
        >>> C = np.array([[-0.5]])
        >>> D = np.array([[0.5]])
        >>> mu = 1.0
        >>> x = np.array([0.0, 0.5, 1.0, 2.0, 5.0])
        >>> cdf = map_m1ps_cdf_respt(C, D, mu, x)
    """
    C = np.asarray(C, dtype=np.float64)
    D = np.asarray(D, dtype=np.float64)
    x = np.asarray(x, dtype=np.float64).ravel()

    M = C.shape[0]
    if C.shape != (M, M):
        raise ValueError("C must be square")
    if D.shape != (M, M):
        raise ValueError("C and D must have same dimensions")
    if mu <= 0:
        raise ValueError("Service rate mu must be positive")

    # Compute stationary distribution pi of underlying Markov chain
    pi = _compute_stationary_distribution(C, D)

    # Compute mean arrival rate lambda
    e = np.ones((M, 1))
    lam = (pi @ D @ e).item()

    # Check stability condition
    rho = lam / mu
    if rho >= 1.0:
        raise ValueError(f"System is unstable (rho = {rho} >= 1)")

    # Compute R matrix
    R = _compute_r_matrix(C, D, mu)

    # Compute pi_0 = pi * (I - R)
    I = np.eye(M)
    pi0 = pi @ (I - R)

    # Determine N(epsilon) - queue length truncation
    N_epsilon = _determine_n_epsilon(pi0, R, D, lam, epsilon)

    # Compute theta (uniformization parameter)
    theta = np.max(np.abs(np.diag(C)))

    # Allocate result array
    result = np.zeros(len(x))

    # Main computation loop for each x value
    for idx, x_val in enumerate(x):
        # Determine K_max for uniformization
        theta_plus_mu = theta + mu
        mean_val = theta_plus_mu * x_val

        if mean_val > 0:
            L = max(0, int(np.floor(mean_val - 10 * np.sqrt(mean_val))))
            K_max = int(np.ceil(mean_val + 10 * np.sqrt(mean_val)))
        else:
            L = 0
            K_max = 0

        # Compute h_{n,k} recursively
        h = _compute_h_recursive(C, D, mu, N_epsilon, K_max, theta)

        # Compute W_bar(x) = (1/lambda) * sum over n of pi0 * R^n * D * sum over k of Poisson * h_{n,k}
        W_bar = 0.0

        Rpower = np.eye(M)
        for n in range(N_epsilon + 1):
            # Compute weight = pi0 * R^n * D
            weight = pi0 @ Rpower @ D

            # Compute sum over k
            sum_k = np.zeros((M, 1))
            for k in range(L, K_max + 1):
                poisson_term = _poisson_pmf(k, mean_val)
                if k <= len(h[n]) - 1:
                    sum_k = sum_k + poisson_term * h[n][k]

            W_bar += (1.0 / lam) * (weight @ sum_k).item()

            Rpower = Rpower @ R

        result[idx] = W_bar

    return result


__all__ = [
    'map_m1ps_cdf_respt',
]
