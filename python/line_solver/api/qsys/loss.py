"""
Loss System Queue Analysis.

Native Python implementations for analyzing loss systems (finite buffer
queues where customers are rejected when buffer is full).

Key functions:
    qsys_mm1k_loss: M/M/1/K loss probability
    qsys_mg1k_loss: M/G/1/K loss probability (Niu-Cooper formula)
    qsys_mg1k_loss_mgs: M/G/1/K loss with MacGregor Smith method

References:
    Original MATLAB: matlab/src/api/qsys/qsys_*k_loss.m
    Niu-Cooper, "Transform-Free Analysis of M/G/1/K and Related Queues", 1993
"""

import numpy as np
from typing import Tuple, Callable, Optional
from scipy.integrate import quad


def qsys_mm1k_loss(lambda_val: float, mu: float, K: int) -> Tuple[float, float]:
    """
    Compute loss probability for M/M/1/K queue.

    Uses the closed-form formula for the M/M/1/K loss system where
    customers are rejected when the buffer (capacity K) is full.

    Args:
        lambda_val: Arrival rate
        mu: Service rate
        K: Buffer capacity (including customer in service)

    Returns:
        Tuple of (lossprob, rho):
            lossprob: Probability that an arriving customer is rejected
            rho: Offered load (lambda/mu)

    References:
        Original MATLAB: matlab/src/api/qsys/qsys_mm1k_loss.m
    """
    rho = lambda_val / mu

    if abs(rho - 1.0) < 1e-10:
        # Special case: rho = 1
        lossprob = 1.0 / (K + 1)
    else:
        lossprob = (1 - rho) / (1 - rho ** (K + 1)) * rho ** K

    return lossprob, rho


def qsys_mg1k_loss(lambda_val: float, service_pdf: Callable[[float], float],
                   K: int, max_t: float = 100.0) -> Tuple[float, float, float]:
    """
    Compute loss probability for M/G/1/K queue.

    Uses the Niu-Cooper transform-free analysis method for computing
    the stationary distribution and loss probability of M/G/1/K queues.

    Args:
        lambda_val: Arrival rate
        service_pdf: Probability density function of service time f(t)
        K: Buffer capacity
        max_t: Maximum integration time (default: 100)

    Returns:
        Tuple of (sigma, rho, lossprob):
            sigma: Normalizing constant term
            rho: Offered load
            lossprob: Probability of loss

    References:
        Original MATLAB: matlab/src/api/qsys/qsys_mg1k_loss.m
        Niu-Cooper, "Transform-Free Analysis of M/G/1/K", 1993
    """
    # Compute mean service time
    mean_service, _ = quad(lambda t: t * service_pdf(t), 0, max_t)
    rho = lambda_val * mean_service

    # Compute the stationary probabilities using Niu-Cooper method
    # This is a simplified implementation

    # Probability of n customers at departure epochs
    pi = np.zeros(K + 1)
    pi[0] = 1.0

    # Iterate to find stationary distribution
    for _ in range(1000):
        pi_new = np.zeros(K + 1)

        # pi_0 contribution
        for n in range(K + 1):
            # Probability of n arrivals during service
            a_n = _poisson_arrivals(lambda_val, n, service_pdf, max_t)
            pi_new[min(n, K)] += pi[0] * a_n

        # pi_j contributions for j > 0
        for j in range(1, K + 1):
            for n in range(K - j + 2):
                a_n = _poisson_arrivals(lambda_val, n, service_pdf, max_t)
                new_state = min(j + n - 1, K)
                if new_state >= 0:
                    pi_new[new_state] += pi[j] * a_n

        # Normalize
        total = np.sum(pi_new)
        if total > 0:
            pi_new /= total

        if np.linalg.norm(pi_new - pi) < 1e-10:
            break

        pi = pi_new

    # Loss probability is related to probability of being in state K
    # at arrival epochs (PASTA)
    sigma = pi[K]
    lossprob = sigma * rho / (1 + sigma * (rho - 1)) if rho > 0 else 0.0

    return sigma, rho, lossprob


def _poisson_arrivals(lambda_val: float, n: int,
                      service_pdf: Callable[[float], float],
                      max_t: float) -> float:
    """Compute probability of n Poisson arrivals during service."""
    from scipy.special import factorial

    def integrand(t):
        return (lambda_val * t) ** n * np.exp(-lambda_val * t) / factorial(n) * service_pdf(t)

    result, _ = quad(integrand, 0, max_t)
    return result


def qsys_mg1k_loss_mgs(lambda_val: float, mu: float, cs: float,
                       K: int) -> Tuple[float, float]:
    """
    Compute loss probability for M/G/1/K using MacGregor Smith approximation.

    Uses the MacGregor Smith approximation which provides good accuracy
    for moderate coefficient of variation.

    Args:
        lambda_val: Arrival rate
        mu: Service rate
        cs: Coefficient of variation of service time
        K: Buffer capacity

    Returns:
        Tuple of (lossprob, rho):
            lossprob: Probability of loss
            rho: Offered load

    References:
        Original MATLAB: matlab/src/api/qsys/qsys_mg1k_loss_mgs.m
        MacGregor Smith, "Queueing with capacity and performance"
    """
    rho = lambda_val / mu

    if rho >= 1:
        # For rho >= 1, use limiting approximation
        lossprob = 1.0 - 1.0 / K
        return lossprob, rho

    # M/M/1/K loss as base
    lossprob_mm1k, _ = qsys_mm1k_loss(lambda_val, mu, K)

    # Correction factor based on coefficient of variation
    if abs(cs - 1.0) < 1e-10:
        lossprob = lossprob_mm1k
    else:
        # MacGregor Smith correction
        # For cs < 1 (more regular), loss decreases
        # For cs > 1 (more bursty), loss increases
        correction = 1 + (cs ** 2 - 1) / 2 * (rho / (1 - rho + rho ** K))
        correction = max(0.1, min(10.0, correction))  # Bound correction
        lossprob = lossprob_mm1k * correction

    lossprob = min(1.0, max(0.0, lossprob))

    return lossprob, rho


def qsys_mxm1(lambda_vec: np.ndarray, mu: float,
              batch_sizes: Optional[np.ndarray] = None
              ) -> Tuple[float, float, float]:
    """
    Analyze M^X/M/1 queue (batch arrivals).

    Computes performance metrics for a queue with batch Poisson arrivals
    and exponential service.

    Args:
        lambda_vec: Vector of arrival rates for each batch size
                   (lambda_vec[k] = rate of batches of size k+1)
        mu: Service rate
        batch_sizes: Optional array of batch sizes (default: 1, 2, 3, ...)

    Returns:
        Tuple of (L, W, rho):
            L: Mean number in system
            W: Mean response time
            rho: Utilization

    References:
        Original MATLAB: matlab/src/api/qsys/qsys_mxm1.m
    """
    lambda_vec = np.asarray(lambda_vec, dtype=float).flatten()

    if batch_sizes is None:
        batch_sizes = np.arange(1, len(lambda_vec) + 1)

    # Total arrival rate of batches
    lambda_total = np.sum(lambda_vec)

    # Mean batch size
    if lambda_total > 0:
        E_X = np.sum(batch_sizes * lambda_vec) / lambda_total
        E_X2 = np.sum(batch_sizes ** 2 * lambda_vec) / lambda_total
    else:
        E_X = 1.0
        E_X2 = 1.0

    # Effective arrival rate of customers
    lambda_eff = lambda_total * E_X

    # Utilization
    rho = lambda_eff / mu

    if rho >= 1:
        raise ValueError(f"System is unstable: utilization rho = {rho:.4g} >= 1")

    # Mean number in system (batch M/M/1 formula)
    # L = rho + (lambda * E[X^2]) / (2 * mu * (1 - rho))
    L = rho + (lambda_total * E_X2) / (2 * mu * (1 - rho))

    # Mean response time (Little's law)
    W = L / lambda_eff

    return L, W, rho


__all__ = [
    'qsys_mm1k_loss',
    'qsys_mg1k_loss',
    'qsys_mg1k_loss_mgs',
    'qsys_mxm1',
]
