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


def qsys_mg1k_loss_mgs(lambda_val: float, mu: float, mu_scv: float,
                       K: int) -> Tuple[float, float]:
    """
    Compute loss probability for M/G/1/K using MacGregor Smith approximation.

    Matches MATLAB qsys_mg1k_loss_mgs.m exactly.

    Args:
        lambda_val: Arrival rate
        mu: Service rate
        mu_scv: Squared coefficient of variation of service time
        K: Buffer capacity

    Returns:
        Tuple of (lossprob, rho):
            lossprob: Probability of loss
            rho: Offered load

    References:
        Original MATLAB: matlab/src/api/qsys/qsys_mg1k_loss_mgs.m
        J. MacGregor Smith, "Optimal Design and Performance Modelling of M/G/1/K Queueing Systems"
    """
    rho = lambda_val / mu
    s = np.sqrt(mu_scv)
    sqrt_rho = np.sqrt(rho)

    lossprob_num = rho**((sqrt_rho * s**2 - sqrt_rho + 2 * K) / (2 + sqrt_rho * s**2 - sqrt_rho)) * (rho - 1)
    lossprob_den = rho**(2 * (1 + sqrt_rho * s**2 - sqrt_rho + K) / (2 + sqrt_rho * s**2 - sqrt_rho)) - 1
    lossprob = lossprob_num / lossprob_den

    return lossprob, rho


def qsys_mxm1(lambda_batch: float, mu: float,
              E_X_or_batch_sizes, E_X2_or_pmf,
              mode: Optional[str] = None) -> Tuple[float, float, float, float]:
    """
    Analyze MX/M/1 queue with batch arrivals.

    Matches MATLAB qsys_mxm1.m exactly.

    Three input formats:
        1. Moment-based: qsys_mxm1(lambda_batch, mu, E_X, E_X2)
        2. PMF-based:    qsys_mxm1(lambda_batch, mu, batch_sizes, pmf)
        3. Variance:     qsys_mxm1(lambda_batch, mu, E_X, Var_X, 'variance')

    Args:
        lambda_batch: Batch arrival rate
        mu: Service rate
        E_X_or_batch_sizes: Mean batch size (scalar) or array of batch sizes
        E_X2_or_pmf: Second moment of batch size, PMF, or variance
        mode: Optional 'variance' flag for variance-based input

    Returns:
        Tuple of (W, Wq, U, Q):
            W: Mean time in system
            Wq: Mean waiting time in queue
            U: Server utilization
            Q: Mean queue length (including service)

    References:
        Original MATLAB: matlab/src/api/qsys/qsys_mxm1.m
    """
    if mode is not None and mode.lower() == 'variance':
        # Format 3: Variance-based
        E_X = float(E_X_or_batch_sizes)
        Var_X = float(E_X2_or_pmf)
        E_X2 = Var_X + E_X**2
    elif hasattr(E_X_or_batch_sizes, '__len__') and len(np.asarray(E_X_or_batch_sizes).flatten()) > 1:
        # Format 2: PMF-based
        batch_sizes = np.asarray(E_X_or_batch_sizes, dtype=float).flatten()
        pmf = np.asarray(E_X2_or_pmf, dtype=float).flatten()
        if len(batch_sizes) != len(pmf):
            raise ValueError("Batch sizes and PMF must have the same length")
        pmf = pmf / np.sum(pmf)
        E_X = np.sum(batch_sizes * pmf)
        E_X2 = np.sum(batch_sizes**2 * pmf)
    else:
        # Format 1: Moment-based (default)
        E_X = float(E_X_or_batch_sizes)
        E_X2 = float(E_X2_or_pmf)

    # Compute effective job arrival rate
    lambda_eff = lambda_batch * E_X

    # Compute utilization
    rho = lambda_eff / mu

    if rho >= 1:
        raise ValueError(f"System is unstable: rho = {rho:.6f} >= 1")

    # Mean waiting time in queue (MATLAB formula):
    # Wq = rho/(mu*(1-rho)) + (E[X^2] - E[X])/(2*mu*E[X]*(1-rho))
    Wq = rho / (mu * (1 - rho)) + (E_X2 - E_X) / (2 * mu * E_X * (1 - rho))

    # Mean time in system
    W = Wq + 1 / mu

    # Server utilization
    U = rho

    # Mean queue length (Little's Law)
    Q = lambda_eff * W

    return W, Wq, U, Q


__all__ = [
    'qsys_mm1k_loss',
    'qsys_mg1k_loss',
    'qsys_mg1k_loss_mgs',
    'qsys_mxm1',
]
