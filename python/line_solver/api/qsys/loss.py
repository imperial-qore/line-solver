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
                   K: int, max_t: Optional[float] = None) -> Tuple[float, float, float]:
    """
    Exact M/G/1/K loss probability via the Markov chain embedded at
    service-start epochs (transform-free analysis in the spirit of
    Niu-Cooper).

    State: number of customers waiting in the queue immediately after a
    service start, q in {0,...,K-2} (capacity K includes the job in service;
    just after a departure at most K-1 jobs remain, one of which enters
    service). With a_j = P(j Poisson arrivals during a service time)::

        q=0 : if no arrival occurs during the service the system empties and
              the next service starts with the next arrival (q'=0), so both
              a_0 and a_1 lead to q'=0 and j>=2 arrivals lead to q'=j-1;
        q>=1: q' = q-1+j, with arrivals beyond the free capacity lost
              (aggregated in the last column).

    The loss probability follows from the renewal-reward argument::

        E[cycle] = E[S] + sigma_0*a_0/lambda,  lambda_eff = 1/E[cycle],
        P_loss = 1 - lambda_eff/lambda = 1 - 1/(rho + sigma_0*a_0)

    where sigma is the stationary distribution at service-start epochs.

    Args:
        lambda_val: Arrival rate
        service_pdf: Probability density function of service time f(t)
        K: Buffer capacity (including the customer in service)
        max_t: Maximum integration time (default: smallest horizon covering
            the service-time distribution mass to within 1e-10)

    Returns:
        Tuple of (sigma0, rho, lossprob):
            sigma0: Stationary probability of an empty queue at service-start
                epochs

            rho: Offered load
            lossprob: Probability of loss

    References:
        Original MATLAB: matlab/src/api/qsys/qsys_mg1k_loss.m
        Niu-Cooper, "Transform-Free Analysis of M/G/1/K", 1993
    """
    from ..mc.dtmc import dtmc_solve, dtmc_makestochastic

    if max_t is None:
        # see _kb/03-api-layer.md for rationale
        max_t = 1.0 / lambda_val
        for _ in range(60):
            mass, _ = quad(service_pdf, 0, max_t)
            if mass >= 1.0 - 1e-10:
                break
            max_t *= 2.0

    # Compute mean service time
    mean_service, _ = quad(lambda t: t * service_pdf(t), 0, max_t)
    rho = lambda_val * mean_service

    # Arrival probabilities a_j = P(j arrivals during a service time)
    a = np.zeros(max(K - 1, 2))
    for j in range(1001):
        aj = _poisson_arrivals(lambda_val, j, service_pdf, max_t)
        if j >= len(a):
            a = np.append(a, 0.0)
        a[j] = aj
        if aj < 1e-12:
            break

    # Embedded chain at service-start epochs, states q=0..K-2
    n = K - 1
    P = np.zeros((n, n))

    # row 0 (q=0): idle period after an empty departure epoch
    P[0, 0] = a[0] + a[1]
    for i in range(1, K - 2):
        P[0, i] = a[i + 1]
    P[0, n - 1] = 1.0 - np.sum(P[0, :n - 1])

    # row 1 (q=1): q' = number of arrivals during the service (capped)
    if n >= 2:
        for i in range(K - 2):
            P[1, i] = a[i]
        P[1, n - 1] = 1.0 - np.sum(P[1, :n - 1])

    # rows j>=2 (q=j): q' = q-1+arrivals (capped)
    for j in range(2, n):
        for i in range(j - 1, K - 2):
            P[j, i] = a[i - j + 1]
        P[j, n - 1] = 1.0 - np.sum(P[j, :n - 1])

    P = dtmc_makestochastic(P)
    sigma = dtmc_solve(P)

    sigma0 = float(sigma[0])
    lossprob = 1.0 - 1.0 / (sigma0 * a[0] + rho)

    return sigma0, rho, lossprob


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
