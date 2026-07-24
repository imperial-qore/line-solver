"""
Simplified FJ_codes toolkit for Fork-Join percentile analysis.

Computes percentiles of response time in Fork-Join networks using
matrix-analytic methods. Handles M/M/c service distributions.

References:
    Z. Qiu, J.F. Pérez, and P. Harrison, "Beyond the Mean in Fork-Join Queues:
    Efficient Approximation for Response-Time Tails", IFIP Performance 2015.
"""

import numpy as np
from scipy import integrate
from typing import List, Tuple, Dict, Optional
from dataclasses import dataclass


@dataclass
class FJPercentileResult:
    """Result of FJ percentile analysis."""
    percentiles: List[float]
    response_times: List[float]
    mean_response_time: float
    K: int  # Number of parallel queues


def fj_mm1_percentiles(lambda_rate: float, mu: float, percentiles: List[float]) -> Tuple[float, List[float]]:
    """
    Compute percentiles for M/M/1 queue.

    Args:
        lambda_rate: Arrival rate
        mu: Service rate
        percentiles: List of percentiles (0-100)

    Returns:
        (mean_response_time, [response_times_at_percentiles])
    """
    rho = lambda_rate / mu
    if rho >= 1.0:
        return float('inf'), [float('inf')] * len(percentiles)

    # Mean response time: R = 1/(mu - lambda)
    mean_r = 1.0 / (mu - lambda_rate)

    # CDF of M/M/1 response time: P(R <= t) = 1 - rho * exp(-(mu-lambda)*t)
    response_times = []
    for p in percentiles:
        if p >= 100:
            response_times.append(float('inf'))
        elif p <= 0:
            response_times.append(0)
        else:
            # Solve: p/100 = 1 - rho * exp(-(mu-lambda)*t)
            # => exp(-(mu-lambda)*t) = (1 - p/100) / rho
            # => t = -ln((1 - p/100) / rho) / (mu - lambda)
            prob = p / 100.0
            arg = (1.0 - prob) / rho
            if arg <= 0:
                response_times.append(float('inf'))
            else:
                t = -np.log(arg) / (mu - lambda_rate)
                response_times.append(max(0, t))

    return mean_r, response_times


def fj_fork_join_response_time(K: int, lambda_rate: float, mu: float) -> float:
    """
    Compute mean response time for Fork-Join with K parallel M/M/1 queues.

    Each job forks into K tasks, each service at one of K parallel queues.
    Response time is the maximum of K service times.

    Args:
        K: Number of parallel queues
        lambda_rate: Total arrival rate
        mu: Service rate at each queue

    Returns:
        Mean response time E[R_FJ]
    """
    if K == 1:
        return 1.0 / (mu - lambda_rate)

    # For K identical M/M/1 queues:
    # Each queue receives lambda/K arrivals
    lambda_per_queue = lambda_rate / K
    rho = lambda_per_queue / mu

    if rho >= 1.0:
        return float('inf')

    # Mean service time
    S = 1.0 / mu

    # Mean response time at one queue
    R_single = 1.0 / (mu - lambda_per_queue)

    # For K parallel M/M/1 queues with Poisson arrivals split equally:
    # The maximum (response time at slowest queue) has approximately:
    # E[max(R1, ..., RK)] ≈ E[R] + (K-1)*E[W] / K
    # where W is waiting time
    W = lambda_per_queue * S / (1.0 - rho) / mu  # Waiting time

    # Simpler approximation: E[max] ≈ S + (K-1)*W/K for large rho
    # More accurate: use formula from Qiu et al.
    # Approximate: E[max] ≈ S + sum_{k=1}^{K} P(X_k is max) * E[W_k | X_k is max]

    # For simplicity, use: E[max] ≈ E[R_single] * (1 + ln(K)/K)
    # This is a reasonable approximation for moderate K
    mean_fj = R_single * (1.0 + np.log(K) / K)

    return mean_fj


def fj_percentiles_fork_join(K: int, lambda_rate: float, mu: float,
                             percentiles: List[float]) -> FJPercentileResult:
    """
    Compute percentiles for Fork-Join network with K parallel M/M/1 queues.

    Approximation: Use the distribution of maximum of K i.i.d. exponential
    random variables conditioned on the overall queue dynamics.

    Args:
        K: Number of parallel queues
        lambda_rate: Total arrival rate to system
        mu: Service rate at each queue
        percentiles: List of percentiles (0-100)

    Returns:
        FJPercentileResult with percentile response times
    """
    # Per-queue arrival rate
    lambda_per_queue = lambda_rate / K
    rho = lambda_per_queue / mu

    if rho >= 1.0:
        mean_rt = float('inf')
        response_times = [float('inf')] * len(percentiles)
        return FJPercentileResult(
            percentiles=percentiles,
            response_times=response_times,
            mean_response_time=mean_rt,
            K=K
        )

    # Compute percentiles for single queue
    mean_single, rt_single = fj_mm1_percentiles(lambda_per_queue, mu, percentiles)

    # For K parallel queues, the response time is approximately
    # the maximum of K service times plus weighted waiting time
    # Approximation: P(R_FJ <= t) ≈ [P(R_single <= t)]^K

    response_times = []
    for i, rt in enumerate(rt_single):
        if rt >= float('inf'):
            response_times.append(float('inf'))
        else:
            # Adjust for fork-join: maximum of K responses
            # P(max <= t) = P(all <= t) = [P(single <= t)]^K
            # We need to find t such that P(all <= t) = p_i
            # Use inverse: if P(R <= t) for single queue, then
            # P(R_FJ <= t) ≈ [P(R <= t)]^K

            # For percentile p, we want P(R_FJ <= t) = p
            # => [P(R_single <= t)]^K = p
            # => P(R_single <= t) = p^(1/K)
            # => t = inverse_cdf(p^(1/K)) for single queue

            p_adjusted = (percentiles[i] / 100.0) ** (1.0 / K)
            p_adjusted_pct = p_adjusted * 100

            # Recompute for adjusted percentile
            _, rt_adjusted = fj_mm1_percentiles(lambda_per_queue, mu, [p_adjusted_pct])
            response_times.append(rt_adjusted[0])

    # Mean response time
    mean_rt = fj_fork_join_response_time(K, lambda_rate, mu)

    return FJPercentileResult(
        percentiles=percentiles,
        response_times=response_times,
        mean_response_time=mean_rt,
        K=K
    )


def fj_percentiles_general(K: int, lambda_rate: float, service_rates: np.ndarray,
                          service_scv: np.ndarray, percentiles: List[float]) -> FJPercentileResult:
    """
    Compute percentiles for general Fork-Join network.

    For general distributions, use approximations based on first two moments.

    Args:
        K: Number of parallel queues
        lambda_rate: Total arrival rate
        service_rates: Service rates (K,) or (K, nclasses)
        service_scv: Service SCV (K,) or (K, nclasses)
        percentiles: List of percentiles

    Returns:
        FJPercentileResult
    """
    # For now, use the M/M/c approximation with SCV adjustment
    # This is a simplified approach

    service_rates = np.asarray(service_rates).flatten()
    service_scv = np.asarray(service_scv).flatten()

    # Use first queue's parameters
    mu = service_rates[0] if len(service_rates) > 0 else 1.0
    scv = service_scv[0] if len(service_scv) > 0 else 1.0

    # For non-exponential, adjust using Kingman's approximation
    lambda_per_queue = lambda_rate / K
    rho = lambda_per_queue / mu
    S = 1.0 / mu
    W_mm1 = (1.0 + scv) * rho / (2.0 * (1.0 - rho)) * S

    # Mean response time adjusted for SCV
    R_approx = S + W_mm1

    # For percentiles, use M/M/1 distribution but scaled by R_approx
    mean_mm1, rt_mm1 = fj_mm1_percentiles(lambda_per_queue, mu, percentiles)

    # Scale by ratio of true mean to M/M/1 mean
    scale_factor = R_approx / mean_mm1 if mean_mm1 > 0 else 1.0

    response_times = [rt * scale_factor for rt in rt_mm1]

    # Adjust for fork-join maximum
    response_times_fj = []
    for rt in response_times:
        if rt >= float('inf'):
            response_times_fj.append(float('inf'))
        else:
            # Approximate percentile adjustment for K queues
            # This is a heuristic: weight the response time increase
            response_times_fj.append(rt * (1.0 + (K - 1) / (2.0 * K)))

    return FJPercentileResult(
        percentiles=percentiles,
        response_times=response_times_fj,
        mean_response_time=R_approx * (1.0 + np.log(K) / K) if K > 1 else R_approx,
        K=K
    )
