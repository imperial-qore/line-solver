"""
Akyildiz-Bolch AMVA method for multi-server BCMP networks.

This module implements the Akyildiz-Bolch linearizer algorithm for solving
queueing networks with multiple servers and various scheduling strategies.

References:
    Akyildiz, I.F. and Bolch, G., "Mean Value Analysis Approximation for
    Multiple Server Queueing Networks", Performance Evaluation, 1988.
"""

import numpy as np
from typing import Tuple, Dict, Optional, List
from dataclasses import dataclass
from enum import IntEnum


class SchedStrategy(IntEnum):
    """Scheduling strategies for queueing stations."""
    FCFS = 0
    LCFS = 1
    INF = 2  # Infinite server (delay)
    PS = 3   # Processor sharing


@dataclass
class AbAmvaResult:
    """Result from Akyildiz-Bolch AMVA algorithm."""
    QN: np.ndarray      # Queue lengths (M x K)
    UN: np.ndarray      # Utilization (M x K)
    RN: np.ndarray      # Residence times (M x K)
    CN: np.ndarray      # Cycle times (1 x K)
    XN: np.ndarray      # Throughput (1 x K)
    totiter: int        # Total iterations


def pfqn_ab_amva(
    D: np.ndarray,
    N: np.ndarray,
    V: np.ndarray,
    nservers: np.ndarray,
    sched: np.ndarray,
    fcfs_schmidt: bool = False,
    marginal_prob_method: str = 'ab'
) -> AbAmvaResult:
    """
    Akyildiz-Bolch AMVA method for multi-server BCMP networks.

    Args:
        D: Service time matrix (M x K)
        N: Population vector (1 x K)
        V: Visit ratio matrix (M x K)
        nservers: Number of servers at each station (M x 1)
        sched: Scheduling strategies for each station (M x 1)
        fcfs_schmidt: Whether to use Schmidt formula for FCFS stations
        marginal_prob_method: Method for marginal probability ('ab' or 'scat')

    Returns:
        AbAmvaResult containing queue lengths, utilization, residence times,
        cycle times, throughput, and iteration count.
    """
    D = np.atleast_2d(D)
    N = np.atleast_1d(N).astype(float)
    V = np.atleast_2d(V)
    nservers = np.atleast_1d(nservers).astype(float)
    sched = np.atleast_1d(sched)

    M, K = D.shape

    QN, UN, RN, CN, XN, totiter = _ab_linearizer(
        K, M, N, nservers, sched, V, D, fcfs_schmidt, marginal_prob_method
    )

    return AbAmvaResult(QN=QN, UN=UN, RN=RN, CN=CN, XN=XN, totiter=totiter)


def _ab_linearizer(
    K: int,
    M: int,
    population: np.ndarray,
    nservers: np.ndarray,
    sched_type: np.ndarray,
    v: np.ndarray,
    s: np.ndarray,
    fcfs_schmidt: bool,
    marginal_prob_method: str
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    """Akyildiz-Bolch linearizer method for multi-server BCMP networks."""

    # Initialize queue length matrix
    L = np.zeros((M, K))
    for i in range(M):
        for r in range(K):
            L[i, r] = population[r] / M

    # Initialize L_without_r matrices
    l_without_r: List[np.ndarray] = []
    for r in range(K):
        l_wr = np.zeros((M, K))
        for i in range(M):
            for t in range(K):
                if r == t:
                    l_wr[i, t] = (population[r] - 1) / M
                else:
                    l_wr[i, t] = L[i, r]
        l_without_r.append(l_wr)

    # Fractional changes matrix (initialized to 0)
    D_frac = np.zeros((M, K, K))

    # STEP 1: Apply core at full population
    L_updated, _, _, _, _, _ = _pfqn_ab_core(
        K, M, population, nservers, sched_type, v, s, 100, D_frac, L,
        fcfs_schmidt, marginal_prob_method
    )

    # STEP 2: Apply core at N-e_k populations
    for r in range(K):
        population_without_c = population.copy()
        population_without_c[r] = population[r] - 1

        l_without_c = np.zeros((M, K))
        for j in range(M):
            for c in range(K):
                l_without_c[j, c] = l_without_r[c][j, c]

        ret_q, _, _, _, _, _ = _pfqn_ab_core(
            K, M, population_without_c, nservers, sched_type, v, s, 100,
            D_frac, l_without_c, fcfs_schmidt, marginal_prob_method
        )

        for j in range(M):
            for c in range(K):
                l_without_r[c][j, r] = ret_q[j, c]

    # STEP 3: Compute estimates of F_mk(N) and F_mk(N-e_j)
    for i in range(M):
        for r in range(K):
            F_ir = L_updated[i, r] / population[r] if population[r] > 0 else 0
            for t in range(K):
                if r == t:
                    divisor = population[r] - 1
                else:
                    divisor = population[r]

                if divisor != 0:
                    F_irt = l_without_r[r][i, t] / divisor
                else:
                    F_irt = 0

                if np.isnan(F_irt):
                    F_irt = 0

                D_frac[i, r, t] = F_irt - F_ir

    # STEP 4: Apply core at full population using L values from step 1 and D values from step 3
    QN, UN, RN, CN, XN, totiter = _pfqn_ab_core(
        K, M, population, nservers, sched_type, v, s, 100, D_frac, L_updated,
        fcfs_schmidt, marginal_prob_method
    )

    return QN, UN, RN, CN, XN, totiter


def _pfqn_ab_core(
    K: int,
    M: int,
    population: np.ndarray,
    nservers: np.ndarray,
    sched_type: np.ndarray,
    v: np.ndarray,
    s: np.ndarray,
    maxiter: int,
    D_frac: np.ndarray,
    l_in: np.ndarray,
    fcfs_schmidt: bool,
    marginal_prob_method: str
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    """Akyildiz-Bolch core method for multi-server BCMP networks."""

    L = l_in.copy()
    tol = 1 / (4000 + 16 * np.sum(population))

    totiter = 0
    W = np.zeros((M, K))

    while totiter < maxiter:
        F = np.zeros((M, K))
        l_without_j = np.zeros((M, K, K))

        for i in range(M):
            for r in range(K):
                if population[r] > 0:
                    F[i, r] = L[i, r] / population[r]
                else:
                    F[i, r] = 0

        for i in range(M):
            for r in range(K):
                for t in range(K):
                    if r == t:
                        scalar = population[r] - 1
                    else:
                        scalar = population[r]
                    l_without_j[i, r, t] = scalar * (F[i, r] + D_frac[i, r, t])

        for i in range(M):
            for r in range(K):
                if sched_type[i] == SchedStrategy.INF:
                    W[i, r] = s[i, r]
                elif nservers[i] == 1:
                    total_queue_length = 0
                    for c in range(K):
                        total_queue_length += l_without_j[i, c, r]
                    W[i, r] = s[i, r] * (1 + total_queue_length)
                elif fcfs_schmidt and sched_type[i] == SchedStrategy.FCFS:
                    # Schmidt formula for FCFS (simplified)
                    wait_time = 0
                    nvec = _pprod_init(population)
                    while nvec is not None:
                        if nvec[r] > 0:
                            bcn = _get_bcn_for_ab(s, i, r, nvec, K, int(nservers[i]))
                            prob = _get_marginal_prob(
                                _oner(nvec, r), _oner(population, r),
                                population, l_in[i, r], K
                            )
                            wait_time += bcn * prob
                        nvec = _pprod_next(nvec, population)
                    if wait_time <= 1e-3:
                        wait_time = 0
                    W[i, r] = wait_time
                else:
                    queue_length = 0
                    for j in range(K):
                        queue_length += l_without_j[i, j, r]

                    num_servers = int(nservers[i])
                    multi_server_weighted_ql = 0

                    if num_servers > 1:
                        population_without_r = population.copy()
                        population_without_r[r] = population_without_r[r] - 1
                        marginal_probs = _find_marginal_probs(
                            queue_length, num_servers, population_without_r,
                            r, marginal_prob_method
                        )

                        for j in range(num_servers - 1):
                            prob_j = marginal_probs.get(j, 0)
                            multi_server_weighted_ql += prob_j * (num_servers - j)

                    wait_time = (s[i, r] / num_servers) * (1 + queue_length + multi_server_weighted_ql)
                    W[i, r] = wait_time

        # Calculate cycle time for each class
        CN = np.zeros(K)
        for r in range(K):
            cycle_time = 0
            for i in range(M):
                cycle_time += v[i, r] * W[i, r]
            CN[r] = cycle_time

        # Calculate queue length L_ir = N_r * W_ir / C_r
        iteration_queue_length = np.zeros((M, K))
        for i in range(M):
            for r in range(K):
                if CN[r] > 0:
                    queue_length = population[r] * (v[i, r] * W[i, r] / CN[r])
                else:
                    queue_length = 0
                iteration_queue_length[i, r] = queue_length

        # Check convergence
        max_difference = 0
        for i in range(M):
            for r in range(K):
                if population[r] > 0:
                    difference = abs(L[i, r] - iteration_queue_length[i, r]) / population[r]
                    max_difference = max(max_difference, difference)

        totiter += 1
        L = iteration_queue_length

        if max_difference < tol:
            break

    # Calculate throughput
    XN = np.zeros(K)
    for r in range(K):
        if W[0, r] > 0:
            XN[r] = L[0, r] / W[0, r]
        else:
            XN[r] = 0

    # Calculate utilization
    UN = np.zeros((M, K))
    for i in range(M):
        for r in range(K):
            if s[i, r] > 0:
                if sched_type[i] == SchedStrategy.INF:
                    UN[i, r] = XN[r] * s[i, r]
                else:
                    UN[i, r] = (XN[r] * s[i, r]) / nservers[i]
            else:
                UN[i, r] = 0

    QN = L

    return QN, UN, W, CN, XN, totiter


def _get_bcn_for_ab(
    D: np.ndarray, i: int, c: int, nvec: np.ndarray, K: int, ns: int
) -> float:
    """Calculate Bcn for AB algorithm."""
    bcn = D[i, c]
    if np.sum(nvec) > 1:
        eps_val = 1e-12
        sum_val = 0
        for t in range(K):
            sum_val += nvec[t] * D[i, t]
        bcn = bcn + (max(0, np.sum(nvec) - ns) / max(ns * (np.sum(nvec) - 1), eps_val) * (sum_val - D[i, c]))
    return bcn


def _get_marginal_prob(
    n: np.ndarray, k: np.ndarray, k_pop: np.ndarray, l_jr: float, R: int
) -> float:
    """Compute marginal probability using binomial distribution."""
    from scipy.special import comb

    prob = 1.0
    for r in range(R):
        frac = l_jr / k_pop[r] if k_pop[r] > 0 else 0
        if frac != 0 and k_pop[r] > 0 and 0 <= n[r] <= k_pop[r]:
            term1 = comb(int(k_pop[r]), int(n[r]), exact=True)
            term2 = frac ** n[r]
            term3 = (1 - frac) ** (k_pop[r] - n[r])
            prob = prob * (term1 * term2 * term3)
    return prob


def _find_marginal_probs(
    avg_jobs: float,
    num_servers: int,
    population: np.ndarray,
    class_idx: int,
    marginal_prob_method: str
) -> Dict[int, float]:
    """Finds marginal probabilities using specified method."""

    marginal_probs: Dict[int, float] = {}

    if marginal_prob_method == 'scat':
        floor_val = int(np.floor(avg_jobs))
        ceil_val = floor_val + 1
        marginal_probs[floor_val] = ceil_val - avg_jobs
        marginal_probs[ceil_val] = avg_jobs - floor_val
        return marginal_probs

    # AB method
    ALPHA = 45.0
    BETA = 0.7

    w = _weight_fun(population, ALPHA, BETA)

    floor_val = int(np.floor(avg_jobs))
    ceiling = floor_val + 1
    max_val = min((2 * floor_val) + 1, num_servers - 2)

    for j in range(max_val + 1):
        if j <= floor_val:
            l_dist = floor_val - j
            lower_val = floor_val - l_dist
            upper_val = ceiling + l_dist
            if l_dist > 25:
                prob = 0
            else:
                if floor_val < population[class_idx]:
                    if upper_val != lower_val:
                        prob = w[floor_val, l_dist] * ((upper_val - avg_jobs) / (upper_val - lower_val))
                    else:
                        prob = 0
                else:
                    prob = 0
            marginal_probs[j] = prob
        else:
            u_dist = j - ceiling
            if u_dist > 25:
                marginal_probs[j] = 0
            elif j > population[class_idx] - 1 and u_dist < 25:
                existing_prob = marginal_probs.get(int(population[class_idx] - 1), 0)
                mp_floor_udist = marginal_probs.get(floor_val - u_dist, 0)
                new_prob = existing_prob + (w[floor_val, u_dist] - mp_floor_udist)
                marginal_probs[int(population[class_idx] - 1)] = new_prob
            else:
                mp_floor_udist = marginal_probs.get(floor_val - u_dist, 0)
                new_prob = w[floor_val, u_dist] - mp_floor_udist
                marginal_probs[j] = new_prob

    return marginal_probs


def _weight_fun(population: np.ndarray, alpha: float, beta: float) -> np.ndarray:
    """Computes weight function for marginal probability calculation."""

    max_class_population = int(np.max(population))

    # Calculate scaling function PR
    scaling_fun = np.zeros(max_class_population + 1)
    if max_class_population >= 1:
        scaling_fun[1] = alpha
        for n in range(2, max_class_population + 1):
            scaling_fun[n] = beta * scaling_fun[n - 1]

    # Calculate weight function W
    w = np.zeros((max_class_population + 1, max_class_population + 1))
    w[0, 0] = 1.0

    for l in range(1, max_class_population + 1):
        for j in range(l):
            w[l, j] = w[l - 1, j] - (w[l - 1, j] * scaling_fun[l]) / 100.0
        sum_val = 0
        for j in range(l):
            sum_val += w[l, j]
        w[l, l] = 1 - sum_val

    return w


def _pprod_init(population: np.ndarray) -> np.ndarray:
    """Initialize population product iterator."""
    return np.zeros_like(population)


def _pprod_next(nvec: np.ndarray, population: np.ndarray) -> Optional[np.ndarray]:
    """Get next population vector in product iteration."""
    K = len(nvec)
    result = nvec.copy()

    for i in range(K):
        if result[i] < population[i]:
            result[i] += 1
            return result
        result[i] = 0

    return None


def _oner(vec: np.ndarray, idx: int) -> np.ndarray:
    """Decrement vector at given index by 1."""
    result = vec.copy()
    result[idx] = max(0, result[idx] - 1)
    return result


def pfqn_ab_core(
    K: int,
    M: int,
    population: np.ndarray,
    nservers: np.ndarray,
    sched_type: np.ndarray,
    v: np.ndarray,
    s: np.ndarray,
    maxiter: int,
    D_frac: np.ndarray,
    l_in: np.ndarray,
    fcfs_schmidt: bool = False,
    marginal_prob_method: str = 'ab'
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    """
    Akyildiz-Bolch core method for multi-server BCMP networks.

    Public wrapper for the internal core algorithm of the Akyildiz-Bolch
    linearizer method.

    Args:
        K: Number of classes
        M: Number of stations
        population: Population vector (K,)
        nservers: Number of servers at each station (M,)
        sched_type: Scheduling strategies for each station (M,)
        v: Visit ratio matrix (M x K)
        s: Service time matrix (M x K)
        maxiter: Maximum iterations
        D_frac: Fractional changes matrix (M x K x K)
        l_in: Initial queue length matrix (M x K)
        fcfs_schmidt: Whether to use Schmidt formula for FCFS stations
        marginal_prob_method: Method for marginal probability ('ab' or 'scat')

    Returns:
        Tuple of (QN, UN, RN, CN, XN, totiter):
            QN: Queue lengths (M x K)
            UN: Utilization (M x K)
            RN: Residence times (M x K)
            CN: Cycle times (K,)
            XN: Throughput (K,)
            totiter: Total iterations

    References:
        Akyildiz, I.F. and Bolch, G., "Mean Value Analysis Approximation for
        Multiple Server Queueing Networks", Performance Evaluation, 1988.
    """
    return _pfqn_ab_core(
        K, M, population, nservers, sched_type, v, s, maxiter, D_frac, l_in,
        fcfs_schmidt, marginal_prob_method
    )


__all__ = [
    'pfqn_ab_amva',
    'pfqn_ab_core',
    'AbAmvaResult',
    'SchedStrategy',
]
