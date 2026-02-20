"""
MAP/BMAP/1 Queue Solver using GI/M/1-type analysis with ETAQA.

Solves a MAP/BMAP/1 queue where:
- Arrivals follow a Markovian Arrival Process (MAP)
- Service follows a Batch Markovian Arrival Process (BMAP) for batch service
- Single server

The queue is modeled as a GI/M/1-type Markov chain because:
- MAP arrivals increase level by exactly 1
- BMAP service can decrease level by 1, 2, 3, ... (batch sizes)

References:
    Riska, A., & Smirni, E. (2003). ETAQA: An Efficient Technique for the
    Analysis of QBD-Processes by Aggregation. Performance Evaluation, 54(2):151-177.
"""

import numpy as np
from dataclasses import dataclass
from typing import List, Tuple
import sys
import scipy.linalg as la

from .map_analysis import map_piq


@dataclass
class MAPBMAP1Result:
    """Result of MAP/BMAP/1 queue analysis."""
    mean_queue_length: float
    """Mean queue length E[N]"""

    utilization: float
    """Server utilization rho"""

    mean_response_time: float
    """Mean response time E[R]"""

    throughput: float
    """Throughput (arrival rate)"""

    pi: np.ndarray
    """Aggregated stationary probabilities [pi0, pi1, piStar]"""

    R: np.ndarray
    """R matrix"""

    mean_batch_size: float
    """Mean batch size of service"""


def solver_mam_map_bmap_1(C0: np.ndarray, C1: np.ndarray,
                          D: List[np.ndarray]) -> MAPBMAP1Result:
    """
    Solve a MAP/BMAP/1 queue using GI/M/1-type matrix-analytic methods.

    The MAP is specified by matrices (C0, C1) where:
    - C0: transitions without arrivals
    - C1: transitions triggering arrivals

    The BMAP for service is specified by matrices {D0, D1, D2, ..., DK} where:
    - D0: transitions without service completions
    - Dk: transitions triggering batch service of k customers (k >= 1)

    The GI/M/1-type structure for MAP/BMAP/1 is:
    ```
          B1  A0  0   0   ...
          B2  A1  A0  0   ...
      Q = B3  A2  A1  A0  ...
          ...
    ```

    Where:
      A0 = C1 otimes I_ms           (MAP arrival, level +1)
      A1 = C0 otimes I_ms + I_ma otimes D0  (phase changes, level 0)
      A_{k+1} = I_ma otimes D_k     (batch size k service, level -k)

    Args:
        C0: MAP matrix for transitions without arrivals
        C1: MAP matrix for arrivals
        D: BMAP service matrices as list [D0, D1, D2, ..., DK]

    Returns:
        MAPBMAP1Result with performance metrics
    """
    if len(D) < 1:
        raise ValueError("BMAP must have at least D0 matrix")
    if len(D) < 2:
        raise ValueError("BMAP must have at least D0 and D1 matrices")

    C0 = np.asarray(C0, dtype=float)
    C1 = np.asarray(C1, dtype=float)
    D = [np.asarray(Dk, dtype=float) for Dk in D]

    K = len(D) - 1  # Maximum batch size
    ma = C0.shape[0]  # Number of MAP phases
    ms = D[0].shape[0]  # Number of BMAP phases
    m = ma * ms  # Combined phases per level

    # Validate dimensions
    if C0.shape != (ma, ma):
        raise ValueError(f"C0 must be {ma}x{ma}")
    if C1.shape != (ma, ma):
        raise ValueError(f"C1 must be {ma}x{ma}")
    for i, Dk in enumerate(D):
        if Dk.shape != (ms, ms):
            raise ValueError(f"All BMAP matrices must be {ms}x{ms}")

    # Compute arrival rate from MAP
    piC = map_piq(C0, C1)
    eA = np.ones(ma)
    lambda_arr = piC @ C1 @ eA

    # Compute total service rate from BMAP
    D1_total = np.zeros((ms, ms))
    for k in range(1, K + 1):
        D1_total = D1_total + D[k]

    piD = map_piq(D[0], D1_total)
    eS = np.ones(ms)

    mu_total = 0.0
    for k in range(1, K + 1):
        rate_k = piD @ D[k] @ eS
        mu_total += k * rate_k

    # Mean batch size
    batch_rate = piD @ D1_total @ eS
    mean_batch_size = mu_total / batch_rate if batch_rate > 0 else 0.0

    # Utilization
    rho = lambda_arr / mu_total

    if rho >= 1.0:
        print(f"Warning: System is unstable (rho = {rho} >= 1). Results may be invalid.",
              file=sys.stderr)

    # Construct A matrix for GI/M/1-type (repeating part)
    # A = [A0; A1; A2; ...; A_{K+1}] with (K+2) blocks of size m x m
    A = np.zeros((m * (K + 2), m))

    # A0 = C1 \otimes I_ms (MAP arrival, level +1)
    I_ms = np.eye(ms)
    A0 = np.kron(C1, I_ms)
    A[0:m, :] = A0

    # A1 = C0 \otimes I_ms + I_ma \otimes D0 (phase changes, level 0)
    I_ma = np.eye(ma)
    A1 = np.kron(C0, I_ms) + np.kron(I_ma, D[0])
    A[m:2*m, :] = A1

    # A_{k+1} = I_ma \otimes D_k for k = 1, ..., K (batch service of size k)
    for k in range(1, K + 1):
        Ak = np.kron(I_ma, D[k])
        A[(k+1)*m:(k+2)*m, :] = Ak

    # Construct B matrix for boundary levels
    B = np.zeros((m * (K + 2), m))

    # B1: Level 0 transitions (no service from empty queue)
    B1 = np.kron(C0, I_ms) + np.kron(I_ma, D[0])
    # Add all service transitions back to level 0
    for k in range(1, K + 1):
        B1 = B1 + np.kron(I_ma, D[k])
    B[0:m, :] = B1

    # B_{j+1} for j = 1, ..., K: Transitions to level 0 from level j
    for j in range(1, K + 1):
        Bj = np.zeros((m, m))
        for k in range(j, K + 1):
            Bj = Bj + np.kron(I_ma, D[k])
        B[m + (j-1)*m : m + j*m, :] = Bj

    # Compute R matrix using GI/M/1 algorithm
    R = _gim1_r_etaqa(A)

    # Compute stationary probabilities using GI/M/1 ETAQA
    pi = _gim1_pi_etaqa(B, A, R, A0)

    # Compute mean queue length
    mean_queue_length = _gim1_qlen_etaqa(B, A, R, pi, 1)

    # Performance metrics
    mean_response_time = mean_queue_length / lambda_arr if lambda_arr > 0 else 0.0

    return MAPBMAP1Result(
        mean_queue_length=mean_queue_length,
        utilization=rho,
        mean_response_time=mean_response_time,
        throughput=lambda_arr,
        pi=pi,
        R=R,
        mean_batch_size=mean_batch_size
    )


def _gim1_r_etaqa(A: np.ndarray) -> np.ndarray:
    """
    Compute R matrix for GI/M/1-type Markov chains.

    R is the minimal nonnegative solution to:
      R = A0 + R*A1 + R^2*A2 + ... + R^K*A_K
    """
    m = A.shape[1]
    num_blocks = A.shape[0] // m

    # Extract A blocks
    blocks = [A[i*m:(i+1)*m, :] for i in range(num_blocks)]

    # Check if continuous time and uniformize if needed
    A1 = blocks[1]
    min_diag = np.min(np.diag(A1))

    if min_diag < 0:
        lamb = -min_diag
        uniformized = []
        for i, blk in enumerate(blocks):
            scaled = blk / lamb
            if i == 1:
                scaled = scaled + np.eye(m)
            uniformized.append(scaled)
    else:
        uniformized = blocks

    # Functional iteration for R
    R = np.zeros((m, m))
    max_iter = 200
    tol = 1e-14

    for _ in range(max_iter):
        # R_new = A0 + A1*R + A2*R^2 + ...
        R_new = uniformized[0].copy()
        R_pow = R.copy()
        for k in range(1, num_blocks):
            R_new = R_new + uniformized[k] @ R_pow
            R_pow = R_pow @ R

        diff = la.norm(R_new - R, ord='fro')
        if diff < tol:
            return R_new
        R = R_new

    return R


def _gim1_pi_etaqa(B: np.ndarray, A: np.ndarray, R: np.ndarray,
                   B0: np.ndarray) -> np.ndarray:
    """
    Compute aggregated stationary probabilities for GI/M/1-type Markov chain.
    """
    m = R.shape[0]

    # Compute g = stationary distribution of R
    g = _stat(R)

    # Compute pi0 using boundary equations
    I_min_R = np.eye(m) - R
    try:
        I_min_R_inv = la.inv(I_min_R)
    except la.LinAlgError:
        I_min_R_inv = la.pinv(I_min_R)

    # Extract B1
    B1 = B[0:m, :]

    # Compute pi0 from B1 + B0 * R * (I-R)^{-1}
    temp = B1 + B0 @ R @ I_min_R_inv
    pi0 = _stat(temp)

    # Normalize
    norm = (pi0 @ I_min_R_inv @ np.ones(m))
    pi0_norm = pi0 / norm if norm > 0 else pi0

    # pi1 = pi0 * R
    pi1 = pi0_norm @ R

    # piStar = pi1 * R
    pi_star = pi1 @ R

    # Concatenate
    result = np.concatenate([pi0_norm, pi1, pi_star])

    return result


def _gim1_qlen_etaqa(B: np.ndarray, A: np.ndarray, R: np.ndarray,
                     pi: np.ndarray, n: int) -> float:
    """
    Compute queue length moments for GI/M/1-type Markov chain.
    """
    m = R.shape[0]

    # Extract probability components
    pi0 = pi[0:m]
    pi1 = pi[m:2*m]
    pi_star = pi[2*m:3*m]

    if n == 1:
        # Mean queue length
        mean = np.sum(pi1) + 2.0 * np.sum(pi_star)
        return mean

    # Higher moments - approximate
    moment = np.sum(pi1)
    pi_star_sum = np.sum(pi_star)
    for k in range(2, 100):
        level_prob = pi_star_sum * (1.0 - pi_star_sum) ** (k - 2)
        moment += (k ** n) * level_prob
        if level_prob < 1e-15:
            break

    return moment


def _stat(A: np.ndarray) -> np.ndarray:
    """
    Compute stationary distribution of a stochastic/generator matrix.
    """
    n = A.shape[0]

    # Solve pi * (A - I) = 0 with sum(pi) = 1
    B = np.hstack([A.T, np.ones((n, 1))])
    y = np.zeros(n + 1)
    y[-1] = 1.0

    pi = la.lstsq(B.T, y)[0]

    # Ensure non-negative and normalized
    pi = np.maximum(pi, 0)
    s = np.sum(pi)
    if s > 0:
        pi = pi / s

    return pi


__all__ = [
    'MAPBMAP1Result',
    'solver_mam_map_bmap_1',
]
