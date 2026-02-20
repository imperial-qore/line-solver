"""
BMAP/MAP/1 Queue Solver using M/G/1-type analysis with ETAQA.

Solves a BMAP/MAP/1 queue where:
- Arrivals follow a Batch Markovian Arrival Process (BMAP)
- Service follows a Markovian Arrival Process (MAP) for single service
- Single server

The queue is modeled as an M/G/1-type Markov chain because:
- BMAP arrivals can increase level by 1, 2, 3, ... (batch sizes)
- MAP service decreases level by exactly 1

References:
    Riska, A., & Smirni, E. (2003). ETAQA: An Efficient Technique for the
    Analysis of QBD-Processes by Aggregation. Performance Evaluation, 54(2):151-177.
"""

import numpy as np
from dataclasses import dataclass
from typing import List, Tuple
import sys

from ..smc import mg1_g_etaqa, mg1_pi_etaqa, mg1_qlen_etaqa
from .map_analysis import map_piq


@dataclass
class BMAPMAP1Result:
    """Result of BMAP/MAP/1 queue analysis."""
    mean_queue_length: float
    """Mean queue length E[N]"""

    utilization: float
    """Server utilization rho"""

    mean_response_time: float
    """Mean response time E[R]"""

    throughput: float
    """Throughput (total customer arrival rate)"""

    pi: np.ndarray
    """Aggregated stationary probabilities [pi0, pi1, piStar]"""

    G: np.ndarray
    """G matrix"""

    mean_batch_size: float
    """Mean batch size of arrivals"""


def solver_mam_bmap_map_1(D: List[np.ndarray], S0: np.ndarray,
                          S1: np.ndarray) -> BMAPMAP1Result:
    """
    Solve a BMAP/MAP/1 queue using M/G/1-type matrix-analytic methods.

    The BMAP is specified by matrices {D0, D1, D2, ..., DK} where:
    - D0: transitions without arrivals
    - Dk: transitions triggering batch size k arrivals (k >= 1)

    The MAP for service is specified by matrices (S0, S1) where:
    - S0: transitions without service completions
    - S1: transitions triggering service completions

    The M/G/1-type structure for BMAP/MAP/1 is:
    ```
        A0 = I_ma otimes S1           (service completion, level -1)
        A1 = D0 otimes I_ms + I_ma otimes S0  (phase changes, level 0)
        A_{k+1} = D_k otimes I_ms     (batch arrival size k, level +k)
    ```

    Args:
        D: BMAP matrices as list [D0, D1, D2, ..., DK]
        S0: MAP service matrix for transitions without service completions
        S1: MAP service matrix for service completions

    Returns:
        BMAPMAP1Result with performance metrics
    """
    if len(D) < 1:
        raise ValueError("BMAP must have at least D0 matrix")
    if len(D) < 2:
        raise ValueError("BMAP must have at least D0 and D1 matrices")

    D = [np.asarray(Dk, dtype=float) for Dk in D]
    S0 = np.asarray(S0, dtype=float)
    S1 = np.asarray(S1, dtype=float)

    K = len(D) - 1  # Maximum batch size
    ma = D[0].shape[0]  # Number of BMAP phases
    ms = S0.shape[0]  # Number of MAP service phases
    m = ma * ms  # Combined phases per level

    # Validate dimensions
    for i, Dk in enumerate(D):
        if Dk.shape != (ma, ma):
            raise ValueError(f"All BMAP matrices must be {ma}x{ma}")
    if S0.shape != (ms, ms):
        raise ValueError(f"S0 must be {ms}x{ms}")
    if S1.shape != (ms, ms):
        raise ValueError(f"S1 must be {ms}x{ms}")

    # Compute arrival rate from BMAP
    D0 = D[0]
    D1_total = np.zeros((ma, ma))
    for k in range(1, K + 1):
        D1_total = D1_total + D[k]

    pi_bmap = map_piq(D0, D1_total)
    e_ma = np.ones(ma)

    # Total customer arrival rate: sum_k (k * pi * Dk * e)
    lambda_total = 0.0
    batch_rate = 0.0
    for k in range(1, K + 1):
        rate_k = pi_bmap @ D[k] @ e_ma
        lambda_total += k * rate_k
        batch_rate += rate_k

    # Mean batch size
    mean_batch_size = lambda_total / batch_rate if batch_rate > 0 else 0.0

    # Compute service rate from MAP
    pi_map = map_piq(S0, S1)
    e_ms = np.ones(ms)
    mu = pi_map @ S1 @ e_ms

    # Utilization
    rho = lambda_total / mu

    if rho >= 1.0:
        print(f"Warning: System is unstable (rho = {rho} >= 1). Results may be invalid.",
              file=sys.stderr)

    # Construct M/G/1-type matrices using Kronecker products
    I_ma = np.eye(ma)
    I_ms = np.eye(ms)

    # A = [A0, A1, A2, ..., A_{K+1}] for levels >= 1
    A = np.zeros((m, m * (K + 2)))

    # A0 = I_ma \otimes S1 (service completion)
    A0 = np.kron(I_ma, S1)
    A[:, 0:m] = A0

    # A1 = D0 \otimes I_ms + I_ma \otimes S0 (phase changes only)
    A1 = np.kron(D0, I_ms) + np.kron(I_ma, S0)
    A[:, m:2*m] = A1

    # A_{k+1} = D_k \otimes I_ms for k >= 1 (batch arrivals)
    for k in range(1, K + 1):
        Ak = np.kron(D[k], I_ms)
        A[:, (k+1)*m:(k+2)*m] = Ak

    # B = [B0, B1, B2, ..., BK] for level 0 (empty queue)
    B = np.zeros((m, m * (K + 1)))

    # B0: At level 0, no service (add S1 back as self-loop)
    B0 = np.kron(D0, I_ms) + np.kron(I_ma, S0 + S1)
    B[:, 0:m] = B0

    # B_k = D_k \otimes I_ms for k >= 1 (batch arrivals from empty queue)
    for k in range(1, K + 1):
        Bk = np.kron(D[k], I_ms)
        B[:, k*m:(k+1)*m] = Bk

    # Compute G matrix using ETAQA
    G = mg1_g_etaqa(A)

    # Compute stationary probabilities using ETAQA
    pi = mg1_pi_etaqa(B, A, G)

    # Compute mean queue length
    mean_queue_length = mg1_qlen_etaqa(B, A, pi, 1)

    # Performance metrics
    mean_response_time = mean_queue_length / lambda_total if lambda_total > 0 else 0.0

    return BMAPMAP1Result(
        mean_queue_length=mean_queue_length,
        utilization=rho,
        mean_response_time=mean_response_time,
        throughput=lambda_total,
        pi=pi,
        G=G,
        mean_batch_size=mean_batch_size
    )


__all__ = [
    'BMAPMAP1Result',
    'solver_mam_bmap_map_1',
]
