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

from line_solver.lib.thirdparty.smc import (
    gim1_r_etaqa, gim1_pi_etaqa, gim1_qlen_etaqa)
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

    The GI/M/1-type structure for MAP/BMAP/1 is::

          B1  A0  0   0   ...
          B2  A1  A0  0   ...
      Q = B3  A2  A1  A0  ...
          ...

    Where::

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

    # Compute R matrix using GI/M/1 ETAQA (GIM1_R_ETAQA -> GIM1_R, dual + FI)
    R = gim1_r_etaqa(A)

    # Compute stationary probabilities using GI/M/1 ETAQA
    pi = gim1_pi_etaqa(B, A, R, A0)

    # Compute mean queue length. The reference asks for the FIRST moment only
    # on this side, and reproduces GIM1_qlen_ETAQA's scalar-A(3) defect, which
    # can make the reported mean negative for more than one phase.
    mean_queue_length = gim1_qlen_etaqa(B, A, R, pi, 1, A0)

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


__all__ = [
    'MAPBMAP1Result',
    'solver_mam_map_bmap_1',
]
