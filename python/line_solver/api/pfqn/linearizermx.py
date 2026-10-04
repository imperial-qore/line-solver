"""
Linearizer for mixed open/closed queueing networks.

Native Python implementation of the mixed open/closed Linearizer algorithm.

For multiserver networks (max_servers > 1), uses pfqn_linearizerms
(Conway/De Souza-Muntz algorithm) to match MATLAB's pfqn_linearizermx.

References:
    MATLAB: matlab/src/api/pfqn/pfqn_linearizermx.m
"""

import numpy as np
from typing import Tuple, List, Optional

from .linearizer import pfqn_linearizer, pfqn_gflinearizer, pfqn_egflinearizer
from .linearizerms import pfqn_linearizerms


def pfqn_linearizermx(
    lambda_arr: np.ndarray,
    L: np.ndarray,
    N: np.ndarray,
    Z: np.ndarray,
    nservers: np.ndarray,
    sched_type: List[str],
    tol: float = 1e-8,
    maxiter: int = 1000,
    method: str = 'egflin',
    QN0: np.ndarray = None
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    """
    Linearizer for mixed open/closed queueing networks.

    This function extends the linearizer algorithm to handle networks with
    both open classes (with external arrivals) and closed classes (with
    fixed populations).

    Args:
        lambda_arr: Arrival rate vector (R,). For closed classes, should be 0 or inf.
        L: Service demand matrix (M x R)
        N: Population vector (R,). Inf for open classes, finite for closed.
        Z: Think time vector (R,) or matrix
        nservers: Number of servers per station (M,)
        sched_type: Scheduling strategy per station (list of strings)
        tol: Convergence tolerance
        maxiter: Maximum iterations
        method: Linearizer variant ('lin', 'gflin', 'egflin')

    Returns:
        QN: Mean queue lengths (M x R)
        UN: Utilization (M x R)
        WN: Waiting times (M x R)
        TN: Throughputs (M x R)
        CN: Cycle times (1 x R)
        XN: System throughput (R,)
        totiter: Total iterations

    References:
        MATLAB: matlab/src/api/pfqn/pfqn_linearizermx.m
    """
    L = np.atleast_2d(L).astype(float)
    N = np.atleast_1d(N).astype(float)
    lambda_arr = np.atleast_1d(lambda_arr).astype(float)

    M, R = L.shape

    if Z is None:
        Z = np.zeros(R)
    else:
        Z = np.atleast_1d(Z).astype(float)
        if Z.ndim > 1:
            Z = np.sum(Z, axis=0)

    # Replace NaN with 0
    lambda_arr = np.nan_to_num(lambda_arr, nan=0.0)
    L = np.nan_to_num(L, nan=0.0)
    Z = np.nan_to_num(Z, nan=0.0)

    # Identify open and closed classes
    open_classes = np.where(np.isinf(N))[0]
    closed_classes = np.array([r for r in range(R) if r not in open_classes])

    # Initialize outputs
    XN = np.zeros(R)
    UN = np.zeros((M, R))
    WN = np.zeros((M, R))
    QN = np.zeros((M, R))
    CN = np.zeros(R)
    TN = np.zeros((M, R))

    # Process open classes: U = lambda * L, X = lambda
    # Reference: MATLAB pfqn_linearizermx.m lines 62-67
    nservers_arr = np.atleast_1d(nservers).astype(float) if nservers is not None else np.ones(M)
    if len(nservers_arr) < M:
        nservers_arr = np.concatenate([nservers_arr, np.ones(M - len(nservers_arr))])

    for r in open_classes:
        for ist in range(M):
            UN[ist, r] = lambda_arr[r] * L[ist, r]
        XN[r] = lambda_arr[r]

    # Total utilization from open classes
    UNt = np.sum(UN, axis=1)

    # see _kb/03-api-layer.md for rationale
    if len(closed_classes) > 0:
        Dc = np.zeros((M, len(closed_classes)))
        for idx, r in enumerate(closed_classes):
            Dc[:, idx] = L[:, r] / (1.0 - UNt)

        N_closed = N[closed_classes]
        Z_closed = Z[closed_classes]

        # initial closed-class queue lengths (warm start) aligned to the algorithm layout
        if QN0 is None:
            QN0c = None
        else:
            QN0a = np.asarray(QN0, dtype=float)
            if QN0a.shape == (M, len(N)):
                QN0c = QN0a[:, closed_classes]
            elif QN0a.shape == (M, len(closed_classes)):
                QN0c = QN0a
            else:
                QN0c = None

        # Call appropriate linearizer for closed classes. MATLAB tests
        # max(nservers(:))==1 without dropping the infinite entries, so an
        # INF-scheduled station routes to pfqn_linearizerms; the single-server
        # linearizers ignore `type` and would queue at it.
        max_servers = np.max(nservers_arr) if nservers_arr.size > 0 else 1

        # Reference: MATLAB pfqn_linearizermx.m lines 78-96
        if max_servers == 1:
            if method == 'gflin':
                lin_alpha = 2.0
                result = pfqn_gflinearizer(Dc, N_closed, Z_closed, sched_type, tol, maxiter, lin_alpha, QN0=QN0c)
            elif method == 'egflin':
                # MATLAB: alphaM = zeros(1,R); for r=closedClasses, alphaM(r) = ...
                # MATLAB passes full R-sized vector but only first len(closedClasses) elements are used
                alphaM = np.zeros(len(closed_classes))
                for idx, r in enumerate(closed_classes):
                    alphaM[idx] = 0.6 + 1.4 * np.exp(-8.0 * np.exp(-0.8 * N[r]))
                result = pfqn_egflinearizer(Dc, N_closed, Z_closed, sched_type, tol, maxiter, alphaM, QN0=QN0c)
            else:  # 'lin' or default
                result = pfqn_linearizer(Dc, N_closed, Z_closed, sched_type, tol, maxiter, QN0=QN0c)
            QNc, UNc, WNc, TNc, CNc, XNc, totiter = result
        else:
            QNc, UNc, WNc, CNc, XNc, totiter = pfqn_linearizerms(
                Dc, N_closed, Z_closed, nservers_arr, sched_type, tol, maxiter, QN0=QN0c
            )

        # Map closed class results back
        # Reference: MATLAB pfqn_linearizermx.m lines 98-102
        XNc_flat = XNc.flatten() if XNc.ndim > 1 else XNc
        CNc_flat = CNc.flatten() if hasattr(CNc, 'flatten') else CNc

        for idx, r in enumerate(closed_classes):
            XN[r] = XNc_flat[idx] if idx < len(XNc_flat) else 0.0
            CN[r] = CNc_flat[idx] if idx < len(CNc_flat) else 0.0
            for ist in range(M):
                if idx < QNc.shape[1]:
                    QN[ist, r] = QNc[ist, idx]
                    WN[ist, r] = WNc[ist, idx]

        # Recompute closed class utilization using original demands (not Dc)
        # Reference: MATLAB pfqn_linearizermx.m lines 104-108
        for ist in range(M):
            for idx, r in enumerate(closed_classes):
                UN[ist, r] = XN[r] * L[ist, r]
    else:
        QNc = np.empty((M, 0))
        totiter = 0

    # Compute open class response times using closed class queue lengths
    # Reference: MATLAB pfqn_linearizermx.m lines 109-118
    for ist in range(M):
        for r in open_classes:
            if QNc.size == 0:
                WN[ist, r] = L[ist, r] / (1.0 - UNt[ist])
            else:
                WN[ist, r] = L[ist, r] * (1.0 + np.sum(QNc[ist, :])) / (1.0 - UNt[ist])
            QN[ist, r] = WN[ist, r] * XN[r]

    # Cycle times for open classes
    # Reference: MATLAB pfqn_linearizermx.m line 119
    CN[open_classes] = np.sum(WN[:, open_classes], axis=0)

    return QN, UN, WN, TN, CN, XN, totiter
