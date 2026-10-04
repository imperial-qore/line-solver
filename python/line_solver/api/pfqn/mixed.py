"""
Mixed MVA for Product-Form Queueing Networks.

Implements MVA algorithms for networks with both open and closed classes.
"""

import numpy as np
from typing import Tuple, Optional


def pfqn_mvamx(
    lambda_arr: np.ndarray,
    L: np.ndarray,
    N: np.ndarray,
    Z: np.ndarray,
    mi: Optional[np.ndarray] = None
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float]:
    """
    Exact MVA for mixed open/closed single-server networks.

    Handles networks with both open classes (infinite population, external arrivals)
    and closed classes (fixed population, no external arrivals).

    Args:
        lambda_arr: Arrival rate vector (R,). Use 0 for closed classes.
        L: Service demand matrix (M x R).
        N: Population vector (R,). Use np.inf for open classes.
        Z: Think time vector (R,).
        mi: Queue replication factors (M,) (default: ones).

    Returns:
        Tuple of (XN, QN, UN, CN, lGN):
            XN: Class throughputs (1 x R).
            QN: Mean queue lengths (M x R).
            UN: Utilizations (M x R).
            CN: Cycle times (M x R).
            lGN: Log normalizing constant.
    """
    from .mva import pfqn_mva

    L = np.atleast_2d(np.asarray(L, dtype=float))
    N = np.asarray(N, dtype=float).ravel()
    Z = np.asarray(Z, dtype=float).ravel()
    lambda_arr = np.asarray(lambda_arr, dtype=float).ravel()

    M, R = L.shape

    if mi is None:
        mi = np.ones(M)
    else:
        mi = np.asarray(mi, dtype=float).ravel()

    # Check that arrival rate is not specified on closed classes
    for r in range(R):
        if lambda_arr[r] > 0 and N[r] > 0 and np.isfinite(N[r]):
            raise ValueError("Arrival rate cannot be specified on closed classes.")

    # Identify open and closed classes
    open_classes = np.where(np.isinf(N))[0]
    closed_classes = np.array([r for r in range(R) if r not in open_classes])

    # Initialize outputs
    XN = np.zeros((1, R))
    UN = np.zeros((M, R))
    CN = np.zeros((M, R))
    QN = np.zeros((M, R))

    # Process open classes first
    for r in open_classes:
        for i in range(M):
            UN[i, r] = lambda_arr[r] * L[i, r]
        XN[0, r] = lambda_arr[r]

    # Total utilization from open classes
    UNt = np.sum(UN, axis=1)  # M x 1

    if len(Z) == 0:
        Z = np.zeros(R)

    # Process closed classes using adjusted demands
    if len(closed_classes) > 0:
        # Adjust demands for open class utilization
        Dc = np.zeros((M, len(closed_classes)))
        for j, c in enumerate(closed_classes):
            for i in range(M):
                if UNt[i] < 1:
                    Dc[i, j] = L[i, c] / (1 - UNt[i])
                else:
                    Dc[i, j] = L[i, c] * 1e10  # Very high demand if saturated

        N_closed = N[closed_classes].astype(int)
        Z_closed = Z[closed_classes]

        # Compute closed class metrics using standard MVA
        XNc, CNc, QNc, UNc, RNc, TNc, ANc = pfqn_mva(Dc, N_closed, Z_closed, mi)

        # Copy results back
        XN[0, closed_classes] = XNc.ravel()
        QN[:, closed_classes] = QNc
        CN[:, closed_classes] = RNc

        # Compute utilizations for closed classes
        for i in range(M):
            for c in closed_classes:
                UN[i, c] = XN[0, c] * L[i, c]

        # Compute open class response times and queue lengths
        # Use tolerance for saturation check to handle floating-point precision issues
        SAT_TOL = 1e-10
        for i in range(M):
            for r in open_classes:
                if UNt[i] < 1 - SAT_TOL:
                    if QNc.size > 0:
                        CN[i, r] = L[i, r] * (1 + np.sum(QNc[i, :])) / (1 - UNt[i])
                    else:
                        CN[i, r] = L[i, r] / (1 - UNt[i])
                else:
                    CN[i, r] = np.inf
                QN[i, r] = CN[i, r] * XN[0, r]

        # Compute log normalizing constant
        from .nc import pfqn_ca
        _, lGN = pfqn_ca(Dc, N_closed, Z_closed)
    else:
        lGN = 0.0
        # Only open classes - compute response times
        # Use tolerance for saturation check to handle floating-point precision issues
        SAT_TOL = 1e-10
        for i in range(M):
            for r in open_classes:
                if UNt[i] < 1 - SAT_TOL:
                    CN[i, r] = L[i, r] / (1 - UNt[i])
                else:
                    CN[i, r] = np.inf
                QN[i, r] = CN[i, r] * XN[0, r]

    return XN, QN, UN, CN, lGN


__all__ = [
    'pfqn_mvamx',
]
