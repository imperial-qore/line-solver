"""Neuse-Chandy SCAT (Self-Correcting Approximation Technique) approximate MVA.

SCAT shares the Linearizer fixed point: it carries the mean queue lengths at the
target population N and at the R reduced populations N-e_s, and corrects the
Bard-Schweitzer proportionality assumption with the fraction difference

    Delta(i,r,s) = Q(i,r|N-e_s)/(N-e_s)_r - Q(i,r|N)/N_r,

held fixed while an inner MVA fixed point is iterated. It differs from
Linearizer in that this correction is refreshed ONCE: SCAT stops after the first
pass, where Linearizer performs the fixed three passes of Chandy and Neuse
(1982), Sec. 4. Cost is therefore about one third of Linearizer's, and accuracy
sits between Bard-Schweitzer (the Delta=0 special case) and Linearizer.

SCAT's second departure from Linearizer, fitting a probability mass function
centred on the mean queue length at queue-dependent centres instead of
propagating the MVA distribution recursion (Krzesinski and Greyling 1984,
Sec. 4), does not arise here: this entry point covers single-server and delay
stations only, exactly as pfqn_linearizer does. That mass function is available
separately as the 'scat' marginal rule of pfqn_ab_amva.

Reference: D. Neuse, K. M. Chandy, "SCAT: A Heuristic Algorithm for Queueing
Network Models of Computing Systems", ACM SIGMETRICS Perform. Eval. Rev. 10(3),
1981.
"""

import numpy as np
from typing import List, Optional, Tuple

from .linearizer import pfqn_egflinearizer

__all__ = ['pfqn_scat']


def pfqn_scat(
    L: np.ndarray,
    N: np.ndarray,
    Z: np.ndarray,
    sched_type: Optional[List[str]] = None,
    tol: float = 1e-8,
    maxiter: int = 1000,
    QN0: np.ndarray = None
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    """SCAT approximate MVA.

    Args:
        L: Demand matrix (M x R)
        N: Population vector (R,)
        Z: Think time vector (R,)
        sched_type: Scheduling strategy per station; accepted for interface
            parity with pfqn_linearizer, but the residence-time recursion is
            discipline-independent
        tol: Convergence tolerance
        maxiter: Maximum inner iterations
        QN0: (M x R) warm start for the Bard-Schweitzer initialization

    Returns:
        Same tuple as pfqn_linearizer: (Q, U, W, T, C, X, iterations), where T
        is the throughput PER REFERENCE VISIT.
    """
    L = np.atleast_2d(np.asarray(L, dtype=float))
    N = np.asarray(N, dtype=float).ravel()
    Z = np.asarray(Z, dtype=float).ravel()
    alpha = np.ones(len(N))
    # npasses=1 is what separates SCAT from Linearizer: one Delta refresh, not three
    return pfqn_egflinearizer(L, N, Z, sched_type, tol, maxiter, alpha,
                              QN0=QN0, npasses=1)
