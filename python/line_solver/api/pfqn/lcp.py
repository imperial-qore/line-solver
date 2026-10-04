"""
Bard Large Customer Population (LCP) approximate MVA, Bard (1979).

Demand matrices follow the LINE convention L[station, class]; the cited papers
index them the other way round as D_ck.
"""

from typing import Tuple

import numpy as np

from .utils import _amva_is_fcfs, _amva_prep

__all__ = ['pfqn_lcp']


def pfqn_lcp(L, N, Z=None, tol: float = 1e-6, maxiter: int = 1000,
             QN0=None, type_sched=None
             ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    """Bard Large Customer Population (LCP) approximate MVA.

    Y. Bard, "Some extensions to multiclass queueing network analysis", in
    Performance of Computer Systems, North-Holland, 1979. The arrival-instant
    queue length is estimated by the time-averaged one WITHOUT removing the
    arriving customer,

        A_k^(c)(N) = Q_k(N - 1_c) ~= Q_k(N) = sum_s Q_ks(N),

    since with a large population one customer less cannot change the mean
    queue lengths appreciably. Setting the Bard-Schweitzer proportional term
    Q_kc(N)/N_c to zero recovers this algorithm, so LCP is uniformly more
    pessimistic than pfqn_bs and is inaccurate at small populations.

    Returns (XN, QN, UN, RN, it), matching pfqn_bs.
    """
    from ...lang.base import SchedStrategy

    L, N, Z, M, R = _amva_prep(L, N, Z)
    QN = np.tile(N, (M, 1)) / M if QN0 is None else np.asarray(QN0, dtype=np.float64).copy()
    if type_sched is None:
        type_sched = [SchedStrategy.PS] * M

    CN = np.zeros((M, R))
    XN = np.zeros(R)
    UN = np.zeros((M, R))
    it = 1
    for it in range(1, maxiter + 1):
        QN_old = QN.copy()
        for r in range(R):
            if N[r] == 0:
                XN[r] = 0.0
                CN[:, r] = 0.0
                QN[:, r] = 0.0
                UN[:, r] = 0.0
                continue
            for ist in range(M):
                CN[ist, r] = L[ist, r]
                if L[ist, r] == 0:
                    continue
                for s in range(R):
                    # the arriving customer is NOT removed: no (N-1)/N factor
                    if s != r and _amva_is_fcfs(type_sched[ist]):
                        CN[ist, r] += L[ist, s] * QN[ist, s]
                    else:
                        CN[ist, r] += L[ist, r] * QN[ist, s]
            CN_sum = float(np.sum(CN[:, r]))
            XN[r] = N[r] / (Z[r] + CN_sum) if (Z[r] + CN_sum) > 0 else 0.0
        for r in range(R):
            QN[:, r] = XN[r] * CN[:, r]
            UN[:, r] = XN[r] * L[:, r]
        with np.errstate(divide='ignore', invalid='ignore'):
            rel = np.abs(1 - QN / QN_old)
            rel = np.nan_to_num(rel, nan=0.0, posinf=0.0, neginf=0.0)
        if np.max(rel) < tol:
            break

    RN = np.zeros((M, R))
    for r in range(R):
        if XN[r] > 0:
            RN[:, r] = QN[:, r] / XN[r]
    return XN.reshape(1, -1), QN, UN, RN, it
