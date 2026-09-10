"""
Chow Second Approximation (SA) approximate MVA, Chow (1983).

Demand matrices follow the LINE convention L[station, class]; the cited papers
index them the other way round as D_ck.
"""

from typing import Tuple

import numpy as np

from .lcp import pfqn_lcp
from .utils import _amva_is_fcfs, _amva_prep

__all__ = ['pfqn_chow']


def pfqn_chow(L, N, Z=None, tol: float = 1e-6, maxiter: int = 1000,
              QN0=None, type_sched=None, variant: str = 'forward'
              ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    """Chow Second Approximation (SA) approximate MVA.

    W.-M. Chow, "Approximations for large scale closed queueing networks",
    Perform. Eval. 3(1), 1983. The arrival-instant queue length is written
    exactly as

        A_k^(c)(N) = Q_k(N - 1_c) = Q_k(N) (1 + theta_ck),
        theta_ck   = [Q_k(N - 1_c) - Q_k(N)] / Q_k(N),

    and the theta-terms are estimated ONCE, off the Bard LCP solution, before
    the fixed point is run. ``variant='forward'`` uses Qhat(N + 1_c) and
    ``'backward'`` uses Qhat(N - 1_c); Chow reports the forward form to be the
    more accurate, so it is the default. Setting every theta to zero recovers
    pfqn_lcp.

    Returns (XN, QN, UN, RN, it).
    """
    from ...lang.base import SchedStrategy

    L, N, Z, M, R = _amva_prep(L, N, Z)
    if type_sched is None:
        type_sched = [SchedStrategy.PS] * M

    _, Qlcp, _, _, _ = pfqn_lcp(L, N, Z, tol, maxiter, QN0, type_sched)
    Qtot = Qlcp.sum(axis=1)
    theta = np.zeros((M, R))
    for r in range(R):
        if N[r] == 0:
            continue
        Nalt = N.copy()
        if variant == 'backward':
            Nalt[r] -= 1
            _, Qalt, _, _, _ = pfqn_lcp(L, Nalt, Z, tol, maxiter, QN0, type_sched)
            base = Qtot
            delta = Qalt.sum(axis=1) - Qtot
        else:
            Nalt[r] += 1
            _, Qalt, _, _, _ = pfqn_lcp(L, Nalt, Z, tol, maxiter, QN0, type_sched)
            base = Qalt.sum(axis=1)
            delta = Qtot - base
        nz = base > 0
        theta[nz, r] = delta[nz] / base[nz]

    QN = np.tile(N, (M, 1)) / M if QN0 is None else np.asarray(QN0, dtype=np.float64).copy()
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
                    if s != r and _amva_is_fcfs(type_sched[ist]):
                        CN[ist, r] += L[ist, s] * QN[ist, s] * (1 + theta[ist, r])
                    else:
                        CN[ist, r] += L[ist, r] * QN[ist, s] * (1 + theta[ist, r])
                # a theta below -1 would make the arrival-instant queue negative
                CN[ist, r] = max(CN[ist, r], L[ist, r])
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
