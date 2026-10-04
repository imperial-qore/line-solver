"""
Hsieh-Lam Proportional Approximation Methods (PAMB/PAMI/PAMT), Hsieh and
Lam (1988).

Demand matrices follow the LINE convention L[station, class]; the cited papers
index them the other way round as D_ck.
"""

from typing import Tuple

import numpy as np

from .utils import _amva_prep

__all__ = ['pfqn_pam']


def pfqn_pam(L, N, Z=None, variant: str = 'pamb'
             ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Hsieh-Lam Proportional Approximation Methods (PAMB/PAMI/PAMT).

    C. T. Hsieh, S. S. Lam, "PAM - A noniterative approximate solution method
    for closed multichain queueing networks", ACM SIGMETRICS Perform. Eval.
    Rev. 16(1), 1988. The three variants are NONITERATIVE: the queue lengths
    are seeded by the proportion of a class demand that falls at each centre,

        E_ck = D_ck / sum_i D_ci,   Q_ck(N) = E_ck N_c,

    and the MVA equations are then unrolled a fixed number of times. 'pamb'
    applies the last MVA step; 'pami' additionally scales a class down wherever
    it would drive a centre past full utilization; 'pamt' seeds at
    N - 1_i - 1_j and applies the last TWO MVA steps before that capping. The
    seed spreads the whole class population over the queueing centres and
    ignores Z, exactly as published: PAM buys speed, not accuracy.

    Returns (XN, QN, UN, RN).
    """
    L, N, Z, M, R = _amva_prep(L, N, Z)
    variant = (variant or 'pamb').lower()

    E = np.zeros((M, R))
    Ltot = L.sum(axis=0)
    for r in range(R):
        if Ltot[r] > 0:
            E[:, r] = L[:, r] / Ltot[r]
    Q = E * np.tile(N, (M, 1))

    RN = np.zeros((M, R))
    XN = np.zeros(R)
    if variant == 'pamt':
        for i in range(R):
            Qmi = np.zeros((M, R))
            for j in range(R):
                # Q_ck(N - 1_i - 1_j) = Q_ck(N) - E_ck [(c==i) + (c==j)]
                Qij = Q.copy()
                Qij[:, i] -= E[:, i]
                Qij[:, j] -= E[:, j]
                Rj = L[:, j] * (1 + Qij.sum(axis=1))
                nj = N[j] - (1.0 if i == j else 0.0)
                Xj = nj / (Rj.sum() + Z[j]) if nj > 0 else 0.0
                Qmi[:, j] = Xj * Rj
            RN[:, i] = L[:, i] * (1 + Qmi.sum(axis=1))
            if N[i] > 0:
                XN[i] = N[i] / (RN[:, i].sum() + Z[i])
    else:
        for r in range(R):
            # Q_jk(N - 1_r) = Q_jk(N) - E_jk [j == r]
            Qmr = Q.copy()
            Qmr[:, r] -= E[:, r]
            RN[:, r] = L[:, r] * (1 + Qmr.sum(axis=1))
            if N[r] > 0:
                XN[r] = N[r] / (RN[:, r].sum() + Z[r])

    if variant in ('pami', 'pamt'):
        # scale a class down when it would drive a centre it visits past U = 1
        U = (L * np.tile(XN, (M, 1))).sum(axis=1)
        for r in range(R):
            visited = L[:, r] != 0
            if np.any(visited):
                S = float(np.max(U[visited]))
                if S > 1:
                    XN[r] = XN[r] / S

    QN = np.tile(XN, (M, 1)) * RN
    UN = np.tile(XN, (M, 1)) * L
    return XN.reshape(1, -1), QN, UN, RN
