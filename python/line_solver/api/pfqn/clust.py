"""
de Souza e Silva-Lavenberg-Muntz Clustering Approximation (CA), de Souza e
Silva, Lavenberg and Muntz (1986).

Demand matrices follow the LINE convention L[station, class]; the cited papers
index them the other way round as D_ck.
"""

from typing import Tuple

import numpy as np

from .pam import pfqn_pam
from .utils import _amva_prep

__all__ = ['pfqn_clust']


def pfqn_clust(L, N, Z=None, subnets=None, localclasses=None, inner: str = 'lin',
               tol: float = 1e-6, maxiter: int = 1000
               ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    """de Souza e Silva-Lavenberg-Muntz Clustering Approximation (CA).

    E. de Souza e Silva, S. S. Lavenberg, R. R. Muntz, "A clustering
    approximation technique for queueing network models with a large number of
    chains", IEEE Trans. Computers C-35(5), 1986. The network is covered by
    subnetworks whose union is the whole network but which need not be
    disjoint. Every class visiting a subnetwork S is either LOCAL to S, and is
    then solved inside it, or FOREIGN, and is then seen only through the
    utilization it leaves behind. Each subnetwork is solved by an ordinary
    approximate MVA algorithm with two replacements: the complement of S is
    collapsed into a per-class delay P_c and the foreign classes into a
    per-centre utilization U_k,

        X_c(N) = N_c / (sum_{k in S} R_ck(N) + Z_c + P_c),
        Q_k(N) = [sum_{c in LC(S)} R_ck(N) X_c(N) + U_k] / (1 - U_k).

    Choosing the PE algorithm for every subnetwork reproduces global PE
    exactly, so the useful setting is Linearizer inside, PE outside.

    When no decomposition is supplied the criterion of the paper is applied
    automatically: the cheap PAMB estimate of the centre utilizations is taken,
    every class is attached to the centre where it loads the most, classes
    sharing that centre form one cluster, and the subnetwork of a cluster is
    the set of centres its classes visit.

    Returns (XN, QN, UN, RN, it).
    """
    L, N, Z, M, R = _amva_prep(L, N, Z)

    X0, Q0, _, _ = pfqn_pam(L, N, Z, 'pamb')
    XN = X0.flatten().copy()
    QN = Q0.copy()

    if not subnets or not localclasses:
        U0 = L * np.tile(XN, (M, 1))
        bottleneck = np.zeros(R, dtype=int)
        for r in range(R):
            if np.any(L[:, r] > 0):
                bottleneck[r] = int(np.argmax(U0[:, r]))
        subnets = []
        localclasses = []
        for centre in sorted(set(bottleneck.tolist())):
            cls = [r for r in range(R) if bottleneck[r] == centre]
            localclasses.append(cls)
            st = [i for i in range(M) if np.any(L[i, cls] > 0)]
            subnets.append(st if st else [centre])
        covered = set()
        for st in subnets:
            covered.update(st)
        missing = [i for i in range(M) if i not in covered]
        if missing:
            subnets.append(missing)
            localclasses.append([])

    owner = [-1] * R
    for g, lc in enumerate(localclasses):
        for c in lc:
            owner[c] = g

    it = 1
    for it in range(1, maxiter + 1):
        QN_old = QN.copy()
        Qk = QN.sum(axis=1)
        for g, S in enumerate(subnets):
            LC = localclasses[g]
            if not LC or not S:
                continue
            outside = [i for i in range(M) if i not in set(S)]
            FC = [r for r in range(R)
                  if r not in set(LC) and np.any(L[np.ix_(S, [r])] > 0)]
            # per-class delay in the complement of S
            P = np.zeros(len(LC))
            for a, c in enumerate(LC):
                if N[c] <= 0:
                    continue
                for ist in outside:
                    P[a] += L[ist, c] * (1 + Qk[ist]) / (1 + L[ist, c] * XN[c] / N[c])
            # utilization left in S by the foreign classes
            Uk = np.zeros(len(S))
            for b, ist in enumerate(S):
                for i in FC:
                    if N[i] > 0:
                        Uk[b] += L[ist, i] * XN[i] / (1 + L[ist, i] * XN[i] / N[i])
            Uk = np.minimum(Uk, 1 - 1e-8)

            Xs, Qs = _subnet_solve(L[np.ix_(S, LC)], N[LC], Z[LC] + P, Uk, inner, tol, maxiter)
            XN[LC] = Xs
            QN[np.ix_(S, LC)] = Qs
            # the local classes still hold jobs outside S
            for c in LC:
                nc = max(N[c], np.finfo(float).eps)
                for ist in outside:
                    QN[ist, c] = XN[c] * L[ist, c] * (1 + Qk[ist]) / (1 + L[ist, c] * XN[c] / nc)
        for r in range(R):
            if owner[r] < 0:
                QN[:, r] = XN[r] * L[:, r]
        if np.max(np.abs(QN - QN_old)) < tol:
            break

    UN = np.tile(XN, (M, 1)) * L
    with np.errstate(divide='ignore', invalid='ignore'):
        RN = QN / np.tile(XN, (M, 1))
    RN[:, N == 0] = 0.0
    RN = np.nan_to_num(RN, nan=0.0, posinf=0.0, neginf=0.0)
    return XN.reshape(1, -1), QN, UN, RN, it


def _subnet_solve(L, N, Z, Uk, inner, tol, maxiter):
    """Approximate MVA restricted to the local classes of one subnetwork."""
    M, R = L.shape
    X = np.zeros(R)
    Q = np.tile(N, (M, 1)) / max(M, 1)
    if M == 0 or R == 0:
        return X, Q
    if inner == 'bs':
        for _ in range(maxiter):
            Qold = Q.copy()
            for r in range(R):
                if N[r] <= 0:
                    continue
                # PE arrival-instant local queue, inflated by the foreign share
                A = (Q.sum(axis=1) - Q[:, r] / N[r] + Uk) / (1 - Uk)
                W = L[:, r] * (1 + A)
                X[r] = N[r] / (Z[r] + W.sum())
                Q[:, r] = X[r] * W
            if np.max(np.abs(Q - Qold)) < tol:
                break
        return X, Q

    Qs = [Q.copy() for _ in range(R + 1)]
    Delta = np.zeros((M, R, R))
    for _ in range(3):
        for s in range(R + 1):
            Ns = N.copy()
            if s > 0:
                Ns[s - 1] -= 1
            Qs[s] = _subnet_core(L, Ns, Z, Uk, Qs[s], Delta, tol, maxiter)
        for r in range(R):
            for s in range(1, R + 1):
                ns = N[r] - (1.0 if r == s - 1 else 0.0)
                if N[r] > 0 and ns > 0:
                    Delta[:, r, s - 1] = Qs[s][:, r] / ns - Qs[0][:, r] / N[r]
                elif N[r] > 0:
                    Delta[:, r, s - 1] = -Qs[0][:, r] / N[r]
                else:
                    Delta[:, r, s - 1] = 0.0
    Qs[0] = _subnet_core(L, N, Z, Uk, Qs[0], Delta, tol, maxiter)
    X, Q = _subnet_forward(L, N, Z, Uk, Qs[0], Delta)
    return X, Q


def _subnet_core(L, N, Z, Uk, Q, Delta, tol, maxiter):
    for _ in range(maxiter):
        Qold = Q
        _, Q = _subnet_forward(L, N, Z, Uk, Q, Delta)
        if np.max(np.abs(Q - Qold)) < tol:
            break
    return Q


def _subnet_forward(L, N, Z, Uk, Q, Delta):
    M, R = L.shape
    Qout = np.zeros((M, R))
    X = np.zeros(R)
    for r in range(R):
        if N[r] <= 0:
            continue
        # Linearizer estimate of the local queue at N - 1_r, inflated by the
        # foreign share
        Qm = np.zeros(M)
        for s in range(R):
            ns = N[s] - (1.0 if s == r else 0.0)
            if N[s] > 0 and ns > 0:
                Qm += ns * (Q[:, s] / N[s] + Delta[:, s, r])
        A = (np.maximum(Qm, 0) + Uk) / (1 - Uk)
        W = L[:, r] * (1 + A)
        X[r] = N[r] / (Z[r] + W.sum())
        Qout[:, r] = X[r] * W
    return X, Qout
