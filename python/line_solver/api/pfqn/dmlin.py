"""
de Souza e Silva-Muntz Improved Linearizer (IL), de Souza e Silva and
Muntz (1990).

Demand matrices follow the LINE convention L[station, class]; the cited papers
index them the other way round as D_ck.
"""

from typing import Tuple

import numpy as np

from .mva import pfqn_bs
from .utils import _amva_prep

__all__ = ['pfqn_dmlin']


def pfqn_dmlin(L, N, Z=None, type_sched=None, tol: float = 1e-8,
               maxiter: int = 1000, QN0=None, npasses: int = 3
               ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray,
                          np.ndarray, int]:
    """de Souza e Silva-Muntz Improved Linearizer (IL).

    E. de Souza e Silva, R. R. Muntz, "A note on the computational cost of the
    Linearizer algorithm for queueing networks", IEEE Trans. Computers 39(6),
    1990. Linearizer evaluates the arrival-instant queue length as

        A_k^(c)(n) = sum_i (n_i - delta_c^(i)) [Q_ik(n)/n_i + Delta^(i)_ck],

    re-summing the C Delta-terms at every Core iteration, at every one of the
    C+1 populations: O(K C^3) per refresh pass. IL splits that sum into the
    part that moves with the Core iterate and the part that does not,

        A_k^(c)(n)     = sum_i (n_i - delta_c^(i)) Q_ik(n)/n_i + xi_ck(n),
        xi_ck(N)       = sum_i (N_i - delta_c^(i)) Delta^(i)_ck,
        xi_ck(N - 1_j) = xi_ck(N) - Delta^(j)_ck,

    so the C K aggregates xi are computed ONCE per refresh pass and each Core
    iteration then costs O(K C) instead of O(K C^2). Because the split is an
    identity and not an approximation, the fixed point is the one Linearizer
    reaches: pfqn_dmlin and pfqn_linearizer agree to round-off.

    ``type_sched`` is accepted for signature parity and unused: the Linearizer
    family in LINE treats every station as single-server PS.

    Returns (Q, U, W, T, C, X, totiter), matching pfqn_linearizer.
    """
    del type_sched
    L, N, Z, M, R = _amva_prep(L, N, Z)

    if M == 0 or not np.any(L):
        X = np.where(Z > 0, N / np.where(Z > 0, Z, 1), 0.0)
        U = np.tile(X, (M, 1)) * L
        return (np.zeros((M, R)), U, np.zeros((M, R)), np.zeros((M, R)),
                np.zeros(R), X.reshape(1, -1), 0)

    def _oner(v, s):
        out = v.copy()
        if s > 0:
            out[s - 1] -= 1
        return out

    def _core(N1, Qin, xi, budget):
        Q = Qin.copy()
        W = L.copy()
        T = np.zeros(R)
        iters = 0
        while True:
            Qlast = Q.copy()
            A = np.zeros((M, R))
            for c in range(R):
                acc = np.zeros(M)
                for r in range(R):
                    if N1[r] > 0:
                        nr = N1[r] - (1.0 if r == c else 0.0)
                        if nr > 0:
                            acc += nr * Q[:, r] / N1[r]
                A[:, c] = acc + xi[:, c]
            W = L * (1 + A)
            for r in range(R):
                T[r] = N1[r] / (Z[r] + W[:, r].sum()) if N1[r] > 0 else 0.0
                Q[:, r] = T[r] * W[:, r]
            converged = np.linalg.norm(Q - Qlast) < tol or iters > budget
            iters += 1
            if converged:
                break
        return Q, W, T, iters

    # Initialize, as Linearizer does, from Bard-Schweitzer at every population
    Qs = [None] * (R + 1)
    for s in range(R + 1):
        N1 = _oner(N, s)
        if QN0 is None:
            _, q, _, _, _ = pfqn_bs(L, N1, Z)
        else:
            _, q, _, _, _ = pfqn_bs(L, N1, Z, tol, maxiter, np.asarray(QN0, dtype=np.float64).copy())
        Qs[s] = q.copy()

    Delta = np.zeros((M, R, R))     # Delta[i, r, c] = Delta^(r)_c at station i
    xi = np.zeros((M, R))

    totiter = 0
    for _ in range(npasses):
        for s in range(R + 1):
            N1 = _oner(N, s)
            # xi at population N - 1_s, exactly; s == 0 leaves xi at N
            xis = xi if s == 0 else xi - Delta[:, s - 1, :]
            Qs[s], _, _, iters = _core(N1, Qs[s], xis, maxiter - totiter)
            totiter += iters
        for r in range(R):
            if N[r] == 1:
                Qs[r + 1][:, r] = 0.0
            for s in range(1, R + 1):
                ns = N[r] - (1.0 if r == s - 1 else 0.0)
                if N[r] > 0 and ns > 0:
                    Delta[:, r, s - 1] = Qs[s][:, r] / ns - Qs[0][:, r] / N[r]
                elif N[r] > 0:
                    Delta[:, r, s - 1] = -Qs[0][:, r] / N[r]
                else:
                    Delta[:, r, s - 1] = 0.0
        for c in range(R):
            acc = np.zeros(M)
            for r in range(R):
                w = N[r] - (1.0 if r == c else 0.0)
                if w > 0:
                    acc += w * Delta[:, r, c]
            xi[:, c] = acc

    Q, W, X, iters = _core(N, Qs[0], xi, maxiter - totiter)
    totiter += iters
    U = np.tile(X, (M, 1)) * L
    T = np.tile(X, (M, 1))
    with np.errstate(divide='ignore', invalid='ignore'):
        C = np.where(X > 0, N / np.where(X > 0, X, 1) - Z, 0.0)
    return Q, U, W, T, C.reshape(1, -1), X.reshape(1, -1), totiter
