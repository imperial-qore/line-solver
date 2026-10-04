# -*- coding: utf-8 -*-
"""
First- and second-order (Brownian) level-dependent, multi-regime Markovian
fluid queues.

Ported from BUTools-family fluid tools (G. Horvath / Kankaya-Akar):
  * SecondOrderLevelDependentFluidSolve  - matrix-exponential building blocks
  * LevelDependentFluidStationaryDistr   - stationary pdf/cdf
  * LevelDependentFluidStationaryMean    - stationary mean fluid level
  * multiregime                          - multi-regime feedback fluid queue

The generator, drift and (optionally) variance change at threshold fluid
levels, yielding a piecewise-homogeneous structure. Setting the variance
cells to zero reduces the model to first order.
"""
import numpy as np
import numpy.matlib as ml
import scipy.linalg as la
from scipy.linalg.lapack import get_lapack_funcs

from .qbd import QBDFundamentalMatrices



__all__ = ["SecondOrderLevelDependentFluidSolve", "LevelDependentFluidStationaryDistr", "LevelDependentFluidStationaryMean", "multiregime"]

def _schur_realparts(T):
    """Real part of each eigenvalue by diagonal position of a real Schur form."""
    n = T.shape[0]
    rp = np.zeros(n)
    i = 0
    while i < n:
        if i < n - 1 and abs(T[i + 1, i]) > 1e-12 * (abs(T[i, i]) + abs(T[i + 1, i + 1]) + 1e-300):
            rp[i] = rp[i + 1] = T[i, i]
            i += 2
        else:
            rp[i] = T[i, i]
            i += 1
    return rp


def _ordered_schur(A, tol=1e-10):
    """
    Real Schur decomposition of A reordered so that the eigenvalues appear
    grouped as [zero real-part, negative real-part, positive real-part], the
    ordering used by MATLAB's ordschur in the multi-regime solver.
    Returns (Z, D, zeroeig, negeig, poseig).
    """
    T, Z = la.schur(np.asarray(A, dtype=float), output='real')
    trsen = get_lapack_funcs('trsen', (T,))
    # pass 1: bring zero real-part eigenvalues to the top
    rp = _schur_realparts(T)
    sel1 = (np.abs(rp) < tol).astype(np.int32)
    T, Z = trsen(sel1, T, Z, job='N')[:2]
    # pass 2: bring {zero, negative} to the top -> [zero, neg, pos]
    rp = _schur_realparts(T)
    sel2 = ((np.abs(rp) < tol) | (rp < -tol)).astype(np.int32)
    T, Z = trsen(sel2, T, Z, job='N')[:2]
    rp = _schur_realparts(T)
    zeroeig = int(np.sum(np.abs(rp) < tol))
    negeig = int(np.sum(rp < -tol))
    poseig = int(np.sum(rp > tol))
    return Z, T, zeroeig, negeig, poseig


def _nullvec(KA):
    """Left null vector of KA (CRPSolve algorithm, without the generator check)."""
    M = np.array(KA, dtype=float)
    M[:, 0] = 1.0
    m = np.zeros(M.shape[0])
    m[0] = 1.0
    return la.solve(M.T, m)


def _integExp(KA, L):
    """int_0^L expm(KA u) du, robust when KA has a zero eigenvalue (deflation)."""
    KA = np.asarray(KA)
    n = KA.shape[0]
    l = _nullvec(KA)
    r = _nullvec(KA.T)
    l = l / (l @ r)
    rl = np.outer(r, l)
    return la.inv(-(KA - rl)) @ (np.eye(n) - la.expm((KA - rl) * L)) + rl * (L + np.exp(-L) - 1)


def _integExp2(KA, KB_, L):
    KA = np.asarray(KA)
    KB_ = np.asarray(KB_)
    eKA = np.min(np.abs(la.eigvals(KA))) if KA.shape[0] > 0 else np.inf
    eKB = np.min(np.abs(la.eigvals(KB_))) if KB_.shape[0] > 0 else np.inf
    if eKA > eKB:
        KAi = la.inv(-KA) @ (np.eye(KA.shape[0]) - la.expm(KA * L))
        KBi = _integExp(KB_, L)
    else:
        KAi = _integExp(KA, L)
        KBi = la.inv(-KB_) @ (np.eye(KB_.shape[0]) - la.expm(KB_ * L))
    return KAi, KBi


def _dg(M):
    return np.diag(np.asarray(M))


def _diag(v):
    return ml.matrix(np.diagflat(np.asarray(v).flatten()))


def _blkdiag(*mats):
    return ml.matrix(la.block_diag(*[np.asarray(m) for m in mats]))


def _sub(M, ri, ci):
    M = np.asarray(M)
    if len(ri) == 0 or len(ci) == 0:
        return ml.zeros((len(ri), len(ci)))
    return ml.matrix(M[np.ix_(ri, ci)])


def _safemax(*vals):
    """max over a set of arrays/scalars, ignoring empty; floor handled by caller."""
    acc = []
    for v in vals:
        a = np.asarray(v, dtype=float).flatten()
        if a.size > 0:
            acc.append(np.max(a))
    if not acc:
        return -np.inf
    return max(acc)


def SecondOrderLevelDependentFluidSolve(Q, R, S, T, boundaryL=None, boundaryU=None, Qt=None, prec=1e-14):
    """
    Solves a multi-regime first/second-order Markovian fluid queue. Q, R, S are
    lists (per regime) of the generator, diagonal drift and diagonal variance
    matrices; T is the vector of thresholds. Returns
    (masses, iniF, KF, cloF, iniB, KB, cloB).
    """
    K = len(T)
    N = np.asarray(Q[0]).shape[0]
    Q = [ml.matrix(Qk, dtype=float) for Qk in Q]
    R = [ml.matrix(Rk, dtype=float) for Rk in R]
    S = [ml.matrix(Sk, dtype=float) for Sk in S]

    if boundaryL is None:
        boundaryL = np.zeros(N)
    if boundaryU is None:
        boundaryU = boundaryL
    boundaryL = np.asarray(boundaryL, dtype=float).flatten()
    boundaryU = np.asarray(boundaryU, dtype=float).flatten()

    if Qt is None or (isinstance(Qt, (list, tuple)) and len(Qt) == 0):
        Qt = [Q[k] for k in range(K)] + [Q[K - 1]]
    elif len(Qt) == 1:
        Qt = [ml.matrix(Qt[0], dtype=float) for _ in range(K + 1)]
    else:
        Qt = [ml.matrix(x, dtype=float) for x in Qt]

    Tarr = np.concatenate(([0.0], np.asarray(T, dtype=float)))

    KF = [None] * K
    KB = [None] * K
    cloF = [None] * K
    cloB = [None] * K
    Np = np.zeros(K, dtype=int)
    Nn = np.zeros(K, dtype=int)
    Ns = np.zeros(K, dtype=int)
    NbF = np.zeros(K, dtype=int)
    NbB = np.zeros(K, dtype=int)
    vixp = [None] * K
    vixn = [None] * K
    vix0 = [None] * K
    vixs = [None] * K

    for k in range(K):
        S[k] = S[k] / 2.0
        ix = np.arange(N)
        dR = _dg(R[k])
        dS = _dg(S[k])
        ix0 = ix[(np.abs(dR) <= prec) & (dS <= prec)]
        ixn0 = np.setdiff1d(ix, ix0)

        Q00 = _sub(Q[k], ix0, ix0)
        Q0n = _sub(Q[k], ix0, ixn0)
        Qn0 = _sub(Q[k], ixn0, ix0)
        Qnn = _sub(Q[k], ixn0, ixn0)

        if len(ix0) > 0:
            Qv = Qnn + Qn0 * la.inv(-Q00) * Q0n
        else:
            Qv = Qnn
        Rv = _sub(R[k], ixn0, ixn0)
        Sv = _sub(S[k], ixn0, ixn0)
        Nnz = Qv.shape[0]

        ixl = np.arange(Nnz)
        dRv = _dg(Rv)
        dSv = _dg(Sv)
        ixp = ixl[(dRv > prec) & (dSv <= prec)]
        ixn = ixl[(dRv < -prec) & (dSv <= prec)]
        ixs = ixl[dSv > prec]
        Np[k] = len(ixp)
        Nn[k] = len(ixn)
        Ns[k] = len(ixs)

        # ----- FORWARD -----
        c1 = _safemax(-_dg(_sub(Qv, ixp, ixp)) / _dg(_sub(Rv, ixp, ixp))) if len(ixp) > 0 else -np.inf
        if len(ixs) > 0:
            RvS = _dg(_sub(Rv, ixs, ixs))
            SvS = _dg(_sub(Sv, ixs, ixs))
            QvS = _dg(_sub(Qv, ixs, ixs))
            discr = RvS ** 2 - 2 * (2 * SvS) * QvS
            cand = (discr > 0) * ((-RvS + np.sqrt(np.maximum(discr, 0))) / (2 * SvS))
            c2 = _safemax(cand)
        else:
            c2 = -np.inf
        c = max(c1, c2, 1.0)

        ixbF = np.concatenate((ixs, ixp)).astype(int)
        NbF[k] = len(ixbF)
        Bm = _blkdiag(c * _sub(Sv, ixbF, ixbF), -_sub(Rv, ixn, ixn))
        Lm = ml.matrix(np.block([
            [np.asarray(-_sub(Rv, ixbF, ixbF) - 2 * c * _sub(Sv, ixbF, ixbF)), np.zeros((NbF[k], Nn[k]))],
            [np.asarray(_sub(Qv, ixn, ixbF) / c), np.asarray(_sub(Qv, ixn, ixn) / c + _sub(Rv, ixn, ixn))]]))
        Fm = ml.matrix(np.block([
            [np.asarray(_sub(Qv, ixbF, ixbF) / c + c * _sub(Sv, ixbF, ixbF) + _sub(Rv, ixbF, ixbF)), np.asarray(_sub(Qv, ixbF, ixn) / c)],
            [np.zeros((Nn[k], NbF[k] + Nn[k]))]]))
        QBDR = ml.matrix(QBDFundamentalMatrices(Bm, Lm, Fm, "R", prec))
        KF[k] = (QBDR[:NbF[k], :NbF[k]] - ml.eye(NbF[k])) * c
        PsiF = QBDR[:NbF[k], NbF[k]:]
        clovF = ml.zeros((NbF[k], NbF[k] + Nn[k]))
        clovF[:, ixbF] = ml.eye(NbF[k])
        if len(ixn) > 0:
            clovF[:, ixn] = PsiF
        cloF[k] = ml.zeros((NbF[k], N))
        cloF[k][:, ixn0] = clovF
        if len(ix0) > 0:
            cloF[k][:, ix0] = clovF * Qn0 * la.inv(-Q00)

        # ----- BACKWARD -----
        c1 = _safemax(-_dg(_sub(Qv, ixn, ixn)) / (-_dg(_sub(Rv, ixn, ixn)))) if len(ixn) > 0 else -np.inf
        if len(ixs) > 0:
            RvS = _dg(_sub(Rv, ixs, ixs))
            SvS = _dg(_sub(Sv, ixs, ixs))
            QvS = _dg(_sub(Qv, ixs, ixs))
            discr = RvS ** 2 - 2 * (2 * SvS) * QvS
            cand = (discr > 0) * ((RvS + np.sqrt(np.maximum(discr, 0))) / (2 * SvS))
            c2 = _safemax(cand)
        else:
            c2 = -np.inf
        c = max(c1, c2, 1.0)

        ixbB = np.concatenate((ixs, ixn)).astype(int)
        NbB[k] = len(ixbB)
        Bm = _blkdiag(c * _sub(Sv, ixbB, ixbB), _sub(Rv, ixp, ixp))
        Lm = ml.matrix(np.block([
            [np.asarray(_sub(Rv, ixbB, ixbB) - 2 * c * _sub(Sv, ixbB, ixbB)), np.zeros((NbB[k], Np[k]))],
            [np.asarray(_sub(Qv, ixp, ixbB) / c), np.asarray(_sub(Qv, ixp, ixp) / c - _sub(Rv, ixp, ixp))]]))
        Fm = ml.matrix(np.block([
            [np.asarray(_sub(Qv, ixbB, ixbB) / c + c * _sub(Sv, ixbB, ixbB) - _sub(Rv, ixbB, ixbB)), np.asarray(_sub(Qv, ixbB, ixp) / c)],
            [np.zeros((Np[k], NbB[k] + Np[k]))]]))
        QBDR = ml.matrix(QBDFundamentalMatrices(Bm, Lm, Fm, "R", prec))
        KB[k] = (QBDR[:NbB[k], :NbB[k]] - ml.eye(NbB[k])) * c
        PsiB = QBDR[:NbB[k], NbB[k]:]
        clovB = ml.zeros((NbB[k], NbB[k] + Np[k]))
        clovB[:, ixbB] = ml.eye(NbB[k])
        if len(ixp) > 0:
            clovB[:, ixp] = PsiB
        cloB[k] = ml.zeros((NbB[k], N))
        cloB[k][:, ixn0] = clovB
        if len(ix0) > 0:
            cloB[k][:, ix0] = clovB * Qn0 * la.inv(-Q00)

        vixp[k] = ixn0[ixp]
        vixn[k] = ixn0[ixn]
        vix0[k] = ix0
        vixs[k] = ixn0[ixs]

    # ---- Boundary equations ----
    Neqns = (K + 1) * N + int(np.sum(Np)) + int(np.sum(Nn)) + 2 * int(np.sum(Ns))
    M = np.zeros((Neqns, Neqns))

    pos = [1]
    for k in range(K):
        pos += [N, int(NbF[k]), int(NbB[k])]
    pp = np.cumsum(pos)  # 1-based block-start positions

    def rows(i, length):
        s = pp[i] - 1
        return slice(s, s + length)

    # flux conservation
    i = 0  # MATLAB i=1 -> python 0 (pp index)
    for k in range(0, K + 1):
        M[rows(i, N), k * N:(k + 1) * N] = -np.asarray(Qt[k])
        if k > 0:
            M[rows(i - 2, NbF[k - 1]), k * N:(k + 1) * N] = np.asarray(
                la.expm(np.asarray(KF[k - 1]) * (Tarr[k] - Tarr[k - 1])) @ (-np.asarray(cloF[k - 1]) @ np.asarray(R[k - 1]) + np.asarray(KF[k - 1]) @ np.asarray(cloF[k - 1]) @ np.asarray(S[k - 1])))
            M[rows(i - 1, NbB[k - 1]), k * N:(k + 1) * N] = np.asarray(
                -np.asarray(cloB[k - 1]) @ np.asarray(R[k - 1]) - np.asarray(KB[k - 1]) @ np.asarray(cloB[k - 1]) @ np.asarray(S[k - 1]))
        if k < K:
            M[rows(i + 1, NbF[k]), k * N:(k + 1) * N] = np.asarray(
                np.asarray(cloF[k]) @ np.asarray(R[k]) - np.asarray(KF[k]) @ np.asarray(cloF[k]) @ np.asarray(S[k]))
            M[rows(i + 2, NbB[k]), k * N:(k + 1) * N] = np.asarray(
                la.expm(np.asarray(KB[k]) * (Tarr[k + 1] - Tarr[k])) @ (np.asarray(cloB[k]) @ np.asarray(R[k]) + np.asarray(KB[k]) @ np.asarray(cloB[k]) @ np.asarray(S[k])))
        i += 3

    col = (K + 1) * N + 1  # 1-based running column
    ix = np.arange(N)
    i = 0
    for k in range(0, K + 1):
        if k == 0:
            ixr0 = np.intersect1d(ix[boundaryL == 0], vixs[0])
            Nr0 = len(ixr0)
            ms0 = np.zeros((N, Np[0] + Nr0))
            sel = np.concatenate((vixp[0], ixr0)).astype(int)
            ms0[sel, :] = np.eye(Np[0] + Nr0)
            M[rows(i, N), col - 1:col - 1 + Np[0] + Nr0] = ms0
            col += Np[0] + Nr0
            ixa0 = np.intersect1d(ix[boundaryL == 1], vixs[0])
            Na0 = len(ixa0)
            pdfF = np.asarray(cloF[0])
            pdfB = np.asarray(la.expm(np.asarray(KB[0]) * Tarr[1]) @ np.asarray(cloB[0]))
            if Na0 > 0:
                M[rows(i + 1, NbF[0]), col - 1:col - 1 + Na0] = pdfF[:, ixa0]
                M[rows(i + 2, NbB[0]), col - 1:col - 1 + Na0] = pdfB[:, ixa0]
            col += Na0
        elif k == K:
            ixrB = np.intersect1d(ix[boundaryU == 0], vixs[K - 1])
            NrB = len(ixrB)
            ms0 = np.zeros((N, Nn[K - 1] + NrB))
            sel = np.concatenate((vixn[K - 1], ixrB)).astype(int)
            ms0[sel, :] = np.eye(Nn[K - 1] + NrB)
            M[rows(i, N), col - 1:col - 1 + Nn[K - 1] + NrB] = ms0
            col += Nn[K - 1] + NrB
            ixaB = np.intersect1d(ix[boundaryU == 1], vixs[K - 1])
            NaB = len(ixaB)
            pdfF = np.asarray(la.expm(np.asarray(KF[K - 1]) * (Tarr[K] - Tarr[K - 1])) @ np.asarray(cloF[K - 1]))
            pdfB = np.asarray(cloB[K - 1])
            if NaB > 0:
                M[rows(i - 2, NbF[K - 1]), col - 1:col - 1 + NaB] = pdfF[:, ixaB]
                M[rows(i - 1, NbB[K - 1]), col - 1:col - 1 + NaB] = pdfB[:, ixaB]
            col += NaB
        else:
            st0 = np.setdiff1d(ix, np.union1d(np.intersect1d(vixp[k - 1], vixn[k]), np.union1d(vix0[k - 1], vix0[k])))
            N0 = len(st0)
            ms0 = np.zeros((N, N0))
            ms0[st0, :] = np.eye(N0)
            M[rows(i, N), col - 1:col - 1 + N0] = ms0
            col += N0
            sts = np.setdiff1d(np.union1d(vixs[k - 1], vixs[k]), np.union1d(vixn[k], vixp[k - 1]))
            Nss = len(sts)
            sqrtSk = np.sqrt(np.asarray(S[k - 1]))
            sqrtSk1 = np.sqrt(np.asarray(S[k]))
            BelowF = np.asarray(la.expm(np.asarray(KF[k - 1]) * (Tarr[k] - Tarr[k - 1])) @ (-np.asarray(cloF[k - 1])) @ sqrtSk)
            BelowB = np.asarray(-np.asarray(cloB[k - 1]) @ sqrtSk)
            AboveF = np.asarray(np.asarray(cloF[k]) @ sqrtSk1)
            AboveB = np.asarray(la.expm(np.asarray(KB[k]) * (Tarr[k + 1] - Tarr[k])) @ np.asarray(cloB[k]) @ sqrtSk1)
            if Nss > 0:
                M[rows(i - 2, NbF[k - 1]), col - 1:col - 1 + Nss] = BelowF[:, sts]
                M[rows(i - 1, NbB[k - 1]), col - 1:col - 1 + Nss] = BelowB[:, sts]
                M[rows(i + 1, NbF[k]), col - 1:col - 1 + Nss] = AboveF[:, sts]
                M[rows(i + 2, NbB[k]), col - 1:col - 1 + Nss] = AboveB[:, sts]
            col += Nss
        i += 3

    # normalizing condition
    h = np.ones(N)
    for k in range(K):
        sumKF, sumKB = _integExp2(KF[k], KB[k], Tarr[k + 1] - Tarr[k])
        h = np.concatenate((h,
                            np.asarray(sumKF @ np.asarray(cloF[k])).sum(axis=1).flatten(),
                            np.asarray(sumKB @ np.asarray(cloB[k])).sum(axis=1).flatten(),
                            np.ones(N)))

    M[:, 0] = h
    b = la.solve(M.T, np.concatenate(([1.0], np.zeros(len(h) - 1))))

    masses = [None] * (K + 1)
    iniF = [None] * K
    iniB = [None] * K
    masses[0] = ml.matrix(b[:N].reshape(1, -1))
    i = 1  # 1-based pp index for the per-regime blocks (MATLAB i=2 -> pp index 2 -> python cumsum idx 1)
    for k in range(K):
        iniF[k] = ml.matrix(b[pp[i] - 1:pp[i] - 1 + NbF[k]].reshape(1, -1))
        iniB[k] = ml.matrix(b[pp[i + 1] - 1:pp[i + 1] - 1 + NbB[k]].reshape(1, -1))
        masses[k + 1] = ml.matrix(b[pp[i + 2] - 1:pp[i + 2] - 1 + N].reshape(1, -1))
        i += 3

    return masses, iniF, KF, cloF, iniB, KB, cloB


def LevelDependentFluidStationaryDistr(masses, iniF, KF, cloF, iniB, KB, cloB, T, what, points):
    """
    Stationary distribution ('pdf','pdfd','cdf','cdfm') of a first/second-order
    level-dependent fluid queue at the requested points.
    """
    K = len(T)
    Tarr = np.concatenate(([0.0], np.asarray(T, dtype=float)))
    N = np.asarray(masses[0]).size
    cummulate = what in ('cdf', 'cdfm')
    res = []
    for p in np.atleast_1d(points):
        pres = np.zeros(N)
        k = 0
        while k < K and p >= Tarr[k]:
            if cummulate:
                if k > 0:
                    sumKF, sumKB = _integExp2(KF[k - 1], KB[k - 1], Tarr[k] - Tarr[k - 1])
                    val = np.asarray(iniF[k - 1]) @ sumKF @ np.asarray(cloF[k - 1]) + np.asarray(iniB[k - 1]) @ sumKB @ np.asarray(cloB[k - 1])
                    pres = pres + np.asarray(val).flatten()
                if p > Tarr[k] or what == 'cdfm':
                    pres = pres + np.asarray(masses[k]).flatten()
            k += 1
        if k == K and p == Tarr[k] and what == 'cdfm':
            pres = pres + np.asarray(masses[k]).flatten()
        ki = k - 1  # regime index (0-based)
        prem = p - Tarr[ki]
        Tk = Tarr[ki + 1] - Tarr[ki]
        if what == 'pdf':
            pres = np.asarray(np.asarray(iniF[ki]) @ la.expm(np.asarray(KF[ki]) * prem) @ np.asarray(cloF[ki])
                              + np.asarray(iniB[ki]) @ la.expm(np.asarray(KB[ki]) * (Tk - prem)) @ np.asarray(cloB[ki])).flatten()
        elif what == 'pdfd':
            pres = np.asarray(np.asarray(iniF[ki]) @ np.asarray(KF[ki]) @ la.expm(np.asarray(KF[ki]) * prem) @ np.asarray(cloF[ki])
                              - np.asarray(iniB[ki]) @ np.asarray(KB[ki]) @ la.expm(np.asarray(KB[ki]) * (Tk - prem)) @ np.asarray(cloB[ki])).flatten()
        elif what in ('cdf', 'cdfm'):
            sumKF, sumKB = _integExp2(KF[ki], KB[ki], prem)
            pres = pres + np.asarray(np.asarray(iniF[ki]) @ sumKF @ np.asarray(cloF[ki])
                                     + np.asarray(iniB[ki]) @ la.expm(np.asarray(KB[ki]) * (Tk - prem)) @ sumKB @ np.asarray(cloB[ki])).flatten()
        res.append(pres)
    return np.array(res)


def LevelDependentFluidStationaryMean(masses, iniF, KF, cloF, iniB, KB, cloB, T):
    """
    Stationary mean fluid level E[X] of a first/second-order level-dependent
    fluid queue (closed form).
    """
    K = len(T)
    Tarr = np.concatenate(([0.0], np.asarray(T, dtype=float)))
    N = np.asarray(masses[0]).size
    h = np.ones((N, 1))

    def expIntMoments(Mmat, L):
        Mmat = np.asarray(Mmat)
        n = Mmat.shape[0]
        Zn = np.zeros((n, n))
        In = np.eye(n)
        A = np.block([[Mmat, In, Zn], [Zn, Zn, In], [Zn, Zn, Zn]])
        W = la.expm(A * L)
        J0 = W[:n, n:2 * n]
        W13 = W[:n, 2 * n:3 * n]
        J1 = L * J0 - W13
        return J0, J1

    res = 0.0
    for j in range(K + 1):
        res += Tarr[j] * float(np.sum(masses[j]))
    for k in range(K):
        Tk = Tarr[k + 1] - Tarr[k]
        J0F, J1F = expIntMoments(KF[k], Tk)
        J0B, J1B = expIntMoments(KB[k], Tk)
        res += float(np.asarray(iniF[k]) @ (Tarr[k] * J0F + J1F) @ np.asarray(cloF[k]) @ h)
        res += float(np.asarray(iniB[k]) @ (Tarr[k + 1] * J0B - J1B) @ np.asarray(cloB[k]) @ h)
    return res


def multiregime(Q, R, Qt, Rt, T, pdfpoints, cdfpoints):
    """
    Multi-regime feedback fluid queue (Kankaya & Akar). Q, R are lists (per
    regime k=1..K) of the generator and diagonal drift-rate VECTORS; Qt, Rt are
    lists (per boundary k=0..K) of boundary generators and rate vectors; T is
    the vector of regime thresholds. Returns (pdf, pdfd, cdf, cdfm).
    """
    K = len(R)
    N = len(np.asarray(R[0]).flatten())
    R = [np.asarray(Rk, dtype=float).flatten() for Rk in R]
    Q = [ml.matrix(Qk, dtype=float) for Qk in Q]
    Tarr = np.concatenate(([0.0], np.asarray(T, dtype=float)))

    if len(Q) == 1:
        Q = [Q[0] for _ in range(K)]
    if Qt is None or len(Qt) == 0:
        Qt = [Q[0]] + [Q[k] for k in range(K)]
    elif len(Qt) == 1:
        Qt = [ml.matrix(Qt[0], dtype=float) for _ in range(K + 1)]
    else:
        Qt = [ml.matrix(x, dtype=float) for x in Qt]
    if Rt is None or len(Rt) == 0:
        Rt = [R[0]] + [R[k] for k in range(K)]
    else:
        Rt = [np.asarray(x, dtype=float).flatten() for x in Rt]

    An = [None] * K; Ap = [None] * K
    L0 = [None] * K; Ln = [None] * K; Lp = [None] * K
    M0 = [None] * K; MT = [None] * K; Mi = [None] * K
    zeig = [0] * K; neig = [0] * K; peig = [0] * K
    Nnz = []
    ix = np.arange(N)
    for k in range(K):
        zix = ix[R[k] == 0]
        nzix = ix[R[k] != 0]
        Nn = len(nzix)
        Qk = np.asarray(Q[k])
        if len(zix) > 0:
            Qnk = Qk[np.ix_(nzix, nzix)] + Qk[np.ix_(nzix, zix)] @ la.inv(-Qk[np.ix_(zix, zix)]) @ Qk[np.ix_(zix, nzix)]
        else:
            Qnk = Qk[np.ix_(nzix, nzix)]
        A = Qnk @ np.diag(1.0 / R[k][nzix])
        Z, D, zeroeig, negeig, poseig = _ordered_schur(A)
        zeig[k], neig[k], peig[k] = zeroeig, negeig, poseig

        D22 = D[zeroeig:, zeroeig:]
        D12 = D[:zeroeig, zeroeig:]
        if zeroeig > 0:
            X1 = la.solve_sylvester(np.zeros((zeroeig, zeroeig)), -D22, D12)
        else:
            X1 = np.zeros((0, Nn - zeroeig))
        Dnn = D[zeroeig:zeroeig + negeig, zeroeig:zeroeig + negeig]
        Dpp = D[zeroeig + negeig:, zeroeig + negeig:]
        Dnp = D[zeroeig:zeroeig + negeig, zeroeig + negeig:]
        if negeig > 0 and poseig > 0:
            X2 = la.solve_sylvester(Dnn, -Dpp, Dnp)
        else:
            X2 = np.zeros((negeig, poseig))
        blk1 = np.block([[np.eye(zeroeig), -X1],
                         [np.zeros((Nn - zeroeig, zeroeig)), np.eye(Nn - zeroeig)]])
        blk2 = np.block([[np.eye(zeroeig), np.zeros((zeroeig, Nn - zeroeig))],
                         [np.zeros((negeig, zeroeig)), np.eye(negeig), -X2],
                         [np.zeros((poseig, zeroeig + negeig)), np.eye(poseig)]])
        Y = Z @ blk1 @ blk2
        iY = la.inv(Y)
        iY0 = iY[:zeroeig, :]
        iYn = iY[zeroeig:zeroeig + negeig, :]
        iYp = iY[zeroeig + negeig:, :]
        At = iY @ A @ Y
        An[k] = At[zeroeig:zeroeig + negeig, zeroeig:zeroeig + negeig]
        Ap[k] = At[zeroeig + negeig:, zeroeig + negeig:]

        def buildL(iYx):
            L = np.zeros((iYx.shape[0], N))
            L[:, nzix] = iYx
            if len(zix) > 0:
                L[:, zix] = iYx @ Qk[np.ix_(nzix, zix)] @ la.inv(-Qk[np.ix_(zix, zix)])
            return L
        L0[k], Ln[k], Lp[k] = buildL(iY0), buildL(iYn), buildL(iYp)
        Tk = Tarr[k + 1] - Tarr[k]
        eAnT = la.expm(An[k] * Tk) if negeig > 0 else np.zeros((0, 0))
        eApT = la.expm(-Ap[k] * Tk) if poseig > 0 else np.zeros((0, 0))
        M0[k] = np.vstack([L0[k], Ln[k], eApT @ Lp[k]])
        MT[k] = np.vstack([L0[k], eAnT @ Ln[k], Lp[k]])
        MiN = la.inv(-An[k]) @ (np.eye(negeig) - eAnT) @ Ln[k] if negeig > 0 else np.zeros((0, N))
        MiP = la.inv(Ap[k]) @ (np.eye(poseig) - eApT) @ Lp[k] if poseig > 0 else np.zeros((0, N))
        Mi[k] = np.vstack([Tk * L0[k], MiN, MiP])
        Nnz.append(Nn)

    d = (K + 1) * N
    Neqns = d + int(np.sum(Nnz))
    M = np.zeros((Neqns, Neqns))

    def dstart(kidx):
        return d + int(np.sum(Nnz[:kidx]))

    # eq. (12)
    p = 0
    M[0:N, p:p + N] = -np.asarray(Qt[0])
    M[dstart(0):dstart(0) + Nnz[0], p:p + N] = M0[0] @ np.diag(R[0])
    p += N
    # eq. (13)
    for k in range(1, K):
        M[k * N:(k + 1) * N, p:p + N] = -np.asarray(Qt[k])
        M[dstart(k):dstart(k) + Nnz[k], p:p + N] = M0[k] @ np.diag(R[k])
        M[dstart(k - 1):dstart(k - 1) + Nnz[k - 1], p:p + N] = -MT[k - 1] @ np.diag(R[k - 1])
        p += N
    # eq. (16)
    M[K * N:(K + 1) * N, p:p + N] = -np.asarray(Qt[K])
    M[dstart(K - 1):dstart(K - 1) + Nnz[K - 1], p:p + N] = -MT[K - 1] @ np.diag(R[K - 1])
    p += N
    # eq. (8)
    for m in range(N):
        if R[0][m] > 0:
            M[:, p] = 0; M[m, p] = 1; p += 1
    # eq. (11)
    for m in range(N):
        if R[K - 1][m] < 0:
            M[:, p] = 0; M[K * N + m, p] = 1; p += 1
    # eq. (9)
    for k in range(1, K):
        for m in range(N):
            if (R[k - 1][m] > 0 and R[k][m] > 0) or (R[k - 1][m] < 0 and R[k][m] < 0):
                M[:, p] = 0; M[k * N + m, p] = 1; p += 1
    # eq. (10)
    for k in range(1, K):
        for m in range(N):
            if (R[k - 1][m] < 0 and R[k][m] > 0) and Rt[k][m] != 0:
                M[:, p] = 0; M[k * N + m, p] = 1; p += 1
    # eq. (14)
    for k in range(1, K):
        for m in range(N):
            if R[k - 1][m] < 0 and Rt[k][m] >= 0:
                M[:, p] = 0
                M[dstart(k - 1):dstart(k - 1) + Nnz[k - 1], p] = MT[k - 1][:, m]
                p += 1
    # eq. (15)
    for k in range(1, K):
        for m in range(N):
            if R[k][m] > 0 and Rt[k][m] <= 0:
                M[:, p] = 0
                M[dstart(k):dstart(k) + Nnz[k], p] = M0[k][:, m]
                p += 1

    # normalization
    M[0:d, 0] = 1.0
    for k in range(K):
        M[dstart(k):dstart(k) + Nnz[k], 0] = Mi[k].sum(axis=1)
    rhs = np.zeros(Neqns); rhs[0] = 1.0
    sol = la.solve(M.T, rhs)

    masses = [None] * (K + 1)
    a0 = [None] * K; an = [None] * K; ap = [None] * K
    masses[0] = sol[0:N]
    for k in range(K):
        masses[k + 1] = sol[k * N:(k + 1) * N] if k < K else None
    for k in range(K):
        masses[k + 1] = sol[(k + 1) * N:(k + 2) * N]
        avec = sol[dstart(k):dstart(k) + Nnz[k]]
        a0[k] = avec[:zeig[k]]
        an[k] = avec[zeig[k]:zeig[k] + neig[k]]
        ap[k] = avec[zeig[k] + neig[k]:zeig[k] + neig[k] + peig[k]]

    def _regime(p):
        k = 0
        while k < K - 1 and p >= Tarr[k + 1]:
            k += 1
        return k  # 0-based regime index

    pdf = []
    for pt in np.atleast_1d(pdfpoints):
        k = _regime(pt)
        presp = a0[k] @ L0[k]
        if neig[k] > 0:
            presp = presp + an[k] @ la.expm(An[k] * (pt - Tarr[k])) @ Ln[k]
        if peig[k] > 0:
            presp = presp + ap[k] @ la.expm(-Ap[k] * (Tarr[k + 1] - pt)) @ Lp[k]
        pdf.append(np.asarray(presp).flatten())
    pdf = np.array(pdf)

    pdfd = []
    for pt in np.atleast_1d(pdfpoints):
        k = _regime(pt)
        presp = np.zeros(N)
        if neig[k] > 0:
            presp = presp + an[k] @ An[k] @ la.expm(An[k] * (pt - Tarr[k])) @ Ln[k]
        if peig[k] > 0:
            presp = presp + ap[k] @ Ap[k] @ la.expm(-Ap[k] * (Tarr[k + 1] - pt)) @ Lp[k]
        pdfd.append(np.asarray(presp).flatten())
    pdfd = np.array(pdfd)

    cdf = []; cdfm = []
    for c in np.atleast_1d(cdfpoints):
        cres = np.zeros(N); cresm = np.zeros(N)
        k = 0  # loop counter (1-based regime accumulation); MATLAB T(k+1) == Tarr[k]
        while k < K and c >= Tarr[k]:
            if k > 0:
                avec = np.concatenate((a0[k - 1], an[k - 1], ap[k - 1]))
                cres = cres + avec @ Mi[k - 1]
                cresm = cresm + avec @ Mi[k - 1]
            cresm = cresm + masses[k]
            if c > Tarr[k]:
                cres = cres + masses[k]
            k += 1
        if k == K and c == Tarr[k]:
            cresm = cresm + masses[k]
        kk = k - 1  # 0-based regime containing c
        crem = c - Tarr[kk]
        Tk = Tarr[kk + 1] - Tarr[kk]
        val = a0[kk] @ L0[kk] * crem
        if neig[kk] > 0:
            val = val + an[kk] @ la.inv(-An[kk]) @ (np.eye(neig[kk]) - la.expm(An[kk] * crem)) @ Ln[kk]
        if peig[kk] > 0:
            val = val + ap[kk] @ la.inv(-Ap[kk]) @ (la.expm(-Ap[kk] * Tk) - la.expm(-Ap[kk] * (Tk - crem))) @ Lp[kk]
        val = np.asarray(val).flatten()
        cres = cres + val
        cresm = cresm + val
        cdf.append(cres); cdfm.append(cresm)
    return pdf, pdfd, np.array(cdf), np.array(cdfm)
