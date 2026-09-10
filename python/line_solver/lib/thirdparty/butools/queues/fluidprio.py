# -*- coding: utf-8 -*-
"""
FluidPrioQueue: performance measures of a continuous-time fluid priority queue.

Ported from BUTools-family fluid tools (G. Horvath).

References
----------
G. Horvath, "Efficient analysis of the MMAP[K]/PH[K]/1 priority queue",
European Journal of Operational Research, 246(1):128-139, 2015.
"""
import math

import numpy as np
import numpy.matlib as ml
import scipy.linalg as la

from ..mam.fluid import GeneralFluidSolve, FluidFundamentalMatrices, FluidStationaryDistr
from ..mc.stst import CTMCSolve
from .flufluqueue import FluFluQueue



__all__ = ["FluidPrioQueue"]

def _dg(M):
    return np.diag(np.asarray(M))


def _diag(v):
    return ml.matrix(np.diagflat(np.asarray(v).flatten()))


def _lyap(A, B, C):
    """Solve A*X + X*B + C = 0."""
    return ml.matrix(la.solve_sylvester(np.asarray(A), np.asarray(B), -np.asarray(C)))


def FluidPrioQueue(Q, R, d, *args, prec=1e-14, erlMaxOrder=200, classes=None):
    """
    Returns performance measures of a continuous-time fluid priority queue.

    Parameters
    ----------
    Q : (N,N) generator of the modulating Markov chain
    R : (K,N) per-class input fluid rates in the background states
    d : constant fluid service rate (scalar, positive)
    *args : measure/points pairs, e.g. 'flMoms', n, 'stDistr', points
    prec, erlMaxOrder, classes : options

    Returns
    -------
    A list with one entry per requested performance measure (each entry a
    matrix whose columns belong to the various job types). A single item is
    returned directly (not wrapped in a list).
    """
    Q = ml.matrix(Q, dtype=float)
    R = ml.matrix(R, dtype=float)
    K = R.shape[0]
    N = Q.shape[0]
    if classes is None:
        classes = list(range(1, K + 1))

    # -------- auxiliary routines (closures over prec / erlMaxOrder) --------
    def DReward(Qm, Rm, n):
        Qm = np.asarray(Qm, dtype=float)
        Rm = np.asarray(Rm, dtype=float)
        NQ = Qm.shape[0]
        if NQ == 0:
            return ml.zeros((0, 0))
        ix = np.arange(NQ)
        dR = np.diag(Rm)
        ixz = ix[np.abs(dR) <= prec]
        ixp = np.concatenate((ix[dR > prec], ix[dR < -prec])).astype(int)
        Nz = len(ixz)
        Np_ = len(ixp)
        Per = np.zeros((NQ, NQ))
        for iw in range(Nz):
            Per[iw, ixz[iw]] = 1.0
        for iw in range(Np_):
            Per[Nz + iw, ixp[iw]] = 1.0
        iPer = la.inv(Per)
        Rp = Rm[np.ix_(ixp, ixp)]
        Qpp = Qm[np.ix_(ixp, ixp)]
        Qpz = Qm[np.ix_(ixp, ixz)]
        Qzp = Qm[np.ix_(ixz, ixp)]
        Qzz = Qm[np.ix_(ixz, ixz)]
        iRp = la.inv(Rp)
        if Nz > 0:
            inner = iRp @ (-Qpp - Qpz @ la.inv(-Qzz) @ Qzp)
        else:
            inner = iRp @ (-Qpp)
        dXvn = (-1) ** n * math.factorial(n) * np.linalg.matrix_power(la.inv(inner), n + 1) @ iRp
        if Nz > 0:
            iQzz = la.inv(-Qzz)
            drpar = np.block([
                [iQzz @ Qzp @ dXvn @ Qpz @ iQzz, iQzz @ Qzp @ dXvn],
                [dXvn @ Qpz @ iQzz, dXvn]])
        else:
            drpar = dXvn
        return ml.matrix(iPer @ drpar @ Per)

    def _partition(C):
        C = np.asarray(C)
        dC = np.diag(C)
        ix = np.arange(C.shape[0])
        ixz = ix[np.abs(dC) <= prec]
        ixp = ix[dC > prec]
        ixn = ix[dC < -prec]
        return ixz, ixp, ixn

    def _perm(NF, ixz, ixp, ixn):
        Per = np.zeros((NF, NF))
        for i in range(len(ixz)):
            Per[i, ixz[i]] = 1.0
        for i in range(len(ixp)):
            Per[len(ixz) + i, ixp[i]] = 1.0
        for i in range(len(ixn)):
            Per[len(ixz) + len(ixp) + i, ixn[i]] = 1.0
        return Per

    def _blocks(PF, Nz, Np_, Nn):
        # returns Fzz,Fpz,Fmz,Fzp,Fpp,Fmp,Fzm,Fpm,Fmm (row,col order first letter=row)
        r = [slice(0, Nz), slice(Nz, Nz + Np_), slice(Nz + Np_, Nz + Np_ + Nn)]
        def B(a, b):
            return ml.matrix(PF[r[a], r[b]])
        Fzz, Fpz, Fmz = B(0, 0), B(1, 0), B(2, 0)
        Fzp, Fpp, Fmp = B(0, 1), B(1, 1), B(2, 1)
        Fzm, Fpm, Fmm = B(0, 2), B(1, 2), B(2, 2)
        return Fzz, Fpz, Fmz, Fzp, Fpp, Fmp, Fzm, Fpm, Fmm

    def BusyPeriodRewardMoms(F, C, D, numOfMoms):
        F = ml.matrix(F); C = ml.matrix(C); D = ml.matrix(D)
        NF = F.shape[0]
        ixz, ixp, ixn = _partition(C)
        Nz, Np_, Nn = len(ixz), len(ixp), len(ixn)
        Per = _perm(NF, ixz, ixp, ixn)
        iPer = la.inv(Per)
        PF = np.asarray(Per @ np.asarray(F) @ iPer)
        Fzz, Fpz, Fmz, Fzp, Fpp, Fmp, Fzm, Fpm, Fmm = _blocks(PF, Nz, Np_, Nn)
        Cm = ml.matrix(np.asarray(C)[np.ix_(ixn, ixn)]) if Nn > 0 else ml.zeros((0, 0))
        Cp = ml.matrix(np.asarray(C)[np.ix_(ixp, ixp)]) if Np_ > 0 else ml.zeros((0, 0))
        Dm = ml.matrix(np.asarray(D)[np.ix_(ixn, ixn)]) if Nn > 0 else ml.zeros((0, 0))
        Dp = ml.matrix(np.asarray(D)[np.ix_(ixp, ixp)]) if Np_ > 0 else ml.zeros((0, 0))
        Dz = ml.matrix(np.asarray(D)[np.ix_(ixz, ixz)]) if Nz > 0 else ml.zeros((0, 0))

        iFzz = la.inv(-Fzz) if Nz > 0 else ml.zeros((0, 0))
        iCp = la.inv(Cp) if Np_ > 0 else ml.zeros((0, 0))
        iCm = la.inv(-Cm) if Nn > 0 else ml.zeros((0, 0))

        Fppd = [None] * (numOfMoms + 1)
        Fpmd = [None] * (numOfMoms + 1)
        Fmpd = [None] * (numOfMoms + 1)
        Fmmd = [None] * (numOfMoms + 1)
        Fppd[0] = iCp * (Fpp + Fpz * iFzz * Fzp)
        Fpmd[0] = iCp * (Fpm + Fpz * iFzz * Fzm)
        Fmpd[0] = iCm * (Fmp + Fmz * iFzz * Fzp)
        Fmmd[0] = iCm * (Fmm + Fmz * iFzz * Fzm)
        for i in range(1, numOfMoms + 1):
            dr = DReward(Fzz, Dz, i)
            Fppd[i] = iCp * Fpz * dr * Fzp
            Fpmd[i] = iCp * Fpz * dr * Fzm
            Fmpd[i] = iCm * Fmz * dr * Fzp
            Fmmd[i] = iCm * Fmz * dr * Fzm
            if i == 1:
                Fppd[i] = Fppd[i] - iCp * Dp
                Fmmd[i] = Fmmd[i] - iCm * Dm
        Psi = ml.matrix(FluidFundamentalMatrices(Fppd[0], Fpmd[0], Fmpd[0], Fmmd[0], "P", prec))
        BPM = [None] * (numOfMoms + 1)
        BPM[0] = Psi
        for i in range(1, numOfMoms + 1):
            X = -Psi * Fmpd[i] * Psi + Fpmd[i]
            for m in range(0, i):
                X = X + math.comb(i, m) * ((Fppd[i - m] + Psi * Fmpd[i - m]) * BPM[m] + BPM[m] * (Fmmd[i - m] + Fmpd[i - m] * Psi))
            for l in range(1, i):
                for m in range(1, i - l + 1):
                    X = X + math.comb(i, l) * math.comb(i - l, m) * BPM[l] * Fmpd[i - l - m] * BPM[m]
            BPM[i] = _lyap(Fppd[0] + Psi * Fmpd[0], Fmmd[0] + Fmpd[0] * Psi, X)
        # re-order states back to original ordering
        out = []
        for i in range(len(BPM)):
            emb = np.zeros((NF, NF))
            emb[Nz:Nz + Np_, Nz + Np_:Nz + Np_ + Nn] = np.asarray(BPM[i])
            out.append(ml.matrix(iPer @ emb @ Per))
        return out

    def BusyPeriodRewardDistr(F, C, D, t):
        F = ml.matrix(F); C = ml.matrix(C); D = ml.matrix(D)
        NF = F.shape[0]
        ixz, ixp, ixn = _partition(C)
        Nz, Np_, Nn = len(ixz), len(ixp), len(ixn)
        Per = _perm(NF, ixz, ixp, ixn)
        iPer = la.inv(Per)
        PF = np.asarray(Per @ np.asarray(F) @ iPer)
        Fzz, Fpz, Fmz, Fzp, Fpp, Fmp, Fzm, Fpm, Fmm = _blocks(PF, Nz, Np_, Nn)
        Cm = ml.matrix(np.asarray(C)[np.ix_(ixn, ixn)]) if Nn > 0 else ml.zeros((0, 0))
        Cp = ml.matrix(np.asarray(C)[np.ix_(ixp, ixp)]) if Np_ > 0 else ml.zeros((0, 0))
        Dm = ml.matrix(np.asarray(D)[np.ix_(ixn, ixn)]) if Nn > 0 else ml.zeros((0, 0))
        Dp = ml.matrix(np.asarray(D)[np.ix_(ixp, ixp)]) if Np_ > 0 else ml.zeros((0, 0))
        Dz = ml.matrix(np.asarray(D)[np.ix_(ixz, ixz)]) if Nz > 0 else ml.zeros((0, 0))
        iCp = la.inv(Cp) if Np_ > 0 else ml.zeros((0, 0))
        iCm = la.inv(-Cm) if Nn > 0 else ml.zeros((0, 0))

        L = erlMaxOrder
        nu = L / t
        Z = la.inv(nu * Dz - Fzz) if Nz > 0 else ml.zeros((0, 0))
        AFpp = iCp * (Fpp - nu * Dp + Fpz * Z * Fzp)
        AFpm = iCp * (Fpm + Fpz * Z * Fzm)
        AFmp = iCm * (Fmp + Fmz * Z * Fzp)
        AFmm = iCm * (Fmm - nu * Dm + Fmz * Z * Fzm)
        Psie = ml.matrix(FluidFundamentalMatrices(AFpp, AFpm, AFmp, AFmm, "P", prec))
        Pn = [Psie]
        pr = Psie.copy()
        AM = AFpp + Psie * AFmp
        BM = AFmm + AFmp * Psie
        nuDz_Z = (nu * Dz) * Z if Nz > 0 else ml.zeros((0, 0))
        for n in range(1, L):
            CM = iCp * (nu * Dp) * Pn[n - 1] + Pn[n - 1] * (iCm * (nu * Dm))
            for i in range(1, n):
                CM = CM + Pn[i] * (iCm * Fmp) * Pn[n - i]
            if Nz > 0:
                nuDzZ_n = np.linalg.matrix_power(np.asarray(nuDz_Z), n)
                CM = CM + iCp * Fpz * Z * ml.matrix(nuDzZ_n) * Fzm - Psie * iCm * Fmz * Z * ml.matrix(nuDzZ_n) * Fzp * Psie
                for i in range(0, n):
                    nuDzZ_ni = ml.matrix(np.linalg.matrix_power(np.asarray(nuDz_Z), n - i))
                    CM = CM + Pn[i] * iCm * Fmz * Z * nuDzZ_ni * (Fzm + Fzp * Psie)
                    CM = CM + (iCp * Fpz + Psie * iCm * Fmz) * Z * nuDzZ_ni * Fzp * Pn[i]
                for i in range(1, n):
                    for j in range(1, n - i + 1):
                        nuDzZ_nij = ml.matrix(np.linalg.matrix_power(np.asarray(nuDz_Z), n - i - j))
                        CM = CM + Pn[i] * iCm * Fmz * Z * nuDzZ_nij * Fzp * Pn[j]
            PM = _lyap(AM, BM, CM)
            Pn.append(PM)
            pr = pr + PM
        # re-order
        embed = lambda Mx: ml.matrix(iPer @ (_embed(Mx, NF, Nz, Np_, Nn)) @ Per)
        pr_full = embed(pr)
        Pn_full = [embed(P) for P in Pn]
        return pr_full, Pn_full

    def _embed(Mx, NF, Nz, Np_, Nn):
        emb = np.zeros((NF, NF))
        emb[Nz:Nz + Np_, Nz + Np_:Nz + Np_ + Nn] = np.asarray(Mx)
        return emb

    # -------- preparation --------
    pi = ml.matrix(CTMCSolve(Q))
    lambda_ = np.asarray(pi * R.T).flatten()   # per-class arrival rate

    Ret = []
    for k in classes:                          # k is 1-based class index
        kk = k - 1                             # 0-based
        # workload process for classes of same-or-higher priority
        Rsum_hi = np.asarray(R[kk:, :]).sum(axis=0)      # sum over classes k..K
        mass0, ini, Km, clo = GeneralFluidSolve(Q, _diag(Rsum_hi) / d - ml.eye(N), [], prec)
        mass0 = ml.matrix(mass0); ini = ml.matrix(ini); Km = ml.matrix(Km); clo = ml.matrix(clo)
        KN = Km.shape[0]
        clok = clo * _diag(np.asarray(R[kk, :]).flatten()) / lambda_[kk]

        Delta = _diag(np.asarray(la.inv(-Km) @ np.asarray(clok).sum(axis=1)).flatten())
        K0 = la.inv(Delta) * Km * Delta
        K1 = la.inv(Delta) * clok
        kappa = ini * Delta
        Km, clok, ini = K0, K1, kappa

        Rk = _diag(np.asarray(R[kk, :]).flatten())
        if k < K:
            Rsum_lo = np.asarray(R[kk + 1:, :]).sum(axis=0)
            F = ml.matrix(np.block([[np.asarray(Km), np.asarray(clok)], [np.zeros((N, KN)), np.asarray(Q)]]))
            Clo = ml.matrix(np.block([[np.eye(KN), np.zeros((KN, N))], [np.zeros((N, KN)), np.asarray(_diag(Rsum_lo) / d - ml.eye(N))]]))
            for opt, val in _pairs(args):
                if opt == 'stMoms':
                    D = ml.matrix(np.block([[np.zeros((KN, KN + N))], [np.zeros((N, KN)), np.eye(N)]]))
                    Tmp = BusyPeriodRewardMoms(F, Clo, D, val)
                    inis = ml.matrix(np.hstack((np.asarray(ini), np.zeros((1, N)))))
                    stMoms = np.array([(-1) ** i * float(np.sum(inis * Tmp[i])) for i in range(1, len(Tmp))])
                    Ret.append(stMoms)
                elif opt == 'flMoms':
                    D = ml.matrix(np.block([[np.zeros((KN, KN + N))], [np.zeros((N, KN)), np.asarray(Rk)]]))
                    FLDPn = BusyPeriodRewardMoms(F, Clo, D, val)
                    FLDPn = [ml.matrix(np.asarray(X)[:, KN:]) for X in FLDPn]
                    inis = ml.matrix(np.hstack((np.asarray(ini), np.zeros((1, N)))))
                    FLDv = []
                    for i in range(len(FLDPn)):
                        Xi = inis * FLDPn[i]
                        if i == 0:
                            Xi = Xi + mass0 * Rk / lambda_[kk]
                        FLDv.append(Xi)
                    fldMoms = np.zeros(val)
                    FLPn = [pi]
                    iTerm = la.inv(ml.ones((N, 1)) * pi - Q)
                    for n in range(1, val + 1):
                        sumP = float(np.sum(FLDv[n])) + n * float((-FLDv[n - 1] + FLPn[n - 1] * Rk / lambda_[kk]) * iTerm * R[kk, :].T)
                        P = sumP * pi + n * (-FLPn[n - 1] * Rk + FLDv[n - 1] * lambda_[kk]) * iTerm
                        FLPn.append(P)
                        fldMoms[n - 1] = (-1) ** n * float(np.sum(P))
                    Ret.append(fldMoms)
                elif opt == 'stDistr':
                    D = ml.matrix(np.block([[np.zeros((KN, KN + N))], [np.zeros((N, KN)), np.eye(N)]]))
                    inis = ml.matrix(np.hstack((np.asarray(ini), np.zeros((1, N)))))
                    pts = np.atleast_1d(val)
                    res = np.zeros(len(pts))
                    for xi, x in enumerate(pts):
                        Tmp, _ = BusyPeriodRewardDistr(F, Clo, D, x)
                        res[xi] = float(np.sum(mass0 * Rk / lambda_[kk])) + float(np.sum(inis * Tmp))
                    Ret.append(res)
                elif opt == 'flDistr':
                    D = ml.matrix(np.block([[np.zeros((KN, KN + N))], [np.zeros((N, KN)), np.asarray(Rk)]]))
                    inis = ml.matrix(np.hstack((np.asarray(ini), np.zeros((1, N)))))
                    pts = np.atleast_1d(val)
                    res = np.zeros(len(pts))
                    for xi, x in enumerate(pts):
                        nu = erlMaxOrder / x
                        _, Psix = BusyPeriodRewardDistr(F, Clo, D, x)
                        Psiy = lambda_[kk] * nu * (mass0 * Rk / lambda_[kk] + inis * ml.matrix(np.asarray(Psix[0])[:, KN:])) * la.inv(nu * Rk - Q)
                        for i in range(1, len(Psix)):
                            Psiy = nu * (lambda_[kk] * inis * ml.matrix(np.asarray(Psix[i])[:, KN:]) + Psiy * Rk) * la.inv(nu * Rk - Q)
                        res[xi] = float(np.sum(Psiy))
                    Ret.append(res)
                else:
                    raise ValueError("FluidPrioQueue: Unknown parameter " + str(opt))
        else:  # k == K, lowest priority
            for opt, val in _pairs(args):
                if opt == 'stMoms':
                    stMoms = np.array([math.factorial(i) * float(np.sum(ini * np.linalg.matrix_power(np.asarray(la.inv(-Km)), i + 1) @ np.asarray(clok))) for i in range(1, val + 1)])
                    Ret.append(stMoms)
                elif opt == 'stDistr':
                    pts = np.atleast_1d(val)
                    res = np.zeros(len(pts))
                    for xi, x in enumerate(pts):
                        res[xi] = float(np.sum(mass0 * Rk / lambda_[kk])) + float(np.sum(ini * la.inv(-Km) * (ml.eye(KN) - la.expm(np.asarray(Km) * x)) * clok))
                    Ret.append(res)
                elif opt == 'flMoms':
                    r = FluFluQueue(np.asarray(Q), np.asarray(Rk), np.array([[0.0]]), np.array([[float(d)]]), False, numFluidMoments=val)
                    Ret.append(np.asarray(r.fluidMoments).flatten())
                elif opt == 'flDistr':
                    pts = np.atleast_1d(val)
                    r = FluFluQueue(np.asarray(Q), np.asarray(Rk), np.array([[0.0]]), np.array([[float(d)]]), False, numFluidMoments=1)
                    y = FluidStationaryDistr(ml.matrix(r.fluid_mass0), ml.matrix(r.fluid_ini), ml.matrix(r.fluid_K), ml.matrix(r.fluid_clo), pts)
                    Ret.append(np.asarray(y).sum(axis=1).flatten())
                else:
                    raise ValueError("FluidPrioQueue: Unknown parameter " + str(opt))

    if len(Ret) == 1:
        return Ret[0]
    return Ret


def _pairs(args):
    """Yield (option, value) pairs from a flat argument list, skipping option keywords already consumed."""
    i = 0
    while i < len(args):
        opt = args[i]
        val = args[i + 1]
        yield opt, val
        i += 2
