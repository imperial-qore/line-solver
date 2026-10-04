# -*- coding: utf-8 -*-
"""
Sojourn-time distribution of Markov-modulated fluid queues, returned as a
matrix-exponential (ME) or phase-type (PH) representation.

Ported from BUTools-family fluid tools (G. Horvath):
  * FluidQueueSTD  - fluid queue with input rate matrix Rin, output rate Rout
  * FluFluSTD      - fluid queue with fluid-modulated arrival and service

References
----------
Horvath G, Telek M, "Sojourn times in fluid queues with independent and
dependent input and output processes", Performance Evaluation 79:160-181, 2014.
"""
import numpy as np
import numpy.matlib as ml
import scipy.linalg as la

from ..mam.fluid import GeneralFluidSolve
from ..mc.stst import CTMCSolve
from ..reptrans.transform_ones import TransformToOnes



__all__ = ["FluidQueueSTD", "FluFluSTD"]

def _diag(v):
    """Diagonal matrix from a vector (row or column)."""
    return ml.matrix(np.diagflat(np.asarray(v).flatten()))


def FluidQueueSTD(Q, Rin, Rout, Q0=None, transToPH=False):
    """
    Sojourn-time distribution of a fluid queue with input rate matrix Rin and
    output (service) rate matrix Rout, modulated by generator Q.

    Returns (alpha, A): an ME (transToPH=False) or PH (transToPH=True)
    representation of the sojourn time.
    """
    Q = ml.matrix(Q)
    Rin = ml.matrix(Rin)
    Rout = ml.matrix(Rout)
    N = Q.shape[0]

    if Q0 is None:
        mass0, ini, K, clo = GeneralFluidSolve(Q, Rin - Rout)
    else:
        mass0, ini, K, clo = GeneralFluidSolve(Q, Rin - Rout, ml.matrix(Q0))
    mass0 = ml.matrix(mass0)
    ini = ml.matrix(ini)
    K = ml.matrix(K)
    clo = ml.matrix(clo)
    nk = K.shape[0]

    iniKi = ml.matrix(la.solve(K.T, -ini.T)).T          # ini*inv(-K)
    lambda_ = float(np.sum(mass0 * Rin + iniKi * clo * Rin))

    if transToPH:
        Delta = _diag(iniKi / lambda_)
        alpha = ml.matrix(np.asarray(clo * Rin).reshape(1, N * nk, order='F')) * ml.matrix(np.kron(ml.eye(N), Delta))
        A = ml.matrix(np.kron(Rout, la.inv(Delta) * K.T * Delta)) + ml.matrix(np.kron(Q, ml.eye(nk)))
    else:
        B = TransformToOnes(ml.matrix(np.asarray(la.inv(-K) * clo * Rin).reshape(N * nk, 1, order='F')))
        Bi = la.inv(B)
        alpha = ml.matrix(np.kron(ml.ones((1, N)), ini / lambda_)) * Bi
        A = B * (ml.matrix(np.kron(Q.T, ml.eye(nk))) + ml.matrix(np.kron(Rout, K))) * Bi
    return alpha, A


def FluFluSTD(Qin, Rin, Qout, Rout, srv0stop, transToPH=False):
    """
    Sojourn-time distribution of a fluid queue in which both the arrival and the
    service processes are Markov-modulated fluid flows. If srv0stop is True the
    service stops while the server fluid level is zero.

    Returns (alpha, A): an ME (transToPH=False) or PH (transToPH=True)
    representation of the sojourn time.
    """
    Qin = ml.matrix(Qin)
    Rin = ml.matrix(Rin)
    Qout = ml.matrix(Qout)
    Rout = ml.matrix(Rout)
    Iin = ml.eye(Qin.shape[0])
    Iout = ml.eye(Qout.shape[0])

    Rh = ml.matrix(np.kron(Rin, Iout)) - ml.matrix(np.kron(Iin, Rout))
    Qh = ml.matrix(np.kron(Qin, Rout)) + ml.matrix(np.kron(Rin, Qout))
    massh, inih, Kh, cloh = GeneralFluidSolve(Qh, Rh)
    inih = ml.matrix(inih)
    Kh = ml.matrix(Kh)
    cloh = ml.matrix(cloh)

    lambda_ = float(np.sum(ml.matrix(CTMCSolve(Qin)) * Rin))
    mu = float(np.sum(ml.matrix(CTMCSolve(Qout)) * Rout))

    if transToPH:
        Delta = _diag(ml.matrix(la.solve(Kh.T, -inih.T)))   # diag(inih*inv(-Kh))
        A = la.inv(Delta) * Kh.T * Delta
        if not srv0stop:
            alpha = np.asarray((Delta * cloh * ml.matrix(np.kron(Rin, Iout)) / lambda_).sum(axis=1)).reshape(1, -1)
        else:
            alpha = np.asarray((Delta * cloh * ml.matrix(np.kron(Rin, Rout)) / lambda_ / mu).sum(axis=1)).reshape(1, -1)
        alpha = ml.matrix(alpha)
    else:
        if not srv0stop:
            B = TransformToOnes((cloh * ml.matrix(np.kron(Rin, Iout)) / lambda_).sum(axis=1))
        else:
            B = TransformToOnes((cloh * ml.matrix(np.kron(Rin, Rout)) / lambda_ / mu).sum(axis=1))
        iB = la.inv(B)
        A = B * Kh * iB
        alpha = inih * la.inv(-Kh) * iB
    return alpha, A
