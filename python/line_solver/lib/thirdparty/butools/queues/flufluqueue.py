# -*- coding: utf-8 -*-
"""
Fluid queue with independent input and output processes.

Standalone implementation that does not require external butools package.

Port from:


References:
    Horvath G, Telek M, "Sojourn times in fluid queues with independent and
    dependent input and output processes", PERFORMANCE EVALUATION 79: pp. 160-181, 2014.
"""

import numpy as np
import numpy.matlib as ml
import scipy.linalg as la
from dataclasses import dataclass
from typing import Optional, Tuple


@dataclass
class FluFluResult:
    """
    Result class for FluFluQueue containing ME distribution parameters
    and computed performance measures.
    """
    # Fluid level distribution parameters (mass0, ini, K, clo from GeneralFluidSolve)
    fluid_mass0: Optional[np.ndarray] = None
    fluid_ini: Optional[np.ndarray] = None
    fluid_K: Optional[np.ndarray] = None
    fluid_clo: Optional[np.ndarray] = None
    # Sojourn time distribution parameters
    sojourn_mass0: Optional[np.ndarray] = None
    sojourn_ini: Optional[np.ndarray] = None
    sojourn_K: Optional[np.ndarray] = None
    sojourn_clo: Optional[np.ndarray] = None
    # Fluid level moments (if requested)
    fluidMoments: Optional[np.ndarray] = None
    # Sojourn time moments (if requested)
    sojournMoments: Optional[np.ndarray] = None
    # Mean input rate (lambda)
    lambda_rate: float = 0.0
    # Mean output rate (mu)
    mu_rate: float = 0.0


def _ctmc_solve(Q: np.ndarray) -> np.ndarray:
    """
    Computes the stationary solution of a continuous time Markov chain.

    Parameters
    ----------
    Q : matrix, shape (M, M)
        The generator matrix of the Markov chain

    Returns
    -------
    pi : row vector, shape (1, M)
        The vector that satisfies pi * Q = 0, sum(pi) = 1
    """
    M = np.array(Q)
    M[:, 0] = np.ones(M.shape[0])
    m = np.zeros(M.shape[0])
    m[0] = 1.0
    return ml.matrix(la.solve(M.T, m))


def _diag(v: np.ndarray) -> np.ndarray:
    """Create diagonal matrix from vector."""
    return np.diag(np.asarray(v).flatten())


def _linsolve(A: np.ndarray, b: np.ndarray) -> np.ndarray:
    """
    Solve linear system Ax = b (if b is column vector) or xA = b (if b is row vector).
    Uses QR decomposition to handle non-square matrices.
    """
    A = np.asarray(A)
    b = np.asarray(b)

    if b.shape[0] == 1:
        # Row vector: solve xA = b, equivalent to A^T x^T = b^T
        x = _linsolve(np.conj(A.T), np.conj(b.T))
        return np.conj(x.T)
    elif len(b.shape) == 1 or b.shape[1] == 1:
        # Column vector: solve Ax = b using QR
        Q, R = la.qr(A)
        N = A.shape[1]
        b_flat = np.asarray(np.conj(Q.T) @ b).flatten()[:N]
        return ml.matrix(la.solve(R[:N, :N], b_flat)).T
    else:
        raise ValueError("b must be a row or column vector")


def _fluid_fundamental_matrices(
    Fpp: np.ndarray,
    Fpm: np.ndarray,
    Fmp: np.ndarray,
    Fmm: np.ndarray,
    matrices: str = "PKU",
    precision: float = 1e-14,
    maxNumIt: int = 50,
    method: str = "ADDA"
) -> Tuple:
    """
    Returns the fundamental matrices corresponding to the canonical Markov fluid model.

    Uses ADDA (Alternating-Directional Doubling Algorithm) by default.
    """
    Fpp = np.asmatrix(Fpp)
    Fpm = np.asmatrix(Fpm)
    Fmp = np.asmatrix(Fmp)
    Fmm = np.asmatrix(Fmm)

    # ADDA algorithm
    A = ml.matrix(-Fpp)
    B = Fpm
    C = Fmp
    D = ml.matrix(-Fmm)

    gamma1 = np.max(np.diag(A))
    gamma2 = np.max(np.diag(D))

    sA = A.shape[0]
    sD = D.shape[0]
    IA = ml.eye(sA)
    ID = ml.eye(sD)

    A = A + gamma2 * IA
    D = D + gamma1 * ID

    Dginv = D.I
    Vginv = (D - C * A.I * B).I
    Wginv = (A - B * Dginv * C).I

    Eg = ID - (gamma1 + gamma2) * Vginv
    Fg = IA - (gamma1 + gamma2) * Wginv
    Gg = (gamma1 + gamma2) * Dginv * C * Wginv
    Hg = (gamma1 + gamma2) * Wginv * B * Dginv

    diff = 1.0
    numit = 0

    while diff > precision and numit < maxNumIt:
        Vginv = Eg * (ID - Gg * Hg).I
        Wginv = Fg * (IA - Hg * Gg).I
        Gg = Gg + Vginv * Gg * Fg
        Hg = Hg + Wginv * Hg * Eg
        Eg = Vginv * Eg
        Fg = Wginv * Fg
        neg = la.norm(Eg, 1)
        nfg = la.norm(Fg, 1)
        # ADDA scaling
        eta = np.sqrt(nfg / neg) if neg > 0 else 1.0
        Eg = Eg * eta
        Fg = Fg / eta
        diff = neg * nfg
        numit += 1

    Psi = ml.matrix(Hg)

    ret = []
    for M in matrices:
        if M == "P":
            ret.append(Psi)
        elif M == "K":
            ret.append(Fpp + Psi * Fmp)
        elif M == "U":
            ret.append(Fmm + Fmp * Psi)

    if len(ret) == 1:
        return ret[0]
    else:
        return tuple(ret)


def _general_fluid_solve(
    Q: np.ndarray,
    R: np.ndarray,
    Q0=None,
    prec: float = 1e-14
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Returns the parameters of the matrix-exponentially distributed stationary
    distribution of a general Markovian fluid model.

    Parameters
    ----------
    Q : matrix, shape (N, N)
        The generator of the background Markov chain
    R : diagonal matrix, shape (N, N)
        The diagonal matrix of the fluid rates associated with the different states
    Q0 : matrix, shape (N, N), optional
        The generator at level 0. If not provided, Q0=Q is assumed.
    prec : float, optional
        Numerical precision (default 1e-14)

    Returns
    -------
    mass0 : matrix, shape (1, N)
        The stationary probability vector of zero level
    ini : matrix, shape (1, Np)
        The initial vector of the stationary density
    K : matrix, shape (Np, Np)
        The matrix parameter of the stationary density
    clo : matrix, shape (Np, N)
        The closing matrix of the stationary density
    """
    Q = np.asmatrix(Q)
    R = np.asmatrix(R)

    N = Q.shape[0]

    # Partition states by fluid rate sign: zero (ixz), positive (ixp), negative (ixn)
    ix = np.arange(N)
    Rdiag = np.diag(R)
    ixz = ix[np.abs(Rdiag) <= prec]
    ixp = ix[Rdiag > prec]
    ixn = ix[Rdiag < -prec]

    Nz = len(ixz)
    Np = len(ixp)
    Nn = len(ixn)

    # Permutation matrix
    P = ml.zeros((N, N))
    for i in range(Nz):
        P[i, ixz[i]] = 1
    for i in range(Np):
        P[Nz + i, ixp[i]] = 1
    for i in range(Nn):
        P[Nz + Np + i, ixn[i]] = 1

    iP = P.I
    Qv = P * Q * iP
    Rv = P * R * iP

    # Censor out zero-rate states
    if Nz > 0:
        iQv00 = la.pinv(-Qv[:Nz, :Nz])
        Qbar = Qv[Nz:, Nz:] + Qv[Nz:, :Nz] * iQv00 * Qv[:Nz, Nz:]
    else:
        Qbar = Qv

    absRi = _diag(np.abs(1.0 / np.diag(Rv[Nz:, Nz:])))
    Qz = absRi @ Qbar

    # Get fundamental matrices
    Psi, K, U = _fluid_fundamental_matrices(
        Qz[:Np, :Np],
        Qz[:Np, Np:],
        Qz[Np:, :Np],
        Qz[Np:, Np:],
        "PKU",
        prec
    )

    # Closing matrix
    Pm = np.hstack((ml.eye(Np), Psi)) @ absRi
    iCn = absRi[Np:, Np:]
    iCp = absRi[:Np, :Np]

    if Nz > 0:
        clo = np.hstack((
            (iCp @ Qv[Nz:Nz+Np, :Nz] + Psi @ iCn @ Qv[Nz+Np:, :Nz]) @ iQv00,
            Pm
        ))
    else:
        clo = Pm

    if Q0 is None or (isinstance(Q0, (list, np.ndarray)) and len(Q0) == 0):
        # Regular boundary behavior
        clo = clo @ P  # Go back to original state ordering

        # Calculate boundary vector
        if Nz > 0:
            Ua = (iCn @ Qv[Nz+Np:, :Nz] @ iQv00 @ ml.ones((Nz, 1)) +
                  iCn @ ml.ones((Nn, 1)) +
                  Qz[Np:, :Np] @ la.inv(-K) @ clo @ ml.ones((Nz+Np+Nn, 1)))
        else:
            Ua = iCn @ ml.ones((Nn, 1)) + Qz[Np:, :Np] @ la.inv(-K) @ clo @ ml.ones((Np+Nn, 1))

        pm = _linsolve(
            np.hstack((U, Ua)).T,
            np.hstack((ml.zeros((1, Nn)), ml.ones((1, 1)))).T
        ).T
        pm = ml.matrix(pm)

        # Create result
        if Nz > 0:
            mass0 = np.hstack((
                pm @ iCn @ Qv[Nz+Np:, :Nz] @ iQv00,
                ml.zeros((1, Np)),
                pm @ iCn
            )) @ P
        else:
            mass0 = np.hstack((ml.zeros((1, Np)), pm @ iCn)) @ P

        ini = pm @ Qz[Np:, :Np]
    else:
        # Custom boundary behavior at Q0
        Q0 = np.asmatrix(Q0)
        Q0v = P @ Q0 @ iP
        M = np.vstack((-clo @ Rv, Q0v[Nz+Np:, :], Q0v[:Nz, :]))
        Ma = np.vstack((np.sum(la.inv(-K) @ clo, axis=1), ml.ones((Nz+Nn, 1))))
        sol = _linsolve(
            np.hstack((M, Ma)).T,
            np.hstack((ml.zeros((1, N)), ml.ones((1, 1)))).T
        ).T
        sol = ml.matrix(sol)
        ini = sol[:, :Np]
        clo = clo @ P
        mass0 = np.hstack((sol[:, Np+Nn:], ml.zeros((1, Np)), sol[:, Np:Np+Nn])) @ P

    return mass0, ini, K, clo


def FluFluQueue(
    Qin: np.ndarray,
    Rin: np.ndarray,
    Qout: np.ndarray,
    Rout: np.ndarray,
    srv0stop: bool = True,
    numFluidMoments: int = 0,
    numSojournMoments: int = 0,
    prec: float = 1e-14
) -> FluFluResult:
    """
    Returns various performance measures of a fluid queue with independent
    fluid arrival and service processes.

    Two types of boundary behavior are available:
    - If srv0stop=False, the output process evolves continuously even if the queue is empty.
    - If srv0stop=True, the output process slows down if there is less fluid in the queue
      than it can serve. If the queue is empty and the fluid input rate is zero, the output
      process freezes until fluid arrives.

    Parameters
    ----------
    Qin : matrix, shape (Nin, Nin)
        Generator of the background Markov chain for the input process
    Rin : matrix, shape (Nin, Nin)
        Diagonal matrix of fluid input rates
    Qout : matrix, shape (Nout, Nout)
        Generator of the background Markov chain for the output process
    Rout : matrix, shape (Nout, Nout)
        Diagonal matrix of fluid output rates
    srv0stop : bool, optional
        If True, service slows down when queue is nearly empty. Default is True.
    numFluidMoments : int, optional
        Number of fluid level moments to compute (0 to skip). Default is 0.
    numSojournMoments : int, optional
        Number of sojourn time moments to compute (0 to skip). Default is 0.
    prec : float, optional
        Numerical precision. Default is 1e-14.

    Returns
    -------
    FluFluResult
        Result containing requested performance measures.
    """
    Qin = np.asmatrix(Qin)
    Rin = np.asmatrix(Rin)
    Qout = np.asmatrix(Qout)
    Rout = np.asmatrix(Rout)

    Nin = Qin.shape[0]
    Nout = Qout.shape[0]

    Iin = np.eye(Nin)
    Iout = np.eye(Nout)

    # Compute mean input and output rates for normalization
    piIn = _ctmc_solve(Qin)
    piOut = _ctmc_solve(Qout)
    lambda_rate = float(np.sum(np.asarray(piIn) * np.diag(Rin)))
    mu_rate = float(np.sum(np.asarray(piOut) * np.diag(Rout)))

    result = FluFluResult(lambda_rate=lambda_rate, mu_rate=mu_rate)

    # Compute fluid level distribution if needed
    if numFluidMoments > 0:
        # Q = kron(Qin, Iout) + kron(Iin, Qout)
        Q = np.kron(Qin, Iout) + np.kron(Iin, Qout)

        # R = kron(Rin, Iout) - kron(Iin, Rout)
        R = np.kron(Rin, Iout) - np.kron(Iin, Rout)

        # Check for degenerate case: all rates non-positive (stable system, no fluid buildup)
        Rdiag = np.diag(R)
        if np.all(Rdiag <= prec):
            # Stable system: service always >= arrival, fluid level is essentially 0
            result.fluidMoments = np.zeros(numFluidMoments)
            result.fluid_mass0 = np.ones((1, R.shape[0])) / R.shape[0]
            result.fluid_ini = np.zeros((1, 0))
            result.fluid_K = np.zeros((0, 0))
            result.fluid_clo = np.zeros((0, R.shape[0]))
        else:
            # Q0 depends on srv0stop
            if srv0stop:
                # Q0 = kron(Qin, Iout) + kron(Rin, pinv(Rout)*Qout)
                RoutPinv = la.pinv(Rout)
                Q0 = np.kron(Qin, Iout) + np.kron(Rin, RoutPinv @ Qout)
            else:
                Q0 = None

            mass0, ini, K, clo = _general_fluid_solve(Q, R, Q0, prec)
            result.fluid_mass0 = mass0
            result.fluid_ini = ini
            result.fluid_K = K
            result.fluid_clo = clo

            # Compute fluid level moments: E[X^m] = m! * ini * (-K)^{-(m+1)} * clo * ones
            negK = -K
            invNegK = la.inv(negK)
            onesN = np.ones((clo.shape[1], 1))

            fluidMoments = np.zeros(numFluidMoments)
            invKPower = np.asarray(invNegK).copy()  # (-K)^{-1}

            for m in range(1, numFluidMoments + 1):
                # (-K)^{-(m+1)} = (-K)^{-m} * (-K)^{-1}
                invKPower = invKPower @ np.asarray(invNegK)
                moment = _factorial(m) * np.sum(np.asarray(ini) @ invKPower @ np.asarray(clo) @ onesN)
                fluidMoments[m - 1] = moment

            result.fluidMoments = fluidMoments

    # Compute sojourn time distribution if needed
    if numSojournMoments > 0:
        # For sojourn time, use different formulation:
        # Rh = kron(Rin, Iout) - kron(Iin, Rout)
        # Qh = kron(Qin, Rout) + kron(Rin, Qout)
        Rh = np.kron(Rin, Iout) - np.kron(Iin, Rout)
        Qh = np.kron(Qin, Rout) + np.kron(Rin, Qout)

        # Check for degenerate case: all rates non-positive
        Rhdiag = np.diag(Rh)
        if np.all(Rhdiag <= prec):
            # Degenerate case: no fluid buildup, sojourn time = service time
            # E[W] = 1/mu for simple service
            result.sojournMoments = np.array([1.0 / mu_rate ** m for m in range(1, numSojournMoments + 1)])
            result.sojourn_mass0 = np.ones((1, Rh.shape[0])) / Rh.shape[0]
            result.sojourn_ini = np.zeros((1, 0))
            result.sojourn_K = np.zeros((0, 0))
            result.sojourn_clo = np.zeros((0, Rh.shape[0]))
        else:
            mass0h, inih, Kh, cloh = _general_fluid_solve(Qh, Rh, None, prec)
            result.sojourn_mass0 = mass0h
            result.sojourn_ini = inih
            result.sojourn_K = Kh
            result.sojourn_clo = cloh

            # kclo depends on srv0stop
            if srv0stop:
                # kclo = cloh * kron(Rin, Rout) / (lambda * mu)
                kclo = np.asarray(cloh) @ np.kron(Rin, Rout) / (lambda_rate * mu_rate)
            else:
                # kclo = cloh * kron(Rin, Iout) / lambda
                kclo = np.asarray(cloh) @ np.kron(Rin, Iout) / lambda_rate

            negKh = -Kh
            invNegKh = la.inv(negKh)

            sojournMoments = np.zeros(numSojournMoments)
            invKhPower = np.asarray(invNegKh).copy()  # (-Kh)^{-1}

            for m in range(1, numSojournMoments + 1):
                # (-Kh)^{-(m+1)} = (-Kh)^{-m} * (-Kh)^{-1}
                invKhPower = invKhPower @ np.asarray(invNegKh)
                moment = _factorial(m) * np.sum(np.asarray(inih) @ invKhPower @ kclo)
                sojournMoments[m - 1] = moment

            result.sojournMoments = sojournMoments

    return result


def _factorial(n: int) -> float:
    """Computes factorial of n."""
    result = 1.0
    for i in range(2, n + 1):
        result *= i
    return result


__all__ = ['FluFluQueue', 'FluFluResult']
