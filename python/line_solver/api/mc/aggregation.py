"""
CTMC Aggregation-Disaggregation Methods.

Native Python implementations for near-completely decomposable (NCD)
Markov chain analysis using aggregation-disaggregation techniques.

Key algorithms:
    ctmc_courtois: Courtois decomposition
    ctmc_kms: Koury-McAllister-Stewart aggregation-disaggregation
    ctmc_takahashi: Takahashi's aggregation-disaggregation
    ctmc_multi: Multigrid aggregation-disaggregation

References:
    Original MATLAB: matlab/src/api/mc/ctmc_*.m
    Courtois, "Decomposability: Queueing and Computer System Applications", 1977
"""

import numpy as np
from numpy.linalg import LinAlgError
from scipy import linalg
from typing import Dict, Any, Optional, List, Tuple
from dataclasses import dataclass

from .ctmc import GMRES_MIN_STATES, ctmc_makeinfgen, ctmc_solve_reducible
from .dtmc import dtmc_solve, dtmc_makestochastic


@dataclass
class CourtoisResult:
    """Result of Courtois decomposition."""
    p: np.ndarray  # Approximate steady-state probability vector
    Qperm: np.ndarray  # Q reordered according to macrostates
    Qdec: np.ndarray  # Infinitesimal generator for macrostates
    eps: float  # Nearly-complete decomposability (NCD) index
    epsMAX: float  # Max acceptable value for eps
    P: np.ndarray  # Probability matrix from randomization
    B: np.ndarray  # Part of P not modelled by decomposition
    q: float  # Randomization coefficient


@dataclass
class KMSResult:
    """Result of KMS aggregation-disaggregation."""
    p: np.ndarray  # Estimated steady-state probability vector
    p_1: np.ndarray  # Previous iteration
    Qperm: np.ndarray  # Permuted Q matrix
    eps: float  # NCD index
    epsMAX: float  # Max acceptable eps
    pcourt: np.ndarray  # Courtois solution


@dataclass
class TakahashiResult:
    """Result of Takahashi aggregation-disaggregation."""
    p: np.ndarray  # Estimated steady-state probability vector
    p_1: np.ndarray  # Previous iteration
    pcourt: np.ndarray  # Courtois solution
    Qperm: np.ndarray  # Permuted Q matrix
    eps: float  # NCD index
    epsMAX: float  # Max acceptable eps


def _ctmc_randomization(Q: np.ndarray, q: Optional[float] = None) -> Tuple[np.ndarray, float]:
    """
    Uniformization (randomization) of a CTMC.

    Converts CTMC generator to uniformized DTMC transition matrix.

    Args:
        Q: Infinitesimal generator matrix
        q: Uniformization rate. If None, derived as 1.05*max|Q|, matching the
           default of MATLAB ctmc_courtois.m. It must be a valid uniformization
           rate (q >= max|q_ii|) or P = I + Q/q is not stochastic and the
           dtmc_makestochastic call below silently patches a negative diagonal.

    Returns:
        Tuple of (P, q) where P is the uniformized transition matrix
    """
    Q = np.asarray(Q, dtype=np.float64)

    if q is None:
        # Deterministic. The former default max|Q| + np.random.rand() drew an
        # unseeded random rate, making every caller irreproducible.
        q = 1.05 * np.max(np.abs(Q))

    P = Q / q + np.eye(Q.shape[0])
    P = dtmc_makestochastic(P)

    return P, q


def ctmc_courtois(Q: np.ndarray, MS: List[List[int]],
                  q: Optional[float] = None) -> CourtoisResult:
    """
    Courtois decomposition for near-completely decomposable CTMCs.

    Decomposes a large CTMC into macrostates and computes approximate
    steady-state probabilities using hierarchical aggregation.

    Args:
        Q: Infinitesimal generator matrix
        MS: List where MS[i] is the list of state indices in macrostate i
        q: Randomization coefficient (optional)

    Returns:
        CourtoisResult with approximate solution and diagnostics

    References:
        Original MATLAB: matlab/src/api/mc/ctmc_courtois.m
        Courtois, "Decomposability: Queueing and Computer System Applications", 1977
    """
    Q = np.asarray(Q, dtype=np.float64)
    nStates = Q.shape[0]
    nMacroStates = len(MS)

    # Rearrange generator according to macrostates
    v = []
    for states in MS:
        v.extend(states)
    v = np.array(v, dtype=int)

    Qperm = Q[np.ix_(v, v)]
    Qdec = Qperm.copy()

    # Zero out off-diagonal blocks
    procRows = 0
    for i in range(nMacroStates):
        block_size = len(MS[i])
        if procRows > 0:
            Qdec[procRows:procRows + block_size, :procRows] = 0
        Qdec[procRows:procRows + block_size, procRows + block_size:] = 0
        procRows += block_size

    # Make each diagonal block a valid generator
    Qdec = ctmc_makeinfgen(Qdec)

    # see _kb/03-api-layer.md for rationale
    if q is None:
        q = 1.05 * np.max(np.abs(Qperm))
    P, q = _ctmc_randomization(Qperm, q)

    # Compute block-diagonal part A
    A = P.copy()
    procRows = 0
    for i in range(nMacroStates):
        block_size = len(MS[i])
        if procRows > 0:
            A[procRows:procRows + block_size, :procRows] = 0
        A[procRows:procRows + block_size, procRows + block_size:] = 0
        procRows += block_size

    B = P - A
    # see _kb/03-api-layer.md for rationale
    eps = np.max(np.sum(B, axis=1))

    # Compute epsMAX
    procRows = 0
    A_norm = A.copy()
    for i in range(nMacroStates):
        block_size = len(MS[i])
        for j in range(block_size):
            pos = j
            row_idx = procRows + j
            col_idx = procRows + pos
            other_cols = list(range(procRows, procRows + block_size))
            other_cols.remove(col_idx)
            A_norm[row_idx, col_idx] = 1 - np.sum(A_norm[row_idx, other_cols])
        procRows += block_size

    eigMS = np.zeros(nMacroStates)
    procRows = 0
    for i in range(nMacroStates):
        block_size = len(MS[i])
        block = A_norm[procRows:procRows + block_size, procRows:procRows + block_size]
        e = np.sort(np.abs(linalg.eigvals(block)))
        if len(e) > 1:
            eigMS[i] = e[-2]  # Second largest eigenvalue
        procRows += block_size

    epsMAX = (1 - np.max(eigMS)) / 2

    # Compute microprobabilities
    pmicro = np.zeros(nStates)
    procRows = 0
    for i in range(nMacroStates):
        block_size = len(MS[i])
        Qmicrostate = Qdec[procRows:procRows + block_size, procRows:procRows + block_size]
        pmicro[procRows:procRows + block_size] = ctmc_solve_reducible(Qmicrostate)
        procRows += block_size

    # Compute macroprobabilities
    G = np.zeros((nMacroStates, nMacroStates))
    procRows = 0
    for i in range(nMacroStates):
        procCols = 0
        for j in range(nMacroStates):
            if i != j:
                for iState in range(len(MS[i])):
                    row_idx = procRows + iState
                    col_slice = slice(procCols, procCols + len(MS[j]))
                    G[i, j] += pmicro[row_idx] * np.sum(P[row_idx, col_slice])
            procCols += len(MS[j])
        procRows += len(MS[i])

    for i in range(nMacroStates):
        G[i, i] = 1 - np.sum(G[i, :])

    pMacro = dtmc_solve(G)

    # Combine micro and macro probabilities
    p_perm = np.zeros(nStates)
    procRows = 0
    for i in range(nMacroStates):
        block_size = len(MS[i])
        p_perm[procRows:procRows + block_size] = pMacro[i] * pmicro[procRows:procRows + block_size]
        procRows += block_size

    # Reorder back to original state ordering
    p = np.zeros(nStates)
    for i, orig_idx in enumerate(v):
        p[orig_idx] = p_perm[i]

    return CourtoisResult(
        p=p,
        Qperm=Qperm,
        Qdec=Qdec,
        eps=eps,
        epsMAX=epsMAX,
        P=P,
        B=B,
        q=q
    )


def ctmc_kms(Q: np.ndarray, MS: List[List[int]],
             numSteps: int = 10) -> KMSResult:
    """
    Koury-McAllister-Stewart aggregation-disaggregation method.

    Iteratively refines the Courtois decomposition solution using
    aggregation and disaggregation steps.

    Args:
        Q: Infinitesimal generator matrix
        MS: List where MS[i] is the list of state indices in macrostate i
        numSteps: Number of iterative steps (default: 10)

    Returns:
        KMSResult with refined solution

    References:
        Original MATLAB: matlab/src/api/mc/ctmc_kms.m
        Koury, McAllister, Stewart, "Iterative Methods for Computing
        Stationary Distributions of Nearly Completely Decomposable
        Markov Chains", 1984
    """
    Q = np.asarray(Q, dtype=np.float64)
    nStates = Q.shape[0]
    nMacroStates = len(MS)

    # Start from Courtois decomposition solution
    courtois = ctmc_courtois(Q, MS)
    pcourt = courtois.p
    Qperm = courtois.Qperm
    eps = courtois.eps
    epsMAX = courtois.epsMAX
    P = courtois.P

    # Initial solution
    pn = pcourt.copy()

    # Reorder pn according to macrostate ordering
    v = []
    for states in MS:
        v.extend(states)
    v = np.array(v, dtype=int)
    pn_perm = pn[v]
    p_1 = pn_perm.copy()  # defined even when numSteps == 0

    # Main iteration loop
    for n in range(numSteps):
        pn_1 = pn_perm.copy()
        p_1 = pn_1.copy()

        # Aggregation step
        pcondn_1 = pn_1.copy()
        procRows = 0
        for I in range(nMacroStates):
            block_size = len(MS[I])
            block_sum = np.sum(pn_1[procRows:procRows + block_size])
            if block_sum > 0:
                pcondn_1[procRows:procRows + block_size] = (
                    pn_1[procRows:procRows + block_size] / block_sum
                )
            procRows += block_size

        G = np.zeros((nMacroStates, nMacroStates))
        procCols = 0
        for I in range(nMacroStates):
            procRows = 0
            for J in range(nMacroStates):
                block_J = len(MS[J])
                block_I = len(MS[I])
                pcond_J = pcondn_1[procRows:procRows + block_J]
                P_block = P[procRows:procRows + block_J, procCols:procCols + block_I]
                G[I, J] = pcond_J @ P_block @ np.ones(block_I)
                procRows += block_J
            procCols += len(MS[I])

        w = dtmc_solve(G.T)

        # Disaggregation step
        z = pcondn_1.copy()
        L = np.zeros((nStates, nStates))
        D = np.zeros((nStates, nStates))
        U = np.zeros((nStates, nStates))

        procRows = 0
        zn = np.zeros(nStates)
        for I in range(nMacroStates):
            block_I = len(MS[I])
            zn[procRows:procRows + block_I] = w[I] * z[procRows:procRows + block_I]

            procCols = 0
            for J in range(nMacroStates):
                block_J = len(MS[J])
                if I > J:
                    L[procRows:procRows + block_I, procCols:procCols + block_J] = (
                        P[procRows:procRows + block_I, procCols:procCols + block_J]
                    )
                if I == J:
                    D[procRows:procRows + block_I, procCols:procCols + block_J] = (
                        np.eye(block_I) -
                        P[procRows:procRows + block_I, procCols:procCols + block_J]
                    )
                if I < J:
                    U[procRows:procRows + block_I, procCols:procCols + block_J] = (
                        P[procRows:procRows + block_I, procCols:procCols + block_J]
                    )
                procCols += block_J
            procRows += block_I

        # see _kb/03-api-layer.md for rationale
        M = D - U
        rhs = (zn @ L).T
        pn_perm = None
        # Above the dispatch threshold the dense solve is what limits the model
        # rather than the arithmetic. The direct solve remains the fallback.
        if M.shape[0] > GMRES_MIN_STATES:
            from .gmres import ctmc_gmres
            x_gmres, gflag, _, _ = ctmc_gmres(M.T, rhs)
            if gflag == 0:
                pn_perm = x_gmres.T
        if pn_perm is None:
            pn_perm = linalg.solve(M.T, rhs).T
        pn_perm = pn_perm / np.sum(pn_perm)

    # Reorder back to original state ordering
    p = np.zeros(nStates)
    for i, orig_idx in enumerate(v):
        p[orig_idx] = pn_perm[i]

    # Map p_1 back to the original ordering too, matching MATLAB ctmc_kms.m
    p_1_out = np.zeros(nStates)
    for i, orig_idx in enumerate(v):
        p_1_out[orig_idx] = p_1[i]

    return KMSResult(
        p=p,
        p_1=p_1_out,
        Qperm=Qperm,
        eps=eps,
        epsMAX=epsMAX,
        pcourt=pcourt
    )


def ctmc_takahashi(Q: np.ndarray, MS: List[List[int]],
                   numSteps: int = 10) -> TakahashiResult:
    """
    Takahashi's aggregation-disaggregation method.

    Iteratively refines the Courtois decomposition solution using
    a different aggregation-disaggregation scheme.

    Args:
        Q: Infinitesimal generator matrix
        MS: List where MS[i] is the list of state indices in macrostate i
        numSteps: Number of iterative steps (default: 10)

    Returns:
        TakahashiResult with refined solution

    References:
        Original MATLAB: matlab/src/api/mc/ctmc_takahashi.m
        Takahashi, "A Lumping Method for Numerical Calculations of
        Stationary Distributions of Markov Chains", 1975
    """
    Q = np.asarray(Q, dtype=np.float64)
    nStates = Q.shape[0]
    nMacroStates = len(MS)

    # Start from Courtois decomposition
    courtois = ctmc_courtois(Q, MS)
    pcourt = courtois.p
    Qperm = courtois.Qperm
    eps = courtois.eps
    epsMAX = courtois.epsMAX

    # see _kb/03-api-layer.md for rationale
    P, _ = _ctmc_randomization(Q, 1.05 * np.max(np.abs(Q)))

    # Reorder according to macrostates
    v = []
    for states in MS:
        v.extend(states)
    v = np.array(v, dtype=int)

    pn = pcourt.copy()

    # Main loop
    for n in range(numSteps):
        pn_1 = pn.copy()
        p_1 = pn_1.copy()

        # Aggregation step
        G = np.zeros((nMacroStates, nMacroStates))
        procRows = 0
        for I in range(nMacroStates):
            S = np.sum(pn_1[MS[I]])
            procCols = 0
            for J in range(nMacroStates):
                if I != J:
                    for i in MS[I]:
                        for j in MS[J]:
                            if S > 1e-14:
                                G[I, J] += P[i, j] * pn_1[i] / S
                procCols += len(MS[J])
            procRows += len(MS[I])

        for i in range(nMacroStates):
            G[i, i] = 1 - np.sum(G[i, :])

        gamma = dtmc_solve(G)

        # Disaggregation step
        GI = np.zeros((nMacroStates, nStates))
        for I in range(nMacroStates):
            S = np.sum(pn_1[MS[I]])
            if S > 1e-14:
                for j in range(nStates):
                    for i in MS[I]:
                        GI[I, j] += P[i, j] * pn_1[i] / S

        for I in range(nMacroStates):
            block_size = len(MS[I])
            A = np.eye(block_size)
            b = np.zeros(block_size)

            for i_local, i in enumerate(MS[I]):
                for j_local, j in enumerate(MS[I]):
                    A[i_local, j_local] -= P[j, i]

                for K in range(nMacroStates):
                    if K != I:
                        b[i_local] += gamma[K] * GI[K, i]

            try:
                pn_block = None
                if A.shape[0] > GMRES_MIN_STATES:
                    from .gmres import ctmc_gmres
                    x_gmres, gflag, _, _ = ctmc_gmres(A, b)
                    if gflag == 0:
                        pn_block = x_gmres
                if pn_block is None:
                    pn_block = linalg.solve(A, b)
                for i_local, i in enumerate(MS[I]):
                    pn[i] = pn_block[i_local]
            except LinAlgError:
                pass

        # Normalize
        if np.sum(pn) > 0:
            pn /= np.sum(pn)

    return TakahashiResult(
        p=pn,
        p_1=p_1,
        pcourt=pcourt,
        Qperm=Qperm,
        eps=eps,
        epsMAX=epsMAX
    )


def ctmc_multi(Q: np.ndarray, MS: List[List[int]],
               MSS: List[List[int]]) -> TakahashiResult:
    """
    Multigrid aggregation-disaggregation method.

    Two-level hierarchical decomposition using nested macrostates.

    Args:
        Q: Infinitesimal generator matrix
        MS: List where MS[i] is the list of state indices in macrostate i
        MSS: List where MSS[i] is the list of macrostate indices in macro-macrostate i

    Returns:
        TakahashiResult with multigrid solution

    References:
        Original MATLAB: matlab/src/api/mc/ctmc_multi.m
    """
    Q = np.asarray(Q, dtype=np.float64)
    nStates = Q.shape[0]
    nMacroStates = len(MS)

    # Rearrange generator according to macrostates
    v = []
    for states in MS:
        v.extend(states)
    v = np.array(v, dtype=int)

    Qperm = Q[np.ix_(v, v)]
    Qdec = Qperm.copy()

    # Zero out off-diagonal blocks
    procRows = 0
    for i in range(nMacroStates):
        block_size = len(MS[i])
        if procRows > 0:
            Qdec[procRows:procRows + block_size, :procRows] = 0
        Qdec[procRows:procRows + block_size, procRows + block_size:] = 0
        procRows += block_size

    Qdec = ctmc_makeinfgen(Qdec)

    # Apply randomization
    q = 1.05 * np.max(np.abs(Qperm))
    P, q = _ctmc_randomization(Qperm, q)

    # Compute block-diagonal part A
    A = P.copy()
    procRows = 0
    for i in range(nMacroStates):
        block_size = len(MS[i])
        if procRows > 0:
            A[procRows:procRows + block_size, :procRows] = 0
        A[procRows:procRows + block_size, procRows + block_size:] = 0
        procRows += block_size

    B = P - A
    # NCD degree of coupling is the induced infinity norm of the coupling
    # matrix, i.e. the largest ROW sum (see ctmc_courtois).
    eps = np.max(np.sum(B, axis=1))
    epsMAX = 0.0

    # Compute microprobabilities
    pmicro = np.zeros(nStates)
    procRows = 0
    for i in range(nMacroStates):
        block_size = len(MS[i])
        Qmicrostate = Qdec[procRows:procRows + block_size, procRows:procRows + block_size]
        pmicro[procRows:procRows + block_size] = ctmc_solve_reducible(Qmicrostate)
        procRows += block_size

    # Compute macroprobabilities
    G = np.zeros((nMacroStates, nMacroStates))
    procRows = 0
    for i in range(nMacroStates):
        procCols = 0
        for j in range(nMacroStates):
            if i != j:
                for iState in range(len(MS[i])):
                    row_idx = procRows + iState
                    col_slice = slice(procCols, procCols + len(MS[j]))
                    G[i, j] += pmicro[row_idx] * np.sum(P[row_idx, col_slice])
            procCols += len(MS[j])
        procRows += len(MS[i])

    for i in range(nMacroStates):
        G[i, i] = 1 - np.sum(G[i, :])

    # Use Courtois on the macro level with MSS
    pMacro_result = ctmc_courtois(G, MSS)
    pMacro = pMacro_result.p

    # Combine micro and macro probabilities
    p_perm = np.zeros(nStates)
    procRows = 0
    for i in range(nMacroStates):
        block_size = len(MS[i])
        p_perm[procRows:procRows + block_size] = pMacro[i] * pmicro[procRows:procRows + block_size]
        procRows += block_size

    # Reorder back to original state ordering
    p = np.zeros(nStates)
    for i, orig_idx in enumerate(v):
        p[orig_idx] = p_perm[i]

    pcourt = ctmc_courtois(Q, MS).p

    return TakahashiResult(
        p=p,
        p_1=np.zeros(nStates),
        pcourt=pcourt,
        Qperm=Qperm,
        eps=eps,
        epsMAX=epsMAX
    )


__all__ = [
    'CourtoisResult',
    'KMSResult',
    'TakahashiResult',
    'ctmc_courtois',
    'ctmc_kms',
    'ctmc_takahashi',
    'ctmc_multi',
]
