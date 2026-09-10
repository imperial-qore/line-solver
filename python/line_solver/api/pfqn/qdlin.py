"""QD-LIN: the Linearizer arm of AMVA-LD, on a plain demand matrix.

Array-level twin of what ``SolverMVA`` computes for ``method='qdlin'``: the
Linearizer of Chandy and Neuse, Commun. ACM 25(2), 1982, run inside the
queue-dependent AMVA framework of Casale, Perez and Wang (IFIP PERFORMANCE
2015), so the load-dependent term ``g_k`` is evaluated at the CORRECTED
arrival-instant queue rather than at the plain one.

This is a transcription of ``solver_amvald`` (MATLAB
``matlab/src/solvers/MVA/solver_amvald.m`` and its forward evaluation,
mirrored by ``line_solver.api.solvers.mva.amvald``) restricted to the domain a
demand matrix describes: closed classes only, one chain per class, unit visits,
PS queueing stations and one optional delay carrying ``Z``. Within that domain
it reproduces ``SolverMVA(model, 'qdlin')`` to machine precision, which is what
this kernel is for; it is NOT an independent re-derivation of the method.

TWO PROPERTIES OF THE REFERENCE ARE REPRODUCED DELIBERATELY, not inherited by
accident, and a caller comparing against a textbook Linearizer will see both:

  1. THE GAMMA CORRECTION IS CLASS-AGGREGATE, STORED IN SLICE 0. ``solver_amvald``
     allocates the (K, M, K) per-class Linearizer array for ``qdlin`` but writes
     the class-aggregate correction into it with a two-subscript assignment,
     ``gamma(s,k) = sum_r Q_s(k,r)/(Nt-1) - sum_r Q(k,r)/Nt``, which MATLAB
     linear-indexes to ``(s,k,1)``. Slices 1..K-1 stay zero while every reader
     indexes gamma per class. The correction that reaches the residence time is
     therefore ``N_0*gamma(r,k,0) - [r==0]*gamma(r,k,0)``: the aggregate
     correction scaled by the population of CHAIN 0 alone, with the self term
     removed only for chain 0. It coincides with the queue-dependent AMVA form
     ``(Nt-1)*gamma_agg`` iff K = 1, so single-chain models are unaffected and
     multichain ones are not. ``method='lin'`` takes the per-class form instead.
  2. A SINGLE-SERVER STATION STILL CARRIES A SOFTMIN TERM. The multiserver
     factor is ``pfqn_lldfun(1 + arrival-instant total, None, nservers)``, whose
     softmin at c = 1 is not exactly 1, so ``qdlin`` does not reduce to a
     textbook single-server AMVA even when every station has one server.

MU AND NSERVERS ARE DIFFERENT MECHANISMS, unlike in :func:`pfqn_qdamva`, which
folds the multiserver curve into ``mu``. Here ``mu`` is ``sn.lldscaling``, an
interpolated rate multiplier per station, and ``nservers`` is the server count
feeding the softmin term. A c-server station is ``nservers[k] = c``, NOT a
``mu`` row of ``minimum(1..smax, c)``; passing the latter reproduces
``Queue.setLoadDependence``, which is a different station.

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

from typing import Optional, Tuple

import numpy as np

from .utils import pfqn_lldfun

__all__ = ['pfqn_qdlin']

_OMICRON = 0.5  # under-relaxation parameter of solver_amvald


def _forward(ST, srv, isdelay, mu, gamma, Qin, Nin, K, wtol):
    """One forward evaluation, solver_amvald_forward restricted to PS/INF.

    Returns the waiting times W (Ms x K) and the effective service times STeff.
    """
    Ms = ST.shape[0]
    nnz = np.where(Nin > 0)[0]
    Nt_in = float(np.sum(Nin))
    delta_in = (Nt_in - 1.0) / Nt_in if Nt_in > 0 else 1.0
    dcl = np.where(Nin > 0, (Nin - 1.0) / np.where(Nin > 0, Nin, 1.0), 1.0)

    # arrival-instant queue lengths, class-aggregate and per class. The row sum
    # is taken one station at a time, as the reference does: a 2-D axis
    # reduction rounds differently from a 1-D one and the fixed point amplifies
    # the ulp into a tolerance-sized gap.
    interp = np.zeros(Ms)
    totArvl = np.zeros((Ms, K))
    for k in range(Ms):
        sumQk = float(np.sum(Qin[k, nnz]))
        interp[k] = delta_in * sumQk
        for r in nnz:
            totArvl[k, r] = dcl[r] * Qin[k, r] + sumQk - Qin[k, r]

    # lld term, evaluated at the gamma-corrected arrival-instant queue
    lldterm = np.ones((Ms, K))
    if nnz.size > 0:
        for r in nnz:
            gcorr = Nin[nnz] @ gamma[r, :, nnz] - gamma[r, :, r]
            lldterm[:, r] = pfqn_lldfun(1.0 + interp + gcorr, mu)
    else:
        lldterm[:] = pfqn_lldfun(1.0 + interp, mu)[:, None]

    # multiserver term; config 'default' leaves PS on the softmin arm
    if nnz.size > 0 and Nt_in > 0:
        g = np.zeros((nnz.size, Ms))
        for r in nnz:
            g = g + ((Nt_in - 1.0) / Nt_in) * Nin[r] * gamma[nnz, :, r]
        msterm = pfqn_lldfun(1.0 + interp + np.mean(g, axis=0), None, srv)
    else:
        msterm = pfqn_lldfun(1.0 + interp, None, srv)

    STeff = np.zeros((Ms, K))
    for r in nnz:
        STeff[:, r] = ST[:, r] * lldterm[:, r] * msterm

    W = np.zeros((Ms, K))
    for r in nnz:
        for k in range(Ms):
            if isdelay[k]:
                W[k, r] = STeff[k, r]
            else:
                corr = Nin[nnz] @ gamma[r, k, nnz] - gamma[r, k, r]
                # the floor is the reference's, see the module header
                W[k, r] = STeff[k, r] * max(wtol, 1.0 + totArvl[k, r] + corr)
    return W, STeff


def pfqn_qdlin(L: np.ndarray,
               N: np.ndarray,
               Z: Optional[np.ndarray] = None,
               mu: Optional[np.ndarray] = None,
               nservers: Optional[np.ndarray] = None,
               tol: float = 1e-6,
               maxiter: int = 1000,
               wtol: float = 1e-4
               ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    """QD-LIN on a closed multiclass network, matching SolverMVA's 'qdlin'.

    Args:
        L: (M x R) service demand matrix, queueing stations only.
        N: (R,) population vector, finite.
        Z: (R,) think time vector; a delay station carrying it is appended to
            the station list when any entry is positive, exactly as the
            equivalent ``Network`` would hold one. None means no think time.
        mu: (M x smax) load-dependent rate multipliers, ``sn.lldscaling``. A
            station whose row is constant is skipped by ``pfqn_lldfun``, so an
            ordinary station is a row of ones or simply None.
        nservers: (M,) server counts; None means one server everywhere.
        tol: convergence tolerance on the queue lengths. Defaults to LINE's own
            ``iter_tol``, so the kernel matches SolverMVA as called by default.
        maxiter: iteration budget, LINE's ``iter_max``. The outer sweep and each
            inner sweep are capped at sqrt(maxiter) and the total number of
            forward evaluations at min(maxiter, 10000), as in solver_amvald.
        wtol: floor on the AMVA wait factor, LINE's ``options.tol``. This is a
            DIFFERENT knob from the convergence tolerance and keeps its own
            default: SolverMVA passes iter_tol through to the fixed point but
            never sets options.tol, so the floor stays at the lineDefaults 1e-4
            while the fixed point converges to 1e-6. The floor is load-bearing
            for qdlin, whose class-aggregate correction drives the wait factor
            negative at a lightly loaded station.

    Returns:
        Q: (M x R) mean queue lengths at the queueing stations.
        U: (M x R) per-class utilizations.
        R: (M x R) per-class residence times.
        X: (1 x R) per-class throughputs.
        C: (1 x R) per-class cycle times, think time included.
        iter: number of forward evaluations performed.
    """
    L = np.atleast_2d(np.asarray(L, dtype=float))
    M, K = L.shape
    N = np.asarray(N, dtype=float).reshape(-1)
    if N.size != K:
        raise ValueError('pfqn_qdlin: the population vector must have one entry per class')
    if np.any(np.isinf(N)):
        raise ValueError('pfqn_qdlin: an infinite population is not supported, closed classes only')
    if Z is None:
        Z = np.zeros(K)
    Z = np.asarray(Z, dtype=float).reshape(-1)
    if Z.size != K:
        raise ValueError('pfqn_qdlin: the think-time vector must have one entry per class')
    if nservers is None:
        nservers = np.ones(M)
    nservers = np.asarray(nservers, dtype=float).reshape(-1)
    if nservers.size != M:
        raise ValueError('pfqn_qdlin: the server-count vector must have one entry per station')
    if mu is not None:
        mu = np.atleast_2d(np.asarray(mu, dtype=float))

    Nt = float(np.sum(N))
    if Nt <= 0:
        z = np.zeros((M, K))
        return z, z.copy(), z.copy(), np.zeros((1, K)), np.zeros((1, K)), 0

    # station list: the delay, when there is one, then the queueing stations
    has_delay = bool(np.any(Z > 0))
    if has_delay:
        ST = np.vstack([Z[None, :], L])
        srv = np.concatenate([[np.inf], nservers])
        isdelay = np.zeros(M + 1, dtype=bool)
        isdelay[0] = True
        mu_full = None if mu is None else np.vstack([np.ones((1, mu.shape[1])), mu])
    else:
        ST = L.copy()
        srv = nservers.copy()
        isdelay = np.zeros(M, dtype=bool)
        mu_full = mu
    Ms = ST.shape[0]

    # balanced initialization, as in solver_amvald
    Q = np.ones((Ms, K))
    Q = Q / np.sum(Q, axis=0, keepdims=True) * N.reshape(1, -1)
    Q[:, N == 0] = 0.0
    with np.errstate(divide='ignore'):
        X = 1.0 / np.sum(ST, axis=0)
    X[~np.isfinite(X)] = 0.0
    nnzclasses = np.where(N > 0)[0]
    U = np.zeros((Ms, K))
    for k in range(Ms):
        for r in nnzclasses:
            U[k, r] = ST[k, r] * X[r] if np.isinf(srv[k]) else ST[k, r] * X[r] / srv[k]

    V = np.ones((Ms, K))  # unit visits; sliced per class to match the reference's ddot
    gamma = np.zeros((K, Ms, K))
    T = np.zeros((Ms, K))
    C = np.zeros(K)
    STeff = np.zeros((Ms, K))

    max_sweep = np.sqrt(maxiter)
    max_totiter = min(maxiter, 10000)
    totiter = 0
    outer_iter = 0
    Q_outer_1 = Q + np.inf

    while (outer_iter < 2 or np.max(np.abs(Q - Q_outer_1)) > tol) \
            and outer_iter < max_sweep and totiter <= max_totiter:
        outer_iter += 1
        Q_outer_1 = Q.copy()
        X_outer_1 = X.copy()

        # Linearizer recursion: one sweep at each reduced population N - 1_s
        for s in range(K):
            if N[s] <= 0:
                continue
            N_s = N.copy()
            N_s[s] -= 1.0
            scale = (Nt - 1.0) / Nt
            Q_s, X_s = Q * scale, X * scale

            iter_s = 0
            Q_s_1 = Q_s + np.inf
            while (iter_s < 2 or np.max(np.abs(Q_s - Q_s_1)) > tol) and iter_s <= max_sweep:
                iter_s += 1
                Q_s_1, X_s_1 = Q_s.copy(), X_s.copy()

                W_s, _ = _forward(ST, srv, isdelay, mu_full, gamma, Q_s_1, N_s, K, wtol)
                totiter += 1
                if totiter >= max_totiter:
                    break

                for r in nnzclasses:
                    if np.sum(W_s[:, r]) == 0:
                        X_s[r] = 0.0
                    elif N_s[r] == 0:
                        X_s[r] = 0.0
                    else:
                        Cs = float(np.dot(V[:, r], W_s[:, r]))
                        if Cs > 1e-14:
                            X_s[r] = _OMICRON * N_s[r] / Cs + (1 - _OMICRON) * X_s_1[r]
                        else:
                            X_s[r] = X_s_1[r]
                    for k in range(Ms):
                        Q_s[k, r] = _OMICRON * X_s[r] * W_s[k, r] + (1 - _OMICRON) * Q_s_1[k, r]

            # class-aggregate correction into slice 0, see the module header
            if Nt > 1:
                for k in range(Ms):
                    gamma[s, k, 0] = np.sum(Q_s_1[k, :]) / (Nt - 1.0) - np.sum(Q_outer_1[k, :]) / Nt
            else:
                gamma[s, :, 0] = 0.0

            if totiter >= max_totiter:
                break
        if totiter >= max_totiter:
            break

        # sweep at the full population N
        inner_iter = 0
        Q_1 = Q + np.inf
        while (inner_iter < 2 or np.max(np.abs(Q - Q_1)) > tol) and inner_iter <= max_sweep:
            inner_iter += 1
            Q_1, X_1, U_1 = Q.copy(), X.copy(), U.copy()

            W, STeff = _forward(ST, srv, isdelay, mu_full, gamma, Q_1, N, K, wtol)
            totiter += 1
            if totiter >= max_totiter:
                break

            for r in nnzclasses:
                if np.sum(W[:, r]) == 0:
                    X[r] = 0.0
                elif N[r] == 0:
                    X[r] = 0.0
                    C[r] = 0.0
                else:
                    C[r] = float(np.dot(V[:, r], W[:, r]))
                    if C[r] > 1e-14:
                        X[r] = _OMICRON * N[r] / C[r] + (1 - _OMICRON) * X_1[r]
                    else:
                        X[r] = X_1[r]
                for k in range(Ms):
                    Q[k, r] = _OMICRON * X[r] * W[k, r] + (1 - _OMICRON) * Q_1[k, r]
                    T[k, r] = X[r]
                    U[k, r] = _OMICRON * STeff[k, r] * X[r] + (1 - _OMICRON) * U_1[k, r]

    # Utilization capping, as in solver_amvald: a queueing station whose class
    # utilizations sum above one has them renormalized in proportion to STeff.
    # Delay stations are exempt.
    for k in range(Ms):
        if isdelay[k]:
            continue
        U_sum = float(np.sum(U[k, :]))
        if U_sum > 1:
            denom = float(np.sum(STeff[k, :] * X))
            if denom > 0:
                for r in range(K):
                    if STeff[k, r] > 0:
                        U[k, r] = min(1.0, U_sum) * STeff[k, r] * X[r] / denom

    # WHICH UTILIZATION SolverMVA REPORTS DEPENDS ON THE MODEL. Its analyzer
    # forwards the iterated Uchain to sn_deaggregate_chain_results ONLY under
    # lld, cd or jd scaling; with none of those the deaggregation recomputes
    # T*S/c from the NOMINAL demand instead, and the two differ by the iteration
    # residual. Reproduced here on the same test, mu being the only one of the
    # three a demand matrix can carry.
    if mu is None:
        for k in range(Ms):
            for r in nnzclasses:
                U[k, r] = ST[k, r] * X[r] if np.isinf(srv[k]) else ST[k, r] * X[r] / srv[k]

    R = np.zeros((Ms, K))
    nz = T > 0
    R[nz] = Q[nz] / T[nz]

    # A class with no jobs keeps its 1/sum(ST) SEED in X unless it is cleared:
    # the sweeps only ever write the classes in nnzclasses, so the initial value
    # would otherwise be reported as that class's throughput. Q, U and C are
    # already zero there because they are written in the same loops.
    Xout = np.zeros(K)
    Xout[nnzclasses] = X[nnzclasses]

    keep = ~isdelay
    return (Q[keep, :], U[keep, :], R[keep, :],
            Xout.reshape(1, K), C.reshape(1, K), totiter)
