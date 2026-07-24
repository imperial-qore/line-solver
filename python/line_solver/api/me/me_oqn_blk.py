"""
Maximum Entropy algorithm for single-class open queueing networks with
finite buffers, loss and transfer blocking.

Extends me_oqn to open networks in which a station has a finite buffer. Two
per-station policies are supported. Under loss a job that finds the
destination full is discarded, so each station is a censored GE/GE/c/0;N
queue and the network is the ME decomposition of Kouvatsos (1994),
Section 4. Under transfer blocking a job that completes service at station i
and finds the destination j full is held in i's server, which cannot serve
anyone else until j has room (blocking after service, BAS).

Transfer blocking is not work conserving, so a product-form approximation
cannot be applied to the network as it stands. Following Tahilramani,
Manjunath and Bose (1999) the network is first made work conserving by
inserting a GE/GE/inf holding node on every routing pair with a
finite-buffer destination. The holding node absorbs the blocked job,
releasing station i's server; the delay it introduces is the residual life
of the minimum of the c_j service times in progress at j, inflated
geometrically because the released job may find j full again. Station i's
own service time is inflated by the same blocking probability so that the
jobs queued behind the blocked one still see the server as busy. The
expanded network is work conserving and is solved node by node with the
censored ME queue of me_gegecn, iterating over the blocking probabilities
and the first two moments of the flows until they converge.
"""

import warnings
from typing import Any, Dict, Optional, Tuple

import numpy as np

from .me_gegecn import me_gegecn, me_gegecn_pb
from .me_oqn import _ge_gec_mql

RULE_LOSS = 0
RULE_BAS = 1


def me_oqn_blk(M: int, lambda0, Ca0, mu, Cs, P, c, N, blockrule=None,
               options: Optional[Dict[str, Any]] = None
               ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray,
                          np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    """Solves a single-class open network with finite buffers.

    Args:
        M: number of stations
        lambda0: external arrival rates, shape (M,)
        Ca0: external interarrival scvs, shape (M,)
        mu: service rates, shape (M,)
        Cs: service scvs, shape (M,)
        P: routing probabilities, shape (M, M); a row sum below one sends the
            residual flow out of the network
        c: servers per station, shape (M,), inf marking an infinite server
        N: buffer capacity per station in jobs, in service included, shape
            (M,), inf marking an unbounded buffer
        blockrule: policy at each finite buffer, RULE_LOSS or RULE_BAS
        options: dict with keys 'tol' (1e-6), 'maxiter' (1000), 'verbose'
            (False) and 'damping' (0.5), the last applied to the blocking
            probabilities as in the relaxation scheme of the source

    Returns:
        (Q, W, T, U, Ca, Cd, PBa, lam, iter): mean queue lengths with the
        jobs held blocked included, response times, carried throughputs,
        utilizations with a blocked server NOT counted as busy, offered
        interarrival scvs, interdeparture scvs, blocking probabilities,
        offered arrival rates and the iteration count.
    """
    opts = dict(options or {})
    tol = opts.get('tol', 1e-6)
    maxiter = int(opts.get('maxiter', 1000))
    verbose = bool(opts.get('verbose', False))
    damping = opts.get('damping', 0.5)

    lambda0 = np.asarray(lambda0, dtype=float).reshape(M)
    Ca0 = np.asarray(Ca0, dtype=float).reshape(M)
    mu = np.asarray(mu, dtype=float).reshape(M)
    Cs = np.asarray(Cs, dtype=float).reshape(M)
    P = np.asarray(P, dtype=float).reshape(M, M)
    c = np.asarray(c, dtype=float).reshape(M)
    N = np.asarray(N, dtype=float).reshape(M)
    if blockrule is None:
        blockrule = np.zeros(M, dtype=int)
    blockrule = np.asarray(blockrule, dtype=int).reshape(M)

    finite_buf = np.isfinite(N) & np.isfinite(c)
    bas = finite_buf & (blockrule == RULE_BAS)

    for i in range(M):
        if finite_buf[i] and Cs[i] < 1 - 1e-12:
            raise ValueError('MEM with finite buffers requires a service scv of at least 1 at station %d: '
                             'the GE distribution is not defined for scv < 1.' % (i + 1))
        if lambda0[i] > 0 and Ca0[i] < 1 - 1e-12:
            raise ValueError('MEM with finite buffers requires an external interarrival scv of at least 1 '
                             'at station %d.' % (i + 1))

    # see _kb/03-api-layer.md for rationale
    Pf = P.copy()
    muf = mu.copy()
    Csf = Cs.copy()
    for i in range(M):
        pii = P[i, i]
        if pii > 0:
            muf[i] = mu[i] * (1 - pii)
            Csf[i] = pii + (1 - pii) * Cs[i]
            Pf[i, :] = P[i, :] / (1 - pii)
            Pf[i, i] = 0.0

    # see _kb/03-api-layer.md for rationale
    sigma_s = 2.0 / (Csf + 1.0)
    mu_res = c * muf * sigma_s

    Ca = np.ones(M)
    Cd = Csf.copy()
    PBs = np.zeros((M, M))
    PBe = np.zeros(M)
    PBh = np.zeros((M, M))
    PBa = np.zeros(M)
    Q = np.zeros(M)
    U = np.zeros(M)
    T = np.zeros(M)
    lam = np.zeros(M)
    ca_stream = np.ones((M, M))
    att_int = np.zeros((M, M))
    att_ext = np.zeros(M)
    delta = np.inf
    it = 0

    for it in range(1, maxiter + 1):
        ca_old = Ca.copy()
        pbs_old = PBs.copy()
        pbe_old = PBe.copy()

        # see _kb/03-api-layer.md for rationale
        pbf = np.zeros(M)
        for i in range(M):
            for j in range(M):
                if Pf[i, j] > 0 and bas[j]:
                    pbf[i] += Pf[i, j] * PBs[i, j]
        if np.any(pbf >= 1 - 1e-9):
            raise ValueError('MEM transfer-blocking fixed point saturates: a station is blocked with '
                             'probability one. The network has no stable operating point under BAS.')
        mu_eff = muf * (1 - pbf)
        cs_eff = pbf + Csf * (1 - pbf)

        # see _kb/03-api-layer.md for rationale
        A = np.zeros((M, M))
        b = lambda0 * (1 - PBe)
        for j in range(M):
            for i in range(M):
                if Pf[i, j] > 0:
                    if finite_buf[j] and not bas[j]:
                        A[i, j] = Pf[i, j] * (1 - PBs[i, j])
                    else:
                        A[i, j] = Pf[i, j]
        T = np.linalg.solve(np.eye(M) - A.T, b)
        T[T < 0] = 0.0

        # see _kb/03-api-layer.md for rationale
        att_ext = lambda0.copy()
        att_int[:] = 0.0
        for j in range(M):
            for i in range(M):
                if Pf[i, j] > 0:
                    if finite_buf[j] and bas[j]:
                        att_int[i, j] = T[i] * Pf[i, j] / max(1 - PBs[i, j], 1e-12)
                    else:
                        att_int[i, j] = T[i] * Pf[i, j]
        lam = att_ext + np.sum(att_int, axis=0)
        for j in range(M):
            if lam[j] > 0:
                PBa[j] = (att_ext[j] * PBe[j] + float(np.sum(att_int[:, j] * PBs[:, j]))) / lam[j]
            else:
                PBa[j] = 0.0

        # see _kb/03-api-layer.md for rationale
        for j in range(M):
            if lam[j] <= 0:
                continue
            sum_inv = 0.0
            if att_ext[j] > 0:
                sum_inv += (att_ext[j] / lam[j]) / (Ca0[j] + 1)
            for i in range(M):
                ca_stream[i, j] = 1.0
                if att_int[i, j] > 0:
                    ca_stream[i, j] = 1 - Pf[i, j] + Pf[i, j] * Cd[i]
                    sum_inv += (att_int[i, j] / lam[j]) / (ca_stream[i, j] + 1)
            if sum_inv > 0:
                Ca[j] = -1 + 1 / sum_inv

        # see _kb/03-api-layer.md for rationale
        pbe_new = np.zeros(M)
        pbs_new = np.zeros((M, M))
        pbh_new = np.zeros((M, M))
        for j in range(M):
            if np.isinf(c[j]):
                # Infinite server: no queueing and no blocking
                Q[j] = lam[j] / mu_eff[j] if mu_eff[j] > 0 else 0.0
                U[j] = Q[j]
                Cd[j] = Ca[j]
                continue
            if not finite_buf[j]:
                # Unbounded buffer: the infinite-capacity GE building blocks
                rho = lam[j] / (c[j] * mu_eff[j]) if mu_eff[j] > 0 else 0.0
                if rho >= 1:
                    Q[j] = np.inf
                    U[j] = 1.0
                    Cd[j] = cs_eff[j]
                elif c[j] == 1:
                    Q[j] = rho * (Ca[j] + 1) / 2 + rho ** 2 * (cs_eff[j] + Ca[j]) / (2 * (1 - rho))
                    U[j] = rho
                    Cd[j] = rho ** 2 * cs_eff[j] + (1 - rho) * Ca[j] + rho * (1 - rho)
                else:
                    Q[j] = _ge_gec_mql(lam[j], Ca[j], mu_eff[j], cs_eff[j], int(c[j]))
                    U[j] = rho
                    Cd[j] = rho ** 2 * cs_eff[j] + (1 - rho) * Ca[j] + rho * (1 - rho)
                continue
            cj = int(c[j])
            nj = int(N[j])
            pj, lj, uj, _, _ = me_gegecn(lam[j], Ca[j], mu_eff[j], cs_eff[j], cj, 0, nj)
            Q[j] = lj
            U[j] = uj
            # see _kb/03-api-layer.md for rationale
            Cd[j] = uj ** 2 * cs_eff[j] + (1 - uj) * Ca[j] + uj * (1 - uj)
            if att_ext[j] > 0:
                pbe_new[j] = me_gegecn_pb(pj, 0, nj, cj, cs_eff[j], Ca0[j])
            for i in range(M):
                if att_int[i, j] > 0:
                    pbs_new[i, j] = me_gegecn_pb(pj, 0, nj, cj, cs_eff[j], ca_stream[i, j])
                    if bas[j]:
                        # see _kb/03-api-layer.md for rationale
                        q = Pf[i, j] * PBs[i, j]
                        ca_h = 1 - q + q * Cd[i]
                        pbh_new[i, j] = me_gegecn_pb(pj, 0, nj, cj, cs_eff[j], ca_h)

        # Relaxation on the blocking probabilities
        PBe = (1 - damping) * PBe + damping * pbe_new
        PBs = (1 - damping) * PBs + damping * pbs_new
        PBh = (1 - damping) * PBh + damping * pbh_new

        delta = max(float(np.max(np.abs(Ca - ca_old))),
                    float(np.max(np.abs(PBs - pbs_old))),
                    float(np.max(np.abs(PBe - pbe_old))))
        if verbose:
            print('Iteration %d: max delta = %e' % (it, delta))
        if delta < tol:
            break

    if it == maxiter and delta >= tol:
        warnings.warn('me_oqn_blk did not converge within %d iterations (delta=%e)' % (maxiter, delta))

    # see _kb/03-api-layer.md for rationale
    l_hold = np.zeros((M, M))
    for i in range(M):
        for j in range(M):
            if bas[j] and Pf[i, j] > 0 and PBs[i, j] > 0:
                rate_h = T[i] * Pf[i, j] * PBs[i, j]
                mu_h = mu_res[j] * (1 - PBh[i, j])
                if mu_h > 0:
                    l_hold[i, j] = rate_h / mu_h
    Q = Q + np.sum(l_hold, axis=1)

    # see _kb/03-api-layer.md for rationale
    for i in range(M):
        if muf[i] > 0:
            U[i] = T[i] / muf[i] if np.isinf(c[i]) else T[i] / (c[i] * muf[i])
        else:
            U[i] = 0.0

    W = np.zeros(M)
    for i in range(M):
        if T[i] > 0:
            W[i] = Q[i] / T[i]
    return Q, W, T, U, Ca, Cd, PBa, lam, it
