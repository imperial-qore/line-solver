"""
Maximum Entropy Method for Closed Queueing Networks.

Implements the two-stage ME algorithm from Kouvatsos (1994) "Entropy
Maximisation and Queueing Network Models", Section 3.3, for closed
multiclass networks of G/G/1 and G/G/inf queues.
"""

import math
import warnings
from typing import Any, Dict, Optional, Tuple

import numpy as np


def me_cqn(
    M: int,
    R: int,
    N,
    mu,
    Cs,
    P,
    c: Optional[np.ndarray] = None,
    refstat: Optional[np.ndarray] = None,
    insens: Optional[np.ndarray] = None,
    options: Optional[Dict[str, Any]] = None
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray,
           np.ndarray, np.ndarray, np.ndarray, int]:
    """
    Maximum Entropy algorithm for Closed Queueing Networks.

    Stage 1 solves a pseudo-open network (no external arrivals) subject to
    job flow conservation and the population constraints sum_i L[i,r]=N[r],
    using the GE-type fixed point of the open algorithm (Section 3.2) on
    class-composed flows. Stage 2 builds the ME product-form solution (3.8)
    from the Lagrangian coefficients of Stage 1, computes the normalising
    constant Z(N) by a multiclass convolution, and iterates the flow (work
    rate) equations until the class throughputs implied by the closed ME
    solution agree with those used to parametrise the building blocks.

    Args:
        M: Number of queues (stations).
        R: Number of job classes.
        N: Class populations [R vector].
        mu: Service rates [M x R matrix].
        Cs: Service scv [M x R matrix].
        P: Routing probability matrix [M x M x R], P[j,i,r] = p_ji,r.
        c: Optional servers per queue [M vector]; np.inf marks an IS queue;
           finite values must be 1 (default: all ones).
        refstat: Optional reference station per class [R vector, 0-based];
           negative entries select the first station visited by the class.
        options: Optional dict with tol / maxiter / verbose fields.

    Returns:
        Tuple of (L, W, Ca, Cd, lambda_arr, rho, X, iters): closed mean
        queue lengths, response times, pseudo-open arrival/departure scvs,
        per-station class throughputs, closed utilizations, class
        throughputs at the reference stations, iteration count.

    Reference:
        Kouvatsos (1994), Section 3.3.
    """
    if options is None:
        options = {}
    tol = options.get('tol', 1e-6)
    maxiter = options.get('maxiter', 1000)
    verbose = options.get('verbose', False)

    N = np.asarray(N, dtype=float).ravel()
    mu = np.atleast_2d(np.asarray(mu, dtype=float))
    Cs = np.atleast_2d(np.asarray(Cs, dtype=float))
    P = np.asarray(P, dtype=float)
    if c is None:
        c = np.ones(M)
    else:
        c = np.asarray(c, dtype=float).ravel()
    if refstat is None:
        refstat = -np.ones(R, dtype=int)
    else:
        refstat = np.asarray(refstat, dtype=int).ravel()
    if insens is None:
        insens = np.zeros(M, dtype=bool)
    else:
        insens = np.asarray(insens, dtype=bool).ravel()

    # Feedback correction (as in the open algorithm)
    P_eff = P.copy()
    mu_eff = mu.copy()
    Cs_eff = Cs.copy()
    selfp = np.zeros((M, R))
    for i in range(M):
        for r in range(R):
            pii = P[i, i, r]
            if pii > 0:
                selfp[i, r] = pii
                mu_eff[i, r] = mu[i, r] * (1.0 - pii)
                Cs_eff[i, r] = pii + (1.0 - pii) * Cs[i, r]
                P_eff[i, :, r] = P[i, :, r] / (1.0 - pii)
                P_eff[i, i, r] = 0.0

    # Visit ratios from the original routing (visit-inclusive), normalised
    # at the reference station of each class
    V = np.zeros((M, R))
    for r in range(R):
        A = np.eye(M) - P[:, :, r].T
        ref = refstat[r]
        if ref < 0:
            ref = int(np.argmax(mu[:, r] > 0))
        A[ref, :] = 0.0
        A[ref, ref] = 1.0
        b = np.zeros(M)
        b[ref] = 1.0
        V[:, r] = np.linalg.solve(A, b)
        V[np.abs(V[:, r]) < 1e-14, r] = 0.0
        refstat[r] = ref

    # Stage 1: pseudo-open network, find X such that sum_i L[i,r] = N[r]
    X = np.zeros(R)
    for r in range(R):
        capr = np.inf
        for i in range(M):
            if not np.isinf(c[i]) and V[i, r] > 0 and mu[i, r] > 0:
                capr = min(capr, mu[i, r] / V[i, r])
        if np.isinf(capr):  # IS-only class
            capr = 1.0
        X[r] = 0.5 * capr / R

    Ca = np.ones((M, R))
    iters = 0
    Lpo = np.zeros((M, R))
    rho_po = np.zeros((M, R))
    Cd = np.ones((M, R))
    lam = np.zeros((M, R))
    # see _kb/03-api-layer.md for rationale
    maxit1 = min(maxiter, 100)
    for it1 in range(maxit1):
        iters += 1
        X = _capacity_cap(X, V, mu, c, M, R)
        lam = V * X[np.newaxis, :]
        Lpo, Ca, Cd, rho_po = _pseudoopen(M, R, lam, mu, mu_eff, Cs_eff,
                                          P_eff, selfp, c, insens, Ca, tol, maxiter)
        Ltot = np.sum(Lpo, axis=0)
        err1 = 0.0
        for r in range(R):
            if N[r] > 0 and Ltot[r] > 0:
                err1 = max(err1, abs(Ltot[r] - N[r]) / N[r])
        if verbose:
            print(f"Stage 1 iteration {it1 + 1}: max population error = {err1:e}")
        if err1 < tol:
            break
        Xold = X.copy()
        for r in range(R):
            if Ltot[r] > 0:
                fac = min(max((N[r] / Ltot[r]) ** 0.5, 0.25), 4.0)  # clamped step
                X[r] = 0.5 * X[r] + 0.5 * X[r] * fac                # damped update
        # Stall guard: the stability cap can bind before the population
        # target is met (bottleneck saturation); stop when X no longer moves
        Xcap = _capacity_cap(X.copy(), V, mu, c, M, R)
        if np.max(np.abs(Xcap - Xold) / np.maximum(Xold, 1e-12)) < tol:
            break

    # Stage 2: closed ME solution by convolution, iterated on the flow
    # (work rate) equations
    sz = (N + 1).astype(int)
    Dec = _lattice(sz, R)
    PIdx = Dec.shape[0]
    rad = np.ones(R, dtype=int)
    for r in range(1, R):
        rad[r] = rad[r - 1] * sz[r - 1]
    L = Lpo.copy()
    rho = rho_po.copy()
    err2 = np.inf
    for it2 in range(maxiter):
        iters += 1
        F = _coefficients(M, R, PIdx, Dec, sz, Lpo, rho_po, lam,
                          mu_eff, Cs_eff, Ca, c, selfp)
        L, U = _convolve(M, R, N.astype(int), PIdx, Dec, rad, F)
        # Utilization split by pseudo-open per-class load; implied
        # throughputs from the work rate theorem with visit weights
        rho = np.zeros((M, R))
        Xhat = np.zeros(R)
        for r in range(R):
            num = 0.0
            den = 0.0
            for i in range(M):
                if lam[i, r] > 0:
                    if np.isinf(c[i]):
                        rho[i, r] = L[i, r]
                        # IS work rate: lambda_eff = L*mu_eff, revisits add 1/(1-p)
                        num += L[i, r] * mu_eff[i, r] / (1.0 - selfp[i, r])
                    else:
                        rho_i = np.sum(rho_po[i, :])
                        if rho_i > 0:
                            rho[i, r] = U[i] * rho_po[i, r] / rho_i
                        num += rho[i, r] * mu[i, r]
                    den += V[i, r]
            if den > 0:
                Xhat[r] = num / den
        err2 = 0.0
        for r in range(R):
            if X[r] > 0:
                err2 = max(err2, abs(Xhat[r] - X[r]) / X[r])
        if verbose:
            print(f"Stage 2 iteration {it2 + 1}: max flow error = {err2:e}")
        if err2 < tol:
            break
        Xold = X.copy()
        X = 0.5 * X + 0.5 * Xhat
        X = _capacity_cap(X, V, mu, c, M, R)
        # Stall guard: when the stability cap binds, X stops moving even
        # though the residual flow error stays above tolerance
        if np.max(np.abs(X - Xold) / np.maximum(Xold, 1e-12)) < tol:
            break
        lam = V * X[np.newaxis, :]
        Lpo, Ca, Cd, rho_po = _pseudoopen(M, R, lam, mu, mu_eff, Cs_eff,
                                          P_eff, selfp, c, insens, Ca, tol, maxiter)

    if it2 == maxiter - 1 and err2 >= tol:
        warnings.warn(
            f"Did not converge within {maxiter} iterations (flow error={err2:e})",
            RuntimeWarning)

    # Response times by Little's law on the visit-inclusive throughputs
    lam = V * X[np.newaxis, :]
    W = np.zeros((M, R))
    for i in range(M):
        for r in range(R):
            if lam[i, r] > 0:
                W[i, r] = L[i, r] / lam[i, r]

    return L, W, Ca, Cd, lam, rho, X, iters


def _capacity_cap(X, V, mu, c, M, R):
    """Scales the class throughputs uniformly so that every single-server
    queue in the pseudo-open network remains stable."""
    maxrho = 0.0
    for i in range(M):
        if not np.isinf(c[i]):
            rho_i = 0.0
            for r in range(R):
                if V[i, r] > 0 and mu[i, r] > 0:
                    rho_i += X[r] * V[i, r] / mu[i, r]
            maxrho = max(maxrho, rho_i)
    if maxrho >= 0.999:
        X = X * (0.999 / maxrho)
    return X


def _pseudoopen(M, R, lam, mu, mu_eff, Cs_eff, P_eff, selfp, c, insens, Ca, tol, maxiter):
    """GE-type fixed point of the open algorithm (Section 3.2) on the
    pseudo-open network: no external arrivals, flows given by lam. The flow
    scvs are computed on the class-composed (aggregate) streams and
    disaggregated per class by thinning, following the class composition
    and disaggregation principle of the closed ME algorithm."""
    lambda_eff = lam * (1.0 - selfp)
    rho = np.zeros((M, R))
    for i in range(M):
        for r in range(R):
            if mu[i, r] > 0:
                if np.isinf(c[i]):
                    rho[i, r] = lambda_eff[i, r] / mu_eff[i, r]
                else:
                    rho[i, r] = lam[i, r] / mu[i, r]
    # Class composition per station: aggregate flow, service process
    # moments and flow-weighted aggregate routing
    lam_a = np.sum(lambda_eff, axis=1)
    mu_a = np.zeros(M)
    Cs_a = np.ones(M)
    for i in range(M):
        if lam_a[i] > 0:
            ES = 0.0
            ES2 = 0.0
            for u in range(R):
                if lambda_eff[i, u] > 0 and mu_eff[i, u] > 0:
                    wu = lambda_eff[i, u] / lam_a[i]
                    ES += wu / mu_eff[i, u]
                    ES2 += wu * (Cs_eff[i, u] + 1.0) / mu_eff[i, u] ** 2
            if ES > 0:
                mu_a[i] = 1.0 / ES
                Cs_a[i] = ES2 / ES ** 2 - 1.0
    Pa = np.zeros((M, M))
    for j in range(M):
        if lam_a[j] > 0:
            for i in range(M):
                num = 0.0
                for r in range(R):
                    if lambda_eff[j, r] > 0:
                        num += lambda_eff[j, r] * P_eff[j, i, r]
                Pa[j, i] = num / lam_a[j]
    # Fixed point on the aggregate arrival scvs
    Ca_a = np.ones(M)
    for i in range(M):
        if lam_a[i] > 0:
            active = np.where(lambda_eff[i, :] > 0)[0]
            if active.size:
                wr = active[0]
                Ca_a[i] = 1.0 + (Ca[i, wr] - 1.0) * lam_a[i] / max(lambda_eff[i, wr], np.finfo(float).tiny)
    Cd_a = np.ones(M)
    L_a = np.zeros(M)
    for _ in range(maxiter):
        Ca_old = Ca_a.copy()
        for i in range(M):
            if lam_a[i] <= 0:
                continue
            rho_i = np.sum(rho[i, :])
            if np.isinf(c[i]):
                # GE/GE/inf: L = lambda/mu, departures inherit the arrival scv
                L_a[i] = lam_a[i] / mu_a[i]
                Cd_a[i] = Ca_a[i]
            elif rho_i < 1:
                if insens[i]:
                    # Insensitive disciplines (PS, LCFS-PR): product-form mql
                    L_a[i] = rho_i / (1.0 - rho_i)
                else:
                    # Single-class GE/GE/1 mql, eq. (3.6)
                    L_a[i] = (rho_i * (Ca_a[i] + 1.0) / 2.0
                              + rho_i ** 2 * (Ca_a[i] + Cs_a[i]) / (2.0 * (1.0 - rho_i)))
                Cd_a[i] = 2.0 * L_a[i] * (1.0 - rho_i) + Ca_a[i] * (1.0 - 2.0 * rho_i)
        # GE-type merging, eq. (3.7) with lambda_o = 0, on aggregate flows
        for i in range(M):
            if lam_a[i] > 0:
                sum_inv = 0.0
                for j in range(M):
                    if Pa[j, i] > 0 and lam_a[j] > 0:
                        Cdji = 1.0 + Pa[j, i] * (Cd_a[j] - 1.0)
                        sum_inv += (lam_a[j] * Pa[j, i] / lam_a[i]) / (Cdji + 1.0)
                if sum_inv > 0:
                    Ca_a[i] = -1.0 + 1.0 / sum_inv
        if np.max(np.abs(Ca_a - Ca_old)) < tol:
            break
    # Disaggregation: per-class arrival scvs by thinning of the composed
    # stream, then per-class mean queue lengths (Section 3.1.1)
    L = np.zeros((M, R))
    Cd = np.ones((M, R))
    for i in range(M):
        rho_i = np.sum(rho[i, :])
        for r in range(R):
            if lambda_eff[i, r] > 0:
                pr = lambda_eff[i, r] / lam_a[i]
                Ca[i, r] = 1.0 + pr * (Ca_a[i] - 1.0)
                Cd[i, r] = 1.0 + pr * (Cd_a[i] - 1.0)
        if np.isinf(c[i]):
            for r in range(R):
                if lambda_eff[i, r] > 0 and mu_eff[i, r] > 0:
                    L[i, r] = lambda_eff[i, r] / mu_eff[i, r]
        elif rho_i < 1:
            if insens[i]:
                # Insensitive disciplines (PS, LCFS-PR): product-form mql
                for r in range(R):
                    if lambda_eff[i, r] > 0 and mu_eff[i, r] > 0:
                        L[i, r] = rho[i, r] / (1.0 - rho_i)
            else:
                resid = 0.0
                for u in range(R):
                    if lambda_eff[i, u] > 0 and mu_eff[i, u] > 0:
                        resid += lambda_eff[i, u] * (Cs_eff[i, u] + Ca[i, u]) / mu_eff[i, u] ** 2
                for r in range(R):
                    if lambda_eff[i, r] > 0 and mu_eff[i, r] > 0:
                        L[i, r] = (rho[i, r] * (Ca[i, r] + 1.0) / 2.0
                                   + lambda_eff[i, r] * resid / (2.0 * (1.0 - rho_i)))
    return L, Ca, Cd, rho


def _lattice(sz, R):
    """Enumerates the population lattice {0..N[0]} x ... x {0..N[R-1]} in
    mixed radix order; Dec[p,:] is the population vector of index p."""
    PIdx = int(np.prod(sz))
    Dec = np.zeros((PIdx, R), dtype=int)
    for p in range(PIdx):
        q = p
        for r in range(R):
            Dec[p, r] = q % sz[r]
            q //= sz[r]
    return Dec


def _coefficients(M, R, PIdx, Dec, sz, Lpo, rho_po, lam, mu_eff, Cs_eff, Ca, c, selfp):
    """Auxiliary functions f_i(n) of the ME solution (3.8): the right-hand
    sides of (3.2) and (3.4) with the (1-rho) factor removed, evaluated from
    the Stage 1 Lagrangian coefficients. Each f_i is rescaled by its maximum
    for numerical stability (per-station constants cancel in the marginals)."""
    F = np.zeros((PIdx, M))
    lambda_eff = lam * (1.0 - selfp)
    for i in range(M):
        if np.isinf(c[i]):
            # GE/GE/inf: f(n) = prod_r prod_{k=1}^{n_r} g_r(k)
            logg = []
            for r in range(R):
                lg = np.zeros(sz[r] - 1) if sz[r] > 1 else np.zeros(0)
                for j in range(1, sz[r]):
                    if lambda_eff[i, r] > 0 and mu_eff[i, r] > 0:
                        gj = ((lambda_eff[i, r] * (1.0 + Cs_eff[i, r])
                               + (j - 1) * mu_eff[i, r] * (Ca[i, r] - 1.0))
                              / (j * mu_eff[i, r] * (Ca[i, r] + Cs_eff[i, r])))
                        lg[j - 1] = math.log(max(gj, 0.0)) if gj > 0 else -np.inf
                    else:
                        lg[j - 1] = -np.inf
                logg.append(lg)
            for p in range(PIdx):
                n = Dec[p, :]
                val = 0.0
                for r in range(R):
                    for j in range(1, n[r] + 1):
                        val += logg[r][j - 1]
                F[p, i] = math.exp(val) if val > -np.inf else 0.0
            F[0, i] = 1.0
        else:
            # see _kb/03-api-layer.md for rationale
            rho_i = np.sum(rho_po[i, :])
            Li = np.sum(Lpo[i, :])
            x = np.zeros(R)
            gx = np.zeros(R)
            if Li > 0 and rho_i < 1:
                for r in range(R):
                    if lam[i, r] > 0:
                        x[r] = max(Lpo[i, r] - rho_po[i, r], 0.0) / Li
                        gx[r] = rho_po[i, r] * rho_i / ((1.0 - rho_i) * Li)
            for p in range(PIdx):
                n = Dec[p, :]
                ntot = int(np.sum(n))
                if ntot == 0:
                    F[p, i] = 1.0
                    continue
                if any(n[r] > 0 and lam[i, r] <= 0 for r in range(R)):
                    F[p, i] = 0.0  # class not visiting this station
                    continue
                logmult = math.lgamma(ntot) - sum(math.lgamma(n[r] + 1) for r in range(R))
                tot = 0.0
                for r in range(R):
                    if n[r] > 0 and gx[r] > 0:
                        lterm = math.log(n[r]) + math.log(gx[r])
                        ok = True
                        for s in range(R):
                            es = n[s] - 1 if s == r else n[s]
                            if es > 0:
                                if x[s] > 0:
                                    lterm += es * math.log(x[s])
                                else:
                                    ok = False
                                    break
                        if ok:
                            tot += math.exp(logmult + lterm)
                F[p, i] = tot
        fmax = np.max(F[:, i])
        if fmax > 0:
            F[:, i] = F[:, i] / fmax
        F[0, i] = max(F[0, i], np.finfo(float).tiny)
    return F


def _convpair(G, f, PIdx, Dec, rad, Nvec):
    """Convolution of a partial normalising constant with one station term
    over the population lattice."""
    G2 = np.zeros(PIdx)
    for p in range(PIdx):
        if f[p] == 0:
            continue
        n = Dec[p, :]
        for q in range(PIdx):
            if G[q] == 0:
                continue
            t = n + Dec[q, :]
            if np.all(t <= Nvec):
                G2[int(np.dot(t, rad))] += f[p] * G[q]
    return G2


def _convolve(M, R, N, PIdx, Dec, rad, F):
    """Computes the normalising constant by convolving the f_i over the
    population lattice, and the per-station marginals by prefix/suffix
    convolutions; returns closed mean queue lengths and busy probabilities."""
    G0 = np.zeros(PIdx)
    G0[0] = 1.0
    Gpre = [G0]
    for k in range(M):
        Gpre.append(_convpair(Gpre[k], F[:, k], PIdx, Dec, rad, N))
    Gsuf = [None] * (M + 1)
    Gsuf[M] = G0
    for k in range(M - 1, -1, -1):
        Gsuf[k] = _convpair(Gsuf[k + 1], F[:, k], PIdx, Dec, rad, N)
    Z = Gpre[M][PIdx - 1]
    L = np.zeros((M, R))
    U = np.zeros(M)
    for i in range(M):
        Grest = _convpair(Gpre[i], Gsuf[i + 1], PIdx, Dec, rad, N)
        for p in range(PIdx):
            if F[p, i] > 0:
                n = Dec[p, :]
                q = int(np.dot(N - n, rad))
                pin = F[p, i] * Grest[q] / Z
                if p > 0:
                    U[i] += pin
                for r in range(R):
                    if n[r] > 0:
                        L[i, r] += n[r] * pin
    return L, U


__all__ = ['me_cqn']
