"""
Birman-Kogan asymptotic evaluation of closed networks with many stations.

Native Python implementation of the saddle point expansion with bottleneck
detection (Propositions 1 and 3 with Algorithm 1), the van der Waerden uniform
expansion for a single chain, and the load concealment reduction of a multichain
network to single chain problems (Algorithm 2).

Key functions:
    pfqn_bk: saddle point normalizing constant with bottleneck detection
    pfqn_bkue: uniform (van der Waerden) expansion, single chain
    pfqn_bklc: load concealment reduction, Algorithm 2

References:
    Original MATLAB: matlab/src/api/pfqn/pfqn_bk.m, pfqn_bkue.m, pfqn_bklc.m
    A. Birman, Y. Kogan, "Asymptotic evaluation of closed queueing networks
    with many stations", Communications in Statistics. Stochastic Models
    8(3):543-563, 1992
"""

import numpy as np
from scipy.special import erfc, erfcx, gammaln
from typing import Optional, Tuple

from .mva import pfqn_mva_single_class

FINE_TOL = 1e-8


def _multiplicity(L: np.ndarray) -> np.ndarray:
    """Number of stations sharing each demand row, up to relative rounding."""
    M = L.shape[0]
    mult = np.ones(M, dtype=int)
    order = np.lexsort(tuple(L[:, j] for j in range(L.shape[1] - 1, -1, -1)))
    Ls = L[order, :]
    start = 0
    for i in range(1, M + 1):
        same = False
        if i < M:
            scale = max(1.0, np.max(np.abs(Ls[i, :])), np.max(np.abs(Ls[i - 1, :])))
            same = np.max(np.abs(Ls[i, :] - Ls[i - 1, :])) <= FINE_TOL * scale
        if not same:
            mult[order[start:i]] = i - start
            start = i
    return mult


def _psi(z, Lg, N, Z):
    f = float(np.dot(Z, z) - np.dot(N, np.log(z)))
    if Lg.shape[0] > 0:
        f -= float(np.sum(np.log(1.0 - Lg @ z)))
    return f


def _grad(z, Lg, N, Z):
    g = Z - N / z
    if Lg.shape[0] > 0:
        g = g + (1.0 / (1.0 - Lg @ z)) @ Lg
    return g


def _hessian(z, Lg, N):
    H = np.diag(N / z ** 2)
    if Lg.shape[0] > 0:
        d = 1.0 / (1.0 - Lg @ z)
        H = H + Lg.T @ ((d ** 2)[:, None] * Lg)
    return H


def _init(Lg, N, Z, mu):
    den = Z.copy()
    if Lg.shape[0] > 0:
        den = den + np.sum(Lg, axis=0)
    z = N / np.maximum(den, FINE_TOL)
    z = np.minimum(z, 0.99 * mu)
    if Lg.shape[0] > 0:
        for _ in range(200):
            if np.max(Lg @ z) < 0.9:
                break
            z = z * 0.7
    return z


def _minimize(Lg, N, Z, mu):
    """Minimize psi over {z>0, Lg z<1, z<=mu}: Algorithm 1 as an active set method."""
    R = len(N)
    on_bound = np.zeros(R, dtype=bool)
    z = _init(Lg, N, Z, mu)
    for _ in range(R + 1):
        z[on_bound] = mu[on_bound]
        free = np.where(~on_bound)[0]
        if free.size == 0:
            break
        for _ in range(500):
            g = _grad(z, Lg, N, Z)
            if np.linalg.norm(g[free]) <= 1e-12 * max(1.0, float(np.sum(N))):
                break
            H = _hessian(z, Lg, N)
            dz = np.zeros(R)
            dz[free] = -np.linalg.solve(H[np.ix_(free, free)], g[free])
            alpha = 1.0
            while True:
                zt = z.copy()
                zt[free] = z[free] + alpha * dz[free]
                feasible = np.all(zt[free] > 0) and np.all(zt[free] <= mu[free])
                if feasible and (Lg.shape[0] == 0 or np.max(Lg @ zt) < 1.0):
                    break
                alpha /= 2
                if alpha < 1e-14:
                    break
            if alpha < 1e-14:
                break
            z[free] = z[free] + alpha * dz[free]
        g = _grad(z, Lg, N, Z)
        newly = free[(z[free] >= mu[free] * (1 - 1e-9)) & (g[free] < 0)]
        if newly.size == 0:
            break
        on_bound[newly] = True
    z[on_bound] = mu[on_bound]
    return z, on_bound


def pfqn_bk(L: np.ndarray, N: np.ndarray, Z: Optional[np.ndarray] = None):
    """
    Birman-Kogan saddle point normalizing constant with bottleneck detection.

    Stations that serve a single chain and appear only once (the paper's
    dedicated single servers) stay outside the exponent as O(1) algebraic
    factors, so their poles may be crossed by the saddle point; Algorithm 1
    detects those chains and pins their coordinate on the pole. The remaining
    stations are the paper's large groups of identical stations.

    Args:
        L: Service demand matrix (M x R)
        N: Population vector (R,)
        Z: Think time vector (R,), optional

    Returns:
        Tuple of (G, lG, X, U, A, B):
            G, lG: normalizing constant and its logarithm
            X: chain throughputs (the saddle point coordinates)
            U: utilizations (M x R)
            A: chains whose dedicated station is not saturated (eq. 29)
            B: chains whose dedicated station is a bottleneck (eq. 30)
    """
    L = np.atleast_2d(np.asarray(L, dtype=float))
    N = np.asarray(N, dtype=float).flatten()
    M, R = L.shape
    if Z is None:
        Z = np.zeros(R)
    else:
        Z = np.asarray(Z, dtype=float)
        Z = Z.sum(axis=0).flatten() if Z.ndim > 1 else Z.flatten()
    if L.size == 0 or N.size == 0 or np.sum(N) == 0:
        return 1.0, 0.0, np.zeros(R), np.zeros((M, R)), np.arange(R), np.array([], dtype=int)
    if np.any(N == 0) and np.any(N > 0):
        keep = N > 0
        idx = np.where(keep)[0]
        _, lG, Xk, Uk, Ak, Bk = pfqn_bk(L[:, keep], N[keep], Z[keep])
        X = np.zeros(R)
        U = np.zeros((M, R))
        X[idx] = Xk
        U[:, idx] = Uk
        return float(np.exp(lG)), lG, X, U, idx[Ak], idx[Bk]

    Lq = L[np.sum(L, axis=1) > 0, :]
    Mq = Lq.shape[0]

    # Dedicated station of each chain: single chain, no identical twin, no think
    # time, and only when the model holds a group of identical stations, since
    # it is against M_j >> 1 replicas that a lone station is an O(1) factor.
    mu = np.full(R, np.inf)
    pole_row = -np.ones(R, dtype=int)
    is_pole = np.zeros(Mq, dtype=bool)
    if Mq > 0:
        mult = _multiplicity(Lq)
        if np.any(mult > 1):
            for ist in range(Mq):
                nz = np.where(Lq[ist, :] > 0)[0]
                if nz.size != 1 or mult[ist] > 1:
                    continue
                r = int(nz[0])
                if Z[r] > FINE_TOL:
                    continue
                if 1.0 / Lq[ist, r] < mu[r]:
                    mu[r] = 1.0 / Lq[ist, r]
                    pole_row[r] = ist
        for r in range(R):
            if pole_row[r] >= 0:
                is_pole[pole_row[r]] = True
    Lg = Lq[~is_pole, :] if Mq > 0 else np.zeros((0, R))

    z, on_bound = _minimize(Lg, N, Z, mu)
    A = np.where(~on_bound)[0]
    B = np.where(on_bound)[0]
    X = z.copy()
    U = L * X[None, :]
    for r in B:
        U[:, r] = np.minimum(U[:, r], 1.0)

    psi0 = _psi(z, Lg, N, Z)
    if A.size == 0:  # eq. (25): the residues carry everything
        lG = psi0
    else:
        H = _hessian(z, Lg, N)[np.ix_(A, A)]
        sign, logdetH = np.linalg.slogdet(H)
        if sign <= 0:
            logdetH = float(np.sum(np.log(np.abs(np.linalg.eigvals(H)))))
        lG = psi0 - 0.5 * A.size * np.log(2 * np.pi) - 0.5 * logdetH - float(np.sum(np.log(z[A])))
        for r in A:
            if np.isfinite(mu[r]):
                lG -= np.log(1.0 - z[r] / mu[r])
    return float(np.exp(lG)), float(lG), X, U, A, B


def _h1(z, D, N, Z):
    return Z * z - N * np.log(z) - float(np.sum(np.log(1.0 - D * z)))


def _h1d1(z, D, N, Z):
    return Z - N / z + float(np.sum(D / (1.0 - D * z)))


def _h1d2(z, D, N):
    return N / z ** 2 + float(np.sum(D ** 2 / (1.0 - D * z) ** 2))


def _h1d3(z, D, N):
    return -2 * N / z ** 3 + 2 * float(np.sum(D ** 3 / (1.0 - D * z) ** 3))


def _saddle1(D, N, Z):
    if D.size == 0:
        return N / Z
    hi = 1.0 / np.max(D)
    z = 0.5 * hi
    for _ in range(200):
        g = _h1d1(z, D, N, Z)
        if abs(g) <= 1e-14 * max(1.0, N):
            break
        dz = -g / _h1d2(z, D, N)
        alpha = 1.0
        while z + alpha * dz <= 0 or z + alpha * dz >= hi:
            alpha /= 2
            if alpha < 1e-14:
                break
        if alpha < 1e-14:
            break
        z = z + alpha * dz
    return z


def pfqn_bkue(L: np.ndarray, N: float, Z: float = 0.0) -> Tuple[float, float]:
    """
    Birman-Kogan uniform (van der Waerden) expansion for a single chain.

    The plain saddle point loses accuracy once the saddle approaches the
    dominant pole of the integrand, which is the regime where the station
    holding that pole saturates. The uniform expansion keeps the pole and the
    saddle in one formula through the complementary error function.

    Args:
        L: Service demand vector (M,), single class
        N: Population (scalar)
        Z: Think time (scalar)

    Returns:
        Tuple of (G, lG)
    """
    L = np.asarray(L, dtype=float).flatten()
    L = L[L > 0]
    N = float(np.sum(np.asarray(N, dtype=float)))
    Z = float(np.sum(np.asarray(Z, dtype=float)))
    if N == 0:
        return 1.0, 0.0
    if L.size == 0:
        lG = N * np.log(Z) - gammaln(N + 1)
        return float(np.exp(lG)), float(lG)
    ipole = int(np.argmax(L))
    dmax = L[ipole]
    tolL = FINE_TOL * max(1.0, dmax)
    has_group = L.size > 1 and bool(np.any(np.abs(np.diff(np.sort(L))) <= tolL))
    has_pole = int(np.count_nonzero(np.abs(L - dmax) <= tolL)) == 1 and has_group
    if has_pole:
        D = np.delete(L, ipole)
        zp = 1.0 / dmax
    else:
        D = L
        zp = np.inf
    z0 = _saddle1(D, N, Z)
    h2 = _h1d2(z0, D, N)
    h3 = _h1d3(z0, D, N)
    if not np.isfinite(zp):
        # No pole to keep out of the exponent: the expansion degenerates to the
        # plain saddle point, and the third derivative term goes with the pole
        # it corrects.
        lG = _h1(z0, D, N, Z) - np.log(z0) - 0.5 * np.log(2 * np.pi * h2)
        return float(np.exp(lG)), float(lG)
    t2 = (1.0 / z0 + h3 / (6 * h2)) / np.sqrt(2 * np.pi * h2)
    b2 = _h1(zp, D, N, Z) - _h1(z0, D, N, Z)  # the paper's M*b^2, always >= 0
    b2 = max(b2, 0.0)
    if zp >= z0:  # saddle before the pole
        lG = _h1(z0, D, N, Z) + np.log(0.5 * erfcx(np.sqrt(b2)) + t2)
    else:  # the pole has been crossed and its residue leads
        lG = _h1(zp, D, N, Z) + np.log(1.0 - 0.5 * erfc(np.sqrt(b2)) + t2 * np.exp(-b2))
    return float(np.exp(lG)), float(lG)


def _ue_chain(D, N, Z):
    """Throughput and queue lengths of a single chain from the uniform expansion."""
    D = np.asarray(D, dtype=float).flatten()
    Q = np.zeros(D.size)
    X = 0.0
    lg_prev = 0.0
    for n in range(1, int(N) + 1):
        _, lg_n = pfqn_bkue(D, n, Z)
        X = float(np.exp(lg_prev - lg_n))
        Q = D * X * (1.0 + Q)
        lg_prev = lg_n
    return X, Q


def pfqn_bklc(L: np.ndarray, N: np.ndarray, Z: Optional[np.ndarray] = None,
                method: str = 'mva', tol: float = 1e-10, maxiter: int = 1000):
    """
    Birman-Kogan load concealment algorithm (Algorithm 2).

    Chain l is solved on its own with every station slowed by the residual
    capacity the other chains leave it, A_i = 1 - sum_{k != l} L(i,k) X_k.
    Sweeping the chains in Gauss-Seidel order and iterating to a fixed point is
    the load concealment algorithm.

    Args:
        L: Service demand matrix (M x R)
        N: Population vector (R,)
        Z: Think time vector (R,), optional
        method: single chain solver, 'mva' (default) or 'ue'
        tol: convergence tolerance on the throughputs
        maxiter: maximum number of sweeps

    Returns:
        Tuple of (X, Q, U, it)
    """
    L = np.atleast_2d(np.asarray(L, dtype=float))
    N = np.asarray(N, dtype=float).flatten()
    M, R = L.shape
    if Z is None:
        Z = np.zeros(R)
    else:
        Z = np.asarray(Z, dtype=float)
        Z = Z.sum(axis=0).flatten() if Z.ndim > 1 else Z.flatten()
    X = np.zeros(R)
    Q = np.zeros((M, R))
    U = np.zeros((M, R))
    it = 0
    if L.size == 0 or np.sum(N) == 0:
        return X, Q, U, it
    # Step 1: the saddle point utilizations of Corollary 1 seed the iteration
    try:
        _, _, X, _, _, _ = pfqn_bk(L, N, Z)
        X = np.asarray(X, dtype=float).flatten().copy()
    except Exception:
        X = np.zeros(R)
    X[~np.isfinite(X)] = 0.0
    X[X < 0] = 0.0
    for r in range(R):
        if X[r] == 0 and N[r] > 0:
            X[r] = N[r] / (Z[r] + float(np.sum(L[:, r])))
        cap = float(np.max(L[:, r]))
        if cap > 0:
            X[r] = min(X[r], 1.0 / cap)

    for it in range(1, maxiter + 1):
        Xold = X.copy()
        for l in range(R):
            if N[l] == 0:
                X[l] = 0.0
                Q[:, l] = 0.0
                continue
            # Step 2a: residual capacity left to chain l at every station
            A = 1.0 - (L @ X - L[:, l] * X[l])
            A = np.maximum(A, FINE_TOL)
            D = L[:, l] / A
            if method == 'ue':
                Xl, Ql = _ue_chain(D, N[l], Z[l])
            else:
                res = pfqn_mva_single_class(int(N[l]), D, float(Z[l]))
                Xl = float(res['X'])
                Ql = np.asarray(res['Q'], dtype=float).flatten()
            X[l] = Xl
            Q[:, l] = Ql
        if np.max(np.abs(X - Xold)) <= tol * max(1.0, float(np.max(np.abs(X)))):
            break
    U = L * X[None, :]
    return X, Q, U, it
