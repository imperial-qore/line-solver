"""
Marie's iterative aggregation-decomposition (Marie 1979/1980) for closed
networks with FCFS non-exponential (Coxian) service. Native-Python port of the
MATLAB pfqn_marie single-class path: each station is analyzed in isolation as a
lambda(n)/Cox/1 queue whose conditional throughputs mu_i(n) drive the exact
load-dependent product-form aggregate (pfqn_mvald), iterated to a fixed point.
Reduces to exact product form for exponential service.

Multiclass (R>1) is handled by marie_multi: a QD-AMVA aggregate with class-
dependent scaling beta_{i,r}(nvec) = (Coxian isolation conditional throughput)/
(exponential isolation conditional throughput), so beta==1 recovers standard
FCFS AMVA and the cd-scaling carries only the non-exponential correction. beta
is supplied by a multiclass lambda_r/Cox/1 isolation CTMC over the joint
population lattice, iterated to a fixed point on X (Baynat-Dallery isolation).
"""

import numpy as np
from .mvald import pfqn_mvald

_COARSE = 1e-3


def _cox_fit(mean, scv):
    """Coxian phase representation matching Coxian.fitMeanAndSCV: returns
    (rates, compl) arrays. scv==1 exponential; scv<0.5 Erlang; 0.5<scv<1 and
    scv>1 two-phase Coxian."""
    if abs(scv - 1.0) <= _COARSE:
        return np.array([1.0 / mean]), np.array([1.0])
    if scv > 0.5 + _COARSE and scv < 1.0 - _COARSE:
        mu = np.zeros(2)
        mu[0] = 2.0 / mean / (1.0 + np.sqrt(1.0 + 2.0 * (scv - 1.0)))
        mu[1] = 2.0 / mean / (1.0 - np.sqrt(1.0 + 2.0 * (scv - 1.0)))
        phi = np.array([0.0, 1.0])
        return mu, phi
    if scv <= 0.5 + _COARSE:
        n = int(np.ceil(1.0 / scv))
        lam = n / mean
        mu = lam * np.ones(n)
        phi = np.zeros(n)
        phi[-1] = 1.0
        return mu, phi
    # scv > 1
    mu = np.zeros(2)
    mu[0] = 2.0 / mean
    mu[1] = mu[0] / (2.0 * scv)
    phi = np.array([1.0 - mu[1] / mu[0], 1.0])
    return mu, phi


def _isol_condtput(lam, ph_rate, ph_compl, N, m):
    """Stationary analysis of a lambda(n)/Cox/1(-m) isolation queue; returns
    mu(n) = conditional departure rate given n present, n=1..N (length-N array,
    index n-1)."""
    P = ph_rate.size
    S = 1 + N * P

    def idx(n, k):  # n>=1, k in 1..P
        return 1 + (n - 1) * P + (k - 1)  # 0-based; state 0 = empty

    Gq = np.zeros((S, S))
    Gq[0, idx(1, 1)] += lam[0]
    for n in range(1, N + 1):
        sc = min(n, m)
        for k in range(1, P + 1):
            r = idx(n, k)
            if n < N:
                Gq[r, idx(n + 1, k)] += lam[n]  # lam[n] = arrival rate at n present
            compl = ph_rate[k - 1] * ph_compl[k - 1] * sc
            adv = ph_rate[k - 1] * (1.0 - ph_compl[k - 1]) * sc
            if adv > 0 and k < P:
                Gq[r, idx(n, k + 1)] += adv
            if compl > 0:
                if n > 1:
                    Gq[r, idx(n - 1, 1)] += compl
                else:
                    Gq[r, 0] += compl
    Gq = Gq - np.diag(Gq.sum(axis=1))
    # solve p*Gq = 0, sum(p)=1  ->  [Gq'; 1...1] p = [0;1]
    A = np.vstack([Gq.T, np.ones(S)])
    b = np.concatenate([np.zeros(S), [1.0]])
    p, *_ = np.linalg.lstsq(A, b, rcond=None)
    muvec = np.zeros(N)
    for n in range(1, N + 1):
        Pn = 0.0
        dep = 0.0
        for k in range(1, P + 1):
            pk = p[idx(n, k)]
            Pn += pk
            dep += pk * ph_rate[k - 1] * ph_compl[k - 1] * min(n, m)
        if Pn > 0:
            muvec[n - 1] = dep / Pn
        else:
            muvec[n - 1] = min(n, m) / np.sum(1.0 / ph_rate)
    return muvec


def pfqn_marie(L, N, Z=None, scv=None, tol=1e-8, maxiter=1000, nservers=None):
    """Marie's method (single-class). L is M x 1 demands, N scalar population,
    Z scalar think time, scv per-station SCV (M,). Returns
    (X, Q, U, C, it, mu): X scalar chain throughput, Q/U/C per-station (M,),
    it iterations, mu the converged LD multiplier matrix (M x N).

    Multiclass (L with >1 column) raises NotImplementedError."""
    L = np.asarray(L, dtype=float)
    if L.ndim == 2 and L.shape[1] > 1:
        return marie_multi(L, N, Z, scv, tol=tol, maxiter=maxiter)
    L = L.ravel()
    M = L.size
    if Z is None:
        Z = 0.0
    Z = float(np.sum(Z))
    if scv is None:
        scv = np.ones(M)
    scv = np.asarray(scv, dtype=float).ravel()
    if nservers is None:
        nservers = np.ones(M)
    nservers = np.asarray(nservers, dtype=float).ravel()
    if nservers.size == 1:
        nservers = nservers * np.ones(M)

    ph_rate = [None] * M
    ph_compl = [None] * M
    for i in range(M):
        r, p = _cox_fit(L[i], scv[i])
        ph_rate[i] = r
        ph_compl[i] = p

    mu = np.zeros((M, N))
    for i in range(M):
        mu[i, :] = np.minimum(np.arange(1, N + 1), nservers[i])

    X = 0.0
    Q = np.zeros(M)
    U = np.zeros(M)
    C = np.zeros(M)
    it = 0
    while it < maxiter:
        it += 1
        XN, QN, UN, CN, _, _, pi = pfqn_mvald(L.reshape(-1, 1), np.array([N]),
                                              np.array([Z]), mu)
        Pg = np.atleast_2d(pi)
        mu_new = mu.copy()
        for i in range(M):
            Pi = Pg[i, 0:N + 1]
            lam = np.zeros(N)
            for n in range(0, N):
                pn = Pi[n]
                if pn > 0:
                    lam[n] = (mu[i, n] / L[i]) * Pi[n + 1] / pn
                else:
                    lam[n] = 0.0
            muabs = _isol_condtput(lam, ph_rate[i], ph_compl[i], N,
                                   nservers[i])
            mu_new[i, :] = muabs * L[i]
        delta = np.max(np.abs(mu_new - mu))
        mu = mu_new
        X = float(np.ravel(XN)[0])
        Q = np.ravel(QN).astype(float)
        U = np.ravel(UN).astype(float)
        C = np.ravel(CN).astype(float)
        if delta < tol:
            break
    return X, Q, U, C, it, mu


def _ndlininterp(A, x, Nvec):
    """Multilinear interpolation of ND array A (shape Nvec+1) at real point x,
    clamped to the box [0, Nvec]."""
    R = len(Nvec)
    x = np.minimum(np.maximum(np.asarray(x, dtype=float).ravel(), 0.0),
                   np.asarray(Nvec, dtype=float))
    lo = np.floor(x).astype(int)
    hi = np.minimum(lo + 1, np.asarray(Nvec, dtype=int))
    fr = x - lo
    v = 0.0
    for mask in range(2 ** R):
        w = 1.0
        sub = np.zeros(R, dtype=int)
        for dbit in range(R):
            if (mask >> dbit) & 1:
                sub[dbit] = hi[dbit]
                w *= fr[dbit]
            else:
                sub[dbit] = lo[dbit]
                w *= (1.0 - fr[dbit])
        if w == 0.0:
            continue
        v += w * A[tuple(sub)]
    return v


def _cdscale_eval(nv, muCox, muExp, Nvec):
    """beta_{r}(nv) = muCox_r(nv)/muExp_r(nv), clamped to [1e-3,1e3]; beta==1
    when the isolation throughputs coincide (exponential service)."""
    R = len(Nvec)
    be = np.ones(R)
    for r in range(R):
        num = _ndlininterp(muCox[r], nv, Nvec)
        den = _ndlininterp(muExp[r], nv, Nvec)
        if den > 0 and num > 0 and np.isfinite(num) and np.isfinite(den):
            be[r] = num / den
    return np.minimum(np.maximum(be, 1e-3), 1e3)


def _isol_mc(lam, ph_rate_row, ph_compl_row, Nvec):
    """Stationary analysis of a multiclass lambda_r/Cox/1 FCFS queue in
    isolation over the joint per-class population box [0..Nvec]. The head-of-
    line job is tracked as (class, phase); on a departure the next head class is
    drawn in random order (prob n_c/sum(n)). Returns mumat[r], an ND array over
    the box giving the conditional class-r throughput mu_r(nvec)."""
    R = len(lam)
    Pc = np.array([len(ph_rate_row[r]) for r in range(R)], dtype=int)
    boxsz = tuple(int(n) + 1 for n in Nvec)
    npops = int(np.prod(boxsz))

    # Enumerate states: id 0 reserved for the empty station.
    ids = [(-1, 0, 0)]                 # (popLinear, head, phase); head 0 => empty
    key2id = {}
    for p in range(npops):
        nvec = np.array(np.unravel_index(p, boxsz), dtype=int)
        if nvec.sum() == 0:
            continue
        for c in range(R):
            if nvec[c] > 0:
                for k in range(Pc[c]):
                    key2id[(p, c, k)] = len(ids)
                    ids.append((p, c, k))
    S = len(ids)

    def poplin(nvec):
        return int(np.ravel_multi_index(tuple(int(x) for x in nvec), boxsz))

    I = []
    J = []
    V = []

    def addrate(a, b, rate):
        I.append(a); J.append(b); V.append(rate)

    for s in range(S):
        p, c, k = ids[s]
        if c == 0:                    # empty station (head sentinel 0)
            nvec = np.zeros(R, dtype=int)
        else:
            nvec = np.array(np.unravel_index(p, boxsz), dtype=int)
        # arrivals
        for r in range(R):
            if nvec[r] < Nvec[r] and lam[r] > 0:
                nnew = nvec.copy(); nnew[r] += 1
                if c == 0:
                    addrate(s, key2id[(poplin(nnew), r, 0)], lam[r])
                else:
                    addrate(s, key2id[(poplin(nnew), c, k)], lam[r])
        if c == 0:
            continue
        rate = ph_rate_row[c][k]
        compl = rate * ph_compl_row[c][k]
        adv = rate * (1.0 - ph_compl_row[c][k])
        if adv > 0 and k < Pc[c] - 1:
            addrate(s, key2id[(p, c, k + 1)], adv)
        if compl > 0:
            nnew = nvec.copy(); nnew[c] -= 1
            if nnew.sum() == 0:
                addrate(s, 0, compl)
            else:
                tot = int(nnew.sum())
                for cp in range(R):
                    if nnew[cp] > 0:
                        addrate(s, key2id[(poplin(nnew), cp, 0)],
                                compl * nnew[cp] / tot)

    Gq = np.zeros((S, S))
    for a, b, rate in zip(I, J, V):
        Gq[a, b] += rate
    Gq = Gq - np.diag(Gq.sum(axis=1))

    A = np.vstack([Gq.T, np.ones(S)])
    b = np.concatenate([np.zeros(S), [1.0]])
    pvec, *_ = np.linalg.lstsq(A, b, rcond=None)

    Ppop = np.zeros(boxsz)
    dep = [np.zeros(boxsz) for _ in range(R)]
    for s in range(S):
        p, c, k = ids[s]
        if c == 0:
            continue
        sub = np.unravel_index(p, boxsz)
        Ppop[sub] += pvec[s]
        dep[c][sub] += pvec[s] * ph_rate_row[c][k] * ph_compl_row[c][k]
    mumat = []
    for r in range(R):
        mr = np.zeros(boxsz)
        pos = Ppop > 0
        mr[pos] = dep[r][pos] / Ppop[pos]
        mumat.append(mr)
    return mumat


def _amva_qd(L, N, Z, cds, tol=1e-9, maxiter=5000):
    """Multiclass Schweitzer AMVA with class-dependent service-rate scaling. The
    effective class-r demand at station i is L[i,r]/beta_{i,r}(nvec_arrival),
    beta supplied by cds[i] evaluated at the arrival-instant per-class
    population."""
    M, R = L.shape
    N = np.asarray(N, dtype=float).ravel()
    Q = np.tile(N, (M, 1)) / max(M, 1)
    W = np.zeros((M, R))
    X = np.zeros(R)
    Qprev = Q + 1.0
    it = 0
    while np.max(np.abs(Q - Qprev)) > tol and it < maxiter:
        it += 1
        Qprev = Q.copy()
        for r in range(R):
            for i in range(M):
                nv = Q[i, :].copy()
                if N[r] > 0:
                    nv[r] = Q[i, r] * (N[r] - 1.0) / N[r]
                be = cds[i](nv)
                Leff = L[i, :] / be
                W[i, r] = Leff[r] + np.sum(Leff * nv)
            denom = Z[r] + np.sum(W[:, r])
            X[r] = N[r] / denom if denom > 0 else 0.0
            for i in range(M):
                Q[i, r] = X[r] * W[i, r]
    C = W.copy()
    U = np.zeros((M, R))
    for r in range(R):
        for i in range(M):
            U[i, r] = X[r] * L[i, r]
    return X, Q, U, C


def marie_multi(L, N, Z=None, scv=None, tol=1e-8, maxiter=1000):
    """Marie's method for multiclass FCFS Coxian closed networks (QD-AMVA
    aggregate with class-dependent scaling). Returns (X, Q, U, C, it, mu) where
    X is (R,), Q/U/C are (M,R), and mu is the list of converged cd-scaling
    callables. Reduces to exact MVA for class-independent exponential service."""
    from .mva import pfqn_mva
    L = np.asarray(L, dtype=float)
    M, R = L.shape
    N = np.asarray(N, dtype=float).ravel()
    if Z is None:
        Z = np.zeros(R)
    Z = np.asarray(Z, dtype=float).ravel()
    if Z.size == 1 and R > 1:
        Z = np.full(R, float(Z))
    if scv is None:
        scv = np.ones((M, R))
    scv = np.asarray(scv, dtype=float)
    if scv.shape != (M, R):
        scv = np.broadcast_to(scv, (M, R)).copy()

    # Exact product-form dispatch: class-independent exponential FCFS is genuine
    # BCMP -> exact MVA.
    isPF = bool(np.all(scv == 1.0))
    if isPF:
        for i in range(M):
            if L[i, :].max() - L[i, :].min() > 1e-12:
                isPF = False
                break
    if isPF:
        XN, CN, QN, UN, RN, TN, AN = pfqn_mva(L, N, Z)
        return (np.ravel(XN).astype(float), np.asarray(QN, float),
                np.asarray(UN, float), np.asarray(CN, float), 0, [])

    Nint = np.round(N).astype(int)
    # Per-station per-class Coxian phases, plus an exponential reference.
    phR = [[None] * R for _ in range(M)]
    phP = [[None] * R for _ in range(M)]
    eR = [[None] * R for _ in range(M)]
    eP = [[None] * R for _ in range(M)]
    for i in range(M):
        for r in range(R):
            mir, pir = _cox_fit(L[i, r], scv[i, r])
            phR[i][r] = np.asarray(mir, dtype=float)
            phP[i][r] = np.asarray(pir, dtype=float)
            eR[i][r] = np.array([1.0 / L[i, r]])
            eP[i][r] = np.array([1.0])

    cds = [(lambda nv: np.ones(R)) for _ in range(M)]

    Xprev = np.full(R, np.inf)
    X = np.zeros(R); Q = np.zeros((M, R)); U = np.zeros((M, R))
    C = np.zeros((M, R))
    it = 0
    while it < maxiter:
        it += 1
        X, Q, U, C = _amva_qd(L, Nint, Z, cds)
        for i in range(M):
            muCox = _isol_mc(X, phR[i], phP[i], Nint)
            muExp = _isol_mc(X, eR[i], eP[i], Nint)

            def make(muC, muE):
                return lambda nv: _cdscale_eval(nv, muC, muE, Nint)
            cds[i] = make(muCox, muExp)
        if np.max(np.abs(X - Xprev)) < tol:
            break
        Xprev = X.copy()
    return X, Q, U, C, it, cds


__all__ = ['pfqn_marie', 'marie_multi']
