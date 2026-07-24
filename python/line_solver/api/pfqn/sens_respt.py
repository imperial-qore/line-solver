"""
Exact moments of the sojourn time of a job at FCFS multiserver centers of a
closed product-form queueing network.

Native Python implementation (no JPype / JVM dependency). Mirrors the MATLAB
reference ``pfqn_sens_respt.m``.

References:
    J. C. Strelen, "Moment Analysis for Closed Queuing Networks and its
    Linearizer", Performance Evaluation 11:127-142, 1990 (Theorem 4.1 with
    equations (4.1)-(4.5) and Remarks 4.2-4.3).
"""

from math import factorial

import numpy as np


class PfqnSensRespt:
    """Result container for :func:`pfqn_sens_respt`.

    Attributes
    ----------
    X, Q, U : np.ndarray
        Base measures at population ``N``. ``X`` is (1 x R), ``Q`` and ``U`` are
        (M x R).
    W : np.ndarray (M, R)
        ``W[i, l] = E[W_(i,l)]``, the mean sojourn time per visit of a class-l
        job at station ``i``. Zero where class ``l`` does not visit ``i``.
    WM : np.ndarray (M, R, tmax)
        ``WM[i, l, t-1] = E[W_(i,l)^t]``.
    WVar : np.ndarray (M, R)
        ``Var[W_(i,l)] = E[W^2] - E[W]^2``. Requires ``tmax >= 2``.
    WSkew : np.ndarray (M, R)
        Skewness of ``W_(i,l)``. Requires ``tmax >= 3``; NaN if the variance is
        zero.
    m : np.ndarray (M,)
        ``E[Q_i]`` at population ``N``, the total queue length.
    Var : np.ndarray (M,)
        ``Var[Q_i]`` at population ``N``.
    p : np.ndarray (M, max(b))
        ``p[i, j] = P[Q_i = j]`` at population ``N``, for ``j = 0..b_i-1``.
        These are the only marginal probabilities the b-server recursion needs,
        so the matrix is RAGGED: row ``i`` is meaningful only up to column
        ``b_i`` and is zero-padded out to ``max(b)``. A padded entry is not
        ``P[Q_i = j]``; it is simply not computed. Read row ``i`` as
        ``p[i, :b[i]]``.
    Wresid : np.ndarray (M, R)
        The residence time ``w_i(l)`` of the MVA recursion. The identity
        ``W[i, l] = Wresid[i, l] / V[i, l]`` is an independent check of the
        ``t = 1`` case of (4.5).
    """

    def __init__(self, X, Q, U, m, Var, p, W, WM, Wresid, WVar, WSkew):
        self.X = X
        self.Q = Q
        self.U = U
        self.m = m
        self.Var = Var
        self.p = p
        self.W = W
        self.WM = WM
        self.Wresid = Wresid
        self.WVar = WVar
        self.WSkew = WSkew


def _jpow(j, tau):
    """``j**tau`` with the convention ``j^0 = 1``, so that ``0^0 = 1`` as the
    reference states."""
    if tau == 0:
        return 1.0
    return float(j) ** tau


def _strelen_a(b, mu, tmax):
    """Coefficients ``a_{t,tau}(0)`` of Remark 4.3 of the reference. They depend
    only on the number of servers ``b`` and on the per-server rate ``mu``, not on
    the network."""
    a = np.zeros((3, 4))
    a[0, 0] = (1.0 - b) / (b * mu)
    a[0, 1] = 1.0 / (b * mu)
    if tmax >= 2:
        a[1, 0] = (2.0 - b - b ** 2) / (b ** 2 * mu ** 2)
        a[1, 1] = 3.0 / (b ** 2 * mu ** 2)
        a[1, 2] = 1.0 / (b ** 2 * mu ** 2)
    if tmax >= 3:
        a[2, 0] = (6.0 - 5.0 * b + 3.0 * b ** 2 - 4.0 * b ** 3) / (b ** 3 * mu ** 3)
        a[2, 1] = (11.0 - 3.0 * b + 3.0 * b ** 2) / (b ** 3 * mu ** 3)
        a[2, 2] = 6.0 / (b ** 3 * mu ** 3)
        a[2, 3] = 1.0 / (b ** 3 * mu ** 3)
    return a


def _pack(X, Q, U, m, Var, p, W, WM, Wresid, tmax) -> PfqnSensRespt:
    M, R = Q.shape
    WVar = np.zeros((M, R))
    WSkew = np.zeros((M, R))
    if tmax >= 2:
        for i in range(M):
            for l in range(R):
                WVar[i, l] = WM[i, l, 1] - WM[i, l, 0] ** 2
    if tmax >= 3:
        for i in range(M):
            for l in range(R):
                mu3 = (WM[i, l, 2] - 3.0 * WM[i, l, 0] * WM[i, l, 1]
                       + 2.0 * WM[i, l, 0] ** 3)
                if WVar[i, l] > 0:
                    WSkew[i, l] = mu3 / WVar[i, l] ** 1.5
                else:
                    WSkew[i, l] = np.nan
    return PfqnSensRespt(X, Q, U, m, Var, p, W, WM, Wresid, WVar, WSkew)


def pfqn_sens_respt(S: np.ndarray, V: np.ndarray, N: np.ndarray,
                    Z: np.ndarray = None, b: np.ndarray = None,
                    tmax: int = 3) -> PfqnSensRespt:
    """Exact raw moments ``E[W_(i,l)^t]``, ``t = 1..tmax``, of the sojourn time of
    a class-l job at an FCFS b-server center ``i`` of a closed product-form
    queueing network, together with the variance of that sojourn time.

    This is Theorem 4.1 of the reference. Its mechanism is the arrival theorem of
    Lavenberg-Reiser and Sevcik-Mitrani: a class-l job arriving at center ``i``
    finds ``j`` jobs already there with probability ``p_i(j, N - 1_l)``.
    Conditioning the sojourn time on ``j`` and inverting the Laplace transform of
    the conditional density gives::

        E[W_(i,l)^t] = t!/mu^t + sum_{tau=0..t} a_(t,tau)(0) E[Qt_i^tau]
                       - sum_{j=0..b-1} p_i(j,N-1_l) sum_{tau=0..t} a_(t,tau)(0) j^tau

    where ``mu = 1/S(i)`` is the rate of each of the ``b`` servers, ``Qt_i`` is
    the total queue length at center ``i`` at population ``N - 1_l`` (so its
    moments are those of :func:`pfqn_sens_mom` evaluated one job down in class
    ``l``), and the coefficients ``a_(t,tau)(0)`` depend only on ``b`` and
    ``mu``, not on the network (Remark 4.3 of the reference). The moments
    ``E[Qt_i^tau]`` up to ``tau = 3`` need the second derivative of the MVA
    recursion, so this routine carries a second-order forward-mode pass exactly
    as :func:`pfqn_sens_mom` does, but over the b-server recursion (4.1)-(4.2)
    rather than the single-server one.

    For ``b = 1`` the coefficients ``a_(t,0)(0)`` vanish identically and the
    double-sum correction disappears, so no marginal probabilities are needed
    (Remark 4.2 of the reference); the routine still evaluates the general
    expression, which reduces to that case on its own.

    Only FCFS centers are covered. The reference is explicit that the
    sojourn-time distribution at PS and LCFS centers is in general not known, so
    no analogue exists there. FCFS in a BCMP network further requires the service
    time to be exponential and class-independent, which is why this routine takes
    a per-station service time ``S(i)`` and a separate visit-ratio matrix ``V``
    rather than a demand matrix: the sojourn time is per visit, so the per-visit
    rate ``mu = 1/S(i)`` must be known and cannot be recovered from the demand
    ``L(i,l) = S(i)*V(i,l)`` alone.

    Parameters
    ----------
    S : (M,) array
        Service time at each station, common to all classes.
    V : (M, R) array
        Visit ratio matrix. The demand is ``L[i, r] = S[i] * V[i, r]``.
    N : (R,) array
        Population per class. Must be finite (closed populations only).
    Z : (R,) array, optional
        Think time per class (default zeros).
    b : (M,) array, optional
        Number of servers at each station (default ones).
    tmax : int, optional
        Highest sojourn-time moment to return, 1..3 (default 3). The
        coefficients ``a_(t,tau)(0)`` are tabulated in the reference up to
        ``t = 3``.

    Returns
    -------
    PfqnSensRespt

    Notes
    -----
    Restricted to closed populations. Load-dependent rates are not covered here;
    the b-server dependence is the only state dependence, and it is carried
    exactly by (4.1)-(4.2).
    """
    V = np.atleast_2d(np.asarray(V, dtype=np.float64))
    S = np.asarray(S, dtype=np.float64).flatten()
    Nf = np.asarray(N, dtype=np.float64).flatten()
    if np.any(np.isinf(Nf)):
        raise ValueError("pfqn_sens_respt requires a closed population")
    N = np.ceil(Nf).astype(int)
    R = len(N)
    M = V.shape[0]
    if V.shape[1] != R:
        raise ValueError("visit matrix and population vector have different "
                         "number of classes")
    if Z is None:
        Z = np.zeros(R)
    else:
        Z = np.asarray(Z, dtype=np.float64).flatten()
    if b is None:
        b = np.ones(M, dtype=int)
    else:
        b = np.round(np.asarray(b, dtype=np.float64).flatten()).astype(int)
    if tmax < 1 or tmax > 3:
        raise ValueError("tmax must be 1, 2 or 3: the coefficients a_{t,tau}(0) "
                         "are tabulated in the reference up to order three.")
    if np.any(b < 1):
        raise ValueError("the number of servers must be at least one at every "
                         "station")
    if np.any(S <= 0):
        raise ValueError("every FCFS station must have a strictly positive "
                         "service time")

    rho = np.zeros((M, R))
    for i in range(M):
        for l in range(R):
            rho[i, l] = S[i] * V[i, l]
    mu = 1.0 / S
    bmax = int(np.max(b))

    X = np.zeros((1, R))
    Q = np.zeros((M, R))
    U = np.zeros((M, R))
    m = np.zeros(M)
    W = np.zeros((M, R))
    WM = np.zeros((M, R, tmax))
    Wresid = np.zeros((M, R))
    pN = np.zeros((M, bmax))

    if not np.any(N > 0):
        return _pack(X, Q, U, m, np.zeros(M), pN, W, WM, Wresid, tmax)

    # ---- population lattice ----------------------------------------------
    prods = np.zeros(max(R - 1, 0))
    for w in range(R - 1):
        prods[w] = np.prod(np.ones(R - w - 1) + N[w + 1:])
    first_non_empty = R - 1
    while N[first_non_empty] == 0:
        first_non_empty -= 1
    totpop = int(np.prod(N + 1))
    ctr = totpop

    # see _kb/03-api-layer.md for rationale
    Mrow = np.zeros((totpop, M))
    D1m = np.zeros((totpop, M, M))
    D2m = np.zeros((totpop, M, M))
    Prow = np.zeros((totpop, M, bmax))
    D1p = np.zeros((totpop, M, bmax, M))
    D2p = np.zeros((totpop, M, bmax, M))
    Prow[0, :, 0] = 1.0   # empty population: every station holds zero jobs

    currentpop = 1
    n = np.zeros(R, dtype=int)
    n[first_non_empty] = 1
    rows = np.zeros(R, dtype=int)

    while ctr > 0:
        hnvec = currentpop
        # ---- residence times, eq. (4.1), and their first two derivatives ---
        wv = np.zeros((M, R))
        d1w = np.zeros((M, R, M))
        d2w = np.zeros((M, R, M))
        for s in range(R):
            if n[s] > 0:
                n[s] -= 1
                pos = int(n[R - 1])
                for w in range(R - 1):
                    pos += int(n[w] * prods[w])
                n[s] += 1
            else:
                pos = 0
            rows[s] = pos
            if n[s] == 0:
                continue   # w and every derivative stay zero, as does X[s]
            for i in range(M):
                # bracket = 1 + m_i(n-e_s) + sum_{j=0}^{b_i-2} (b_i-1-j) p_i(j,n-e_s)
                brk = 1.0 + Mrow[pos, i]
                for j in range(b[i] - 1):
                    brk += (b[i] - 1 - j) * Prow[pos, i, j]
                wv[i, s] = (rho[i, s] / b[i]) * brk
                for h in range(M):
                    dbrk = D1m[pos, i, h]
                    d2brk = D2m[pos, i, h]
                    for j in range(b[i] - 1):
                        dbrk += (b[i] - 1 - j) * D1p[pos, i, j, h]
                        d2brk += (b[i] - 1 - j) * D2p[pos, i, j, h]
                    # w = y_i * (rho/b) * brk
                    if i == h:
                        d1w[i, s, h] = (rho[i, s] / b[i]) * (brk + dbrk)
                        d2w[i, s, h] = (rho[i, s] / b[i]) * (2.0 * dbrk + d2brk)
                    else:
                        d1w[i, s, h] = (rho[i, s] / b[i]) * dbrk
                        d2w[i, s, h] = (rho[i, s] / b[i]) * d2brk

        # ---- throughputs and their derivatives -----------------------------
        lam = np.zeros(R)
        d1lam = np.zeros((R, M))
        d2lam = np.zeros((R, M))
        for s in range(R):
            if n[s] == 0:
                continue
            den = Z[s] + np.sum(wv[:, s])
            lam[s] = n[s] / den
            for h in range(M):
                dden = np.sum(d1w[:, s, h])
                d2den = np.sum(d2w[:, s, h])
                d1lam[s, h] = -n[s] * dden / den ** 2
                d2lam[s, h] = (-n[s] * d2den / den ** 2
                               + 2.0 * n[s] * dden ** 2 / den ** 3)

        # ---- mean queue lengths --------------------------------------------
        for i in range(M):
            acc = 0.0
            for s in range(R):
                if n[s] == 0:
                    continue
                acc += lam[s] * wv[i, s]
            Mrow[hnvec, i] = acc
            for h in range(M):
                d1acc = 0.0
                d2acc = 0.0
                for s in range(R):
                    if n[s] == 0:
                        continue
                    d1acc += d1lam[s, h] * wv[i, s] + lam[s] * d1w[i, s, h]
                    d2acc += (d2lam[s, h] * wv[i, s]
                              + 2.0 * d1lam[s, h] * d1w[i, s, h]
                              + lam[s] * d2w[i, s, h])
                D1m[hnvec, i, h] = d1acc
                D2m[hnvec, i, h] = d2acc

        # ---- marginal probabilities, eq. (4.2), and their derivatives -------
        nc = int(np.sum(n))
        for i in range(M):
            # p_i(j,n) = (1/j) sum_l lam(l) * rho_i(l)*y_i * p_i(j-1, n-e_l)
            for j in range(1, b[i]):
                if j > nc:
                    Prow[hnvec, i, j] = 0.0
                    continue
                acc = 0.0
                for l in range(R):
                    if n[l] == 0:
                        continue
                    acc += lam[l] * rho[i, l] * Prow[rows[l], i, j - 1]
                Prow[hnvec, i, j] = acc / j
                for h in range(M):
                    d1acc = 0.0
                    d2acc = 0.0
                    for l in range(R):
                        if n[l] == 0:
                            continue
                        pprev = Prow[rows[l], i, j - 1]
                        d1prev = D1p[rows[l], i, j - 1, h]
                        d2prev = D2p[rows[l], i, j - 1, h]
                        # g = rho * u * v with u = lam, v = y_i * pprev
                        if i == h:
                            v1 = pprev + d1prev
                            v2 = 2.0 * d1prev + d2prev
                        else:
                            v1 = d1prev
                            v2 = d2prev
                        d1acc += rho[i, l] * (d1lam[l, h] * pprev + lam[l] * v1)
                        d2acc += rho[i, l] * (d2lam[l, h] * pprev
                                              + 2.0 * d1lam[l, h] * v1
                                              + lam[l] * v2)
                    D1p[hnvec, i, j, h] = d1acc / j
                    D2p[hnvec, i, j, h] = d2acc / j
            # u_i = sum_l lam(l) * rho_i(l)*y_i  (mean number of busy servers)
            ui = 0.0
            d1ui = np.zeros(M)
            d2ui = np.zeros(M)
            for l in range(R):
                if n[l] == 0:
                    continue
                ui += lam[l] * rho[i, l]
                for h in range(M):
                    if i == h:
                        d1ui[h] += rho[i, l] * (d1lam[l, h] + lam[l])
                        d2ui[h] += rho[i, l] * (d2lam[l, h] + 2.0 * d1lam[l, h])
                    else:
                        d1ui[h] += rho[i, l] * d1lam[l, h]
                        d2ui[h] += rho[i, l] * d2lam[l, h]
            # p_i(0,n) = 1 - (1/b)(u_i + sum_{j=1}^{b-1} (b-j) p_i(j,n))
            acc0 = ui
            for j in range(1, b[i]):
                acc0 += (b[i] - j) * Prow[hnvec, i, j]
            Prow[hnvec, i, 0] = 1.0 - acc0 / b[i]
            for h in range(M):
                d1acc0 = d1ui[h]
                d2acc0 = d2ui[h]
                for j in range(1, b[i]):
                    d1acc0 += (b[i] - j) * D1p[hnvec, i, j, h]
                    d2acc0 += (b[i] - j) * D2p[hnvec, i, j, h]
                D1p[hnvec, i, 0, h] = -d1acc0 / b[i]
                D2p[hnvec, i, 0, h] = -d2acc0 / b[i]

        # keep the measures of the last (full) population
        X[0, :] = lam
        for i in range(M):
            for s in range(R):
                Wresid[i, s] = wv[i, s]
                Q[i, s] = lam[s] * wv[i, s]
                U[i, s] = lam[s] * rho[i, s]

        # ---- odometer advance -----------------------------------------------
        s = R - 1
        while s >= 0 and (n[s] == N[s] or s > first_non_empty):
            s -= 1
        if s < 0:
            break
        n[s] += 1
        for i in range(s + 1, R):
            n[i] = 0
        ctr -= 1
        currentpop += 1

    lastrow = currentpop
    for i in range(M):
        m[i] = Mrow[lastrow, i]
        for j in range(bmax):
            pN[i, j] = Prow[lastrow, i, j]
    Var = np.zeros(M)
    for i in range(M):
        Var[i] = D1m[lastrow, i, i]

    # ---- sojourn-time moments, eq. (4.5) ---------------------------------
    # index of N - e_l on the lattice
    rowsN = np.zeros(R, dtype=int)
    for l in range(R):
        if N[l] > 0:
            nn = N.copy()
            nn[l] -= 1
            pos = int(nn[R - 1])
            for w in range(R - 1):
                pos += int(nn[w] * prods[w])
            rowsN[l] = pos

    for i in range(M):
        for l in range(R):
            if N[l] == 0 or V[i, l] <= 0:
                continue
            rl = rowsN[l]
            # moments of the queue length seen by an arriving class-l job, i.e.
            # of the total queue at station i at population N - e_l, from (3.2)
            mt = Mrow[rl, i]
            d1t = D1m[rl, i, i]
            d2t = D2m[rl, i, i]
            EQ = np.zeros(4)               # EQ[tau] = E[Qt_i^tau]
            EQ[0] = 1.0
            EQ[1] = mt
            EQ[2] = d1t + mt ** 2
            EQ[3] = d2t + (1.0 + 3.0 * mt) * d1t + mt ** 3
            acoef = _strelen_a(b[i], mu[i], tmax)
            for t in range(1, tmax + 1):
                val = factorial(t) / mu[i] ** t
                for tau in range(t + 1):
                    val += acoef[t - 1, tau] * EQ[tau]
                # correction over the states in which a server is idle
                for j in range(b[i]):
                    inner = 0.0
                    for tau in range(t + 1):
                        inner += acoef[t - 1, tau] * _jpow(j, tau)
                    val -= Prow[rl, i, j] * inner
                WM[i, l, t - 1] = val
            W[i, l] = WM[i, l, 0]

    return _pack(X, Q, U, m, Var, pN, W, WM, Wresid, tmax)


__all__ = [
    'pfqn_sens_respt',
    'PfqnSensRespt',
]
