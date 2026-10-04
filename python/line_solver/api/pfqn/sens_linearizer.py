"""
Approximate higher moments of the queue lengths of a closed product-form
queueing network, by differentiating the Linearizer fixed point. Polynomial in
the population, unlike the exact ``pfqn_sens_mom``.

Native Python implementation (no JPype / JVM dependency). Mirrors the MATLAB
reference ``pfqn_sens_linearizer.m``.

References:
    J. C. Strelen, "Moment Analysis for Closed Queuing Networks and its
    Linearizer", Performance Evaluation 11:127-142, 1990 (Section 5, equations
    (5.1)-(5.8), the CORE-2 and LINEARIZER-2 algorithms).
"""

import numpy as np


class PfqnSensLinearizer:
    """Result container for :func:`pfqn_sens_linearizer`.

    Attributes
    ----------
    X, Q, U, W : np.ndarray
        Approximate base measures. ``X`` is (1 x R); ``Q``, ``U`` and ``W`` are
        (M x R).
    m : np.ndarray (M,)
        Approximate ``E[Q_i]``, the total queue length at station ``i``.
    dm : np.ndarray (M, M)
        ``dm[i, h] = x_h dm_i/dx_h``, the scaled first derivative.
    d2m : np.ndarray (M,)
        ``d2m[i] = x_i^2 d^2m_i/dx_i^2``.
    Var, Cov, M2, M3, Skew : np.ndarray
        The moments of (3.2), formed exactly as in :func:`pfqn_sens_mom` but
        from the approximate derivatives. ``Var``, ``M2``, ``M3`` and ``Skew``
        are (M,); ``Cov`` is (M x M).
    CovAsym : float
        Raw asymmetry of ``Cov`` before symmetrization. Unlike the exact
        routines, this is NOT expected to sit at roundoff: the Linearizer fixed
        point does not enforce the symmetry that the product form guarantees, so
        this is a useful measure of the approximation error.
    iter : int
        Total CORE iterations performed.
    """

    def __init__(self, X, Q, U, W, m, dm, d2m, Var, Cov, M2, M3, Skew,
                 CovAsym, iter):
        self.X = X
        self.Q = Q
        self.U = U
        self.W = W
        self.m = m
        self.dm = dm
        self.d2m = d2m
        self.Var = Var
        self.Cov = Cov
        self.M2 = M2
        self.M3 = M3
        self.Skew = Skew
        self.CovAsym = CovAsym
        self.iter = iter


def _pack(X, Q, U, W, m, dm, d2m, it) -> PfqnSensLinearizer:
    """Form the moments of equation (3.2) from the approximate derivatives.

    ``Cov(i,j) = x_j dm_i/dx_j``. The product form makes this symmetric, but the
    Linearizer fixed point does not enforce that, so the raw asymmetry is a
    genuine error indicator here rather than a roundoff residual.
    """
    M = Q.shape[0]
    Cov = np.array(dm, dtype=np.float64)
    if Cov.size == 0:
        CovAsym = 0.0
    else:
        CovAsym = float(np.max(np.abs(Cov - Cov.T)))
    Cov = (Cov + Cov.T) / 2.0

    Var = np.zeros(M)
    M2 = np.zeros(M)
    M3 = np.zeros(M)
    Skew = np.zeros(M)
    for i in range(M):
        Var[i] = dm[i, i]
        M2[i] = dm[i, i] + m[i] ** 2
        M3[i] = d2m[i] + (1.0 + 3.0 * m[i]) * dm[i, i] + m[i] ** 3
        mu3 = M3[i] - 3.0 * m[i] * M2[i] + 2.0 * m[i] ** 3
        if Var[i] > 0:
            Skew[i] = mu3 / Var[i] ** 1.5
        else:
            Skew[i] = np.nan
    return PfqnSensLinearizer(X, Q, U, W, m, dm, d2m, Var, Cov, M2, M3, Skew,
                              CovAsym, it)


def _fractions(m, dm, d2m, n):
    """(5.1) and (5.5): ``v_i(l) = m_i(l)/n(l)``, and likewise for the
    derivatives."""
    M, R = m.shape
    v = np.zeros((M, R))
    dv = np.zeros((M, R, M))
    d2v = np.zeros((M, R, M))
    for l in range(R):
        if n[l] <= 0:
            continue
        v[:, l] = m[:, l] / n[l]
        dv[:, l, :] = dm[:, l, :] / n[l]
        d2v[:, l, :] = d2m[:, l, :] / n[l]
    return v, dv, d2v


def _core2(L, Z, n, delta, ddelta, d2delta, m, dm, d2m, tol, maxiter):
    """CORE-2 of the reference, extended to second derivatives.

    Iterates (5.1), (5.2) and the MVA equations (1.4)-(1.5) together with their
    first and second derivatives (5.5)-(5.6) and (3.4), until the mean queue
    lengths and the variances both stop moving.
    """
    M, R = L.shape
    nc = float(np.sum(n))
    if tol is None:
        tolm = 1.0 / (4000.0 + 16.0 * nc)     # termination test of the reference
    else:
        tolm = tol
    tolv = 1e-3

    lam = np.zeros(R)
    w = np.zeros((M, R))
    varprev = np.zeros(M)
    it = 0
    for iter_ in range(1, maxiter + 1):
        it = iter_
        mprev = m.copy()

        # ---- (5.1)-(5.2) and (5.5)-(5.6): queue lengths one job down -------
        # mtot(i,l) estimates sum_l2 m_i^{(n - e_l)}(l2)
        v, dv, d2v = _fractions(m, dm, d2m, n)
        mtot = np.zeros((M, R))
        dmtot = np.zeros((M, R, M))
        d2mtot = np.zeros((M, R, M))
        for l in range(R):
            if n[l] <= 0:
                continue
            for i in range(M):
                acc = 0.0
                dacc = np.zeros(M)
                d2acc = np.zeros(M)
                for l2 in range(R):
                    cnt = n[l2] - (1 if l2 == l else 0)      # (n - 1_l)_{l2}
                    if cnt <= 0:
                        continue
                    acc += cnt * (v[i, l2] + delta[i, l, l2])
                    for h in range(M):
                        dacc[h] += cnt * (dv[i, l2, h] + ddelta[i, l, l2, h])
                        d2acc[h] += cnt * (d2v[i, l2, h] + d2delta[i, l, l2, h])
                mtot[i, l] = acc
                for h in range(M):
                    dmtot[i, l, h] = dacc[h]
                    d2mtot[i, l, h] = d2acc[h]

        # ---- MVA (1.4)-(1.5) and its derivatives (3.4) ---------------------
        # w_i(l) = y_i * L(i,l) * (1 + mtot(i,l))
        w = np.zeros((M, R))
        dw = np.zeros((M, R, M))
        d2w = np.zeros((M, R, M))
        for l in range(R):
            if n[l] <= 0:
                continue
            for i in range(M):
                A = 1.0 + mtot[i, l]
                w[i, l] = L[i, l] * A
                for h in range(M):
                    dA = dmtot[i, l, h]
                    d2A = d2mtot[i, l, h]
                    if i == h:
                        dw[i, l, h] = L[i, l] * (A + dA)
                        d2w[i, l, h] = L[i, l] * (2.0 * dA + d2A)
                    else:
                        dw[i, l, h] = L[i, l] * dA
                        d2w[i, l, h] = L[i, l] * d2A
        lam = np.zeros(R)
        dlam = np.zeros((R, M))
        d2lam = np.zeros((R, M))
        for l in range(R):
            if n[l] <= 0:
                continue
            den = Z[l] + np.sum(w[:, l])
            lam[l] = n[l] / den
            for h in range(M):
                dden = np.sum(dw[:, l, h])
                d2den = np.sum(d2w[:, l, h])
                dlam[l, h] = -n[l] * dden / den ** 2
                d2lam[l, h] = (-n[l] * d2den / den ** 2
                               + 2.0 * n[l] * dden ** 2 / den ** 3)
        m = np.zeros((M, R))
        dm = np.zeros((M, R, M))
        d2m = np.zeros((M, R, M))
        for l in range(R):
            if n[l] <= 0:
                continue
            for i in range(M):
                m[i, l] = lam[l] * w[i, l]
                for h in range(M):
                    dm[i, l, h] = dlam[l, h] * w[i, l] + lam[l] * dw[i, l, h]
                    d2m[i, l, h] = (d2lam[l, h] * w[i, l]
                                    + 2.0 * dlam[l, h] * dw[i, l, h]
                                    + lam[l] * d2w[i, l, h])

        # ---- termination test of the reference -----------------------------
        dev = 0.0
        for l in range(R):
            if n[l] <= 0:
                continue
            dev = max(dev, float(np.max(np.abs(m[:, l] - mprev[:, l]))) / n[l])
        varnow = np.zeros(M)
        for i in range(M):
            varnow[i] = np.sum(dm[i, :, i])
        sv = float(np.sum(varnow))
        if sv > 0:
            vdev = float(np.max(np.abs(varnow - varprev))) / sv
        else:
            vdev = 0.0
        varprev = varnow
        if dev <= tolm and vdev <= tolv:
            break

    return m, dm, d2m, lam, w, it


def pfqn_sens_linearizer(L: np.ndarray, N: np.ndarray, Z: np.ndarray = None,
                         tol: float = None,
                         maxiter: int = 200) -> PfqnSensLinearizer:
    """Approximate moments E[Q_i], Var[Q_i], Cov[Q_i,Q_j], E[Q_i^2] and E[Q_i^3]
    of the per-station total queue lengths of a closed product-form queueing
    network, by the LINEARIZER-2 / LINEARIZER-3 algorithms of the reference
    (Section 5).

    Motivation. The exact moment analysis of :func:`pfqn_sens_mom` evaluates the
    MVA recursion on the whole population lattice, so it costs ``O(prod(N+1))``
    and is unusable once the populations are large. The Linearizer replaces that
    lattice by a fixed point over a handful of populations, and the reference
    observes that the same trick applies to the derivatives: differentiate the
    Linearizer equations, append the differentiated equations to the originals,
    and iterate all of them together. This routine does that, carrying both the
    first and the second derivative, so it returns everything (3.2) needs,
    including the third moment. Carrying only the first derivative is the
    reference's LINEARIZER-2; carrying the second as well is its LINEARIZER-3.

    The approximation. CORE (equations (5.1)-(5.2)) estimates the queue lengths
    at population ``n - 1_l`` from those at ``n`` by::

        v_i(l)          = m_i^(n)(l) / n(l)
        m_i^(n-1_l')(l) = (n - 1_l')_l * ( v_i(l) + delta_i(l',l) )

    and substitutes them into the exact MVA equations. Setting the delta terms to
    zero gives Bard-Schweitzer; Linearizer instead estimates them from (5.3),
    ``delta_i^(N)(l',l) = v_i^(N-1_l')(l) - v_i^(N)(l)``, by running CORE at each
    of the ``N - 1_l`` populations, and holds them fixed across populations (the
    heuristic (5.4)). Differentiating (5.1)-(5.3) gives (5.5)-(5.8), which are
    carried alongside.

    Accuracy. The reference reports, over 51 networks including 34 stress cases,
    relative errors below 2.1% on E[Q], 4.1% on E[Q^2] and 6.2% on E[Q^3]. The
    tests measure the error against the exact :func:`pfqn_sens_mom` on models
    small enough for both, and assert bands of that order rather than machine
    precision: this routine is an approximation and is expected to disagree with
    the exact answer.

    Parameters
    ----------
    L : (M, R) array
        Service demand matrix, ``L[i, r] = visits_ir / rate_ir``.
    N : (R,) array
        Population per class. Must be finite (closed populations only).
    Z : (R,) array, optional
        Think time per class (default zeros).
    tol : float, optional
        Convergence tolerance of the CORE fixed point on the mean queue lengths.
        Default: the test of the reference, ``1/(4000+16*sum(n))``, which is also
        applied to the variances at 1e-3.
    maxiter : int, optional
        Maximum CORE iterations (default 200).

    Returns
    -------
    PfqnSensLinearizer

    Notes
    -----
    Single-server stations plus an optional delay ``Z``, matching
    :func:`pfqn_linearizer`.

    The exact counterpart is :func:`pfqn_sens_mom`; the per-class exact second
    moments are in :func:`pfqn_sens_mva`.
    """
    L = np.atleast_2d(np.asarray(L, dtype=np.float64))
    Nf = np.asarray(N, dtype=np.float64).flatten()
    if np.any(np.isinf(Nf)):
        raise ValueError("pfqn_sens_linearizer requires a closed population")
    N = np.ceil(Nf).astype(int)
    R = len(N)
    if L.ndim == 1:
        L = L.reshape(-1, 1) if R == 1 else L.reshape(1, -1)
    M = L.shape[0]
    if L.shape[1] != R:
        raise ValueError("demand matrix and population vector have different "
                         "number of classes")
    if Z is None:
        Z = np.zeros(R)
    else:
        Z = np.asarray(Z, dtype=np.float64).flatten()

    if not np.any(N > 0):
        return _pack(np.zeros((1, R)), np.zeros((M, R)), np.zeros((M, R)),
                     np.zeros((M, R)), np.zeros(M), np.zeros((M, M)),
                     np.zeros(M), 0)

    # ---- Linearizer state -------------------------------------------------
    # pops: index 0 = N, index 1+l = N - e_l
    npops = 1 + R
    pv = np.zeros((npops, R), dtype=int)
    pv[0, :] = N
    for l in range(R):
        pv[1 + l, :] = N
        if N[l] > 0:
            pv[1 + l, l] = N[l] - 1

    # mE[p][i,l] = estimate of m_i^{(pv[p,:])}(l), and its derivatives w.r.t. y_h
    mE = []
    dmE = []
    d2mE = []
    for p in range(npops):
        mp = np.zeros((M, R))
        for l in range(R):
            mp[:, l] = pv[p, l] / M          # initialization of the reference
        mE.append(mp)
        dmE.append(np.zeros((M, R, M)))
        d2mE.append(np.zeros((M, R, M)))
    # delta[i,lp,l], indexed by the removed class lp and the class l
    delta = np.zeros((M, R, R))
    ddelta = np.zeros((M, R, R, M))
    d2delta = np.zeros((M, R, R, M))

    totiter = 0
    Xf = np.zeros((1, R))
    Wf = np.zeros((M, R))
    dmf = np.zeros((M, M))
    d2mf = np.zeros(M)

    for outer in range(1, 4):
        # ---- Step 1: CORE at the full population ------------------------
        mE[0], dmE[0], d2mE[0], lam, Wf, it = _core2(
            L, Z, N, delta, ddelta, d2delta, mE[0], dmE[0], d2mE[0], tol, maxiter)
        Xf = lam.reshape(1, R)
        totiter += it

        if outer < 3:
            # ---- Step 2: CORE at each of the N - e_l populations --------
            for l in range(R):
                if N[l] == 0:
                    continue
                mE[1 + l], dmE[1 + l], d2mE[1 + l], _, _, it = _core2(
                    L, Z, pv[1 + l, :], delta, ddelta, d2delta,
                    mE[1 + l], dmE[1 + l], d2mE[1 + l], tol, maxiter)
                totiter += it

            # ---- Step 3: refresh delta from (5.1) and (5.3) -------------
            vN, dvN, d2vN = _fractions(mE[0], dmE[0], d2mE[0], N)
            for lp in range(R):
                if N[lp] == 0:
                    continue
                vL, dvL, d2vL = _fractions(mE[1 + lp], dmE[1 + lp], d2mE[1 + lp],
                                           pv[1 + lp, :])
                for i in range(M):
                    for l in range(R):
                        delta[i, lp, l] = vL[i, l] - vN[i, l]
                        for h in range(M):
                            ddelta[i, lp, l, h] = dvL[i, l, h] - dvN[i, l, h]
                            d2delta[i, lp, l, h] = d2vL[i, l, h] - d2vN[i, l, h]

    # ---- final measures ---------------------------------------------------
    Qf = mE[0]
    m = np.sum(Qf, axis=1)
    for i in range(M):
        for h in range(M):
            dmf[i, h] = np.sum(dmE[0][i, :, h])
        d2mf[i] = np.sum(d2mE[0][i, :, i])
    U = np.zeros((M, R))
    for i in range(M):
        for r in range(R):
            U[i, r] = Xf[0, r] * L[i, r]

    return _pack(Xf, Qf, U, Wf, m, dmf, d2mf, totiter)


__all__ = [
    'pfqn_sens_linearizer',
    'PfqnSensLinearizer',
]
