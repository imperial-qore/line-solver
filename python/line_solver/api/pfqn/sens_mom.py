"""
Exact higher moments (up to order three) of the per-station total queue lengths
of a closed product-form queueing network, by second-order differentiation of
the MVA recursion.

Native Python implementation (no JPype / JVM dependency). Mirrors the MATLAB
reference ``pfqn_sens_mom.m``.

References:
    J. C. Strelen, "Moment Analysis for Closed Queuing Networks and its
    Linearizer", Performance Evaluation 11:127-142, 1990 (Theorems 2.1, 3.1,
    3.2, 3.5 and equation (3.2)).
"""

import numpy as np


class PfqnSensMom:
    """Result container for :func:`pfqn_sens_mom`.

    Attributes
    ----------
    X, Q, U, R : np.ndarray
        Base MVA measures, identical entry by entry to ``pfqn_mva(L, N, Z, mi)``.
        ``X`` is (1 x R), ``Q``, ``U`` and ``R`` are (M x R).
    m : np.ndarray (M, G)
        ``m[i, g] = E[Q_(i,g)]``, the mean queue length of group ``g`` at station
        ``i``. COLLAPSED to (M,) in the default single-group case, where it is
        the per-station total.
    dm : np.ndarray (M, G, M, G)
        ``dm[i, g, j, g2]``, the scaled first derivative of ``m_(i,g)`` with
        respect to the parameter of ``(j, g2)``. Collapsed to (M, M) when G == 1.
    d2m : np.ndarray (M, G)
        ``d2m[i, g]``, the scaled pure second derivative with respect to the
        parameter of ``(i, g)``. Only that entry is needed by (3.2); the mixed
        second derivatives are not required for moments of a single ``Q_(i,g)``
        and would cost an extra factor ``M*G`` to carry. Collapsed to (M,) when
        G == 1.
    Var : np.ndarray (M, G)
        ``Var[Q_(i,g)]``. Collapsed to (M,) when G == 1.
    Cov : np.ndarray (M, G, M, G)
        ``Cov[i, g, j, g2] = Cov[Q_(i,g), Q_(j,g2)]``. Collapsed to (M, M) when
        G == 1, where the group index carries no information.
    M2 : np.ndarray (M, G)
        ``E[Q_(i,g)^2]``. Collapsed to (M,) when G == 1.
    M3 : np.ndarray (M, G)
        ``E[Q_(i,g)^3]``. Collapsed to (M,) when G == 1.
    Skew : np.ndarray (M, G)
        Skewness of ``Q_(i,g)``, i.e. the third central moment divided by
        ``Var^(3/2)``. NaN where the variance is zero (a deterministic queue
        length). Collapsed to (M,) when G == 1.
    CovAsym : float
        ``max |x_j dm_i/dx_j - x_i dm_j/dx_i|`` before symmetrization of ``Cov``.
        The two are distinct expressions that must agree, so this is a live
        residual of the recursion; expect roundoff.
    """

    def __init__(self, X, Q, U, R, m, dm, d2m, Var, Cov, M2, M3, Skew, CovAsym):
        self.X = X
        self.Q = Q
        self.U = U
        self.R = R
        self.m = m
        self.dm = dm
        self.d2m = d2m
        self.Var = Var
        self.Cov = Cov
        self.M2 = M2
        self.M3 = M3
        self.Skew = Skew
        self.CovAsym = CovAsym


def _pack(X, Q, U, C, m, dm, d2m, G) -> PfqnSensMom:
    """Symmetrize the scaled first derivative into a covariance matrix and form
    the moments of equation (3.2).

    ``Cov((i,g),(j,g'))`` and its transpose are distinct expressions for the same
    quantity, obtained along different derivative tracks of the recursion; the
    raw disagreement is reported in ``CovAsym`` rather than discarded, and only
    then is ``Cov`` symmetrized so that it is a covariance matrix.

    ``m``, ``dm`` and ``d2m`` arrive shaped (M, G), (M, G, M, G) and (M, G). When
    G == 1 every returned field is collapsed along the group axis, so that the
    default (per-station total) call keeps the shapes it has always had.

    The (M, G) -> M*G flattening is row-major here and column-major in MATLAB,
    but the SAME relabeling is applied to both axes of the square matrix, so the
    transpose, ``CovAsym`` and the symmetrization are all invariant under it and
    the reshaped-back 4-D array is identical.
    """
    M = Q.shape[0]
    flat = np.array(dm, dtype=np.float64).reshape(M * G, M * G)
    if flat.size == 0:
        CovAsym = 0.0
    else:
        CovAsym = float(np.max(np.abs(flat - flat.T)))
    flat = (flat + flat.T) / 2.0

    Var = np.zeros((M, G))
    M2 = np.zeros((M, G))
    M3 = np.zeros((M, G))
    Skew = np.zeros((M, G))
    for i in range(M):
        for g in range(G):
            d1 = dm[i, g, i, g]
            Var[i, g] = d1                                          # (3.2)
            M2[i, g] = d1 + m[i, g] ** 2                            # (3.2)
            M3[i, g] = (d2m[i, g] + (1.0 + 3.0 * m[i, g]) * d1
                        + m[i, g] ** 3)                             # (3.2)
            # third central moment mu3 = E[Q^3] - 3 m E[Q^2] + 2 m^3
            mu3 = M3[i, g] - 3.0 * m[i, g] * M2[i, g] + 2.0 * m[i, g] ** 3
            if Var[i, g] > 0:
                Skew[i, g] = mu3 / Var[i, g] ** 1.5
            else:
                Skew[i, g] = np.nan

    if G == 1:
        # the group index carries no information here; collapse it so the
        # default call keeps the per-station-total shapes
        Cov = flat.reshape(M, M)
        dm_out = np.array(dm, dtype=np.float64).reshape(M, M)
        m_out = m.reshape(M)
        d2m_out = d2m.reshape(M)
        Var = Var.reshape(M)
        M2 = M2.reshape(M)
        M3 = M3.reshape(M)
        Skew = Skew.reshape(M)
    else:
        Cov = flat.reshape(M, G, M, G)
        dm_out = np.array(dm, dtype=np.float64)
        m_out = m
        d2m_out = d2m
    return PfqnSensMom(X, Q, U, C, m_out, dm_out, d2m_out, Var, Cov, M2, M3,
                       Skew, CovAsym)


def pfqn_sens_mom(L: np.ndarray, N: np.ndarray, Z: np.ndarray = None,
                  mi: np.ndarray = None,
                  groups: np.ndarray = None) -> PfqnSensMom:
    """Exact moments E[Q_i], E[Q_i^2], E[Q_i^3] and the covariances Cov[Q_i,Q_j]
    of the TOTAL queue lengths ``Q_i = sum_r n(i,r)`` of a closed product-form
    (BCMP) queueing network.

    The method is the moment analysis of Strelen. Its Theorem 3.1 states that one
    further factor ``Q_i`` in a moment costs one differentiation with respect to
    ``x_i``, the reciprocal of the capacity of station ``i``, i.e. a parameter
    that scales the service times of ALL classes at station ``i``::

        E[Q_i^j] = m_i E[Q_i^(j-1)] + x_i d/dx_i E[Q_i^(j-1)],   Q_i^0 = 1

    Iterating from ``E[Q_i^0] = 1`` gives, with ``m_i = E[Q_i]``, equation (3.2)::

        Var[Q_i]     = x_i dm_i/dx_i
        Cov[Q_i,Q_j] = x_j dm_i/dx_j = x_i dm_j/dx_i
        E[Q_i^2]     = x_i dm_i/dx_i + m_i^2
        E[Q_i^3]     = x_i^2 d^2m_i/dx_i^2 + (x_i + 3 x_i m_i) dm_i/dx_i + m_i^3

    so the third moment requires the SECOND derivative of the MVA recursion,
    which is what this routine adds over :func:`pfqn_sens_mva` and
    :func:`pfqn_sens` (both first order only). The derivatives are obtained by
    second-order forward-mode differentiation of the Reiser-Lavenberg recursion,
    i.e. by carrying, for each parameter, the value together with its first and
    second derivative along the population lattice (Theorem 3.2 for one class,
    Theorem 3.5 for several).

    The parameter need not scale a whole column. Theorem 1 of Akyildiz and
    Strelen states the same recursion for a parameter that scales the service
    times of an arbitrary class subset ``T`` at station ``i``, and the moments it
    generates are then those of ``Q_(i,T) = sum_(r in T) n(i,r)``. ``groups``
    supplies that subset structure: it partitions the classes, and the routine
    reports the moments of each group's queue length at each station. The three
    useful settings are::

        groups = ones(R)   the whole column: per-station TOTALS (default, this
                           is Strelen's x_i)
        groups = 1..R      one class per group: PER-CLASS moments, so that even
                           the third moment is per class
        groups = chain(r)  one group per chain: PER-CHAIN moments

    Strelen states only the first; the generalization is Akyildiz and Strelen's
    ``T``. The per-class setting reproduces :func:`pfqn_sens_mva`'s second
    moments exactly.

    Parameters
    ----------
    L : (M, R) array
        Service demand matrix, ``L[i, r] = visits_ir / rate_ir``.
    N : (R,) array
        Population per class. Must be finite (closed populations only).
    Z : (R,) array, optional
        Think time per class (default zeros).
    mi : (M,) array, optional
        Server multiplicity per station (default ones).
    groups : (R,) array, optional
        Class-to-group map, a partition of the classes into ``G = max(groups)``
        groups labelled consecutively 1..G with no empty group. Default
        ``ones(R)``, i.e. one group holding every class, the per-station total.
        Labels are 1-based, as in MATLAB.

    Returns
    -------
    PfqnSensMom

    Notes
    -----
    Restricted to closed populations, as is the moment analysis of the
    reference. Mixed and load-dependent second moments are in
    :func:`pfqn_sens_mvaldmx`.

    Moments of the sojourn times at FCFS centers are built on top of these
    queue-length moments by :func:`pfqn_sens_respt`, following Theorem 4.1.

    The exact recursion costs ``O(prod(N+1))`` lattice points;
    :func:`pfqn_sens_linearizer` approximates the same quantities in polynomial
    time.
    """
    L = np.atleast_2d(np.asarray(L, dtype=np.float64))
    Nf = np.asarray(N, dtype=np.float64).flatten()
    if np.any(np.isinf(Nf)):
        raise ValueError("pfqn_sens_mom requires a closed population")
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
    if mi is None:
        mi = np.ones(M)
    else:
        mi = np.asarray(mi, dtype=np.float64).flatten()
    if groups is None:
        groups = np.ones(R, dtype=int)
    else:
        groups = np.round(np.asarray(groups, dtype=np.float64).flatten()
                          ).astype(int)
    if len(groups) != R or np.any(groups < 1):
        raise ValueError("groups must be a (1 x R) vector of group labels "
                         "starting at 1")
    G = int(np.max(groups))
    if not np.array_equal(np.unique(groups), np.arange(1, G + 1)):
        raise ValueError("groups must label the classes consecutively from 1 to "
                         "max(groups), with no empty group")
    # 0-based group index for internal use; the API is 1-based like MATLAB
    g_of = groups - 1

    X = np.zeros((1, R))
    Q = np.zeros((M, R))
    U = np.zeros((M, R))
    C = np.zeros((M, R))
    m = np.zeros((M, G))
    dm = np.zeros((M, G, M, G))
    d2m = np.zeros((M, G))

    if not np.any(N > 0):
        return _pack(X, Q, U, C, m, dm, d2m, G)

    # ---- population-lattice odometer, identical to pfqn_mva ---------------
    prods = np.zeros(max(R - 1, 0))
    for w in range(R - 1):
        prods[w] = np.prod(np.ones(R - w - 1) + N[w + 1:])
    first_non_empty = R - 1
    while N[first_non_empty] == 0:
        first_non_empty -= 1
    totpop = int(np.prod(N + 1))
    ctr = totpop

    # see _kb/03-api-layer.md for rationale
    P = M * G
    pidx = np.zeros((M, G), dtype=int)
    for h in range(M):
        for g in range(G):
            pidx[h, g] = h * G + g
    Qtot = np.zeros((totpop, M))
    D1 = np.zeros((totpop, M, P))
    D2 = np.zeros((totpop, M, P))
    Qg = np.zeros((M, G))
    D1Qg = np.zeros((M, G, P))
    D2Qg = np.zeros((M, G, P))
    currentpop = 1
    n = np.zeros(R, dtype=int)
    n[first_non_empty] = 1

    Cs = np.zeros(M)
    dCs = np.zeros((M, P))
    d2Cs = np.zeros((M, P))

    while ctr > 0:
        # the group accumulators describe one population only
        Qg[:] = 0.0
        D1Qg[:] = 0.0
        D2Qg[:] = 0.0
        for s in range(R):
            if n[s] > 0:
                n[s] -= 1
                pos = int(n[R - 1])
                for w in range(R - 1):
                    pos += int(n[w] * prods[w])
                n[s] += 1
            else:
                # see _kb/03-api-layer.md for rationale
                pos = 0
            row = pos

            # see _kb/03-api-layer.md for rationale
            gs = g_of[s]
            CNtot = 0.0
            dCNtot = np.zeros(P)
            d2CNtot = np.zeros(P)
            for i in range(M):
                A = mi[i] + Qtot[row, i]
                Cs[i] = L[i, s] * A
                C[i, s] = Cs[i]
                CNtot += Cs[i]
                own = pidx[i, gs]
                for p in range(P):
                    dA = D1[row, i, p]
                    d2A = D2[row, i, p]
                    if p == own:
                        dCs[i, p] = L[i, s] * (A + dA)
                        d2Cs[i, p] = L[i, s] * (2.0 * dA + d2A)
                    else:
                        dCs[i, p] = L[i, s] * dA
                        d2Cs[i, p] = L[i, s] * d2A
                    dCNtot[p] += dCs[i, p]
                    d2CNtot[p] += d2Cs[i, p]

            # see _kb/03-api-layer.md for rationale
            den = Z[s] + CNtot
            X[0, s] = n[s] / den if den > 0 else 0.0
            dX = np.zeros(P)
            d2X = np.zeros(P)
            if den > 0:
                for p in range(P):
                    dX[p] = -n[s] * dCNtot[p] / den ** 2
                    d2X[p] = (-n[s] * d2CNtot[p] / den ** 2
                              + 2.0 * n[s] * dCNtot[p] ** 2 / den ** 3)

            # ---- queue lengths ---------------------------------------------
            # Q = X*C,  dQ = dX*C + X*dC,  d2Q = d2X*C + 2*dX*dC + X*d2C
            for i in range(M):
                Q[i, s] = X[0, s] * Cs[i]
                Qtot[currentpop, i] += Q[i, s]
                Qg[i, gs] += Q[i, s]
                for p in range(P):
                    dQ = dX[p] * Cs[i] + X[0, s] * dCs[i, p]
                    d2Q = (d2X[p] * Cs[i] + 2.0 * dX[p] * dCs[i, p]
                           + X[0, s] * d2Cs[i, p])
                    D1[currentpop, i, p] += dQ
                    D2[currentpop, i, p] += d2Q
                    D1Qg[i, gs, p] += dQ
                    D2Qg[i, gs, p] += d2Q

        # ---- odometer advance ---------------------------------------------
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

    # ---- utilization -----------------------------------------------------
    for i in range(M):
        for r in range(R):
            U[i, r] = X[0, r] * L[i, r]

    # see _kb/03-api-layer.md for rationale
    for i in range(M):
        for g in range(G):
            m[i, g] = Qg[i, g]
            for j in range(M):
                for g2 in range(G):
                    dm[i, g, j, g2] = D1Qg[i, g, pidx[j, g2]]
            d2m[i, g] = D2Qg[i, g, pidx[i, g]]

    return _pack(X, Q, U, C, m, dm, d2m, G)


__all__ = [
    'pfqn_sens_mom',
    'PfqnSensMom',
]
