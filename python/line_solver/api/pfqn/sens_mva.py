"""
Exact per-station queue-length variances and covariances for closed
product-form queueing networks, computed by an MVA-like moment recursion that
does not require the full sensitivity Jacobian.

Native Python implementation (no JPype / JVM dependency). Mirrors the MATLAB
reference ``pfqn_sens_mva.m``.

References:
    E. de Souza e Silva and R. R. Muntz, "Simple Relationships Among Moments of
    Queue Lengths in Product Form Queueing Networks", IEEE Trans. Computers
    37(9):1125-1129, 1988 (Theorems 1-2 and Corollary 1).
    I. F. Akyildiz and J. C. Strelen, "Moment Analysis for Load-Dependent Mixed
    Product Form Queueing Networks", IEEE Trans. Communications 39(6):828-832,
    1991 (Theorem 1, whose k=2 case is the identity
    Cov[n(i,r),n(j,s)] = theta(i,r) dQ(j,s)/dtheta(i,r)).
"""

import numpy as np


class PfqnSensMva:
    """Result container for :func:`pfqn_sens_mva`.

    Attributes
    ----------
    X, Q, U, R : np.ndarray
        Base MVA measures, identical entry by entry to ``pfqn_mva(L, N, Z, mi)``.
        ``X`` is (1 x R) system throughput per class, ``Q`` (M x R) mean queue
        length, ``U`` (M x R) utilization, ``R`` (M x R) residence time per
        visit.
    QCov : np.ndarray (M x R x R)
        ``QCov[i, r, s] = Cov[n(i,r), n(i,s)]``, the queue-length covariance of
        classes ``r`` and ``s`` at station ``i``. Symmetric in ``(r, s)``.
    QVar : np.ndarray (M x R)
        ``QVar[i, r] = QCov[i, r, r] = Var[n(i,r)]``.
    QTotVar : np.ndarray (M,)
        ``QTotVar[i] = Var[sum_r n(i,r)]``, i.e. ``sum_{r,s} QCov[i, r, s]``.
        This is Theorem 3 of de Souza e Silva and Muntz, obtained here without a
        capacity derivative.
    QCovAsym : float
        ``max |W(r,s) - W(s,r)|`` over the covariance entries before
        symmetrization. The two triangles come from differentiating two
        different classes' MVA equations, so this is an independent residual of
        the recursion and should sit at roundoff; a large value signals a bug.
    """

    def __init__(self, X, Q, U, R, QCov, QVar, QTotVar, QCovAsym):
        self.X = X
        self.Q = Q
        self.U = U
        self.R = R
        self.QCov = QCov
        self.QVar = QVar
        self.QTotVar = QTotVar
        self.QCovAsym = QCovAsym


def _pack(X, Q, U, C, QCovRaw) -> PfqnSensMva:
    """Symmetrize the raw covariance block and derive QVar / QTotVar.

    The recursion obtains ``W(k,j;t,j)`` by differentiating class ``t``'s MVA
    equation and ``W(t,j;k,j)`` by differentiating class ``k``'s, so the two
    triangles are numerically distinct expressions that agree only up to
    roundoff. Averaging them keeps QCov exactly symmetric, as a covariance
    matrix must be. The raw discrepancy is reported in ``QCovAsym`` rather than
    discarded, so that a genuine disagreement cannot hide behind the averaging.
    """
    M, R = Q.shape
    QCov = (QCovRaw + np.transpose(QCovRaw, (0, 2, 1))) / 2.0
    if QCovRaw.size == 0:
        QCovAsym = 0.0
    else:
        QCovAsym = float(np.max(np.abs(QCovRaw - np.transpose(QCovRaw, (0, 2, 1)))))
    QVar = np.zeros((M, R))
    QTotVar = np.zeros(M)
    for i in range(M):
        for r in range(R):
            QVar[i, r] = QCov[i, r, r]
        QTotVar[i] = np.sum(QCov[i, :, :])
    return PfqnSensMva(X, Q, U, C, QCov, QVar, QTotVar, QCovAsym)


def pfqn_sens_mva(L: np.ndarray, N: np.ndarray, Z: np.ndarray = None,
                 mi: np.ndarray = None) -> PfqnSensMva:
    """Exact second moments of the queue lengths of a closed product-form network.

    The moments are obtained by an MVA-type recursion evaluated on the same
    population lattice as :func:`pfqn_mva`, so no derivative of the model is
    ever formed and the cost is O(M*R^2) per lattice point rather than the
    O(M^2*R^2) of the differentiated-MVA kernel used by :func:`pfqn_sens`.

    The recursion is obtained by differentiating the Reiser-Lavenberg MVA
    equation ``Q(j,v|N) = X(v|N) L(j,v) (mi(j) + Qtot(j|N-e_v))`` with respect
    to the visit ratio ``theta(i,k)`` of class ``k`` at station ``i`` and
    rescaling. Writing ``W(k,i;v,j|N) = Cov[n(i,k), n(j,v)]`` at population N,

        W(k,i;v,j|N) = Q(j,v|N) (Q(i,k|N-e_v) - Q(i,k|N))
                     + [i==j & k==v] Q(j,v|N)
                     + X(v|N) L(j,v) sum_t W(k,i;t,j|N-e_v)

    with ``W(.|0) = 0``. This routine evaluates the same-station case ``i==j``,
    which is self-contained: the inner sum then only involves same-station
    terms, so a single scalar ``Ssum(j,k|N) = sum_t W(k,j;t,j|N)`` carried along
    the lattice closes the recursion. The cross-station case ``i != j`` is not
    self-contained (it couples every station pair) and costs as much as the full
    Jacobian, so it is left to :func:`pfqn_sens`.

    Setting ``i == j`` and ``mi == 1`` reproduces Corollary 1 of de Souza e
    Silva and Muntz, i.e. its equations (2.9a) for the variance and (2.10) for
    the covariance; the station multiplicity ``mi`` cancels identically because
    ``X(v|N) L(j,v) (mi(j) + Qtot(j|N-e_v)) = Q(j,v|N)`` is the MVA equation for
    any ``mi``. The equivalent statement for an infinite-server station,
    equation (2.9b), is recovered automatically because LINE folds the delay
    into the think time ``Z``, which enters only through ``X(v|N)`` and carries
    no queue-length moment of its own.

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

    Returns
    -------
    PfqnSensMva

    Notes
    -----
    Restricted to closed populations. Mixed and load-dependent models are
    handled by :func:`pfqn_sens_mvaldmx`.

    For a station of multiplicity ``mi[i] > 1``, which LINE treats as ``mi[i]``
    identical replicas sharing the demand row ``L[i, :]``, the moments returned
    are those of the aggregate queue length over the replicas.
    """
    L = np.atleast_2d(np.asarray(L, dtype=np.float64))
    Nf = np.asarray(N, dtype=np.float64).flatten()
    if np.any(np.isinf(Nf)):
        raise ValueError("pfqn_sens_mva requires a closed population; "
                         "use pfqn_sens_mvaldmx for mixed models")
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

    X = np.zeros((1, R))
    Q = np.zeros((M, R))
    U = np.zeros((M, R))
    C = np.zeros((M, R))
    QCov = np.zeros((M, R, R))

    if not np.any(N > 0):
        return _pack(X, Q, U, C, QCov)

    # see _kb/03-api-layer.md for rationale
    prods = np.zeros(max(R - 1, 0))
    for w in range(R - 1):
        prods[w] = np.prod(np.ones(R - w - 1) + N[w + 1:])
    first_non_empty = R - 1
    while N[first_non_empty] == 0:
        first_non_empty -= 1
    totpop = int(np.prod(N + 1))

    Qtot = np.zeros((totpop, M))          # Qtot[m, i]    = sum_r Q(i,r) at pop m
    Qcls = np.zeros((totpop, M, R))       # Qcls[m, i, r] = Q(i,r) at pop m
    Xall = np.zeros((totpop, R))          # Xall[m, r]    = X(r) at pop m
    Ssum = np.zeros((totpop, M, R))       # Ssum[m, j, k] = sum_t Cov[n(j,k),n(j,t)]

    n = np.zeros(R, dtype=int)
    n[first_non_empty] = 1
    rows = np.zeros(R, dtype=int)         # rows[s] = lattice index of n - e_s
    currentpop = 1
    ctr = totpop

    while ctr > 0:
        # ---- mean value analysis step at population n --------------------
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
            rows[s] = pos
            CNtot = 0.0
            for i in range(M):
                C[i, s] = L[i, s] * (mi[i] + Qtot[pos, i])
                CNtot += C[i, s]
            den = Z[s] + CNtot
            X[0, s] = n[s] / den if den > 0 else 0.0
            Xall[currentpop, s] = X[0, s]
            for i in range(M):
                Q[i, s] = X[0, s] * C[i, s]
                Qcls[currentpop, i, s] = Q[i, s]
                Qtot[currentpop, i] += Q[i, s]

        # see _kb/03-api-layer.md for rationale
        for j in range(M):
            for k in range(R):
                Qjk = Qcls[currentpop, j, k]
                sk = 0.0
                for t in range(R):
                    Qjt = Qcls[currentpop, j, t]
                    wkt = Qjt * (Qcls[rows[t], j, k] - Qjk)
                    if k == t:
                        wkt += Qjt
                    wkt += Xall[currentpop, t] * L[j, t] * Ssum[rows[t], j, k]
                    QCov[j, k, t] = wkt
                    sk += wkt
                Ssum[currentpop, j, k] = sk

        # ---- odometer advance --------------------------------------------
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

    return _pack(X, Q, U, C, QCov)


__all__ = [
    'pfqn_sens_mva',
    'PfqnSensMva',
]
