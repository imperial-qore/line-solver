"""
Exact analytic performance sensitivities for closed product-form queueing
networks. Dispatches to a CoMoM-backed kernel for the single-station repairman
model and to differentiated Mean Value Analysis otherwise.

Native Python implementation (no JPype / JVM dependency). Mirrors the MATLAB
reference ``pfqn_sens.m`` and the JAR ``Pfqn_sens``.

References:
    Z. Liu and P. Nain, "Sensitivity Results in Open, Closed and Mixed
    Product-Form Queueing Networks", INRIA RR-1144, 1989.
    X.-R. Cao and D.-J. Ma, "Performance sensitivity formulae, algorithms and
    estimates for closed queueing networks with exponential servers",
    Performance Evaluation 26:181-199, 1996.
    G. Casale, "CoMoM: Efficient Class-Oriented Evaluation of Multiclass
    Performance Models", IEEE TSE 2011.
"""

from typing import Dict, List
import numpy as np

_FINE_TOL = 1e-8


class PfqnSens:
    """Result container for :func:`pfqn_sens`.

    Attributes
    ----------
    X, Q, U, R : np.ndarray
        Base MVA measures. ``X`` is (1 x R) system throughput per class, ``Q``
        (M x R) mean queue length, ``U`` (M x R) utilization, ``R`` (M x R)
        residence time per visit.
    params : list of dict
        One entry per differentiation parameter ``p`` with keys ``type``
        ('L' or 'Z'), ``station`` (i, or -1 for Z) and ``jobclass`` (r).
    dX : np.ndarray (R x P)
    dQ, dU, dR : np.ndarray (M x R x P)
        Derivative of each base measure w.r.t. parameter ``p``.
    QCov : np.ndarray (M x R x M x R)
        Queue-length covariance, ``QCov[i, r, j, s] = Cov[n_ir, n_js]
        = D_js dQ_ir/dD_js``. The same-station blocks (i == j) come from
        :func:`pfqn_sens_mva`; the cross-station ones are read off the Jacobian.
    QVar : np.ndarray (M x R)
        Queue-length variance, ``QVar[i, r] = QCov[i, r, i, r]``.
    QTotVar : np.ndarray (M,)
        Variance of the total queue length per station, ``QTotVar[i] =
        Var[sum_r n_ir]``.
    QCovAsym : float
        Roundoff-level residual of the moment recursion, see
        :func:`pfqn_sens_mva`.
    """

    def __init__(self, X, Q, U, R, params, dX, dQ, dU, dR, QCov=None, QVar=None,
                 QTotVar=None, QCovAsym=None):
        self.X = X
        self.Q = Q
        self.U = U
        self.R = R
        self.params = params
        self.dX = dX
        self.dQ = dQ
        self.dU = dU
        self.dR = dR
        self.QCov = QCov
        self.QVar = QVar
        self.QTotVar = QTotVar
        self.QCovAsym = QCovAsym


def pfqn_sens(L: np.ndarray, N: np.ndarray, Z: np.ndarray = None,
              mi: np.ndarray = None) -> PfqnSens:
    """Exact derivatives of {X, Q, U, R} w.r.t. demands L and think times Z.

    Derivatives are analytic (exact to machine precision), not finite
    differences. The CoMoM-backed kernel is selected for the repairman model
    (a single single-server queue, M=1, plus a think-time delay where every
    populated class has positive think time and demand); it is polynomial in
    the number of classes R. Every other model uses forward-mode
    differentiation of the exact MVA recursion.

    Parameters
    ----------
    L : (M, R) array
        Service demand matrix, ``L[i, r] = visits_ir / rate_ir``.
    N : (R,) array
        Population per class.
    Z : (R,) array, optional
        Think time per class (default zeros).
    mi : (M,) array, optional
        Station residence multiplicity (default ones).

    Returns
    -------
    PfqnSens
    """
    L = np.asarray(L, dtype=np.float64)
    N = np.ceil(np.asarray(N, dtype=np.float64).flatten()).astype(int)
    R = len(N)
    if L.ndim == 1:
        L = L.reshape(-1, 1) if R == 1 else L.reshape(1, -1)
    M = L.shape[0]
    if L.shape[1] != R:
        raise ValueError("demand matrix columns must match population size")
    if Z is None:
        Z = np.zeros(R)
    else:
        Z = np.asarray(Z, dtype=np.float64).flatten()
    if mi is None:
        mi = np.ones(M)
    else:
        mi = np.asarray(mi, dtype=np.float64).flatten()

    # see _kb/03-api-layer.md for rationale
    populated = N > 0
    use_comom = (M == 1 and np.all(mi == 1) and np.any(populated)
                 and np.all(Z[populated] > _FINE_TOL)
                 and np.all(L[0, populated] > _FINE_TOL))
    if use_comom:
        sens = _sens_comom(L, N, Z)
    else:
        sens = _sens_mva(L, N, Z, mi)
    _attach_second_moments(sens, L, N, Z, mi)
    return sens


# =========================================================================
# CoMoM-backed kernel (M=1 repairman model)
# =========================================================================
class _ReplNC:
    """Replicated-model normalizing-constant moments via memoized CoMoM."""

    def __init__(self, D, Z):
        from .comom import pfqn_comomrm
        self._comomrm = pfqn_comomrm
        self.D = D
        self.Z = Z
        self.R = len(D)
        self._cache = {}

    def lgm(self, m, n):
        """log normalizing constant of the m-replica model at population n."""
        n = np.asarray(n, dtype=int)
        nz = n > 0
        if not np.any(nz):
            return 0.0
        key = (m, tuple(int(x) for x in n))
        cached = self._cache.get(key)
        if cached is not None:
            return cached
        # strip zero-population classes (they leave the NC unchanged)
        Ls = self.D[nz].reshape(1, -1)
        Ns = n[nz].astype(float)
        Zs = self.Z[nz]
        lg = self._comomrm(Ls, Ns, Zs, m, _FINE_TOL).lG
        self._cache[key] = lg
        return lg

    def qmean(self, n, s):
        """mean class-s queue at population n (single station, m=1 model)."""
        if n[s] < 1:
            return 0.0
        nm = np.array(n, dtype=int)
        nm[s] -= 1
        return self.D[s] * np.exp(self.lgm(2, nm) - self.lgm(1, n))

    def xput(self, m, n, s):
        """class-s throughput G_m(n-1_s)/G_m(n) in the m-replica model."""
        if n[s] < 1:
            return 0.0
        nm = np.array(n, dtype=int)
        nm[s] -= 1
        return np.exp(self.lgm(m, nm) - self.lgm(m, n))

    def qplus(self, n, s):
        """Q^{+1}_{1,s}(n): class-s queue at one replica of the doubled station."""
        if n[s] < 1:
            return 0.0
        nm = np.array(n, dtype=int)
        nm[s] -= 1
        return self.D[s] * np.exp(self.lgm(3, nm) - self.lgm(2, n))


def _sens_comom(L: np.ndarray, N: np.ndarray, Z: np.ndarray) -> PfqnSens:
    R = len(N)
    D = L[0, :].astype(float)
    Nn = np.asarray(N, dtype=int)

    # parameter list: L(0,r) first, then Z(r) (mirrors _sens_mva ordering)
    P = 2 * R
    params: List[Dict] = [None] * P
    pL = np.zeros(R, dtype=int)
    pZ = np.zeros(R, dtype=int)
    for r in range(R):
        pL[r] = r
        params[pL[r]] = {'type': 'L', 'station': 0, 'jobclass': r}
        pZ[r] = R + r
        params[pZ[r]] = {'type': 'Z', 'station': -1, 'jobclass': r}

    X = np.zeros((1, R))
    Q = np.zeros((1, R))
    U = np.zeros((1, R))
    C = np.zeros((1, R))
    dX = np.zeros((R, P))
    dQ = np.zeros((1, R, P))
    dU = np.zeros((1, R, P))
    dC = np.zeros((1, R, P))

    if not np.any(Nn > 0):
        return PfqnSens(X, Q, U, C, params, dX, dQ, dU, dC)

    nc = _ReplNC(D, np.asarray(Z, dtype=float))

    # base measures
    for r in range(R):
        if Nn[r] >= 1:
            X[0, r] = nc.xput(1, Nn, r)
    for s in range(R):
        Q[0, s] = nc.qmean(Nn, s)
    for r in range(R):
        U[0, r] = X[0, r] * D[r]
        if X[0, r] > 0:
            C[0, r] = Q[0, r] / X[0, r]

    # Jacobian
    for r in range(R):
        if Nn[r] < 1:
            continue   # empty class: X=Q=0, all derivatives 0
        Nr = np.array(Nn, dtype=int)
        Nr[r] -= 1
        Xr = X[0, r]
        Qr = Q[0, r]
        for s in range(R):
            drs = 1.0 if r == s else 0.0

            # L(0,s) parameter
            p = pL[s]
            Vrs = Qr * (drs + 2.0 * nc.qplus(Nr, s) - Q[0, s])   # Cov[n_r, n_s]
            dQ_L = Vrs / D[s]
            dX_L = Xr * (nc.qmean(Nr, s) - Q[0, s]) / D[s]
            dU_L = dX_L * D[r] + Xr * drs
            dQ[0, r, p] = dQ_L
            dX[r, p] = dX_L
            dU[0, r, p] = dU_L
            if Xr > 0:
                dC[0, r, p] = (dQ_L * Xr - Qr * dX_L) / (Xr * Xr)

            # Z(s) parameter: d log G_m(n)/dZ_s = G_m(n-1_s)/G_m(n)
            p = pZ[s]
            dX_Z = Xr * (nc.xput(1, Nr, s) - X[0, s])
            dQ_Z = Qr * (nc.xput(2, Nr, s) - X[0, s])
            dU_Z = dX_Z * D[r]
            dQ[0, r, p] = dQ_Z
            dX[r, p] = dX_Z
            dU[0, r, p] = dU_Z
            if Xr > 0:
                dC[0, r, p] = (dQ_Z * Xr - Qr * dX_Z) / (Xr * Xr)

    return PfqnSens(X, Q, U, C, params, dX, dQ, dU, dC)


# =========================================================================
# Differentiated-MVA kernel (general model)
# =========================================================================
def _sens_mva(L: np.ndarray, N: np.ndarray, Z: np.ndarray,
              mi: np.ndarray) -> PfqnSens:
    M, R = L.shape

    # ---- parameter list: all L[i, r], then all Z[r] ---------------------
    P = M * R + R
    params: List[Dict] = []
    pL = np.zeros((M, R), dtype=int)
    pZ = np.zeros(R, dtype=int)
    p = 0
    for i in range(M):
        for r in range(R):
            params.append({'type': 'L', 'station': i, 'jobclass': r})
            pL[i, r] = p
            p += 1
    for r in range(R):
        params.append({'type': 'Z', 'station': -1, 'jobclass': r})
        pZ[r] = p
        p += 1

    X = np.zeros((1, R))
    Q = np.zeros((M, R))
    U = np.zeros((M, R))
    C = np.zeros((M, R))
    dX = np.zeros((R, P))
    dQ = np.zeros((M, R, P))
    dU = np.zeros((M, R, P))
    dC = np.zeros((M, R, P))

    if not np.any(N > 0):
        return PfqnSens(X, Q, U, C, params, dX, dQ, dU, dC)

    # ---- population-lattice odometer, identical to pfqn_mva -------------
    prods = np.zeros(R - 1)
    for w in range(R - 1):
        prods[w] = np.prod(np.ones(R - w - 1) + N[w + 1:])
    first_non_empty = R - 1
    while first_non_empty >= 0 and N[first_non_empty] == 0:
        first_non_empty -= 1
    totpop = int(np.prod(N + 1))
    Qtot = np.zeros((totpop, M))
    Qtotd = np.zeros((totpop, M, P))

    n = np.zeros(R, dtype=int)
    n[first_non_empty] = 1
    currentpop = 1
    ctr = totpop

    while ctr > 0:
        for s in range(R):
            if n[s] > 0:
                n[s] -= 1
                pos = int(n[R - 1])
                for w in range(R - 1):
                    pos += int(n[w] * prods[w])
                n[s] += 1
            else:
                pos = 0
            base = mi + Qtot[pos, :]                        # (M,)
            C[:, s] = L[:, s] * base
            Cd = L[:, s][:, None] * Qtotd[pos, :, :]        # (M, P)
            for i in range(M):
                Cd[i, pL[i, s]] += base[i]                  # dL[i, s] indicator
            CNtot = C[:, s].sum()
            CNtotd = Cd.sum(axis=0)                          # (P,)
            den = Z[s] + CNtot
            X[0, s] = n[s] / den if den > 0 else 0.0
            Xd = -n[s] * CNtotd / (den * den) if den > 0 else np.zeros(P)
            if den > 0:
                Xd[pZ[s]] -= n[s] / (den * den)             # dZ[s] indicator
            dX[s, :] = Xd
            for i in range(M):
                Q[i, s] = X[0, s] * C[i, s]
                Qd = Xd * C[i, s] + X[0, s] * Cd[i, :]
                dQ[i, s, :] = Qd
                dC[i, s, :] = Cd[i, :]
                Qtot[currentpop, i] += Q[i, s]
                Qtotd[currentpop, i, :] += Qd
        # advance odometer (identical to pfqn_mva)
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

    # ---- utilization and its derivatives -------------------------------
    for i in range(M):
        for r in range(R):
            U[i, r] = X[0, r] * L[i, r]
            Ud = dX[r, :] * L[i, r]
            Ud[pL[i, r]] += X[0, r]                          # dL[i, r] indicator
            dU[i, r, :] = Ud

    return PfqnSens(X, Q, U, C, params, dX, dQ, dU, dC)


# see _kb/03-api-layer.md for rationale
def _attach_second_moments(sens: PfqnSens, L: np.ndarray, N: np.ndarray,
                           Z: np.ndarray, mi: np.ndarray):
    from .sens_mva import pfqn_sens_mva

    M, R = L.shape
    mom = pfqn_sens_mva(L, N, Z, mi)
    pLidx = np.zeros((M, R), dtype=int)
    hasL = np.zeros((M, R), dtype=bool)
    for p, pr in enumerate(sens.params):
        if pr['type'] == 'L':
            pLidx[pr['station'], pr['jobclass']] = p
            hasL[pr['station'], pr['jobclass']] = True
    QCov = np.zeros((M, R, M, R))
    for i in range(M):
        for r in range(R):
            for j in range(M):
                for s in range(R):
                    if i == j:
                        QCov[i, r, j, s] = mom.QCov[i, r, s]
                    elif hasL[j, s]:
                        QCov[i, r, j, s] = L[j, s] * sens.dQ[i, r, pLidx[j, s]]
    sens.QCov = QCov
    sens.QVar = mom.QVar
    sens.QTotVar = mom.QTotVar
    sens.QCovAsym = mom.QCovAsym
