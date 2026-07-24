"""
Exact queue-length variances and covariances for mixed open/closed product-form
queueing networks with limited load dependence, together with the effective
capacity terms of the mixed load-dependent MVA and their exact derivatives with
respect to the open-class load.

Native Python implementation (no JPype / JVM dependency). Mirrors the MATLAB
references ``pfqn_sens_ldmx_ec.m`` and ``pfqn_sens_mvaldmx.m``.

References:
    I. F. Akyildiz and J. C. Strelen, "Moment Analysis for Load-Dependent Mixed
    Product Form Queueing Networks", IEEE Trans. Communications 39(6):828-832,
    1991.
    The closed load-independent case reduces to E. de Souza e Silva and
    R. R. Muntz, IEEE Trans. Computers 37(9):1125-1129, 1988, which
    :func:`pfqn_sens_mva` implements directly.
"""

from typing import Optional, Tuple
import numpy as np

from .schmidt import pprod, hashpop
from .utils import oner

_EPS = np.finfo(float).eps


def pfqn_sens_ldmx_ec(lam: np.ndarray, D: np.ndarray, mu: np.ndarray
                     ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray,
                                np.ndarray, np.ndarray, np.ndarray]:
    """Effective capacity terms of the mixed load-dependent MVA and their exact
    derivatives with respect to the open-class load.

    Computes the terms ``EC``, ``E`` and ``Eprime`` of the mixed load-dependent
    MVA of Bruell-Balbo-Afshari, exactly as :func:`pfqn_ldmx_ec` does, and
    additionally their exact analytic derivatives with respect to the open-class
    load ``Lo[i]`` of each station.

    ``Lo[i] = sum_r lambda[r] D[i, r]`` is the only channel through which a
    service demand enters ``E``, ``Eprime`` and ``EC``: the load-dependent rates
    ``mu`` are independent of the demands. Station ``i``'s terms depend on
    ``Lo[i]`` alone, so a single derivative per station is enough, and the chain
    rule then yields the derivative with respect to any demand-scaling
    parameter. This is the factorization behind equations (19), (21) and
    (24)-(31) of Akyildiz and Strelen, which are reproduced here term by term.

    Parameters
    ----------
    lam : (R,) array
        Arrival rate vector. Zero for closed classes.
    D : (M, R) array
        Service demand matrix.
    mu : (M, Nt) array
        Load-dependent rate matrix, limited load dependence.

    Returns
    -------
    EC : (M, Nt) array
        Effective capacity matrix.
    E : (M, 1 + Nt) array
        E-function values, ``E[i, n]`` holding the MATLAB ``E(i, 1+n)``.
    Eprime : (M, 1 + Nt) array
        E-prime function values, indexed as ``E``.
    Lo : (M,) array
        Open class load vector.
    dEC : (M, Nt) array
        ``dEC[i, n] = dEC(i, n)/dLo(i)``.
    dE : (M, 1 + Nt) array
        ``dE[i, n] = dE(i, n)/dLo(i)``.
    dEprime : (M, 1 + Nt) array
        ``dEprime[i, n] = dEprime(i, n)/dLo(i)``.
    """
    mu = np.atleast_2d(np.asarray(mu, dtype=float))
    D = np.atleast_2d(np.asarray(D, dtype=float))
    lam = np.asarray(lam, dtype=float).flatten()

    M, Nt = mu.shape

    Lo = np.zeros(M)
    for ist in range(M):
        Lo[ist] = np.dot(lam, D[ist, :])

    # limited load dependence level: first index at which mu attains its final
    # (saturated) value, stored 1-based to mirror the MATLAB recursion indices
    b = np.zeros(M, dtype=int)
    for ist in range(M):
        b[ist] = Nt
        for k in range(Nt):
            if mu[ist, k] == mu[ist, -1]:
                b[ist] = k + 1
                break
    max_b = int(np.max(b))

    # extend mu past Nt with its saturated value, as the recursion reads C at
    # indices up to n + n0 + 1
    mu_ext = np.zeros((M, Nt + 1 + max_b))
    mu_ext[:, :Nt] = mu
    mu_ext[:, Nt:] = mu[:, -1].reshape(-1, 1)
    C = 1.0 / mu_ext

    EC = np.zeros((M, Nt))
    E = np.zeros((M, 1 + Nt))
    Eprime = np.zeros((M, 1 + Nt))
    dEC = np.zeros((M, Nt))
    dE = np.zeros((M, 1 + Nt))
    dEprime = np.zeros((M, 1 + Nt))

    for ist in range(M):
        bi = int(b[ist])
        Cb = C[ist, bi - 1]
        Lo_i = Lo[ist]
        den = 1.0 - Lo_i * Cb   # geometric tail factor of the limited load dep.

        E1 = np.zeros(1 + Nt)
        dE1 = np.zeros(1 + Nt)
        E2 = np.zeros(1 + Nt)
        dE2 = np.zeros(1 + Nt)
        E3 = np.zeros(1 + Nt)
        dE3 = np.zeros(1 + Nt)
        E2prime = np.zeros(1 + Nt)
        dE2prime = np.zeros(1 + Nt)
        ncol = 1 + max(0, bi - 2)
        F2 = np.zeros((1 + Nt, ncol))
        dF2 = np.zeros((1 + Nt, ncol))
        F3 = np.zeros((1 + Nt, ncol))
        dF3 = np.zeros((1 + Nt, ncol))
        F2prime = np.zeros((1 + Nt, ncol))
        dF2prime = np.zeros((1 + Nt, ncol))

        for n in range(Nt + 1):
            if n >= bi:
                # E(n) = 1/den^(n+1)  =>  dE/dLo = (n+1)*Cb/den^(n+2)
                E[ist, n] = 1.0 / den ** (n + 1)
                dE[ist, n] = (n + 1) * Cb / den ** (n + 2)
                Eprime[ist, n] = Cb * E[ist, n]
                dEprime[ist, n] = Cb * dE[ist, n]
            else:  # n <= bi-1
                # ---- E1 and its derivative, eq. (25)-(26) -----------------
                if n == 0:
                    E1[n] = 1.0 / den
                    dE1[n] = Cb / den ** 2
                    for j in range(bi - 1):
                        E1[n] *= C[ist, j] / Cb
                        dE1[n] *= C[ist, j] / Cb
                else:
                    fac = Cb / C[ist, n - 1]
                    E1[n] = (1.0 / den) * fac * E1[n - 1]
                    dE1[n] = (Cb / den ** 2) * fac * E1[n - 1] \
                        + (1.0 / den) * fac * dE1[n - 1]

                # ---- F2 and its derivative, eq. (27)-(28) -----------------
                for n0 in range(bi - 1):
                    if n0 == 0:
                        F2[n, n0] = 1.0
                        dF2[n, n0] = 0.0
                    else:
                        coef = (n + n0) / n0 * C[ist, n + n0 - 1]
                        F2[n, n0] = coef * Lo_i * F2[n, n0 - 1]
                        dF2[n, n0] = coef * (F2[n, n0 - 1] + Lo_i * dF2[n, n0 - 1])
                E2[n] = np.sum(F2[n, :max(0, bi - 1)])
                dE2[n] = np.sum(dF2[n, :max(0, bi - 1)])

                # ---- F3 and its derivative, eq. (29)-(30) -----------------
                for n0 in range(bi - 1):
                    if n == 0 and n0 == 0:
                        F3[n, n0] = 1.0
                        for j in range(bi - 1):
                            F3[n, n0] *= C[ist, j] / Cb
                        dF3[n, n0] = 0.0
                    elif n > 0 and n0 == 0:
                        fac = Cb / C[ist, n - 1]
                        F3[n, n0] = fac * F3[n - 1, 0]
                        dF3[n, n0] = fac * dF3[n - 1, 0]
                    else:
                        coef = (n + n0) / n0 * Cb
                        F3[n, n0] = coef * Lo_i * F3[n, n0 - 1]
                        dF3[n, n0] = coef * (F3[n, n0 - 1] + Lo_i * dF3[n, n0 - 1])
                E3[n] = np.sum(F3[n, :max(0, bi - 1)])
                dE3[n] = np.sum(dF3[n, :max(0, bi - 1)])

                # ---- F2prime and its derivative ---------------------------
                for n0 in range(bi - 1):
                    if n0 == 0:
                        F2prime[n, n0] = C[ist, n]
                        dF2prime[n, n0] = 0.0
                    else:
                        coef = (n + n0) / n0 * C[ist, n + n0]
                        F2prime[n, n0] = coef * Lo_i * F2prime[n, n0 - 1]
                        dF2prime[n, n0] = coef * (F2prime[n, n0 - 1]
                                                  + Lo_i * dF2prime[n, n0 - 1])
                E2prime[n] = np.sum(F2prime[n, :max(0, bi - 1)])
                dE2prime[n] = np.sum(dF2prime[n, :max(0, bi - 1)])

                # E = E1 + E2 - E3, eq. (23)-(24)
                E[ist, n] = E1[n] + E2[n] - E3[n]
                dE[ist, n] = dE1[n] + dE2[n] - dE3[n]
                if n < bi - 1:
                    Eprime[ist, n] = Cb * E1[n] + E2prime[n] - Cb * E3[n]
                    dEprime[ist, n] = Cb * dE1[n] + dE2prime[n] - Cb * dE3[n]
                else:  # n >= bi-1
                    Eprime[ist, n] = Cb * E[ist, n]
                    dEprime[ist, n] = Cb * dE[ist, n]

        # EC(n) = C(n)*E(n)/E(n-1), eq. (19); quotient rule for the derivative
        for n in range(1, Nt + 1):
            EC[ist, n - 1] = C[ist, n - 1] * E[ist, n] / E[ist, n - 1]
            dEC[ist, n - 1] = C[ist, n - 1] * (dE[ist, n] * E[ist, n - 1]
                                               - E[ist, n] * dE[ist, n - 1]) \
                / E[ist, n - 1] ** 2

    return EC, E, Eprime, Lo, dEC, dE, dEprime


class PfqnSensMvaldmx:
    """Result container for :func:`pfqn_sens_mvaldmx`.

    Attributes
    ----------
    X, Q, U, R : np.ndarray
        Base measures, identical to
        ``pfqn_mvaldmx(lambda, D, N, Z, mu, S)``. ``X`` is (R,) throughput per
        class, ``Q`` (M x R) mean queue length, ``U`` (M x R) utilization,
        ``R`` (M x R) residence time.
    QCov : np.ndarray (M x R x R)
        ``QCov[i, r, s] = Cov[n(i,r), n(i,s)]``, the same-station block.
    QCovFull : np.ndarray (M x R x M x R)
        ``QCovFull[i, r, j, s] = Cov[n(i,r), n(j,s)]``.
    QVar : np.ndarray (M x R)
        ``QVar[i, r] = Var[n(i,r)]``.
    QTotVar : np.ndarray (M,)
        ``QTotVar[i] = Var[sum_r n(i,r)]``.
    QCovAsym : float
        ``max |QCovFull[i, r, j, s] - QCovFull[j, s, i, r]|`` before
        symmetrization. The two entries are produced by differentiating two
        different classes' equations, so this residual is an independent check
        of the recursion and should sit at roundoff.
    """

    def __init__(self, X, Q, U, R, QCov, QCovFull, QVar, QTotVar, QCovAsym):
        self.X = X
        self.Q = Q
        self.U = U
        self.R = R
        self.QCov = QCov
        self.QCovFull = QCovFull
        self.QVar = QVar
        self.QTotVar = QTotVar
        self.QCovAsym = QCovAsym


def pfqn_sens_mvaldmx(lam: np.ndarray, D: np.ndarray, N: np.ndarray,
                     Z: np.ndarray, mu: Optional[np.ndarray] = None,
                     S: Optional[np.ndarray] = None) -> PfqnSensMvaldmx:
    """Exact second moments of the queue lengths of a mixed open/closed
    product-form network with limited load-dependent service rates.

    This is the load-dependent and mixed counterpart of :func:`pfqn_sens_mva`,
    which is restricted to closed load-independent models.

    The method is the moment analysis of Akyildiz and Strelen. Their Theorem 1,
    equation (11), states that multiplying a queue-length moment by one further
    factor ``Q_jT`` costs one derivative with respect to a parameter ``y_j``
    that scales the service demands ``s_ir`` of the classes ``r`` in ``T`` at
    station ``j``:

        E[Q_jT^k ...] = d/dy_j E[Q_jT^(k-1) ...]|_{y_j=1}
                        + nbar_jT E[Q_jT^(k-1) ...]

    Taking ``k=2`` and ``T={s}`` gives the second moment, hence

        Cov[n(i,r), n(j,s)] = d nbar(i,r) / dy_(j,s) |_{y=1}

    which is evaluated here by forward-mode differentiation of the mixed
    load-dependent MVA of Bruell-Balbo-Afshari, i.e. of exactly the recursion
    implemented by :func:`pfqn_mvaldmx`. The differentiated equations are (13)
    for the residence times, (15)-(17) for the conditional marginal
    probabilities, (18) for the throughputs, (19) and (24)-(31) for the
    effective capacities (delegated to :func:`pfqn_sens_ldmx_ec`), (32) for the
    closed-class queue lengths and (33) for the open-class ones.

    Because a demand-scaling parameter perturbs the whole network through the
    closed-class throughputs, the derivatives must be propagated for every
    parameter, so the cross-station covariances come out at no extra cost and
    are returned in ``QCovFull``. This is unlike :func:`pfqn_sens_mva`, whose
    cheaper same-station recursion cannot reach them.

    Parameters
    ----------
    lam : (R,) array
        Arrival rate vector. Must be zero on closed classes.
    D : (M, R) array
        Service demand matrix.
    N : (R,) array
        Population vector. ``inf`` entries denote open classes.
    Z : (R,) array
        Think time vector.
    mu : (M, sum(N)) array, optional
        Load-dependent rate matrix, limited load dependence (default ones).
    S : (M,) array, optional
        Number of servers per station. Accepted for signature compatibility with
        :func:`pfqn_mvaldmx`, which likewise does not read it: the multiserver
        behaviour is carried entirely by the rates ``mu``.

    Returns
    -------
    PfqnSensMvaldmx

    Notes
    -----
    Open classes are supported: an open class contributes to the load ``Lo[i]``
    that drives the effective capacities, and equation (21) supplies the
    corresponding ``dLo/dy``.

    The moments of an open class are those of its queue length at a station,
    which is finite even though its population is infinite.
    """
    D = np.atleast_2d(np.asarray(D, dtype=float))
    N = np.asarray(N, dtype=float).flatten()
    Z = np.asarray(Z, dtype=float).flatten()
    lam = np.asarray(lam, dtype=float).flatten()
    M, R = D.shape

    NCtot_all = int(np.sum(N[np.isfinite(N)]))
    if mu is None:
        mu = np.ones((M, NCtot_all))
        S = np.ones(M)
    if S is None:
        S = np.ones(M)
    mu = np.atleast_2d(np.asarray(mu, dtype=float))
    if mu.shape[1] < NCtot_all:
        raise ValueError("pfqn_sens_mvaldmx requires to specify the load-dependent "
                         "rates with one job more than the maximum closed "
                         "population.")
    nz = np.where(lam != 0)[0]
    if np.any((N[nz] > 0) & np.isfinite(N[nz])):
        raise ValueError("Arrival rate cannot be specified on closed classes.")

    openClasses = np.where(np.isinf(N))[0]
    closedClasses = np.array([r for r in range(R) if r not in openClasses],
                             dtype=int)
    C = len(closedClasses)
    if C == 0:
        raise ValueError("pfqn_sens_mvaldmx requires at least one closed class; "
                         "use the open-class formulas directly otherwise.")

    # up to sum(N)+1, limited load dependence
    mu = np.hstack([mu, mu[:, -1:]])
    EC, E, Eprime, _, dEC_dLo = pfqn_sens_ldmx_ec(lam, D, mu)[:5]

    Dc = D[:, closedClasses]
    Nc = N[closedClasses].astype(int)
    Zc = Z[closedClasses]
    NCtot = int(np.sum(Nc))

    # ---- parameter list: y(j,r) multiplies the demand D(j,r) -------------
    P = M * R
    pidx = np.zeros((M, R), dtype=int)
    pj = np.zeros(P, dtype=int)
    pr = np.zeros(P, dtype=int)
    p = 0
    for j in range(M):
        for r in range(R):
            pidx[j, r] = p
            pj[p] = j
            pr[p] = r
            p += 1
    # see _kb/03-api-layer.md for rationale
    dLo = np.zeros((M, P))
    for q in range(P):
        if np.isinf(N[pr[q]]):
            dLo[pj[q], q] = lam[pr[q]] * D[pj[q], pr[q]]

    # ---- population recursion -------------------------------------------
    prods = np.zeros(C, dtype=int)
    for r in range(C):
        prods[r] = int(np.prod(Nc[:r] + 1))
    NT = int(np.prod(1 + Nc))
    Pc = np.zeros((M, 1 + NCtot, NT))
    dPc = np.zeros((M, 1 + NCtot, NT, P))
    x = np.zeros((C, NT))
    dx = np.zeros((C, NT, P))
    w = np.zeros((M, C, NT))
    dw = np.zeros((M, C, NT, P))

    nvec = pprod(Nc)
    h0 = hashpop(nvec, Nc, C, prods) - 1
    for ist in range(M):
        Pc[ist, 0, h0] = 1.0                                  # eq. (16)

    while np.all(nvec >= 0):
        hnvec = hashpop(nvec, Nc, C, prods) - 1
        nc = int(np.sum(nvec))

        # ---- residence times, eq. (12) and its derivative eq. (13) -------
        for ist in range(M):
            for c in range(C):
                if nvec[c] > 0:
                    hnvec_c = hashpop(oner(nvec, c), Nc, C, prods) - 1
                    cls = closedClasses[c]
                    acc = 0.0
                    dacc = np.zeros(P)
                    for n in range(1, nc + 1):
                        Pprev = Pc[ist, n - 1, hnvec_c]
                        acc += n * EC[ist, n - 1] * Pprev
                        dacc += n * (dEC_dLo[ist, n - 1] * dLo[ist, :] * Pprev
                                     + EC[ist, n - 1] * dPc[ist, n - 1, hnvec_c, :])
                    w[ist, c, hnvec] = Dc[ist, c] * acc
                    dwq = Dc[ist, c] * dacc
                    # d(D*y)/dy = D
                    dwq[pidx[ist, cls]] += Dc[ist, c] * acc
                    dw[ist, c, hnvec, :] = dwq

        # ---- throughputs, eq. (18) --------------------------------------
        for c in range(C):
            den = Zc[c] + np.sum(w[:, c, hnvec])
            x[c, hnvec] = nvec[c] / den if den > 0 else 0.0
            if nvec[c] > 0 and den > 0:
                sdw = np.sum(dw[:, c, hnvec, :], axis=0)
                dx[c, hnvec, :] = -nvec[c] / den ** 2 * sdw

        # ---- conditional marginal probabilities, eq. (14)-(15) -----------
        for ist in range(M):
            for n in range(1, nc + 1):
                for c in range(C):
                    if nvec[c] > 0:
                        hnvec_c = hashpop(oner(nvec, c), Nc, C, prods) - 1
                        cls = closedClasses[c]
                        Pprev = Pc[ist, n - 1, hnvec_c]
                        Pc[ist, n, hnvec] += Dc[ist, c] * EC[ist, n - 1] \
                            * x[c, hnvec] * Pprev
                        dt = Dc[ist, c] * (
                            dEC_dLo[ist, n - 1] * dLo[ist, :] * x[c, hnvec] * Pprev
                            + EC[ist, n - 1] * dx[c, hnvec, :] * Pprev
                            + EC[ist, n - 1] * x[c, hnvec] * dPc[ist, n - 1, hnvec_c, :])
                        dt[pidx[ist, cls]] += Dc[ist, c] * EC[ist, n - 1] \
                            * x[c, hnvec] * Pprev
                        dPc[ist, n, hnvec, :] += dt
            # see _kb/03-api-layer.md for rationale
            Pc[ist, 0, hnvec] = max(_EPS, 1.0 - np.sum(Pc[ist, 1:nc + 1, hnvec]))
            dPc[ist, 0, hnvec, :] = -np.sum(dPc[ist, 1:nc + 1, hnvec, :], axis=0)

        nvec = pprod(nvec, Nc)

    # ---- measures and their derivatives at the full population -----------
    hnvec = hashpop(Nc, Nc, C, prods) - 1
    XN = np.zeros(R)
    QN = np.zeros((M, R))
    UN = np.zeros((M, R))
    CN = np.zeros((M, R))
    dQN = np.zeros((M, R, P))

    # closed classes, eq. (32)
    for c in range(C):
        cls = closedClasses[c]
        XN[cls] = x[c, hnvec]
        hnvec_c = hashpop(oner(Nc, c), Nc, C, prods) - 1
        for ist in range(M):
            CN[ist, cls] = w[ist, c, hnvec]
            QN[ist, cls] = XN[cls] * CN[ist, cls]
            dQN[ist, cls, :] = dx[c, hnvec, :] * w[ist, c, hnvec] \
                + x[c, hnvec] * dw[ist, c, hnvec, :]
            uacc = 0.0
            for n in range(1, NCtot + 1):
                uacc += Dc[ist, c] * x[c, hnvec] * Eprime[ist, n - 1] \
                    / E[ist, n - 1] * Pc[ist, n - 1, hnvec_c]
            UN[ist, cls] = uacc

    # open classes, eq. (33)
    for r in openClasses:
        XN[r] = lam[r]
        for ist in range(M):
            acc = 0.0
            dacc = np.zeros(P)
            for n in range(NCtot + 1):
                Pn = Pc[ist, n, hnvec]
                acc += (n + 1) * EC[ist, n] * Pn
                dacc += (n + 1) * (dEC_dLo[ist, n] * dLo[ist, :] * Pn
                                   + EC[ist, n] * dPc[ist, n, hnvec, :])
            QN[ist, r] = lam[r] * D[ist, r] * acc
            CN[ist, r] = QN[ist, r] / lam[r] if lam[r] > 0 else 0.0
            dq = lam[r] * D[ist, r] * dacc
            dq[pidx[ist, r]] += lam[r] * D[ist, r] * acc
            dQN[ist, r, :] = dq
            uacc = 0.0
            for n in range(NCtot + 1):
                uacc += lam[r] * Eprime[ist, n + 1] / E[ist, n + 1] \
                    * Pc[ist, n, hnvec]
            UN[ist, r] = uacc

    # ---- moments ---------------------------------------------------------
    # Cov[n(i,r),n(j,s)] = d nbar(i,r) / dy_(j,s)
    QCovRaw = np.zeros((M, R, M, R))
    for i in range(M):
        for r in range(R):
            for j in range(M):
                for s in range(R):
                    QCovRaw[i, r, j, s] = dQN[i, r, pidx[j, s]]
    QCovFull = (QCovRaw + np.transpose(QCovRaw, (2, 3, 0, 1))) / 2.0
    if QCovRaw.size == 0:
        QCovAsym = 0.0
    else:
        QCovAsym = float(np.max(np.abs(QCovRaw - np.transpose(QCovRaw, (2, 3, 0, 1)))))

    QCov = np.zeros((M, R, R))
    QVar = np.zeros((M, R))
    QTotVar = np.zeros(M)
    for i in range(M):
        for r in range(R):
            for s in range(R):
                QCov[i, r, s] = QCovFull[i, r, i, s]
            QVar[i, r] = QCov[i, r, r]
        QTotVar[i] = np.sum(QCov[i, :, :])

    return PfqnSensMvaldmx(XN, QN, UN, CN, QCov, QCovFull, QVar, QTotVar,
                           QCovAsym)


__all__ = [
    'pfqn_sens_ldmx_ec',
    'pfqn_sens_mvaldmx',
    'PfqnSensMvaldmx',
]
