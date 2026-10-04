"""
Sojourn-time moments at the processor-sharing station of the closed
terminal-driven system of Mitra and Morrison (1983).

Native Python implementation (no JPype / JVM dependency). Mirrors the MATLAB
reference ``pfqn_respt_ps_moments.m``.

References:
    D. Mitra, J. A. Morrison, "Asymptotic Expansions of Moments of the Waiting
    Time in Closed and Open Processor-Sharing Systems with Multiple Job
    Classes", Adv. Appl. Prob. 15(4):813-839, 1983, Propositions 3 and 6.
"""

import numpy as np
from scipy.sparse import coo_matrix, eye as sparse_eye
from scipy.sparse.linalg import spsolve
from scipy.special import gammaln

AUTO_MAX = 4096        # state-space size below which auto goes exact
EXACT_MAX = 65536      # hard bound on an explicitly requested exact solve


class PfqnResptPsMoments:
    """Result container for :func:`pfqn_respt_ps_moments`.

    Attributes
    ----------
    method : list of str
        Per-class route taken: 'exact', 'asymptotic', 'unavailable' or 'none'.
    c0, c1 : np.ndarray (R,)
        Coefficients of the expansion ``E[W^2] ~ c0 + c1/expansionParam``, NaN
        on the exact route.
    alpha : np.ndarray (R,)
        Per-class unutilized fraction ``1 - sum_r lambda_r/q_r`` of the CPU in
        the corresponding open system.
    nstates : np.ndarray (R,)
        Size of the exact state space that the tagged class would need.
    expansionParam : float
        The large parameter ``Nexp = max_r Z(r)/S(r)``.
    """

    def __init__(self, R):
        self.method = ['none'] * R
        self.c0 = np.full(R, np.nan)
        self.c1 = np.full(R, np.nan)
        self.alpha = np.full(R, np.nan)
        self.nstates = np.full(R, np.nan)
        self.expansionParam = np.nan


def pfqn_respt_ps_moments(S, N, Z, method='auto'):
    """Sojourn-time moments at the PS station of a closed terminal-driven system.

    The system is a bank of terminals in series with a single processor-sharing
    CPU, with class-dependent exponential think times (mean ``Z[r]``) and
    class-dependent exponential service times (mean ``S[r]``), and ``N[r]`` jobs
    of class r cycling between the two.

    Two routes to the moments are implemented, both from Mitra and Morrison
    (1983):

    ``'exact'``
        solves the linear system ``c'[A - q_J I] = -pi'B`` of Proposition 3 on
        the state space ``{n : 0 <= n <= K}``, K being the population vector
        with the tagged class decremented by one. The moments are then
        ``E[W_J] = sum_n c(n)`` and ``(q_J/2) E[W_J^2] = sum_n (n'1+1) c(n)``.
        Exact to solver precision, at the cost of a linear solve of dimension
        ``prod_r (K[r]+1)``.

    ``'asymptotic'``
        evaluates the two leading terms of the asymptotic expansion in inverse
        powers of the large parameter ``Nexp = max_r Z[r]/S[r]``,
        ``E[W_J^2] ~ c0 + c1/Nexp``, of Proposition 6. The cost is a linear
        system of dimension R, the number of classes, and is therefore
        independent of the populations. Note that the expansion parameter is
        the think-to-service ratio and NOT the population, so a model with
        short think times is expanded in a small parameter no matter how many
        jobs it holds.

    ``'auto'``
        (default) takes the exact route when the state space has at most
        ``AUTO_MAX`` states and the asymptotic route otherwise.

    The asymptotic route requires the normal-usage condition ``alpha > 0``.
    Where it fails and the exact route is not affordable, the entry of W and W2
    is NaN and the result records 'unavailable'; asking for 'asymptotic'
    explicitly in that regime raises rather than returning a blank.

    Parameters
    ----------
    S : array_like (R,)
        Per-class mean service times at the PS station, positive.
    N : array_like (R,)
        Per-class populations, non-negative integers.
    Z : array_like (R,)
        Per-class mean think times, positive where ``N > 0``.
    method : str
        'auto' (default), 'exact' or 'asymptotic'.

    Returns
    -------
    W : np.ndarray (R,)
        Per-class mean sojourn times at the PS station.
    W2 : np.ndarray (R,)
        Per-class second moments of the sojourn time.
    out : PfqnResptPsMoments
        Route taken and expansion diagnostics.

    A class with ``N[r] = 0`` has no sojourn time and its entries are NaN.

    See also
    --------
    qsys_mm1_ps : the open counterpart, exact in closed form.
    """
    method = str(method).strip().lower()
    if method not in ('auto', 'exact', 'asymptotic'):
        raise ValueError("method must be one of auto, exact, asymptotic")

    S = np.asarray(S, dtype=float).flatten()
    N = np.asarray(N, dtype=float).flatten()
    Z = np.asarray(Z, dtype=float).flatten()
    R = S.size
    if N.size != R or Z.size != R:
        raise ValueError("S, N and Z must have the same number of classes")
    if not np.all(np.isfinite(S)) or np.any(S <= 0):
        raise ValueError("S must be finite and positive")
    if not np.all(np.isfinite(N)) or np.any(N < 0) or np.any(N != np.round(N)):
        raise ValueError("N must contain non-negative integers")
    act = np.flatnonzero(N > 0)
    if not np.all(np.isfinite(Z[act])) or np.any(Z[act] <= 0):
        raise ValueError("Z must be finite and positive for every populated class")

    W = np.full(R, np.nan)
    W2 = np.full(R, np.nan)
    out = PfqnResptPsMoments(R)
    if act.size == 0:
        return W, W2, out

    qa = 1.0 / S[act]
    pa = 1.0 / Z[act]
    out.expansionParam = float(np.max(qa / pa))

    for jj, J in enumerate(act):
        K = N[act].astype(int).copy()
        K[jj] -= 1
        ns = int(np.prod(K + 1))
        lam = pa * K
        alpha = 1.0 - float(np.sum(lam / qa))
        out.alpha[J] = alpha
        out.nstates[J] = ns
        use_exact = method == 'exact' or (method == 'auto' and ns <= AUTO_MAX)
        if use_exact:
            if ns > EXACT_MAX:
                raise ValueError(
                    "the exact route needs a linear solve of dimension %d, above "
                    "the bound of %d; use method = 'asymptotic'" % (ns, EXACT_MAX))
            W[J], W2[J] = _exact_moments(pa, qa, K, jj)
            out.method[J] = 'exact'
            continue
        if alpha <= 0:
            if method == 'asymptotic':
                raise ValueError(
                    "the asymptotic expansion needs normal usage alpha > 0, but "
                    "class %d gives alpha = %.6f" % (J, alpha))
            out.method[J] = 'unavailable'
            continue
        W[J], W2[J], out.c0[J], out.c1[J] = _asymptotic_moments(pa, qa, K, jj)
        out.method[J] = 'asymptotic'
    return W, W2, out


def _exact_moments(p, q, K, J):
    """Proposition 3: the moments follow from c, the solution of
    ``c'[A - q_J I] = -pi'B``, with A the generator-like operator of equation
    (26) and B the diagonal operator ``B(n,n) = n'1+1``."""
    R = K.size
    dims = K + 1
    ns = int(np.prod(dims))
    stride = np.concatenate(([1], np.cumprod(dims)[:-1])).astype(int)
    lin = np.arange(ns)
    states = np.zeros((ns, R), dtype=int)
    res = lin.copy()
    for j in range(R):
        states[:, j] = res % dims[j]
        res = res // dims[j]
    tot = states.sum(axis=1)

    # stationary law (15), in logs so that large populations do not overflow
    r = p / q
    logpi = gammaln(tot + 1.0)
    for j in range(R):
        nj = states[:, j]
        logpi = (logpi + gammaln(K[j] + 1.0) - gammaln(nj + 1.0)
                 - gammaln(K[j] - nj + 1.0))
        if r[j] > 0:
            logpi = logpi + nj * np.log(r[j])
        else:
            logpi = np.where(nj > 0, -np.inf, logpi)
    logpi = logpi - np.max(logpi)
    pin = np.exp(logpi)
    pin = pin / pin.sum()

    rows = []
    cols = []
    vals = []
    diagv = np.zeros(ns)
    for j in range(R):
        nj = states[:, j]
        dn = nj >= 1
        if np.any(dn):
            rows.append(lin[dn] - stride[j])
            cols.append(lin[dn])
            vals.append(p[j] * (K[j] - nj[dn] + 1) * tot[dn])
        up = nj <= K[j] - 1
        if np.any(up):
            rows.append(lin[up] + stride[j])
            cols.append(lin[up])
            vals.append((nj[up] + 1) * q[j])
        diagv = diagv - (p[j] * (K[j] - nj) * (tot + 1) + nj * q[j])
    rows.append(lin)
    cols.append(lin)
    vals.append(diagv)
    A = coo_matrix((np.concatenate(vals),
                    (np.concatenate(rows), np.concatenate(cols))),
                   shape=(ns, ns)).tocsc()

    M = (A - q[J] * sparse_eye(ns, format='csc')).T.tocsc()
    c = spsolve(M, -(tot + 1.0) * pin)
    W = float(np.sum(c))
    W2 = float(2.0 / q[J] * np.sum((tot + 1.0) * c))
    return W, W2


def _asymptotic_moments(p, q, K, J):
    """Proposition 6: the two leading terms of the expansion in 1/Nexp. Equation
    numbers below are those of Mitra and Morrison (1983)."""
    lam = p * K
    alpha = 1.0 - float(np.sum(lam / q))
    Nexp = float(np.max(q / p))                       # (50)
    Gam = Nexp * p / q                                # (51)
    beta = K / Nexp                                   # (51)
    qJ = q[J]

    den = 1.0 - float(np.sum(lam / (q + qJ)))
    F10 = (-1.0 / (alpha ** 2 * qJ)) * (
        1.0 - float(np.sum(lam * (q - qJ) / (q * (q + qJ))))) / den      # (110)
    c0 = -2.0 / qJ * F10

    bg2 = float(np.sum(beta * Gam ** 2))
    f1 = lam / (q + qJ) * (F10 - 2.0 / (alpha ** 2 * q))                 # (113iii)
    S2j = 6.0 / alpha ** 4 * (alpha * beta * Gam ** 2
                              + 2.0 * bg2 * beta * Gam)                  # (113i)
    S2js = 3.0 / alpha ** 3 * np.outer(beta * Gam, beta * Gam)           # (113ii)

    R = K.size
    Amat = np.eye(R)
    rhs = np.zeros(R)
    for j in range(R):
        for s in range(R):
            d = q[j] + q[s] + qJ
            Amat[j, s] -= lam[j] / d
            Amat[j, j] -= lam[s] / d
            rhs[j] += S2js[j, s] / d
        rhs[j] -= f1[j]
    F2 = np.linalg.solve(Amat, rhs)                                      # (112)

    f10 = -3.0 / (alpha ** 3 * qJ) * bg2                                 # (98)
    F20 = (float(np.sum((2.0 * Gam * q * F2 + S2j) / (q + qJ))) - f10) / den   # (111)
    c1 = -2.0 / qJ * F20 + c0 / alpha ** 2 * bg2                         # (114ii)

    W = 1.0 / (alpha * qJ) * (1.0 - 2.0 / Nexp * bg2 / alpha ** 2)       # (68)
    W2 = c0 + c1 / Nexp                                                  # (114i)
    return W, W2, c0, c1
