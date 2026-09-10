"""
Fork-join formulas of A. Thomasian, "Analysis of Fork/Join and Related Queueing
Systems", ACM Computing Surveys 47(2), Article 17, 2014.

Native Python port of matlab/src/api/fj/fj_qgb.m, fj_amva.m, fj_respt_closed.m,
fj_xmax_het.m, fj_lst_max_het.m, fj_xmax_moments_het.m, fj_xmax_hz.m,
fj_xmax_hz_het.m, fj_char_max_discrete.m, fj_char_max_blom.m, fj_cox_fit.m,
fj_xmax_coxian.m, fj_dispersion.m, fj_delay_opt.m, fj_respt_nosplit.m,
fj_respt_bulk.m, fj_ism_green.m, fj_tsm_capacity.m, fj_serialization.m and
fj_dag_makespan.m.
"""

import math
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.optimize import linprog
from scipy.special import erfcinv


def fj_harmonic(K: int) -> float:
    """K-th harmonic number, the expected maximum of K unit-rate exponentials."""
    if K < 1:
        raise ValueError(f"K must be a positive integer. Got K={K}.")
    return float(sum(1.0 / k for k in range(1, K + 1)))


# ---------------------------------------------------------------- closed models

def fj_qgb(D: Sequence[float], P: Sequence[int], M: int,
           Z: float = 0.0) -> Tuple[np.ndarray, np.ndarray]:
    """
    Geometric bound on the queue length of each fork-join subnetwork, Eq. (70).

        y_n(M) = D_n M / (Z + sum_j D_j H_{P_j} + Dmax M)
        Q_n(M) = H_{P_n} [ y_n/(1-y_n) - y_n^(M+1)/(1-y_n) ]

    The harmonic weights are what distinguishes this from the ordinary geometric
    bound of pfqn_qzgblow: a P-way subnetwork inflates its own demand by H_P in
    the denominator and its queue length by H_P in the numerator, and setting
    every fork degree to one recovers that bound exactly.
    """
    D = np.asarray(D, dtype=float).ravel()
    P = np.asarray(P, dtype=int).ravel()
    if D.size != P.size:
        raise ValueError(f"D and P must have the same number of elements. Got {D.size} and {P.size}.")
    if D.size == 0:
        raise ValueError("At least one subnetwork is required.")
    if np.any(D < 0):
        raise ValueError("Service demands must be non-negative.")
    if np.any(P < 1):
        raise ValueError("Fork degrees must be positive integers.")
    if M < 1:
        raise ValueError(f"M must be a positive integer. Got M={M}.")
    if Z < 0:
        raise ValueError("The think time must be non-negative.")

    H = np.array([fj_harmonic(int(p)) for p in P])
    Dtot = float(np.sum(D * H))
    Dmax = float(np.max(D))

    y = D * M / (Z + Dtot + Dmax * M)
    Q = np.empty_like(y)
    for n in range(D.size):
        if y[n] < 1.0:
            Q[n] = H[n] * (y[n] / (1 - y[n]) - y[n] ** (M + 1) / (1 - y[n]))
        else:
            # Degenerate ratio: the bound collapses onto the full population
            Q[n] = M
    return Q, y


def fj_amva(D: Sequence[float], P: Sequence[int], M: int,
            Z: float = 0.0) -> Tuple[np.ndarray, np.ndarray, float, np.ndarray]:
    """
    Mean value analysis of a closed network of fork-join subnetworks,
    Eqs. (68)-(69).

        R_n(m) = D_n [ H_{P_n} + Q_n(m-1) ]
        X(m)   = m / (Z + sum_n R_n(m))
        Q_n(m) = X(m) R_n(m)

    With every fork degree equal to one this is the exact single-class mean value
    analysis, because H_1 = 1; above that it is an approximation whose
    per-subnetwork residence time is an upper bound in the sense of Varki.
    """
    D = np.asarray(D, dtype=float).ravel()
    P = np.asarray(P, dtype=int).ravel()
    if D.size != P.size:
        raise ValueError(f"D and P must have the same number of elements. Got {D.size} and {P.size}.")
    if D.size == 0:
        raise ValueError("At least one subnetwork is required.")
    if np.any(D < 0):
        raise ValueError("Service demands must be non-negative.")
    if np.any(P < 1):
        raise ValueError("Fork degrees must be positive integers.")
    if M < 1:
        raise ValueError(f"M must be a positive integer. Got M={M}.")
    if Z < 0:
        raise ValueError("The think time must be non-negative.")

    H = np.array([fj_harmonic(int(p)) for p in P])
    Q = np.zeros(D.size)
    R = np.zeros(D.size)
    X = 0.0
    for m in range(1, M + 1):
        R = D * (H + Q)
        Rtot = float(np.sum(R))
        if not Rtot > 0:
            raise ValueError("The total residence time vanished; every demand is zero.")
        X = m / (Z + Rtot)
        Q = X * R
    # Each subnetwork holds P(n) queues sharing the demand equally
    U = X * D / P
    return R, Q, X, U


def fj_respt_closed(K: int, x: float, M: int,
                    A: Optional[float] = None) -> Tuple[float, bool]:
    """
    Varki bound on the residence time of a closed fork-join subnetwork,
    R_{P_K}(M) <= x [ H_K + A ], Eq. (67).

    With A omitted the subnetwork is the whole closed network, so every other job
    is inside it and A = M-1; that bound is tight at K = 2 (Theorem 4.1).
    """
    if K < 1:
        raise ValueError(f"K must be a positive integer. Got K={K}.")
    if not x > 0:
        raise ValueError(f"The mean service time must be positive. Got x={x}.")
    if M < 1:
        raise ValueError(f"M must be a positive integer. Got M={M}.")
    isolated = A is None
    if isolated:
        A = M - 1
    if A < 0:
        raise ValueError(f"The arrival-instant queue length must be non-negative. Got A={A}.")
    return x * (fj_harmonic(K) + A), bool(isolated and K == 2)


# ------------------------------------------------------------------- order stats

def fj_xmax_het(lam: Sequence[float], n: int = 1) -> float:
    """
    Exact n-th moment of the maximum of K heterogeneous exponentials, by
    inclusion-exclusion on the survival function:

        E[Y^n] = sum over the nonempty subsets S of (-1)^(|S|+1) n! / (sum_S l_i)^n

    At n = 1 and K = 2 this collapses to 1/l1 + 1/l2 - 1/(l1+l2), and for equal
    rates to H_K/lambda.
    """
    lam = np.asarray(lam, dtype=float).ravel()
    K = lam.size
    if K < 1:
        raise ValueError("At least one rate is required.")
    if np.any(lam <= 0):
        raise ValueError("All exponential rates must be positive.")
    if n < 1:
        raise ValueError(f"The moment order must be a positive integer. Got n={n}.")
    if K > 24:
        raise ValueError(
            f"Inclusion-exclusion over K={K} rates needs 2^K terms; use fj_xmax_moments_het.")
    nfact = float(math.factorial(n))
    acc = 0.0
    for mask in range(1, 1 << K):
        rate = 0.0
        card = 0
        for i in range(K):
            if mask & (1 << i):
                rate += lam[i]
                card += 1
        term = nfact / rate ** n
        acc += term if card % 2 == 1 else -term
    return acc


def fj_lst_max_het(lam: Sequence[float], s: float) -> float:
    """
    Laplace-Stieltjes transform of the maximum of heterogeneous exponentials,
    by the Harrison-Zertal recurrence, Eq. (29):

        ( s + sum_{j<=m} l_j ) L*_m = sum_{j<=m} l_j L*_{m-1}(lam \\ j)

    anchored at L*_0 = 1 because the maximum of an empty collection is zero.
    """
    lam = np.asarray(lam, dtype=float).ravel()
    K = lam.size
    if K < 1:
        raise ValueError("At least one rate is required.")
    if np.any(lam <= 0):
        raise ValueError("All exponential rates must be positive.")
    if s < 0:
        raise ValueError("The transform argument must be non-negative.")
    if K > 22:
        raise ValueError(f"The recurrence enumerates 2^K sub-collections; K={K} is too large.")
    nmask = 1 << K
    tab = np.zeros(nmask)
    tab[0] = 1.0
    for mask in range(1, nmask):
        tot = 0.0
        acc = 0.0
        for j in range(K):
            if mask & (1 << j):
                tot += lam[j]
                acc += lam[j] * tab[mask ^ (1 << j)]
        tab[mask] = acc / (s + tot)
    return float(tab[nmask - 1])


def fj_xmax_moments_het(lam: Sequence[float], n: int = 1) -> np.ndarray:
    """
    Moments of orders 1..n of the maximum of heterogeneous exponentials, by the
    n-th derivative of the transform recurrence at the origin:

        M_m(n) = [ n M_m(n-1) + sum_j l_j M_{m-1}(lam \\ j, n) ] / sum_j l_j

    Eq. (30) of the survey prints the second sum WITHOUT the l_j weight; that
    form is not the derivative of Eq. (29) and misses the textbook two-variable
    answer, so the weight is restored here. fj_xmax_het is the independent
    inclusion-exclusion check.
    """
    lam = np.asarray(lam, dtype=float).ravel()
    K = lam.size
    if K < 1:
        raise ValueError("At least one rate is required.")
    if np.any(lam <= 0):
        raise ValueError("All exponential rates must be positive.")
    if n < 1:
        raise ValueError(f"The moment order must be a positive integer. Got n={n}.")
    if K > 22:
        raise ValueError(f"The recurrence enumerates 2^K sub-collections; K={K} is too large.")
    nmask = 1 << K
    tab = np.zeros((nmask, n + 1))
    tab[:, 0] = 1.0
    for order in range(1, n + 1):
        # The empty sub-collection has a zero maximum, so all its moments vanish
        tab[0, order] = 0.0
        for mask in range(1, nmask):
            tot = 0.0
            acc = order * tab[mask, order - 1]
            for j in range(K):
                if mask & (1 << j):
                    tot += lam[j]
                    acc += lam[j] * tab[mask ^ (1 << j), order]
            tab[mask, order] = acc / tot
    return tab[nmask - 1, 1:]


def fj_xmax_hz(m1: float, m2: float, K: int) -> Tuple[float, float]:
    """
    Harrison-Zertal closed form for K identically distributed branches,
    X_K^max ~ m1 + ( m2/(2 m1) ) ( H_K - 1 ).

    The correction is the equilibrium mean of the branch law scaled by H_K - 1:
    one branch, plus the residual work still owed by the branches that finish
    later. It is exact for the exponential and reduces to m1 at K = 1.
    """
    if not m1 > 0:
        raise ValueError(f"The branch mean must be positive. Got m1={m1}.")
    if m2 < m1 ** 2:
        raise ValueError(f"The second moment {m2} is below m1^2={m1 ** 2}.")
    if K < 1:
        raise ValueError(f"K must be a positive integer. Got K={K}.")
    resid = m2 / (2 * m1)
    return m1 + resid * (fj_harmonic(K) - 1), resid


def fj_xmax_hz_het(m1: Sequence[float], m2: Sequence[float],
                   cdf: Sequence[Callable[[np.ndarray], np.ndarray]],
                   tol: float = 1e-10, npoints: int = 2001) -> float:
    """
    Harrison-Zertal recurrence for independent but not identically distributed
    branches, Eq. (46):

        I(S) = (1/|S|) sum_{i in S} [ I(S \\ i) + (m2_i/(2 m1_i)) L*_{S\\i}(1/m1_i) ]

    anchored at I({i}) = m1_i. The transform of the maximum over a
    sub-collection is recovered from the product of the distribution functions by
    composite Simpson quadrature on a horizon widened until the product is within
    TOL of one.
    """
    m1 = np.asarray(m1, dtype=float).ravel()
    m2 = np.asarray(m2, dtype=float).ravel()
    K = m1.size
    if m2.size != K or len(cdf) != K:
        raise ValueError("m1, m2 and cdf must have the same length.")
    if K < 1:
        raise ValueError("At least one branch is required.")
    if np.any(m1 <= 0):
        raise ValueError("All branch means must be positive.")
    if np.any(m2 < m1 ** 2):
        raise ValueError("Some second moment is below the square of its mean.")
    if K > 14:
        raise ValueError(f"The recurrence enumerates 2^K sub-collections with a quadrature each; K={K}.")
    if npoints % 2 == 0:
        npoints += 1

    alpha = 1.0 / m1
    resid = m2 / (2 * m1)

    # Horizon: widen until every branch is essentially complete
    U = 8 * float(np.max(m1))
    for _ in range(60):
        prodF = 1.0
        for j in range(K):
            prodF *= float(cdf[j](np.array([U]))[0])
        if 1 - prodF < tol:
            break
        U *= 2

    t = np.linspace(0.0, U, npoints)
    h = t[1] - t[0]
    w = np.ones(npoints)
    w[1:-1:2] = 4
    w[2:-2:2] = 2
    Fvals = [np.asarray(cdf[j](t), dtype=float) for j in range(K)]

    def lst_max(mask: int, s: float) -> float:
        if mask == 0:
            return 1.0
        g = np.exp(-s * t)
        for j in range(K):
            if mask & (1 << j):
                g = g * Fvals[j]
        return float(s * (h / 3) * np.sum(w * g))

    nmask = 1 << K
    Ival = np.zeros(nmask)
    for mask in range(1, nmask):
        members = [i for i in range(K) if mask & (1 << i)]
        if len(members) == 1:
            Ival[mask] = m1[members[0]]
            continue
        acc = 0.0
        for i in members:
            rest = mask ^ (1 << i)
            acc += Ival[rest] + resid[i] * lst_max(rest, float(alpha[i]))
        Ival[mask] = acc / len(members)
    return float(Ival[nmask - 1])


def fj_char_max_discrete(K: int, dist_type: str, param: float) -> Tuple[float, int, float]:
    """
    Gravey's characteristic maximum for a lattice law, Eq. (48). With m_K the
    smallest integer at which P(X > m_K) <= 1/K,

        M_K = m_K + K sum_{k >= m_K} P(X > k),

    which upper bounds the expected maximum of K i.i.d. copies at O(1) instead of
    the alternating binomial sum. Returns (M_K, m_K, exact).
    """
    if K < 1:
        raise ValueError(f"K must be a positive integer. Got K={K}.")
    kind = dist_type.lower()
    if kind in ("geometric", "geom"):
        p = param
        if not 0 < p < 1:
            raise ValueError(f"The geometric parameter must lie in (0,1). Got p={p}.")
        mK = max(0, int(math.ceil(-math.log(K) / math.log(p) - 1e-12)))
        MK = mK + K * p ** (mK + 1) / (1 - p)
        # Exact maximum by inclusion-exclusion on the geometric tail
        exact = 0.0
        for k in range(1, K + 1):
            term = math.comb(K, k) * p ** k / (1 - p ** k)
            exact += term if k % 2 == 1 else -term
        return MK, mK, exact
    if kind in ("poisson", "pois"):
        theta = param
        if not theta > 0:
            raise ValueError(f"The Poisson mean must be positive. Got theta={theta}.")
        kmax = int(math.ceil(theta + 12 * math.sqrt(theta) + 40))
        k = np.arange(kmax + 1)
        with np.errstate(divide="ignore"):
            logpmf = -theta + k * np.log(theta) - np.array([math.lgamma(v + 1) for v in k])
        pmf = np.exp(logpmf)
        cdf = np.minimum(np.cumsum(pmf), 1.0)
        tail = 1.0 - cdf
        idx = np.nonzero(tail <= 1.0 / K)[0]
        if idx.size == 0:
            raise ValueError(f"The Poisson lattice truncation at {kmax} never reached a tail of 1/K.")
        mK = int(idx[0])
        tail_prev = 1.0 if mK == 0 else float(tail[mK - 1])
        MK = mK * (1 - K * float(tail[mK])) + K * theta * tail_prev
        exact = float(np.sum(1.0 - cdf ** K))
        return MK, mK, exact
    raise ValueError(f'Unsupported discrete distribution "{dist_type}". Use geometric or poisson.')


def fj_char_max_blom(K: int, Finv: Optional[Callable[[float], float]] = None,
                     alpha: float = 0.4886,
                     beta: float = 0.3140) -> Tuple[float, float, float]:
    """
    Blom-corrected plotting position for the characteristic maximum,
    m_K = F^-1( (K - alpha) / (K - alpha - beta + 1) ), which for
    alpha = beta = 0 falls back on the naive K/(K+1).

    For the standard normal the position is bracketed for K >= 5 by
    sqrt(2 ln K - ln ln K - 3) < m_K < sqrt(2 ln K - ln ln K); those two are
    returned as lo and hi, NaN below K = 5.
    """
    if K < 1:
        raise ValueError(f"K must be a positive integer. Got K={K}.")
    den = K - alpha - beta + 1
    if not den > 0:
        raise ValueError(f"The Blom offsets leave a non-positive denominator at K={K}.")
    q = (K - alpha) / den
    if not 0 < q < 1:
        raise ValueError(f"The Blom plotting position {q} fell outside (0,1).")
    # sqrt(2) erfinv(2q-1) == -sqrt(2) erfcinv(2q)
    mK = float(Finv(q)) if Finv is not None else float(-math.sqrt(2.0) * erfcinv(2 * q))
    lo = math.nan
    hi = math.nan
    if K >= 5:
        z = 2 * math.log(K) - math.log(math.log(K))
        if z > 3:
            lo = math.sqrt(z - 3)
        if z > 0:
            hi = math.sqrt(z)
    return mK, lo, hi


def fj_cox_fit(m1: float, c2: float) -> Tuple[float, float, float, int, int]:
    """
    Marie's balanced-stage fit of a two-stage Coxian law to a target mean and
    squared coefficient of variation. Requiring 1/mu1 = q/mu2 closes the system
    of two moment equations in three unknowns and gives

        mu1 = 2 mu,   q = 1/(2 c2),   mu2 = mu/c2,

    which needs q <= 1, hence c2 >= 0.5. The admissible Erlang stage count for
    the same target is bracketed by ceil(1/c2) <= k <= floor(1/c2) + 1.
    """
    if not m1 > 0:
        raise ValueError(f"The target mean must be positive. Got m1={m1}.")
    if c2 < 0.5:
        raise ValueError(
            f"The balanced-stage Coxian fit needs c2 >= 0.5. Got c2={c2}; "
            f"use an Erlang with {max(1, math.ceil(1 / c2))} stages instead.")
    mu = 1.0 / m1
    kmin = max(1, int(math.ceil(1 / c2 - 1e-12)))
    kmax = max(kmin, int(math.floor(1 / c2 + 1e-12)) + 1)
    return 2 * mu, mu / c2, 1.0 / (2 * c2), kmin, kmax


def fj_xmax_coxian(K: int, mu1: float, mu2: float, q: float) -> Tuple[float, float, float]:
    """
    Exact expected maximum of K i.i.d. two-stage Coxian branches. The survival
    function is a two-term exponential mixture, so expanding 1 - (1-S)^K
    binomially and integrating term by term is closed:

        E[Y_K] = sum_j (-1)^(j+1) C(K,j) sum_i C(j,i) A^(j-i) B^i / ((j-i) mu1 + i mu2)

    At coincident stage rates the mixture degenerates into
    S(t) = (1 + q mu t) exp(-mu t), handled by the same expansion with the
    polynomial integrals. Returns (Xmax, m1, c2).
    """
    if K < 1:
        raise ValueError(f"K must be a positive integer. Got K={K}.")
    if not mu1 > 0 or not mu2 > 0:
        raise ValueError(f"Both stage rates must be positive. Got mu1={mu1}, mu2={mu2}.")
    if not 0 <= q <= 1:
        raise ValueError(f"The branching probability must lie in [0,1]. Got q={q}.")
    if K > 60:
        raise ValueError(f"The binomial expansion loses precision past K=60. Got K={K}.")

    m1 = 1 / mu1 + q / mu2
    var1 = 1 / mu1 ** 2 + q * (2 - q) / mu2 ** 2
    c2 = var1 / m1 ** 2

    acc = 0.0
    if abs(mu2 - mu1) > 1e-12 * max(mu1, mu2):
        A = (1 - q) + q * mu2 / (mu2 - mu1)
        B = -q * mu1 / (mu2 - mu1)
        for j in range(1, K + 1):
            inner = 0.0
            for i in range(j + 1):
                rate = (j - i) * mu1 + i * mu2
                inner += math.comb(j, i) * A ** (j - i) * B ** i / rate
            term = math.comb(K, j) * inner
            acc += term if j % 2 == 1 else -term
    else:
        # Coincident stage rates: S(t) = (1 + q mu t) exp(-mu t)
        mu = mu1
        for j in range(1, K + 1):
            inner = 0.0
            jmu = j * mu
            for i in range(j + 1):
                inner += math.comb(j, i) * (q * mu) ** i * math.factorial(i) / jmu ** (i + 1)
            term = math.comb(K, j) * inner
            acc += term if j % 2 == 1 else -term
    return acc, m1, c2


# -------------------------------------------------------------------- dispersion

def _erlang_cdf(t: np.ndarray, k: int, mu: float) -> np.ndarray:
    """Erlang-k distribution function, zero on the negative half line."""
    t = np.asarray(t, dtype=float)
    out = np.zeros_like(t)
    pos = t > 0
    tp = t[pos]
    acc = np.zeros_like(tp)
    term = np.ones_like(tp)
    for j in range(k):
        if j > 0:
            term = term * (mu * tp) / j
        acc = acc + term
    out[pos] = 1 - np.exp(-mu * tp) * acc
    return out


def fj_dispersion(shape: Sequence[int], rate: Sequence[float],
                  d: Optional[Sequence[float]] = None, tol: float = 1e-10,
                  npoints: int = 4001) -> Tuple[float, float, float]:
    """
    Mean subtask dispersion of a split-merge system with shifted Erlang branches,
    Eqs. (24)-(25):

        E[D_d] = int_0^inf [ 1 - prod_i F_i(x-d_i) - prod_i (1-F_i(x-d_i)) ] dx

    That integrand is non-negative and vanishes at both ends; the difference of
    the two products printed in the survey is not the dispersion and can go
    negative. Returns (Edisp, Emax, Emin).
    """
    shape = np.asarray(shape, dtype=int).ravel()
    rate = np.asarray(rate, dtype=float).ravel()
    N = shape.size
    if d is None:
        d = np.zeros(N)
    d = np.asarray(d, dtype=float).ravel()
    if rate.size != N or d.size != N:
        raise ValueError("shape, rate and d must have the same length.")
    if N < 1:
        raise ValueError("At least one branch is required.")
    if np.any(shape < 1):
        raise ValueError("Erlang stage counts must be positive integers.")
    if np.any(rate <= 0):
        raise ValueError("Erlang stage rates must be positive.")
    if np.any(d < 0):
        raise ValueError("Delays must be non-negative.")
    if npoints % 2 == 0:
        npoints += 1

    U = float(np.max(d) + 8 * np.max(shape / rate))
    for _ in range(60):
        prodF = 1.0
        for i in range(N):
            prodF *= float(_erlang_cdf(np.array([U - d[i]]), int(shape[i]), float(rate[i]))[0])
        if 1 - prodF < tol:
            break
        U *= 2

    x = np.linspace(0.0, U, npoints)
    h = x[1] - x[0]
    w = np.ones(npoints)
    w[1:-1:2] = 4
    w[2:-2:2] = 2

    Fprod = np.ones(npoints)
    Sprod = np.ones(npoints)
    for i in range(N):
        Fi = _erlang_cdf(x - d[i], int(shape[i]), float(rate[i]))
        Fprod = Fprod * Fi
        Sprod = Sprod * (1 - Fi)

    Emax = float((h / 3) * np.sum(w * (1 - Fprod)))
    Emin = float((h / 3) * np.sum(w * Sprod))
    return Emax - Emin, Emax, Emin


def fj_delay_opt(shape: Sequence[int], rate: Sequence[float], maxsweeps: int = 40,
                 dtol: float = 1e-8, npoints: int = 2001) -> Tuple[np.ndarray, float, float]:
    """
    The deterministic delays that minimise the mean dispersion. Holding back a
    fast branch costs little at the last completion and buys a great deal at the
    first, so the minimiser is generally interior and strictly positive on every
    branch but the slowest.

    Cyclic coordinate descent with a golden section line search on each
    coordinate: deterministic, derivative-free, and the same sequence of
    evaluations in all four codebases. Adding a constant to every delay shifts
    both order statistics equally, so the representative with min(d) = 0 is
    returned.
    """
    shape = np.asarray(shape, dtype=int).ravel()
    rate = np.asarray(rate, dtype=float).ravel()
    N = shape.size
    if rate.size != N:
        raise ValueError("shape and rate must have the same length.")
    if N < 1:
        raise ValueError("At least one branch is required.")
    d = np.zeros(N)
    if N < 2:
        Edisp, Emax, _ = fj_dispersion(shape, rate, d, npoints=npoints)
        return d, Edisp, Emax

    means = shape / rate
    ub = float(np.max(means) + 8 * np.max(np.sqrt(shape) / rate))
    invphi = (math.sqrt(5.0) - 1) / 2

    def obj(dv: np.ndarray) -> float:
        return fj_dispersion(shape, rate, dv, npoints=npoints)[0]

    fcur = obj(d)
    for _ in range(maxsweeps):
        fprev = fcur
        for i in range(N):
            a, b = 0.0, ub
            c = b - invphi * (b - a)
            dd = a + invphi * (b - a)
            probe = d.copy()
            probe[i] = c
            fc = obj(probe)
            probe[i] = dd
            fd = obj(probe)
            for _ in range(60):
                if fc < fd:
                    b, dd, fd = dd, c, fc
                    c = b - invphi * (b - a)
                    probe[i] = c
                    fc = obj(probe)
                else:
                    a, c, fc = c, dd, fd
                    dd = a + invphi * (b - a)
                    probe[i] = dd
                    fd = obj(probe)
                if (b - a) <= dtol * max(1.0, ub):
                    break
            d[i] = c if fc < fd else dd
        d = d - np.min(d)
        fcur = obj(d)
        if abs(fprev - fcur) <= dtol * max(1.0, abs(fprev)):
            break

    Edisp, Emax, _ = fj_dispersion(shape, rate, d, npoints=npoints)
    return d, Edisp, Emax


# ------------------------------------------------------- parallel-processing models

def fj_respt_nosplit(K: int, lam: float, mu: float) -> Tuple[float, float]:
    """
    Distributed no splitting: a job of K tasks goes in one piece to a single
    server chosen uniformly among the K, which is an M/E_K/1 queue and reduces to
    R = [ K - (K-1) rho/2 ] / (mu - lambda). At K = 1 it is the M/M/1 response
    time. Returns (R, rho).
    """
    if K < 1:
        raise ValueError(f"K must be a positive integer. Got K={K}.")
    if not lam > 0:
        raise ValueError(f"The arrival rate must be positive. Got lambda={lam}.")
    if not mu > 0:
        raise ValueError(f"The service rate must be positive. Got mu={mu}.")
    rho = lam / mu
    if rho >= 1:
        raise ValueError(f"System is unstable: rho = lambda/mu = {rho:.4f} >= 1.")
    return (K - (K - 1) * rho / 2) / (mu - lam), rho


def fj_respt_bulk(K: int, lam: float, mu: float, c: int,
                  nmax: Optional[int] = None) -> Tuple[float, float, float, np.ndarray]:
    """
    Centralized splitting as an M[K]/M/c bulk arrival system, solved by
    truncating the level chain. The request response time is the completion of
    the LAST of the K tasks: by PASTA the batch finds n tasks in system, its last
    task is the (n+K)-th in line, and under first come first served with c
    exponential servers it starts after max(0, n+K-c) departures, each an
    exponential of rate c mu. Returns (Rreq, Rtask, Q, p).
    """
    if K < 1:
        raise ValueError(f"K must be a positive integer. Got K={K}.")
    if c < 1:
        raise ValueError(f"c must be a positive integer. Got c={c}.")
    if not lam > 0 or not mu > 0:
        raise ValueError(f"lambda and mu must be positive. Got lambda={lam}, mu={mu}.")
    rho = lam * K / (c * mu)
    if rho >= 1:
        raise ValueError(f"System is unstable: rho = lambda*K/(c*mu) = {rho:.4f} >= 1.")
    if nmax is None:
        nmax = int(max(200, math.ceil(K + 40 * c / (1 - rho))))
    ns = nmax + 1

    Q = np.zeros((ns, ns))
    for i in range(ns):
        if i > 0:
            srv = min(i, c) * mu
            Q[i, i - 1] += srv
            Q[i, i] -= srv
        j = i + K
        if j < ns:
            Q[i, j] += lam
            Q[i, i] -= lam

    A = np.vstack([Q.T, np.ones((1, ns))])
    b = np.zeros(ns + 1)
    b[-1] = 1.0
    p, *_ = np.linalg.lstsq(A, b, rcond=None)
    p = np.maximum(p, 0.0)
    p = p / np.sum(p)

    n = np.arange(ns)
    Qmean = float(np.sum(n * p))
    Rtask = Qmean / (lam * K)
    wait = np.maximum(0, n + K - c) / (c * mu)
    Rreq = float(np.sum(p * (wait + 1 / mu)))
    return Rreq, Rtask, Qmean, p


def fj_ism_green(lam: float, mu: float, s: int,
                 c: Sequence[float]) -> Tuple[float, float, Dict[str, object]]:
    """
    Green's independent server model: a customer needs j servers at once with
    probability c(j) and releases them asynchronously as each of its j tasks
    completes at rate mu, so its own service is the maximum of j exponentials.

    E[B] is the j-th order statistic of s exponentials, because all s servers are
    busy whenever a customer enters service during a queueing period; E[D] is the
    initial delay of the customer that starts one; and the waiting-time transform
    of Eq. (61) factors into the equilibrium transform of D and the
    Pollaczek-Khinchine transform with service B, so the means add.

    Eq. (65) of the survey prints the inner sum of E[D] as starting at 1/(s mu)
    even though only i servers are busy; it is started at 1/(i mu) here, which is
    what the accompanying text prescribes and what makes E[D] reduce to E[B] at
    i = s. Returns (W, R, out).
    """
    c = np.asarray(c, dtype=float).ravel()
    if s < 1:
        raise ValueError(f"s must be a positive integer. Got s={s}.")
    if c.size != s:
        raise ValueError(f"c must have s={s} entries, one per server requirement. Got {c.size}.")
    if np.any(c < 0):
        raise ValueError("The server-requirement probabilities must be non-negative.")
    if abs(float(np.sum(c)) - 1) > 1e-9:
        raise ValueError(f"The server-requirement probabilities must sum to one. Got {np.sum(c)}.")
    if not lam > 0 or not mu > 0:
        raise ValueError(f"lambda and mu must be positive. Got lambda={lam}, mu={mu}.")

    EB = 0.0
    EB2 = 0.0
    for j in range(1, s + 1):
        stages = 1.0 / ((s - np.arange(j)) * mu)
        mj = float(np.sum(stages))
        vj = float(np.sum(stages ** 2))
        EB += c[j - 1] * mj
        EB2 += c[j - 1] * (vj + mj ** 2)

    rho = lam * EB
    if rho >= 1:
        raise ValueError(f"System is unstable: rho = lambda*E[B] = {rho:.4f} >= 1.")

    n = s + 1
    T = np.zeros((n, n))
    for i in range(n):
        tot = lam + i * mu
        if i > 0:
            T[i, i - 1] = i * mu / tot
        for j in range(1, s - i + 1):
            T[i, i + j] += lam * c[j - 1] / tot

    V = np.linalg.solve(np.eye(n) - T, np.eye(n))
    # A nonqueue period starts with all s servers busy
    v = V[s, :]
    hold = 1.0 / (lam + np.arange(n) * mu)
    EQbar = float(np.sum(v * hold))
    q = (v * hold) / EQbar

    pd = 0.0
    for i in range(n):
        free = s - i
        if free < s:
            pd += q[i] * float(np.sum(c[free:s]))
    if not pd > 0:
        raise ValueError(f"No arrival can ever be delayed; the model degenerates to M/M/{s}.")

    ED = 0.0
    ED2 = 0.0
    for i in range(1, s + 1):
        for k in range(1, i + 1):
            j = s - i + k
            if j < 1 or j > s:
                continue
            wgt = q[i] * c[j - 1] / pd
            if wgt == 0:
                continue
            stages = 1.0 / ((i - np.arange(k)) * mu)
            mk = float(np.sum(stages))
            vk = float(np.sum(stages ** 2))
            ED += wgt * mk
            ED2 += wgt * (vk + mk ** 2)

    EQ = ED / (1 - rho)
    pq = EQ / (EQ + EQbar)
    pi0 = (1 - rho) / (1 - lam * (EB - ED))
    W = (1 - pi0) * (ED2 / (2 * ED) + lam * EB2 / (2 * (1 - rho)))
    ES = float(sum(c[j - 1] * fj_harmonic(j) / mu for j in range(1, s + 1)))

    out = {"EB": EB, "EB2": EB2, "ED": ED, "ED2": ED2, "EQ": EQ, "EQbar": EQbar,
           "pq": pq, "pd": pd, "pi0": pi0, "rho": rho, "q": q, "ES": ES}
    return W, W + ES, out


def _tsm_states(r: np.ndarray, s: int) -> List[List[int]]:
    """Every multiset of jobs whose total server demand is at most s, minus the empty one."""
    K = r.size
    out: List[List[int]] = []
    stack = [0] * K

    def rec(k: int, left: int) -> None:
        if k == K:
            if any(v > 0 for v in stack):
                out.append(list(stack))
            return
        for n in range(left // int(r[k]) + 1):
            stack[k] = n
            rec(k + 1, left - n * int(r[k]))
        stack[k] = 0

    rec(0, s)
    return out


def fj_tsm_capacity(s: int, f: Sequence[float], r: Sequence[int],
                    x: Sequence[float]) -> Tuple[float, float, float, Dict[str, object]]:
    """
    Saturation throughput of the team service model, Eqs. (57)-(58). The apparent
    rate Lambda_max = s / sum_k f(k) r(k) x(k) is attainable only when the
    scheduler can pack jobs into execution states that leave no server idle; the
    attainable capacity is the optimum of

        max Lambda s.t. sum_j p_j n(j,k)/x(k) = Lambda f(k), sum_j p_j = 1, p >= 0

    over the feasible execution states. For the two-server two-class case with
    r = (1,2), strict first come first served cannot pack at all and reaches only

        lambda_FCFS = 2 mu1 mu2 / (f1^2 mu2 + 2 f2^2 mu1 + 2 f1 f2 (mu1+mu2)).

    Returns (Lmax, Llp, Lfcfs, detail); Lfcfs is NaN outside that special case.
    """
    f = np.asarray(f, dtype=float).ravel()
    r = np.asarray(r, dtype=int).ravel()
    x = np.asarray(x, dtype=float).ravel()
    K = f.size
    if s < 1:
        raise ValueError(f"s must be a positive integer. Got s={s}.")
    if r.size != K or x.size != K:
        raise ValueError("f, r and x must have the same length.")
    if np.any(f < 0) or abs(float(np.sum(f)) - 1) > 1e-9:
        raise ValueError(f"The class frequencies must be non-negative and sum to one. Got {np.sum(f)}.")
    if np.any(r < 1) or np.any(r > s):
        raise ValueError(f"The server requirements must be integers in 1..{s}.")
    if np.any(x <= 0):
        raise ValueError("The mean service times must be positive.")

    Lmax = s / float(np.sum(f * r * x))

    states = _tsm_states(r, s)
    ns = len(states)
    if ns == 0:
        raise ValueError("No feasible execution state.")
    S = np.array(states, dtype=float)

    active = np.nonzero(f > 0)[0]
    Aeq = np.zeros((active.size + 1, ns + 1))
    beq = np.zeros(active.size + 1)
    for a, k in enumerate(active):
        Aeq[a, :ns] = S[:, k] / x[k]
        Aeq[a, ns] = -f[k]
    Aeq[-1, :ns] = 1.0
    beq[-1] = 1.0

    cost = np.zeros(ns + 1)
    cost[ns] = -1.0
    res = linprog(cost, A_eq=Aeq, b_eq=beq, bounds=[(0, None)] * (ns + 1), method="highs")
    if not res.success:
        raise ValueError(f"The capacity linear program did not solve: {res.message}")
    Llp = float(res.x[ns])

    Lfcfs = math.nan
    if s == 2 and K == 2 and sorted(r.tolist()) == [1, 2]:
        a = 0 if r[0] == 1 else 1
        b = 1 - a
        f1, f2 = float(f[a]), float(f[b])
        mu1, mu2 = 1 / float(x[a]), 1 / float(x[b])
        Lfcfs = 2 * mu1 * mu2 / (f1 ** 2 * mu2 + 2 * f2 ** 2 * mu1 + 2 * f1 * f2 * (mu1 + mu2))

    return Lmax, Llp, Lfcfs, {"states": states, "prob": res.x[:ns]}


def fj_serialization(Rs: Sequence[float], R0: float, M: int,
                     alpha: float = 0.5) -> Tuple[np.ndarray, np.ndarray, float]:
    """
    Blocking probability and pseudoserver delay of serialization phases:

        P_s(M) = 1 - [ 1 - R_s(M)/R(M) ]^(M-1),   R(M) = R_0 + sum_s R_s(M),

    with the delay charged at the pseudoserver equal to alpha R_s(M); alpha = 1/2
    is the value for an arrival uniform in a lightly utilized phase, the regime in
    which the approximation is stated. Returns (P, delay, Rtot).
    """
    Rs = np.asarray(Rs, dtype=float).ravel()
    if Rs.size < 1:
        raise ValueError("At least one serialization phase must be supplied.")
    if np.any(Rs < 0):
        raise ValueError("The residence times inside the serialization phases must be non-negative.")
    if R0 < 0:
        raise ValueError(f"The nonserialized residence time must be non-negative. Got R0={R0}.")
    if M < 1:
        raise ValueError(f"M must be a positive integer. Got M={M}.")
    if not 0 <= alpha <= 1:
        raise ValueError(f"alpha must lie in [0,1]. Got alpha={alpha}.")

    R = R0 + float(np.sum(Rs))
    if not R > 0:
        raise ValueError("The total residence time vanished; every phase has zero demand.")
    P = 1 - (1 - Rs / R) ** (M - 1)
    delay = P * (alpha * Rs)
    return P, delay, R + float(np.sum(delay))


def fj_dag_makespan(pred: np.ndarray,
                    rate: np.ndarray) -> Tuple[float, np.ndarray, np.ndarray, np.ndarray]:
    """
    Makespan of a task system with precedence constraints, Eq. (66) and
    Section 7.4. The chain whose state is the SET of completed tasks is acyclic,
    so it is swept level by level: task i among the k eligible completes at rate
    rate[i, k-1], the state is held for 1/T(S), and

        p(R) = sum p(S) b(S,R),  D(R) = M(R) p(R) + sum b(S,R) D(S),

    started at p(empty) = 1. Making the rate depend on the concurrency is what
    couples the task system to the queueing network underneath it: two tasks
    sharing a processor each run slower than either would alone.

    Returns (C, I, Cend, E).
    """
    pred = np.asarray(pred)
    n = pred.shape[0]
    if pred.shape[1] != n:
        raise ValueError(f"pred must be square. Got {pred.shape}.")
    if n < 1:
        raise ValueError("At least one task is required.")
    if n > 20:
        raise ValueError(f"The completed-set sweep enumerates 2^n states; n={n} is too large.")
    rate = np.asarray(rate, dtype=float)
    if rate.ndim == 1:
        if rate.size != n:
            raise ValueError(f"The rate vector must have one entry per task. Got {rate.size}.")
        rate = np.tile(rate.reshape(-1, 1), (1, n))
    elif rate.shape != (n, n):
        raise ValueError(f"The rate matrix must be {n}x{n}. Got {rate.shape}.")
    if np.any(rate <= 0):
        raise ValueError("All completion rates must be positive.")

    predmask = np.zeros(n, dtype=int)
    indeg = np.zeros(n, dtype=int)
    for j in range(n):
        for i in range(n):
            if pred[i, j]:
                predmask[j] |= (1 << i)
                indeg[j] += 1
    seen = np.zeros(n, dtype=bool)
    remaining = n
    for _ in range(n):
        pick = -1
        for i in range(n):
            if not seen[i] and indeg[i] == 0:
                pick = i
                break
        if pick < 0:
            break
        seen[pick] = True
        indeg[pick] = -1
        for j in range(n):
            if pred[pick, j] and indeg[j] > 0:
                indeg[j] -= 1
        remaining -= 1
    if remaining > 0:
        raise ValueError("The precedence relation contains a cycle.")

    nmask = 1 << n
    eligmask = np.zeros(nmask, dtype=int)
    closed = np.zeros(nmask, dtype=bool)
    for mask in range(nmask):
        ok = True
        em = 0
        for i in range(n):
            bit = 1 << i
            if mask & bit:
                # A completed task must have all of its predecessors completed
                if (predmask[i] & mask) != predmask[i]:
                    ok = False
                    break
            elif (predmask[i] & mask) == predmask[i]:
                em |= bit
        closed[mask] = ok
        if ok:
            eligmask[mask] = em

    p = np.zeros(nmask)
    D = np.zeros(nmask)
    p[0] = 1.0
    I = np.zeros(n)
    Cend = np.zeros(n)

    for mask in range(nmask):
        if not closed[mask]:
            continue
        em = int(eligmask[mask])
        elig = [i for i in range(n) if em & (1 << i)]
        k = len(elig)
        if k == 0:
            continue
        Ttot = float(sum(rate[i, k - 1] for i in elig))
        # Holding time of this state, weighted by the probability of reaching it
        D[mask] += p[mask] / Ttot
        for i in elig:
            b = rate[i, k - 1] / Ttot
            nxt = mask | (1 << i)
            contrib = b * D[mask]
            p[nxt] += p[mask] * b
            D[nxt] += contrib
            # Task i completes on this transition
            Cend[i] += contrib
            # Tasks that first become eligible on this transition start on it
            fresh = int(eligmask[nxt]) & ~em
            for j in range(n):
                if fresh & (1 << j):
                    I[j] += contrib

    return float(D[nmask - 1]), I, Cend, Cend - I
