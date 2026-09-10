"""
Hierarchical and multiserver/load-dependent bound methods for closed
product-form queueing networks. Native-Python ports (no JVM) of the MATLAB
pfqn_{pbh,cbh,pbk,bjbk,mcub,ssd,sib,ldbcmp,scb} bound functions used by
SolverBA, plus the class-aggregation error and class-count bounds
pfqn_{scbgap,usumbound,minclasses}.

References:
- Eager-Sevcik 1983 (PBH), Dowdy et al. 1984 (CBH), Casale-Muntz-Serazzi 2008
  (iterative PB(k)/BJB(k)), Kerola 1986 (multiclass composite upper bound),
  Suri-Dallery 1986 (multiserver disaggregation), Srinivasan 1985 (SIB),
  Anselmi-Cremonesi 2008 (LD-BCMP closed-open equivalence),
  Dowdy-Carlson-Krantz-Tripathi 1992 (single-class bounds of multiclass
  networks, and the class-count bound of their Section 4.7).
"""

from typing import Tuple

import numpy as np
from scipy.special import gammaln
from math import factorial

from .utils import _amva_prep


# ----------------------------- PBH (Eager-Sevcik) --------------------------

def _pbh_residence(L, N, Z, level, side):
    """Per-station residence-time vector of the level-`level` PBH bound."""
    L = np.asarray(L, dtype=float).ravel()
    K = L.size
    b = int(np.argmax(L))
    level = min(level, N)
    n0 = N - level
    if side != 'opt':
        # Pessimistic start, carried in QUEUE LENGTHS: all n0 customers at the
        # bottleneck (eq 13). Seeding a residence instead makes the assumed
        # population n0*Rb/(Z+Rtot) < n0 once Z > 0, so the pessimism is
        # diluted and the resulting Xlo stops being a bound (violated exact on
        # 14% of random delay models, worst 43%).
        Q = np.zeros(K)
        Q[b] = n0
        Rk = L * (1.0 + Q)
        for n in range(n0 + 1, N + 1):
            Rk = L * (1.0 + Q)
            Q = (n / (Z + np.sum(Rk))) * Rk
        return Rk
    Rk = np.ones(K) * max(n0 * L[b] - Z, np.sum(L)) / K
    if n0 == 0:
        Rk = np.zeros(K)
    for n in range(n0 + 1, N + 1):
        Rtot = np.sum(Rk)
        if n == 1 or (Z + Rtot) == 0:
            Rk = L.copy()
        else:
            Rk = L * (1.0 + (n - 1) * Rk / (Z + Rtot))
    return Rk


def pfqn_pbh(L, N, Z=0.0, level=1):
    """Performance Bound Hierarchy (Eager-Sevcik 1983), single-class.

    Returns (Xlo, Xhi, Qlo, Qhi). Level-`level` throughput/queue bounds;
    level 1 (Z=0) equals the BJB optimistic bound, and the bracket tightens
    to exact MVA as level -> N.
    """
    L = np.asarray(L, dtype=float).ravel()
    if Z is None:
        Z = 0.0
    if level is None:
        level = 1
    Lmax = np.max(L)
    Ro = _pbh_residence(L, N, Z, level, 'opt')
    Rp = _pbh_residence(L, N, Z, level, 'pess')
    RoC = max(np.sum(Ro), max(N * Lmax - Z, np.sum(L)))
    Xhi = min(1.0 / Lmax, N / (Z + RoC))
    Xlo = N / (Z + np.sum(Rp))
    Qlo = Xlo * Ro
    Qhi = Xhi * Rp
    return Xlo, Xhi, Qlo, Qhi


def pfqn_pbk(L, N, Z=0.0, k=1):
    """Iterative PB(k) proportional bounds (Eager-Sevcik / CMS08). Backed by
    the PBH recursion at level k. Returns (Xlo, Xhi)."""
    if Z is None:
        Z = 0.0
    if k is None:
        k = 1
    Xlo, Xhi, _, _ = pfqn_pbh(L, N, Z, k)
    return Xlo, Xhi


def pfqn_bjbk(L, N, Z=0.0, k=1):
    """Iterative BJB(k) balanced job bounds (CMS08). BJB(1) recovers the
    noniterative balanced job bound. Returns (Xlo, Xhi)."""
    if Z is None:
        Z = 0.0
    if k is None:
        k = 1
    Xlo, Xhi, _, _ = pfqn_pbh(L, N, Z, k)
    return Xlo, Xhi


# ----------------------------- CBH (Dowdy et al.) --------------------------

def _cbh_hier(L, N, Z, c, side):
    L = np.asarray(L, dtype=float).ravel()
    M = L.size
    Rc = np.sum(L[:c])
    Lbc = np.max(L[:c])
    Lac = np.mean(L[:c])
    e = np.zeros(N + 1)
    e[0] = 1.0
    for i in range(1, N + 1):
        if side == 'upper':
            B = i / (Rc + (i - 1) * Lac)
        else:
            B = i / (Rc + (i - 1) * Lbc)
        e[i] = e[i - 1] / B
    if c == 1:
        e = L[0] ** np.arange(N + 1)
    g = e.copy()
    for m in range(c, M):  # servers c+1..M (0-based c..M-1)
        for n in range(1, N + 1):
            g[n] = g[n] + L[m] * g[n - 1]
    if Z > 0:
        # Poisson weight Z^j/j! through logs: the naive ratio overflows for
        # j >~ 171 in double, and j runs to the POPULATION here.
        _j = np.arange(N + 1, dtype=float)
        gd = np.exp(_j * np.log(Z) - gammaln(_j + 1.0))
        gfull = np.zeros(N + 1)
        for n in range(N + 1):
            gfull[n] = sum(g[j] * gd[n - j] for j in range(n + 1))
        g = gfull
    return g[N - 1] / g[N]


def pfqn_cbh(L, N, Z=0.0, level=2):
    """Convolutional Bound Hierarchy (Dowdy et al. 1984), single-class.

    `level` exactly-convolved servers (1..M); the bracket tightens
    monotonically and equals exact at level M. Returns (Xlo, Xhi).
    """
    L = np.asarray(L, dtype=float).ravel()
    if Z is None:
        Z = 0.0
    if level is None:
        level = 2
    M = L.size
    level = max(1, min(level, M))
    c = max(1, M - level)
    Xlo = _cbh_hier(L, N, Z, c, 'lower')
    Xhi = _cbh_hier(L, N, Z, c, 'upper')
    return Xlo, Xhi


# ----------------------------- MCUB (Kerola) -------------------------------

def pfqn_mcub(L, N, Z=None):
    """Multiclass Composite Upper Bound (Kerola 1986). L is M x R.

    Returns (Xub, Xlb): Xub the per-class composite UPPER bound (eqs 13-16),
    Xlb the per-class multiclass Balanced Job Bounds LOWER bound (eq 10) that
    seeds it. Both are 1 x R arrays.
    """
    L = np.atleast_2d(np.asarray(L, dtype=float))
    M, R = L.shape
    N = np.asarray(N, dtype=float).ravel()
    if Z is None:
        Z = np.zeros(R)
    Z = np.asarray(Z, dtype=float).ravel()
    Ntot = np.sum(N)
    R0 = np.sum(L, axis=0)
    Lb = np.max(L, axis=0)
    Xlb = N / (R0 + Z + (Ntot - 1) * Lb)
    Xub = np.zeros(R)
    for r in range(R):
        Uoth = np.zeros(M)
        for s in range(R):
            if s != r:
                Uoth = Uoth + Xlb[s] * L[:, s]
        Ucub = 1.0 - Uoth
        dev = np.full(M, np.inf)
        for k in range(M):
            if L[k, r] > 0:
                dev[k] = Ucub[k] / L[k, r]
        Xub[r] = np.min(dev)
    return Xub, Xlb


# ----------------------------- SSD (Suri-Dallery) --------------------------

def pfqn_ssd(L, N, Z=0.0, nservers=None):
    """Server-Station Disaggregation bounds (Suri-Dallery 1986, Thm 5),
    single-class multiserver. Returns (Xlo, Xhi).

    With Z>0 the queueing terms carry the terminal-workload correction of
    Lazowska et al. 1984, Table 5.2; adding Z without it is not a bound.
    """
    L = np.asarray(L, dtype=float).ravel()
    K = L.size
    if Z is None:
        Z = 0.0
    if nservers is None:
        C = np.ones(K)
    else:
        C = np.asarray(nservers, dtype=float).ravel()
        if C.size == 1:
            C = C * np.ones(K)
    Rl = np.sum(L)
    Yl = np.max(L / C)
    Ru = np.sum(L / C)
    Yu = Ru / K
    b = int(np.argmax(L / C))
    Xlo = N / (Rl + Z + (N - 1) * Yl / (1 + Z / (N * Rl)))
    Xhi = min(N / (Ru + Z + (N - 1) * Yu / (1 + Z / Ru)), C[b] / L[b], N / (Rl + Z))
    return Xlo, Xhi


# ----------------------------- SIB (Srinivasan) ----------------------------

def pfqn_sib(L, N, Z=0.0, level=3):
    """Successively Improving Bounds (Srinivasan 1985), single-class, Z=0 only.

    Returns (Xlo, Xhi, Wlo, Whi). Raises ValueError for Z>0 (delay needs the
    Section-3.2 demand substitution, not yet implemented).
    """
    L = np.asarray(L, dtype=float).ravel()
    if Z is None:
        Z = 0.0
    Z = float(np.sum(Z))
    if level is None:
        level = 3
    level = max(1, int(round(level)))
    if Z > 0:
        raise ValueError("pfqn_sib supports Z=0 only (delay needs the "
                         "Section-3.2 demand substitution, not yet implemented).")
    Lsum = np.sum(L)
    rho = L / Lsum
    rho_u = np.max(rho)
    imax = level + 3
    S = np.array([np.sum(rho ** i) for i in range(1, imax + 1)])  # S[i-1]=S_i
    S2 = S[1]

    def Sv(i):  # 1-indexed access S_i
        return S[i - 1]

    alpha = np.zeros(level + 1)  # alpha[k] = alpha_k
    alpha[0] = S2
    for i in range(1, level + 1):
        acc = 0.0
        for j in range(0, i):
            acc += Sv(i + 1 - j) * alpha[j]
        alpha[i] = Sv(i + 2) - acc

    def phi_u1(K):
        if K <= 0:
            return 0.0
        if K == 1:
            return S2
        eta = (K - 1) / K
        T1 = (K - 1) * rho_u - 1
        return 0.5 / eta * (T1 + np.sqrt(T1 ** 2 + 4 * (K - 1) * S2))

    def sigma(NN, i):
        s = 0.0
        if i <= 0:
            return 0.0
        Dbar = 1 + phi_u1(NN - 2)
        pnum = 1.0
        for j in range(1, i + 1):
            pnum *= (NN - 1 - (j - 1))
            s += (rho_u * Sv(j + 1) - Sv(j + 2)) * pnum / Dbar ** j
        return s

    def betaL(NN, i):
        b = 0.0
        if i <= 0:
            return 0.0
        for j in range(1, i):
            p = 1.0
            for m in range(2, j + 1):
                p *= (NN - m) / (1 + phi_u1(NN - m))
            b += alpha[j] * p
        p = 1.0
        for m in range(2, i + 1):
            p *= (NN - m) / (1 + phi_u1(NN - m))
        Nim2 = NN - 1 - i - 1
        corr = 1 + (alpha[i] / alpha[i - 1]) * Nim2 / (1 + Nim2 * alpha[0])
        b += alpha[i] * p * corr
        return b

    NN = N - 1
    phi_lo = (N - 1) * S2
    T1s2 = (N - 1) * rho_u - 1
    phi_hi = 0.5 * (T1s2 + np.sqrt(T1s2 ** 2 + 4 * (N - 1) * S2))

    if N >= 3:
        eta = (N - 2) / (N - 1)
        T1u = (N - 2) * rho_u - 1
        su = sigma(NN, level - 1)
        phi_u_n = 0.5 / eta * (T1u + np.sqrt(max(0.0, T1u ** 2 + 4 * (N - 2) * (S2 - su))))
        phi_hi = min(phi_hi, phi_u_n)
        T1l = (N - 2) * S2 - 1
        bl = betaL(NN, level - 1)
        # eq (3.23) divides by 2*eta, the SAME constant eq (3.22) applies as
        # 0.5/eta above; the Greek eta on the scan was read as the level index n
        phi_l_n = (T1l + np.sqrt(max(0.0, T1l ** 2 + 4 * (N - 2) * (S2 + (N - 2) * bl)))) / (2 * eta)
        phi_lo = max(phi_lo, phi_l_n)

    phi_lo = max(0.0, phi_lo)
    if phi_hi < phi_lo:
        phi_hi = phi_lo
    Wlo = Lsum * (1 + phi_lo) + Z
    Whi = Lsum * (1 + phi_hi) + Z
    Xlo = N / Whi
    Xhi = N / Wlo
    return Xlo, Xhi, Wlo, Whi


# ----------------------------- LD-BCMP (Anselmi-Cremonesi) -----------------

def pfqn_ldbcmp(L, N, Z=0.0, c=None, tol=1e-10):
    """Anselmi-Cremonesi (2008) lower throughput bound for closed single-class
    BCMP networks with load-dependent stations. Returns (Xlo, Rhi, Qhat).

    c[i]=0 marks a fixed-rate (LI) station; c[i]>0 a Heffes LD station with open
    queue (c[i]+1)*rho/(1-rho). Bottleneck assumed fixed-rate. NaN if N < Qhat.
    """
    L = np.asarray(L, dtype=float).ravel()
    M = L.size
    if Z is None:
        Z = 0.0
    Z = float(np.sum(Z))
    if c is None:
        c = np.zeros(M)
    c = np.asarray(c, dtype=float).ravel()
    Dstar = L
    Dm = np.max(Dstar)
    isbott = np.abs(Dstar - Dm) <= 1e-12 * Dm
    bmax = int(np.sum(isbott))
    lam = 1.0 / Dm
    Qhat = 0.0
    for i in range(M):
        if isbott[i]:
            continue
        rho_i = lam * Dstar[i]
        if rho_i >= 1:
            return np.nan, np.nan, np.nan
        Qhat += (c[i] + 1) * rho_i / (1 - rho_i)
    Qhat += lam * Z
    if N - Qhat < 0:
        return np.nan, np.nan, Qhat
    a = N - Qhat
    Xprime = 0.0
    Xlo = 0.0
    for _ in range(10000):
        Xprev = Xlo
        denom = Dm * (bmax + N - Qhat) - bmax * (Dm * Xprime) ** N * Dm
        Xlo = a / denom
        Xprime = Xlo
        if Xprev > 0 and abs(Xprev - Xlo) / Xprev <= tol:
            break
    Rhi = N / Xlo
    return Xlo, Rhi, Qhat


# ------------------- SCB (Dowdy-Carlson-Krantz-Tripathi) -------------------

def pfqn_scb(L, N):
    """Single-class bounds of multiclass networks (Dowdy et al. 1992, JACM 39(1)).

    Returns (Xlo, Xhi, Ulo, Uhi), a bracket on the total throughput and on the
    per-device utilizations of the UNKNOWN multiclass system whose single-class
    counterpart has demand vector L at population N.

    SEMANTICS DIFFER FROM EVERY OTHER pfqn_* BOUND. aba/bjb/gb/... bracket the
    exact solution OF THE GIVEN MODEL; this brackets the multiclass system that
    the given single-class model aggregates. The lower side is therefore the
    EXACT single-class solution, not an approximation of it.

    Theorem 2 / Corollary 2: aggregating an R-class model into its single-class
    counterpart can only understate performance, U_k,1 <= U_k,R and X_1 <= X_R,
    and Corollary 1 makes the utilization ratio uniform, U_k,R/U_k,1 = X_R/X_1
    for every k. Theorem 3 (their Expression 3) caps the relative throughput
    error at (m-1)/(N+m-1), m = min(N,K), independently of the demands. The
    single-server capacity U_k,R <= 1 caps the same ratio at 1/(X_1*max(L)),
    tight on the paper's own worst case, so both are applied. L holds queueing
    stations only: Theorem 3 rests on the delay-free balanced-network
    throughput, so a delay station is not admitted.
    """
    L = np.asarray(L, dtype=float).ravel()
    K = L.size
    if K == 0:
        raise ValueError("pfqn_scb requires at least one queueing station.")
    N = int(round(N))
    if N < 1:
        raise ValueError("pfqn_scb requires N >= 1.")
    # Exact single-class MVA at Z=0. This IS the lower bound (Theorem 2), so it
    # is computed exactly rather than bounded: a bounded X1 would not bracket X_R.
    Q = np.zeros(K)
    X1 = 0.0
    for n in range(1, N + 1):
        Rk = L * (1.0 + Q)
        X1 = n / np.sum(Rk)
        Q = X1 * Rk
    U1 = X1 * L
    m = min(N, K)
    ratio = (N + m - 1) / float(N)          # Theorem 3, Expression (3)
    Dmax = float(np.max(L))
    if X1 * Dmax > 0:
        # U_k,R <= 1 with the uniform ratio of Corollary 1. Tight at the worst case.
        ratio = min(ratio, 1.0 / (X1 * Dmax))
    return float(X1), float(X1 * ratio), U1, U1 * ratio


def pfqn_scbgap(N, K, r=None, undominated=False):
    """Maximum relative throughput error of merging r of N classes (Dowdy 1992).

    Demand-free bound on the relative throughput error incurred when r of the N
    single-customer classes of a closed product-form network are merged into
    one class. With r = N (the default) this is the full single-class
    aggregation error of their Theorem 3, at most 50%; with r < N it is the
    partial-aggregation error of their Theorem 4. The bound never reads the
    demands, so it can be attached as a certified error bar to any result
    computed on merged chains.

    General case, dominating classes allowed (Expression 4, and with r = N
    Expression 3): e = (min(r,K)-1)/(r+min(r,K)-1). Undominated case, every
    customer placing the same total demand (Theorem 5 and its comment (3),
    which lifts the N = R restriction): e = r(r-1)/(min(N,K)(2r-1)), valid for
    r <= K only, smaller than the general case by the factor r/min(N,K) and
    equal to it at r = K. THE DOMAIN IS NOT COSMETIC: Theorem 5 gives each of
    its R classes a dedicated device, so r never exceeds K there, and comment
    (3) states the generalization for r < K. Evaluated at r > K the expression
    climbs past the general bound and past the 50% cap of Theorem 3, i.e. it
    stops being a bound, so r > K is refused rather than returned.
    """
    N = int(round(N)); K = int(round(K))
    if r is None:
        r = N
    r = int(round(r))
    if N < 1 or K < 1:
        raise ValueError("pfqn_scbgap requires N >= 1 and K >= 1.")
    if r < 1 or r > N:
        raise ValueError("pfqn_scbgap requires 1 <= r <= N (r=%d, N=%d)." % (r, N))
    if r == 1:
        return 0.0                          # merging one class changes nothing
    if undominated:
        if r > K:
            raise ValueError(
                "The undominated (Theorem 5) form is defined for r <= K only "
                "(r=%d, K=%d); beyond it the expression exceeds the general "
                "bound and the 50%% cap." % (r, K))
        return r * (r - 1) / float(min(N, K) * (2 * r - 1))   # Thm 5, comment (3)
    m = min(r, K)
    return (m - 1) / float(r + m - 1)       # Expression (4); r=N gives (3)


def pfqn_usumbound(R, K, N):
    """Upper bound on sum_k U_k in a closed R-class network (Dowdy 1992, Thm 6).

    sum_k U_k,R <= (H-1) + (K-H+1)(N-H+1)/(K+N-2H+1), H = min(R,K). Demand-free
    and nondecreasing in R, which is what makes it invertible into a lower
    bound on the number of necessary classes; see pfqn_minclasses. At
    R >= min(N,K) it reaches min(N,K), the trivial one-busy-server-per-device
    cap. The paper's worked case is K = 2, N = 3, R = 1, giving 2N/(N+1) = 1.5.
    """
    R = int(round(R)); K = int(round(K)); N = int(round(N))
    if N < 1 or K < 1:
        raise ValueError("pfqn_usumbound requires N >= 1 and K >= 1.")
    if R < 1 or R > N:
        raise ValueError("pfqn_usumbound requires 1 <= R <= N (R=%d, N=%d)." % (R, N))
    H = min(R, K)
    return (H - 1) + (K - H + 1) * (N - H + 1) / float(K + N - 2 * H + 1)


def pfqn_minclasses(Usum, K, N):
    """Lower bound on the class count from a measured utilization sum (Dowdy 1992).

    Smallest number of customer classes R consistent with an observed sum of
    device utilizations, obtained by inverting the demand-free Expression (6)
    bound of pfqn_usumbound, which is nondecreasing in R. Only measured
    quantities are needed -- the utilizations, the device count and the
    population -- so the answer is available BEFORE any class-specific demand
    has been characterized. An upper bound on R is meaningless (extra classes
    can always be introduced by splitting) and none is returned.

    Returns NaN when Usum exceeds min(N,K) and so is unattainable by ANY class
    structure, which signals a measurement error rather than a workload needing
    more classes. The paper's example: K = 2, N = 3, Usum = 1.6 -> 2, since a
    single class admits at most 2N/(N+1) = 1.5.
    """
    K = int(round(K)); N = int(round(N))
    if N < 1 or K < 1:
        raise ValueError("pfqn_minclasses requires N >= 1 and K >= 1.")
    if Usum < 0:
        raise ValueError("pfqn_minclasses requires a nonnegative utilization sum.")
    tol = 1e-12 * max(1.0, abs(Usum))
    for R in range(1, N + 1):
        if pfqn_usumbound(R, K, N) >= Usum - tol:
            return float(R)
    return float('nan')


def pfqn_looping(L, N, Z=None, tol: float = 1e-6, maxiter: int = 1000
                 ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    """Eager Looping approximate MVA bounds.

    D. L. Eager, "Bounding Algorithms for Queueing Network Models of Computer
    Systems", Ph.D. thesis, Tech. Rept. CSRG-156, University of Toronto, 1984.
    Looping supplies the initial pessimistic and optimistic estimates that the
    multiple-class performance bound hierarchy starts from, so it carries a
    pair of bounds rather than a single fixed point. It is built on the
    convolution identity of Zahorjan (1980)

        Q_jk(N - 1_c) = [X_j^{+k}(N - 1_c) / X_j(N)] Q_jk(N),

    with X_j^{+k}(N - 1_c) estimated from the level-0 multiple-class PBH upper
    bound B_j and X_j(N) from the optimistic response time R_j^(opt). A HEAP
    H_j is the class-j congestion the current queue-length lower bounds have
    not accounted for; it is charged back at the pessimistic inflation factor
    V_c = max_k D_ck or the optimistic one L_c = min_k D_ck. The level-0
    multiple-class PBH bounds on the mean response time are

        J_j(n) = sum_k D_jk,  B_j(n) = sum_k D_jk + (sum(n) - 1) max_k D_jk,

    i.e. an arriving customer queues behind nobody, respectively behind every
    other customer in the network at its own worst centre.

    Returns (Xlo, Xup, QN, RN, it) with Xlo the pessimistic and Xup the
    optimistic throughput bound.
    """
    L, N, Z, M, R = _amva_prep(L, N, Z)

    Dtot = L.sum(axis=0)
    Vpess = L.max(axis=0) if M > 0 else np.zeros(R)
    Lopt = L.min(axis=0) if M > 0 else np.zeros(R)
    Ntot = float(N.sum())
    Jbnd = Dtot
    Bm = Dtot + max(Ntot - 2, 0) * Vpess
    beta = np.eye(R)

    Qm = np.zeros((M, R, R))            # Qm[k, j, c] = Q_jk(N - 1_c)
    for c in range(R):
        for j in range(R):
            Qm[:, j, c] = max(N[j] - beta[c, j], 0) / M if M > 0 else 0.0
    Hopt = np.zeros((R, R))
    Hpess = np.zeros((R, R))

    QN = np.zeros((M, R))
    RN = np.zeros((M, R))
    XN = np.zeros(R)
    Rc = np.zeros(R)
    Rpess = np.zeros(R)
    Ropt = np.zeros(R)
    it = 1
    for it in range(1, maxiter + 1):
        QN_old = QN.copy()
        for c in range(R):
            if N[c] == 0:
                RN[:, c] = 0.0
                XN[c] = 0.0
                Rc[c] = 0.0
                Rpess[c] = 0.0
                continue
            Qk = Qm[:, :, c].sum(axis=1)
            RN[:, c] = L[:, c] * (1 + Qk)
            Rc[c] = RN[:, c].sum()
            Rpess[c] = Rc[c] + Vpess[c] * Hpess[:, c].sum()
            XN[c] = N[c] / (Z[c] + Rpess[c])
        for c in range(R):
            if N[c] == 0:
                Ropt[c] = 0.0
                continue
            sat = -np.inf
            for ist in range(M):
                den = 1 - (float(np.sum(XN * L[ist, :])) - XN[c] * L[ist, c])
                if den > 0:
                    sat = max(sat, L[ist, c] * N[c] / den - Z[c])
            heaped = Rc[c] + Lopt[c] * Hopt[:, c].sum()
            # an optimistic bound can never exceed the pessimistic one
            Ropt[c] = min(max(sat, heaped, Dtot[c]), Rpess[c])
        for c in range(R):
            QN[:, c] = XN[c] * RN[:, c]
        for c in range(R):
            for j in range(R):
                nj = N[j] - beta[c, j]
                if N[j] <= 0 or nj <= 0:
                    Qm[:, j, c] = 0.0
                else:
                    Qm[:, j, c] = (nj / N[j]) * ((Z[j] + Ropt[j]) / (Z[j] + Bm[j])) * QN[:, j]
                qsum = float(Qm[:, j, c].sum())
                if nj > 0:
                    Hopt[j, c] = max(0.0, Jbnd[j] / (Z[j] + Jbnd[j]) * nj - qsum)
                    Hpess[j, c] = max(0.0, Bm[j] / (Z[j] + Bm[j]) * nj - qsum)
                else:
                    Hopt[j, c] = 0.0
                    Hpess[j, c] = 0.0
        nz = N > 0
        if not np.any(nz) or (it > 1 and np.max(np.abs(QN[:, nz] - QN_old[:, nz])) < tol):
            break

    Xup = np.zeros(R)
    for c in range(R):
        if N[c] > 0:
            Xup[c] = N[c] / (Z[c] + Ropt[c])
    return XN.reshape(1, -1), Xup.reshape(1, -1), QN, RN, it


__all__ = [
    'pfqn_looping',
    'pfqn_pbh', 'pfqn_pbk', 'pfqn_bjbk', 'pfqn_cbh',
    'pfqn_mcub', 'pfqn_ssd', 'pfqn_sib', 'pfqn_ldbcmp',
    'pfqn_scb', 'pfqn_scbgap', 'pfqn_usumbound', 'pfqn_minclasses',
]
