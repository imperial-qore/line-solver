"""
Hierarchical and multiserver/load-dependent bound methods for closed
product-form queueing networks. Native-Python ports (no JVM) of the MATLAB
pfqn_{pbh,cbh,pbk,bjbk,mcub,ssd,sib,ldbcmp} bound functions used by SolverBA.

References:
- Eager-Sevcik 1983 (PBH), Dowdy et al. 1984 (CBH), Casale-Muntz-Serazzi 2008
  (iterative PB(k)/BJB(k)), Kerola 1986 (multiclass composite upper bound),
  Suri-Dallery 1986 (multiserver disaggregation), Srinivasan 1985 (SIB),
  Anselmi-Cremonesi 2008 (LD-BCMP closed-open equivalence).
"""

import numpy as np
from math import factorial


# ----------------------------- PBH (Eager-Sevcik) --------------------------

def _pbh_residence(L, N, Z, level, side):
    """Per-station residence-time vector of the level-`level` PBH bound."""
    L = np.asarray(L, dtype=float).ravel()
    K = L.size
    b = int(np.argmax(L))
    level = min(level, N)
    n0 = N - level
    if side == 'opt':
        Rk = np.ones(K) * max(n0 * L[b] - Z, np.sum(L)) / K
    else:  # 'pess'
        Rk = np.zeros(K)
        Rk[b] = n0
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
        gd = np.array([Z ** j / factorial(j) for j in range(N + 1)])
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
    single-class multiserver. Returns (Xlo, Xhi)."""
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
    Xlo = N / (Rl + Z + (N - 1) * Yl)
    Xhi = min(N / (Ru + Z + (N - 1) * Yu), C[b] / L[b], N / (Rl + Z))
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
        phi_l_n = (T1l + np.sqrt(max(0.0, T1l ** 2 + 4 * (N - 2) * (S2 + (N - 2) * bl)))) / (2 * level)
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


__all__ = [
    'pfqn_pbh', 'pfqn_pbk', 'pfqn_bjbk', 'pfqn_cbh',
    'pfqn_mcub', 'pfqn_ssd', 'pfqn_sib', 'pfqn_ldbcmp',
]
