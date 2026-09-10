"""
Asymptotic Methods for Normalizing Constants.

Implements asymptotic approximation methods for computing normalizing
constants in closed product-form queueing networks, including:
- Logistic Expansion (LE)
- Cubature Methods (Grundmann-Moeller rules)
"""

import numpy as np
from typing import Tuple, Optional
from scipy.special import gammaln

from ...constants import GlobalConstants

_TINY = np.finfo(float).tiny

def _factln(n) -> np.ndarray:
    """Log factorial using gamma function."""
    n = np.asarray(n, dtype=float)
    return gammaln(1 + n)

def _multinomialln(n: np.ndarray) -> float:
    """Log multinomial coefficient."""
    n = np.asarray(n, dtype=float)
    return float(gammaln(1 + np.sum(n)) - np.sum(gammaln(1 + n)))

def _allbut(y: np.ndarray, idx: int) -> np.ndarray:
    """Return array without element at index idx."""
    return np.delete(y, idx)

def _pfqn_le_fpi(L: np.ndarray, N: np.ndarray) -> np.ndarray:
    """Fixed-point iteration to find mode location (no think time)."""
    M, R = L.shape

    u = np.ones(M) / M
    u_prev = np.full(M, np.inf)

    max_iter = 1000
    for _ in range(max_iter):
        if np.linalg.norm(u - u_prev, 1) < 1e-10:
            break
        u_prev = u.copy()

        for i in range(M):
            u[i] = 1.0 / (np.sum(N) + M)
            for r in range(R):
                denom = np.dot(u_prev, L[:, r])
                if denom > 0:
                    u[i] += N[r] / (np.sum(N) + M) * L[i, r] * u_prev[i] / denom

    return u

def _pfqn_le_fpiZ(
    L: np.ndarray, N: np.ndarray, Z: np.ndarray
) -> Tuple[np.ndarray, float]:
    """Fixed-point iteration to find mode location (with think time)."""
    M, R = L.shape

    eta = np.sum(N) + M
    u = np.ones(M) / M
    # Note: eq. (35) in the SIGMETRICS 2017 paper has a spurious +1 in the v
    # equation; the correct stationary point is v = eta - sum_r xi_r*Z_r.
    v = eta
    u_prev = np.full(M, np.inf)

    max_iter = 1000
    for _ in range(max_iter):
        if np.linalg.norm(u - u_prev, 1) < 1e-10:
            break
        u_prev = u.copy()
        v_prev = v

        for i in range(M):
            u[i] = 1.0 / eta
            for r in range(R):
                denom = Z[r] + v * np.dot(u_prev, L[:, r])
                if denom > 0:
                    u[i] += (N[r] / eta) * (Z[r] + v * L[i, r]) * u_prev[i] / denom

        # Compute xi and update v
        xi = np.zeros(R)
        for r in range(R):
            denom = Z[r] + v_prev * np.dot(u_prev, L[:, r])
            if denom > 0:
                xi[r] = N[r] / denom

        v = eta - np.sum(xi * Z)

    return u, v

def _pfqn_le_hessian(L: np.ndarray, N: np.ndarray, u: np.ndarray) -> np.ndarray:
    """Compute Hessian matrix (no think time case)."""
    M, R = L.shape

    Ntot = np.sum(N)
    hu = np.zeros((M - 1, M - 1))

    for i in range(M - 1):
        for j in range(M - 1):
            if i != j:
                hu[i, j] = -(Ntot + M) * u[i] * u[j]
                for r in range(R):
                    denom = np.dot(u, L[:, r]) ** 2
                    if denom > 0:
                        hu[i, j] += N[r] * L[i, r] * L[j, r] * u[i] * u[j] / denom
            else:
                u_others = _allbut(u, i)
                hu[i, j] = (Ntot + M) * u[i] * np.sum(u_others)
                for r in range(R):
                    L_others = _allbut(L[:, r], i)
                    denom = np.dot(u, L[:, r]) ** 2
                    if denom > 0:
                        hu[i, j] -= N[r] * L[i, r] * u[i] * np.dot(u_others, L_others) / denom

    return hu

def _pfqn_le_hessianZ(
    L: np.ndarray, N: np.ndarray, Z: np.ndarray, u: np.ndarray, v: float
) -> np.ndarray:
    """Compute Hessian matrix (with think time case)."""
    K, R = L.shape

    Ntot = np.sum(N)

    # Compute csi
    csi = np.zeros(R)
    for r in range(R):
        denom = Z[r] + v * np.dot(u, L[:, r])
        if denom > 0:
            csi[r] = N[r] / denom

    # Compute Lhat
    Lhat = np.zeros((K, R))
    for k in range(K):
        for r in range(R):
            Lhat[k, r] = Z[r] + v * L[k, r]

    eta = Ntot + K
    A = np.zeros((K, K))

    # Off-diagonal elements
    for i in range(K):
        for j in range(K):
            if i != j:
                A[i, j] = -eta * u[i] * u[j]
                for r in range(R):
                    if N[r] > 0:
                        A[i, j] += csi[r] ** 2 * Lhat[i, r] * Lhat[j, r] * u[i] * u[j] / N[r]

    # Diagonal elements
    for i in range(K):
        row_sum = np.sum(_allbut(A[i, :], i))
        A[i, i] = -row_sum

    # Reduce to (K-1) x (K-1)
    A_reduced = A[: K - 1, : K - 1]

    # Add extra element for v
    A_full = np.zeros((K, K))
    A_full[: K - 1, : K - 1] = A_reduced

    A_full[K - 1, K - 1] = 1.0
    for r in range(R):
        if N[r] > 0:
            A_full[K - 1, K - 1] -= (csi[r] ** 2 / N[r]) * Z[r] * np.dot(u, L[:, r])
    A_full[K - 1, K - 1] *= v

    for i in range(K - 1):
        val = 0.0
        for r in range(R):
            if N[r] > 0:
                val += v * u[i] * (
                    (csi[r] ** 2 / N[r]) * Lhat[i, r] * np.dot(u, L[:, r]) - csi[r] * L[i, r]
                )
        A_full[i, K - 1] = val
        A_full[K - 1, i] = val

    return A_full

def pfqn_le(
    L: np.ndarray, N: np.ndarray, Z: Optional[np.ndarray] = None
) -> Tuple[float, float]:
    """
    Logistic Expansion (LE) asymptotic approximation for normalizing constant.

    Provides an asymptotic estimate of the normalizing constant for closed
    product-form queueing networks. Useful for large populations where exact
    methods become computationally expensive.

    Args:
        L: Service demand matrix (M x R).
        N: Population vector (R,).
        Z: Think time vector (R,). Optional.

    Returns:
        Tuple of (Gn, lGn):
            Gn: Estimated normalizing constant.
            lGn: Logarithm of normalizing constant.

    Reference:
        G. Casale. "Accelerating performance inference over closed systems by
        asymptotic methods." ACM SIGMETRICS 2017.
    """
    L = np.atleast_2d(np.asarray(L, dtype=float))
    N = np.asarray(N, dtype=float).ravel()

    M, R = L.shape

    # Handle empty or trivial cases
    if L.size == 0 or N.size == 0 or np.sum(N) == 0 or np.sum(L) < 1e-4:
        # Per-class Z, as MATLAB's sum(Z,1) is, and an empty class contributes 0
        # rather than 0*log(0). Z may also be absent altogether on this branch.
        Zt = np.zeros(N.size) if Z is None else np.asarray(Z, dtype=float).ravel()
        lGn = -np.sum(_factln(N))
        for r in range(N.size):
            if N[r] > 0:
                lGn += N[r] * np.log(Zt[r] if r < Zt.size else 0.0)
        Gn = np.exp(lGn)
        return float(Gn), float(lGn)

    if Z is None or np.sum(Z) < GlobalConstants.Zero:
        # Case without think time
        umax = _pfqn_le_fpi(L, N)
        A = _pfqn_le_hessian(L, N, umax)

        S = 0.0
        for r in range(R):
            term = np.dot(umax, L[:, r])
            if term > 0:
                S += N[r] * np.log(term)

        det_A = np.linalg.det(A)
        if det_A <= 0:
            det_A = 1e-100  # Fallback for numerical issues

        # Cas17 eq.(34) as published; pfqn_ble adds the eps->0 bias correction.
        lGn = (
            _multinomialln(np.append(N, M - 1))
            + _factln(np.array([M - 1]))[0]
            + (M - 1) * np.log(np.sqrt(2 * np.pi))
            - np.log(np.sqrt(det_A))
            + np.sum(np.log(np.maximum(umax, 1e-100)))
            + S
        )
    else:
        # Case with think time
        Z = np.asarray(Z, dtype=float).ravel()
        umax, vmax = _pfqn_le_fpiZ(L, N, Z)
        A = _pfqn_le_hessianZ(L, N, Z, umax, vmax)

        S = 0.0
        for r in range(R):
            term = Z[r] + vmax * np.dot(umax, L[:, r])
            if term > 0:
                S += N[r] * np.log(term)

        det_A = np.linalg.det(A)
        if det_A <= 0:
            det_A = 1e-100

        lGn = (
            -np.sum(_factln(N))
            - vmax
            + M * np.log(max(vmax, 1e-100))
            + M * np.log(np.sqrt(2 * np.pi))
            - np.log(np.sqrt(det_A))
            + np.sum(np.log(np.maximum(umax, 1e-100)))
            + S
        )

    Gn = np.exp(lGn)
    return float(Gn), float(lGn)

def pfqn_ble(
    L: np.ndarray, N: np.ndarray, Z: np.ndarray = None
) -> Tuple[float, float]:
    """Logistic expansion with the eps->0 bias correction (BLE).

    Cas17 Theorem 4.1 holds for eps >= eps_N > 0, where the K(1+eps*N) self-looping
    populations make the integrand concentrate. Evaluated at eps->0, as pfqn_le does,
    the curvature at the saddle tends to 1 rather than growing with N, so Laplace's
    method has no asymptotic regime there and carries an O(1) relative bias of
    e/sqrt(2*pi) PER LAPLACED DIRECTION. The count is the exponent on sqrt(2*pi) in
    the branch taken: M-1 with Z=0, where the radial integral is exact as
    Gamma(N+M), and M with Z>0, where the radius is Laplaced too. Measured over the
    1562 models of the Cas17 dataset (Zenodo 546873, sec5.3.1, sigma=100) the Z>0
    deficit is M to within 0.01 units. The published expansion is NOT in error; the
    correction is EMPIRICAL and is not part of Cas17. See _kb/03-api-layer.md.

    Args:
        L: Service demand matrix (MxR).
        N: Population vector (1xR).
        Z: Think time vector (1xR), optional.

    Returns:
        (Gn, lGn), the normalizing constant and its logarithm.
    """
    L = np.atleast_2d(np.asarray(L, dtype=float))
    M = L.shape[0]
    N = np.atleast_1d(np.asarray(N, dtype=float))
    _, lGn = pfqn_le(L, N, Z)
    if L.size == 0 or N.size == 0 or np.sum(N) == 0 or np.sum(L) < 1e-4:
        # Degenerate branch: the delay term is exact, no Laplace step to correct.
        return float(np.exp(lGn)), float(lGn)
    # Same predicate pfqn_le branches on, so the count always matches the branch.
    no_delay = Z is None or np.sum(Z) < GlobalConstants.Zero
    n_gauss = M - 1 if no_delay else M
    lGn = lGn + n_gauss * (1.0 - np.log(2 * np.pi) / 2)
    return float(np.exp(lGn)), float(lGn)


def pfqn_lekt_route(L: np.ndarray, N: np.ndarray, Z: np.ndarray = None) -> str:
    """The side pfqn_lekt computes on: 'kt' when R <= M or a class self-loops
    (one nonzero demand and no think time, which pfqn_kt extracts exactly), 'le'
    otherwise. The KT side is an R-dimensional convex solve and an R x R
    determinant, the LE side an M-dimensional fixed point and an (M-1) x (M-1) one."""
    L = np.atleast_2d(np.asarray(L, dtype=float))
    M, R = L.shape
    Z = np.zeros(R) if Z is None else np.asarray(Z, dtype=float).flatten()
    selfloop = False
    if R > 1:
        selfloop = bool(np.any((np.count_nonzero(L, axis=0) == 1) & (Z == 0)))
    return 'kt' if (R <= M or selfloop) else 'le'


def pfqn_lekt(
    L: np.ndarray, N: np.ndarray, Z: np.ndarray = None
) -> Tuple[float, float]:
    """The common corrected asymptotic expansion (LE-KT), computed on the cheaper side.

    The corrected logistic expansion (pfqn_ble) and the corrected Knessl-Tier
    expansion (pfqn_bkt) are ONE estimator, evaluated in M-1 and in R dimensions.
    With a think time their stationary points are one point in dual coordinates,
    xi_r = N_r/(Z_r + v u'L_r) being the class throughputs of the LE fixed point and
    v u_k = 1/(1-U_k) the M/M/1 factor of the KT saddle, and Sylvester's identity
    exchanges the R x R Hessian determinant for the M x M one, after which every 2 pi
    cancels; they agree to the accuracy of the two saddle-point solvers (~1e-7 nats,
    1e-14 with polished saddles). Without a think time the LE branch integrates the
    radius exactly as Gamma(N+M) while KT Laplaces it, so the two differ by the
    constant (1-log(2 pi)/2) - r(N+M), r the Stirling remainder of a Gamma direction;
    the common estimator is defined as the KT value, and the LE side here carries
    M(1-log(2 pi)/2) - r(N+M) rather than pfqn_ble's (M-1)(1-log(2 pi)/2). The route
    is pfqn_lekt_route. See _kb/03-api-layer.md.

    Args:
        L: Service demand matrix (MxR).
        N: Population vector (1xR).
        Z: Think time vector (1xR), optional.

    Returns:
        (Gn, lGn), the normalizing constant and its logarithm.
    """
    L = np.atleast_2d(np.asarray(L, dtype=float))
    M, R = L.shape
    N = np.atleast_1d(np.asarray(N, dtype=float))
    Z = np.zeros(R) if Z is None else np.asarray(Z, dtype=float).flatten()
    if pfqn_lekt_route(L, N, Z) == 'kt':
        from .kt import pfqn_bkt
        return pfqn_bkt(L, N, Z)
    G, lGn = pfqn_ble(L, N, Z)
    if L.size == 0 or N.size == 0 or np.sum(N) == 0 or np.sum(L) < 1e-4:
        return G, lGn  # pfqn_ble's degenerate branch: the delay term is exact
    if np.sum(Z) >= GlobalConstants.Zero:
        return G, lGn
    # the Z = 0 branch of pfqn_ble counts M-1 directions; the common estimator
    # carries M kappa - r(N+M)
    eta = float(np.sum(N) + M)
    kappa = 1.0 - np.log(2 * np.pi) / 2
    r = gammaln(eta) - (eta - 0.5) * np.log(eta) + eta - 0.5 * np.log(2 * np.pi)
    lGn = float(lGn + kappa - r)
    return float(np.exp(lGn)), lGn


def _grnmol(
    f, V: np.ndarray, s: int, tol: float = 1e-8
) -> Tuple[np.ndarray, int]:
    """
    Grundmann-Moeller cubature rule for simplex integration.

    Reference:
        "Invariant Integration Formulas for the N-Simplex by Combinatorial Methods",
        A. Grundmann and H. M. Moller, SIAM J Numer. Anal. 15(1978), pp. 282-290.
    """
    n = V.shape[0]
    Q = np.zeros(s + 1)
    Qv = np.zeros(s + 1)
    import math
    Vol = 1.0 / math.factorial(n)
    nv = 0
    d = 0

    while True:
        m = n + 2 * d + 1
        al = np.ones(n, dtype=float)
        alz = 2 * d + 1
        Qs = 0.0

        while True:
            x = V @ np.append([alz], al) / m
            Qs += f(x)
            nv += 1

            for j in range(n):
                alz -= 2
                if alz > 0:
                    al[j] += 2
                    break
                alz += al[j] + 1
                al[j] = 1

            if alz == 2 * d + 1:
                break

        d += 1
        Qv[d - 1] = Vol * Qs

        Q[d - 1] = 0
        p = 2.0 / np.prod(np.arange(n + 1, m + 1) * 2.0)
        for i in range(1, d + 1):
            Q[d - 1] += ((m + 2 - 2 * i) ** (2 * d - 1)) * p * Qv[d - i]
            p = -p * (m + 1 - i) / i

        if d > s or (d > 1 and abs(Q[d - 1] - Q[d - 2]) < tol * abs(Q[d - 2])):
            return Q[:d], nv

    return Q, nv

# quadrature points of v in the think-time branch of pfqn_cub; matches
# pfqn_cub.m and Pfqn_cub.java, and pfqn_nc prices CUB against it
CUB_THINK_STEPS = 10000

# integrand-evaluation budget above which pfqn_nc prefers le over cub
CUB_MAX_EVALS = 10000000


def pfqn_cub_evals(M: int, order: int, Z: Optional[np.ndarray] = None,
                   atol: float = 1e-8) -> int:
    """Number of integrand evaluations pfqn_cub performs at this order."""
    from math import comb
    n = M - 1
    nodes = sum(comb(n + 2 * d, n) for d in range(order + 1))
    has_think = Z is not None and float(np.sum(Z)) >= atol
    return nodes * (CUB_THINK_STEPS if has_think else 1)


def pfqn_cub(
    L: np.ndarray,
    N: np.ndarray,
    Z: Optional[np.ndarray] = None,
    order: Optional[int] = None,
    atol: float = 1e-8,
) -> Tuple[float, float]:
    """
    Cubature method for normalizing constant using Grundmann-Moeller rules.

    Uses numerical integration over simplices to compute the normalizing
    constant exactly (for sufficient order) or approximately.

    Args:
        L: Service demand matrix (M x R).
        N: Population vector (R,).
        Z: Think time vector (R,). Optional.
        order: Degree of cubature rule (default: ceil((sum(N)-1)/2)).
        atol: Absolute tolerance (default: 1e-8).

    Returns:
        Tuple of (Gn, lGn):
            Gn: Estimated normalizing constant.
            lGn: Logarithm of normalizing constant.

    Reference:
        G. Casale. "Accelerating performance inference over closed systems by
        asymptotic methods." ACM SIGMETRICS 2017.
    """
    L = np.atleast_2d(np.asarray(L, dtype=float))
    N = np.asarray(N, dtype=float).ravel()

    M, R = L.shape

    # Handle empty or trivial cases
    if L.size == 0 or N.size == 0 or np.sum(N) == 0:
        return 1.0, 0.0

    if order is None:
        order = int(np.ceil((np.sum(N) - 1) / 2))

    if Z is None or np.sum(Z) < atol:
        # Case without think time
        Nt = np.sum(N)
        beta = N / Nt

        V = np.eye(M - 1, M)  # Simplex vertices

        def f(x):
            # Complete the simplex point
            x_full = np.append(x, 1 - np.sum(x))
            Lx = x_full @ L
            h = beta * np.log(np.maximum(Lx, 1e-100))
            return np.prod(np.exp(Nt * h))

        Q, _ = _grnmol(f, V, order, atol)
        if len(Q) > 0:
            Gn = Q[-1] * np.exp(
                gammaln(1 + np.sum(N) + M - 1) - np.sum(gammaln(1 + N))
            )
        else:
            Gn = 1.0
    else:
        # Case with think time - numerical integration over v
        Z = np.asarray(Z, dtype=float).ravel()
        steps = CUB_THINK_STEPS
        Nt = np.sum(N)
        beta = N / Nt
        Gn = 0.0
        vmax = Nt * 10
        dv = vmax / steps

        V = np.eye(M - 1, M)

        for v_val in np.arange(0, vmax + dv, dv):
            Lv = L * v_val + np.tile(Z, (M, 1))

            def f(x):
                x_full = np.append(x, 1 - np.sum(x))
                Lx = x_full @ Lv
                h = beta * np.log(np.maximum(Lx, 1e-100))
                return np.exp(np.sum(Nt * h))

            Q, _ = _grnmol(f, V, order, atol)
            if len(Q) > 0:
                dG = np.exp(-v_val) * (v_val ** (M - 1)) * Q[-1] * dv
                Gn += dG

                if v_val > 0 and Gn > 0 and dG / Gn < atol:
                    break

        Gn *= np.exp(-np.sum(_factln(N)))

    lGn = np.log(max(Gn, 1e-300))
    return float(Gn), float(lGn)

def _logmeanexp(x: np.ndarray) -> float:
    """
    Compute log(mean(exp(x))) in a numerically stable way.

    Uses the log-sum-exp trick to avoid overflow/underflow.
    """
    x = np.asarray(x, dtype=float).ravel()
    if len(x) == 0:
        return -np.inf
    x_max = np.max(x)
    if np.isinf(x_max):
        return x_max
    return x_max + np.log(np.mean(np.exp(x - x_max)))

def pfqn_mci(
    D: np.ndarray,
    N: np.ndarray,
    Z: Optional[np.ndarray] = None,
    I: int = 100000,
    variant: str = 'imci'
) -> Tuple[float, float, np.ndarray]:
    """
    Monte Carlo Integration (MCI) for normalizing constant estimation.

    Provides a Monte Carlo estimate of the normalizing constant for closed
    product-form queueing networks.

    Args:
        D: Service demand matrix (M x R).
        N: Population vector (R,).
        Z: Think time vector (R,). Optional, defaults to zeros.
        I: Number of samples (default: 100000).
        variant: MCI variant - 'mci', 'imci' (improved), 'amci', 'lhsmci' or
            'rm' (repairman). Default: 'imci'. 'amci' and 'lhsmci' use the
            'imci' tilt and differ only in how the uniforms are drawn.
            'amci' draws ANTITHETIC pairs (u, 1-u). This does NOT reliably
            reduce variance here: the tilted integrand is not monotone in the
            exponential draws (the tilt term -(1-gamma)V decreases while the
            N log(VD+Z) term increases), so the pair correlation is not
            systematically negative; measured variance ratios against 'imci'
            range from 0.54 to 1.6 across models. It is kept because it is the
            Ross-Wang construction, not because it is the better default.
            'lhsmci' stratifies each coordinate by Latin hypercube sampling,
            which IS reliably variance-reducing on the same models (ratios 0.0
            to 0.48, exact quadrature in the limit of one station) at
            O(I log I) extra cost.

    Returns:
        Tuple of (G, lG, lZ):
            G: Estimated normalizing constant.
            lG: Logarithm of normalizing constant.
            lZ: Individual random sample log values.

    Reference:
        Implementation based on MonteQueue methodology.
    """
    from .mva import pfqn_bs

    D = np.atleast_2d(np.asarray(D, dtype=float))
    N = np.asarray(N, dtype=float).ravel()

    M, R = D.shape

    if Z is None:
        Z = np.zeros(R)
    else:
        Z = np.asarray(Z, dtype=float).ravel()

    # Handle empty or trivial cases
    if D.size == 0 or np.sum(D) < 1e-4:
        lGn = -np.sum(_factln(N)) + np.sum(N * np.log(np.sum(Z)))
        G = np.exp(lGn)
        return float(G), float(lGn), np.array([])

    # Compute throughput estimate using balanced system
    if variant.lower() in ('imci', 'amci', 'lhsmci'):
        # Improved MCI tilt; the three differ only in how the uniforms are drawn
        XN, _, _, _, _ = pfqn_bs(D, N, Z)
        tput = XN.ravel()
        util = D @ tput
        gamma = np.maximum(0.01, 1 - util)
    elif variant.lower() == 'mci':
        # Original MCI
        XN, _, _, _, _ = pfqn_bs(D, N, Z)
        tput = XN.ravel()
        util = D @ tput
        gamma = np.zeros(M)
        for i in range(M):
            if util[i] > 0.9:
                gamma[i] = 1.0 / np.sqrt(max(N))
            else:
                gamma[i] = 1.0 - util[i]
    elif variant.lower() == 'rm':
        # Repairman problem
        tput = N / (np.sum(D, axis=0) + Z + np.max(D, axis=0) * (np.sum(N) - 1))
        util = D @ tput
        gamma = np.zeros(M)
        for i in range(M):
            if util[i] > 0.9:
                gamma[i] = 1.0 / np.sqrt(max(N))
            else:
                gamma[i] = 1.0 - util[i]
    else:
        raise ValueError(f"Unknown variant: {variant}. Use 'mci', 'imci', 'amci', "
                         f"'lhsmci', or 'rm'.")

    # Ensure gamma is positive
    gamma = np.maximum(gamma, 1e-6)

    # Compute log factorials
    logfact = np.array([np.sum(np.log(np.arange(1, int(N[r]) + 1))) if N[r] > 0 else 0.0
                        for r in range(R)])

    # Uniform sampling with importance sampling
    if variant.lower() == 'amci':
        # Antithetic pairs (u, 1-u), the Ross-Wang construction
        Ih = (I + 1) // 2
        U = np.random.rand(Ih, M)
        VL = np.log(np.vstack([U, 1.0 - U]))[:I, :]
    elif variant.lower() == 'lhsmci':
        # Latin hypercube: one sample per stratum in every coordinate, so no
        # region of the tilted density is over- or under-sampled by chance.
        U = np.empty((I, M))
        for m in range(M):
            U[:, m] = (np.random.permutation(I) + np.random.rand(I)) / I
        VL = np.log(U)
    else:
        VL = np.log(np.random.rand(I, M))
    V = (-1.0 / gamma).reshape(1, -1) * VL

    ZI = np.tile(Z, (I, 1))

    # Importance sampling formula
    # lZ = -(ones(1,M) - gamma) * V' - sum(log(gamma)) - sum(logfact) + N*log(V*D+ZI)'
    term1 = -np.sum((1 - gamma) * V, axis=1)  # Shape: (I,)
    term2 = -np.sum(np.log(gamma))
    term3 = -np.sum(logfact)
    VD_plus_ZI = V @ D + ZI  # Shape: (I, R)
    term4 = np.sum(N * np.log(np.maximum(VD_plus_ZI, 1e-300)), axis=1)  # Shape: (I,)

    lZ = term1 + term2 + term3 + term4

    # Compute log of mean
    lG = _logmeanexp(lZ)

    if np.isinf(lG):
        lG = np.max(lZ)

    G = np.exp(lG)

    return float(G), float(lG), lZ

def pfqn_grnmol(L: np.ndarray, N: np.ndarray) -> Tuple[float, float]:
    """
    Normalizing constant using Grundmann-Moeller quadrature.

    Computes the normalizing constant for closed product-form queueing
    networks using Grundmann-Moeller cubature rules on simplices.

    This is an exact method that uses polynomial quadrature to compute
    the normalizing constant integral representation.

    Args:
        L: Service demand matrix (M x R).
        N: Population vector (R,).

    Returns:
        Tuple of (G, lG):
            G: Normalizing constant.
            lG: Logarithm of normalizing constant.

    Reference:
        Grundmann, A. and Moller, H.M. "Invariant Integration Formulas for
        the N-Simplex by Combinatorial Methods", SIAM J Numer. Anal. 15 (1978),
        pp. 282-290.
    """
    L = np.atleast_2d(np.asarray(L, dtype=float))
    N = np.asarray(N, dtype=float).ravel()

    M, R = L.shape

    # Handle trivial cases
    if L.size == 0 or N.size == 0 or np.sum(N) == 0:
        return 1.0, 0.0

    G = 0.0
    S = int(np.ceil((np.sum(N) - 1) / 2))

    H = np.zeros(1 + S)
    c = np.zeros(1 + S)
    w = np.zeros(1 + S)

    for i in range(S + 1):
        c[i] = 2 * (S - i) + M
        # w(1+i) = 2^-(2*S) * (-1)^i * c(1+i)^(2*S+1) / factorial(i) / factorial(i+c(1+i))
        log_w = (-(2 * S) * np.log(2) +
                 (2 * S + 1) * np.log(max(c[i], 1e-300)) -
                 _factln(np.array([i]))[0] -
                 _factln(np.array([i + c[i]]))[0])
        w[i] = ((-1) ** i) * np.exp(log_w)

        # Enumerate simplex states
        s_iter = S - i
        if s_iter < 0:
            continue

        # Generate all partitions of s_iter into M parts
        from scipy.special import comb
        num_states = int(comb(s_iter + M - 1, M - 1))

        if num_states == 0:
            H[i] = 1.0  # Only one state: all zeros
        else:
            H[i] = 0.0
            bvec = np.zeros(M, dtype=int)
            bvec[0] = s_iter

            for _ in range(num_states):
                # Compute term: prod((((2*bvec+1)/c(i))*L).^N)
                scale = (2 * bvec + 1) / c[i]
                Lscaled = scale.reshape(-1, 1) * L  # (M x R)
                # Sum over stations for each class, then raise to power N
                Lsum = np.sum(Lscaled, axis=0)  # (R,)
                log_term = np.sum(N * np.log(np.maximum(Lsum, 1e-300)))
                H[i] += np.exp(log_term)

                # Next partition
                bvec = _next_partition(bvec, s_iter)
                if bvec is None:
                    break

        G += w[i] * H[i]

    # Multiply by factorial(sum(N)+M-1) / prod(factorial(N))
    lG_coeff = _factln(np.array([np.sum(N) + M - 1]))[0] - np.sum(_factln(N))
    G = G * np.exp(lG_coeff)

    lG = np.log(max(abs(G), 1e-300))
    if G < 0:
        G = 0.0
        lG = -np.inf

    return float(G), float(lG)

def _next_partition(bvec: np.ndarray, total: int) -> Optional[np.ndarray]:
    """
    Generate next partition of total into M non-negative integers.

    Uses reverse lexicographic order.

    Args:
        bvec: Current partition
        total: Sum constraint

    Returns:
        Next partition, or None if exhausted
    """
    M = len(bvec)
    if M <= 1:
        return None

    # Find rightmost position that can be decremented
    for i in range(M - 2, -1, -1):
        if bvec[i] > 0:
            bvec[i] -= 1
            # Compute remaining
            remaining = total - np.sum(bvec[:i+1])
            # Put remaining in position i+1
            bvec[i+1] = remaining
            # Zero out positions after i+1
            for j in range(i + 2, M):
                bvec[j] = 0
            return bvec

    return None

def pfqn_le_fpi(L: np.ndarray, N: np.ndarray) -> np.ndarray:
    """
    Fixed-point iteration to find mode location (no think time).

    Public wrapper for the internal _pfqn_le_fpi function.

    Args:
        L: Service demand matrix (M x R).
        N: Population vector (R,).

    Returns:
        Mode location vector u (M,).
    """
    return _pfqn_le_fpi(L, N)

def pfqn_le_fpiZ(
    L: np.ndarray, N: np.ndarray, Z: np.ndarray
) -> Tuple[np.ndarray, float]:
    """
    Fixed-point iteration to find mode location (with think time).

    Public wrapper for the internal _pfqn_le_fpiZ function.

    Args:
        L: Service demand matrix (M x R).
        N: Population vector (R,).
        Z: Think time vector (R,).

    Returns:
        Tuple of (u, v):
            u: Mode location vector (M,).
            v: Scale factor.
    """
    return _pfqn_le_fpiZ(L, N, Z)

def pfqn_le_hessian(L: np.ndarray, N: np.ndarray, u: np.ndarray) -> np.ndarray:
    """
    Compute Hessian matrix (no think time case).

    Public wrapper for the internal _pfqn_le_hessian function.

    Args:
        L: Service demand matrix (M x R).
        N: Population vector (R,).
        u: Mode location vector (M,).

    Returns:
        Hessian matrix (M-1 x M-1).
    """
    return _pfqn_le_hessian(L, N, u)

def pfqn_le_hessianZ(
    L: np.ndarray, N: np.ndarray, Z: np.ndarray, u: np.ndarray, v: float
) -> np.ndarray:
    """
    Compute Hessian matrix (with think time case).

    Public wrapper for the internal _pfqn_le_hessianZ function.

    Args:
        L: Service demand matrix (M x R).
        N: Population vector (R,).
        Z: Think time vector (R,).
        u: Mode location vector (M,).
        v: Scale factor.

    Returns:
        Hessian matrix (M x M).
    """
    return _pfqn_le_hessianZ(L, N, Z, u, v)

# ---------------------------------------------------------------------------
# Dirichlet closure and adaptive Gauss-Hermite quadrature of the McKenna-Mitra
# integral.  Both reuse the LE mode and Hessian and differ only in what is
# wrapped around them; see _kb/03-api-layer.md.
# ---------------------------------------------------------------------------

def _golub_welsch(diag_off: np.ndarray, mu0: float):
    """Nodes and weights of the Gauss rule whose Jacobi off-diagonal is given."""
    n = len(diag_off) + 1
    J = np.diag(diag_off, 1) + np.diag(diag_off, -1)
    lam, V = np.linalg.eigh(J)
    return lam, mu0 * V[0] ** 2

def _gausslegendre(n: int):
    """N-point Gauss-Legendre rule on [-1,1]."""
    k = np.arange(1, n, dtype=float)
    return _golub_welsch(k / np.sqrt(4 * k * k - 1), 2.0)

def _gausshermite(q: int):
    """Q-point Gauss-Hermite rule of the probabilists' weight exp(-z^2/2)."""
    k = np.arange(1, q, dtype=float)
    return _golub_welsch(np.sqrt(k), float(np.sqrt(2 * np.pi)))

_GL64 = _gausslegendre(64)

def _logdet(A: np.ndarray) -> float:
    """Log-determinant of a positive-definite matrix (0 for the empty matrix)."""
    if A.size == 0:
        return 0.0
    try:
        return float(2.0 * np.sum(np.log(np.diag(np.linalg.cholesky(A)))))
    except np.linalg.LinAlgError:
        sign, ld = np.linalg.slogdet(A)
        return float(ld)

def _radial_logf(t: float, c: np.ndarray, N: np.ndarray, Z: np.ndarray, M: int) -> float:
    """Log-integrand of the radial integral in t = log v, Jacobian included."""
    v = np.exp(t)
    return float(-v + M * t + N @ np.log(np.maximum(Z + v * c, _TINY)))

def _radial(c: np.ndarray, N: np.ndarray, Z: np.ndarray, M: int):
    """log J(c) and the moments of the tilted law of the radius.

    J(c) = int_0^inf exp(-v) v^(M-1) prod_r (Z_r + v c_r)^N_r dv, with
    G_r = E[T_r], vbar = E[v] and Lam = cov(T) - diag(E[T^2]/N) = grad^2_c log J,
    where T_r(v) = N_r v / (Z_r + v c_r).  Quadrature runs in t = log v, where the
    integrand is bounded at both ends, over two Gauss-Legendre panels meeting at
    the mode.
    """
    vg, wg = _GL64
    R = N.size
    t = np.log(N.sum() + M)
    for _ in range(200):
        v = np.exp(t)
        d = np.maximum(Z + v * c, _TINY)
        F1 = -v + M + float(N @ ((v * c) / d))
        F2 = -v + float(N @ ((v * c) * Z / d ** 2))
        if F2 > -1e-300:
            break
        step = max(min(-F1 / F2, 2.0), -2.0)
        if abs(step) < 1e-13:
            t += step
            break
        t += step
    v = np.exp(t)
    d = np.maximum(Z + v * c, _TINY)
    F2 = -v + float(N @ ((v * c) * Z / d ** 2))
    sig = 1.0 / np.sqrt(-F2) if F2 < -1e-300 else 1.0
    Fm = _radial_logf(t, c, N, Z, M)
    # Widen each half-window until the log-integrand has fallen by 60 nats, so the
    # discarded tails are below 1e-26 in relative terms.
    a = min(12.0 * sig, t + 745.0)
    for _ in range(60):
        if t - a <= -745.0 or _radial_logf(t - a, c, N, Z, M) < Fm - 60.0:
            break
        a = min(1.6 * a, t + 745.0)
    b = 12.0 * sig
    for _ in range(60):
        if _radial_logf(t + b, c, N, Z, M) < Fm - 60.0:
            break
        b *= 1.6
    tt = np.concatenate([0.5 * a * vg + (t - 0.5 * a), 0.5 * b * vg + (t + 0.5 * b)])
    W = np.concatenate([0.5 * a * wg, 0.5 * b * wg])
    vv = np.exp(tt)
    D = np.maximum(Z + np.outer(vv, c), _TINY)
    Fv = -vv + M * tt + np.log(D) @ N
    mx = Fv.max()
    e = W * np.exp(Fv - mx)
    se = e.sum()
    lJ = mx + np.log(se)
    p = e / se
    T = (vv[:, None] * N[None, :]) / D
    G = p @ T
    vbar = float(p @ vv)
    ET2 = T.T @ (p[:, None] * T)
    dg = np.zeros(R)
    nz = N > 0
    dg[nz] = np.diag(ET2)[nz] / N[nz]
    Lam = ET2 - np.outer(G, G) - np.diag(dg)
    return float(lJ), G, vbar, 0.5 * (Lam + Lam.T)

def _simplex_mode(L: np.ndarray, N: np.ndarray, Z: np.ndarray):
    """Mode and curvature of h(w) = log J(L'x(w)) + sum_i log x_i.

    The fixed point x = (1 + x*(L@G))/vbar is the Z>0 analogue of pfqn_le_fpi:
    integrating by parts gives sum_i x_i (L@G)_i = vbar - M, so the update is
    normalised by construction.  At Z=0 it reduces to pfqn_le_fpi.  The term in
    the second derivative of x(w) drops at the mode against sum_i x_i == 1.
    """
    M = L.shape[0]
    x = _pfqn_le_fpiZ(L, N, Z)[0]
    x_1 = np.full(M, np.inf)
    it = 0
    while np.abs(x - x_1).sum() > 1e-11 and it < 10000:
        x_1 = x
        _, G, vbar, _ = _radial(x_1 @ L, N, Z, M)
        x = (1.0 + x_1 * (L @ G)) / vbar
        x = x / x.sum()
        it += 1
    lJ, _, _, Lam = _radial(x @ L, N, Z, M)
    P = L @ Lam @ L.T - np.diag(1.0 / x ** 2)
    Jm = (np.diag(x) - np.outer(x, x))[:, : M - 1]
    A = -(Jm.T @ P @ Jm)
    A = 0.5 * (A + A.T)
    return x, A, _logdet(A), float(lJ + np.log(x).sum())

def _softmax_gauge(w: np.ndarray) -> np.ndarray:
    a = np.concatenate([np.asarray(w, dtype=float).ravel(), [0.0]])
    e = np.exp(a - a.max())
    return e / e.sum()

def pfqn_aghq(
    L: np.ndarray, N: np.ndarray, Z: Optional[np.ndarray] = None, q: int = 3
) -> Tuple[float, float]:
    """
    Adaptive Gauss-Hermite quadrature of the McKenna-Mitra integral.

    Rescaling the simplex integral by the LE mode and curvature, w = w* + A^-1/2 z,
    and applying the q-node probabilists' Gauss-Hermite rule in each of the M-1
    directions gives a convergent rule whose q=1 member is pfqn_le itself (single
    node at the mode, weight sqrt(2*pi)), to the tolerance of the shared fixed
    point.  Cost is q^(M-1) evaluations, which confines the method to small M.

    A tensor rule is not invariant to the choice of A^-1/2: any B with B B' =
    inv(A) is admissible and they place the nodes differently.  The principal-axis
    frame from the eigendecomposition is used, as in the reference results; where
    two curvatures are close to equal the frame is close to arbitrary and two
    valid rules can part company well above their own error, converging back
    together as q grows.  Do not compare across codebases node by node.

    With Z>0 the radius is integrated numerically and the rule is
    applied to the M-1 simplex directions, so every node costs one radial
    quadrature.  q=1 there is LE with an exact radius, NOT pfqn_le's own Z>0
    branch.

    Args:
        L: Service demand matrix (M x R).
        N: Population vector (R,).
        Z: Think time vector (R,). Optional.
        q: Nodes per simplex direction (default 3).

    Returns:
        Tuple of (Gn, lGn).

    Reference:
        J. McKenna, D. Mitra. "Integral representations and asymptotic expansions
        for closed Markovian queueing networks: normal usage." BSTJ 61(5), 1982.
        G. Casale. "Accelerating performance inference over closed systems by
        asymptotic methods." ACM SIGMETRICS 2017.
    """
    L = np.atleast_2d(np.asarray(L, dtype=float))
    N = np.asarray(N, dtype=float).ravel()
    M, R = L.shape
    Z = np.zeros(N.size) if Z is None else np.asarray(Z, dtype=float).ravel()
    q = 3 if q is None else int(q)

    if L.size == 0 or N.size == 0 or np.sum(N) == 0 or np.sum(L) < 1e-4:
        lGn = -np.sum(_factln(N))
        if np.sum(Z) > 0:
            lGn += float(np.sum(N[N > 0] * np.log(Z[N > 0])))
        return float(np.exp(lGn)), float(lGn)

    if np.sum(Z) < GlobalConstants.Zero:
        umax = _pfqn_le_fpi(L, N)
        A = _pfqn_le_hessian(L, N, umax)
        ld = _logdet(A)
        h0 = float(N @ np.log(umax @ L)) + float(np.log(umax).sum())
        w0 = np.log(umax[: M - 1] / umax[M - 1]) if M > 1 else np.zeros(0)
        hfun = lambda w: (float(N @ np.log(_softmax_gauge(w) @ L))
                          + float(np.log(_softmax_gauge(w)).sum()))
        lacc = _aghq_rule(hfun, w0, h0, A, q, M - 1)
        lGn = (_multinomialln(np.append(N, M - 1)) + _factln(np.array([M - 1]))[0]
               + h0 + lacc - 0.5 * ld)
    else:
        xmax, A, ld, h0 = _simplex_mode(L, N, Z)
        w0 = np.log(xmax[: M - 1] / xmax[M - 1]) if M > 1 else np.zeros(0)

        def hfun(w):
            x = _softmax_gauge(w)
            return _radial(x @ L, N, Z, M)[0] + float(np.log(x).sum())

        lacc = _aghq_rule(hfun, w0, h0, A, q, M - 1)
        lGn = -np.sum(_factln(N)) + h0 + lacc - 0.5 * ld
    return float(np.exp(lGn)), float(lGn)

def _aghq_rule(hfun, w0: np.ndarray, h0: float, A: np.ndarray, q: int, d: int) -> float:
    """Log of the tensor Gauss-Hermite sum, accumulated with a running maximum.

    The det(A)^-1/2 of the rule is applied by the caller.
    """
    if d == 0:
        return 0.0
    nodes = q ** d
    if nodes > 10 ** 7:
        raise ValueError(
            'pfqn_aghq: the tensor rule needs q^(M-1)=%d nodes; reduce q or use pfqn_le.'
            % nodes)
    lam, V = np.linalg.eigh(0.5 * (A + A.T))
    if lam.min() <= 0:
        return float('nan')
    B = V @ np.diag(1.0 / np.sqrt(lam))
    z, wt = _gausshermite(q)
    lwt = np.log(wt)
    idx = np.zeros(d, dtype=int)
    lmax = -np.inf
    s = 0.0
    for _ in range(nodes):
        zz = z[idx]
        lt = float(lwt[idx].sum()) + hfun(w0 + B @ zz) - h0 + 0.5 * float(zz @ zz)
        if lt > lmax:
            s = s * np.exp(lmax - lt) + 1.0
            lmax = lt
        else:
            s += np.exp(lt - lmax)
        for j in range(d - 1, -1, -1):
            idx[j] += 1
            if idx[j] < q:
                break
            idx[j] = 0
    return float(lmax + np.log(s))

__all__ = [
    'pfqn_le',
    'pfqn_aghq',
    'pfqn_cub',
    'pfqn_mci',
    'pfqn_grnmol',
    'pfqn_le_fpi',
    'pfqn_le_fpiZ',
    'pfqn_le_hessian',
    'pfqn_le_hessianZ',
]
