"""
Asymptotic Methods for Normalizing Constants.

Implements asymptotic approximation methods for computing normalizing
constants in closed product-form queueing networks, including:
- Laguerre Expansion (LE)
- Cubature Methods (Grundmann-Moeller rules)
"""

import numpy as np
from typing import Tuple, Optional
from scipy.special import gammaln


def _factln(n: np.ndarray) -> np.ndarray:
    """Log factorial using gamma function."""
    return gammaln(1 + np.asarray(n, dtype=float))


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
    v = eta + 1
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

        v = eta + 1 - np.sum(xi * Z)

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
    Laguerre Expansion (LE) asymptotic approximation for normalizing constant.

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
        if Z is not None:
            Z = np.asarray(Z, dtype=float).ravel()
            sum_Z = np.sum(Z)
            if sum_Z > 0:
                lGn = -np.sum(_factln(N)) + np.sum(N * np.log(sum_Z))
            else:
                lGn = -np.sum(_factln(N))
        else:
            lGn = -np.sum(_factln(N))
        Gn = np.exp(lGn)
        return float(Gn), float(lGn)

    if Z is None or np.sum(Z) < 1e-10:
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
        steps = 1000
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
        variant: MCI variant - 'mci', 'imci' (improved), or 'rm' (repairman).
            Default: 'imci'.

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
    if variant.lower() == 'imci':
        # Improved MCI
        XN, _, _, _, _, _, _ = pfqn_bs(D, N, Z)
        tput = XN.ravel()
        util = D @ tput
        gamma = np.maximum(0.01, 1 - util)
    elif variant.lower() == 'mci':
        # Original MCI
        XN, _, _, _, _, _, _ = pfqn_bs(D, N, Z)
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
        raise ValueError(f"Unknown variant: {variant}. Use 'mci', 'imci', or 'rm'.")

    # Ensure gamma is positive
    gamma = np.maximum(gamma, 1e-6)

    # Compute log factorials
    logfact = np.array([np.sum(np.log(np.arange(1, int(N[r]) + 1))) if N[r] > 0 else 0.0
                        for r in range(R)])

    # Uniform sampling with importance sampling
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


__all__ = [
    'pfqn_le',
    'pfqn_cub',
    'pfqn_mci',
    'pfqn_grnmol',
    'pfqn_le_fpi',
    'pfqn_le_fpiZ',
    'pfqn_le_hessian',
    'pfqn_le_hessianZ',
]
