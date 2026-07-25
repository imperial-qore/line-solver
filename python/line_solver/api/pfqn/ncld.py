"""
Load-dependent Normalizing Constant methods for Product-Form Queueing Networks.

Native Python implementations of methods for computing normalizing constants
in load-dependent queueing networks.

Key functions:
    pfqn_ncld: Main dispatcher for load-dependent NC computation
    pfqn_gld: Generic load-dependent NC
    pfqn_gldsingle: Single-class load-dependent NC
    pfqn_comomrm_ld: COMOM method for load-dependent repairman models

References:
    Casale, G., et al. "LINE: A unified library for queueing network modeling."
"""

import numpy as np
from math import log, exp, log1p, factorial, lgamma
from typing import Tuple, Dict, Optional, Any
from dataclasses import dataclass


# Fine tolerance for numerical comparisons
FINE_TOL = 1e-12
ZERO = 1e-10
NEG_INF = float('-inf')


def _logsumexp2(a: float, b: float) -> float:
    """Pairwise log-sum-exp, stable when either argument is -inf."""
    if a > b:
        if b == NEG_INF:
            return a
        return a + log1p(exp(b - a))
    if a == NEG_INF:
        return b
    return b + log1p(exp(a - b))


def _factln(n: float) -> float:
    """Compute log(n!) using log-gamma function."""
    if n <= 0:
        return 0.0
    return lgamma(n + 1)


def _factln_array(arr: np.ndarray) -> np.ndarray:
    """Compute log(n!) element-wise for an array."""
    from scipy.special import gammaln
    arr = np.asarray(arr, dtype=float)
    result = np.zeros_like(arr)
    mask = arr > 0
    result[mask] = gammaln(arr[mask] + 1)
    return result


@dataclass
class PfqnNcResult:
    """Result of normalizing constant computation."""
    G: float
    lG: float
    method: str = "default"


@dataclass
class PfqnComomrmLdResult:
    """Result of COMOM load-dependent computation."""
    G: float
    lG: float
    prob: np.ndarray


def pfqn_mushift(mu: np.ndarray, k: int) -> np.ndarray:
    """
    Shift a load-dependent scaling vector by one position.

    Used in recursive normalizing constant computations.

    Args:
        mu: Load-dependent scalings matrix (M x N)
        k: Row index to shift

    Returns:
        Shifted mu matrix (M x N-1)
    """
    mu = np.atleast_2d(np.asarray(mu, dtype=float))
    M, N = mu.shape

    if N <= 1:
        return np.zeros((M, 0))

    mushift = mu[:, :-1].copy()
    mushift[k, :] = mu[k, 1:]

    return mushift


def pfqn_gldsingle(L: np.ndarray, N: np.ndarray, mu: np.ndarray,
                   options: Optional[Dict[str, Any]] = None) -> PfqnNcResult:
    """
    Compute normalizing constant for single-class load-dependent model.

    Auxiliary function used by pfqn_gld to compute the normalizing constant
    in a single-class load-dependent model using dynamic programming.

    Args:
        L: Service demands at all stations (M x 1)
        N: Number of jobs (scalar or 1x1 array)
        mu: Load-dependent scaling factors (M x Ntot)
        options: Solver options (unused, for API compatibility)

    Returns:
        PfqnNcResult with G (normalizing constant) and lG (log)

    Raises:
        RuntimeError: If multiclass model is detected
    """
    L_arr = np.asarray(L)
    L_dtype = complex if np.iscomplexobj(L_arr) else float
    L = np.atleast_2d(np.asarray(L, dtype=L_dtype))
    N = np.asarray(N, dtype=float).flatten()
    mu = np.atleast_2d(np.asarray(mu, dtype=float))

    M = L.shape[0]
    R = L.shape[1]

    if R > 1:
        raise RuntimeError("pfqn_gldsingle: multiclass model detected. "
                          "pfqn_gldsingle is for single class models.")

    N_val = int(np.ceil(N[0]))

    if N_val <= 0:
        return PfqnNcResult(G=1.0, lG=0.0)

    # see _kb/03-api-layer.md for rationale
    use_log = (not np.iscomplexobj(L) and not np.iscomplexobj(mu)
               and bool(np.all(np.real(L) >= 0)) and bool(np.all(mu > 0)))

    if use_log:
        # see _kb/03-api-layer.md for rationale
        lg = {}
        # lg[(0, n, 1)] stays -inf for n>=1: no station can hold n>=1 jobs.
        for n in range(1, N_val + 1):
            lg[(0, n, 1)] = NEG_INF

        with np.errstate(divide='ignore'):
            lL_all = np.log(np.real(L[:, 0]))   # -inf where the demand is zero
            lmu_all = np.log(mu)               # +inf where the rate is infinite

        for m in range(1, M + 1):
            for tm in range(1, N_val + 2):
                lg[(m, 0, tm)] = 0.0           # log(1): zero jobs

            lL = float(lL_all[m - 1])
            for n in range(1, N_val + 1):
                for tm in range(1, N_val - n + 2):
                    a = lg.get((m - 1, n, 1), NEG_INF)
                    lg_curr = lg.get((m, n - 1, tm + 1), NEG_INF)

                    mu_idx = tm - 1  # 0-indexed
                    if mu_idx < mu.shape[1]:
                        lmu_val = float(lmu_all[m - 1, mu_idx])
                    else:
                        lmu_val = 0.0          # mu = 1
                    b = lL + lg_curr - lmu_val
                    lg[(m, n, tm)] = _logsumexp2(a, b)

        lG = lg.get((M, N_val, 1), NEG_INF)
        return PfqnNcResult(G=float(np.exp(lG)), lG=lG)

    # Use dictionary for sparse storage with tuple keys
    # g[(m, n, tm)] maps to the value (complex if L is complex)
    g = {}

    # Initialize boundary conditions: g(0, n, 1) = 0 for n=1:N
    for n in range(1, N_val + 1):
        g[(0, n, 1)] = 0.0

    for m in range(1, M + 1):
        # Initialize boundary conditions: g(m, 0, tm) = 1 for tm=1:(N+1)
        for tm in range(1, N_val + 2):
            g[(m, 0, tm)] = 1.0

        for n in range(1, N_val + 1):
            for tm in range(1, N_val - n + 2):
                g_prev = g.get((m - 1, n, 1), 0.0)
                g_curr = g.get((m, n - 1, tm + 1), 0.0)

                # Get mu value safely
                mu_idx = tm - 1  # 0-indexed
                if mu_idx < mu.shape[1]:
                    mu_val = mu[m - 1, mu_idx]
                else:
                    mu_val = 1.0

                # MATLAB divides by mu even when negative, so we should too
                if mu_val != 0:
                    g[(m, n, tm)] = g_prev + L[m - 1, 0] * g_curr / mu_val
                else:
                    g[(m, n, tm)] = g_prev

    G = g.get((M, N_val, 1), 0.0)
    if abs(G) > 0:
        # Use complex log to handle negative G values (MATLAB's log does this)
        # The caller should use np.real() if they need only the real part
        import cmath
        lG = cmath.log(G)
    else:
        lG = NEG_INF

    return PfqnNcResult(G=G, lG=lG)


def pfqn_gld(L: np.ndarray, N: np.ndarray, mu: np.ndarray,
             options: Optional[Dict[str, Any]] = None) -> PfqnNcResult:
    """
    Compute normalizing constant of a load-dependent closed queueing network.

    Uses the generalized convolution algorithm for computing normalizing
    constants in load-dependent closed queueing networks.

    Args:
        L: Service demands at all stations (M x R)
        N: Number of jobs for each class (1 x R)
        mu: Load-dependent scalings (M x Ntot)
        options: Solver options

    Returns:
        PfqnNcResult with G (normalizing constant) and lG (log)
    """
    from .nc import pfqn_nc

    L_arr = np.asarray(L)
    L_dtype = complex if np.iscomplexobj(L_arr) else float
    L = np.atleast_2d(np.asarray(L, dtype=L_dtype))
    N = np.asarray(N, dtype=float).flatten()
    mu = np.atleast_2d(np.asarray(mu, dtype=float)) if mu is not None else None

    if options is None:
        options = {'tol': 1e-6, 'method': 'default'}

    M = L.shape[0]
    R = L.shape[1]
    Ntot = int(np.ceil(np.sum(N)))

    # Validate dimensions
    if len(N) != R:
        # Dimension mismatch between L columns and N elements
        # This can happen with chain aggregation - use the minimum
        R = min(R, len(N))
        L = L[:, :R]

    # Handle single station case
    if M == 1:
        N_tmp = []
        L_tmp = []
        for i in range(R):
            if abs(L[0, i]) > FINE_TOL:
                N_tmp.append(N[i])
                L_tmp.append(np.log(L[0, i]))

        if len(N_tmp) == 0:
            return PfqnNcResult(G=1.0, lG=0.0)

        N_tmp = np.array(N_tmp)
        L_tmp = np.array(L_tmp)

        # Ensure mu has enough columns
        if mu is not None:
            if Ntot >= mu.shape[1]:
                mu_row = mu[0, :].copy()
            else:
                mu_row = mu[0, :Ntot].copy()
        else:
            mu_row = np.ones(Ntot)

        # Compute log of mu values
        with np.errstate(divide='ignore'):
            log_mu = np.log(mu_row)
            log_mu = np.where(np.isfinite(log_mu), log_mu, 0.0)

        lG = (_factln(np.sum(N_tmp)) - np.sum(_factln_array(N_tmp)) +
              np.dot(N_tmp, L_tmp) - np.sum(log_mu[:Ntot]))
        G = np.exp(lG) if np.real(lG) > -700 else 0.0

        return PfqnNcResult(G=G, lG=lG)

    # Handle single-class case
    if R == 1:
        return pfqn_gldsingle(L, N, mu, options)

    # Handle empty L
    if L.size == 0 or np.sum(L) < FINE_TOL:
        return PfqnNcResult(G=0.0, lG=NEG_INF)

    # Initialize mu if None
    if mu is None:
        mu = np.ones((M, Ntot))

    # Check if load-dependent
    is_load_dep = False
    is_inf_server = np.zeros(M, dtype=bool)

    for i in range(M):
        mu_row = mu[i, :Ntot] if Ntot <= mu.shape[1] else np.concatenate([mu[i, :], np.ones(Ntot - mu.shape[1])])

        # Check if delay station (mu = [1, 2, 3, ...])
        expected_delay = np.arange(1, Ntot + 1, dtype=float)
        if len(mu_row) >= Ntot:
            is_delay = np.allclose(mu_row[:Ntot], expected_delay[:Ntot], atol=FINE_TOL)
        else:
            is_delay = False

        # Check if single server (mu = [1, 1, 1, ...])
        is_single = np.allclose(mu_row, 1.0, atol=FINE_TOL)

        if is_single:
            is_inf_server[i] = False
        elif is_delay:
            is_inf_server[i] = True
        else:
            is_inf_server[i] = False
            is_load_dep = True

    # If not load-dependent, use standard NC
    if not is_load_dep:
        Lli = L[~is_inf_server, :] if np.any(~is_inf_server) else np.zeros((1, R))
        Zli = L[is_inf_server, :] if np.any(is_inf_server) else np.zeros((1, R))

        if Lli.size == 0 or Lli.shape[0] == 0:
            Lli = np.zeros((1, R))
        if Zli.size == 0 or Zli.shape[0] == 0:
            Zli = np.zeros((1, R))

        Z_sum = np.sum(Zli, axis=0)
        result = pfqn_nc(Lli, N, Z_sum, method='exact')
        return PfqnNcResult(G=result[0], lG=result[1])

    # Handle zero population
    if Ntot == 0 or (np.abs(np.max(N)) < FINE_TOL and np.abs(np.min(N)) < FINE_TOL):
        return PfqnNcResult(G=1.0, lG=0.0)

    # Single-class case
    if R == 1:
        result = pfqn_gldsingle(L, N, mu, options)
        return result

    # Recursive case: G_M(N) = G_{M-1}(N) + sum_r L[M-1,r]/mu[M-1,0] * G_M(N-e_r)
    G = pfqn_gld(L[:-1, :], N, mu[:-1, :], options).G

    for r in range(R):
        if N[r] > FINE_TOL:
            N_1 = N.copy()
            N_1[r] -= 1

            mu_shifted = pfqn_mushift(mu, M - 1)

            if mu[M - 1, 0] > 0:
                G += (L[M - 1, r] / mu[M - 1, 0]) * pfqn_gld(L, N_1, mu_shifted, options).G

    if np.iscomplex(G):
        lG = np.log(G) if abs(G) > 0 else NEG_INF
    else:
        G = float(np.real(G))
        lG = log(G) if G > 0 else NEG_INF

    return PfqnNcResult(G=G, lG=lG)


def pfqn_comomrm_ld(L: np.ndarray, N: np.ndarray, Z: np.ndarray,
                    mu: np.ndarray, options: Optional[Dict[str, Any]] = None
                    ) -> PfqnComomrmLdResult:
    """
    Run the COMOM normalizing constant method on a load-dependent repairman model.

    Implements the Class-Oriented Method of Moments (COMOM) for computing normalizing
    constants in load-dependent repairman queueing models.

    Args:
        L: Service demands at all stations (M x R)
        N: Number of jobs for each class (1 x R)
        Z: Think times for each class (1 x R)
        mu: Load-dependent scalings (M x Ntot)
        options: Solver options

    Returns:
        PfqnComomrmLdResult with G, lG, and marginal probabilities
    """
    from .nc import pfqn_ca

    L = np.atleast_2d(np.asarray(L, dtype=float)).copy()
    N = np.asarray(N, dtype=float).flatten().copy()
    Z = np.asarray(Z, dtype=float).flatten().copy()
    mu = np.atleast_2d(np.asarray(mu, dtype=float)).copy()

    if options is None:
        options = {'tol': 1e-6}
    atol = options.get('tol', 1e-6)

    N = np.ceil(N)
    M = L.shape[0]
    R = L.shape[1]
    Nt = int(np.sum(N))

    # Sum Z across rows if 2D
    if Z.ndim > 1:
        Z = np.sum(Z, axis=0)

    # Handle case where Z is negligible
    if np.sum(Z) < ZERO:
        # see _kb/03-api-layer.md for rationale
        OneToMuCols = np.arange(1, mu.shape[1] + 1, dtype=float)

        zset = []
        non_zset = []
        for i in range(M):
            # Compare the FULL mu row with expected delay pattern [1, 2, 3, ..., mu_cols]
            if np.linalg.norm(mu[i, :] - OneToMuCols) < atol:
                zset.append(i)
            else:
                non_zset.append(i)

        if len(zset) > 0:
            Z = np.sum(L[zset, :], axis=0)
            L = L[non_zset, :] if len(non_zset) > 0 else np.zeros((0, R))
            mu = mu[non_zset, :] if len(non_zset) > 0 else np.zeros((0, mu.shape[1]))

    M = L.shape[0]

    # Handle negligible demands
    if np.sum(L) < ZERO:
        G, lG = pfqn_ca(L, N, Z)
        prob = np.zeros(Nt + 1)
        prob[Nt] = 1.0
        return PfqnComomrmLdResult(G=G, lG=lG, prob=prob)

    # Sanitize inputs
    lG0 = 0.0

    # Remove classes with zero demands and zero think times
    non_zero_classes = []
    for r in range(R):
        if np.sum(L[:, r]) >= atol or (len(Z) > r and Z[r] >= atol):
            non_zero_classes.append(r)
        else:
            if N[r] > 0:
                # Handle zero-demand classes
                pass

    if len(non_zero_classes) < R:
        L = L[:, non_zero_classes]
        N = N[non_zero_classes]
        if len(Z) > 0:
            Z = Z[non_zero_classes]
        R = len(non_zero_classes)

    # Handle empty cases
    if Z.size == 0 or np.sum(Z) < ZERO:
        if L.size == 0 or np.sum(L) < ZERO:
            prob = np.zeros(Nt + 1)
            prob[0] = 1.0
            return PfqnComomrmLdResult(G=exp(lG0), lG=lG0, prob=prob)
        Z = np.zeros(R) if R > 0 else np.zeros(1)
    elif L.size == 0 or np.sum(L) < ZERO:
        L = np.zeros((1, R)) if R > 0 else np.zeros((1, 1))

    M = L.shape[0]

    if M == 0:
        prob = np.zeros(Nt + 1)
        prob[Nt] = 1.0
        return PfqnComomrmLdResult(G=exp(lG0), lG=lG0, prob=prob)

    if M != 1:
        raise ValueError("pfqn_comomrm_ld: The solver accepts at most a single queueing station.")

    # COMOM algorithm
    h = np.zeros(Nt + 1)
    h[Nt] = 1.0
    scale = np.zeros(Nt)
    nt = 0

    for r in range(R):
        # Build transition matrix Tr
        Tr = np.eye(Nt + 1) * Z[r]
        for i in range(Nt):
            mu_idx = Nt - i - 1
            if mu_idx < mu.shape[1]:
                mu_val = mu[0, mu_idx]
            else:
                mu_val = 1.0

            if mu_val > 0:
                Tr[i, i + 1] = L[0, r] * (Nt - i) / mu_val

        nr = 0
        while nr < N[r]:
            hT = Tr / (1.0 + nr)
            h = hT @ h

            scale[nt] = np.abs(np.sum(np.sort(h)))
            h = np.abs(h)
            if scale[nt] > 0:
                h = h / scale[nt]
            nt += 1
            nr += 1

    # Compute final result
    with np.errstate(divide='ignore'):
        log_scale = np.log(scale)
        log_scale = np.where(np.isfinite(log_scale), log_scale, 0.0)

    lG = lG0 + np.sum(log_scale)
    G = exp(lG) if lG > -700 else 0.0

    prob = h[::-1]
    if G > 0:
        prob = prob / G
    prob = prob / np.sum(prob) if np.sum(prob) > 0 else prob

    return PfqnComomrmLdResult(G=G, lG=lG, prob=prob)


def pfqn_ncld(L: np.ndarray, N: np.ndarray, Z: np.ndarray,
              mu: np.ndarray, options: Optional[Dict[str, Any]] = None
              ) -> PfqnNcResult:
    """
    Main method to compute normalizing constant of a load-dependent model.

    Provides the main entry point for computing normalizing constants in
    load-dependent queueing networks with automatic method selection and
    preprocessing.

    Args:
        L: Service demands at all stations (M x R)
        N: Number of jobs for each class (1 x R)
        Z: Think times for each class (1 x R)
        mu: Load-dependent scalings (M x Ntot)
        options: Solver options with keys:
            - method: 'default', 'exact', 'rd', 'comomld', etc.
            - tol: Numerical tolerance

    Returns:
        PfqnNcResult with G (normalizing constant), lG (log), and method used
    """
    from .nc import pfqn_ca

    L = np.atleast_2d(np.asarray(L, dtype=float)).copy()
    N = np.asarray(N, dtype=float).flatten().copy()
    Z = np.asarray(Z, dtype=float).flatten().copy()
    mu = np.atleast_2d(np.asarray(mu, dtype=float)).copy()

    if options is None:
        options = {'method': 'default', 'tol': 1e-6}
    method = options.get('method', 'default')
    tol = options.get('tol', 1e-6)

    lG = np.nan
    G = np.nan

    Ntot = int(np.ceil(np.sum(N)))

    # Ensure mu has enough columns
    if Ntot > mu.shape[1]:
        # Extend mu with last column values
        extra_cols = Ntot - mu.shape[1]
        mu_extended = np.zeros((mu.shape[0], Ntot))
        mu_extended[:, :mu.shape[1]] = mu
        for i in range(mu.shape[1], Ntot):
            mu_extended[:, i] = mu[:, -1]
        mu = mu_extended
    elif Ntot < mu.shape[1]:
        mu = mu[:, :Ntot]

    # Remove classes with zero population
    L_new = []
    N_new = []
    Z_new = []
    for i in range(len(N)):
        if np.abs(N[i]) >= FINE_TOL:
            L_new.append(L[:, i])
            N_new.append(N[i])
            if i < len(Z):
                Z_new.append(Z[i])
            else:
                Z_new.append(0.0)

    if len(N_new) == 0:
        return PfqnNcResult(G=1.0, lG=0.0, method=method)

    L_new = np.column_stack(L_new) if len(L_new) > 0 else np.zeros((L.shape[0], 1))
    N_new = np.array(N_new)
    Z_new = np.array(Z_new)
    R = len(N_new)

    # Scaling for numerical stability
    scalevec = np.ones(R)
    for r in range(R):
        max_L = np.max(L_new[:, r]) if L_new.shape[0] > 0 else 0
        max_Z = Z_new[r] if r < len(Z_new) else 0
        scalevec[r] = max(max_L, max_Z, FINE_TOL)

    L_new = L_new / scalevec
    Z_new = Z_new / scalevec

    # Compute demand statistics
    Lsum = np.sum(L_new, axis=1)
    Lmax = np.max(L_new, axis=1)

    # Filter stations with non-zero demands
    dem_stations = []
    for i in range(L_new.shape[0]):
        with np.errstate(divide='ignore', invalid='ignore'):
            ratio = Lmax[i] / Lsum[i]
        if not np.isnan(ratio) and ratio > FINE_TOL:
            dem_stations.append(i)

    if len(dem_stations) > 0:
        L_new = L_new[dem_stations, :]
        mu = mu[dem_stations, :]
    else:
        L_new = np.zeros((0, R))
        mu = np.zeros((0, Ntot))

    M = L_new.shape[0]

    # Check for zero demands with positive population
    flag = False
    for i in range(R):
        L_sum_r = np.sum(L_new[:, i]) if M > 0 else 0
        Z_r = Z_new[i] if i < len(Z_new) else 0
        if np.abs(L_sum_r + Z_r) < FINE_TOL and N_new[i] > FINE_TOL:
            flag = True
            break

    if flag:
        print("pfqn_ncld warning: The model has no positive demands in any class.")
        if Z_new.size == 0 or np.sum(Z_new) < tol:
            lG = 0.0
        else:
            Z_sum = np.sum(Z_new)
            with np.errstate(divide='ignore'):
                log_Z = np.log(np.sum(Z_new))
                log_scale = np.log(scalevec)
            lG = (-np.sum(_factln_array(N_new)) +
                  np.dot(N_new, np.where(np.isfinite(log_Z), log_Z, 0.0) * np.ones(R)) +
                  np.dot(N_new, np.where(np.isfinite(log_scale), log_scale, 0.0)))
        return PfqnNcResult(G=np.nan, lG=lG, method=method)

    # Handle empty or negligible demands
    if L_new.size == 0 or np.sum(L_new) < tol:
        if Z_new.size == 0 or np.sum(Z_new) < tol:
            lG = 0.0
        else:
            with np.errstate(divide='ignore'):
                log_Z_sum = np.log(np.sum(Z_new, axis=0) if Z_new.ndim > 1 else Z_new)
                log_scale = np.log(scalevec)
            log_Z_sum = np.where(np.isfinite(log_Z_sum), log_Z_sum, 0.0)
            log_scale = np.where(np.isfinite(log_scale), log_scale, 0.0)
            lG = (-np.sum(_factln_array(N_new)) +
                  np.dot(N_new, log_Z_sum) + np.dot(N_new, log_scale))
        G = exp(lG) if lG > -700 else 0.0
        return PfqnNcResult(G=G, lG=lG, method=method)

    # Single station with no think times
    if M == 1 and (Z_new.size == 0 or np.sum(Z_new) < tol):
        with np.errstate(divide='ignore'):
            log_L_sum = np.log(np.sum(L_new, axis=0))
            log_scale = np.log(scalevec)
            log_mu = np.log(mu.flatten()[:Ntot]) if mu.size > 0 else np.zeros(Ntot)

        log_L_sum = np.where(np.isfinite(log_L_sum), log_L_sum, 0.0)
        log_scale = np.where(np.isfinite(log_scale), log_scale, 0.0)
        log_mu = np.where(np.isfinite(log_mu), log_mu, 0.0)

        lG = (_factln(np.sum(N_new)) - np.sum(_factln_array(N_new)) +
              np.dot(N_new, log_L_sum) + np.dot(N_new, log_scale) - np.sum(log_mu))
        G = exp(lG) if lG > -700 else 0.0
        return PfqnNcResult(G=G, lG=lG, method=method)

    # Separate zero-demand and nonzero-demand classes
    zero_demand_classes = []
    nonzero_demand_classes = []

    for i in range(R):
        if np.sum(L_new[:, i]) < tol:
            zero_demand_classes.append(i)
        else:
            nonzero_demand_classes.append(i)

    # Compute contribution from zero-demand classes (delay only)
    lGzdem = 0.0
    if len(zero_demand_classes) > 0:
        Zz = Z_new[zero_demand_classes]
        Nz = N_new[zero_demand_classes]
        scalevecz = scalevec[zero_demand_classes]

        if np.sum(Zz) >= tol:
            with np.errstate(divide='ignore'):
                log_Zz = np.log(Zz)
                log_scalevecz = np.log(scalevecz)
            log_Zz = np.where(np.isfinite(log_Zz), log_Zz, 0.0)
            log_scalevecz = np.where(np.isfinite(log_scalevecz), log_scalevecz, 0.0)

            lGzdem = (-np.sum(_factln_array(Nz)) +
                      np.dot(Nz, log_Zz) + np.dot(Nz, log_scalevecz))

    # Extract nonzero demand classes
    if len(nonzero_demand_classes) > 0:
        L_nnz = L_new[:, nonzero_demand_classes]
        N_nnz = N_new[nonzero_demand_classes]
        Z_nnz = Z_new[nonzero_demand_classes]
        scalevec_nnz = scalevec[nonzero_demand_classes]
    else:
        L_nnz = np.zeros((M, 1))
        N_nnz = np.zeros(1)
        Z_nnz = np.zeros(1)
        scalevec_nnz = np.ones(1)

    # Compute normalizing constant for nonzero demand classes
    lGnnzdem = 0.0
    if np.min(N_nnz) >= 0:
        result = _compute_norm_const_ld(L_nnz, N_nnz, Z_nnz, mu, options)
        lGnnzdem = result.lG
        method = result.method

    # Combine results
    with np.errstate(divide='ignore'):
        log_scalevec_nnz = np.log(scalevec_nnz)
    log_scalevec_nnz = np.where(np.isfinite(log_scalevec_nnz), log_scalevec_nnz, 0.0)

    lG = lGnnzdem + lGzdem + np.dot(N_nnz, log_scalevec_nnz)
    G = exp(lG) if lG > -700 else 0.0

    return PfqnNcResult(G=G, lG=lG, method=method)


def _compute_norm_const_ld(L: np.ndarray, N: np.ndarray, Z: np.ndarray,
                           mu: np.ndarray, options: Dict[str, Any]
                           ) -> PfqnNcResult:
    """
    Run a normalizing constant solution method on a load-dependent model.

    Internal function that dispatches to the appropriate algorithm based
    on the method option.

    Args:
        L: Service demands at all stations (M x R)
        N: Number of jobs for each class (1 x R)
        Z: Think times for each class (1 x R)
        mu: Load-dependent scalings (M x Ntot)
        options: Solver options

    Returns:
        PfqnNcResult with G, lG, and method used
    """
    M = L.shape[0]
    R = L.shape[1]
    method = options.get('method', 'default')
    lG = None

    # Ensure N has R elements to match L's columns
    N = np.atleast_1d(N).flatten()
    if len(N) != R:
        # Dimension mismatch - this can happen with chain aggregation
        # Try to handle gracefully by falling back to exact method
        if method not in ['default', 'exact']:
            method = 'exact'

    # see _kb/03-api-layer.md for rationale
    CLW_MAX_CLASSES = 5    # class-count gate (clw cost is exponential in R)
    CLW_MAX_POP = 200      # total-population cap (numerical validity; NaN onset ~450)
    CLW_MAX_COST = 2e7     # contour-point budget (~2s at ~1e7 pts/s, see profiler)
    _lvec = np.full(R, 3.0)
    if R >= 1:
        _lvec[0] = 1.0
    if R >= 2:
        _lvec[1] = 2.0
    if R >= 3:
        _lvec[2] = 2.0
    clw_pred_cost = float(np.prod(2.0 * _lvec * np.asarray(N, dtype=float)))
    if (method == 'default' and M > 1 and R >= 2
            and R <= CLW_MAX_CLASSES and float(np.sum(N)) <= CLW_MAX_POP
            and clw_pred_cost <= CLW_MAX_COST):
        from .nc import pfqn_clw_lld
        Z_row = np.sum(Z, axis=0) if np.ndim(Z) > 1 else Z
        _, lG = pfqn_clw_lld(L, N, Z_row, mu)
        method = "clw"
        lG_real = np.real(lG)
        G = exp(lG_real) if lG_real > -700 else 0.0
        return PfqnNcResult(G=G, lG=lG_real, method=method)

    if method in ['default', 'exact']:
        # Combine L and Z for stations with infinite servers
        if np.sum(Z) < FINE_TOL:
            Lz = L
            muz = mu
        else:
            D = 1  # Z is 1D
            Lz = np.vstack([L, Z.reshape(1, -1)])

            # Create mu for delay stations
            Ntot = mu.shape[1]
            delay_mu = np.arange(1, Ntot + 1, dtype=float).reshape(1, -1)
            muz = np.vstack([mu, delay_mu])

        if R == 1:
            result = pfqn_gldsingle(Lz, N, muz, options)
            lG = result.lG
            method = "exact/gld"
        elif M == 1 and np.max(Z) > 0:
            result = pfqn_comomrm_ld(L, N, Z, muz, options)
            lG = result.lG
            method = "exact/comomld"
        elif M == 1 and np.max(Z) < FINE_TOL:
            # see _kb/03-api-layer.md for rationale
            result = pfqn_comomrm_ld(L, N, np.zeros_like(N), mu, options)
            lG = result.lG
            method = "exact/comomld"
        else:
            result = pfqn_gld(Lz, N, muz, options)
            lG = result.lG
            method = "exact/gld"

    elif method == 'is':
        # see _kb/03-api-layer.md for rationale
        Z_row = np.sum(Z, axis=0) if np.ndim(Z) > 1 else Z
        lG = pfqn_ld_is(L, N, Z_row, mu, options).lG
        method = "is"

    elif method == 'clw':
        # see _kb/03-api-layer.md for rationale
        from .nc import pfqn_clw_lld
        Z_row = np.sum(Z, axis=0) if np.ndim(Z) > 1 else Z
        _, lG = pfqn_clw_lld(L, N, Z_row, mu)
        method = "clw"

    elif method in ('panacea', 'panaceald'):
        # Mitra-McKenna load-dependent PANACEA asymptotic expansion. Delay terms
        # may arrive either in Z or as mu(i,n)=n rows of L, both are recognized
        # by pfqn_panaceald.
        Z_row = np.sum(Z, axis=0) if np.ndim(Z) > 1 else Z
        _, lG = pfqn_panaceald(L, N, Z_row, mu)
        method = "panaceald"
        if np.isnan(lG):
            # normal usage (1 - lambda_i/mu_i(Ntot) > 0 at every queueing
            # center) is the domain of the expansion, not a numerical failure
            raise ValueError(
                "The model is not in normal usage, so the 'panaceald' "
                "asymptotic expansion does not apply. Use 'exact', 'clw' or an "
                "approximate load-dependent method instead.")

    elif method == 'comomld':
        if M <= 1 or np.sum(Z) <= ZERO:
            result = pfqn_comomrm_ld(L, N, Z, mu, options)
            lG = result.lG
        else:
            print("pfqn_ncld warning: Load-dependent CoMoM is available only in "
                  "models with a delay and m identical stations.")
            # Fall back to gld
            if np.sum(Z) < FINE_TOL:
                Lz = L
                muz = mu
            else:
                Lz = np.vstack([L, Z.reshape(1, -1)])
                Ntot = mu.shape[1]
                delay_mu = np.arange(1, Ntot + 1, dtype=float).reshape(1, -1)
                muz = np.vstack([mu, delay_mu])
            result = pfqn_gld(Lz, N, muz, options)
            lG = result.lG
            method = "gld"

    elif method == 'nrl':
        from .laplace import pfqn_nrl
        lG = pfqn_nrl(L, N, Z, alpha=mu)
        method = "nrl"

    elif method == 'nrp':
        from .laplace import pfqn_nrp
        lG = pfqn_nrp(L, N, Z, alpha=mu)
        method = "nrp"

    elif method == 'rd':
        from .rd import pfqn_rd
        result = pfqn_rd(L, N, Z, mu=mu)
        lG = result[0] if isinstance(result, tuple) else result.lGN
        method = "rd"

    else:
        # Default to exact/gld
        if np.sum(Z) < FINE_TOL:
            Lz = L
            muz = mu
        else:
            Lz = np.vstack([L, Z.reshape(1, -1)])
            Ntot = mu.shape[1]
            delay_mu = np.arange(1, Ntot + 1, dtype=float).reshape(1, -1)
            muz = np.vstack([mu, delay_mu])
        result = pfqn_gld(Lz, N, muz, options)
        lG = result.lG
        method = "exact/gld"

    # Handle complex lG by taking real part for comparison
    lG_real = np.real(lG) if lG is not None else NEG_INF
    G = exp(lG_real) if lG_real > -700 else 0.0
    lG = lG_real  # Use real part for result

    return PfqnNcResult(G=G, lG=lG if lG is not None else NEG_INF, method=method)


@dataclass
class PfqnFncResult:
    """Result of functional server scaling computation."""
    mu: np.ndarray
    c: np.ndarray


def pfqn_fnc(alpha: np.ndarray, c: Optional[np.ndarray] = None) -> PfqnFncResult:
    """
    Compute scaling factor of a load-dependent functional server.

    Used to calculate the mean queue length in load-dependent systems by
    computing functional scaling factors from load-dependent service rate
    parameters.

    Args:
        alpha: Load-dependent scalings (M x N)
        c: Scaling constants (1 x M), optional. If None, auto-selected.

    Returns:
        PfqnFncResult with mu (functional server scalings) and c (scaling constants)
    """
    alpha = np.atleast_2d(np.asarray(alpha, dtype=float))
    M = alpha.shape[0]
    N = alpha.shape[1]

    if N == 0:
        # see _kb/03-api-layer.md for rationale
        return PfqnFncResult(mu=np.zeros((M, 0)), c=np.zeros((1, M)))

    if c is None:
        # First try c = 0
        c = np.zeros((1, M))
        result = _pfqn_fnc_with_c(alpha, c)

        if not np.all(np.isfinite(result.mu)):
            # Try c = -0.5
            c = np.full((1, M), -0.5)
            result = _pfqn_fnc_with_c(alpha, c)

        # If still not finite, search for valid c
        dt = 0.0
        while not np.all(np.isfinite(result.mu)):
            dt += 0.05
            # see _kb/03-api-layer.md for rationale
            c = np.full((1, M), -0.5 + dt)
            result = _pfqn_fnc_with_c(alpha, c)
            if (-0.5 + dt) >= 2:
                break

        return result
    else:
        c = np.atleast_2d(np.asarray(c, dtype=float))
        return _pfqn_fnc_with_c(alpha, c)


def _pfqn_fnc_with_c(alpha: np.ndarray, c: np.ndarray) -> PfqnFncResult:
    """
    Compute functional server scalings with specified scaling constant.

    Internal function that performs the actual computation of functional
    server scalings.

    Args:
        alpha: Load-dependent scalings (M x N)
        c: Scaling constants (1 x M)

    Returns:
        PfqnFncResult with mu and c
    """
    alpha = np.atleast_2d(np.asarray(alpha, dtype=float))
    c = np.atleast_2d(np.asarray(c, dtype=float)).flatten()

    M = alpha.shape[0]
    N = alpha.shape[1]

    mu = np.zeros((M, N))

    for i in range(M):
        c_i = c[i] if i < len(c) else 0.0
        mu[i, 0] = alpha[i, 0] / (1 + c_i)

        alphanum = np.zeros((N, N))
        alphaden = np.zeros((N, N))

        for n in range(1, N):
            alphanum[n, 0] = alpha[i, n]
            alphaden[n, 0] = alpha[i, n - 1]

            for k in range(1, n):
                alphanum[n, k] = alphanum[n, k - 1] * alpha[i, n - k]
                alphaden[n, k] = alphaden[n, k - 1] * alpha[i, n - k - 1]

        for n in range(1, N):
            rho = 0.0
            muden = 1.0

            for k in range(n):
                with np.errstate(invalid='ignore'):
                    muden *= mu[i, k]
                if muden != 0 and np.isfinite(muden):
                    rho += (alphanum[n, k] - alphaden[n, k]) / muden

            if muden != 0 and np.isfinite(muden) and (1 - rho) != 0:
                mu[i, n] = (alphanum[n, n - 1] * alpha[i, 0] / muden) / (1 - rho)
            else:
                mu[i, n] = np.inf

    # Clean up non-finite values
    for i in range(M):
        for j in range(N):
            if np.isnan(mu[i, j]) or np.abs(mu[i, j]) > 1e15:
                mu[i, j] = np.inf

    # Replace values after first inf with inf
    for i in range(M):
        if not np.all(np.isfinite(mu[i, :])):
            replace_with_inf = False
            for j in range(N):
                if replace_with_inf:
                    mu[i, j] = np.inf
                elif np.isinf(mu[i, j]):
                    replace_with_inf = True

    return PfqnFncResult(mu=mu, c=c.reshape(1, -1))


def pfqn_ld_is(L, N, Z=None, mu=None, options=None) -> PfqnNcResult:
    """
    Importance-sampling (IS) estimate of the normalizing constant of a closed
    LOAD-DEPENDENT product-form queueing network. Load-dependent counterpart of
    pfqn_pas_is / pfqn_oi_is: the same sample-an-ordering estimator, with the
    order-independent rank rate replaced by the load-dependent capacity.

    Identity. Every product-form station's balance function is the sum, over the
    orderings q of a given per-class count vector n, of an ordered product of a
    per-position factor::

        F_i(n) = |n|!/prod_r(n_r!) * prod_r L(i,r)^{n_r} / prod_{k=1}^{|n|} mu_i(k)
               = sum_{q: |q|=n} prod_{p=1}^{|n|} L(i,q_p) / mu_i(p)

    since the multiset has |n|!/prod_r(n_r!) orderings, each contributing the same
    ordered product. The delay (infinite-server) node is the special case
    mu_Z(k)=k, giving F_Z(n)=prod_r Z_r^{n_r}/n_r!; a single-server queue is
    mu_i(k)=1; a c-server queue is mu_i(k)=min(k,c).

    Consequently, with ell = sum(N) and a "cut vector" splitting an ordering c of
    all ell jobs into S contiguous segments (one per station)::

        G(N) = sum_{c} sum_{cuts} prod_{m=1}^{S} w_m(seg_m),
        w_m(q) = prod_{p=1}^{|q|} L(m,q_p) / mu_m(p)

    because summing over the orderings of each segment independently reproduces
    prod_m F_m(n_m), and each count split is realized exactly once.

    Estimator. An ordering c is drawn by placing, at each step, a uniformly random
    present class; p(c) is the product of the reciprocal branching factors. For the
    sampled c the inner sum over ALL cut vectors is computed exactly by the dynamic
    program A_0(0)=1, A_m(k) = sum_{j<=k} A_{m-1}(j) * w_m(c_{j+1..k}), so
    S(c)=A_S(ell) in O(S*ell^2) time (no cut enumeration). Then
    G = E_{C~p}[S(C)/p(C)] is unbiased, estimated by the sample mean.

    Parameters
    ----------
    L : (M, R) array
        Per-class service demands at the M queueing stations.
    N : (R,) array
        Closed population vector, finite.
    Z : (R,) array, optional
        Aggregated think time (delay) demand; None or zeros if none.
    mu : (M, ell) array or sequence of callables, optional
        Load-dependent capacities; ``mu[i][k-1]`` is the capacity of station i
        holding k jobs. None for the load-independent case mu(i,k)=1 (see
        :func:`pfqn_is`).
    options : dict or options object, optional
        Fields ``samples`` (default 1e4) and ``seed`` (optional).

    Returns
    -------
    PfqnNcResult with ``G`` the IS estimate of the normalizing constant and
    ``lG = log(G)``.

    Examples
    --------
    >>> L = np.array([[0.5, 0.3], [0.2, 0.4]]); N = np.array([3, 2]); Z = np.array([1.0, 1.0])
    >>> mu = np.array([[1, 2, 2, 2, 2], [1, 1, 1, 1, 1]], dtype=float)
    >>> res = pfqn_ld_is(L, N, Z, mu, {'samples': 100000, 'seed': 7})

    See Also
    --------
    pfqn_is, pfqn_ncld, pfqn_nc
    """
    from .pas import _opt

    L = np.asarray(L, dtype=float)
    if L.ndim == 1:
        L = L.reshape(1, -1)
    M, R = L.shape
    N = np.round(np.asarray(N, dtype=float)).astype(int).ravel()
    if N.size != R:
        raise ValueError('L must have as many columns as N has classes.')
    if not np.all(np.isfinite(N)):
        raise ValueError('pfqn_ld_is requires finite (closed) populations.')
    if Z is None:
        Z = np.zeros(R)
    Z = np.asarray(Z, dtype=float)
    if Z.ndim > 1:
        Z = np.sum(Z, axis=0)
    Z = Z.ravel()
    ell = int(np.sum(N))

    nsamples = int(round(_opt(options, 'samples', 10000)))
    seed = _opt(options, 'seed', None)
    rng = np.random.default_rng(seed if seed is not None else None)

    if ell == 0:
        return PfqnNcResult(G=1.0, lG=0.0, method='is')

    # ---- station list: M queues, plus the delay as mu_Z(k)=k ----------------
    # D[m, r] per-class demand of station m; B[m, k-1] its capacity at k jobs.
    has_z = bool(np.any(Z > 0))
    S = M + (1 if has_z else 0)
    D = np.zeros((S, R))
    B = np.ones((S, ell))
    for i in range(M):
        D[i, :] = L[i, :]
        if mu is None:
            B[i, :] = 1.0                      # load-independent single server
        elif callable(mu[i]):
            for k in range(1, ell + 1):
                B[i, k - 1] = mu[i](k)
        else:
            mu_i = np.asarray(mu, dtype=float)[i, :]
            ncol = min(ell, mu_i.size)
            B[i, :ncol] = mu_i[:ncol]
            if mu_i.size < ell:                # extend with the last capacity
                B[i, mu_i.size:ell] = mu_i[-1]
    if has_z:
        D[S - 1, :] = Z
        B[S - 1, :] = np.arange(1, ell + 1)    # delay: mu_Z(k) = k
    if np.any(B <= 0):
        raise ValueError('load-dependent capacities must be strictly positive.')

    acc = 0.0
    for _ in range(nsamples):
        # ---- draw an ordering c (uniformly random present class each step) --
        x = N.copy()
        c = np.zeros(ell, dtype=int)
        logp = 0.0
        for ppos in range(ell):
            avail = np.flatnonzero(x > 0)
            na = avail.size
            pick = int(avail[rng.integers(na)])
            c[ppos] = pick
            logp -= log(na)
            x[pick] -= 1

        # see _kb/03-api-layer.md for rationale
        A = np.zeros(ell + 1)
        A[0] = 1.0
        for m in range(S):
            Anew = np.zeros(ell + 1)
            for j in range(ell + 1):
                if A[j] == 0.0:
                    continue
                Anew[j] += A[j]                       # empty segment
                w = 1.0
                for k in range(j + 1, ell + 1):
                    w *= D[m, c[k - 1]] / B[m, k - j - 1]   # position in segment
                    if w == 0.0:
                        break
                    Anew[k] += A[j] * w
            A = Anew
        acc += A[ell] * exp(-logp)

    G = acc / nsamples
    lG = log(G) if G > 0 else -np.inf
    return PfqnNcResult(G=G, lG=lG, method='is')


def pfqn_panaceald(L: np.ndarray, N: np.ndarray, Z: np.ndarray = None,
                   mu: np.ndarray = None, terms: int = 3) -> Tuple[float, float]:
    """
    PANACEA asymptotic expansion for load-dependent closed networks.

    Mitra-McKenna (JACM 33(3):568-592, 1986) load-dependent PANACEA: the
    expansion coefficients A_n are linear combinations of partition functions
    of a pseudonetwork whose load dependence is the phi(n) transform of the
    original {f(n)}. See _kb/03-api-layer.md (pfqn/ family, panaceald).

    Args:
        L: Service demand matrix (M x R)
        N: Population vector (R,)
        Z: Think time vector (R,) or matrix (D x R), summed over rows
        mu: Load-dependent rate matrix (M x sum(N))
        terms: Number of terms in the normal-usage asymptotic series (1, 2 or 3)

    Returns:
        Tuple (G, lG) - normalizing constant and its log, both NaN when the
        model is not in normal usage
    """
    if terms not in (1, 2, 3):
        raise ValueError("The terms parameter must be 1, 2, or 3 "
                         "(higher-order coefficients are not implemented).")
    L = np.atleast_2d(np.asarray(L, dtype=float))
    N = np.asarray(N, dtype=float).flatten()
    M, R = L.shape
    if Z is None:
        Ztot = np.zeros(R)
    else:
        Z = np.asarray(Z, dtype=float)
        Ztot = np.sum(Z, axis=0) if Z.ndim > 1 else Z.copy()
    Ntot = int(round(float(np.sum(N))))
    if Ntot == 0:
        return 1.0, 0.0
    if mu is None:
        mu = np.ones((M, Ntot))
    else:
        mu = np.atleast_2d(np.asarray(mu, dtype=float))
    if mu.shape[1] < Ntot:
        mu = np.hstack([mu, np.tile(mu[:, -1:], (1, Ntot - mu.shape[1]))])

    # Type-3 (infinite-server) rows are absent from the pseudonetwork and enter
    # only through rho_j0; solver_ncld encodes them as mu(i,n)=n rows of L.
    lattice = np.arange(1, Ntot + 1, dtype=float)
    is_is = np.all(np.abs(mu[:, :Ntot] - lattice) < FINE_TOL, axis=1)
    if np.any(is_is):
        Ztot = Ztot + np.sum(L[is_is, :], axis=0)
    Lq = L[~is_is, :]
    muq = mu[~is_is, :Ntot]
    Mq = Lq.shape[0]

    if np.any((N > 0) & (Ztot <= 0)):
        # no IS center on the route of a populated class: the expansion
        # parameter rho_j0 is undefined and PANACEA does not apply
        return float('nan'), float('nan')

    lGdelay = -float(np.sum([_factln(n) for n in N]))
    for j in range(R):
        if N[j] != 0:
            lGdelay += N[j] * log(Ztot[j])
    if Mq == 0:
        return exp(lGdelay), lGdelay
    if np.any(muq <= 0) or not np.all(np.isfinite(muq)):
        return float('nan'), float('nan')

    r = np.zeros((Mq, R))
    for j in range(R):
        if Ztot[j] > 0:
            r[:, j] = Lq[:, j] / Ztot[j]
    lam = r.dot(N)
    muK = muq[:, Ntot - 1]
    alpha = 1.0 - lam / muK
    if np.min(alpha) <= 0:
        # model is not in normal usage: the {phi(n)} series diverges
        return float('nan'), float('nan')

    # log-partial products log prod_{k=1}^{s} mu_i(k), s=0..Ntot
    lPi = np.hstack([np.zeros((Mq, 1)), np.cumsum(np.log(muq), axis=1)])

    nmax = 2 * (terms - 1)
    lpsi = np.zeros((Mq, nmax + 1))
    for i in range(Mq):
        for n in range(nmax + 1):
            lpsi[i, n] = _logpsi(n, lam[i], lPi[i, :], muK[i], alpha[i], Ntot)

    # load dependence of the pseudonetwork centers:
    # psi_i(n) = psi_i(0) n! / prod_{k=1}^{n} mups_i(k)
    mups = np.ones((Mq, max(1, nmax)))
    for i in range(Mq):
        for n in range(1, nmax + 1):
            mups[i, n - 1] = exp(log(n) + lpsi[i, n - 1] - lpsi[i, n])

    # Expansion coefficients (5.4). The large parameter N cancels identically
    # between beta_j=K_j/N, Gamma=N*r and the 1/N^n scaling, so the demands are
    # taken as r and beta as N.
    A = [1.0, 0.0, 0.0]
    if terms >= 2:
        for j in range(R):
            k = np.zeros(R, dtype=int)
            k[j] = 2
            A[1] -= N[j] * _pseudonet(r, k, mups)
    if terms >= 3:
        for j in range(R):
            k = np.zeros(R, dtype=int)
            k[j] = 3
            A[2] += 2 * N[j] * _pseudonet(r, k, mups)
            k[j] = 4
            A[2] += 3 * N[j] ** 2 * _pseudonet(r, k, mups)
            for s in range(R):
                if s == j:
                    continue
                k2 = np.zeros(R, dtype=int)
                k2[j] = 2
                k2[s] = 2
                A[2] += 0.5 * N[j] * N[s] * _pseudonet(r, k2, mups)
    I = sum(A[:terms])
    if I <= 0:
        return float('nan'), float('nan')

    lG = lGdelay + float(np.sum(lpsi[:, 0])) + log(I)
    if not np.isfinite(lG):
        return float('nan'), float('nan')
    return exp(lG), lG


def _logpsi(n: int, lam: float, lPirow: np.ndarray, muK: float,
            alpha: float, K: int) -> float:
    """
    log of psi(n) = sum_{s>=n} [s!/(s-n)!] lam^(s-n) / prod_k mu(k), the mu-free
    part of the phi(n) transform in eq. (3.7)-(3.8a). The series is split into
    the exact head s<=K and a geometric tail summed in closed form via the
    Vandermonde identity, all terms positive.
    """
    t = []
    for s in range(n, K + 1):
        t.append(_factln(s) - _factln(s - n) + _xlogy(s - n, lam) - lPirow[s])
    T = max(n, K + 1)
    for i in range(n + 1):
        t.append(_factln(n) + _factln(T) - _factln(n - i) - _factln(T - n + i)
                 + _xlogy(T + i - n, lam) + (K - T - i) * log(muK)
                 - (i + 1) * log(alpha) - lPirow[K])
    tmax = max(t)
    if tmax == NEG_INF:
        return NEG_INF
    return tmax + log(sum(exp(v - tmax) for v in t))


def _pseudonet(gam: np.ndarray, k: np.ndarray, mups: np.ndarray) -> float:
    """
    Partition function of the pseudonetwork at population k, normalized so that
    G(0)=1. Populations are at most 2*(terms-1), so a direct load-dependent
    convolution over the population lattice is used.
    """
    nz = [j for j in range(len(k)) if k[j] > 0]
    Mq = gam.shape[0]
    sizes = [int(k[j]) + 1 for j in nz]
    npop = 1
    for sz in sizes:
        npop *= sz

    def idx2vec(idx):
        v = []
        t = idx
        for sz in sizes:
            v.append(t % sz)
            t //= sz
        return v

    def vec2idx(v):
        idx = 0
        mult = 1
        for j, sz in enumerate(sizes):
            idx += mult * v[j]
            mult *= sz
        return idx

    sterm = np.zeros((Mq, npop))
    for i in range(Mq):
        for jdx in range(npop):
            m = idx2vec(jdx)
            sm = int(sum(m))
            v = _factln(sm)
            zero = False
            for jj, cls in enumerate(nz):
                if m[jj] > 0:
                    if gam[i, cls] <= 0:
                        zero = True
                        break
                    v += m[jj] * log(gam[i, cls]) - _factln(m[jj])
            if zero:
                sterm[i, jdx] = 0.0
            else:
                for l in range(sm):
                    v -= log(mups[i, l])
                sterm[i, jdx] = exp(v)

    g = np.zeros(npop)
    g[0] = 1.0
    for i in range(Mq):
        gnew = np.zeros(npop)
        for idx in range(npop):
            nvec = idx2vec(idx)
            acc = 0.0
            for jdx in range(npop):
                m = idx2vec(jdx)
                if all(m[j] <= nvec[j] for j in range(len(sizes))):
                    acc += g[vec2idx([nvec[j] - m[j] for j in range(len(sizes))])] * sterm[i, jdx]
            gnew[idx] = acc
        g = gnew
    return float(g[npop - 1])


def _xlogy(e: int, x: float) -> float:
    """e*log(x) with the convention 0*log(0)=0."""
    if e == 0:
        return 0.0
    return e * log(x) if x > 0 else NEG_INF


__all__ = [
    'pfqn_ncld',
    'pfqn_panaceald',
    'pfqn_gld',
    'pfqn_gldsingle',
    'pfqn_mushift',
    'pfqn_comomrm_ld',
    'pfqn_fnc',
    'pfqn_ld_is',
    'PfqnNcResult',
    'PfqnComomrmLdResult',
    'PfqnFncResult',
]
