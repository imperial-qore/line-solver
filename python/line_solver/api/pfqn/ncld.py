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
from math import log, exp, factorial, lgamma
from typing import Tuple, Dict, Optional, Any
from dataclasses import dataclass


# Fine tolerance for numerical comparisons
FINE_TOL = 1e-12
ZERO = 1e-10
NEG_INF = float('-inf')


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

    Implements the Convolution Method of Moments (COMOM) for computing normalizing
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
        # Identify delay stations (mu = [1, 2, 3, ...])
        # A delay station has mu[j] = j+1 for ALL columns, not just the first Nt
        # This prevents false positives when Nt is small
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
        elif M == 2 and np.max(Z) < FINE_TOL:
            # For M==2 with no Z, pfqn_comomrm_ld assumes one station is a delay
            # Check if exactly one station matches the delay pattern [1, 2, 3, ...]
            # If mu was trimmed to small Ntot, both stations might falsely match
            delay_count = 0
            OneToMuCols = np.arange(1, mu.shape[1] + 1, dtype=float)
            atol = options.get('tol', 1e-6)
            for i in range(M):
                if np.linalg.norm(mu[i, :] - OneToMuCols) < atol:
                    delay_count += 1

            if delay_count == 1:
                result = pfqn_comomrm_ld(L, N, np.zeros_like(N), mu, options)
                lG = result.lG
                method = "exact/comomld"
            else:
                # Either no delays or both stations appear as delays (e.g., small Ntot)
                # Fall back to general solver
                result = pfqn_gld(Lz, N, muz, options)
                lG = result.lG
                method = "exact/gld"
        else:
            result = pfqn_gld(Lz, N, muz, options)
            lG = result.lG
            method = "exact/gld"

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
            c_scalar = -0.5 + dt
            c = np.array([[c_scalar]])
            result = _pfqn_fnc_with_c(alpha, c)
            if c_scalar >= 2:
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


__all__ = [
    'pfqn_ncld',
    'pfqn_gld',
    'pfqn_gldsingle',
    'pfqn_mushift',
    'pfqn_comomrm_ld',
    'pfqn_fnc',
    'PfqnNcResult',
    'PfqnComomrmLdResult',
    'PfqnFncResult',
]
