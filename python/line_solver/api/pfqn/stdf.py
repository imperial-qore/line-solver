"""
Sojourn Time Distribution Functions for PFQN.

Native Python implementation for computing sojourn time distributions
in product-form queueing networks with multiserver FCFS nodes.

References:
    McKenna, J. "Mean Value Analysis for networks with service-time-dependent
    queue disciplines." Performance Evaluation, 1987.
"""

import numpy as np
from typing import List, Dict, Optional, Tuple
from scipy import linalg

from ..mam import (
    map_exponential, map_erlang, map_sumind, map_cdf
)
from .utils import oner, factln
from .mvald import pfqn_mvald
from .ncld import pfqn_comomrm_ld
from .schmidt import pprod


def pfqn_stdf(L: np.ndarray, N: np.ndarray, Z: np.ndarray,
              S: np.ndarray, fcfs_nodes: np.ndarray, rates: np.ndarray,
              tset: np.ndarray) -> Dict[Tuple[int, int], np.ndarray]:
    """
    Compute sojourn time distribution for multiserver FCFS nodes.

    Implements McKenna's method for computing the response time distribution
    at multiserver FCFS stations in closed queueing networks.

    Args:
        L: Load matrix (M x R) - service demands
        N: Population vector (R,) - number of jobs per class
        Z: Think time vector (R,) or matrix (1 x R)
        S: Number of servers at each station (M,)
        fcfs_nodes: Array of FCFS station indices (1-indexed as in MATLAB)
        rates: Service rates matrix (M x R)
        tset: Time points at which to evaluate the distribution

    Returns:
        Dictionary with (station, class) tuples as keys and
        response time distribution arrays as values.
        Each array has shape (len(tset), 2) with:
            - column 0: CDF values
            - column 1: time points

    Examples:
        >>> L = np.array([[0.5, 0.3], [0.2, 0.4]])
        >>> N = np.array([2, 3])
        >>> Z = np.array([1.0, 1.0])
        >>> S = np.array([2, 1])
        >>> fcfs_nodes = np.array([0])  # 0-indexed
        >>> rates = np.array([[1.0, 1.0], [2.0, 2.0]])
        >>> tset = np.linspace(0.1, 5.0, 50)
        >>> RD = pfqn_stdf(L, N, Z, S, fcfs_nodes, rates, tset)
    """
    L = np.atleast_2d(np.asarray(L, dtype=float))
    N = np.asarray(N, dtype=int).flatten()
    Z = np.atleast_1d(np.asarray(Z, dtype=float)).flatten()
    S = np.asarray(S, dtype=int).flatten()
    fcfs_nodes = np.asarray(fcfs_nodes, dtype=int).flatten()
    rates = np.atleast_2d(np.asarray(rates, dtype=float))
    tset = np.asarray(tset, dtype=float).flatten()

    M, R = L.shape
    T = len(tset)

    # Replace zero times with small value to avoid numerical issues
    FINE_TOL = 1e-10
    tset = np.where(tset == 0, FINE_TOL, tset)

    # Build mu matrix (load-dependent service rates)
    mu = np.zeros((M, int(np.sum(N))))
    for k in range(M):
        for n in range(int(np.sum(N))):
            mu[k, n] = min(S[k], n + 1)

    stability_warn_issued = False
    RD = {}

    for k in fcfs_nodes:
        # Check that all rates at this station are equal (product-form requirement)
        if np.ptp(rates[k, :]) > FINE_TOL:
            raise ValueError(
                f"FCFS station {k} has distinct service rates, the model is invalid."
            )

        # Build h_kc(t, n) - sojourn time CDF for queue length n
        hMAP = {}
        hkc = np.zeros((T, 1 + int(np.sum(N))))

        for n in range(int(np.sum(N)) + 1):
            if n < S[k]:
                # Exponential service (no queueing)
                D0, D1 = map_exponential(1.0 / rates[k, 0])
                hMAP[n] = (D0, D1)
            else:
                # Service + waiting in queue
                exp_map = map_exponential(1.0 / rates[k, 0])
                erlang_phases = n - S[k] + 1
                erlang_mean = (n - S[k] + 1) / (S[k] * rates[k, 0])
                erlang = map_erlang(erlang_mean, erlang_phases)
                D0, D1 = map_sumind([exp_map, erlang])
                hMAP[n] = (D0, D1)

            # Compute CDF at all time points
            D0, D1 = hMAP[n]
            cdf_vals = map_cdf(D0, D1, tset)
            hkc[:, n] = cdf_vals

        # Compute response time distribution for each class
        for r in range(R):
            if L[k, r] > FINE_TOL:
                Nr = oner(N, r)

                # Compute normalizing constant for reduced population
                if M == 1:
                    # Use COMOM for single station
                    _, lGr = _pfqn_comomrm_ld_wrapper(L, Nr, Z, mu)
                    is_num_stable = True
                else:
                    # Use load-dependent MVA (returns 7 values)
                    _, _, _, _, lGr, is_num_stable, _ = pfqn_mvald(L, Nr, Z, mu)

                lGr = lGr[-1] if hasattr(lGr, '__len__') else lGr

                if not is_num_stable and not stability_warn_issued:
                    stability_warn_issued = True

                # Initialize result array
                RD[(k, r)] = np.zeros((T, 2))
                RD[(k, r)][:, 1] = tset

                # Use recursive form for load-dependent models (faster)
                Hkrt = np.zeros(T)

                # Compute lGk for network without station k
                other_stations = np.array([i for i in range(M) if i != k])
                if len(other_stations) > 0:
                    if len(other_stations) == 1:
                        _, lGk = _pfqn_comomrm_ld_wrapper(
                            L[other_stations, :], Nr, Z, mu[other_stations, :]
                        )
                    else:
                        _, _, _, _, lGk, _, _ = pfqn_mvald(
                            L[other_stations, :], Nr, Z, mu[other_stations, :]
                        )
                    lGk = lGk[-1] if hasattr(lGk, '__len__') else lGk
                else:
                    lGk = 0.0

                for t_idx in range(T):
                    gammat = mu.copy()

                    # Adjust gamma for time-dependent service
                    for m in range(1, int(np.sum(Nr))):
                        if hkc[t_idx, m] > FINE_TOL:
                            gammat[k, m - 1] = mu[k, m - 1] * hkc[t_idx, m - 1] / hkc[t_idx, m]

                    gammak = _mushift(gammat, k)
                    if gammak.shape[1] > 0:
                        gammak = gammak[:, :int(np.sum(Nr)) - 1]

                    # nvec = 0 case
                    Hkrt[t_idx] = hkc[t_idx, 0] * np.exp(lGk)

                    # nvec >= 1 cases
                    for s in range(R):
                        if Nr[s] > 0:
                            Nr_s = oner(Nr, s)

                            if M == 1:
                                _, lYks_t = _pfqn_comomrm_ld_wrapper(
                                    L, Nr_s, Z, gammak
                                )
                            else:
                                _, _, _, _, lYks_t, _, _ = pfqn_mvald(
                                    L, Nr_s, Z, gammak
                                )

                            lYks_t = lYks_t[-1] if hasattr(lYks_t, '__len__') else lYks_t

                            if gammat[k, 0] > FINE_TOL:
                                Hkrt[t_idx] += (L[k, s] * hkc[t_idx, 0] / gammat[k, 0]) * np.exp(lYks_t)

                # Handle NaN values
                Hkrt = np.where(np.isnan(Hkrt), FINE_TOL, Hkrt)
                lHkrt = np.log(np.maximum(Hkrt, FINE_TOL))

                # Compute final CDF
                RD[(k, r)][:, 0] = np.minimum(1.0, np.exp(lHkrt - lGr))

    return RD


def _mushift(mu: np.ndarray, i: int) -> np.ndarray:
    """
    Shift the service rate vector.

    For station i, shifts rates left by one position.
    For other stations, keeps rates unchanged.

    Args:
        mu: Service rate matrix (M x N)
        i: Station index to shift

    Returns:
        Shifted service rate matrix
    """
    M, N = mu.shape
    if N <= 1:
        return np.zeros((M, 0))

    mushifted = np.zeros((M, N - 1))
    for m in range(M):
        if m == i:
            mushifted[m, :] = mu[m, 1:N]
        else:
            mushifted[m, :] = mu[m, :N - 1]

    return mushifted


def _pfqn_comomrm_ld_wrapper(L: np.ndarray, N: np.ndarray, Z: np.ndarray,
                              mu: np.ndarray) -> Tuple[np.ndarray, float]:
    """
    Wrapper for pfqn_comomrm_ld that handles edge cases.

    Args:
        L: Load matrix
        N: Population vector
        Z: Think time vector
        mu: Service rate matrix

    Returns:
        Tuple of (result, log_normalizing_constant)
    """
    try:
        result = pfqn_comomrm_ld(L, N, Z, mu)
        if hasattr(result, 'lG'):
            return result.G, result.lG
        else:
            # Handle tuple return
            return result[0], result[1]
    except Exception:
        # Fallback for numerical issues
        return np.zeros(len(N)), 0.0


def pfqn_stdf_heur(L: np.ndarray, N: np.ndarray, Z: np.ndarray,
                    S: np.ndarray, fcfs_nodes: np.ndarray, rates: np.ndarray,
                    tset: np.ndarray) -> Dict[Tuple[int, int], np.ndarray]:
    """
    Heuristic sojourn time distribution for multiserver FCFS nodes.

    Implements a variant of McKenna's 1987 method that uses per-class
    service rate weighting for improved accuracy in multiclass networks.

    Unlike pfqn_stdf which uses a simpler Erlang model for the waiting
    component, this heuristic accounts for class-dependent waiting times
    based on the expected queue composition.

    Args:
        L: Load matrix (M x R) - service demands
        N: Population vector (R,) - number of jobs per class
        Z: Think time vector (R,) or matrix (1 x R)
        S: Number of servers at each station (M,)
        fcfs_nodes: Array of FCFS station indices (0-indexed)
        rates: Service rates matrix (M x R)
        tset: Time points at which to evaluate the distribution

    Returns:
        Dictionary with (station, class) tuples as keys and
        response time distribution arrays as values.
        Each array has shape (len(tset), 2) with:
            - column 0: CDF values
            - column 1: time points

    References:
        McKenna, J. "Mean Value Analysis for networks with service-time-dependent
        queue disciplines." Performance Evaluation, 1987.
    """
    from .rd import pfqn_rd

    L = np.atleast_2d(np.asarray(L, dtype=float))
    N = np.asarray(N, dtype=int).flatten()
    Z = np.atleast_1d(np.asarray(Z, dtype=float)).flatten()
    S = np.asarray(S, dtype=int).flatten()
    fcfs_nodes = np.asarray(fcfs_nodes, dtype=int).flatten()
    rates = np.atleast_2d(np.asarray(rates, dtype=float))
    tset = np.asarray(tset, dtype=float).flatten()

    M, R = L.shape
    T = len(tset)

    FINE_TOL = 1e-10
    tset = np.where(tset == 0, FINE_TOL, tset)

    # Build mu matrix
    mu = np.zeros((M, int(np.sum(N))))
    for k in range(M):
        for n in range(int(np.sum(N))):
            mu[k, n] = min(S[k], n + 1)

    stability_warn_issued = False
    RD = {}

    for k in fcfs_nodes:
        for r in range(R):
            if L[k, r] > FINE_TOL:
                Nr = oner(N, r)

                # Compute normalizing constant and queue lengths
                if M == 1:
                    _, lGr = _pfqn_comomrm_ld_wrapper(L, Nr, Z, mu)
                    # Compute Q1 using repairman formula
                    from .utils import pfqn_mu_ms
                    mu_ms = pfqn_mu_ms(int(np.sum(N)), 2, S[k])
                    mu_ms_mat = np.zeros((1, int(np.sum(N))))
                    mu_ms_mat[0, :len(mu_ms)] = mu_ms
                    _, lGra = _pfqn_comomrm_ld_wrapper(L, Nr, Z, mu_ms_mat)
                    Q1 = np.zeros((M, R))
                    for s in range(R):
                        Q1[k, s] = L[k, s] * np.exp(lGra - lGr)
                    is_num_stable = True
                else:
                    _, Q1, _, _, lGr, is_num_stable, _ = pfqn_mvald(L, Nr, Z, mu)

                lGr = lGr[-1] if hasattr(lGr, '__len__') else lGr

                if not is_num_stable and not stability_warn_issued:
                    stability_warn_issued = True

                # Build heuristic h_kc(t, n) with per-class weighting
                hkc = np.zeros((T, 1 + int(np.sum(N))))

                for n in range(int(np.sum(N)) + 1):
                    if n < S[k]:
                        # Exponential service (no queueing)
                        D0, D1 = map_exponential(1.0 / rates[k, r])
                    else:
                        # Service + weighted waiting in queue
                        # The heuristic: weight waiting by class-specific queue lengths
                        C = [map_exponential(1.0 / rates[k, r])]

                        Q1_k_total = np.sum(Q1[k, :])
                        if Q1_k_total > FINE_TOL:
                            for s in range(R):
                                if Q1[k, s] > FINE_TOL:
                                    # Mean waiting time component for class s
                                    mean_wait_s = Q1[k, s] * (n - S[k] + 1) / Q1_k_total / rates[k, s]
                                    C.append(map_exponential(mean_wait_s))
                        else:
                            # Fallback to uniform weighting
                            mean_wait = (n - S[k] + 1) / (S[k] * rates[k, r])
                            erlang = map_erlang(mean_wait, n - S[k] + 1)
                            C = [map_exponential(1.0 / rates[k, r]), erlang]

                        D0, D1 = map_sumind(C)

                    # Compute CDF at all time points
                    cdf_vals = map_cdf(D0, D1, tset)
                    hkc[:, n] = cdf_vals

                # Initialize result array
                RD[(k, r)] = np.zeros((T, 2))
                RD[(k, r)][:, 1] = tset

                # Compute lGk for network without station k
                other_stations = np.array([i for i in range(M) if i != k])
                if len(other_stations) > 0:
                    if len(other_stations) == 1:
                        _, lGk = _pfqn_comomrm_ld_wrapper(
                            L[other_stations, :], Nr, Z, mu[other_stations, :]
                        )
                    else:
                        _, _, _, _, lGk, _, _ = pfqn_mvald(
                            L[other_stations, :], Nr, Z, mu[other_stations, :]
                        )
                    lGk = lGk[-1] if hasattr(lGk, '__len__') else lGk
                else:
                    lGk = 0.0

                # Use recursive form for load-dependent models
                Hkrt = np.zeros(T)

                for t_idx in range(T):
                    gammat = mu.copy()

                    # Adjust gamma for time-dependent service
                    for m in range(1, int(np.sum(Nr)) + 1):
                        if hkc[t_idx, m] > FINE_TOL:
                            gammat[k, m - 1] = mu[k, m - 1] * hkc[t_idx, m - 1] / hkc[t_idx, m]

                    gammak = _mushift(gammat, k)
                    if gammak.shape[1] > 0:
                        gammak = gammak[:, :int(np.sum(Nr))]

                    # nvec = 0 case
                    Hkrt[t_idx] = hkc[t_idx, 0] * np.exp(lGk)

                    # nvec >= 1 cases
                    for s in range(R):
                        if Nr[s] > 0:
                            Nr_s = oner(Nr, s)

                            # Use pfqn_rd for the recursive computation
                            try:
                                lYks_t = pfqn_rd(L, Nr_s, Z, gammak)
                                lYks_t = lYks_t[-1] if hasattr(lYks_t, '__len__') else lYks_t
                            except Exception:
                                # Fallback to mvald
                                if M == 1:
                                    _, lYks_t = _pfqn_comomrm_ld_wrapper(
                                        L, Nr_s, Z, gammak
                                    )
                                else:
                                    _, _, _, _, lYks_t, _, _ = pfqn_mvald(
                                        L, Nr_s, Z, gammak
                                    )
                                lYks_t = lYks_t[-1] if hasattr(lYks_t, '__len__') else lYks_t

                            if gammat[k, 0] > FINE_TOL:
                                Hkrt[t_idx] += (L[k, s] * hkc[t_idx, 0] / gammat[k, 0]) * np.exp(lYks_t)

                # Handle NaN values
                Hkrt = np.where(np.isnan(Hkrt), FINE_TOL, Hkrt)
                lHkrt = np.log(np.maximum(Hkrt, FINE_TOL))

                # Compute final CDF
                RD[(k, r)][:, 0] = np.minimum(1.0, np.exp(lHkrt - lGr))

    return RD


__all__ = [
    'pfqn_stdf',
    'pfqn_stdf_heur',
]
