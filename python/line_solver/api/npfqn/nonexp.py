"""
Non-exponential approximation methods for NPFQN.

Native Python implementations of approximation methods for non-product-form
queueing networks with non-exponential service and inter-arrival time
distributions.

Key functions:
    npfqn_nonexp_approx: Main non-exponential approximation method

References:
    Casale, G., et al. "LINE: A unified library for queueing network modeling."
"""

import numpy as np
from math import exp, pow
from typing import Optional, Dict, Any
from dataclasses import dataclass

from ..sn import NetworkStruct, SchedStrategy


# Fine tolerance for numerical comparisons
FINE_TOL = 1e-12


@dataclass
class NpfqnNonexpApproxResult:
    """Result of non-exponential approximation."""
    ST: np.ndarray       # Updated service times
    gamma: np.ndarray    # Correction factors
    nservers: np.ndarray # Updated server counts
    rho: np.ndarray      # Utilizations
    scva: np.ndarray     # SCV of arrivals
    scvs: np.ndarray     # SCV of service
    eta: np.ndarray      # Eta factors


def npfqn_nonexp_approx(
    method: str,
    sn: NetworkStruct,
    ST: np.ndarray,
    V: Optional[np.ndarray],
    SCV: np.ndarray,
    Tin: np.ndarray,
    Uin: np.ndarray,
    gamma: np.ndarray,
    nservers: np.ndarray
) -> NpfqnNonexpApproxResult:
    """
    Approximates non-product-form queueing networks using the specified method.

    This function adjusts service times and other parameters to account for
    non-exponential service time distributions in product-form analysis.

    Args:
        method: Approximation method ("default", "none", "hmva", "interp")
        sn: Network structure
        ST: Service time matrix (M x K)
        V: Visit ratios matrix (M x K), optional
        SCV: Squared coefficient of variation matrix (M x K)
        Tin: Initial throughput matrix (M x K)
        Uin: Initial utilization matrix (M x K)
        gamma: Gamma correction matrix (M x 1)
        nservers: Number of servers matrix (M x 1)

    Returns:
        NpfqnNonexpApproxResult with updated matrices

    Raises:
        ValueError: If unknown approximation method is specified
    """
    M = sn.nstations

    # Initialize output arrays
    rho = np.zeros(M)
    scva = np.ones(M)
    scvs = np.ones(M)
    eta = np.ones(M)

    # Create local copies to avoid modifying originals
    T = Tin.copy()
    U = Uin.copy()
    ST_local = ST.copy()
    gamma_local = gamma.copy() if gamma.ndim > 0 else np.zeros(M)
    nservers_local = nservers.copy() if nservers.ndim > 0 else np.ones(M)

    # Ensure arrays are 1D for indexing
    if gamma_local.ndim > 1:
        gamma_local = gamma_local.flatten()
    if nservers_local.ndim > 1:
        nservers_local = nservers_local.flatten()

    if method in ["default", "none", "hmva"]:
        # No-op - just return initialized values
        return NpfqnNonexpApproxResult(
            ST=ST_local,
            gamma=gamma_local,
            nservers=nservers_local,
            rho=rho,
            scva=scva,
            scvs=scvs,
            eta=eta
        )

    elif method == "interp":
        # Interpolation method for non-exponential service times
        sched = sn.sched
        stations = sn.stations if hasattr(sn, 'stations') and sn.stations is not None else list(range(M))

        for ist in range(M):
            # Find non-zero classes for this station
            K = ST_local.shape[1]
            nnz_classes = np.zeros(K)

            for j in range(K):
                if np.isfinite(ST_local[ist, j]) and np.isfinite(SCV[ist, j]):
                    nnz_classes[j] = 1

            # Compute utilization for non-zero classes
            for j in range(K):
                if nnz_classes[j] > 0:
                    rho[ist] += U[ist, j]

            # Apply correction for FCFS stations with non-exponential service
            if np.sum(nnz_classes) > 0:
                # Get scheduling strategy for this station
                station_key = stations[ist] if ist < len(stations) else ist
                station_sched = sched.get(station_key, None) if sched else None

                if station_sched == SchedStrategy.FCFS:
                    # Extract non-zero class data
                    nnz_indices = np.where(nnz_classes > 0)[0]
                    n_nnz = len(nnz_indices)

                    ST_nnz = np.array([ST_local[ist, j] for j in nnz_indices])
                    SCV_nnz = np.array([SCV[ist, j] for j in nnz_indices])
                    T_nnz = np.array([T[ist, j] for j in nnz_indices])

                    # Check if correction is needed
                    # (heterogeneous service times or non-exponential service)
                    if (np.max(ST_nnz) - np.min(ST_nnz) > 0 or
                        np.max(SCV_nnz) > 1 + FINE_TOL or
                        np.min(SCV_nnz) < 1 - FINE_TOL):

                        scva[ist] = 1.0

                        # Weighted average SCV of service
                        T_sum = np.sum(T_nnz)
                        if T_sum > 0:
                            scvs[ist] = np.dot(SCV_nnz, T_nnz) / T_sum
                        else:
                            scvs[ist] = 1.0

                        # Compute gamma correction factor
                        gamma_local[ist] = (pow(rho[ist], nservers_local[ist]) + rho[ist]) / 2

                        # Compute eta factor
                        if (scvs[ist] > 1 - 1e-6 and scvs[ist] < 1 + 1e-6 and
                            nservers_local[ist] == 1.0):
                            eta[ist] = rho[ist]
                        else:
                            denom = scvs[ist] + scva[ist] * rho[ist]
                            if denom > 0:
                                eta[ist] = exp(-2 * (1 - rho[ist]) / denom)
                            else:
                                eta[ist] = 1.0

                        # Apply service time correction
                        order = 8
                        ai = pow(rho[ist], order)
                        bi = pow(rho[ist], order)

                        for k in nnz_indices:
                            rates = sn.rates
                            if rates is not None and ist < rates.shape[0] and k < rates.shape[1]:
                                if rates[ist, k] > 0:
                                    # Interpolate between original ST and adjusted ST
                                    term1 = max(0, 1 - ai) * ST_local[ist, k]
                                    term2 = ai * (bi * eta[ist] + max(0, 1 - bi) * gamma_local[ist])
                                    if T_sum > 0:
                                        term2 *= nservers_local[ist] / T_sum
                                    ST_local[ist, k] = term1 + term2

                        # Adjust server count
                        nservers_local[ist] = 1.0

        return NpfqnNonexpApproxResult(
            ST=ST_local,
            gamma=gamma_local,
            nservers=nservers_local,
            rho=rho,
            scva=scva,
            scvs=scvs,
            eta=eta
        )

    else:
        raise ValueError(f"Unknown approximation method: {method}")


__all__ = [
    'npfqn_nonexp_approx',
    'NpfqnNonexpApproxResult',
]
