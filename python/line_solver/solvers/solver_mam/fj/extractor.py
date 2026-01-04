"""
Fork-Join parameter extraction.

Extracts arrival and service parameters from a Fork-Join network
in a format suitable for FJ_codes percentile analysis.
"""

import numpy as np
from typing import Tuple, Optional
from dataclasses import dataclass

from ....api.sn.network_struct import NetworkStruct
from .validator import fj_isfj, classify_distribution


@dataclass
class FJParams:
    """Extracted Fork-Join network parameters."""
    K: int  # Number of parallel queues
    lambda_arr: np.ndarray  # Arrival rates (shape: (nclasses,))
    service_rates: np.ndarray  # Service rates (shape: (K, nclasses))
    service_scv: np.ndarray  # Service SCV (shape: (K, nclasses))
    nclasses: int
    lambda_total: float  # Total arrival rate
    service_type: str  # 'M', 'E', 'H', 'G', 'MAP'
    is_homogeneous: bool  # All queues identical


def fj_extract_params(sn: NetworkStruct) -> Optional[FJParams]:
    """
    Extract Fork-Join parameters from NetworkStruct.

    Args:
        sn: Compiled NetworkStruct (must be valid Fork-Join topology)

    Returns:
        FJParams if valid Fork-Join, None otherwise
    """
    # Validate Fork-Join topology
    result = fj_isfj(sn)
    if not result.is_fj:
        return None

    K = result.K
    queue_indices = result.queue_indices
    nclasses = sn.nclasses

    # Extract arrival rates
    lambda_arr = np.ones(nclasses)
    if hasattr(sn, 'lambda_arr') and sn.lambda_arr is not None:
        sn_lambda = np.asarray(sn.lambda_arr, dtype=np.float64).flatten()
        for k in range(nclasses):
            if k < len(sn_lambda) and sn_lambda[k] > 0:
                lambda_arr[k] = sn_lambda[k]

    lambda_total = np.sum(lambda_arr)

    # Extract service rates and SCV for parallel queues
    service_rates = np.asarray(sn.rates, dtype=np.float64)
    service_scv = np.asarray(sn.scv, dtype=np.float64)

    if service_rates.ndim == 1:
        service_rates = service_rates.reshape(-1, 1)
    if service_scv.ndim == 1:
        service_scv = service_scv.reshape(-1, 1)

    # Extract only the parallel queue rows
    parallel_rates = np.zeros((K, nclasses))
    parallel_scv = np.zeros((K, nclasses))

    for k_idx, station_idx in enumerate(queue_indices):
        if station_idx < service_rates.shape[0]:
            parallel_rates[k_idx, :] = service_rates[station_idx, :nclasses]
        if station_idx < service_scv.shape[0]:
            parallel_scv[k_idx, :] = service_scv[station_idx, :nclasses]

    # Determine service type (for single-class case)
    service_type = 'G'
    if nclasses == 1:
        scv = parallel_scv[0, 0] if parallel_scv.shape[0] > 0 else 1.0
        service_type = classify_distribution(scv, parallel_rates[0, 0] if parallel_rates.shape[0] > 0 else 1.0)

    # Check homogeneity
    is_homogeneous = True
    if K > 1:
        ref_rate = parallel_rates[0, :]
        ref_scv = parallel_scv[0, :]
        for k in range(1, K):
            if not (np.allclose(parallel_rates[k, :], ref_rate, rtol=1e-6) and
                    np.allclose(parallel_scv[k, :], ref_scv, rtol=1e-6)):
                is_homogeneous = False
                break

    return FJParams(
        K=K,
        lambda_arr=lambda_arr,
        service_rates=parallel_rates,
        service_scv=parallel_scv,
        nclasses=nclasses,
        lambda_total=lambda_total,
        service_type=service_type,
        is_homogeneous=is_homogeneous
    )


def fj_convert_to_fj_codes_format(params: FJParams) -> dict:
    """
    Convert FJParams to FJ_codes input format.

    FJ_codes expects:
    - K: number of parallel queues
    - L: arrival MAP (D0, D1)
    - S: service distribution representations
    - SerChoice: service type code

    Args:
        params: Extracted Fork-Join parameters

    Returns:
        Dictionary with FJ_codes format parameters
    """
    # For exponential arrivals, create MAP
    lambda_total = params.lambda_total
    D0 = np.array([[-lambda_total]])
    D1 = np.array([[lambda_total]])

    # Service distribution selection
    # SerChoice: 1=M/M/c, 2=Erlang, 3=HyperExp, 4=General PH
    ser_choice = 1  # Default: Exponential
    if params.service_type == 'E':
        ser_choice = 2
    elif params.service_type == 'H':
        ser_choice = 3
    elif params.service_type == 'G':
        ser_choice = 4

    return {
        'K': params.K,
        'lambda_total': params.lambda_total,
        'D0': D0,
        'D1': D1,
        'service_rates': params.service_rates,
        'service_scv': params.service_scv,
        'service_type': params.service_type,
        'ser_choice': ser_choice,
        'nclasses': params.nclasses,
    }
