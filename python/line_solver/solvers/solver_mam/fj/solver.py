"""
Fork-Join solver integration for SolverMAM.

Detects Fork-Join topologies and applies specialized analysis using
FJ_codes for accurate percentile response time computation.
"""

import numpy as np
import time
from typing import Optional, List, Tuple

from ....api.sn.network_struct import NetworkStruct
from ....api.fj.fj_codes import (
    fj_percentiles_fork_join,
    fj_percentiles_general,
    FJPercentileResult,
)
from .validator import fj_isfj
from .extractor import fj_extract_params


class FJSolver:
    """Specialized solver for Fork-Join networks."""

    def __init__(self, verbose: bool = False):
        """Initialize FJ solver."""
        self.verbose = verbose

    def can_solve(self, sn: NetworkStruct) -> Tuple[bool, Optional[str]]:
        """
        Check if network has valid Fork-Join topology.

        Args:
            sn: NetworkStruct to analyze

        Returns:
            (can_solve, reason_if_not)
        """
        result = fj_isfj(sn)
        return result.is_fj, result.reason

    def compute_percentiles(self, sn: NetworkStruct, percentiles: List[float],
                          mean_response_time: float) -> Optional[FJPercentileResult]:
        """
        Compute response time percentiles for Fork-Join network.

        Args:
            sn: NetworkStruct (must be valid Fork-Join)
            percentiles: List of percentiles to compute (0-100)
            mean_response_time: Mean response time (optional, for validation)

        Returns:
            FJPercentileResult with percentile response times, or None if not FJ
        """
        # Validate topology
        result = fj_isfj(sn)
        if not result.is_fj:
            return None

        # Extract parameters
        params = fj_extract_params(sn)
        if params is None:
            return None

        K = params.K
        lambda_total = params.lambda_total
        service_rates = params.service_rates
        service_scv = params.service_scv

        if self.verbose:
            print(f"FJ Solver: K={K} parallel queues")
            print(f"  λ_total={lambda_total:.4f}")
            print(f"  Service rates: {service_rates[:, 0] if service_rates.ndim > 1 else service_rates}")
            print(f"  Service SCV: {service_scv[:, 0] if service_scv.ndim > 1 else service_scv}")

        # Check if all queues are exponential (M/M/c case)
        is_exponential = True
        if service_scv.size > 0:
            scv_vals = service_scv.flatten()
            for scv in scv_vals:
                if not np.isclose(scv, 1.0, rtol=0.1):
                    is_exponential = False
                    break

        # Compute percentiles
        if is_exponential:
            # Use optimized M/M/c analysis
            mu = service_rates[0, 0] if service_rates.ndim > 1 else service_rates[0]
            fj_result = fj_percentiles_fork_join(K, lambda_total, mu, percentiles)
        else:
            # Use general distribution analysis
            fj_result = fj_percentiles_general(K, lambda_total, service_rates, service_scv, percentiles)

        return fj_result

    def integrate_into_result(self, fj_result: FJPercentileResult) -> dict:
        """
        Convert FJ percentile result to format compatible with SolverMAM.

        Args:
            fj_result: FJPercentileResult from percentile computation

        Returns:
            Dictionary with percentile response times per class
        """
        # For multi-class, replicate percentiles for each class
        # For now, assume single-class or that percentiles apply to all classes

        percentile_dict = {
            'method': 'fj',
            'percentiles': fj_result.percentiles,
            'response_times': {
                0: fj_result.response_times  # Assume class 0 for single-class networks
            },
            'mean_response_time': fj_result.mean_response_time,
            'K': fj_result.K,
        }

        return percentile_dict


def solve_fork_join(sn: NetworkStruct, percentiles: Optional[List[float]] = None,
                   verbose: bool = False) -> Optional[dict]:
    """
    Solve Fork-Join network and compute percentiles.

    Convenience function for direct Fork-Join analysis without full SolverMAM.

    Args:
        sn: NetworkStruct (must be valid Fork-Join)
        percentiles: List of percentiles to compute (default: [50, 90, 95, 99])
        verbose: Enable verbose output

    Returns:
        Dictionary with results, or None if not a valid Fork-Join network
    """
    if percentiles is None:
        percentiles = [50, 90, 95, 99]

    solver = FJSolver(verbose=verbose)

    # Validate
    can_solve, reason = solver.can_solve(sn)
    if not can_solve:
        if verbose:
            print(f"Not a Fork-Join network: {reason}")
        return None

    # Compute metrics (use decomposition method for mean response time)
    # For now, approximate mean based on per-queue analysis
    from ..algorithms.dec_source import DecSourceAlgorithm
    from ..algorithms import MAMResult

    try:
        dec_solver = DecSourceAlgorithm()
        basic_result = dec_solver.solve(sn)

        mean_rt_system = np.sum(basic_result.RN[:, 0]) if basic_result.RN.size > 0 else 1.0
    except:
        mean_rt_system = 1.0  # Fallback

    # Compute percentiles
    fj_result = solver.compute_percentiles(sn, percentiles, mean_rt_system)
    if fj_result is None:
        return None

    # Format result
    result = solver.integrate_into_result(fj_result)
    return result
