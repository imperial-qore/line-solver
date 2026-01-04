"""
Age of Information (AoI) analysis method for SolverFLD.

Integrates aoi-fluid solvers into the FLD solver framework for computing AoI
and PAoI distributions in single-queue systems.
"""

import numpy as np
import time
from typing import Optional, Tuple, Dict, Any
from ..options import SolverFLDOptions, FLDResult
from line_solver.api.sn import NetworkStruct
from line_solver.api.aoi import (
    aoi_is_aoi,
    aoi_extract_params,
    solve_bufferless,
    solve_singlebuffer,
)


class AoIMethod:
    """Age of Information solver method for single-queue systems."""

    def __init__(self, sn: NetworkStruct, options: Optional[SolverFLDOptions] = None):
        """
        Initialize AoI method.

        Parameters
        ----------
        sn : NetworkStruct
            Network structure
        options : SolverFLDOptions, optional
            Solver options (may contain aoi_preemption override)
        """
        self.sn = sn
        self.options = options or SolverFLDOptions()

    def solve(self) -> FLDResult:
        """
        Solve AoI analysis.

        Workflow:
        1. Validate network topology (aoi_is_aoi)
        2. Extract parameters (aoi_extract_params)
        3. Solve using appropriate solver (bufferless or singlebuffer)
        4. Compute standard FLD metrics (QN, UN, RN, TN)
        5. Package results in FLDResult

        Returns
        -------
        result : FLDResult
            Result container with:
            - Standard metrics: QN, UN, RN, TN
            - AoI-specific: aoiResults dict with AoI/PAoI parameters
        """
        start_time = time.time()

        try:
            # Step 1: Validate topology
            is_valid, aoi_info = aoi_is_aoi(self.sn)
            if not is_valid:
                raise ValueError(f"Invalid AoI topology: {aoi_info.get('reason', 'Unknown error')}")

            # Step 2: Extract parameters
            aoi_params, solver_config = aoi_extract_params(self.sn, aoi_info, self.options)

            # Step 3: Solve
            if aoi_info['systemType'] == 'bufferless':
                aoi_result = solve_bufferless(**aoi_params)
            else:
                aoi_result = solve_singlebuffer(**aoi_params)

            if aoi_result['status'] != 'success':
                raise RuntimeError(f"AoI solver failed: {aoi_result.get('error_message', 'Unknown')}")

            # Step 4: Compute standard FLD metrics
            queue_idx = aoi_info['queue_idx']
            lambda_arr = self.sn.lambda_arr[0] if hasattr(self.sn, 'lambda_arr') else 0.5
            mu = self.sn.rates[queue_idx, 0] if (hasattr(self.sn, 'rates') and self.sn.rates is not None) else 1.0
            rho = lambda_arr / mu if mu > 0 else 0.0

            # M/M/1 approximation for QN, UN, RN, TN
            if rho < 1.0:
                QN = np.array([[rho / (1.0 - rho)]])  # Queue length
                UN = np.array([[rho]])  # Utilization
                RN = np.array([[1.0 / (mu - lambda_arr)]])  # Response time
                TN = np.array([[lambda_arr]])  # Throughput
            else:
                # Unstable: return inf
                QN = np.array([[np.inf]])
                UN = np.array([[np.inf]])
                RN = np.array([[np.inf]])
                TN = np.array([[lambda_arr]])

            CN = np.array([[RN[0, 0]]])  # Cycle time (single station)
            XN = np.array([[lambda_arr]])  # System throughput

            # Step 5: Package results
            runtime = time.time() - start_time

            result = FLDResult(
                QN=QN,
                UN=UN,
                RN=RN,
                TN=TN,
                CN=CN,
                XN=XN,
                t=None,
                QNt={},
                UNt={},
                xvec=None,
                iterations=1,
                runtime=runtime,
                method='aoi',
            )

            # Attach AoI-specific results
            result.aoiResults = aoi_result
            result.aoiConfig = solver_config

            return result

        except Exception as e:
            runtime = time.time() - start_time
            raise RuntimeError(f"AoI method failed after {runtime:.3f}s: {str(e)}")


def run_aoi_analysis(sn: NetworkStruct, options: Optional[SolverFLDOptions] = None) -> FLDResult:
    """
    Run AoI analysis on network.

    Convenience function that creates AoIMethod and solves.

    Parameters
    ----------
    sn : NetworkStruct
        Network structure
    options : SolverFLDOptions, optional
        Solver options

    Returns
    -------
    result : FLDResult
        Solution with AoI metrics attached
    """
    method = AoIMethod(sn, options)
    return method.solve()
