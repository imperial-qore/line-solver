"""
Age of Information (AoI) analysis method for SolverFLD.

Integrates aoi-fluid solvers into the FLD solver framework for computing AoI
and PAoI distributions in single-queue systems.
"""

import time
from typing import Optional

import numpy as np

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
            source_station = int(aoi_info['sourceStation'])
            queue_station = int(aoi_info['queueStation'])
            class_idx = int(aoi_info['openClass'])
            m = int(self.sn.nstations)
            k = int(self.sn.nclasses)

            lambda_arr = self._get_rate(source_station, class_idx, fallback=0.5)
            mean_service_time = float(
                np.real_if_close(
                    -aoi_params['sigma'] @ np.linalg.inv(aoi_params['S']) @ np.ones((aoi_params['S'].shape[0], 1))
                ).item()
            )
            mu = 1.0 / mean_service_time if mean_service_time > 0 else 0.0
            rho = lambda_arr / mu if mu > 0 else np.inf

            QN = np.zeros((m, k))
            UN = np.zeros((m, k))
            RN = np.zeros((m, k))
            TN = np.zeros((m, k))
            CN = np.zeros((1, k))
            XN = np.zeros((1, k))

            UN[queue_station, class_idx] = min(1.0, rho) if np.isfinite(rho) else 1.0
            TN[queue_station, class_idx] = lambda_arr if rho < 1.0 else mu
            if 0 <= source_station < m:
                TN[source_station, class_idx] = TN[queue_station, class_idx]
            if rho < 1.0:
                QN[queue_station, class_idx] = rho / (1.0 - rho)
                RN[queue_station, class_idx] = 1.0 / (mu - lambda_arr)
            else:
                QN[queue_station, class_idx] = np.inf
                RN[queue_station, class_idx] = np.inf
            CN[0, class_idx] = RN[queue_station, class_idx]
            XN[0, class_idx] = TN[queue_station, class_idx]

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
                aoiResults=aoi_result,
                aoiConfig=solver_config,
            )

            return result

        except Exception as e:
            runtime = time.time() - start_time
            raise RuntimeError(f"AoI method failed after {runtime:.3f}s: {str(e)}")

    def _get_rate(self, station: int, class_idx: int, fallback: float) -> float:
        if hasattr(self.sn, 'rates') and self.sn.rates is not None:
            rates = np.asarray(self.sn.rates, dtype=float)
            if rates.ndim == 2 and station < rates.shape[0] and class_idx < rates.shape[1]:
                return float(rates[station, class_idx])
        if hasattr(self.sn, 'lambda_arr') and self.sn.lambda_arr is not None:
            lambda_arr = np.asarray(self.sn.lambda_arr, dtype=float).reshape(-1)
            if class_idx < len(lambda_arr):
                return float(lambda_arr[class_idx])
        return float(fallback)


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
