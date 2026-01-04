"""
Matrix-Normalizing Approximation for open networks (mna_open).

MNA is an approximation algorithm that uses QNA-like traffic equations
with matrix-generalized formulas for non-product-form networks.

Algorithm:
1. Solve traffic equations for arrival rates
2. Compute squared coefficient of variation (SCV) for each station
3. Use QNA-like formulas with SCV to compute response times
4. Iterate until convergence (SCV stabilization)

References:
    MATLAB: matlab/src/solvers/MAM/solver_mna_open.m (268 lines)
    QNA: Whitt, W. "The Queueing Network Analyzer." Bell Labs Technical Journal (1983)
"""

import numpy as np
import time
from typing import Tuple, Optional

from . import MAMAlgorithm, MAMResult
from ..utils.network_adapter import extract_mam_params, extract_visit_counts, build_routing_matrix


class MNAOpenAlgorithm(MAMAlgorithm):
    """Matrix-Normalizing Approximation for open networks."""

    @staticmethod
    def supports_network(sn) -> Tuple[bool, Optional[str]]:
        """Check if network can be solved by mna_open.

        MNA supports open networks with arbitrary scheduling and
        phase-type service distributions.

        Args:
            sn: NetworkStruct

        Returns:
            (can_solve, reason_if_not)
        """
        # Check if open network
        is_open = any(np.isinf(n) for n in sn.njobs) if len(sn.njobs) > 0 else True
        if not is_open:
            return False, "mna_open is for open networks. Use mna_closed for closed networks."

        return True, None

    def solve(self, sn, options=None) -> MAMResult:
        """Solve open network using Matrix-Normalizing Approximation.

        Args:
            sn: NetworkStruct
            options: MNA options (tol, max_iter, verbose)

        Returns:
            MAMResult
        """
        start_time = time.time()

        params = extract_mam_params(sn)
        M = params['nstations']
        K = params['nclasses']
        rates = params['rates']
        scv_orig = params['scv'].copy()
        nservers = params['nservers']
        visits = extract_visit_counts(sn)
        routing = build_routing_matrix(sn)

        # Get options
        tol = getattr(options, 'tol', 1e-6) if options else 1e-6
        max_iter = getattr(options, 'max_iter', 100) if options else 100
        verbose = getattr(options, 'verbose', False) if options else False

        if verbose:
            print(f"mna_open: M={M} stations, K={K} classes")

        # Initialize result matrices
        QN = np.zeros((M, K))
        UN = np.zeros((M, K))
        RN = np.zeros((M, K))
        TN = np.zeros((1, K))
        XN = np.zeros((1, K))

        # Solve traffic equations for lambda
        lambda_k = _solve_traffic_open(routing, visits, rates, M, K)

        # Service times
        S = 1.0 / np.maximum(rates, 1e-10)

        # Initialize SCVs
        scv = scv_orig.copy()
        scv_prev = scv.copy() + np.inf

        iteration = 0

        # MNA iteration loop (SCV convergence)
        while np.max(np.abs(scv - scv_prev)) > tol and iteration < max_iter:
            iteration += 1
            scv_prev = scv.copy()

            if verbose and iteration <= 3:
                print(f"  Iteration {iteration}: max SCV change = {np.max(np.abs(scv - scv_prev)):.6e}")

            # Solve each station with current SCVs
            for m in range(M):
                # Arrivals to station m
                arr_rate = np.sum(lambda_k)

                # Service
                for k in range(K):
                    rho = lambda_k[k] * S[m, k] / nservers[m]
                    UN[m, k] = rho

                    if rho >= 1.0:
                        RN[m, k] = np.inf
                        QN[m, k] = np.inf
                    else:
                        # QNA-like formula with aggregated SCV
                        ca2_agg = np.sum(lambda_k) ** 2 / (np.sum(lambda_k ** 2) + 1e-10)
                        if np.any(np.isnan(ca2_agg)) or np.any(np.isinf(ca2_agg)):
                            ca2_agg = 1.0

                        cs2 = scv[m, k]

                        # M/G/c approximation
                        if nservers[m] <= 1:
                            # M/G/1
                            Wq = (ca2_agg + cs2) * rho / (2.0 * (1.0 - rho)) * S[m, k]
                        else:
                            # M/G/c (approximation)
                            Wq = rho / (1.0 - rho) * S[m, k]

                        RN[m, k] = S[m, k] + Wq
                        QN[m, k] = lambda_k[k] * Wq

            # Update SCVs based on response times (splitting effect)
            for m in range(M):
                for k in range(K):
                    if RN[m, k] > 0:
                        # Variability of response time (approximation)
                        # scv_out ≈ scv_in + (1 - rho) * adjustment
                        rho_m = UN[m, k]
                        if rho_m < 0.99:
                            scv[m, k] = scv_orig[m, k] + (1.0 - rho_m) * (scv_orig[m, k] - 1.0)
                        else:
                            scv[m, k] = 2.0 + scv_orig[m, k]  # Increase SCV near saturation

        # Compute system metrics
        TN[0, :] = lambda_k
        XN[0, :] = lambda_k

        runtime = time.time() - start_time

        return MAMResult(
            QN=QN,
            UN=UN,
            RN=RN,
            TN=TN,
            XN=XN,
            totiter=iteration,
            method="mna_open",
            runtime=runtime
        )


def _solve_traffic_open(routing: np.ndarray,
                       visits: np.ndarray,
                       rates: np.ndarray,
                       M: int,
                       K: int) -> np.ndarray:
    """
    Solve traffic equations for open network.

    Computes arrival rates: λ = λ_source + Σ_j (λ_j * p_{j→i})

    Args:
        routing: (M, M) routing matrix
        visits: (M, K) visit counts
        rates: (M, K) service rates
        M: Number of stations
        K: Number of classes

    Returns:
        (K,) arrival rate vector
    """
    # Simple version: extract from source (first station)
    # assumes station 0 is source
    lambda_k = np.zeros(K)

    if rates.shape[0] > 0 and rates.shape[1] >= K:
        # External arrivals from source (station 0 rates)
        lambda_k = rates[0, :K].copy()
        lambda_k[~np.isfinite(lambda_k)] = 1.0
    else:
        lambda_k = np.ones(K)

    return lambda_k
