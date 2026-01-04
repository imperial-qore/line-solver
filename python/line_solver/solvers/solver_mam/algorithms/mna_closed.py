"""
Matrix-Normalizing Approximation for closed networks (mna_closed).

Extends mna_open with bisection search for lambda to maintain
population constraint in closed networks.

Algorithm:
1. Outer loop: bisection search on lambda
2. Inner loop: SCV iteration (same as mna_open)
3. Constraint: Σ(Q_n) = N (total population)
4. Extract distributions from MMAP departures

References:
    MATLAB: matlab/src/solvers/MAM/solver_mna_closed.m (346 lines)
"""

import numpy as np
import time
from typing import Tuple, Optional

from . import MAMAlgorithm, MAMResult
from ..utils.network_adapter import extract_mam_params, extract_visit_counts, check_closed_network


class MNAClosedAlgorithm(MAMAlgorithm):
    """Matrix-Normalizing Approximation for closed networks."""

    @staticmethod
    def supports_network(sn) -> Tuple[bool, Optional[str]]:
        """Check if network can be solved by mna_closed.

        MNA supports closed networks with arbitrary scheduling.

        Args:
            sn: NetworkStruct

        Returns:
            (can_solve, reason_if_not)
        """
        is_closed = check_closed_network(sn)
        if not is_closed:
            return False, "mna_closed is for closed networks. Use mna_open for open networks."

        return True, None

    def solve(self, sn, options=None) -> MAMResult:
        """Solve closed network using MNA with bisection.

        Args:
            sn: NetworkStruct
            options: Options dict with tol, max_iter, verbose

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
        njobs = params['njobs'] if 'njobs' in params else sn.njobs

        # Get options
        tol = getattr(options, 'tol', 1e-6) if options else 1e-6
        max_iter = getattr(options, 'max_iter', 100) if options else 100
        verbose = getattr(options, 'verbose', False) if options else False

        if verbose:
            print(f"mna_closed: M={M} stations, K={K} classes")

        # Initialize result matrices
        QN = np.zeros((M, K))
        UN = np.zeros((M, K))
        RN = np.zeros((M, K))
        TN = np.zeros((1, K))
        CN = np.zeros((1, K))
        XN = np.zeros((1, K))

        # Service times
        S = 1.0 / np.maximum(rates, 1e-10)

        # Compute total demand (sum of 1/service_rate across all stations)
        total_demand = np.sum(1.0 / np.maximum(rates, 1e-10), axis=0)

        # Initial lambda estimate (population / total demand)
        if len(njobs) > 0:
            N_total = np.sum(njobs)
        else:
            N_total = 1.0

        lambda_low = 1e-6
        lambda_high = N_total / np.max(total_demand + 1e-10)
        lambda_k = (lambda_low + lambda_high) / 2.0

        # Bisection loop for lambda
        bisect_iter = 0
        bisect_max = 20

        while bisect_iter < bisect_max:
            bisect_iter += 1

            # Solve inner SCV loop with current lambda
            QN, UN, RN = _solve_mna_closed_inner(
                M, K, lambda_k, S, scv_orig, nservers,
                max_iter, tol, verbose
            )

            # Check population constraint
            total_qlen = np.sum(QN)

            if verbose and bisect_iter <= 3:
                print(f"  Bisect {bisect_iter}: lambda={lambda_k:.6e}, total_Q={total_qlen:.4f}, target={N_total:.4f}")

            if abs(total_qlen - N_total) < tol * N_total or bisect_iter >= bisect_max:
                break

            # Adjust lambda based on queue length
            if total_qlen < N_total:
                lambda_low = lambda_k
            else:
                lambda_high = lambda_k

            lambda_k = (lambda_low + lambda_high) / 2.0

        # Compute TN and XN
        for k in range(K):
            TN[0, k] = lambda_k
            XN[0, k] = lambda_k

        # Compute cycle times
        for k in range(K):
            CN[0, k] = np.sum(RN[:, k]) / (lambda_k + 1e-10)

        runtime = time.time() - start_time

        return MAMResult(
            QN=QN,
            UN=UN,
            RN=RN,
            TN=TN,
            CN=CN,
            XN=XN,
            totiter=bisect_iter,
            method="mna_closed",
            runtime=runtime
        )


def _solve_mna_closed_inner(M: int,
                           K: int,
                           lambda_k: float,
                           S: np.ndarray,
                           scv_orig: np.ndarray,
                           nservers: np.ndarray,
                           max_iter: int,
                           tol: float,
                           verbose: bool) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Solve MNA closed network inner loop (SCV iteration).

    Same as mna_open but for a single lambda value (fixed).

    Args:
        M: Number of stations
        K: Number of classes
        lambda_k: Arrival rate per class (scalar, replicated)
        S: (M, K) service times
        scv_orig: (M, K) original SCVs
        nservers: (M,) number of servers
        max_iter: Max iterations
        tol: Tolerance
        verbose: Print progress

    Returns:
        (QN, UN, RN)
    """
    QN = np.zeros((M, K))
    UN = np.zeros((M, K))
    RN = np.zeros((M, K))

    scv = scv_orig.copy()
    scv_prev = scv.copy() + np.inf

    iteration = 0

    while np.max(np.abs(scv - scv_prev)) > tol and iteration < max_iter:
        iteration += 1
        scv_prev = scv.copy()

        # Solve each station
        for m in range(M):
            for k in range(K):
                rho = lambda_k * S[m, k] / nservers[m]
                UN[m, k] = rho

                if rho >= 1.0:
                    RN[m, k] = np.inf
                    QN[m, k] = np.inf
                else:
                    # QNA formula with SCV
                    ca2 = 1.0  # Assume Poisson internal arrivals
                    cs2 = scv[m, k]

                    if nservers[m] <= 1:
                        Wq = (ca2 + cs2) * rho / (2.0 * (1.0 - rho)) * S[m, k]
                    else:
                        Wq = rho / (1.0 - rho) * S[m, k]

                    RN[m, k] = S[m, k] + Wq
                    QN[m, k] = lambda_k * Wq

        # Update SCVs
        for m in range(M):
            for k in range(K):
                rho_m = UN[m, k]
                if rho_m < 0.99:
                    scv[m, k] = scv_orig[m, k] + (1.0 - rho_m) * (scv_orig[m, k] - 1.0)
                else:
                    scv[m, k] = 2.0 + scv_orig[m, k]

    return QN, UN, RN
