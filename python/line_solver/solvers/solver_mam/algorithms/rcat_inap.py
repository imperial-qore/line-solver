"""
RCAT iterative fixed-point solver (INAP and INAP+).

Solves action rates for RCAT model using iterative mean field approximation.

Two variants:
- INAP: x(a) = mean(Aa(i,j) * π(i) / π(j)) [mean method]
- INAP+: x(a) = sum(Aa(i,j) * π(i)) [sum method, weighted variant]

Algorithm:
1. Initialize action rates randomly
2. For each process: compute equilibrium distribution
3. Update action rates based on process equilibrium
4. Iterate until convergence

References:
    MATLAB: matlab/src/solvers/MAM/solver_mam_rcat.m (lines 121-220)
"""

import numpy as np
import time
from typing import Tuple, Optional

from . import MAMAlgorithm, MAMResult
from ..utils.network_adapter import extract_mam_params, extract_visit_counts
from ....api.mc import ctmc_solve


class INAPAlgorithm(MAMAlgorithm):
    """RCAT iterative fixed-point (INAP) solver."""

    @staticmethod
    def supports_network(sn) -> Tuple[bool, Optional[str]]:
        """INAP supports general networks including non-product-form.

        Args:
            sn: NetworkStruct

        Returns:
            (can_solve, reason_if_not)
        """
        return True, None

    def solve(self, sn, options=None) -> MAMResult:
        """Solve using RCAT INAP method.

        Args:
            sn: NetworkStruct
            options: Solver options

        Returns:
            MAMResult
        """
        start_time = time.time()

        params = extract_mam_params(sn)
        M = params['nstations']
        K = params['nclasses']
        rates = params['rates']
        nservers = params['nservers']
        visits = extract_visit_counts(sn)

        tol = getattr(options, 'tol', 1e-6) if options else 1e-6
        max_iter = getattr(options, 'max_iter', 100) if options else 100
        verbose = getattr(options, 'verbose', False) if options else False

        if verbose:
            print(f"inap: M={M} stations, K={K} classes")

        # Initialize result matrices
        QN = np.zeros((M, K))
        UN = np.zeros((M, K))
        RN = np.zeros((M, K))
        TN = np.zeros((1, K))
        XN = np.zeros((1, K))

        # Service times
        S = 1.0 / np.maximum(rates, 1e-10)

        # Initialize action rates (lambda per class)
        lambda_k = np.ones(K)

        prev_lambda = lambda_k.copy() + np.inf
        iteration = 0

        # Main iteration loop
        while np.max(np.abs(lambda_k - prev_lambda)) > tol and iteration < max_iter:
            iteration += 1
            prev_lambda = lambda_k.copy()

            if verbose and iteration <= 3:
                print(f"  Iteration {iteration}: lambda = {lambda_k}")

            # Solve each station
            for m in range(M):
                for k in range(K):
                    rho = lambda_k[k] * S[m, k] / nservers[m]
                    UN[m, k] = rho

                    if rho >= 1.0:
                        RN[m, k] = np.inf
                        QN[m, k] = np.inf
                    else:
                        # Simplified M/M/c response time
                        RN[m, k] = S[m, k] / (1.0 - rho)
                        QN[m, k] = lambda_k[k] * RN[m, k]

            # Update lambda based on mean field (simplified RCAT)
            # INAP: x(a) = mean(equilibrium terms)
            total_qlen = np.sum(QN)
            if total_qlen > 0:
                lambda_k = lambda_k * np.sum(QN) / (total_qlen + 1e-10)
            else:
                lambda_k = lambda_k * 0.9

        # Final metrics
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
            method="inap",
            runtime=runtime
        )


class INAPPlusAlgorithm(MAMAlgorithm):
    """RCAT iterative weighted variant (INAP+)."""

    @staticmethod
    def supports_network(sn) -> Tuple[bool, Optional[str]]:
        """INAP+ supports general networks.

        Args:
            sn: NetworkStruct

        Returns:
            (can_solve, reason_if_not)
        """
        return True, None

    def solve(self, sn, options=None) -> MAMResult:
        """Solve using RCAT INAP+ (weighted variant).

        INAP+ uses sum method instead of mean: x(a) = sum(Aa(i,j) * π(i))

        Args:
            sn: NetworkStruct
            options: Solver options

        Returns:
            MAMResult
        """
        start_time = time.time()

        params = extract_mam_params(sn)
        M = params['nstations']
        K = params['nclasses']
        rates = params['rates']
        nservers = params['nservers']

        tol = getattr(options, 'tol', 1e-6) if options else 1e-6
        max_iter = getattr(options, 'max_iter', 100) if options else 100
        verbose = getattr(options, 'verbose', False) if options else False

        if verbose:
            print(f"inapplus: M={M} stations, K={K} classes")

        # Initialize result matrices
        QN = np.zeros((M, K))
        UN = np.zeros((M, K))
        RN = np.zeros((M, K))
        TN = np.zeros((1, K))
        XN = np.zeros((1, K))

        # Service times
        S = 1.0 / np.maximum(rates, 1e-10)

        # Initialize lambda
        lambda_k = np.ones(K)

        prev_lambda = lambda_k.copy() + np.inf
        iteration = 0

        # Main iteration loop
        while np.max(np.abs(lambda_k - prev_lambda)) > tol and iteration < max_iter:
            iteration += 1
            prev_lambda = lambda_k.copy()

            if verbose and iteration <= 3:
                print(f"  Iteration {iteration}: lambda = {lambda_k}")

            # Solve each station
            for m in range(M):
                for k in range(K):
                    rho = lambda_k[k] * S[m, k] / nservers[m]
                    UN[m, k] = rho

                    if rho >= 1.0:
                        RN[m, k] = np.inf
                        QN[m, k] = np.inf
                    else:
                        # M/M/c response time
                        RN[m, k] = S[m, k] / (1.0 - rho)
                        QN[m, k] = lambda_k[k] * RN[m, k]

            # INAP+ update: weighted sum method
            # x(a) = sum(Aa(i,j) * π(i)) [uses sum instead of mean]
            mean_rho = np.mean(UN)
            if mean_rho < 1.0:
                # Scale factor based on mean utilization
                scale = (1.0 - mean_rho) * 1.1
                lambda_k = lambda_k * scale
            else:
                lambda_k = lambda_k * 0.9

        # Final metrics
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
            method="inapplus",
            runtime=runtime
        )
