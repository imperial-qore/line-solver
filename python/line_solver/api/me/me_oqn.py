"""
Maximum Entropy Methods for Open Queueing Networks.

Implements the ME algorithm from Kouvatsos (1994) "Entropy Maximisation
and Queueing Network Models", Section 3.2.
"""

import numpy as np
from typing import Tuple, Optional, Dict, Any
import warnings


def me_oqn(
    M: int,
    R: int,
    lambda0: np.ndarray,
    Ca0: np.ndarray,
    mu: np.ndarray,
    Cs: np.ndarray,
    P: np.ndarray,
    options: Optional[Dict[str, Any]] = None
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Maximum Entropy algorithm for Open Queueing Networks.

    Implements the ME algorithm from Kouvatsos (1994) for analyzing
    open queueing networks with general arrival and service processes.

    Args:
        M: Number of queues (stations).
        R: Number of job classes.
        lambda0: External arrival rates [M x R matrix], lambda0[i,r] = λ_oi,r.
        Ca0: External arrival squared coefficient of variation [M x R matrix].
        mu: Service rates [M x R matrix], mu[i,r] = μ_i,r.
        Cs: Service squared coefficient of variation [M x R matrix].
        P: Routing probability matrix [M x M x R], P[j,i,r] = p_ji,r
           (probability class r goes from queue j to queue i).
        options: Optional dict with fields:
            - tol: convergence tolerance (default: 1e-6)
            - maxiter: maximum iterations (default: 1000)
            - verbose: print iteration info (default: False)

    Returns:
        Tuple of (L, W, Ca, Cd, lambda_arr, rho):
            L: Mean queue lengths [M x R matrix].
            W: Mean waiting times [M x R matrix].
            Ca: Arrival scv at each queue [M x R matrix].
            Cd: Departure scv at each queue [M x R matrix].
            lambda_arr: Total arrival rates [M x R matrix].
            rho: Utilizations [M x R matrix].

    Reference:
        Kouvatsos (1994), Equations 3.6 and 3.7.
    """
    # Handle optional arguments
    if options is None:
        options = {}
    tol = options.get('tol', 1e-6)
    maxiter = options.get('maxiter', 1000)
    verbose = options.get('verbose', False)

    # Ensure arrays
    lambda0 = np.atleast_2d(np.asarray(lambda0, dtype=float))
    Ca0 = np.atleast_2d(np.asarray(Ca0, dtype=float))
    mu = np.atleast_2d(np.asarray(mu, dtype=float))
    Cs = np.atleast_2d(np.asarray(Cs, dtype=float))
    P = np.asarray(P, dtype=float)

    # Step 1: Feedback correction
    # If p_ii,r > 0, apply feedback elimination transformation
    P_eff = P.copy()
    mu_eff = mu.copy()
    Cs_eff = Cs.copy()

    for i in range(M):
        for r in range(R):
            pii = P[i, i, r]
            if pii > 0:
                # Feedback correction: adjust service rate and scv
                mu_eff[i, r] = mu[i, r] * (1 - pii)
                # Adjusted service scv accounting for feedback
                Cs_eff[i, r] = Cs[i, r] / (1 - pii) + pii / (1 - pii)
                # Remove self-loop from routing
                P_eff[i, i, r] = 0
                # Renormalize routing probabilities
                row_sum = np.sum(P_eff[i, :, r])
                if row_sum > 0:
                    P_eff[i, :, r] = P_eff[i, :, r] / row_sum * (1 - pii)

    # Step 2: Initialize arrival scv
    Ca = np.ones((M, R))

    # Step 3: Solve job flow balance equations using ORIGINAL P (including feedback)
    # λ_i,r = λ_oi,r + Σ_j λ_j,r * p_ji,r
    lambda_arr = np.zeros((M, R))
    for r in range(R):
        # Build the system (I - P') * lambda = lambda0
        Pr = P[:, :, r].T  # P[j,i,r] -> need sum over j
        A = np.eye(M) - Pr
        lambda_arr[:, r] = np.linalg.solve(A, lambda0[:, r])

    # Compute utilizations using ORIGINAL mu
    rho = np.zeros((M, R))
    for i in range(M):
        for r in range(R):
            if mu[i, r] > 0:
                rho[i, r] = lambda_arr[i, r] / mu[i, r]

    # Check stability
    rho_total = np.sum(rho, axis=1)
    if np.any(rho_total >= 1):
        warnings.warn("Network is unstable (utilization >= 1 at some queues)",
                      RuntimeWarning)

    # Initialize outputs
    L = np.zeros((M, R))
    Cd = np.ones((M, R))

    # Step 4-5: Iterative computation of Ca, Cd, and L
    for iteration in range(maxiter):
        Ca_old = Ca.copy()

        # Step 4: Apply GE-type formulae for mean queue length
        for i in range(M):
            rho_i = np.sum(rho[i, :])
            if rho_i > 0 and rho_i < 1:
                for r in range(R):
                    if rho[i, r] > 0:
                        # Proportion of class r traffic
                        prop_r = rho[i, r] / rho_i
                        # GE/GE/1 mean queue length formula (ME approximation)
                        Ca_agg = Ca[i, r]
                        Cs_agg = Cs_eff[i, r]
                        L_total = rho_i + (rho_i ** 2 * (Ca_agg + Cs_agg)) / (2 * (1 - rho_i))
                        L[i, r] = prop_r * L_total

        # Step 5a: Compute departure scv using equation (3.6)
        for j in range(M):
            rho_j = np.sum(rho[j, :])
            for r in range(R):
                if lambda_arr[j, r] > 0:
                    Cd[j, r] = 2 * L[j, r] * (1 - rho_j) + Ca[j, r] * (1 - 2 * rho_j)
                    Cd[j, r] = max(0, Cd[j, r])

        # Step 5b: Compute arrival scv using equation (3.7)
        for i in range(M):
            for r in range(R):
                if lambda_arr[i, r] > 0:
                    sum_inv = 0.0

                    # Contribution from other queues
                    for j in range(M):
                        pji = P_eff[j, i, r]
                        if pji > 0 and lambda_arr[j, r] > 0:
                            # Thinning formula for departure scv after splitting
                            Cdji = 1 + pji * (Cd[j, r] - 1)
                            weight = (lambda_arr[j, r] * pji) / lambda_arr[i, r]
                            sum_inv += weight / (Cdji + 1)

                    # Contribution from external arrivals
                    if lambda0[i, r] > 0:
                        weight0 = lambda0[i, r] / lambda_arr[i, r]
                        sum_inv += weight0 / (Ca0[i, r] + 1)

                    # Compute new arrival scv
                    if sum_inv > 0:
                        Ca[i, r] = -1 + 1 / sum_inv
                        Ca[i, r] = max(0, Ca[i, r])

        # Check convergence
        delta = np.max(np.abs(Ca - Ca_old))
        if verbose:
            print(f"Iteration {iteration + 1}: max delta = {delta:e}")
        if delta < tol:
            if verbose:
                print(f"Converged after {iteration + 1} iterations")
            break

    if iteration == maxiter - 1 and delta >= tol:
        warnings.warn(
            f"Did not converge within {maxiter} iterations (delta={delta:e})",
            RuntimeWarning
        )

    # Step 6: Compute final statistics
    # Mean waiting time via Little's law: W = L / λ
    W = np.zeros((M, R))
    for i in range(M):
        for r in range(R):
            if lambda_arr[i, r] > 0:
                W[i, r] = L[i, r] / lambda_arr[i, r]

    return L, W, Ca, Cd, lambda_arr, rho


__all__ = ['me_oqn']
