"""
Maximum Entropy Methods for Open Queueing Networks.

Implements the ME algorithm from Kouvatsos (1994) "Entropy Maximisation
and Queueing Network Models", Section 3.2, with the GE/GE/c building block
of Section 3.4 (eq. 3.9) and the GE/GE/inf building block.
"""

import numpy as np
from typing import Tuple, Optional, Dict, Any
import warnings


def _ge_gec_mql(lam: float, Ca: float, mu: float, Cs: float, c: int) -> float:
    """Mean queue length of a stable GE/GE/c/FCFS queue via the exact ME
    solution of Kouvatsos (1994), eq. (3.9)."""
    c = int(round(c))
    alpha2 = 2.0 / (Cs + 1.0)
    alpha1 = 1.0 - alpha2
    beta2 = 2.0 / (Ca + 1.0)
    beta1 = 1.0 - beta2
    lambda2 = beta2 * lam
    mu2 = alpha2 * mu
    g = np.zeros(c)
    for j in range(1, c):
        g[j - 1] = (lambda2 + (j - 1) * mu2 * beta1) * alpha2 / (j * mu2 * (1.0 - alpha1 * beta1))
    g[c - 1] = (lambda2 + (c - 1) * mu2 * beta1) * alpha2 / (lambda2 * alpha1 + c * mu2)
    x = (lambda2 + c * mu2 * beta1) / (lambda2 * alpha1 + c * mu2)
    Gn = np.cumprod(g)
    Z = 1.0 + np.sum(Gn[:c - 1]) + Gn[c - 1] / (1.0 - x)
    S1 = sum((n + 1) * Gn[n] for n in range(c - 1))
    S2 = Gn[c - 1] * (c / (1.0 - x) + x / (1.0 - x) ** 2)
    return (S1 + S2) / Z


def me_oqn(
    M: int,
    R: int,
    lambda0: np.ndarray,
    Ca0: np.ndarray,
    mu: np.ndarray,
    Cs: np.ndarray,
    P: np.ndarray,
    c: Optional[np.ndarray] = None,
    insens: Optional[np.ndarray] = None,
    options: Optional[Dict[str, Any]] = None
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    """
    Maximum Entropy algorithm for Open Queueing Networks.

    Implements the ME algorithm from Kouvatsos (1994) for analyzing
    open queueing networks with general arrival and service processes.

    Args:
        M: Number of queues (stations).
        R: Number of job classes.
        lambda0: External arrival rates [M x R matrix], lambda0[i,r] = lambda_oi,r.
        Ca0: External arrival squared coefficient of variation [M x R matrix].
        mu: Service rates [M x R matrix], mu[i,r] = mu_i,r.
        Cs: Service squared coefficient of variation [M x R matrix].
        P: Routing probability matrix [M x M x R], P[j,i,r] = p_ji,r
           (probability class r goes from queue j to queue i).
        c: Optional servers per queue [M vector]; np.inf marks an
           infinite-server (IS) queue (default: all ones).
        insens: Optional boolean vector [M]; True marks a station with an
           insensitive scheduling discipline (PS, LCFS-PR), solved with the
           product-form mean queue length L_r = rho_r/(1-rho) instead of
           the FCFS GE formula (default: all False).
        options: Optional dict with fields:
            - tol: convergence tolerance (default: 1e-6)
            - maxiter: maximum iterations (default: 1000)
            - verbose: print iteration info (default: False)

    Returns:
        Tuple of (L, W, Ca, Cd, lambda_arr, rho, iters):
            L: Mean queue lengths [M x R matrix].
            W: Mean response times [M x R matrix], W = L / lambda.
            Ca: Arrival scv at each queue [M x R matrix].
            Cd: Departure scv at each queue [M x R matrix].
            lambda_arr: Total arrival rates [M x R matrix], inclusive of
                self-loop revisits (visit-based throughput).
            rho: Utilizations [M x R matrix] (per-server for finite c,
                mean number of busy servers for IS queues).
            iters: Number of iterations until convergence.

    Reference:
        Kouvatsos (1994), equations (3.3), (3.6), (3.7), (3.9) and the
        multiclass GE/GE/1/FCFS mql of Section 3.1.1.
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
    if c is None:
        c = np.ones(M)
    else:
        c = np.asarray(c, dtype=float).ravel()
    if insens is None:
        insens = np.zeros(M, dtype=bool)
    else:
        insens = np.asarray(insens, dtype=bool).ravel()

    # see _kb/03-api-layer.md for rationale
    P_eff = P.copy()
    mu_eff = mu.copy()
    Cs_eff = Cs.copy()
    for i in range(M):
        for r in range(R):
            pii = P[i, i, r]
            if pii > 0:
                mu_eff[i, r] = mu[i, r] * (1.0 - pii)
                Cs_eff[i, r] = pii + (1.0 - pii) * Cs[i, r]
                P_eff[i, :, r] = P[i, :, r] / (1.0 - pii)
                P_eff[i, i, r] = 0.0

    # see _kb/03-api-layer.md for rationale
    lambda_arr = np.zeros((M, R))
    lambda_eff = np.zeros((M, R))
    for r in range(R):
        Pr = P[:, :, r]
        A = np.eye(M) - Pr.T
        lambda_arr[:, r] = np.linalg.solve(A, lambda0[:, r])
        lambda_eff[:, r] = lambda_arr[:, r] * (1.0 - np.diag(Pr))

    # see _kb/03-api-layer.md for rationale
    rho = np.zeros((M, R))
    for i in range(M):
        for r in range(R):
            if mu[i, r] > 0:
                if np.isinf(c[i]):
                    rho[i, r] = lambda_eff[i, r] / mu_eff[i, r]
                else:
                    rho[i, r] = lambda_arr[i, r] / (c[i] * mu[i, r])

    # Stability check (finite-server queues only)
    unstable = np.zeros(M, dtype=bool)
    for i in range(M):
        if not np.isinf(c[i]) and np.sum(rho[i, :]) >= 1:
            unstable[i] = True
    if np.any(unstable):
        warnings.warn("Network is unstable (utilization >= 1 at some queues)",
                      RuntimeWarning)

    # Step 2: Initialize arrival scvs
    Ca = np.ones((M, R))
    Cd = np.ones((M, R))
    L = np.zeros((M, R))

    # Steps 4-5: fixed-point iteration on the arrival scvs
    delta = np.inf
    iters = 0
    for iteration in range(maxiter):
        iters = iteration + 1
        Ca_old = Ca.copy()

        # Step 4: GE-type mean queue length formulae
        for i in range(M):
            rho_i = np.sum(rho[i, :])
            if np.isinf(c[i]):
                # GE/GE/inf queue: L = lambda/mu
                for r in range(R):
                    if lambda_eff[i, r] > 0 and mu_eff[i, r] > 0:
                        L[i, r] = lambda_eff[i, r] / mu_eff[i, r]
            elif unstable[i]:
                for r in range(R):
                    if lambda_eff[i, r] > 0:
                        L[i, r] = np.inf
            elif c[i] == 1:
                if insens[i]:
                    # Insensitive disciplines (PS, LCFS-PR): product-form mql,
                    # exact irrespective of the service distribution
                    for r in range(R):
                        if lambda_eff[i, r] > 0 and mu_eff[i, r] > 0:
                            L[i, r] = rho[i, r] / (1.0 - rho_i)
                else:
                    # see _kb/03-api-layer.md for rationale
                    resid = 0.0
                    for u in range(R):
                        if lambda_eff[i, u] > 0 and mu_eff[i, u] > 0:
                            resid += lambda_eff[i, u] * (Cs_eff[i, u] + Ca[i, u]) / mu_eff[i, u] ** 2
                    for r in range(R):
                        if lambda_eff[i, r] > 0 and mu_eff[i, r] > 0:
                            L[i, r] = (rho[i, r] * (Ca[i, r] + 1.0) / 2.0
                                       + lambda_eff[i, r] * resid / (2.0 * (1.0 - rho_i)))
            else:
                # see _kb/03-api-layer.md for rationale
                lam_a = 0.0
                for u in range(R):
                    if lambda_eff[i, u] > 0 and mu_eff[i, u] > 0:
                        lam_a += lambda_eff[i, u]
                if lam_a > 0:
                    inv_a = 0.0
                    ES = 0.0
                    ES2 = 0.0
                    for u in range(R):
                        if lambda_eff[i, u] > 0 and mu_eff[i, u] > 0:
                            wu = lambda_eff[i, u] / lam_a
                            inv_a += wu / (Ca[i, u] + 1.0)
                            ES += wu / mu_eff[i, u]
                            ES2 += wu * (Cs_eff[i, u] + 1.0) / mu_eff[i, u] ** 2
                    Ca_a = -1.0 + 1.0 / inv_a
                    Cs_a = ES2 / ES ** 2 - 1.0
                    L_a = _ge_gec_mql(lam_a, Ca_a, 1.0 / ES, Cs_a, c[i])
                    Lq_a = L_a - lam_a * ES  # mean waiting-line length
                    for r in range(R):
                        if lambda_eff[i, r] > 0 and mu_eff[i, r] > 0:
                            L[i, r] = c[i] * rho[i, r] + (lambda_eff[i, r] / lam_a) * Lq_a

        # Step 5a: departure scvs
        for j in range(M):
            rho_j = np.sum(rho[j, :])
            for r in range(R):
                if lambda_eff[j, r] > 0:
                    if np.isinf(c[j]):
                        # GE/GE/inf queue: interdeparture scv = interarrival scv
                        Cd[j, r] = Ca[j, r]
                    elif unstable[j]:
                        # Saturated server: departures follow the service process
                        Cd[j, r] = Cs_eff[j, r]
                    elif c[j] == 1:
                        # Eq. (3.6) on the class-r virtual queue, with the
                        # marginal utilization rhohat_r of eq. (3.3)
                        rhohat = rho[j, r] * L[j, r] / (L[j, r] + rho_j - rho[j, r])
                        Cd[j, r] = 2.0 * L[j, r] * (1.0 - rhohat) + Ca[j, r] * (1.0 - 2.0 * rhohat)
                    else:
                        # GE/GE/c interdeparture scv (Section 4.2)
                        Cd[j, r] = (rho_j * (1.0 - rho_j) + (1.0 - rho_j) * Ca[j, r]
                                    + rho_j ** 2 * Cs_eff[j, r])

        # see _kb/03-api-layer.md for rationale
        for i in range(M):
            for r in range(R):
                if lambda_eff[i, r] > 0:
                    sum_inv = 0.0
                    for j in range(M):
                        pji = P_eff[j, i, r]
                        if pji > 0 and lambda_eff[j, r] > 0:
                            Cdji = 1.0 + pji * (Cd[j, r] - 1.0)
                            sum_inv += (lambda_eff[j, r] * pji / lambda_eff[i, r]) / (Cdji + 1.0)
                    if lambda0[i, r] > 0:
                        sum_inv += (lambda0[i, r] / lambda_eff[i, r]) / (Ca0[i, r] + 1.0)
                    if sum_inv > 0:
                        Ca[i, r] = -1.0 + 1.0 / sum_inv

        # Check convergence
        delta = np.max(np.abs(Ca - Ca_old))
        if verbose:
            print(f"Iteration {iters}: max delta = {delta:e}")
        if delta < tol:
            if verbose:
                print(f"Converged after {iters} iterations")
            break

    if iters == maxiter and delta >= tol:
        warnings.warn(
            f"Did not converge within {maxiter} iterations (delta={delta:e})",
            RuntimeWarning
        )

    # Step 6: response times by Little's law on the reported arrival rates
    W = np.zeros((M, R))
    for i in range(M):
        for r in range(R):
            if lambda_arr[i, r] > 0:
                W[i, r] = L[i, r] / lambda_arr[i, r]

    return L, W, Ca, Cd, lambda_arr, rho, iters


__all__ = ['me_oqn']
