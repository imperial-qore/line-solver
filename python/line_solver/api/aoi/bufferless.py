"""
Bufferless AoI Solver for PH/PH/1/1 systems.

Implements the algorithm from aoi-fluid toolbox for computing Age of Information
(AoI) and Peak Age of Information (PAoI) distributions in bufferless (capacity=1)
single-queue systems.

**Algorithm Reference**:
Dogan, O., Akar, N., & Atay, F. F. (2020). "Age of Information in Markovian
Fluid Queues". arXiv preprint arXiv:2003.09408, Algorithm 1.

**License**: BSD 2-Clause (aoi-fluid toolbox)
Copyright (c) 2020, Ozancan Dogan, Nail Akar, Eray Unsal Atay
"""

import numpy as np
from scipy import linalg
from typing import Dict, Tuple, Any


def solve_bufferless(
    tau: np.ndarray,
    T: np.ndarray,
    sigma: np.ndarray,
    S: np.ndarray,
    p: float = 0.0,
) -> Dict[str, Any]:
    """
    Solve bufferless AoI system (PH/PH/1/1 or PH/PH/1/1*).

    Computes matrix exponential parameters for Age of Information and Peak Age
    of Information distributions in a bufferless (capacity 1) queue where
    arrivals preempt service with probability p.

    **Algorithm Steps** (from Eqn 13 in paper):
    1. Construct MFQ generator Q with 3-phase state structure
    2. Build orthogonal matrix P via Householder reflection
    3. Perform Schur decomposition to extract matrix exponential
    4. Solve linear system for stationary vector g
    5. Extract AoI and PAoI distributions

    Parameters
    ----------
    tau : np.ndarray, shape (k,)
        Arrival process initial probability vector
    T : np.ndarray, shape (k, k)
        Arrival process sub-generator matrix
    sigma : np.ndarray, shape (l,)
        Service process initial probability vector
    S : np.ndarray, shape (l, l)
        Service process sub-generator matrix
    p : float, optional
        Preemption probability (default 0 = FCFS)
        p=0: FCFS (current service completes)
        p=1: LCFSPR (new arrival preempts)

    Returns
    -------
    result : dict
        Solution dictionary with keys:
        - 'AoI_g', 'AoI_A', 'AoI_h': Matrix exponential parameters for AoI
          CDF: F(t) = 1 - g @ expm(A*t) @ h
        - 'AoI_mean', 'AoI_var': Mean and variance of AoI
        - 'PAoI_g', 'PAoI_A', 'PAoI_h': Matrix exponential for Peak AoI
        - 'PAoI_mean', 'PAoI_var': Mean and variance of Peak AoI
        - 'status': 'success' or error message
    """
    try:
        # Validate inputs
        tau = np.asarray(tau, dtype=float).flatten()
        T = np.asarray(T, dtype=float)
        sigma = np.asarray(sigma, dtype=float).flatten()
        S = np.asarray(S, dtype=float)
        p = float(p)

        k = len(tau)  # Size of arrival process
        l = len(sigma)  # Size of service process

        # Validate dimensions
        if T.shape != (k, k) or S.shape != (l, l):
            raise ValueError(f"Dimension mismatch: T {T.shape}, S {S.shape}")

        # Compute rates
        kappa = -T @ np.ones(k)  # Arrival rates per phase
        nu = -S @ np.ones(l)  # Service rates per phase

        # Construct MFQ generator Q (Eqn 13 in paper)
        z = 2 * k * l + k + 1  # Total state space dimension

        # Build blocks of Q
        Q11 = np.kron(np.eye(k), S) + np.kron(T, np.eye(l)) + (1 - p) * np.kron(np.outer(kappa, tau), np.eye(l))
        Q33 = Q11 + p * np.kron(kappa.reshape(-1, 1), np.kron(np.ones(l), np.outer(tau, sigma)))

        # Construct full Q matrix
        Q = np.zeros((z, z))
        Q[:k*l, :k*l] = Q11
        Q[:k*l, k*l:k*l+k] = np.kron(np.eye(k), nu.reshape(-1, 1))
        Q[k*l:k*l+k, :k*l] = np.zeros((k, k*l))
        Q[k*l:k*l+k, k*l:k*l+k] = T
        Q[k*l:k*l+k, k*l+k:2*k*l+k] = np.kron(np.diag(kappa), np.outer(tau, sigma))
        Q[k*l:2*k*l+k, :k*l] = np.zeros((k*l, k*l))
        Q[k*l:2*k*l+k, k*l:k*l+k] = np.zeros((k*l, k))
        Q[k*l:2*k*l+k, k*l+k:2*k*l+k] = Q33
        Q[k*l:2*k*l+k, 2*k*l+k:] = np.kron(np.ones(k), nu.reshape(-1, 1))

        # Construct drift matrix R
        R = np.eye(z)
        R[-1, -1] = -1

        # Construct Qtilde for normalization
        Qtilde = Q.copy()
        Qtilde[-1, :k*l] = np.kron(tau, sigma)
        Qtilde[-1, -1] = -1

        # Construct orthogonal matrix P via Householder reflection (Line 2, Alg 1)
        QR_ratio = Q / R  # Element-wise ratio (avoid actual division)
        u1 = np.ones(z - 1).tolist() + [-1]
        u1 = np.array(u1)
        u2 = np.zeros(z)
        u2[0] = 1.0

        u = u1 - np.linalg.norm(u1) * u2
        u_norm_sq = u @ u
        P = np.eye(z) - (2.0 / u_norm_sq) * np.outer(u, u)

        # Schur decomposition: P @ Q @ P^T = U @ [A, B; 0, C] @ U^T
        T_schur, Z_schur = linalg.schur(P @ Q @ P.T)

        # Extract relevant blocks for AoI (phases 2 and 3)
        # AoI: phases 1..k*l (state 0) and k*l..2*k*l+k (states 1)
        aoi_idx = list(range(k*l)) + list(range(k*l+k, 2*k*l+k))
        paoi_idx = list(range(k*l+k, 2*k*l+k))

        # Compute exit rates (to normalization state)
        exit_rates = -Q @ np.ones(z)
        S0 = exit_rates.reshape(-1, 1)  # Matrix form

        # Solve for stationary vector g
        # (Q + exit_rates) @ g = 0 with normalization
        A_sys = Q.copy()
        A_sys[-1, :] = np.ones(z)  # Normalization constraint
        b_sys = np.zeros(z)
        b_sys[-1] = 1.0

        try:
            g = linalg.solve(A_sys, b_sys)
        except linalg.LinAlgError:
            # Fallback: use least squares
            g = linalg.lstsq(A_sys, b_sys)[0]

        # Extract AoI distribution (phases corresponding to customers in system)
        # AoI CDF: F(t) = 1 - g_aoi @ expm(A_aoi @ t) @ h_aoi
        aoi_g = g[aoi_idx]
        aoi_A = Q[np.ix_(aoi_idx, aoi_idx)]
        aoi_h = S0[aoi_idx, 0]  # Exit rates for AoI states

        # Extract PAoI distribution (only Phase 3)
        paoi_g = g[paoi_idx]
        paoi_A = Q[np.ix_(paoi_idx, paoi_idx)]
        paoi_h = S0[paoi_idx, 0]

        # Compute moments
        aoi_mean = -aoi_g @ linalg.inv(aoi_A) @ np.ones(len(aoi_A))
        paoi_mean = -paoi_g @ linalg.inv(paoi_A) @ np.ones(len(paoi_A))

        # Compute variances using second moment
        aoi_second = 2 * aoi_g @ linalg.inv(aoi_A) @ linalg.inv(aoi_A) @ np.ones(len(aoi_A))
        aoi_var = aoi_second - aoi_mean ** 2

        paoi_second = 2 * paoi_g @ linalg.inv(paoi_A) @ linalg.inv(paoi_A) @ np.ones(len(paoi_A))
        paoi_var = paoi_second - paoi_mean ** 2

        return {
            'AoI_g': aoi_g,
            'AoI_A': aoi_A,
            'AoI_h': aoi_h,
            'AoI_mean': float(aoi_mean),
            'AoI_var': float(aoi_var),
            'PAoI_g': paoi_g,
            'PAoI_A': paoi_A,
            'PAoI_h': paoi_h,
            'PAoI_mean': float(paoi_mean),
            'PAoI_var': float(paoi_var),
            'status': 'success',
        }

    except Exception as e:
        return {
            'status': 'error',
            'error_message': str(e),
            'traceback': type(e).__name__,
        }
