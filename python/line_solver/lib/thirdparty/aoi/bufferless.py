"""
Bufferless AoI solver for PH/PH/1/1 systems.

This is a direct NumPy/SciPy port of LINE's MATLAB `solveBufferless`
implementation from the bundled aoi-fluid toolbox.
"""

from __future__ import annotations

from typing import Any, Dict

import numpy as np
from scipy import linalg


def _as_row(vec: np.ndarray) -> np.ndarray:
    return np.asarray(vec, dtype=float).reshape(1, -1)


def _as_matrix(mat: np.ndarray) -> np.ndarray:
    return np.asarray(mat, dtype=float)


def _solve_mldivide_transpose(eqn_matrix: np.ndarray, rhs: np.ndarray) -> np.ndarray:
    """Mirror MATLAB `mldivide(EqnMatrix', RightHandSide')`."""
    sol, _, _, _ = linalg.lstsq(eqn_matrix.T, rhs)
    return sol


def _real_scalar(value: np.ndarray) -> float:
    return float(np.real_if_close(np.asarray(value).item()))


def solve_bufferless(
    tau: np.ndarray,
    T: np.ndarray,
    sigma: np.ndarray,
    S: np.ndarray,
    p: float = 0.0,
) -> Dict[str, Any]:
    """
    Solve bufferless AoI system (PH/PH/1/1 or PH/PH/1/1*).

    Returns matrix-exponential parameters for AoI and Peak AoI plus their
    first two moments.
    """
    try:
        tau = _as_row(tau)
        T = _as_matrix(T)
        sigma = _as_row(sigma)
        S = _as_matrix(S)
        p = float(p)

        k = T.shape[1]
        l = S.shape[1]
        if T.shape != (k, k) or S.shape != (l, l):
            raise ValueError(f"Dimension mismatch: T {T.shape}, S {S.shape}")

        ones_k = np.ones((k, 1))
        ones_l = np.ones((l, 1))
        kappa = -T @ ones_k
        nu = -S @ ones_l

        z = 2 * k * l + k + 1
        a = 1
        b = z - 1

        q11 = (
            np.kron(np.eye(k), S)
            + np.kron(T, np.eye(l))
            + (1.0 - p) * np.kron(np.kron(kappa, tau), np.eye(l))
        )
        q33 = q11 + p * np.kron(kappa, np.kron(np.ones((l, 1)), np.kron(tau, sigma)))

        q = np.block([
            [
                q11,
                np.kron(np.eye(k), nu),
                np.zeros((k * l, k * l)),
                p * np.kron(kappa, np.ones((l, 1))),
            ],
            [
                np.zeros((k, k * l)),
                T,
                np.kron(kappa, np.kron(tau, sigma)),
                np.zeros((k, 1)),
            ],
            [
                np.zeros((k * l, k * l)),
                np.zeros((k * l, k)),
                q33,
                np.kron(np.ones((k, 1)), nu),
            ],
            [np.zeros((1, z))],
        ])

        r = np.eye(z)
        r[-1, -1] = -1.0

        qtilde = q.copy()
        qtilde[-1, :k * l] = np.kron(tau, sigma)
        qtilde[-1, -1] = -1.0

        qr = q @ np.linalg.inv(r)

        u1 = np.vstack([np.ones((z - 1, 1)), [[-1.0]]])
        u2 = np.vstack([[[1.0]], np.zeros((z - 1, 1))])
        u = u1 - np.linalg.norm(u1, ord=2) * u2
        pmat = np.eye(z) - (2.0 * (u @ u.T)) / float((u.T @ u).item())

        ae = pmat.T @ qr @ pmat
        amat = ae[a:, a:]
        hmat = pmat[:, a:].T
        qtildestar = qtilde[b:, :]

        amat_inv = np.linalg.inv(amat)
        eqn_matrix = np.block([
            [hmat @ r, -amat_inv @ hmat @ np.ones((z, 1))],
            [-qtildestar, np.ones((a, 1))],
        ])
        rhs = np.zeros((z + 1, 1))
        rhs[-1, 0] = 1.0
        solution = _solve_mldivide_transpose(eqn_matrix, rhs)

        g = solution[:b, :].T

        selected_aoi = np.vstack([
            np.zeros((k * l, 1)),
            np.ones((k * l + k, 1)),
            [[0.0]],
        ])
        aoi_h = hmat @ selected_aoi
        norm_constant = -(g @ amat_inv @ aoi_h)
        aoi_g = g / norm_constant
        aoi_a = amat

        amat_inv2 = amat_inv @ amat_inv
        amat_inv3 = amat_inv2 @ amat_inv
        aoi_mean = _real_scalar(aoi_g @ amat_inv2 @ aoi_h)
        aoi_var = max(0.0, _real_scalar(-2.0 * aoi_g @ amat_inv3 @ aoi_h - aoi_mean ** 2))

        selected_paoi = np.vstack([
            np.zeros((k * l + k, 1)),
            np.kron(np.ones((k, 1)), nu),
            [[0.0]],
        ])
        paoi_h = hmat @ selected_paoi
        norm_constant = -(g @ amat_inv @ paoi_h)
        paoi_g = g / norm_constant
        paoi_a = amat
        paoi_mean = _real_scalar(paoi_g @ amat_inv2 @ paoi_h)
        paoi_var = max(0.0, _real_scalar(-2.0 * paoi_g @ amat_inv3 @ paoi_h - paoi_mean ** 2))

        return {
            'AoI_g': np.real_if_close(aoi_g).ravel(),
            'AoI_A': np.real_if_close(aoi_a),
            'AoI_h': np.real_if_close(aoi_h).ravel(),
            'AoI_mean': aoi_mean,
            'AoI_var': aoi_var,
            'PAoI_g': np.real_if_close(paoi_g).ravel(),
            'PAoI_A': np.real_if_close(paoi_a),
            'PAoI_h': np.real_if_close(paoi_h).ravel(),
            'PAoI_mean': paoi_mean,
            'PAoI_var': paoi_var,
            'systemType': 'bufferless',
            'preemption': p,
            'status': 'success',
        }
    except Exception as err:
        return {
            'status': 'error',
            'error_message': str(err),
            'traceback': type(err).__name__,
        }
