"""
Single-buffer AoI solver for M/PH/1/2 systems.

This is a direct NumPy/SciPy port of LINE's MATLAB `solveSingleBuffer`
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
    sol, _, _, _ = linalg.lstsq(eqn_matrix.T, rhs)
    return sol


def _real_scalar(value: np.ndarray) -> float:
    return float(np.real_if_close(np.asarray(value).item()))


def solve_singlebuffer(
    lambda_rate: float,
    sigma: np.ndarray,
    S: np.ndarray,
    r: float = 0.0,
) -> Dict[str, Any]:
    """
    Solve single-buffer AoI system (M/PH/1/2 or M/PH/1/2*).

    Returns matrix-exponential parameters for AoI and Peak AoI plus their
    first two moments.
    """
    try:
        lambda_rate = float(lambda_rate)
        sigma = _as_row(sigma)
        S = _as_matrix(S)
        r = float(r)

        l = S.shape[1]
        if S.shape != (l, l):
            raise ValueError(f"Dimension mismatch: S {S.shape}")

        ones_l = np.ones((l, 1))
        nu = -S @ ones_l

        # Part 1: waiting-time MFQ (MATLAB Eqn. 16 path).
        z1 = l + 2
        a1 = 2
        b1 = l

        q1 = np.block([
            [S, nu, np.zeros((l, 1))],
            [np.zeros((1, l)), np.array([[-lambda_rate, lambda_rate]])],
            [np.zeros((1, l + 2))],
        ])

        r1 = np.eye(z1)
        r1[-2, -2] = -1.0
        r1[-1, -1] = -1.0

        qtilde1 = np.block([
            [np.zeros((l, l + 2))],
            [lambda_rate * sigma, np.array([[-lambda_rate, 0.0]])],
            [sigma, np.array([[0.0, -1.0]])],
        ])

        e1 = np.ones((z1, 1))
        pik = linalg.solve((q1 + e1 @ e1.T).T, e1).T
        qr1 = q1 @ np.linalg.inv(r1)
        x_r = r1 @ e1
        x_l = pik
        a1mat = qr1 + (x_r @ x_l) / float((x_l @ x_r).item())

        t_schur, z_schur, _ = linalg.schur(
            a1mat,
            output='real',
            sort=lambda eig: np.real(eig) > 0,
        )
        p1 = z_schur

        ae1 = p1.T @ qr1 @ p1
        amat1 = ae1[a1:, a1:]
        hmat1 = p1[:, a1:].T
        qtildestar1 = qtilde1[b1:, :]

        amat1_inv = np.linalg.inv(amat1)
        eqn_matrix1 = np.block([
            [hmat1 @ r1, -amat1_inv @ hmat1 @ np.ones((z1, 1))],
            [-qtildestar1, np.ones((a1, 1))],
        ])
        rhs1 = np.zeros((z1 + 1, 1))
        rhs1[-1, 0] = 1.0
        solution1 = _solve_mldivide_transpose(eqn_matrix1, rhs1)

        wait_g = solution1[:b1, :].T
        wait_d = solution1[b1:, :].T
        c0 = wait_d[0, 0]

        wait_a = amat1 - r * lambda_rate * np.eye(l)
        wait_h = hmat1 @ np.vstack([np.zeros((l, 1)), [[1.0]], [[r]]])

        wait_a_inv = np.linalg.inv(wait_a)
        n1 = 1.0 / _real_scalar(-(wait_g @ wait_a_inv @ wait_h) + c0)
        wait_g = n1 * wait_g

        mdiag = (-wait_a_inv @ wait_h).ravel()
        mmat = np.diag(mdiag)
        bmat = np.linalg.solve(mmat, wait_a @ mmat)
        beta = wait_g @ mmat
        beta0 = 1.0 - float(np.sum(beta))
        psi = -bmat @ ones_l

        # Part 2: AoI MFQ (MATLAB Eqn. 19 path).
        z2 = 4 * l + 2
        a2 = 1
        b2 = z2 - 1

        q2 = np.block([
            [bmat, np.kron(psi, sigma), np.zeros((l, 2 * l + 2))],
            [
                np.zeros((l, l)),
                S - lambda_rate * np.eye(l),
                lambda_rate * np.eye(l),
                nu,
                np.zeros((l, l + 1)),
            ],
            [
                np.zeros((l, 2 * l)),
                S,
                np.zeros((l, 1)),
                np.kron(nu, sigma),
                np.zeros((l, 1)),
            ],
            [np.zeros((1, 3 * l)), np.array([[-lambda_rate]]), lambda_rate * sigma, np.zeros((1, 1))],
            [np.zeros((l, 3 * l + 1)), S, nu],
            [np.zeros((1, z2))],
        ])

        r2 = np.eye(z2)
        r2[-1, -1] = -1.0

        qtilde2 = q2.copy()
        qtilde2[-1, :l] = beta
        qtilde2[-1, l:2 * l] = beta0 * sigma
        qtilde2[-1, -1] = -1.0

        u1 = np.vstack([np.ones((z2 - 1, 1)), [[-1.0]]])
        u2 = np.vstack([[[1.0]], np.zeros((z2 - 1, 1))])
        u = u1 - np.linalg.norm(u1, ord=2) * u2
        p2 = np.eye(z2) - (2.0 * (u @ u.T)) / float((u.T @ u).item())

        qr2 = q2 @ np.linalg.inv(r2)
        ae2 = p2.T @ qr2 @ p2
        amat2 = ae2[a2:, a2:]
        hmat2 = p2[:, a2:].T
        qtildestar2 = qtilde2[b2:, :]

        amat2_inv = np.linalg.inv(amat2)
        eqn_matrix2 = np.block([
            [hmat2 @ r2, -amat2_inv @ hmat2 @ np.ones((z2, 1))],
            [-qtildestar2, np.ones((a2, 1))],
        ])
        rhs2 = np.zeros((z2 + 1, 1))
        rhs2[-1, 0] = 1.0
        solution2 = _solve_mldivide_transpose(eqn_matrix2, rhs2)

        g = solution2[:b2, :].T

        selected_aoi = np.vstack([
            np.zeros((3 * l, 1)),
            np.ones((l + 1, 1)),
            [[0.0]],
        ])
        aoi_h = hmat2 @ selected_aoi
        norm_constant = -(g @ amat2_inv @ aoi_h)
        aoi_g = g / norm_constant
        aoi_a = amat2

        amat2_inv2 = amat2_inv @ amat2_inv
        amat2_inv3 = amat2_inv2 @ amat2_inv
        aoi_mean = _real_scalar(aoi_g @ amat2_inv2 @ aoi_h)
        aoi_var = max(0.0, _real_scalar(-2.0 * aoi_g @ amat2_inv3 @ aoi_h - aoi_mean ** 2))

        selected_paoi = np.vstack([
            np.zeros((3 * l + 1, 1)),
            nu,
            [[0.0]],
        ])
        paoi_h = hmat2 @ selected_paoi
        norm_constant = -(g @ amat2_inv @ paoi_h)
        paoi_g = g / norm_constant
        paoi_a = amat2
        paoi_mean = _real_scalar(paoi_g @ amat2_inv2 @ paoi_h)
        paoi_var = max(0.0, _real_scalar(-2.0 * paoi_g @ amat2_inv3 @ paoi_h - paoi_mean ** 2))

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
            'systemType': 'singlebuffer',
            'preemption': r,
            'status': 'success',
        }
    except Exception as err:
        return {
            'status': 'error',
            'error_message': str(err),
            'traceback': type(err).__name__,
        }
