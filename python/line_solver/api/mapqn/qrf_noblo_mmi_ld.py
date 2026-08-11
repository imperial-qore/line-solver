"""
QRF No-Blocking MMI with Load-Dependent Rates Approximation.

Solves a constrained NLP with load-dependent transition rates.
Uses mutual information minimization as the objective.

Port of MATLAB qrf_noblo_mmi_ld.m.
"""

import numpy as np
from scipy.optimize import minimize

from .qrf_noblo_common import (
    sub_qrfvar,
    sub_qrfcon_noblo,
    mmi_objective,
    compute_num_vars,
    extract_results,
    extract_mu_v_from_maps,
    build_q_ld,
    affine_constraint_matrices,
    reduce_equalities,
    feasible_start,
)


def qrf_noblo_mmi_ld(MAPs, N, rt, alpha=None):
    """QRF no-blocking MMI approximation with load-dependent rates.

    Args:
        MAPs: List of [D0, D1] pairs per queue
        N: Total population
        rt: Routing matrix [M, M]
        alpha: Load-dependent scaling factors [M, N]. If None, defaults to ones.

    Returns:
        UN: Utilization per queue [M]
        QN: Queue length per queue [M]
    """
    M = len(MAPs)
    K = np.array([MAPs[i][0].shape[0] for i in range(M)], dtype=int)

    if alpha is None:
        alpha = np.ones((M, N))

    mu, v = extract_mu_v_from_maps(MAPs, M, K)

    MR = 1
    BB = np.zeros((1, M))
    F = np.full(M, N, dtype=int)

    # see _kb/03-api-layer.md for rationale
    q = build_q_ld(M, K, mu, v, rt, N, alpha)

    num_vars = compute_num_vars(M, N, K, MR)

    # Bounds
    bounds = [(0.0, 1.0)] * num_vars

    # see _kb/03-api-layer.md for rationale
    Aeq, beq = affine_constraint_matrices(
        lambda z: sub_qrfcon_noblo(z, q, M, MR, BB, F, N, K)[1], num_vars)
    Aub, bub = affine_constraint_matrices(
        lambda z: sub_qrfcon_noblo(z, q, M, MR, BB, F, N, K)[0], num_vars)
    Aeq, beq, _ = reduce_equalities(Aeq, beq)

    # see _kb/03-api-layer.md for rationale
    x0 = feasible_start(Aeq, beq, Aub, bub, num_vars)

    constraints = [{'type': 'eq', 'fun': lambda x: Aeq @ x - beq}]
    if Aub.shape[0] > 0:
        constraints.append({'type': 'ineq', 'fun': lambda x: -(Aub @ x - bub)})

    result = minimize(
        lambda x: mmi_objective(x, M, N, K, F, MR),
        x0,
        method='SLSQP',
        bounds=bounds,
        constraints=constraints,
        options={'maxiter': 100, 'disp': False, 'ftol': 1e-8},
    )

    p2opt, _ = sub_qrfvar(result.x, M, N, K, MR)
    UN, QN = extract_results(p2opt, M, K, F, MR)

    return UN, QN
