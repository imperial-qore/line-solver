"""
QRF No-Blocking MEM (Maximum Entropy Method) Approximation.

Solves a constrained NLP to approximate performance metrics of single-class
closed queueing networks with PH service, using maximum entropy
as the objective.

Port of MATLAB qrf_noblo_mem.m.
"""

import numpy as np

from .qrf_noblo_common import (
    sub_qrfvar,
    sub_qrfcon_noblo,
    mem_objective,
    compute_num_vars,
    extract_results,
    build_q_from_mu_v_rt,
    extract_mu_v_from_maps,
    affine_constraint_matrices,
    reduce_equalities,
    feasible_start,
    qrf_index_map,
    solve_qrf_nlp,
    mem_gradient,
)


def qrf_noblo_mem(MAPs, N, rt):
    """QRF no-blocking MEM approximation.

    Args:
        MAPs: List of [D0, D1] pairs per queue, where D0 is the background
              transition matrix and D1 is the completion rate matrix.
        N: Total population
        rt: Routing matrix [M, M]

    Returns:
        UN: Utilization per queue [M]
        QN: Queue length per queue [M]
    """
    M = len(MAPs)

    # Extract phases per queue
    K = np.array([MAPs[i][0].shape[0] for i in range(M)], dtype=int)

    # Extract mu, v from MAPs
    mu, v = extract_mu_v_from_maps(MAPs, M, K)

    # No-blocking defaults
    MR = 1
    BB = np.zeros((1, M))
    F = np.full(M, N, dtype=int)

    # Build transition rates
    q = build_q_from_mu_v_rt(M, K, mu, v, rt)

    num_vars = compute_num_vars(M, N, K, MR)

    # Bounds: all variables in [0, 1]
    bounds = [(0.0, 1.0)] * num_vars

    # see _kb/03-api-layer.md for rationale
    Aeq, beq = affine_constraint_matrices(
        lambda z: sub_qrfcon_noblo(z, q, M, MR, BB, F, N, K)[1], num_vars)
    Aub, bub = affine_constraint_matrices(
        lambda z: sub_qrfcon_noblo(z, q, M, MR, BB, F, N, K)[0], num_vars)
    Aeq, beq, _ = reduce_equalities(Aeq, beq)

    # see _kb/03-api-layer.md for rationale
    x0 = feasible_start(Aeq, beq, Aub, bub, num_vars)

    idx = qrf_index_map(M, N, K, MR)
    xopt = solve_qrf_nlp(
        lambda x: mem_objective(x, M, N, K, F, MR),
        lambda x: mem_gradient(x, M, N, K, F, MR, idx),
        x0, Aeq, beq, Aub, bub, 'qrf_noblo_mem')

    p2opt, _ = sub_qrfvar(xopt, M, N, K, MR)
    UN, QN = extract_results(p2opt, M, K, F, MR)

    return UN, QN
