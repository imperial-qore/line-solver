"""QRF No-Blocking BETHE approximation.

Minimises the tree-reweighted (Bethe) free entropy over the same no-blocking
polytope qrf_noblo_mmi uses. The polytope, the phase-1 feasible start and the
SLSQP call are shared infrastructure; the objective is the only difference.
See qrf_noblo_common.bethe_objective for what it is and why the uniform edge
weight is lambda = 1/M.

ONE SOLVE, NO RESTARTS. The objective is convex on this polytope, so there is
no second local minimum for a restart to find: the single solve from the
phase-1 point returns the global optimum, and the answer is a property of the
model rather than of the start point.

Port of MATLAB qrf_noblo_bethe.m.
"""

import numpy as np

from .qrf_noblo_common import (
    sub_qrfvar,
    sub_qrfcon_noblo,
    bethe_objective,
    compute_num_vars,
    extract_results,
    build_q_from_mu_v_rt,
    affine_constraint_matrices,
    reduce_equalities,
    feasible_start,
    qrf_index_map,
    solve_qrf_nlp,
    bethe_gradient,
)


def qrf_noblo_bethe(M, MR, K, N, mu, v, rt):
    """QRF no-blocking approximation under the tree-reweighted free entropy.

    Args:
        M: Number of queues
        MR: Ignored -- no blocking means one configuration, by definition
        K: Phases per queue, array of int [M]
        N: Total population
        mu: Completion rates [M, Kmax, Kmax]
        v: Background rates [M, Kmax, Kmax]
        rt: Routing matrix [M, M]

    Returns:
        UN: Utilization per queue [M]
        QN: Queue length per queue [M]
    """
    K = np.asarray(K, dtype=int)

    # No-blocking defaults
    MR = 1
    BB = np.zeros((1, M))
    F = np.full(M, N, dtype=int)

    # Build transition rates
    q = build_q_from_mu_v_rt(M, K, mu, v, rt)

    num_vars = compute_num_vars(M, N, K, MR)

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
        lambda x: bethe_objective(x, M, N, K, F, MR),
        lambda x: bethe_gradient(x, M, N, K, F, MR, idx),
        x0, Aeq, beq, Aub, bub, 'qrf_noblo_bethe')

    p2opt, _ = sub_qrfvar(xopt, M, N, K, MR)
    UN, QN = extract_results(p2opt, M, K, F, MR)

    return UN, QN
