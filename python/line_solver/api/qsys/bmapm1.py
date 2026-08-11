"""
BMAP/M/1 queue analysed by the matrix-analytic (M/G/1-type) method.

Port from: matlab/src/api/qsys/qsys_bmapm1.m

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

from typing import Any, Dict, List, Optional, Sequence

import numpy as np

__all__ = ['qsys_bmapm1']

# Matches MATLAB GlobalConstants.FineTol, which the input validation below uses
# as its shape and consistency tolerance.
_FINE_TOL = 1e-12


def qsys_bmapm1(D: Sequence[np.ndarray], mu: float,
                uniformization: Optional[float] = None,
                max_iter: int = 10000,
                tolerance: float = 1e-12,
                max_level: Optional[int] = None,
                tail_tolerance: float = 1e-10) -> Dict[str, Any]:
    """
    Analyze a single-server queue fed by a batch Markovian arrival process and
    with exponential service of rate mu.

    Args:
        D: list of BMAP matrices [D0, D1, ..., DK]. D0 carries the hidden
            transitions, Dk (k >= 1) the transitions that release a batch of k
            customers.
        mu: exponential service rate.
        uniformization: uniformization constant q used to randomize the
            generator into a discrete-time M/G/1-type chain. It must dominate
            every total outflow rate; by default it is chosen as
            max_i(-D0[i,i]) + mu.
        max_iter: maximum functional iterations for G.
        tolerance: convergence tolerance for G.
        max_level: level truncation used for the queue-length distribution
            (default: adaptive).
        tail_tolerance: relative truncation target for the level distribution.

    Beyond the usual performance measures the result exposes the intermediate
    matrix-analytic quantities themselves, so that the algorithm can be
    inspected and taught rather than only its output:

        theta        - stationary vector of the BMAP phase process, sum_k D_k
        lambda       - mean arrival rate, theta * sum_k k*D_k * e
        rho          - offered load lambda/mu
        q            - uniformization constant actually used
        A0, A1, Bk   - randomized blocks: A0 = (mu/q)I is a service completion
                       (level down by one), A1 = (1/q)(D0 - mu*I) + I keeps the
                       level, Bk[k] = (1/q)D_(k+1) raises the level by k+1
        B0           - boundary local block (1/q)D0 + I, used at level 0 where
                       no service can complete
        A            - A0 + A1 + sum_k Bk[k], the phase process of the chain
        alpha        - stationary vector of A
        G            - minimal non-negative solution of
                       G = A0 + A1*G + sum_k Bk[k]*G^(k+1)
        drift        - alpha*(sum_k k*Bk[k])*e - alpha*A0*e. The queue is
                       stable iff this is strictly negative
        decayRate    - geometric decay rate of the level probabilities,
                       measured as the limiting ratio pi_(n+1)/pi_n. Reported
                       rather than derived from a spectral convention so that
                       it is unambiguous
        levelProb    - level probabilities pi_n as rows (level 0 first)
        pi0          - probability the system is empty (equals 1-rho exactly)

    Example:
        # Example 6.4 of Bolch et al.
        D0 = np.array([[-2, 0.5], [1/3, -3]])
        D1 = np.array([[0.25, 0.5], [1/3, 1.0]])
        D2 = np.array([[0.25, 0.5], [1.0, 1/3]])
        result = qsys_bmapm1([D0, D1, D2], 11)

    See also qsys_mapm1, qsys_mapph1, qsys_bmapphnn_retrial.
    """
    from ..mc.ctmc import ctmc_solve
    from ..mc.dtmc import dtmc_solve
    from ..io.logging import line_warning

    # --- Validate inputs ---
    if not isinstance(D, (list, tuple)) or len(D) < 2:
        raise ValueError('The BMAP must be given as a list [D0,D1,...,DK] with '
                         'at least D0 and D1.')
    D = [np.atleast_2d(np.asarray(Dk, dtype=float)) for Dk in D]
    V = D[0].shape[0]
    for k, Dk in enumerate(D):
        if Dk.shape != (V, V) or not np.all(np.isfinite(Dk)):
            raise ValueError('BMAP matrix D[%d] must be a finite %dx%d matrix.' % (k, V, V))
        if k > 0 and np.any(Dk < -_FINE_TOL):
            raise ValueError('BMAP arrival matrix D[%d] must be non-negative.' % k)
    Dsum = sum(D)
    if np.any(np.abs(Dsum @ np.ones(V)) > np.sqrt(_FINE_TOL)):
        raise ValueError('BMAP matrices are inconsistent: sum_k D_k must have '
                         'zero row sums.')
    mu = float(mu)
    if not np.isfinite(mu) or mu <= 0:
        raise ValueError('The service rate mu must be a finite positive scalar.')

    K = len(D) - 1

    # --- Arrival characterization ---
    theta = np.asarray(ctmc_solve(Dsum), dtype=float).ravel()
    sum_k_Dk = np.zeros((V, V))
    for k in range(1, K + 1):
        sum_k_Dk += k * D[k]
    lam = float(theta @ sum_k_Dk @ np.ones(V))
    rho = lam / mu

    # --- Randomization (uniformization) into a discrete-time M/G/1-type chain ---
    outflow = float(np.max(-np.diag(D[0])) + mu)
    if uniformization is None:
        q = outflow
    else:
        q = float(uniformization)
    if q < outflow - _FINE_TOL:
        raise ValueError('The uniformization constant q = %g does not dominate the '
                         'total outflow rate %g; the randomized chain would have '
                         'negative entries.' % (q, outflow))

    I_V = np.eye(V)
    A0 = (mu / q) * I_V                        # level down by one: service completion
    A1 = (1.0 / q) * (D[0] - mu * I_V) + I_V   # level unchanged
    B0 = (1.0 / q) * D[0] + I_V                # level 0: no service can complete
    Bk = [(1.0 / q) * D[k] for k in range(1, K + 1)]  # level up by k

    A = A0 + A1
    for Bkk in Bk:
        A = A + Bkk
    alpha = np.asarray(dtmc_solve(A), dtype=float).ravel()

    # --- Matrix G: minimal non-negative solution of the M/G/1-type equation ---
    G = np.zeros((V, V))
    converged = False
    last_change = np.inf
    for _ in range(int(max_iter)):
        Gpow = G
        Gnew = A0 + A1 @ G
        for k in range(K):
            Gpow = Gpow @ G          # G^(k+2) after k increments from G^1
            Gnew = Gnew + Bk[k] @ Gpow
        last_change = float(np.max(np.abs(Gnew - G)))
        if last_change < tolerance:
            G = Gnew
            converged = True
            break
        G = Gnew
    if not converged:
        line_warning('qsys_bmapm1',
                     'The functional iteration for G did not converge to %g in %d '
                     'iterations (last change %g). The queue may be unstable.'
                     % (tolerance, int(max_iter), last_change))

    # --- Stability drift ---
    up_drift = np.zeros((V, V))
    for k in range(K):
        up_drift = up_drift + (k + 1) * Bk[k]
    drift = float(alpha @ up_drift @ np.ones(V) - alpha @ A0 @ np.ones(V))

    # see _kb/03-api-layer.md for rationale
    if max_level is not None and max_level > 0:
        level_max = int(round(max_level))
        level_prob = _solve_levels(D, mu, V, K, level_max)
        trunc_error = _level_tail_error(level_prob, level_max)
    else:
        level_max = max(50, int(np.ceil(20.0 / max(1.0 - min(rho, 0.999),
                                                   np.finfo(float).eps))))
        while True:
            level_prob = _solve_levels(D, mu, V, K, level_max)
            trunc_error = _level_tail_error(level_prob, level_max)
            if trunc_error <= tail_tolerance or (2 * level_max + 1) * V > 2e5:
                break
            level_max = 2 * level_max
        if trunc_error > tail_tolerance:
            line_warning('qsys_bmapm1',
                         'The level distribution did not reach the requested accuracy: '
                         'residual %.3e > TailTolerance %.3e at level %d.'
                         % (trunc_error, tail_tolerance, level_max))

    level_mass = level_prob.sum(axis=1)
    # Measured decay rate: the ratio settles geometrically, so read it where the
    # mass is still numerically meaningful rather than at the truncation boundary.
    usable_idx = np.nonzero(level_mass > 1e-12)[0]
    if usable_idx.size == 0 or (usable_idx[-1] + 1) < 3:
        decay_rate = float('nan')
    else:
        usable = int(usable_idx[-1]) + 1        # 1-based, as in MATLAB
        ref = max(2, usable // 2)               # 1-based reference level
        decay_rate = float(level_mass[ref] / level_mass[ref - 1])

    mean_queue_length = float(np.arange(level_prob.shape[0]) @ level_mass)

    return {
        'theta': theta,
        'lambda': lam,
        'rho': rho,
        'q': q,
        'A0': A0,
        'A1': A1,
        'B0': B0,
        'Bk': Bk,
        'A': A,
        'alpha': alpha,
        'G': G,
        'drift': drift,
        'decayRate': decay_rate,
        'levelProb': level_prob,
        'pi0': float(level_mass[0]),
        'meanQueueLength': mean_queue_length,
        'utilization': rho,
        'throughput': lam,
        'truncLevel': level_prob.shape[0] - 1,
        'truncError': trunc_error,
        'analyzer': 'LINE:qsys_bmapm1',
    }


def _level_tail_error(level_prob: np.ndarray, level_max: int) -> float:
    """Relative contribution the truncated tail would add to the mean level."""
    level_mass = level_prob.sum(axis=1)
    mean_level = float(np.arange(level_max + 1) @ level_mass)
    return float(level_max * level_mass[-1] / max(mean_level, np.finfo(float).tiny))


def _solve_levels(D: List[np.ndarray], mu: float, V: int, K: int,
                  level_max: int) -> np.ndarray:
    """
    Level-truncated CTMC generator of the BMAP/M/1 queue and its stationary
    distribution. Level n holds n customers in the system; the phase is the
    BMAP state. Service fires only above level 0.
    """
    from scipy.sparse import coo_matrix, diags
    from scipy.sparse.linalg import spsolve

    total_dim = (level_max + 1) * V
    levels = np.arange(level_max + 1)

    r0, c0 = np.nonzero(D[0])
    v0 = D[0][r0, c0]
    svc = mu * np.eye(V)
    rs, cs = np.nonzero(svc)
    vs = svc[rs, cs]

    parts = []
    # Local blocks: D0 on every level. The service outflow is put back on the
    # diagonal by the row-sum correction below.
    parts.append(_tile_block(r0, c0, v0, levels, levels, V))
    # Service: level n -> n-1 for n >= 1
    sub = levels[levels >= 1]
    parts.append(_tile_block(rs, cs, vs, sub, sub - 1, V))
    # Batch arrivals: level n -> n+k
    for k in range(1, K + 1):
        rk, ck = np.nonzero(D[k])
        if rk.size == 0:
            continue
        up = levels[levels <= level_max - k]
        if up.size == 0:
            continue
        parts.append(_tile_block(rk, ck, D[k][rk, ck], up, up + k, V))

    I = np.concatenate([pp[0] for pp in parts])
    J = np.concatenate([pp[1] for pp in parts])
    X = np.concatenate([pp[2] for pp in parts])

    Q = coo_matrix((X, (I, J)), shape=(total_dim, total_dim)).tocsr()
    row_sum = np.asarray(Q.sum(axis=1)).ravel()
    Q = Q - diags(row_sum, 0, shape=(total_dim, total_dim))

    # Solve pi * Q = 0, pi * e = 1: replace the last column with ones
    Qc = Q.tocoo()
    keep = Qc.col != (total_dim - 1)
    I2 = np.concatenate([Qc.row[keep], np.arange(total_dim)])
    J2 = np.concatenate([Qc.col[keep], np.full(total_dim, total_dim - 1)])
    X2 = np.concatenate([Qc.data[keep], np.ones(total_dim)])
    Amat = coo_matrix((X2, (I2, J2)), shape=(total_dim, total_dim))

    b = np.zeros(total_dim)
    b[-1] = 1.0
    if total_dim > 5000:
        pi_flat = spsolve(Amat.T.tocsc(), b)
    else:
        pi_flat = np.linalg.solve(Amat.toarray().T, b)

    level_prob = np.asarray(pi_flat).reshape(level_max + 1, V)
    level_prob = np.where(level_prob < 0, 0.0, level_prob)
    return level_prob / level_prob.sum()


def _tile_block(r, c, v, row_levels, col_levels, V):
    """Place a V x V block pattern at every (row_levels, col_levels) pair."""
    r = np.asarray(r).ravel()
    c = np.asarray(c).ravel()
    v = np.asarray(v, dtype=float).ravel()
    row_levels = np.asarray(row_levels).ravel()
    col_levels = np.asarray(col_levels).ravel()
    if r.size == 0 or row_levels.size == 0:
        empty_i = np.zeros(0, dtype=int)
        return empty_i, empty_i.copy(), np.zeros(0)
    I = (r[:, None] + row_levels[None, :] * V).ravel()
    J = (c[:, None] + col_levels[None, :] * V).ravel()
    X = np.tile(v[:, None], (1, row_levels.size)).ravel()
    return I, J, X
