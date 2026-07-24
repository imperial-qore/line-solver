"""MAPH(2,m) fitting: marks a second-order acyclic phase-type distribution so
that the class probabilities are matched exactly and the first-order backward
moments as closely as the form allows.

Port of m3a ``maph2m_fit.m``, ``maph2m_fit_multiclass.m``, ``maph2m_fit_mmap.m``
and ``maph2m_fit_trace.m``. An MMAP is a list ``[D0, D1, D11, ..., D1m]``.
"""
from typing import List, Optional, Sequence, Tuple

import numpy as np

from .qp import solve_qp

DEGENTOL = 1e-6
FEASTOL = 1e-8


def _isfeasible(qj: np.ndarray) -> bool:
    return bool(np.min(qj) >= -FEASTOL and np.sum(qj) <= 1 + FEASTOL)


def _fix(qj: np.ndarray) -> np.ndarray:
    out = np.maximum(qj, 0.0)
    return out / np.sum(out)


def maph2m_fit_multiclass(aph: Sequence[np.ndarray], p: Sequence[float],
                          B: Sequence[float],
                          classWeights: Optional[Sequence[float]] = None
                          ) -> Tuple[List[np.ndarray], np.ndarray]:
    """Mark a canonical acyclic APH(2).

    Args:
        aph: underlying APH(2) as [D0, D1], canonical acyclic form
        p: class probabilities
        B: target first-order backward moments
        classWeights: per-class weights in the objective (default: unit)

    Returns:
        (maph, fB): the fitted MAPH and the backward moments it realizes
    """
    D0 = np.asarray(aph[0], dtype=float)
    D1 = np.asarray(aph[1], dtype=float)
    if D0.shape[0] != 2:
        raise ValueError('Underlying APH must be of second-order.')
    if D0[1, 0] != 0:
        raise ValueError('Underlying APH must be acyclic')
    if D1[0, 1] != 0 or D1[1, 1] != 0:
        raise ValueError('Underlying APH must be in canonical acyclic form')

    p = np.asarray(p, dtype=float).ravel()
    B = np.asarray(B, dtype=float).ravel()
    k = p.size
    w = np.ones(k) if classWeights is None else np.asarray(classWeights, float).ravel()

    h1 = -1.0 / D0[0, 0]
    h2 = -1.0 / D0[1, 1]
    r1 = D0[0, 1] * h1

    q = np.zeros((2, k))
    fB = np.zeros(k)

    if abs(1 - r1) < DEGENTOL:
        # degenerate form: one degree of freedom, match the class probabilities
        q[0, :] = p
        q[1, :] = p
    else:
        # q(j,c) = fB(c) * q_b(j,c) + q_0(j,c)
        q_b = np.zeros((2, k))
        q_0 = np.zeros((2, k))
        for c in range(k):
            q_b[0, c] = p[c] * (1.0 / (h2 * (r1 - 1)))
            q_0[0, c] = p[c] * (-(h1 + h2) / (h2 * (r1 - 1)))
            q_b[1, c] = p[c] * (1.0 / (h2 * r1))
            q_0[1, c] = p[c] * (-h1 / (h2 * r1))

        A = np.zeros((4 * k, k))
        b = np.zeros(4 * k)
        for c in range(k):
            for j in range(2):
                row = c * 4 + j * 2
                A[row, c] = q_b[j, c]
                b[row] = 1 - q_0[j, c]
                A[row + 1, c] = -q_b[j, c]
                b[row + 1] = q_0[j, c]

        Aeq = np.zeros((2, k))
        beq = np.ones(2)
        for c in range(k):
            for j in range(2):
                Aeq[j, c] = q_b[j, c]
                beq[j] -= q_0[j, c]

        H = np.zeros((k, k))
        h = np.zeros(k)
        for c in range(k):
            H[c, c] = 2.0 / B[c] ** 2 * w[c]
            h[c] = -2.0 / B[c] * w[c]

        fB = solve_qp(H, h, A, b, Aeq, beq, lb=1e-6, ub=1e6, x0=B.copy())
        for c in range(k):
            for j in range(2):
                q[j, c] = fB[c] * q_b[j, c] + q_0[j, c]

    for j in range(2):
        if not _isfeasible(q[j, :]):
            raise ValueError('Fitting MAPH(2,m): Feasibility could not be restored')
        q[j, :] = _fix(q[j, :])

    maph: List[np.ndarray] = [D0.copy(), D1.copy()]
    for c in range(k):
        mask = np.array([[q[0, c], 0.0], [q[1, c], 0.0]])
        maph.append(D1 * mask)

    if abs(1 - r1) < DEGENTOL:
        from line_solver.api.mam.mmap_ops import mmap_backward_moment
        fB = np.asarray(mmap_backward_moment(maph, [1]), dtype=float).ravel()

    return maph, np.asarray(fB, dtype=float).ravel()


def maph2m_fit(M1: float, M2: float, M3: float, P: Sequence[float],
               B: Sequence[float]) -> List[np.ndarray]:
    """Second-order MAPH[m] matching three moments, the class probabilities
    (exactly) and the first-order backward moments (as closely as possible)."""
    from line_solver.api.mam.aph2_fitting import aph2_fitall, aph2_adjust

    B = np.asarray(B, dtype=float).ravel()
    P = np.asarray(P, dtype=float).ravel()

    aphs = aph2_fitall(M1, M2, M3)
    if not aphs:
        M2a, M3a = aph2_adjust(M1, M2, M3)
        aphs = aph2_fitall(M1, M2a, M3a)
    if not aphs:
        raise ValueError('Fitting MAPH(2,m): feasibility could not be restored')

    best, best_err = None, np.inf
    for aph in aphs:
        try:
            maph, fB = maph2m_fit_multiclass(aph, P, B)
        except (ValueError, RuntimeError):
            continue
        err = float(np.sum((fB / B - 1) ** 2))
        if err < best_err:
            best, best_err = maph, err
    if best is None:
        raise ValueError('Fitting MAPH(2,m): no feasible marking of the APH(2) forms')
    return best


def maph2m_fit_mmap(mmap: Sequence[np.ndarray]) -> List[np.ndarray]:
    """MAPH(2,m) fitting the characteristics of a given MMAP."""
    from line_solver.api.mam.map_analysis import map_moment
    from line_solver.api.mam.mmap_ops import mmap_backward_moment, mmap_pc

    D0 = np.asarray(mmap[0], float)
    D1 = np.asarray(mmap[1], float)
    M1 = map_moment(D0, D1, 1)
    M2 = map_moment(D0, D1, 2)
    M3 = map_moment(D0, D1, 3)
    P = np.asarray(mmap_pc(list(mmap)), dtype=float).ravel()
    B = np.asarray(mmap_backward_moment(list(mmap), [1]), dtype=float).ravel()
    return maph2m_fit(M1, M2, M3, P, B)


def maph2m_fit_trace(T: Sequence[float], A: Sequence[int]) -> List[np.ndarray]:
    """MAPH(2,m) fitting the characteristics of a marked trace.

    Args:
        T: inter-arrival times
        A: class label of each arrival
    """
    T = np.asarray(T, dtype=float).ravel()
    A = np.asarray(A).ravel()
    classes = np.unique(A)
    M1 = float(np.mean(T))
    M2 = float(np.mean(T ** 2))
    M3 = float(np.mean(T ** 3))
    P = np.array([np.mean(A == c) for c in classes], dtype=float)
    # first-order backward moment: mean inter-arrival time preceding a class-c
    # arrival, i.e. E[T_i | class(i) = c] weighted by the class probability
    from line_solver.api.trace.trace_analysis import mtrace_backward_moment
    B = np.asarray(mtrace_backward_moment(T, A, [1]), dtype=float).ravel()
    return maph2m_fit(M1, M2, M3, P, B)
