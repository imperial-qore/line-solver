"""Small convex quadratic program, the Python counterpart of the m3a QP() call.

The m3a fitters all reduce to

    min  0.5 x' H x + h' x
    s.t. A x <= b,  Aeq x = beq,  lb <= x <= ub

with H diagonal and positive, so the problem is convex and tiny (2m variables
at most). SLSQP is used because it handles equalities, inequalities and bounds
directly; the MATLAB reference uses an interior-point convex QP, which reaches
the same optimum on a convex problem.
"""
from typing import Optional

import numpy as np
from scipy.optimize import minimize


class QPFailure(RuntimeError):
    """Raised when the quadratic program cannot be solved, as MATLAB errors."""


def _independent_rows(Aeq: np.ndarray, beq: np.ndarray, tol: float = 1e-10):
    """Keep a maximal independent subset of the equality rows.

    A dependent row is redundant when it is consistent with the ones kept, and
    signals an infeasible system when it is not.
    """
    kept_rows = []
    kept_rhs = []
    basis = np.zeros((0, Aeq.shape[1]))
    for i in range(Aeq.shape[0]):
        row = Aeq[i, :]
        nrm = np.linalg.norm(row)
        if nrm <= tol:
            if abs(beq[i]) > tol * max(1.0, abs(beq[i])):
                raise QPFailure('Quadratic programming solver failed: inconsistent equality')
            continue
        trial = np.vstack([basis, row])
        if np.linalg.matrix_rank(trial, tol=tol * max(1.0, np.linalg.norm(trial))) > basis.shape[0]:
            basis = trial
            kept_rows.append(row)
            kept_rhs.append(beq[i])
        else:
            # dependent: verify consistency against the retained rows
            if kept_rows:
                K = np.array(kept_rows)
                coeff, *_ = np.linalg.lstsq(K.T, row, rcond=None)
                if abs(coeff @ np.array(kept_rhs) - beq[i]) > 1e-8 * max(1.0, abs(beq[i])):
                    raise QPFailure('Quadratic programming solver failed: inconsistent equality')
    if not kept_rows:
        return np.zeros((0, Aeq.shape[1])), np.zeros(0)
    return np.array(kept_rows), np.array(kept_rhs)


def solve_qp(H: np.ndarray, h: np.ndarray,
             A: Optional[np.ndarray] = None, b: Optional[np.ndarray] = None,
             Aeq: Optional[np.ndarray] = None, beq: Optional[np.ndarray] = None,
             lb: float = 1e-6, ub: float = 1e6, x0: Optional[np.ndarray] = None
             ) -> np.ndarray:
    H = np.atleast_2d(np.asarray(H, dtype=float))
    h = np.asarray(h, dtype=float).ravel()
    n = h.size

    def obj(x):
        return 0.5 * x @ H @ x + h @ x

    def jac(x):
        return H @ x + h

    cons = []
    if A is not None and np.size(A) > 0:
        A = np.atleast_2d(np.asarray(A, dtype=float))
        b = np.asarray(b, dtype=float).ravel()
        cons.append({'type': 'ineq',
                     'fun': lambda x, A=A, b=b: b - A @ x,
                     'jac': lambda x, A=A: -A})
    if Aeq is not None and np.size(Aeq) > 0:
        Aeq = np.atleast_2d(np.asarray(Aeq, dtype=float))
        beq = np.asarray(beq, dtype=float).ravel()
        # The m3a equality blocks are frequently rank deficient (the per-phase
        # normalizations coincide), which makes the SLSQP least-squares
        # subproblem singular. MATLAB's interior-point QP tolerates that, so
        # drop the dependent rows here and keep only an independent subset.
        Aeq, beq = _independent_rows(Aeq, beq)
        if Aeq.size:
            cons.append({'type': 'eq',
                         'fun': lambda x, Aeq=Aeq, beq=beq: Aeq @ x - beq,
                         'jac': lambda x, Aeq=Aeq: Aeq})

    bounds = [(lb, ub)] * n
    if x0 is None:
        # unconstrained minimizer of the separable objective, clipped into the box
        with np.errstate(divide='ignore', invalid='ignore'):
            diag = np.diag(H).copy()
            diag[diag == 0] = 1.0
            x0 = np.clip(-h / diag, lb, ub)
    res = minimize(obj, np.asarray(x0, dtype=float), jac=jac, bounds=bounds,
                   constraints=cons, method='SLSQP',
                   options={'maxiter': 3000, 'ftol': 1e-12})
    if not res.success:
        # one retry from the centre of the box before giving up, as the MATLAB
        # solver also restarts internally
        res = minimize(obj, np.full(n, min(max(1.0, lb), ub)), jac=jac,
                       bounds=bounds, constraints=cons, method='SLSQP',
                       options={'maxiter': 5000, 'ftol': 1e-12})
    if not res.success:
        raise QPFailure('Quadratic programming solver failed: %s' % res.message)
    return res.x
