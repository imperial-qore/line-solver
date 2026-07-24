"""
Restarted GMRES with ILUT preconditioning for sparse linear systems.

The iterative counterpart of the direct sparse solve used by ctmc_solve,
intended for generators whose LU fill-in exceeds available memory. The kernel
performs no CTMC-specific processing, so it also serves the stochastic
complementation and aggregation kernels.

Two preparation steps are not optional on a generator. Rows are equilibrated to
unit max norm, so the O(1) normalization row does not mix with rows carrying
rates of a different magnitude. The states are then reordered by reverse
Cuthill-McKee: in the natural ordering of a birth-death chain the unpivoted
elimination has growth factor (mu/lambda)^n, which overflows by a few thousand
states, and a fill-reducing ordering rather than pivoting is what removes it.

Key algorithms:
    ctmc_gmres: Restarted GMRES with an ILUT right preconditioner
    ctmc_gmres_multi: The same, over every column of a right-hand-side block
"""

from typing import Optional, Tuple

import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import reverse_cuthill_mckee
from scipy.sparse.linalg import LinearOperator, gmres, spilu

# see _kb/03-api-layer.md for rationale
GMRES_DEFAULT_TOL = 1e-12

# Default Krylov subspace dimension between restarts.
GMRES_DEFAULT_RESTART = 50

# Relative threshold below which ILUT discards a fill-in entry, and the factor
# bounding the ILUT factors as a multiple of the nonzeros of A.
ILUT_DROP_TOL = 1e-4
ILUT_FILL_FACTOR = 10.0


def _rcm_permutation(A: sp.csc_matrix) -> np.ndarray:
    """
    Reverse Cuthill-McKee ordering of the symmetrized sparsity pattern, as an
    array mapping a new index to the old one.
    """
    pattern = sp.csr_matrix(abs(A) + abs(A).T, dtype=np.int8)
    return reverse_cuthill_mckee(pattern, symmetric_mode=True)


def _build_preconditioner(A: sp.csc_matrix, n: int) -> LinearOperator:
    """
    ILUT right preconditioner, falling back to Jacobi on breakdown.

    The matrix is already reordered by the caller, so no further fill-reducing
    permutation is applied here: NATURAL keeps the three codebases factorizing
    the same matrix in the same order. A breakdown is common on generators
    holding an absorbing or a rate-free state, hence the guarded factorization.
    """
    try:
        ilu = spilu(
            A, drop_tol=ILUT_DROP_TOL, fill_factor=ILUT_FILL_FACTOR,
            permc_spec='NATURAL', diag_pivot_thresh=0.0,
        )
        return LinearOperator((n, n), matvec=ilu.solve)
    except (RuntimeError, ValueError):
        # Jacobi fallback. A zero diagonal entry would make the preconditioner
        # singular, so those rows are left unscaled rather than inverted.
        d = np.asarray(A.diagonal(), dtype=np.float64).copy()
        d[d == 0.0] = 1.0
        dinv = 1.0 / d
        return LinearOperator((n, n), matvec=lambda v: dinv * v)


def ctmc_gmres(
    A,
    b,
    tol: Optional[float] = None,
    restart: Optional[int] = None,
    maxit: Optional[int] = None,
    x0=None,
) -> Tuple[np.ndarray, int, float, int]:
    """
    Solve the sparse nonsymmetric system A*x = b by restarted GMRES.

    Args:
        A: Coefficient matrix, dense or sparse. Converted to CSC internally.
        b: Right-hand side
        tol: Relative residual tolerance (default 1e-12)
        restart: Krylov subspace dimension between restarts (default min(n, 50))
        maxit: Maximum number of restart cycles (default ceil(n/restart))
        x0: Initial guess (default uniform 1/n)

    Returns:
        (x, flag, relres, iter), where flag follows the MATLAB gmres
        convention: 0 converged, 1 iteration limit reached, 2 preconditioner
        ill-conditioned, 3 stagnation or breakdown. Callers must check flag and
        fall back to the direct solve when it is nonzero.
    """
    A = sp.csc_matrix(A, dtype=np.float64)
    n = A.shape[0]
    b = np.asarray(b, dtype=np.float64).reshape(n)

    if tol is None or tol <= 0.0:
        tol = GMRES_DEFAULT_TOL
    if restart is None or restart <= 0:
        restart = min(n, GMRES_DEFAULT_RESTART)
    restart = min(restart, n)
    if maxit is None or maxit <= 0:
        maxit = int(np.ceil(n / restart))
    maxit = max(1, min(maxit, n))

    if x0 is None:
        x0 = np.full(n, 1.0 / n)
    else:
        x0 = np.asarray(x0, dtype=np.float64).reshape(n)

    # see _kb/03-api-layer.md for rationale
    rownorm = np.asarray(abs(A).max(axis=1).todense()).ravel()
    rownorm[rownorm == 0.0] = 1.0
    A = (sp.diags(1.0 / rownorm) @ A).tocsc()
    b = b / rownorm

    perm = _rcm_permutation(A)
    A = A[perm][:, perm].tocsc()
    b = b[perm]
    x0 = x0[perm]

    M = _build_preconditioner(A, n)

    # Count inner iterations so the three codebases report a single comparable
    # scalar; scipy reports only the outer cycle count through info.
    counter = {'n': 0}

    def _callback(_res):
        counter['n'] += 1

    try:
        x, info = gmres(
            A, b, x0=x0, rtol=tol, atol=0.0, restart=restart,
            maxiter=maxit, M=M, callback=_callback, callback_type='pr_norm',
        )
    except (RuntimeError, ValueError):
        # A breakdown inside GMRES is reported as non-convergence rather than
        # propagated, so the caller falls back to the direct solve.
        out = np.empty(n)
        out[perm] = x0
        return out, 3, np.inf, 0

    bnorm = float(np.linalg.norm(b))
    if bnorm == 0.0:
        bnorm = 1.0
    relres = float(np.linalg.norm(b - A @ x) / bnorm)

    if info == 0:
        flag = 0
    elif info > 0:
        flag = 1
    else:
        flag = 3

    if not np.all(np.isfinite(x)):
        flag = 3
        relres = np.inf

    out = np.empty(n)
    out[perm] = x
    return out, flag, relres, counter['n']


def ctmc_gmres_multi(
    A,
    B,
    tol: Optional[float] = None,
    restart: Optional[int] = None,
    maxit: Optional[int] = None,
):
    """
    Solve A*X = B for every column of B, reusing one ILUT factorization across
    all of them and starting each column from the previous solution.

    This is the shape of the stochastic complement, whose right-hand side is a
    whole block of the generator: refactorizing per column would cost more than
    the direct solve it replaces.

    Args:
        A: Coefficient matrix, dense or sparse
        B: Right-hand sides, one per column
        tol: Relative residual tolerance (default 1e-12)
        restart: Krylov subspace dimension between restarts
        maxit: Maximum number of restart cycles

    Returns:
        (X, flag). flag is 0 only if every column converged; on any other value
        X is None and the caller must fall back to the direct solve. Returning a
        partial block would leave that fallback ambiguous.
    """
    A = sp.csc_matrix(A, dtype=np.float64)
    n = A.shape[0]
    B = np.asarray(B, dtype=np.float64)
    if B.ndim == 1:
        B = B.reshape(n, 1)

    if tol is None or tol <= 0.0:
        tol = GMRES_DEFAULT_TOL
    if restart is None or restart <= 0:
        restart = min(n, GMRES_DEFAULT_RESTART)
    restart = min(restart, n)
    if maxit is None or maxit <= 0:
        maxit = int(np.ceil(n / restart))
    maxit = max(1, min(maxit, n))

    # Same preparation as ctmc_gmres, hoisted out of the column loop.
    rownorm = np.asarray(abs(A).max(axis=1).todense()).ravel()
    rownorm[rownorm == 0.0] = 1.0
    A = (sp.diags(1.0 / rownorm) @ A).tocsc()
    B = B / rownorm[:, None]

    perm = _rcm_permutation(A)
    A = A[perm][:, perm].tocsc()
    B = B[perm, :]

    M = _build_preconditioner(A, n)

    Xp = np.zeros((n, B.shape[1]))
    guess = np.full(n, 1.0 / n)
    for c in range(B.shape[1]):
        try:
            x, info = gmres(
                A, B[:, c], x0=guess, rtol=tol, atol=0.0, restart=restart,
                maxiter=maxit, M=M,
            )
        except (RuntimeError, ValueError):
            return None, 3
        if info != 0 or not np.all(np.isfinite(x)):
            return None, 1 if info > 0 else 3
        Xp[:, c] = x
        guess = x

    X = np.zeros_like(Xp)
    X[perm, :] = Xp
    return X, 0
