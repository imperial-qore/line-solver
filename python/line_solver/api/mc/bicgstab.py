"""
Preconditioned stabilized biconjugate gradients for sparse linear systems.

The short-recurrence counterpart of ctmc_gmres: work and storage per iteration
are constant rather than growing with the Krylov dimension, so the method does
not restart and does not lose the optimality that restarting costs GMRES. Where
GMRES(m) stagnates because the useful subspace is wider than m, this converges;
where it does not, GMRES(m) is the more robust of the two, hence the order in
which ctmc_solve tries them.

The equilibration, reverse Cuthill-McKee reordering and ILUT preconditioner are
imported from the GMRES kernel rather than reimplemented, so both methods
factorize the same matrix in the same order and a switch between them cannot
move a reported metric for a reason other than the iteration itself.

Key algorithms:
    ctmc_bicgstab: BiCGSTAB of van der Vorst (1992) with an ILUT preconditioner
    ctmc_bicgstab_multi: The same, over every column of a right-hand-side block
"""

from typing import Optional, Tuple

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import bicgstab

from .gmres import _build_preconditioner, _rcm_permutation

# Linear-solve residual, as in the GMRES kernel. Much tighter than the
# fixed-point tolerance options.iter_tol: switching solve method must not move a
# reported metric.
BICGSTAB_DEFAULT_TOL = 1e-12

# Default cap on complete iterations. BiCGSTAB storage is O(n) regardless of the
# count, so the cap bounds time rather than memory.
BICGSTAB_DEFAULT_MAXIT = 200


def ctmc_bicgstab(
    A,
    b,
    tol: Optional[float] = None,
    maxit: Optional[int] = None,
    x0=None,
) -> Tuple[np.ndarray, int, float, int]:
    """
    Solve the sparse nonsymmetric system A*x = b by preconditioned BiCGSTAB.

    Args:
        A: Coefficient matrix, dense or sparse. Converted to CSC internally.
        b: Right-hand side
        tol: Relative residual tolerance (default 1e-12)
        maxit: Maximum number of complete iterations (default min(n, 200))
        x0: Initial guess (default uniform 1/n)

    Returns:
        (x, flag, relres, iter), where flag follows the MATLAB bicgstab
        convention: 0 converged, 1 iteration limit reached, 2 preconditioner
        ill-conditioned, 3 stagnation, 4 a scalar quantity became too small or
        too large to continue. Callers must check flag and fall back to another
        solve when it is nonzero. iter counts matrix-vector products with A: two
        per complete iteration, which is what makes it comparable with the iter
        of ctmc_gmres and across the four codebases. scipy reports only complete
        iterations, so a solve that converges at a half step is counted here as
        the full pair.
    """
    A = sp.csc_matrix(A, dtype=np.float64)
    n = A.shape[0]
    b = np.asarray(b, dtype=np.float64).reshape(n)

    if tol is None or tol <= 0.0:
        tol = BICGSTAB_DEFAULT_TOL
    if maxit is None or maxit <= 0:
        maxit = min(n, BICGSTAB_DEFAULT_MAXIT)
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

    # scipy fires the callback once per complete iteration, which is two
    # matrix-vector products with A.
    counter = {'n': 0}

    def _callback(_xk):
        counter['n'] += 1

    try:
        x, info = bicgstab(
            A, b, x0=x0, rtol=tol, atol=0.0, maxiter=maxit, M=M,
            callback=_callback,
        )
    except (RuntimeError, ValueError):
        # A breakdown inside BiCGSTAB is reported as non-convergence rather than
        # propagated, so the caller falls back to another solve.
        out = np.empty(n)
        out[perm] = x0
        return out, 4, np.inf, 0

    bnorm = float(np.linalg.norm(b))
    if bnorm == 0.0:
        bnorm = 1.0
    relres = float(np.linalg.norm(b - A @ x) / bnorm)

    if info == 0:
        flag = 0
    elif info > 0:
        flag = 1
    else:
        # A negative info is a breakdown of the underlying Lanczos process, not
        # slow convergence.
        flag = 4

    if not np.all(np.isfinite(x)):
        flag = 4
        relres = np.inf

    out = np.empty(n)
    out[perm] = x
    return out, flag, relres, 2 * counter['n']


def ctmc_bicgstab_multi(
    A,
    B,
    tol: Optional[float] = None,
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
        maxit: Maximum number of complete iterations per column

    Returns:
        (X, flag). flag is 0 only if every column converged; on any other value
        X is None and the caller must fall back to another solve. Returning a
        partial block would leave that fallback ambiguous.
    """
    A = sp.csc_matrix(A, dtype=np.float64)
    n = A.shape[0]
    B = np.asarray(B, dtype=np.float64)
    if B.ndim == 1:
        B = B.reshape(n, 1)

    if tol is None or tol <= 0.0:
        tol = BICGSTAB_DEFAULT_TOL
    if maxit is None or maxit <= 0:
        maxit = min(n, BICGSTAB_DEFAULT_MAXIT)
    maxit = max(1, min(maxit, n))

    # Same preparation as ctmc_bicgstab, hoisted out of the column loop.
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
            x, info = bicgstab(
                A, B[:, c], x0=guess, rtol=tol, atol=0.0, maxiter=maxit, M=M,
            )
        except (RuntimeError, ValueError):
            return None, 4
        if info != 0 or not np.all(np.isfinite(x)):
            return None, 1 if info > 0 else 4
        Xp[:, c] = x
        guess = x

    X = np.zeros_like(Xp)
    X[perm, :] = Xp
    return X, 0
