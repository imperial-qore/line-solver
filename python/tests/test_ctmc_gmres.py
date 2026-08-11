"""
Analytic validation of ctmc_gmres.

The oracle is the M/M/1/K stationary distribution, a truncated geometric, rather
than a recorded baseline: a recorded value cannot catch an error that the
reference implementation shares.
"""

import numpy as np
import pytest
import scipy.sparse as sp

from line_solver.api.mc import ctmc_gmres


def mm1k_generator(lam, mu, K):
    """Generator of an M/M/1/K queue."""
    n = K + 1
    Q = sp.lil_matrix((n, n))
    for i in range(n):
        if i < n - 1:
            Q[i, i + 1] = lam
        if i > 0:
            Q[i, i - 1] = mu
    Q.setdiag(-np.asarray(Q.sum(axis=1)).ravel())
    return Q.tocsc()


def normalized_system(Q):
    """
    Assembles the linear system ctmc_solve poses: the last column of Q is
    replaced by ones to carry the normalization, and the transposed system is
    solved against e_n.
    """
    n = Q.shape[0]
    A = Q.tolil()
    A[:, -1] = 1.0
    b = np.zeros(n)
    b[-1] = 1.0
    return sp.csc_matrix(A.T), b


@pytest.mark.parametrize('rho', [0.7, 1.0])
@pytest.mark.parametrize('K', [10, 20000])
def test_gmres_matches_truncated_geometric(rho, K):
    # K = 20000 is past the size at which a direct factorization is the intended
    # method, and rho < 1 is the case in which an unpreconditioned or naturally
    # ordered elimination overflows.
    A, b = normalized_system(mm1k_generator(rho, 1.0, K))
    x, flag, relres, _ = ctmc_gmres(A, b)
    assert flag == 0, 'gmres flag %d (relres %g)' % (flag, relres)

    pi = rho ** np.arange(K + 1)
    pi = pi / pi.sum()
    assert np.max(np.abs(x - pi)) < 1e-9


def test_gmres_agrees_with_direct_solve():
    # The dispatch in ctmc_solve must be numerically invisible, so the two
    # methods have to agree far below any tolerance a caller would notice.
    from scipy.sparse.linalg import spsolve

    for K, rho in [(3, 0.5), (10, 0.85), (200, 0.85)]:
        A, b = normalized_system(mm1k_generator(rho, 1.0, K))
        xd = spsolve(A, b)
        x, flag, _, _ = ctmc_gmres(A, b)
        assert flag == 0
        assert np.max(np.abs(x - xd)) < 1e-9


def test_gmres_reports_non_convergence():
    # A caller may only trust the answer when the flag is zero. One restart cycle
    # of dimension one cannot converge on a 500-state chain, and the kernel has
    # to say so rather than return the iterate it happens to hold.
    A, b = normalized_system(mm1k_generator(0.9, 1.0, 500))
    _, flag, relres, _ = ctmc_gmres(A, b, tol=1e-12, restart=1, maxit=1)
    assert flag != 0
    assert relres > 1e-12


def test_ctmc_solve_dispatch_is_numerically_invisible():
    # The two methods must agree far below any tolerance a caller would notice,
    # and the answer must not jump as a model grows past the threshold.
    from line_solver.api.mc import ctmc_solve

    for K in [50, 400]:
        Q = mm1k_generator(0.8, 1.0, K).toarray()
        pd = ctmc_solve(Q, method='direct')
        pg = ctmc_solve(Q, method='gmres')
        assert np.max(np.abs(pd - pg)) < 1e-9

        pi = 0.8 ** np.arange(K + 1)
        pi = pi / pi.sum()
        assert np.max(np.abs(pg - pi)) < 1e-9


def test_gmres_survives_zero_diagonal():
    # A generator with no diagonal entries breaks the incomplete factorization.
    # The kernel must fall back rather than propagate the breakdown.
    Z = sp.csc_matrix((np.ones(3), ([0, 1, 2], [1, 2, 0])), shape=(3, 3))
    b = np.zeros(3)
    b[0] = 1.0
    x, _, _, _ = ctmc_gmres(Z, b)
    assert np.all(np.isfinite(x))
