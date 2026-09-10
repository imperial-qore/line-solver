"""
The two Krylov kernels of the CTMC path, against an ANALYTICAL oracle.

ctmc_gmres and ctmc_bicgstab share the same equilibration, reverse Cuthill-McKee
reordering and ILUT preconditioner, so checking one against the other proves
nothing about either: a fault in the shared preparation would cancel. The oracle
here is the stationary vector of M/M/1/K, which is the truncated geometric in
closed form, plus the direct solve on a two-stage tandem lattice.

Twin of the ctmc_gmres/ctmc_bicgstab cases in cpp/tests/test_mc_aggregation.cpp
and of jar/src/test/java/jline/api/mc/CtmcKrylovTest.java.
"""

import os

import numpy as np
import pytest
import scipy.sparse as sp

import line_solver
# The worktree copy of line_solver must be the one under test; a global install
# would silently validate the wrong code.
assert os.path.realpath(__file__).rsplit('/python/', 1)[0] in os.path.realpath(
    line_solver.__file__), (
    "test must import the worktree line_solver, got %s" % line_solver.__file__)

from line_solver.api.mc import (ctmc_bicgstab, ctmc_bicgstab_multi, ctmc_gmres,
                                ctmc_gmres_multi)


def _mm1k(K, lam=0.7, mu=1.0):
    """The system ctmc_solve builds for M/M/1/K, and its closed-form answer."""
    n = K + 1
    main = np.empty(n)
    main[0] = -lam
    main[-1] = -mu
    main[1:-1] = -(lam + mu)
    Q = sp.diags([np.full(n - 1, mu), main, np.full(n - 1, lam)], [-1, 0, 1], format='lil')
    # Replace the last column by the normalization row, then transpose.
    Q[:, -1] = 1.0
    A = Q.tocsc().T.tocsc()
    b = np.zeros(n)
    b[-1] = 1.0
    rho = lam / mu
    exact = rho ** np.arange(n)
    exact /= exact.sum()
    return A, b, exact


@pytest.mark.parametrize('K', [200, 5000])
@pytest.mark.parametrize('kernel', [ctmc_gmres, ctmc_bicgstab])
def test_kernels_reproduce_the_truncated_geometric(kernel, K):
    A, b, exact = _mm1k(K)
    x, flag, relres, iters = kernel(A, b)
    assert flag == 0, 'flag %d, relative residual %g' % (flag, relres)
    assert np.max(np.abs(x - exact)) < 1e-11
    assert iters > 0


def test_bicgstab_reports_matvecs_not_iterations():
    """iter counts matrix-vector products with A, two per complete iteration.

    Counting products rather than iterations is what makes the number comparable
    with ctmc_gmres and across the four codebases, whose iteration bookkeeping
    differs (MATLAB counts half iterations, scipy counts full ones). Only the
    bound is assertable: an implementation that converges at a half step reports
    an odd count.
    """
    maxit = 10
    A, b, _ = _mm1k(200)
    _, flag, _, iters = ctmc_bicgstab(A, b, None, maxit)
    assert flag == 0
    assert 1 <= iters <= 2 * maxit


def _tandem(m, lam=0.5, mu=1.0):
    """Two-stage tandem lattice with m buffer slots per stage."""
    n = m * m
    Q = sp.lil_matrix((n, n))
    idx = lambda i, j: i * m + j
    for i in range(m):
        for j in range(m):
            s = idx(i, j)
            if i < m - 1:
                Q[s, idx(i + 1, j)] += lam
            if i > 0 and j < m - 1:
                Q[s, idx(i - 1, j + 1)] += mu
            if j > 0:
                Q[s, idx(i, j - 1)] += mu
    Q.setdiag(0.0)
    Q.setdiag(-np.asarray(Q.tocsr().sum(axis=1)).ravel())
    Q = Q.tolil()
    Q[:, -1] = 1.0
    A = Q.tocsc().T.tocsc()
    b = np.zeros(n)
    b[-1] = 1.0
    return A, b


@pytest.mark.parametrize('kernel', [ctmc_gmres, ctmc_bicgstab])
def test_kernels_match_the_direct_solve_on_a_tandem_lattice(kernel):
    A, b = _tandem(31)
    reference = np.linalg.solve(A.toarray(), b)
    x, flag, relres, _ = kernel(A, b)
    assert flag == 0, 'flag %d, relative residual %g' % (flag, relres)
    assert np.max(np.abs(x - reference)) < 1e-9


@pytest.mark.parametrize('kernel', [ctmc_gmres_multi, ctmc_bicgstab_multi])
def test_multi_column_kernels_solve_every_column(kernel):
    """The shape of the stochastic complement: one factorization, many columns."""
    A, _ = _tandem(21)
    rng = np.random.default_rng(23000)
    B = rng.standard_normal((A.shape[0], 3))
    X, flag = kernel(A, B, 1e-10)
    assert flag == 0
    assert X is not None
    assert np.max(np.abs(A @ X - B)) < 1e-7
