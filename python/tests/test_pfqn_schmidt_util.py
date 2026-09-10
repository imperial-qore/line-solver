"""
Regression: pfqn_schmidt must return a utilization, not zeros.

pfqn_schmidt declared an (M x R) UN, allocated a `u` accumulator that was never
written, and returned it -- MATLAB even carried the comment "this will return 0".
Python mirrored the omission. The JAR was the only codebase that computed it, via
the Utilization Law U = D*X/c, so the JAR was taken as the reference.

pfqn_schmidt_ext computed U from a throughput vector that had already been
replicated over the stations, so in MATLAB the linear index XN(c) silently read
the WRONG class's throughput whenever M > 1. Python and the JAR index the first
row explicitly and were correct.

Oracle: the Utilization Law U = D*X/c is exact for these product-form networks,
and pfqn_schmidt already returns X and Q, so U is over-determined by its own
outputs.
"""

import numpy as np
import pytest

from line_solver.api.pfqn import pfqn_schmidt, pfqn_schmidt_ext

TOL = 1e-9

D = np.array([[0.4, 0.7], [0.9, 0.3]])
N = np.array([2, 2])
SCHED_FCFS = np.array([0, 0])


def _expected(D, X_row, S):
    M, R = D.shape
    return np.array([[D[m, c] * X_row[c] / S[m] for c in range(R)] for m in range(M)])


@pytest.mark.parametrize('S', [np.array([2, 2]), np.array([1, 1]), np.array([3, 2])])
def test_schmidt_util_follows_utilization_law(S):
    X, Q, U, C = pfqn_schmidt(D, N, S, SCHED_FCFS)
    assert not np.all(U == 0), 'UN must not be all zeros'
    np.testing.assert_allclose(U, _expected(D, X[0, :], S), atol=TOL)


@pytest.mark.parametrize('S', [np.array([2, 2]), np.array([1, 1]), np.array([3, 2])])
def test_schmidt_ext_util_follows_utilization_law(S):
    X, Q, U, C = pfqn_schmidt_ext(D, N, S, SCHED_FCFS)
    assert not np.all(U == 0), 'UN must not be all zeros'
    np.testing.assert_allclose(U, _expected(D, X[0, :], S), atol=TOL)


def test_schmidt_ext_uses_the_right_class_throughput():
    """
    The classes must not share a throughput. With distinct per-class throughputs,
    reading XN by a linear index into the (M x R) replica would give class 2 the
    throughput of class 1; this pins that it does not.
    """
    S = np.array([2, 2])
    X, Q, U, C = pfqn_schmidt_ext(D, N, S, SCHED_FCFS)
    x = X[0, :]
    assert not np.isclose(x[0], x[1]), 'fixture must have distinct class throughputs'
    # U(0,1) must use x[1], not x[0]
    assert U[0, 1] == pytest.approx(D[0, 1] * x[1] / S[0], abs=TOL)
    assert not np.isclose(U[0, 1], D[0, 1] * x[0] / S[0])


def test_schmidt_single_server_is_plain_offered_load():
    """With one server the Utilization Law degenerates to U = X*D."""
    S = np.array([1, 1])
    X, Q, U, C = pfqn_schmidt(D, N, S, SCHED_FCFS)
    np.testing.assert_allclose(U, D * X[0, :][None, :], atol=TOL)
