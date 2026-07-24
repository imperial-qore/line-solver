"""
Contract tests for pfqn_mvams, the general-purpose MVA dispatcher for mixed
networks with multiserver nodes.

pfqn_mvams routes to one of four algorithms depending on whether the model has
open classes and whether any station is multiserver:

    single-server + closed  -> pfqn_mva
    single-server + mixed   -> pfqn_mvamx
    multiserver   + mixed   -> pfqn_mvaldms
    multiserver   + closed  -> pfqn_mvald

Those four do NOT share a return contract of their own: the load-dependent
family (pfqn_mvald) reports a (1 x R) cycle time, while pfqn_mva reports the
(1 x R) cycle time as CN and the (M x R) residence time as RN. pfqn_mvams must
nonetheless present ONE contract on every branch, namely MATLAB pfqn_mvams':
an (M x R) per-station residence time. These tests pin that, because a branch
returning a different QUANTITY under the same name is invisible to any test that
only exercises one branch.
"""

import numpy as np
import pytest

from line_solver.api.pfqn import pfqn_mvams, pfqn_mva, pfqn_mvald

TOL = 1e-9

L = np.array([[0.4, 0.7], [0.9, 0.3]])
Z = np.array([1.0, 0.5])

# (name, lambda, N, S) covering all four dispatch branches
BRANCHES = [
    ('closed single-server (pfqn_mva)', [0.0, 0.0], [2.0, 2.0], [1, 1]),
    ('closed multiserver (pfqn_mvald)', [0.0, 0.0], [2.0, 2.0], [2, 2]),
    ('mixed single-server (pfqn_mvamx)', [0.2, 0.0], [np.inf, 2.0], [1, 1]),
    ('mixed multiserver (pfqn_mvaldms)', [0.2, 0.0], [np.inf, 2.0], [2, 2]),
    ('closed multiserver, empty class', [0.0, 0.0], [3.0, 0.0], [2, 2]),
    ('closed multiserver, 3 servers', [0.0, 0.0], [3.0, 2.0], [3, 2]),
]


@pytest.mark.parametrize('name,lam,N,S', BRANCHES)
def test_residence_time_contract_on_every_branch(name, lam, N, S):
    """CN must be an (M x R) residence time satisfying QN = XN*CN, on every branch."""
    X, Q, U, C, lG = pfqn_mvams(np.array(lam), L, np.array(N), Z,
                                mi=np.ones(2), S=np.array(S))
    X = np.asarray(X).flatten()
    C = np.asarray(C)
    assert C.shape == L.shape, 'CN must be (M x R) on branch: ' + name
    assert not np.any(np.isnan(C)), 'CN must not be NaN on branch: ' + name
    for r in range(L.shape[1]):
        if np.isinf(N[r]) or N[r] == 0:
            continue
        # the defining identity of a residence time
        np.testing.assert_allclose(Q[:, r], X[r] * C[:, r], atol=TOL,
                                   err_msg='QN != XN*CN on branch: ' + name)


def test_single_server_limit_equals_pfqn_mva():
    """With S=1 the dispatcher takes the pfqn_mva branch and must agree with it."""
    X, Q, U, C, lG = pfqn_mvams(np.array([0.0, 0.0]), L, np.array([2.0, 2.0]), Z,
                                mi=np.ones(2), S=np.array([1, 1]))
    ref = pfqn_mva(L, np.array([2, 2]), Z)
    Xr, Qr, RNr = np.asarray(ref[0]).flatten(), np.asarray(ref[2]), np.asarray(ref[4])
    np.testing.assert_allclose(np.asarray(X).flatten(), Xr, atol=TOL)
    np.testing.assert_allclose(Q, Qr, atol=TOL)
    # pfqn_mva's RESIDENCE time is RN, not its CN (which is the cycle time)
    np.testing.assert_allclose(C, RNr, atol=TOL)


def test_multiserver_branch_still_matches_mvald():
    """The fix must not disturb X and Q, which come straight from pfqn_mvald."""
    N = np.array([2, 2])
    mu = np.tile(np.minimum(np.arange(1, N.sum() + 1), 2), (2, 1)).astype(float)
    X, Q, U, C, lG = pfqn_mvams(np.array([0.0, 0.0]), L, N.astype(float), Z,
                                mi=np.ones(2), S=np.array([2, 2]))
    ref = pfqn_mvald(L, N, Z, mu)
    Xr, Qr, Cr = (np.asarray(ref[0]).flatten(), np.asarray(ref[1]),
                  np.asarray(ref[3]).flatten())
    np.testing.assert_allclose(np.asarray(X).flatten(), Xr, atol=TOL)
    np.testing.assert_allclose(Q, Qr, atol=TOL)
    # the residence times must decompose pfqn_mvald's cycle time
    np.testing.assert_allclose(np.asarray(C).sum(axis=0), Cr, atol=TOL)
