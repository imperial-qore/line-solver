"""
Validation of the Section V extension of MVAC to queue-length dependent (QLD)
service centers (Conway, de Souza e Silva and Lavenberg, IEEE TC 38(3), 1989).

Four independent oracles are used, since MVAC-LD is exact and must not merely be
self-consistent:
  1. mu = 1 must reduce EXACTLY to pfqn_mvac (Sections II-IV), a structurally
     different recursion that propagates means rather than distributions.
  2. Load-dependent results must match pfqn_mvald, the classic load-dependent MVA
     population recursion.
  3. Eq. (25) is self-normalizing, so sum_n P_j(n) = 1 identically; a deviation
     means the implementation is broken, not the model.
  4. The marginals must reproduce the means: Q_j = sum_n n P_j(n).
The three chain-resolution paths (part 1 / part 2 / label interchange) are
covered explicitly, as in test_pfqn_mvac.
"""

import numpy as np
import pytest

from line_solver.api.pfqn import pfqn_mvacld, pfqn_mvac, pfqn_mvald

EXACT_TOL = 1e-9


def _ms_rates(M, Nt, c):
    """Rates of a c-server queue, mu(n) = min(n,c)."""
    return np.tile(np.minimum(np.arange(1, Nt + 1), c), (M, 1)).astype(float)


def _check_internal_oracles(L, N, Z, mu):
    L = np.array(L, float)
    N = np.array(N, int)
    Z = np.atleast_2d(np.array(Z, float))
    X, Q, U, C, pij = pfqn_mvacld(L, N, Z, mu)
    Nt = int(N.sum())
    # (3) eq. (25) is self-normalizing
    np.testing.assert_allclose(pij.sum(axis=1), 1.0, atol=EXACT_TOL)
    # (4) the marginals reproduce the means
    np.testing.assert_allclose(Q.sum(axis=1), pij @ np.arange(Nt + 1), atol=EXACT_TOL)
    # utilization is 1-P_j(0), per station
    np.testing.assert_allclose(U, 1.0 - pij[:, 0], atol=EXACT_TOL)
    # every customer is accounted for
    np.testing.assert_allclose(Q.sum() + (X * Z.sum(axis=0)).sum(), N.sum(), atol=1e-8)
    return X, Q, U, C, pij


def _assert_matches_mvald(L, N, Z, mu):
    # every output must match pfqn_mvald, not just X and Q: the LD family
    # contract is X (1xR), Q (MxR), U (Mx1) = 1-P_j(0), C (1xR) cycle time
    X, Q, U, C, pij = _check_internal_oracles(L, N, Z, mu)
    Z = np.atleast_2d(np.array(Z, float))
    ref = pfqn_mvald(np.array(L, float), np.array(N, int), Z.sum(axis=0), mu)
    Xref, Qref, Uref, Cref = (np.asarray(ref[0]).flatten(), np.asarray(ref[1]),
                              np.asarray(ref[2]), np.asarray(ref[3]))
    np.testing.assert_allclose(X, Xref, atol=EXACT_TOL)
    np.testing.assert_allclose(Q, Qref, atol=EXACT_TOL)
    assert U.shape == Uref.flatten().shape, 'U shape must match pfqn_mvald'
    np.testing.assert_allclose(U, Uref.flatten(), atol=EXACT_TOL)
    assert C.shape == Cref.flatten().shape, 'C shape must match pfqn_mvald'
    np.testing.assert_allclose(C, Cref.flatten(), atol=EXACT_TOL)


# ---- oracle 1: mu = 1 must reduce exactly to pfqn_mvac ----

@pytest.mark.parametrize('L,N,Z', [
    ([[0.4, 0.7, 0.2], [0.9, 0.3, 0.5]], [1, 2, 1], [[1.0, 0.5, 2.0]]),   # S=0
    ([[0.4, 0.7, 0.2], [0.9, 0.3, 0.5]], [1, 2, 1], [[0.0, 0.0, 0.0]]),   # S=D
    ([[0.4, 0.7, 0.2], [0.9, 0.3, 0.5]], [2, 1, 2], [[1.0, 0.0, 0.5]]),   # 0<S<D
    ([[0.4, 0.4, 0.2], [0.9, 0.9, 0.5]], [2, 3, 1], [[1.0, 1.0, 0.0]]),   # identical
    ([[0.4], [0.9]], [6], [[1.5]]),                                        # D=1
    ([[0.4, 0.7], [0.9, 0.3]], [2, 2], [[1.0, 0.0], [0.5, 2.0]]),          # two IS
    ([[0.4, 0.7], [0.9, 0.3]], [3, 0], [[1.0, 2.0]]),                      # empty class
])
def test_fixed_rate_reduces_to_mvac(L, N, Z):
    L = np.array(L, float)
    N = np.array(N, int)
    Z = np.array(Z, float)
    mu = np.ones((L.shape[0], max(int(N.sum()), 1)))
    X, Q, U, C, pij = pfqn_mvacld(L, N, Z, mu)
    Xr, Qr, Ur, Cr = pfqn_mvac(L, N, Z)
    np.testing.assert_allclose(X, Xr, atol=EXACT_TOL)
    np.testing.assert_allclose(Q, Qr, atol=EXACT_TOL)
    np.testing.assert_allclose(pij.sum(axis=1), 1.0, atol=EXACT_TOL)
    # C is the LD-family (1xR) cycle time, NOT pfqn_mvac's (MxR) residence time
    assert C.shape == (L.shape[1],)


def test_default_mu_is_fixed_rate():
    L = np.array([[0.4, 0.7], [0.9, 0.3]])
    N = np.array([2, 2])
    Z = np.array([[1.0, 0.5]])
    X, Q, U, C, pij = pfqn_mvacld(L, N, Z)
    Xr, Qr, Ur, Cr = pfqn_mvac(L, N, Z)
    np.testing.assert_allclose(X, Xr, atol=EXACT_TOL)
    np.testing.assert_allclose(Q, Qr, atol=EXACT_TOL)


# ---- oracle 2: genuine load dependence vs pfqn_mvald ----

def test_multiserver_two_servers():
    _assert_matches_mvald([[0.4, 0.7], [0.9, 0.3]], [2, 2], [[1.0, 0.5]],
                          _ms_rates(2, 4, 2))


def test_multiserver_no_delay():
    # S=D: every chain resolved by a label interchange, with load dependence
    _assert_matches_mvald([[0.4, 0.7], [0.9, 0.3]], [2, 3], [[0.0, 0.0]],
                          _ms_rates(2, 5, 3))


def test_delay_emulated_by_rate_n():
    # mu(j,n)=n makes a QLD center an infinite server; c_i = 1 falls out of (21)
    _assert_matches_mvald([[0.4, 0.7], [0.9, 0.3]], [2, 2], [[0.0, 0.0]],
                          np.tile(np.arange(1, 5, dtype=float), (2, 1)))


def test_nonmonotone_rates():
    # nothing in (21)-(25) requires mu to be monotone or concave
    Nt = 4
    mu = np.tile(1.0 + 0.5 * np.sin(np.arange(1, Nt + 1)), (2, 1))
    _assert_matches_mvald([[0.5, 0.2], [0.3, 0.8]], [2, 2], [[1.0, 1.0]], mu)


def test_single_class_multiserver():
    _assert_matches_mvald([[0.4], [0.9]], [5], [[1.5]], _ms_rates(2, 5, 2))


def test_three_classes_multiserver():
    _assert_matches_mvald([[0.3, 0.5, 0.7], [0.6, 0.1, 0.4]], [2, 1, 2],
                          [[1.0, 0.0, 0.5]], _ms_rates(2, 5, 2))


def test_per_station_rates_differ():
    # a fixed-rate queue, a 2-server queue and a delay in one network
    Nt = 4
    mu = np.vstack([np.ones(Nt), np.minimum(np.arange(1, Nt + 1), 2)]).astype(float)
    _assert_matches_mvald([[0.4, 0.7], [0.9, 0.3]], [2, 2], [[1.0, 0.5]], mu)


def test_randomized_against_mvald():
    rng = np.random.default_rng(3)
    checked = 0
    for _ in range(30):
        M = int(rng.integers(1, 3))
        R = int(rng.integers(1, 3))
        L = np.round(10 * rng.random((M, R))) / 5
        N = rng.integers(0, 3, R)
        Z = [np.zeros((1, R)), np.round(10 * rng.random((1, R))) / 5][int(rng.integers(0, 2))]
        if N.sum() == 0:
            continue
        if not np.all(np.any(np.vstack([L, Z]) > 0, axis=0) | (N == 0)):
            continue
        Nt = int(N.sum())
        mu = np.minimum(np.tile(np.arange(1, Nt + 1), (M, 1)),
                        rng.integers(1, 4, (M, 1))).astype(float)
        _assert_matches_mvald(L, N, Z, mu)
        checked += 1
    assert checked > 10, 'too few feasible random cases were generated'


def test_mvald_empty_class_has_no_cycle_time():
    """
    Regression: pfqn_mvald left an empty class's cycle time as an unguarded
    N/XN = 0/0. MATLAB/JAR returned NaN there (Python already guarded), which
    propagated to callers such as pfqn_mvams. An empty class has zero
    throughput, hence no cycle time: all three codebases now report 0, as
    pfqn_dac and pfqn_mvacld do.
    """
    L = np.array([[0.4, 0.7], [0.9, 0.3]])
    N = np.array([3, 0])
    Z = np.array([1.0, 2.0])
    mu = _ms_rates(2, 3, 2)
    ref = pfqn_mvald(L, N, Z, mu)
    C = np.asarray(ref[3]).flatten()
    assert not np.any(np.isnan(C)), 'an empty class must not yield a NaN cycle time'
    assert C[1] == pytest.approx(0.0, abs=EXACT_TOL)
    assert C[0] > 0
    # and it must agree with pfqn_mvacld on every class, empty ones included
    _, _, _, Cld, _ = pfqn_mvacld(L, N, np.atleast_2d(Z), mu)
    np.testing.assert_allclose(Cld, C, atol=EXACT_TOL)


# ---- input validation ----

def test_rejects_bad_rate_matrix():
    L = np.array([[0.4, 0.7], [0.9, 0.3]])
    N = np.array([2, 2])
    with pytest.raises(ValueError):
        pfqn_mvacld(L, N, None, np.ones((3, 4)))  # wrong number of centers
    with pytest.raises(ValueError):
        pfqn_mvacld(L, N, None, np.ones((2, 2)))  # too few populations
    with pytest.raises(ValueError):
        pfqn_mvacld(L, N, None, np.zeros((2, 4)))  # nonpositive rates
