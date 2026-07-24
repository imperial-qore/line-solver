"""
Validation of the MVAC (mean value analysis by chain) algorithm of Conway, de Souza e
Silva and Lavenberg, IEEE Trans. Computers 38(3):432-442, 1989.

MVAC is exact, so it is checked against two independent references: the worked example of
Section III of the paper, whose closed-form values are quoted there, and the classic
population-recursion MVA of pfqn_mva, which computes the same measures by a structurally
unrelated recursion. The three chain-resolution paths of the algorithm (part 1 for chain
K, part 2 for the chains visiting an IS center, and the label interchanges for the chains
visiting only SSFR centers) are covered explicitly, since a network exercising only one of
them would leave the others untested.
"""

import numpy as np
import pytest

from line_solver.api.pfqn import pfqn_mvac, pfqn_mva

EXACT_TOL = 1e-10


def _mva_ref(L, N, Z):
    """Throughput and queue-length from pfqn_mva, which returns (X, C, Q, U, ...)."""
    ret = pfqn_mva(L, N, Z.sum(axis=0))
    return np.asarray(ret[0]).flatten(), np.asarray(ret[2])


def _assert_agrees_with_mva(L, N, Z):
    L = np.array(L, dtype=float)
    N = np.array(N, dtype=int)
    Z = np.atleast_2d(np.array(Z, dtype=float))
    X, Q, U, C = pfqn_mvac(L, N, Z)
    Xref, Qref = _mva_ref(L, N, Z)
    np.testing.assert_allclose(X, Xref, atol=EXACT_TOL)
    np.testing.assert_allclose(Q, Qref, atol=EXACT_TOL)
    np.testing.assert_allclose(U, X * L, atol=EXACT_TOL)
    # the customers must all be accounted for
    total = Q.sum() + (X * Z.sum(axis=0)).sum()
    np.testing.assert_allclose(total, N.sum(), atol=1e-8)


def test_paper_section_iii_example():
    # J = 2 SSFR centers, K = 2 single-customer chains, a11=1, a21=2, a12=2, a22=3.
    # The paper reports lambda_2 = 3/23, L_12 = 8/23, L_22 = 15/23, L_1 = 15/23,
    # L_2 = 31/23.
    X, Q, U, C = pfqn_mvac(np.array([[1.0, 2.0], [2.0, 3.0]]), np.array([1, 1]))
    assert X[1] == pytest.approx(3.0 / 23.0, abs=EXACT_TOL)
    assert Q[0, 1] == pytest.approx(8.0 / 23.0, abs=EXACT_TOL)
    assert Q[1, 1] == pytest.approx(15.0 / 23.0, abs=EXACT_TOL)
    # total mean number of customers at each center, chains 1 and 2 summed
    assert Q[0, :].sum() == pytest.approx(15.0 / 23.0, abs=EXACT_TOL)
    assert Q[1, :].sum() == pytest.approx(31.0 / 23.0, abs=EXACT_TOL)


def test_all_chains_visit_delay():
    # S = 0: all chains but K are resolved by part 2
    _assert_agrees_with_mva([[0.4, 0.7, 0.2], [0.9, 0.3, 0.5]], [1, 2, 1], [[1.0, 0.5, 2.0]])


def test_no_delay():
    # S = D: every chain is resolved by a label interchange
    _assert_agrees_with_mva([[0.4, 0.7, 0.2], [0.9, 0.3, 0.5]], [1, 2, 1], [[0.0, 0.0, 0.0]])


def test_mixed_delay_and_queue_only_chains():
    # 0 < S < D: part 2 and the label interchanges are both used
    _assert_agrees_with_mva([[0.4, 0.7, 0.2], [0.9, 0.3, 0.5]], [2, 1, 2], [[1.0, 0.0, 0.5]])


def test_identical_classes_collapse():
    _assert_agrees_with_mva([[0.4, 0.4, 0.2], [0.9, 0.9, 0.5]], [2, 3, 1], [[1.0, 1.0, 0.0]])


def test_single_class():
    _assert_agrees_with_mva([[0.4], [0.9]], [6], [[1.5]])


def test_delay_only_network():
    _assert_agrees_with_mva([[0.0, 0.0], [0.0, 0.0]], [2, 2], [[1.0, 2.0]])


def test_multiple_delay_centers():
    _assert_agrees_with_mva([[0.4, 0.7], [0.9, 0.3]], [2, 2], [[1.0, 0.0], [0.5, 2.0]])


def test_empty_class():
    X, Q, U, C = pfqn_mvac(np.array([[0.4, 0.7], [0.9, 0.3]]), np.array([3, 0]),
                           np.array([[1.0, 2.0]]))
    assert X[1] == pytest.approx(0.0, abs=EXACT_TOL)
    np.testing.assert_allclose(Q[:, 1], 0.0, atol=EXACT_TOL)
    _assert_agrees_with_mva([[0.4, 0.7], [0.9, 0.3]], [3, 0], [[1.0, 2.0]])


def test_many_chains_few_centers():
    # the regime MVAC targets: few centers, many chains
    _assert_agrees_with_mva([[0.3, 0.5, 0.7, 0.2, 0.9], [0.6, 0.1, 0.4, 0.8, 0.3]],
                            [2, 2, 2, 2, 2], [[1.0, 0.0, 2.0, 0.0, 0.5]])


def test_randomized_against_mva():
    rng = np.random.default_rng(7)
    checked = 0
    for _ in range(60):
        M = int(rng.integers(1, 4))
        R = int(rng.integers(1, 4))
        L = np.round(10 * rng.random((M, R))) / 5
        N = rng.integers(0, 4, R)
        Z = [np.zeros((1, R)),
             np.round(10 * rng.random((1, R))) / 5,
             np.round(10 * rng.random((2, R))) / 5][int(rng.integers(0, 3))]
        if N.sum() == 0:
            continue
        # every populated class must have a nonzero demand somewhere
        if not np.all(np.any(np.vstack([L, Z]) > 0, axis=0) | (N == 0)):
            continue
        _assert_agrees_with_mva(L, N, Z)
        checked += 1
    assert checked > 20, 'too few feasible random cases were generated'


def test_rejects_class_count_mismatch():
    with pytest.raises(ValueError):
        pfqn_mvac(np.array([[0.4, 0.7]]), np.array([1, 1, 1]))
    with pytest.raises(ValueError):
        pfqn_mvac(np.array([[0.4, 0.7]]), np.array([1, 1]), np.array([[1.0]]))


def test_rejects_zero_demand_network():
    with pytest.raises(ValueError):
        pfqn_mvac(np.zeros((2, 2)), np.array([1, 1]))
