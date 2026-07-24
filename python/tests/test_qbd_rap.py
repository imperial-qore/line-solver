"""
Validates qbd_rap, the general QBD with Rational Arrival Process components,
against the published benchmark of N. G. Bean and B. F. Nielsen,
"Quasi-Birth-and-Death Processes with Rational Arrival Process Components",
Stochastic Models, 26(3), 2010, pp. 309-334: the marginal level distribution of
Table 1, the closed form of R, the eigenvalues of R and the stability boundary
gamma = 1/2. Also covers the general (non rank-one) path for G by cross-checking
a MAP/MAP/1 queue against qbd_raprap1, and adds the PH-in-ME-clothing identity
for qbd_raprap1 itself.
"""

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from line_solver.api.mam import (map_lambda, map_mean, qbd_mapmap1,  # noqa: E402
                                 qbd_rap, qbd_raprap1)

# Marginal level distribution of Table 1 of Bean and Nielsen (2010), gamma = 0.25.
TABLE1 = [0.6736, 0.2175, 0.0726, 0.0242, 0.0081, 0.0027, 0.0009, 0.0003, 0.0001]

# Local block A1 = Ca = Cs of the example process.
EXAMPLE_A1 = np.array([[-1.0, 0.0, 0.0],
                       [-2.0 / 3, -1.0, 1.0],
                       [2.0 / 3, -1.0, -1.0]])

# Jump matrix Da of the arrival RAP of the example process.
EXAMPLE_DA = np.array([[14.0 / 5, -9.0 / 10, -9.0 / 10],
                       [26.0 / 15, -8.0 / 15, -8.0 / 15],
                       [58.0 / 15, -19.0 / 15, -19.0 / 15]])

# Jump matrix Ds of the service RAP, the rank-one product (1,2/3,4/3)'*(3,-1,-1).
EXAMPLE_DS = np.array([[1.0], [2.0 / 3], [4.0 / 3]]) @ np.array([[3.0, -1.0, -1.0]])


def example_blocks(g):
    return g * EXAMPLE_DA, EXAMPLE_A1, (1.0 - g) * EXAMPLE_DS


def solve_example(g, num_levels):
    A0, A1, A2 = example_blocks(g)
    return qbd_rap(A0, A1, A2, A0.copy(), g * A1, num_levels)


def example_R(g):
    """Closed form of R given in Section 5 of Bean and Nielsen (2010)."""
    return np.array([
        [-1.0 / 15 * g * (33 + 2 * g) / (g - 1), 0.0, 1.0 / 10 * g * (g + 9) / (g - 1)],
        [-2.0 / 45 * g * (31 + 4 * g) / (g - 1), 0.0, 2.0 / 15 * g * (g + 4) / (g - 1)],
        [-4.0 / 45 * g * (34 + g) / (g - 1), 0.0, 1.0 / 15 * g * (g + 19) / (g - 1)],
    ])


def test_table1_level_distribution():
    res = solve_example(0.25, 8)
    for n in range(9):
        assert abs(res.levelProb[n] - TABLE1[n]) < 5e-5, \
            'Table 1 of Bean and Nielsen (2010), level %d' % n
    # The reported levels 0..8 carry essentially all of the mass.
    assert abs(float(np.sum(res.levelProb)) - 1.0) < 1e-4


def test_R_closed_form_and_eigenvalues():
    g = 0.25
    res = solve_example(g, 8)
    assert np.linalg.norm(res.R - example_R(g)) < 1e-12
    # Eigenvalues of R are (0, -gamma/15, gamma/(1-gamma)).
    ev = np.sort(np.real(np.linalg.eigvals(res.R)))
    expected = np.sort(np.array([0.0, -g / 15.0, g / (1.0 - g)]))
    assert np.max(np.abs(ev - expected)) < 1e-12
    assert abs(res.spr - g / (1.0 - g)) < 1e-12


def test_G_is_rank_one_closed_form():
    g = 0.25
    A0, A1, A2 = example_blocks(g)
    res = solve_example(g, 4)
    G = res.G
    assert np.linalg.norm(A0 @ G @ G + A1 @ G + A2) < 1e-12
    # Rank-one A2 gives G = e*v/(v*e) with v = (3,-1,-1) and v*e = 1.
    expected_G = np.ones((3, 1)) @ np.array([[3.0, -1.0, -1.0]])
    assert np.max(np.abs(G - expected_G)) < 1e-12
    assert np.max(np.abs(res.U - (A1 + A0 @ G))) < 1e-12


def test_stability_boundary_at_gamma_one_half():
    # The queue is stable exactly when gamma < 1/2 (Corollary 8).
    assert solve_example(0.49, 4).QN > 0.0
    with pytest.raises(ValueError):
        solve_example(0.5, 4)
    with pytest.raises(ValueError):
        solve_example(0.55, 4)


def test_mean_queue_length_increases_with_gamma():
    # Figure 1 of the paper: the mean queue length grows monotonically to
    # infinity as gamma approaches 1/2 from below.
    prev = -1.0
    for g in (0.05, 0.15, 0.25, 0.35, 0.45, 0.49):
        qn = solve_example(g, 4).QN
        assert qn > prev
        prev = qn
    assert prev > 20.0


def test_general_G_path_against_mapmap1():
    # A MAP/MAP/1 queue whose service jump matrix has full rank: A2 is not rank
    # one, so the general functional-iteration plus Newton path for G is
    # exercised. The oracle is qbd_mapmap1, which reaches the same answer
    # through cyclic reduction and shares no code with qbd_rap. qbd_raprap1
    # would NOT be an independent oracle here: it now delegates to qbd_rap.
    C0 = np.array([[-0.5]])
    C1 = np.array([[0.5]])
    S0 = np.array([[-2.0, 0.1], [0.2, -3.0]])
    S1 = np.array([[1.9, 0.0], [0.0, 2.8]])
    assert np.linalg.matrix_rank(S1) == 2

    I2 = np.eye(2)
    A0 = np.kron(C1, I2)
    A1 = np.kron(C0, I2) + np.kron(np.eye(1), S0)
    A2 = np.kron(np.eye(1), S1)
    B1 = np.kron(C0, I2)
    res = qbd_rap(A0, A1, A2, A0.copy(), B1, 400)

    QNref, pqueue_ref = qbd_mapmap1((C0, C1), (S0, S1))[1], qbd_mapmap1((C0, C1), (S0, S1))[3]
    assert abs(res.QN - QNref) < 1e-8
    assert abs(res.levelProb[0] - float(np.sum(pqueue_ref[0, :]))) < 1e-9
    G = res.G
    assert np.linalg.norm(A0 @ G @ G + A1 @ G + A2) < 1e-10
    assert np.linalg.norm(G @ np.ones((2, 1)) - np.ones((2, 1))) < 1e-10


def test_non_conservative_blocks_are_rejected():
    A0, A1, A2 = example_blocks(0.25)
    with pytest.raises(ValueError):
        qbd_rap(A0, A1, A2 + np.eye(3))


def me_process(alpha, A):
    """Builds the ME process pair (A, (-A e) alpha) from a representation."""
    A = np.asarray(A, dtype=float)
    alpha = np.asarray(alpha, dtype=float)
    return A, np.outer(-A @ np.ones(A.shape[0]), alpha)


def exp_process(rate):
    return np.array([[-rate]]), np.array([[rate]])


def genuine_me():
    """Genuine (non phase-type) ME: poles -0.5 and -1 +- 2*pi*i."""
    w = 2 * np.pi
    alpha = [0.984694494294579, -0.040430911430916, 0.0557364171363366]
    A = [[-0.5, 0.0, 0.0], [0.0, -1.0, w], [0.0, -w, -1.0]]
    return me_process(alpha, A)


def test_raprap1_production_regression():
    """
    Production regression for qbd_raprap1 after its internals were replaced by
    the Theorem 7 construction of qbd_rap. These are the four models the MAM
    ME/RAP path is validated on; the values are the ones the cyclic-reduction
    implementation produced and must not move. They are the truncated level
    series, not the exact closed form, which is why qbd_raprap1 keeps its own
    truncation rule.
    """
    erlang_as_me = me_process([1.0, 0.0], [[-2.0, 2.0], [0.0, -2.0]])
    hyper_as_me = me_process([0.6, 0.4], [[-2.0, 0.0], [0.0, -0.5]])
    gen = genuine_me()

    for label, arr, svc, qn_expected, levels_expected in (
        ('a', exp_process(0.5), erlang_as_me, 0.874999998681, 27),
        ('b', exp_process(0.5), hyper_as_me, 1.522222217439, 55),
        ('c', exp_process(0.255775446238906), gen, 1.015214675095, 34),
        ('d', gen, exp_process(1.0 / 0.977420), 1.037045549083, 35),
    ):
        _, QN, _, pqueue = qbd_raprap1(arr, svc)[:4]
        assert abs(QN - qn_expected) < 1e-9, '(%s) QN %.12f vs %.12f' % (label, QN, qn_expected)
        assert pqueue.shape[0] == levels_expected, '(%s) truncated level count' % label


def test_raprap1_utilization_equals_rho():
    """
    Utilization of a single-server queue is exactly rho = lambda*E[S], whatever
    the correlation structure. This is a free exact oracle and it is what
    exposed the pre-rewrite defect: cyclic reduction silently returned a G with
    residual 1.19 and ||G*e-e|| = 0.39 for a multi-phase RAP arrival, giving
    UN = 0.4040 against the true 0.5000000733.
    """
    gen = genuine_me()
    svc = exp_process(1.0 / 0.977420)
    _, _, UN, _, _, _, G, B, L, F = qbd_raprap1(gen, svc)
    rho = map_lambda(*gen) * map_mean(*svc)
    assert abs(UN - rho) < 1e-12, 'utilization %.12f vs rho %.12f' % (UN, rho)
    assert np.linalg.norm(F @ G @ G + L @ G + B) < 1e-10, 'G must solve F*G^2 + L*G + B = 0'
    e = np.ones((G.shape[0], 1))
    assert np.linalg.norm(G @ e - e) < 1e-10, 'G*e must equal e'


def test_raprap1_matches_mapmap1_for_ph_in_me_clothing():
    """
    PH-in-ME-clothing identity for qbd_raprap1: an M/ME/1 queue whose
    matrix-exponential service is a similarity transform of an Erlang-2, and
    therefore not a nonnegative representation, must give the same answer as the
    M/Erlang/1 queue solved through the phase-type path of qbd_mapmap1.
    """
    lam = 0.4
    mu = 2.0
    C0 = np.array([[-lam]])
    C1 = np.array([[lam]])
    # Erlang-2 as a phase-type MAP.
    D0 = np.array([[-mu, mu], [0.0, -mu]])
    D1 = np.array([[0.0, 0.0], [mu, 0.0]])
    # Similarity transform with S*e = e: preserves the process but destroys
    # nonnegativity, so (H0,H1) is a genuine ME representation.
    S = np.array([[1.0, 0.0], [-0.5, 1.5]])
    Sinv = np.linalg.inv(S)
    H0 = Sinv @ D0 @ S
    H1 = Sinv @ D1 @ S
    off_diag_neg = any(H0[i, j] < 0.0 for i in range(2) for j in range(2) if i != j)
    assert off_diag_neg or np.any(H1 < 0.0), \
        'the transformed representation must have a negative entry, otherwise it is still a MAP'

    XNph, QNph, UNph = qbd_mapmap1((C0, C1), (D0, D1))[:3]
    XNme, QNme, UNme = qbd_raprap1((C0, C1), (H0, H1))[:3]
    assert abs(XNph - XNme) < 1e-6
    assert abs(UNph - UNme) < 1e-6
    assert abs(QNph - QNme) < 1e-6
