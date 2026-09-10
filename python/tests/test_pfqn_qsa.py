"""Queue-Shift Approximation (pfqn_qsa).

The three stress networks of Sect. 5.2 of Schweitzer, Serazzi and Broglia
(Tools'98, LNCS 1469) are fully specified in print together with the QSA errors
they produce under the paper's own metric eq. (17), so they pin the algorithm
rather than this port: any drift in the shift definition or in the
extrapolation (15) moves them.
"""

import numpy as np
import pytest

from line_solver.api.pfqn import pfqn_mva, pfqn_qsa, pfqn_linearizer
from line_solver.constants import SchedStrategy

TABLE2 = np.array([[16., 50], [19, 28], [0, 42], [12, 25], [39, 29]])
TABLE3 = np.array([[16., 73], [59, 15], [6, 26], [2, 36], [10, 0]])
TABLE4 = np.array([[10., 1, 1], [1, 10, 1], [1, 1, 10]])


def eq17(L, N, Q, U):
    """err(Q) = max |Q_appr - Q_MVA| / K_r, err(U) = max |U_appr - U_MVA|."""
    M = L.shape[0]
    _, _, Qe, Ue, _, _, _ = pfqn_mva(L, N, np.zeros(len(N)))
    eQ = np.max(np.abs(Q - Qe) / np.tile(N, (M, 1)))
    eU = np.max(np.abs(U - Ue))
    return eQ, eU


def test_table2_stress_case():
    # K = (21,29), the case with the largest Linearizer error.
    # Published: QSA err(Q) = 0.001812, err(U) = 0.000639.
    N = np.array([21., 29.])
    Q, U = pfqn_qsa(TABLE2, N, np.zeros(2))[:2]
    eQ, eU = eq17(TABLE2, N, Q, U)
    assert eQ == pytest.approx(0.001812, rel=1e-3)
    assert eU == pytest.approx(0.000639, rel=1e-3)


def test_table3_stress_case():
    # K = (111,89), maximizing the Linearizer-to-QSA error ratio. The published
    # 1e-10 is that paper's own residual tolerance; the Newton solve here
    # reaches machine precision, consistent with its claim of "errors of zero
    # (at least 6 digits)".
    N = np.array([111., 89.])
    Q, U = pfqn_qsa(TABLE3, N, np.zeros(2))[:2]
    eQ, eU = eq17(TABLE3, N, Q, U)
    assert eQ < 4.18e-10
    assert eU < 5.54e-10


def test_table4_three_class_stress_case():
    # K = (2,2,2) on the Chandy-Neuse Example 2 loadings. Published:
    # err(Q) = 0.000185, err(U) = 0.000415. The tiny population is what makes
    # this discriminating: the extrapolation (15) has no 1/Ksum to hide in.
    N = np.array([2., 2., 2.])
    Q, U = pfqn_qsa(TABLE4, N, np.zeros(3))[:2]
    eQ, eU = eq17(TABLE4, N, Q, U)
    assert eQ == pytest.approx(0.000185, rel=1e-2)
    assert eU == pytest.approx(0.000415, rel=1e-2)


def test_beats_linearizer_on_table4():
    # Sect. 5.2: on this model QSA has roughly half the Linearizer error.
    N = np.array([2., 2., 2.])
    Q, U = pfqn_qsa(TABLE4, N, np.zeros(3))[:2]
    eQ, eU = eq17(TABLE4, N, Q, U)
    assert eQ < 0.000318   # the published Linearizer error on the same model
    assert eU < 0.000716


def test_population_is_conserved_and_queues_stay_non_negative():
    N = np.array([21., 29.])
    Q = pfqn_qsa(TABLE2, N, np.zeros(2))[0]
    assert Q.min() >= -1e-12
    assert Q.sum() == pytest.approx(50.0, abs=1e-10)
    for r in range(2):
        assert Q[:, r].sum() == pytest.approx(N[r], abs=1e-10)   # eq. (6)


def test_think_time_is_held_outside_the_queues():
    L = np.array([[1.0, 0.5], [0.7, 1.2], [0.3, 0.9]])
    N = np.array([6., 4.])
    Z = np.array([2., 1.])
    Q, U, W, C, X, _ = pfqn_qsa(L, N, Z)
    assert Q.sum() + float(X @ Z) == pytest.approx(10.0, abs=1e-10)
    _, _, Qe, _, _, _, _ = pfqn_mva(L, N, Z)
    assert np.abs(Q - Qe).max() < 2e-2


def test_delay_centre_declared_through_type_has_no_queueing_term():
    L = np.array([[1.0, 0.5], [0.7, 1.2], [0.3, 0.9]])
    N = np.array([6., 4.])
    type_ = [SchedStrategy.PS, SchedStrategy.INF, SchedStrategy.PS]
    Q, U, W, C, X, _ = pfqn_qsa(L, N, np.zeros(2), type_)
    for r in range(2):
        assert Q[1, r] == pytest.approx(X[r] * L[1, r], rel=1e-12)  # eq. (13b)
        assert W[1, r] == pytest.approx(L[1, r], rel=1e-12)


def test_two_level_variant_is_the_cruder_one():
    # Sect. 3.1 gives eq. (14) an error of O(1/Ksum); on Table 2 (Ksum = 50)
    # that is the order observed, and it must not be mistaken for the default.
    N = np.array([21., 29.])
    Q2 = pfqn_qsa(TABLE2, N, np.zeros(2), None, 1e-10, 100, 2)[0]
    Q3 = pfqn_qsa(TABLE2, N, np.zeros(2), None, 1e-10, 100, 3)[0]
    _, _, Qe, _, _, _, _ = pfqn_mva(TABLE2, N, np.zeros(2))
    e2 = np.max(np.abs(Q2 - Qe) / np.tile(N, (5, 1)))
    e3 = np.max(np.abs(Q3 - Qe) / np.tile(N, (5, 1)))
    assert e3 < e2
    assert 1.0 / 50.0 < e2 < 0.1


def test_exact_when_one_job_per_class():
    # Every arrival sees the other classes only, so the arrival-instant queue
    # length is exact and so is the shift.
    L = np.array([[1.0, 2.0], [3.0, 1.0], [2.0, 2.0]])
    N = np.array([1., 1.])
    Q = pfqn_qsa(L, N, np.zeros(2))[0]
    _, _, Qe, _, _, _, _ = pfqn_mva(L, N, np.zeros(2))
    assert np.abs(Q - Qe).max() < 1e-10


def test_empty_class_contributes_nothing():
    L = np.array([[1.0, 2.0], [3.0, 1.0], [2.0, 2.0]])
    N = np.array([3., 0.])
    Q, U, W, C, X, _ = pfqn_qsa(L, N, np.zeros(2))
    assert np.all(Q[:, 1] == 0.0)
    assert X[1] == 0.0
    assert Q[:, 0].sum() == pytest.approx(3.0, abs=1e-10)


def test_single_class_matches_exact_mva_away_from_saturation():
    L = np.array([[1.], [2.], [3.]])
    N = np.array([200.])
    Q = pfqn_qsa(L, N, np.zeros(1))[0]
    _, _, Qe, _, _, _, _ = pfqn_mva(L, N, np.zeros(1))
    Ql = pfqn_linearizer(L, N, np.zeros(1), [SchedStrategy.PS] * 3, 1e-12, 1000)[0]
    assert np.abs(Q - Qe).max() < 1e-9
    assert np.abs(Q - Qe).max() < np.abs(Ql - Qe).max()
