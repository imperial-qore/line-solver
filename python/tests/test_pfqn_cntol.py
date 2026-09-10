"""Chandy-Neuse population-scaled termination test (pfqn_cntol).

The cutoff 1/(4000 + 16*sum(N)) and the metric max_{i,r}|dQ(i,r)|/N_r it is compared
against are published in K. M. Chandy, D. Neuse, Commun. ACM 25(2):126-134, 1982, p.129
and appendix; LQNS runs the same expression from SchweitzerCommon. It is opt-in here:
the sentinel is tol='cn' (or NaN) and the default tolerances are untouched, so the
figures below double as a guard that the sentinel path is the only thing it changes.

The asserted values are MATLAB's, printed by pfqn_linearizer/pfqn_bs with tol = 'cn'
on the same model.
"""

import math

import numpy as np
import pytest

from line_solver.api.pfqn import is_cntol, pfqn_bs, pfqn_cntol, pfqn_linearizer

L = np.array([[1.0, 2.0], [3.0, 1.0], [0.5, 0.5]])
N = np.array([8.0, 1.0])
Z = np.array([0.0, 0.0])
TYPE = ['PS'] * 3


def test_cutoff_matches_the_published_expression():
    assert pfqn_cntol([10]) == pytest.approx(1.0 / (4000.0 + 160.0), abs=1e-15)
    assert pfqn_cntol(N) == pytest.approx(1.0 / (4000.0 + 16.0 * 9.0), abs=1e-15)
    # Below 0.00025 even at a population of one, as the paper's appendix states.
    assert pfqn_cntol([1]) < 0.00025
    # Decreasing in the population, which is the point of the scaling.
    assert pfqn_cntol([100]) < pfqn_cntol([10])


def test_sentinel_recognition():
    assert is_cntol('cn')
    assert is_cntol('CN')
    assert is_cntol(float('nan'))
    assert not is_cntol(1e-8)
    with pytest.raises(ValueError):
        is_cntol('loose')


def test_linearizer_under_the_chandy_neuse_test_matches_matlab():
    Q, U, W, T, C, X, it = pfqn_linearizer(L, N, Z, TYPE, 'cn', 1000)
    X = np.asarray(X).ravel()
    assert X[0] == pytest.approx(0.304718171165, abs=1e-9)
    assert X[1] == pytest.approx(0.084038879640, abs=1e-9)
    assert Q[0, 0] == pytest.approx(0.555262542605, abs=1e-9)
    assert Q[1, 0] == pytest.approx(7.255702544745, abs=1e-9)
    assert Q[2, 0] == pytest.approx(0.189034912650, abs=1e-9)
    assert Q[0, 1] == pytest.approx(0.251901005767, abs=1e-9)
    assert Q[1, 1] == pytest.approx(0.697720815306, abs=1e-9)
    assert Q[2, 1] == pytest.approx(0.050378178927, abs=1e-9)

    # NaN is the cross-language spelling of the same sentinel.
    Xnan = np.asarray(pfqn_linearizer(L, N, Z, TYPE, float('nan'), 1000)[5]).ravel()
    assert Xnan[0] == pytest.approx(X[0], abs=1e-15)

    # The looser published cutoff must stop earlier than the 1e-8 default, and the two
    # answers must differ: an identical result would mean the sentinel did nothing.
    Xdef = np.asarray(pfqn_linearizer(L, N, Z, TYPE, 1e-8, 1000)[5]).ravel()
    itdef = pfqn_linearizer(L, N, Z, TYPE, 1e-8, 1000)[6]
    assert it < itdef
    assert abs(X[0] - Xdef[0]) > 1e-12
    assert X[0] == pytest.approx(Xdef[0], abs=1e-4)


def test_bard_schweitzer_under_the_chandy_neuse_test_matches_matlab():
    XN, QN, UN, RN, it = pfqn_bs(L, N, Z, 'cn', 1000)
    XN = np.asarray(XN).ravel()
    assert XN[0] == pytest.approx(0.301066244185, abs=1e-9)
    assert XN[1] == pytest.approx(0.083878140075, abs=1e-9)

    XNdef, _, _, _, itdef = pfqn_bs(L, N, Z, 1e-6, 1000)
    XNdef = np.asarray(XNdef).ravel()
    assert it < itdef
    assert XN[0] == pytest.approx(XNdef[0], abs=1e-4)


def test_empty_class_does_not_divide_by_zero():
    # An empty class has N_r = 0, so the published metric would divide by zero there;
    # it must be skipped rather than guarded downstream.
    Nz = np.array([8.0, 0.0])
    XN, QN, _, _, _ = pfqn_bs(L, Nz, Z, 'cn', 1000)
    assert np.all(np.isfinite(np.asarray(XN)))
    assert np.all(np.isfinite(QN))
    assert math.isclose(float(np.asarray(XN).ravel()[1]), 0.0, abs_tol=1e-15)
