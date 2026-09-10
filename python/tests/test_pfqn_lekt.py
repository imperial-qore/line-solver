"""pfqn_lekt computes the estimator pfqn_ble and pfqn_bkt share, on the cheaper side.

With a think time the two corrected expansions are one estimator (dual saddle points,
Sylvester on the Hessians), so the value must not depend on the route. Without one the
LE side carries M kappa - r(N+M) in place of ble's (M-1) kappa, which is exactly what
lands it on the KT value.
"""
import numpy as np
import pytest
from scipy.special import gammaln

from line_solver.api.pfqn.asymptotic import pfqn_le, pfqn_ble, pfqn_lekt, pfqn_lekt_route
from line_solver.api.pfqn.kt import pfqn_bkt

KAPPA = 1.0 - np.log(2 * np.pi) / 2
MODEL_A = np.array([[1.0, 0.5], [0.7, 1.2], [0.3, 0.9]])          # 3 x 2: the KT side
WIDE = MODEL_A.T.copy()                                            # 2 x 3: the LE side


def _r(a):
    return gammaln(a) - (a - 0.5) * np.log(a) + a - 0.5 * np.log(2 * np.pi)


def test_the_route_follows_the_smaller_dimension_and_self_loops():
    assert pfqn_lekt_route(MODEL_A, [2, 3], [1, 0.5]) == 'kt'
    assert pfqn_lekt_route(WIDE, [2, 3, 4], [1, 0.5, 0.2]) == 'le'
    assert pfqn_lekt_route(WIDE, [2, 3, 4]) == 'le'
    Ls = np.array([[1.0, 0.0, 0.4], [0.7, 1.2, 0.0]])           # class 2 self-loops at Z = 0
    assert pfqn_lekt_route(Ls, [3, 2, 2]) == 'kt'
    assert pfqn_lekt_route(Ls, [3, 2, 2], [0, 1, 1]) == 'le'      # ... but not once both single-station classes think


def test_the_kt_side_is_bkt():
    N, Z = np.array([2.0, 3.0]), np.array([1.0, 0.5])
    assert pfqn_lekt(MODEL_A, N, Z)[1] == pfqn_bkt(MODEL_A, N, Z)[1]
    assert pfqn_lekt(MODEL_A, N)[1] == pfqn_bkt(MODEL_A, N)[1]


def test_the_le_side_with_a_think_time_is_ble_and_lands_on_bkt():
    N, Z = np.array([2.0, 3.0, 4.0]), np.array([1.0, 0.5, 0.2])
    v = pfqn_lekt(WIDE, N, Z)[1]
    assert v == pfqn_ble(WIDE, N, Z)[1]
    assert abs(v - pfqn_bkt(WIDE, N, Z)[1]) < 1e-5              # solver floors, not the identity


def test_the_le_side_without_a_think_time_carries_the_common_constant():
    N = np.array([2.0, 3.0, 4.0])
    M = WIDE.shape[0]
    eta = N.sum() + M
    v = pfqn_lekt(WIDE, N)[1]
    assert v == pytest.approx(pfqn_le(WIDE, N)[1] + M * KAPPA - _r(eta), abs=1e-12)
    assert v == pytest.approx(pfqn_ble(WIDE, N)[1] + KAPPA - _r(eta), abs=1e-12)
    assert abs(v - pfqn_bkt(WIDE, N)[1]) < 1e-5                 # Proposition: bkt - LE = M kappa - r(eta)


def test_the_nc_dispatcher_reaches_it():
    from line_solver.api.pfqn.nc import pfqn_nc
    N, Z = np.array([2.0, 3.0]), np.array([1.0, 0.5])
    assert pfqn_nc(MODEL_A, N, Z, method='lekt')[1] == pytest.approx(pfqn_lekt(MODEL_A, N, Z)[1], abs=1e-12)


def test_matches_matlab():
    # MATLAB pfqn_lekt on the transposed modelA (2 x 3, the LE side) and on a self-looping model
    N, Z = np.array([2.0, 3.0, 4.0]), np.array([1.0, 0.5, 0.2])
    assert pfqn_lekt(WIDE, N, Z)[1] == pytest.approx(7.64923112714508, abs=1e-8)
    assert pfqn_lekt(WIDE, N)[1] == pytest.approx(7.09110182626673, abs=1e-8)
    Ls = np.array([[1.0, 0.0, 0.4], [0.7, 1.2, 0.0]])
    assert pfqn_lekt(Ls, np.array([3.0, 2.0, 2.0]))[1] == pytest.approx(2.08635455377189, abs=1e-8)


def test_degenerate_inputs_pass_through():
    assert pfqn_lekt(np.zeros((0, 0)), np.array([]))[1] == 0.0
    assert pfqn_lekt(WIDE, np.array([0.0, 0.0, 0.0]))[1] == pfqn_ble(WIDE, np.array([0.0, 0.0, 0.0]))[1]
