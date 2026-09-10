"""BKT subtracts the exact Stirling remainder of each class direction pfqn_kt Laplaces.

The expansion extracts [u^N] from the generating function by steepest descent and so
carries Stirling's approximation of log(N!) in place of log(N!), one remainder
s(N) = log(N!) - (N log N - N + log(2 pi N)/2) per class. With a think time the
corrected expansion is the SAME estimator as pfqn_ble (LE corrected by M units):
the two saddle points are one point in dual coordinates and Sylvester's identity
exchanges the R x R Hessian determinant for the M x M one. Without a think time
they differ by the single constant kappa - r(N+M), the remainder of the radial
Gamma(N+M) direction that LE integrates exactly and KT does not.
"""
import numpy as np
import pytest
from scipy.special import gammaln

from line_solver.api.pfqn.kt import pfqn_kt, pfqn_bkt, stirling_remainder
from line_solver.api.pfqn.asymptotic import pfqn_le, pfqn_ble
from line_solver.api.pfqn.nc import pfqn_ca

KAPPA = 1.0 - np.log(2 * np.pi) / 2


def _r(a):
    """Stirling remainder of a Gamma(a) direction, r(1) = kappa."""
    return gammaln(a) - (a - 0.5) * np.log(a) + a - 0.5 * np.log(2 * np.pi)


def _demands(M, R, seed):
    return 0.5 + np.random.default_rng(seed).random((M, R))


def test_the_remainder_at_one_is_the_ble_constant():
    assert stirling_remainder(1.0) == pytest.approx(KAPPA, abs=1e-15)
    assert stirling_remainder(200.0) == pytest.approx(1 / (12 * 200.0), rel=2e-4)


@pytest.mark.parametrize("M,R", [(2, 1), (3, 2), (5, 3), (4, 4)])
def test_the_correction_is_the_sum_of_the_class_remainders(M, R):
    L, N = _demands(M, R, 10 * M + R), np.arange(2.0, 2.0 + R)
    for Z in (None, np.zeros(R), np.full(R, 3.0)):
        gap = pfqn_bkt(L, N, Z)[1] - pfqn_kt(L, N, Z)[1]
        assert gap == pytest.approx(-np.sum(stirling_remainder(N)), abs=1e-12)


@pytest.mark.parametrize("M,R", [(3, 2), (4, 3), (6, 2)])
def test_with_a_think_time_bkt_is_ble(M, R):
    """Proposition: for Z > 0 the two corrected expansions coincide; ~1e-7 is the
    accuracy of the two saddle-point solvers, not of the identity."""
    L, N, Z = _demands(M, R, 300 + M), np.full(R, 6.0), np.full(R, 2.5)
    assert abs(pfqn_bkt(L, N, Z)[1] - pfqn_ble(L, N, Z)[1]) < 1e-5


@pytest.mark.parametrize("M,R", [(3, 2), (4, 3), (6, 2)])
def test_without_a_think_time_the_gap_to_le_is_a_known_constant(M, R):
    L, N = _demands(M, R, 400 + M), np.full(R, 6.0)
    eta = N.sum() + M
    gap = pfqn_bkt(L, N)[1] - pfqn_le(L, N)[1]
    assert gap == pytest.approx(M * KAPPA - _r(eta), abs=1e-5)


def test_light_load_against_exact_convolution():
    """At Z = 100 the class directions are Poisson and the remainder is the whole error."""
    L, N, Z = _demands(4, 2, 7), np.array([8.0, 6.0]), np.array([100.0, 100.0])
    exact = pfqn_ca(L, N, Z)[1]
    err_kt = abs(pfqn_kt(L, N, Z)[1] - exact)
    err_pp = abs(pfqn_bkt(L, N, Z)[1] - exact)
    assert err_kt > 0.01                       # about sum_r 1/(12 N_r)
    assert err_pp < 1e-3 and err_pp < err_kt / 20


def test_empty_and_self_looping_classes_carry_no_remainder():
    # class 2 empty: only class 1 is Laplaced
    L, N = _demands(3, 2, 11), np.array([5.0, 0.0])
    assert pfqn_bkt(L, N)[1] - pfqn_kt(L, N)[1] == pytest.approx(-stirling_remainder(5.0), abs=1e-12)
    # class 2 visits one station with no think time: extracted exactly by pfqn_kt
    Ls = np.array([[1.0, 0.0], [0.7, 1.2], [0.3, 0.0]])
    N = np.array([3.0, 2.0])
    assert pfqn_bkt(Ls, N)[1] - pfqn_kt(Ls, N)[1] == pytest.approx(-stirling_remainder(3.0), abs=1e-12)
    # ... but with a think time it is Laplaced like any other
    Z = np.array([0.0, 1.0])
    assert pfqn_bkt(Ls, N, Z)[1] - pfqn_kt(Ls, N, Z)[1] == pytest.approx(-np.sum(stirling_remainder(N)), abs=1e-12)


def test_matches_matlab():
    L = np.array([[1.0, 0.5], [0.7, 1.2], [0.3, 0.9]])
    N, Z = np.array([2.0, 3.0]), np.array([1.0, 0.5])
    # MATLAB pfqn_bkt(L,[2 3],[1 .5]) and pfqn_bkt(L,[2 3])
    assert pfqn_bkt(L, N, Z)[1] == pytest.approx(4.77242756160185, abs=1e-8)
    assert pfqn_bkt(L, N)[1] == pytest.approx(4.13902507926115, abs=1e-8)
    for n, ref in ((10, 14.1321407982845), (20, 26.141862103174), (40, 50.201513735287)):
        assert pfqn_bkt(L, np.array([n, n], float), Z)[1] == pytest.approx(ref, abs=1e-6)
    Ls = np.array([[1.0, 0.0], [0.7, 1.2], [0.3, 0.0]])
    assert pfqn_bkt(Ls, np.array([3.0, 2.0]))[1] == pytest.approx(2.81538680478809, abs=1e-8)
    assert pfqn_bkt(L, np.array([3.0, 0.0]), Z)[1] == pytest.approx(1.97269392277537, abs=1e-8)


def test_the_nc_dispatcher_reaches_it():
    from line_solver.api.pfqn.nc import pfqn_nc
    L = np.array([[1.0, 0.5], [0.7, 1.2], [0.3, 0.9]])
    N, Z = np.array([2.0, 3.0]), np.array([1.0, 0.5])
    assert pfqn_nc(L, N, Z, method='bkt')[1] == pytest.approx(pfqn_bkt(L, N, Z)[1], abs=1e-12)


def test_degenerate_inputs_pass_through():
    assert pfqn_bkt(np.zeros((0, 0)), np.array([]))[1] == 0.0
    assert pfqn_bkt(_demands(2, 2, 1), np.array([0.0, 0.0]))[1] == 0.0
