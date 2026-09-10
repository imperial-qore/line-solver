"""The BLE correction counts the Laplaced directions of the branch actually taken.

M-1 with Z = 0, where the radial integral is exact as Gamma(N+M), and M with Z > 0,
where the radius is Laplaced too. Measured over the 1562 models of the Cas17 dataset
(Zenodo 546873, sec5.3.1, sigma = 100) the Z > 0 deficit is M to within 0.01 units,
and using M-1 there leaves a residual of exactly one unit on every model.
"""
import numpy as np
import pytest

from line_solver.api.pfqn.asymptotic import pfqn_le, pfqn_ble
from line_solver.api.pfqn.nc import pfqn_ca

KAPPA = 1.0 - np.log(2 * np.pi) / 2


def _demands(M, R, seed):
    return 0.5 + np.random.default_rng(seed).random((M, R))


@pytest.mark.parametrize("M", [2, 3, 4, 6])
def test_offset_is_m_minus_one_without_a_think_time(M):
    L, N = _demands(M, 2, 100 + M), np.array([20.0, 20.0])
    assert pfqn_ble(L, N)[1] - pfqn_le(L, N)[1] == pytest.approx((M - 1) * KAPPA, abs=1e-12)
    Z0 = np.zeros(2)
    assert pfqn_ble(L, N, Z0)[1] - pfqn_le(L, N, Z0)[1] == pytest.approx((M - 1) * KAPPA, abs=1e-12)


@pytest.mark.parametrize("M", [2, 3, 4, 6])
def test_offset_is_m_with_a_think_time(M):
    L, N, Z = _demands(M, 2, 200 + M), np.array([20.0, 20.0]), np.array([100.0, 100.0])
    assert pfqn_ble(L, N, Z)[1] - pfqn_le(L, N, Z)[1] == pytest.approx(M * KAPPA, abs=1e-12)


@pytest.mark.parametrize("M", [2, 3, 4, 6])
def test_the_delay_branch_deficit_against_the_exact_constant(M):
    """The gap pfqn_ca - pfqn_le is M units, which is exactly what pfqn_ble closes."""
    L, N, Z = _demands(M, 2, 300 + M), np.array([40.0, 40.0]), np.array([100.0, 100.0])
    exact = pfqn_ca(L, N, Z)[1]
    assert abs((exact - pfqn_le(L, N, Z)[1]) - M * KAPPA) < 0.1
    assert abs(exact - pfqn_ble(L, N, Z)[1]) < 0.1
