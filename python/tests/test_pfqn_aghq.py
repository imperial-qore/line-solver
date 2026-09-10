"""The adaptive Gauss-Hermite rule over the simplex factor of the McKenna-Mitra integral.

At Z = 0 the integrand is homogeneous, so the radius integrates exactly to Gamma(N+M) and
every bit of the error of pfqn_le sits in the closure of the remaining simplex integral.
pfqn_aghq reads the mode and curvature of that closure as a quadrature rule whose q = 1
member IS pfqn_le. With Z > 0 it integrates the radius instead of closing it.
"""
import numpy as np
import pytest

from line_solver.api.pfqn.asymptotic import pfqn_aghq, pfqn_le
from line_solver.api.pfqn.nc import pfqn_ca, pfqn_nc


def _demands(M, R, seed):
    return 0.5 + np.random.default_rng(seed).random((M, R))


@pytest.mark.parametrize("M", [2, 3, 4, 5])
def test_aghq_at_one_node_is_the_logistic_expansion(M):
    """A single node sits at the mode with weight sqrt(2 pi), which is pfqn_le itself."""
    L, N = _demands(M, 2, 10 + M), np.array([6.0, 4.0])
    assert pfqn_aghq(L, N, None, 1)[1] == pytest.approx(pfqn_le(L, N)[1], abs=1e-9)


@pytest.mark.parametrize("M", [2, 3, 4, 5])
def test_the_rule_is_exact_at_one_station(M):
    """M = 1 leaves a point mass on the simplex, so there is nothing left to close."""
    L, N = _demands(1, 2, 20 + M), np.array([4.0, 3.0])
    exact = pfqn_ca(L, N, np.zeros(2))[1]
    assert pfqn_aghq(L, N)[1] == pytest.approx(exact, abs=1e-9)


@pytest.mark.parametrize("q", [3, 5, 7])
def test_the_quadrature_converges_in_q(q):
    L, N = _demands(3, 2, 77), np.array([9.0, 6.0])
    exact = pfqn_ca(L, N, np.zeros(2))[1]
    coarse = abs(pfqn_aghq(L, N, None, 2)[1] - exact)
    assert abs(pfqn_aghq(L, N, None, q)[1] - exact) <= coarse


def test_an_all_zero_think_time_is_the_z_zero_branch():
    """pfqn_nc always passes sum(Z,1), so this is the shape a delay-free model arrives in."""
    L, N = _demands(4, 2, 91), np.array([7.0, 5.0])
    for f in (pfqn_aghq, pfqn_le):
        assert f(L, N, np.zeros(2))[1] == pytest.approx(f(L, N)[1], abs=1e-12)


@pytest.mark.parametrize("method", ["aghq"])
@pytest.mark.parametrize("Z", [np.zeros(2), np.array([1.5, 0.5])])
def test_it_is_reachable_through_the_nc_dispatcher(method, Z):
    L, N = _demands(3, 2, 64), np.array([5.0, 4.0])
    lG = pfqn_nc(L, N, Z, method)[1]
    assert np.isfinite(lG)
    assert abs(lG - pfqn_ca(L, N, Z)[1]) < 0.5
