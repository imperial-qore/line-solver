"""
pfqn_lld must reproduce pfqn_gld BIT FOR BIT.

Saturating the rate shift at the LLD threshold returns the same value, because
past it a further shift leaves the row unchanged over the columns the recursion
can still read; memoising the resulting repeated state changes what is
recomputed, never what is computed. Every terminal case is delegated back to
pfqn_gld on the materialised block, so exact equality is the right assertion.
"""

import numpy as np
import pytest

from line_solver.api.pfqn import pfqn_gld, pfqn_lld


def rates(kind, M, Ntot, rng):
    mu = np.ones((M, Ntot))
    lat = np.arange(1, Ntot + 1, dtype=float)
    for k in range(M):
        if kind == "multiserver":
            mu[k, :] = np.minimum(lat, 2 + (k % 3))
        elif kind == "loadindep":
            mu[k, :] = 1.0
        elif kind == "infserver":
            mu[k, :] = lat
        elif kind == "mixed":
            mu[k, :] = [np.ones(Ntot), np.minimum(lat, 2), lat][k % 3]
        elif kind == "settles_late":
            mu[k, :] = np.where(lat < 4, lat, 4.0 + k)
        elif kind == "arbitrary":
            mu[k, :] = 0.5 + rng.random(Ntot)
        else:
            raise ValueError(kind)
    return mu


@pytest.mark.parametrize(
    "kind",
    ["multiserver", "loadindep", "infserver", "mixed", "settles_late", "arbitrary"],
)
@pytest.mark.parametrize("M", [1, 2, 3, 4])
@pytest.mark.parametrize("R", [1, 2, 3])
@pytest.mark.parametrize("Ntot", [1, 3, 6])
def test_matches_gld_bitwise(kind, M, R, Ntot):
    rng = np.random.default_rng(hash((kind, M, R, Ntot)) % (2**32))
    N = np.full(R, Ntot // R, dtype=float)
    N[0] += Ntot - N.sum()
    L = 1 + 9 * rng.random((M, R))
    mu = rates(kind, M, Ntot, rng)
    assert pfqn_lld(L, N, mu).G == pfqn_gld(L, N, mu).G


def test_null_rates_default_identically():
    L = np.array([[2.0, 1.0], [1.0, 3.0]])
    N = np.array([2.0, 2.0])
    assert pfqn_lld(L, N, None).G == pfqn_gld(L, N, None).G


@pytest.mark.parametrize("Ntot", [3, 6])
def test_delay_row_never_settles(Ntot):
    """s_k = Ntot for a delay: the saturation must be a no-op, not a shortcut."""
    lat = np.arange(1, Ntot + 1, dtype=float)
    mu = np.vstack([np.minimum(lat, 2), np.minimum(lat, 3), lat])
    L = np.array([[2.0, 1.0], [1.0, 3.0], [0.5, 0.5]])
    N = np.array([float(Ntot // 2), float(Ntot - Ntot // 2)])   # integer populations
    assert pfqn_lld(L, N, mu).G == pfqn_gld(L, N, mu).G


def test_zero_population():
    L = np.array([[2.0, 1.0], [1.0, 3.0]])
    N = np.array([0.0, 0.0])
    assert pfqn_lld(L, N, np.ones((2, 1))).G == pfqn_gld(L, N, np.ones((2, 1))).G
