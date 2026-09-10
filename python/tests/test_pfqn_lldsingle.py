"""
pfqn_lldsingle must reproduce pfqn_gldsingle BIT FOR BIT.

The capped sweep performs a subset of the full sweep's operations, never a
rearrangement of them, so the two constants are equal in the last bit and not
merely to a tolerance. Every assertion here is exact equality for that reason:
a tolerance would hide precisely the regression these tests exist to catch.
"""

import numpy as np
import pytest

from line_solver.api.pfqn import pfqn_gldsingle, pfqn_lldsingle


def rates(kind, M, Ntot, rng):
    """Rate matrices spanning every threshold regime the cap has to handle."""
    mu = np.ones((M, Ntot))
    lattice = np.arange(1, Ntot + 1)
    for k in range(M):
        if kind == "loadindep":            # s_k = 1
            mu[k, :] = 1.0
        elif kind == "multiserver":        # s_k = server count
            mu[k, :] = np.minimum(lattice, 1 + (k % 3) + 1)
        elif kind == "infserver":          # never settles, s_k = Ntot
            mu[k, :] = lattice
        elif kind == "mixed":
            mu[k, :] = [np.ones(Ntot), np.minimum(lattice, 2), lattice][k % 3]
        elif kind == "settles_late":       # threshold strictly inside the row
            mu[k, :] = np.where(lattice < 4, lattice, 4.0 + k)
        elif kind == "arbitrary":          # no exploitable structure at all
            mu[k, :] = 0.5 + rng.random(Ntot)
        else:
            raise ValueError(kind)
    return mu


@pytest.mark.parametrize(
    "kind",
    ["loadindep", "multiserver", "infserver", "mixed", "settles_late", "arbitrary"],
)
@pytest.mark.parametrize("M", [1, 2, 5, 9])
@pytest.mark.parametrize("Ntot", [1, 2, 5, 16, 33])
def test_matches_gldsingle_bitwise(kind, M, Ntot):
    rng = np.random.default_rng(hash((kind, M, Ntot)) % (2**32))
    L = 1 + 9 * rng.random((M, 1))
    mu = rates(kind, M, Ntot, rng)
    N = np.array([Ntot], dtype=float)
    assert pfqn_lldsingle(L, N, mu).lG == pfqn_gldsingle(L, N, mu).lG


@pytest.mark.parametrize("Ntot", [5, 16])
def test_zero_demand_station(Ntot):
    """A zero demand sends lL to -inf; the cap must not turn that into a nan."""
    rng = np.random.default_rng(7)
    L = 1 + 9 * rng.random((4, 1))
    L[0, 0] = 0.0
    mu = rates("multiserver", 4, Ntot, rng)
    N = np.array([Ntot], dtype=float)
    assert pfqn_lldsingle(L, N, mu).lG == pfqn_gldsingle(L, N, mu).lG


@pytest.mark.parametrize("Ntot", [5, 16])
def test_negative_demand_uses_linear_branch(Ntot):
    """Negative demands leave the log domain, which is a separate code path."""
    rng = np.random.default_rng(11)
    L = -(1 + rng.random((3, 1)))
    mu = rates("multiserver", 3, Ntot, rng)
    N = np.array([Ntot], dtype=float)
    assert pfqn_lldsingle(L, N, mu).G == pfqn_gldsingle(L, N, mu).G


@pytest.mark.parametrize("Ntot", [5, 16])
def test_complex_demand_uses_linear_branch(Ntot):
    """Complex demands take the same branch and the spelled-out product."""
    rng = np.random.default_rng(13)
    L = rng.random((3, 1)) + 1j * rng.random((3, 1))
    mu = rates("multiserver", 3, Ntot, rng)
    N = np.array([Ntot], dtype=float)
    assert pfqn_lldsingle(L, N, mu).G == pfqn_gldsingle(L, N, mu).G


def test_empty_population():
    """N = 0 is the empty product, as pfqn_gldsingle returns from unrun loops."""
    L = np.array([[1.0], [2.0]])
    mu = np.ones((2, 1))
    res = pfqn_lldsingle(L, np.array([0.0]), mu)
    assert res.G == 1.0 and res.lG == 0.0


def test_multiclass_is_refused():
    L = np.ones((2, 3))
    mu = np.ones((2, 4))
    with pytest.raises(RuntimeError, match="multiclass"):
        pfqn_lldsingle(L, np.array([4.0]), mu)


def test_threshold_scan_does_not_merge_distinct_rates():
    """
    A FALSE tie would be a wrong answer, so the scan must not collapse a row
    whose tail only nearly repeats. Rates one ulp apart stay distinct.
    """
    Ntot = 8
    mu = np.ones((1, Ntot))
    mu[0, :] = 3.0
    mu[0, Ntot - 3] = np.nextafter(3.0, 4.0) * (1 + 1e-9)
    L = np.array([[0.7]])
    N = np.array([float(Ntot)])
    assert pfqn_lldsingle(L, N, mu).lG == pfqn_gldsingle(L, N, mu).lG


@pytest.mark.parametrize("Ntot", [4, 9, 17])
def test_infinite_tail_does_not_swallow_finite_rates(Ntot):
    """
    A row that ends in +inf must not tie its FINITE entries to that tail.

    With tail = inf the tolerance bound eps*max(abs(tail), 1) is itself inf and
    abs(prev - inf) <= inf holds for every finite prev, which would collapse the
    entire row onto alpha = inf and zero the constant. pfqn_rd reaches exactly
    this input, since it maps nan rates to inf before calling in.
    """
    mu = np.empty((2, Ntot))
    lattice = np.arange(1, Ntot + 1)
    mu[0, :] = np.minimum(lattice, 2)
    mu[1, :] = np.where(lattice < Ntot - 1, lattice.astype(float), np.inf)
    L = np.array([[0.7], [0.4]])
    N = np.array([float(Ntot)])
    assert pfqn_lldsingle(L, N, mu).lG == pfqn_gldsingle(L, N, mu).lG


def test_all_infinite_row_still_collapses():
    """An entirely infinite row ties with itself and must reach s_k = 1."""
    Ntot = 8
    mu = np.ones((2, Ntot))
    mu[0, :] = np.minimum(np.arange(1, Ntot + 1), 3)
    mu[1, :] = np.inf
    L = np.array([[0.7], [0.4]])
    N = np.array([float(Ntot)])
    assert pfqn_lldsingle(L, N, mu).lG == pfqn_gldsingle(L, N, mu).lG
