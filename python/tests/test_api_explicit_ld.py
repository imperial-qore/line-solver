"""Validation of pfqn_explicit_ld, the limited load-dependent explicit constant.

Theorem 1 of Casale, Harrison and Ong (Perform. Eval. 2021), carried over the
divided-difference form of Casale (SIGMETRICS 2017), Corollary 3.2. The closed
form is EXACT, so it is checked against oracles rather than for self-consistency:

1. a multi-dimensional convolution over the state space, in exact rational
   arithmetic, on multi-server, fixed-rate, arbitrary limited load-dependent and
   never-settling rate lattices;
2. pfqn_explicit itself on the fixed-rate degeneration mu = 1, which the
   load-dependent route must reproduce to the bit;
3. the near-tie contract: an induced tie that lands two ulps apart is MISSED at
   tol = eps and the loss report has to say so, while any looser tolerance
   merges the pair and Eq. (16) is exact.
"""

import itertools
import math
from fractions import Fraction

import numpy as np
import pytest

from line_solver.api.pfqn import pfqn_explicit, pfqn_explicit_ld
from line_solver.api.pfqn.ncld import pfqn_gld


def _convolve_exact(L, N, mu):
    """G(N) by load-dependent convolution over the state space, in Fractions."""
    L = np.asarray(L, dtype=float)
    N = np.asarray(N, dtype=int)
    mu = np.asarray(mu, dtype=float)
    M, R = L.shape
    states = [tuple(n) for n in itertools.product(*[range(int(x) + 1) for x in N])]

    def station(k, n):
        t = sum(n)
        v = Fraction(math.factorial(t))
        for u in range(1, t + 1):
            v /= Fraction(mu[k][u - 1]).limit_denominator(10 ** 9)
        for r in range(R):
            v *= (Fraction(L[k][r]).limit_denominator(10 ** 9) ** n[r]
                  / Fraction(math.factorial(n[r])))
        return v

    g = {s: Fraction(0) for s in states}
    g[tuple([0] * R)] = Fraction(1)
    for k in range(M):
        ng = {s: Fraction(0) for s in states}
        for s in states:
            for n in itertools.product(*[range(s[r] + 1) for r in range(R)]):
                ng[s] += g[tuple(s[r] - n[r] for r in range(R))] * station(k, n)
        g = ng
    return float(g[tuple(int(x) for x in N)])


def _ms_rates(M, Nt, c):
    """mu(n) = min(n,c) at each of M centers."""
    return np.tile(np.minimum(np.arange(1, Nt + 1), c), (M, 1)).astype(float)


@pytest.mark.parametrize('tag,L,N,mu', [
    ('multiserver R=2',
     [[1.2, 0.7], [0.4, 1.9]], [3, 2],
     [[1, 2, 2, 2, 2], [1, 2, 3, 3, 3]]),
    ('fixed rate R=3',
     [[1.0, 0.7, 0.3], [0.5, 1.3, 0.9], [0.9, 0.4, 1.7]], [2, 1, 2],
     np.ones((3, 5))),
    # rates that neither increase nor follow a multi-server shape, settling at s=3
    ('arbitrary LLD R=1',
     [[1.1], [0.6], [2.3]], [6],
     [[0.5, 1.4, 2.2, 2.2, 2.2, 2.2],
      [1.0, 0.8, 1.7, 1.7, 1.7, 1.7],
      [2.0, 1.1, 0.9, 0.9, 0.9, 0.9]]),
    # an infinite server never settles at a finite s, but s = |N| is admissible
    # because no larger population occurs
    ('never settles R=2',
     [[0.9, 1.4], [1.1, 0.5]], [2, 2],
     [[1, 2, 3, 4], [1, 2, 2, 2]]),
    ('zero column R=2',
     [[0.0, 1.3], [0.9, 0.4]], [2, 2],
     [[1, 2, 2, 2], [1, 2, 2, 2]]),
])
def test_matches_the_convolution(tag, L, N, mu):
    L = np.asarray(L, dtype=float)
    N = np.asarray(N, dtype=float)
    mu = np.asarray(mu, dtype=float)
    ref = _convolve_exact(L, N, mu)
    lG, G, _, loss = pfqn_explicit_ld(L, N, mu)
    assert np.isfinite(lG), tag
    assert loss < 15, tag
    assert G == pytest.approx(ref, rel=1e-8), tag


def test_fixed_rate_degeneration_reproduces_pfqn_explicit():
    L = np.array([[1.0, 0.7], [0.5, 1.3], [0.9, 0.4]])
    N = np.array([3.0, 2.0])
    a = pfqn_explicit(L, N)
    b = pfqn_explicit_ld(L, N, np.ones((3, 5)))
    assert a[2] == b[2]
    assert a[0] == pytest.approx(b[0], rel=0, abs=1e-12)


def test_repeated_scaled_demands_take_equation_sixteen():
    L = np.array([[1.3, 0.8], [1.3, 0.8], [1.3, 0.8]])
    N = np.array([2.0, 2.0])
    mu = _ms_rates(3, 4, 2)
    lG, G, method, _ = pfqn_explicit_ld(L, N, mu)
    assert method == 'repeated'
    assert G == pytest.approx(_convolve_exact(L, N, mu), rel=1e-8)


def test_single_class_skips_the_outer_sum():
    # h_theta(N) is homogeneous of degree N in theta, so the divided difference
    # is the identity at R=1 and the constant must still match the convolution
    L = np.array([[1.4], [0.9], [0.35]])
    N = np.array([7.0])
    mu = _ms_rates(3, 7, 2)
    _, G, _, _ = pfqn_explicit_ld(L, N, mu)
    assert G == pytest.approx(_convolve_exact(L, N, mu), rel=1e-8)


def test_pfqn_gld_agrees_on_a_zero_demand_load_dependent_model():
    """Regression for the pfqn_gld base case (fixed 2026-09-03).

    A class with jobs and no demand at the ONLY station used to be dropped from
    the single-station sum instead of zeroing the constant, and the recursion
    reaches that base case with the full population every time it peels a
    station. So any load-dependent model carrying a zero demand came back wrong.
    Fixed rate was unaffected, which is why it hid for so long.
    """
    N = np.array([2.0, 2.0])
    mu = _ms_rates(2, 4, 2)
    for L, exact in [(np.array([[0.0, 1.3], [0.9, 0.4]]), 0.755325),
                     (np.array([[0.0, 0.0], [0.9, 0.4]]), None)]:
        ref = _convolve_exact(L, N, mu) if exact is None else exact
        assert pfqn_gld(L, N, mu, None).G == pytest.approx(ref, rel=1e-9)
    # the fixed-rate model must be untouched by the fix
    L = np.array([[0.0, 1.3], [0.9, 0.4]])
    assert pfqn_gld(L, N, np.ones((2, 4)), None).G == pytest.approx(
        _convolve_exact(L, N, np.ones((2, 4))), rel=1e-9)


def test_an_ulp_wide_tie_is_missed_at_eps_and_reported():
    # both scaled demands are 1.95 at t=[2 3], but land two ulps apart in doubles
    L = np.array([[0.0, 1.3], [0.9, 0.7]])
    N = np.array([2.0, 3.0])
    mu = _ms_rates(2, 5, 2)
    ref = _convolve_exact(L, N, mu)

    lG, _, method, loss = pfqn_explicit_ld(L, N, mu)
    assert method == 'distinct'
    assert loss > 15
    assert abs(lG - math.log(ref)) > 1.0

    lG, _, method, _ = pfqn_explicit_ld(L, N, mu, tol=1e-12)
    assert method == 'repeated'
    assert lG == pytest.approx(math.log(ref), rel=1e-8)

    # a caller holding a budget is refused rather than warned
    lG, _, _, _ = pfqn_explicit_ld(L, N, mu, maxloss=8.0)
    assert math.isnan(lG)


def test_inadmissible_arguments_are_refused():
    L = np.array([[1.0, 0.5], [0.5, 1.0]])
    N = np.array([2.0, 2.0])
    with pytest.raises(ValueError):
        pfqn_explicit_ld(L, N, np.ones((2, 4)), method='bogus')
    # a rate lattice shorter than the population cannot answer at |N|
    with pytest.raises(ValueError):
        pfqn_explicit_ld(L, N, np.ones((2, 3)))
    # and a rate of zero would divide by zero in phi_k
    bad = np.ones((2, 4))
    bad[0, 2] = 0.0
    with pytest.raises(ValueError):
        pfqn_explicit_ld(L, N, bad)


def test_the_explicit_token_routes_through_pfqn_ncld():
    from line_solver.api.pfqn.ncld import pfqn_ncld
    from line_solver.solvers.solver_nc.solver_nc import SolverNC

    L = np.array([[1.2, 0.7], [0.4, 1.9]])
    N = np.array([3.0, 2.0])
    Z = np.zeros((1, 2))
    mu = np.vstack([np.minimum(np.arange(1, 6), 2),
                    np.minimum(np.arange(1, 6), 3)]).astype(float)

    exact = SolverNC.default_options()
    exact.method = 'exact'
    ref = pfqn_ncld(L, N, Z, mu, exact)

    opts = SolverNC.default_options()
    opts.method = 'divdiff'
    got = pfqn_ncld(L, N, Z, mu, opts)
    assert got.method == 'divdiff.ld/distinct'
    assert got.lG == pytest.approx(ref.lG, rel=1e-8)

    # a think time would have to enter g_sigma, whose closed form covers queues only
    with pytest.raises(ValueError):
        pfqn_ncld(L, N, np.array([[0.5, 0.5]]), mu, opts)
