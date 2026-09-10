"""Boundaries of the SPM cache normalizing constant.

cache_spm returns Ehat(m) = E(m) * prod_l m_l!, the same constant cache_erec
evaluates exactly, so cache_erec is the oracle throughout.

Two boundaries used to be mishandled in every codebase and are covered here:

  m_l == 0     list l has xi(l)=0, a boundary of the Laplace integral rather
               than a direction of it. Left in, the -(1/2) sum_l log xi(l)
               prefactor gained ~+17 per empty list: n=12 at gamma=(.8,.6,.4)
               and m=(0,3,3) returned lZ=24.814 against an exact 9.127.
  n == sum(m)  every item is cached, the multipliers diverge and the
               fixed-point iteration cannot converge. It used to be run
               anyway, and did not terminate.
"""
import math

import numpy as np
import pytest

from line_solver.api.cache import cache_erec, cache_spm
from line_solver.api.cache.spm import cache_prob_spm

UNIFORM = np.tile([[0.8, 0.6, 0.4]], (12, 1))


def skewed(n=12):
    """Non-uniform popularity, monotone across the list index as SPM assumes."""
    g = 0.3 + 2.7 * (np.arange(1, n + 1) / n)
    return np.column_stack([g, 0.6 * g, 0.3 * g])


def exact(gamma, m):
    return math.log(cache_erec(gamma, np.asarray(m, float)))


# --------------------------------------------------------------------------
# m_l == 0: a zero-capacity list must leave the expansion
# --------------------------------------------------------------------------

@pytest.mark.parametrize("gamma", [UNIFORM, skewed()])
@pytest.mark.parametrize("m,keep", [((0, 3, 3), [1, 2]),
                                    ((3, 0, 3), [0, 2]),
                                    ((0, 0, 6), [2]),
                                    ((4, 4, 0), [0, 1])])
def test_zero_capacity_list_equals_dropping_it(gamma, m, keep):
    """Dropping an empty list is exact, so the two calls must agree bit for bit."""
    _, with_empty, xi = cache_spm(gamma, np.array(m, float))
    _, without, _ = cache_spm(gamma[:, keep], np.array([m[l] for l in keep], float))
    assert with_empty == without
    # the dropped lists carry xi = 0, the root of their capacity equation
    for l in range(3):
        assert (xi[l] == 0.0) == (m[l] == 0)


@pytest.mark.parametrize("gamma", [UNIFORM, skewed()])
@pytest.mark.parametrize("m", [(0, 3, 3), (3, 0, 3), (0, 0, 6), (4, 4, 0)])
def test_zero_capacity_list_stays_near_the_exact_constant(gamma, m):
    """The old prefactor blow-up was ~1e8 relative; SPM's own error is ~10%."""
    _, lZ, _ = cache_spm(gamma, np.array(m, float))
    assert lZ == pytest.approx(exact(gamma, m), abs=0.30)


def test_empty_cache_is_the_unit_constant():
    Z, lZ, xi = cache_spm(UNIFORM, np.zeros(3))
    assert Z == 1.0
    assert lZ == 0.0
    assert np.all(xi == 0.0)


# --------------------------------------------------------------------------
# n == sum(m): degenerate saddle, exact fallback, no iteration
# --------------------------------------------------------------------------

@pytest.mark.parametrize("gamma", [UNIFORM, skewed()])
@pytest.mark.parametrize("m", [(4, 4, 4), (12, 0, 0), (6, 5, 1)])
def test_full_cache_returns_the_exact_constant_and_terminates(gamma, m):
    Z, lZ, xi = cache_spm(gamma, np.array(m, float))
    assert lZ == pytest.approx(exact(gamma, m), rel=1e-12)
    assert Z == pytest.approx(cache_erec(gamma, np.array(m, float)), rel=1e-12)
    assert np.all(np.isinf(xi))          # the multipliers' limit, not a number


# --------------------------------------------------------------------------
# interior points must be untouched by the boundary handling
# --------------------------------------------------------------------------

@pytest.mark.parametrize("m,want", [((3, 3, 3), 13.34844416820355),
                                    ((4, 3, 2), 14.04836686545269),
                                    ((2, 2, 2), 10.23839884632498),
                                    ((6, 4, 1), 15.878606855246987)])
def test_interior_saddle_is_unchanged(m, want):
    _, lZ, _ = cache_spm(UNIFORM, np.array(m, float))
    assert lZ == pytest.approx(want, rel=1e-12)


# --------------------------------------------------------------------------
# cache_prob_spm reads xi, so it inherits both boundaries
#
# NOTE, cross-codebase: this native cache_prob_spm returns the per-item MISS
# probability 1/(1+S_k), whereas cache_prob_spm.m and Cache_prob_spm.java
# return the (n, h+1) placement matrix P(item i in list l). That divergence
# predates this file and is not what these tests are about.
# --------------------------------------------------------------------------

@pytest.mark.parametrize("gamma", [UNIFORM, skewed()])
def test_prob_spm_ignores_an_empty_list(gamma):
    """xi=0 on the empty list, so it must not shift any miss probability."""
    with_empty = np.asarray(cache_prob_spm(gamma, np.array([0.0, 3.0, 2.0])))
    without = np.asarray(cache_prob_spm(gamma[:, [1, 2]], np.array([3.0, 2.0])))
    # not bit-for-bit: S = gamma @ xi contracts over one more (zero) column here
    assert np.allclose(with_empty, without, rtol=1e-12, atol=0.0)
    assert np.all((with_empty >= 0.0) & (with_empty <= 1.0))


@pytest.mark.parametrize("gamma", [UNIFORM, skewed()])
def test_prob_spm_at_a_full_cache_misses_nothing(gamma):
    """Every item is cached, so xi diverges and the miss probability is 0."""
    prob = np.asarray(cache_prob_spm(gamma, np.array([4.0, 4.0, 4.0])))
    assert np.all(prob == 0.0)
