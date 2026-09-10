"""Tests for fj_tail_ordstat, the k-of-n (quorum) fork-join tail.

A quorum join fires on the k-th of n siblings, so the request response time is
the k-th ORDER STATISTIC of the branch times and not their maximum. Reading the
maximum there returns the AND-join tail under a quorum's name -- the same number
for every k, which is what the forktail route did before this leaf existed.

The reference values are agreed across MATLAB, the JAR, native python and the C++
port, and were checked against a 400k-sample Monte Carlo of the SAME fitted GE
branches (so the check is of the order-statistic inversion, not of the GE fit):
every k agreed to within 1.2% at the 99th percentile, which is the sampling error
there.
"""

import numpy as np
import pytest

from line_solver.api.fjnative import fj_tail_forktail, fj_tail_ordstat

HOM = dict(ET=2.0, VT=6.0, n=4)
HET = dict(ET=[1.0, 2.0, 4.0], VT=[1.0, 8.0, 40.0])

REF_HOM = {  # (k, p) -> percentile
    (1, 0.90): 0.871483, (1, 0.99): 2.179109,
    (2, 0.90): 2.146981, (2, 0.99): 4.220359,
    (3, 0.90): 4.189345, (3, 0.99): 7.427549,
    (4, 0.90): 8.718015, (4, 0.99): 15.050364,
}
REF_HET = {
    (1, 0.90): 0.980116, (1, 0.99): 2.397990,
    (2, 0.90): 3.006171, (2, 0.99): 7.388476,
    (3, 0.90): 12.421746, (3, 0.99): 29.881069,
}


@pytest.mark.parametrize('key', sorted(REF_HOM))
def test_homogeneous_reference(key):
    k, p = key
    got = fj_tail_ordstat(HOM['ET'], HOM['VT'], HOM['n'], p, k)[0]
    assert got == pytest.approx(REF_HOM[key], rel=1e-5)


@pytest.mark.parametrize('key', sorted(REF_HET))
def test_heterogeneous_reference(key):
    k, p = key
    got = fj_tail_ordstat(HET['ET'], HET['VT'], None, p, k)[0]
    assert got == pytest.approx(REF_HET[key], rel=1e-5)


@pytest.mark.parametrize('p', [0.5, 0.9, 0.99, 0.999])
def test_full_join_is_forktail_exactly(p):
    """k = n must reproduce the AND-join bit for bit, so no existing result moves."""
    a = fj_tail_ordstat(HOM['ET'], HOM['VT'], HOM['n'], p, HOM['n'])[0]
    b = fj_tail_forktail(HOM['ET'], HOM['VT'], HOM['n'], p)[0]
    assert a == b
    c = fj_tail_ordstat(HET['ET'], HET['VT'], None, p, len(HET['ET']))[0]
    d = fj_tail_forktail(HET['ET'], HET['VT'], None, p)[0]
    assert c == d


def test_percentile_grows_with_the_quorum():
    """Waiting for more siblings can only take longer."""
    for p in (0.5, 0.9, 0.99):
        xs = [fj_tail_ordstat(HOM['ET'], HOM['VT'], HOM['n'], p, k)[0]
              for k in range(1, HOM['n'] + 1)]
        assert all(xs[i] < xs[i + 1] for i in range(len(xs) - 1))
        ys = [fj_tail_ordstat(HET['ET'], HET['VT'], None, p, k)[0] for k in (1, 2, 3)]
        assert all(ys[i] < ys[i + 1] for i in range(len(ys) - 1))


def test_k_equals_one_is_the_first_completion():
    """With one branch the k=1 order statistic IS that branch."""
    xp = fj_tail_ordstat(HOM['ET'], HOM['VT'], 1, 0.99, 1)[0]
    ref = fj_tail_forktail(HOM['ET'], HOM['VT'], 1, 0.99)[0]
    assert xp == pytest.approx(ref, rel=1e-12)


def test_matches_monte_carlo_of_the_fitted_branches():
    """The inversion, checked against sampling the very law it fits."""
    rng = np.random.default_rng(23000)
    _, alpha, beta = fj_tail_ordstat(HET['ET'], HET['VT'], None, 0.9, 3)
    ns = 200000
    S = np.empty((3, ns))
    for i in range(3):
        S[i, :] = -beta[i] * np.log(1.0 - rng.random(ns) ** (1.0 / alpha[i]))
    S.sort(axis=0)
    for k in (1, 2, 3):
        for p in (0.9, 0.99):
            got = fj_tail_ordstat(HET['ET'], HET['VT'], None, p, k)[0]
            mc = float(np.quantile(S[k - 1, :], p))
            assert abs(got - mc) / mc < 0.03


def test_quorum_out_of_range_is_refused():
    with pytest.raises(ValueError):
        fj_tail_ordstat(HOM['ET'], HOM['VT'], HOM['n'], 0.99, 0)
    with pytest.raises(ValueError):
        fj_tail_ordstat(HOM['ET'], HOM['VT'], HOM['n'], 0.99, HOM['n'] + 1)
    with pytest.raises(ValueError):
        fj_tail_ordstat(HET['ET'], HET['VT'], None, 0.99, 4)


def test_percentile_out_of_range_is_refused():
    with pytest.raises(ValueError):
        fj_tail_ordstat(HOM['ET'], HOM['VT'], HOM['n'], 0.0, 2)
    with pytest.raises(ValueError):
        fj_tail_ordstat(HOM['ET'], HOM['VT'], HOM['n'], 1.0, 2)
