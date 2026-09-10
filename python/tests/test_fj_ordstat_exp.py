"""
Tests for E[X_(k)] of independent exponential branches, the instant a k-of-n
(quorum) join fires, against the closed forms that hold independently of any
implementation: the minimum is 1/sum(lambda), the i.i.d. case is a partial
harmonic sum, and k = n reproduces the alternating inclusion-exclusion series
the fork-join fixed point used before quorum joins were supported.
"""

from itertools import combinations

import pytest

from line_solver.api.fjnative import fj_ordstat_exp


def test_minimum_is_reciprocal_of_the_rate_sum():
    ri = [1 / 3, 1 / 4, 1 / 2]
    assert fj_ordstat_exp(ri, 1) == pytest.approx(1 / 9, abs=1e-12)


def test_iid_case_is_a_partial_harmonic_sum():
    ri = [1.0] * 4
    assert fj_ordstat_exp(ri, 1) == pytest.approx(1 / 4, abs=1e-12)
    assert fj_ordstat_exp(ri, 2) == pytest.approx(1 / 4 + 1 / 3, abs=1e-12)
    assert fj_ordstat_exp(ri, 3) == pytest.approx(1 / 4 + 1 / 3 + 1 / 2, abs=1e-12)
    assert fj_ordstat_exp(ri, 4) == pytest.approx(1 / 4 + 1 / 3 + 1 / 2 + 1, abs=1e-12)


def test_maximum_matches_the_classical_series():
    """k = n must reproduce E[max] TERM BY TERM: that is what keeps a standard
    join bit-identical to the pre-quorum fork-join fixed point."""
    ri = [0.7, 1.3, 2.9, 0.4]
    lambdai = [1 / x for x in ri]
    expected = 0.0
    for j in range(1, len(ri) + 1):
        term = sum(1.0 / sum(c) for c in combinations(lambdai, j))
        expected += ((-1) ** (j - 1)) * term
    assert fj_ordstat_exp(ri, len(ri)) == pytest.approx(expected, abs=1e-12)


def test_monotone_in_k():
    ri = [0.5, 1.5, 2.5, 3.5, 4.5]
    values = [fj_ordstat_exp(ri, k) for k in range(1, len(ri) + 1)]
    assert all(b > a for a, b in zip(values, values[1:]))


def test_zero_branch_completes_instantly():
    """A branch of zero mean counts toward the quorum at once and never delays
    the join; the MATLAB reference reaches the same values through 1/Inf."""
    assert fj_ordstat_exp([0.0, 2.0], 1) == pytest.approx(0.0, abs=1e-14)
    assert fj_ordstat_exp([0.0, 2.0], 2) == pytest.approx(2.0, abs=1e-14)


def test_degenerate_branch_sets():
    assert fj_ordstat_exp([3.0], 1) == pytest.approx(3.0, abs=1e-14)
    assert fj_ordstat_exp([], 1) == 0.0


def test_quorum_out_of_range_raises():
    with pytest.raises(ValueError):
        fj_ordstat_exp([1.0, 2.0], 0)
    with pytest.raises(ValueError):
        fj_ordstat_exp([1.0, 2.0], 3)
