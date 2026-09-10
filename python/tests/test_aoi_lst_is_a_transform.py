"""The AoI LST returned by the FCFS routines must itself BE an LST.

Nothing asserted that until 2026-08-19, which is how BOTH rows shipped an
expression that was not one, in all four codebases:

  M/GI/1  returned (lambda*H*(s))/(s + lambda - lambda*H*(s)) * W*(s), whose
          first factor diverges as s -> 0 (the denominator vanishes like
          s*(1+rho)), so A*(0) was +inf and A*(s) exceeded 1 for small s.
  GI/M/1  returned (mu*sigma(s))/(s + mu - mu*sigma(s)) * D*(s), which is
          sigma/(1-sigma) at the origin rather than 1.

The mean and peak AoI columns were correct throughout, so every existing
assertion passed over both. See BUGS.md (Known open items, lstAoI) and
_kb/03-api-layer.md (AoI family).
"""

import numpy as np
import pytest

from line_solver.lib.thirdparty.aoi.analytical import aoi_fcfs_mgi1, aoi_fcfs_gim1


def _erlang_lst(k, rate):
    return lambda s: (rate / (rate + s)) ** k


# M/E2/1: lambda 0.6, Erlang(2, 4) service (mean 0.5, second moment 0.375).
MGI1 = dict(lambd=0.6, H_lst=_erlang_lst(2, 4.0), E_H=0.5, E_H2=0.375)
# E2/M/1: Erlang(2, 1.2) interarrivals, exponential service at rate 1.
GIM1 = dict(Y_lst=_erlang_lst(2, 1.2), mu=1.0, E_Y=2.0 / 1.2, E_Y2=6.0 / 1.44)


def _rows():
    m, lst_m, _ = aoi_fcfs_mgi1(**MGI1)
    g, lst_g, _ = aoi_fcfs_gim1(**GIM1)
    return [("M/E2/1", m, lst_m), ("E2/M/1", g, lst_g)]


@pytest.mark.parametrize("name,mean,lst", _rows())
def test_aoi_lst_is_unity_at_the_origin(name, mean, lst):
    """A*(0) = 1. This is the check that fails outright on the old expressions."""
    assert float(np.real(lst(0.0))) == pytest.approx(1.0, abs=1e-9)


@pytest.mark.parametrize("name,mean,lst", _rows())
def test_aoi_lst_is_in_range_and_nonincreasing(name, mean, lst):
    """0 < A*(s) <= 1 and nonincreasing -- the "above 1 at small s" signature."""
    prev = 1.0
    for s in (0.05, 0.2, 0.5, 1.0, 2.0):
        v = float(np.real(lst(s)))
        assert 0.0 < v <= 1.0 + 1e-12, "A*(%g) = %g is outside (0,1]" % (s, v)
        assert v <= prev + 1e-12, "A*(%g) = %g rose above A*(previous)" % (s, v)
        prev = v


@pytest.mark.parametrize("name,mean,lst", _rows())
def test_aoi_lst_slope_at_zero_is_the_mean(name, mean, lst):
    """-A*'(0) = E[A], tying the transform to the mean column beside it."""
    h = 1e-5
    slope = -(float(np.real(lst(h))) - float(np.real(lst(0.0)))) / h
    assert slope == pytest.approx(mean, rel=1e-2)


def test_aoi_lst_matches_a_sample_path():
    """Against simulation of the same two queues.

    Values are the cycle average of exp(-s*age) over departure intervals, from
    4e6 (M/E2/1) and 3e6 (E2/M/1) cycles of a Lindley recursion. The old forms
    gave 1.3062 and 0.23056 at the first point of each row.
    """
    _, lst_m, _ = aoi_fcfs_mgi1(**MGI1)
    for s, simulated in ((0.3, 0.569475), (1.0, 0.228510), (2.0, 0.093143)):
        assert float(np.real(lst_m(s))) == pytest.approx(simulated, abs=5e-4)

    _, lst_g, _ = aoi_fcfs_gim1(**GIM1)
    for s, simulated in ((0.2, 0.593321), (0.8, 0.198780), (2.0, 0.054780)):
        assert float(np.real(lst_g(s))) == pytest.approx(simulated, abs=5e-4)


def test_aoi_lst_accepts_an_array():
    """The array path is a second branch in each routine and carried the same bug."""
    _, lst_m, _ = aoi_fcfs_mgi1(**MGI1)
    got = np.real(lst_m(np.array([0.3, 1.0, 2.0])))
    assert got.shape == (3,)
    np.testing.assert_allclose(got, [0.569354, 0.228429, 0.093102], atol=1e-6)
