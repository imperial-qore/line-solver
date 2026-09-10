"""
The six approximate MVA algorithms of Chapter 2 of H. Wang, "Approximate MVA
Algorithms for Solving Queueing Network Models", M.Sc. thesis, University of
Toronto, 1997, that LINE lacked until 2026-08-03: Bard LCP, Chow SA, the
Hsieh-Lam PAM family, Eager Looping, dSeS-Lavenberg-Muntz Clustering and the
dSeS-Muntz Improved Linearizer.

The reference values are the MATLAB, JAR and C++ ports on the same input; all
four agree to the printed digits. Two of the checks pin PROPERTIES rather than
numbers, and those are the ones that catch a bad transcription:

- dmlin MUST equal linearizer. IL is a cost reduction, not an approximation;
  survey eq. (2.50) drops the Delta^(j) correction and breaks it.
- looping MUST bracket the exact solution, and lcp MUST be worse than bs,
  which is the accuracy ordering of Table 2.1.
"""

import numpy as np
import pytest

from line_solver.api.pfqn import (pfqn_bs, pfqn_chow, pfqn_clust, pfqn_dmlin,
                                  pfqn_lcp, pfqn_linearizer, pfqn_looping,
                                  pfqn_mva, pfqn_pam)
from line_solver.lang.base import SchedStrategy

L = np.array([[0.10, 0.30], [0.20, 0.05], [0.40, 0.15]])
N = np.array([5.0, 4.0])
Z = np.array([1.0, 2.0])


def _exact():
    return np.asarray(pfqn_mva(L, N, Z)[0]).flatten()


def _xerr(X):
    Xex = _exact()
    return float(np.max(np.abs(np.asarray(X).flatten() - Xex) / Xex))


def test_lcp_matches_the_other_ports():
    X = np.asarray(pfqn_lcp(L, N, Z)[0]).flatten()
    assert X == pytest.approx([1.502598, 1.188360], rel=1e-6)


def test_lcp_is_less_accurate_than_bard_schweitzer():
    # Table 2.1 ranks LCP below the PE algorithm it seeded: dropping the
    # (N-1)/N factor can only over-count the queue the arriving customer sees.
    assert _xerr(pfqn_lcp(L, N, Z)[0]) > _xerr(pfqn_bs(L, N, Z)[0])


def test_chow_matches_the_other_ports():
    X = np.asarray(pfqn_chow(L, N, Z)[0]).flatten()
    assert X == pytest.approx([1.721975, 1.250086], rel=1e-6)


def test_chow_beats_bard_schweitzer():
    # Rank 3 against PE's rank 4; a theta-correction of zero would show here.
    assert _xerr(pfqn_chow(L, N, Z)[0]) < _xerr(pfqn_bs(L, N, Z)[0])


def test_chow_backward_estimator_is_a_different_answer():
    fwd = np.asarray(pfqn_chow(L, N, Z)[0]).flatten()
    bwd = np.asarray(pfqn_chow(L, N, Z, variant='backward')[0]).flatten()
    assert not np.allclose(fwd, bwd, rtol=1e-9)


def test_pam_variants():
    b = np.asarray(pfqn_pam(L, N, Z, 'pamb')[0]).flatten()
    assert b == pytest.approx([1.351351, 1.024515], rel=1e-6)
    # No centre is overloaded here, so the PAMI capping is inert and the two
    # must coincide; a capping applied unconditionally would not.
    i = np.asarray(pfqn_pam(L, N, Z, 'pami')[0]).flatten()
    assert i == pytest.approx(b, rel=1e-12)
    # PAMT unrolls one MVA step more, which moves the answer.
    t = np.asarray(pfqn_pam(L, N, Z, 'pamt')[0]).flatten()
    assert t == pytest.approx([1.672982, 1.197762], rel=1e-6)


def test_dmlin_is_linearizer():
    # The defining claim of de Souza e Silva and Muntz (1990): same fixed
    # point, lower cost. Transcribing survey eq. (2.50) literally breaks this.
    d = pfqn_dmlin(L, N, Z)
    l = pfqn_linearizer(L, N, Z, [SchedStrategy.PS] * L.shape[0])
    assert np.asarray(d[5]).flatten() == pytest.approx(np.asarray(l[5]).flatten(), rel=1e-9)
    assert np.asarray(d[0]) == pytest.approx(np.asarray(l[0]), rel=1e-9)
    assert np.asarray(d[5]).flatten()[0] == pytest.approx(1.790568, rel=1e-6)


def test_clust_matches_the_other_ports():
    X = np.asarray(pfqn_clust(L, N, Z)[0]).flatten()
    assert X == pytest.approx([1.799264, 1.241524], rel=1e-6)


def test_clust_with_pe_inside_is_a_different_algorithm():
    # Linearizer inside is the only setting that buys anything: CA with the PE
    # algorithm inside every subnetwork IS global PE, per the paper.
    lin = np.asarray(pfqn_clust(L, N, Z, inner='lin')[0]).flatten()
    pe = np.asarray(pfqn_clust(L, N, Z, inner='bs')[0]).flatten()
    assert not np.allclose(lin, pe, rtol=1e-9)


def test_looping_matches_the_other_ports():
    lo, up = pfqn_looping(L, N, Z)[:2]
    assert np.asarray(lo).flatten() == pytest.approx([1.341502, 0.987254], rel=1e-6)
    assert np.asarray(up).flatten() == pytest.approx([2.129780, 1.455485], rel=1e-6)


def test_looping_brackets_the_exact_solution():
    # What makes it a bound rather than an estimate.
    lo, up = pfqn_looping(L, N, Z)[:2]
    Xex = _exact()
    assert np.all(np.asarray(lo).flatten() <= Xex + 1e-9)
    assert np.all(np.asarray(up).flatten() >= Xex - 1e-9)
