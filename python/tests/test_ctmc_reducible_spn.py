"""Reducible SPN with a transient SCC feeding two disjoint recurrent classes.

Three tokens sit in P1. T1 (rate 1) moves all three into cycle A (P2 <-> P4),
T2 (rate 2) moves all three into cycle B (P3 <-> P5). Nothing ever returns a
token to P1, so the choice is irreversible and the chain decomposes into one
transient SCC (the initial marking) and two BSCCs.

The reference is CLOSED FORM, not a recorded output: the branch is the T1/T2
race, P(A) = 1/3 and P(B) = 2/3, and each cycle is a closed two-place
exponential cycle whose marginal is geometric in the rate ratio. An independent
external GSPN tool reproduces these values to five decimals.

This model is the regression for three defects that were each silent:
  - a zero-visit station had its capacity zeroed, deleting the declared initial
    marking from the state space, so the transient SCC never existed;
  - the SCC support graph was built with a sign test, which shatters a
    matrix-exponential generator;
  - a chain reducible within a SINGLE weak component was handed to the plain
    solve as a singular system, which returned NaN.

See _kb/11-conventions-and-gotchas.md.
"""
import numpy as np
import pytest
from scipy.sparse import csc_matrix
from scipy.sparse.csgraph import connected_components

from line_solver import (Network, Place, Transition, ClosedClass, Exp,
                         SolverCTMC, GlobalConstants)

# exact: P(A) = 1/3, P(B) = 2/3; cycle A marginal ratio 4/3, cycle B ratio 6/5
EXACT = {'P2': 108.0 / 175.0, 'P3': 772.0 / 671.0,
         'P4': 67.0 / 175.0, 'P5': 570.0 / 671.0}
TOL = 1e-6


def _model():
    m = Network('spn_twobscc')
    P1 = Place(m, 'P1'); P2 = Place(m, 'P2'); P3 = Place(m, 'P3')
    P4 = Place(m, 'P4'); P5 = Place(m, 'P5')
    T1 = Transition(m, 'T1'); T2 = Transition(m, 'T2')
    TA = Transition(m, 'TA'); TA2 = Transition(m, 'TA2')
    TB = Transition(m, 'TB'); TB2 = Transition(m, 'TB2')
    jc = ClosedClass(m, 'C', 3, P1, 0)

    def mode(T, name, rate, pin, nin, pout, nout):
        md = T.addMode(name)
        T.setDistribution(md, Exp(rate))
        T.setEnablingConditions(md, jc, pin, nin)
        T.setFiringOutcome(md, jc, pout, nout)

    mode(T1, 'mA', 1.0, P1, 3, P2, 3)
    mode(T2, 'mB', 2.0, P1, 3, P3, 3)
    mode(TA, 'a1', 3.0, P2, 1, P4, 1)
    mode(TA2, 'a2', 4.0, P4, 1, P2, 1)
    mode(TB, 'b1', 5.0, P3, 1, P5, 1)
    mode(TB2, 'b2', 6.0, P5, 1, P3, 1)

    R = m.initRoutingMatrix()
    R.set(jc, jc, P1, T1, 0.5); R.set(jc, jc, P1, T2, 0.5)
    R.set(jc, jc, T1, P2, 1.0); R.set(jc, jc, T2, P3, 1.0)
    R.set(jc, jc, P2, TA, 1.0); R.set(jc, jc, TA, P4, 1.0)
    R.set(jc, jc, P4, TA2, 1.0); R.set(jc, jc, TA2, P2, 1.0)
    R.set(jc, jc, P3, TB, 1.0); R.set(jc, jc, TB, P5, 1.0)
    R.set(jc, jc, P5, TB2, 1.0); R.set(jc, jc, TB2, P3, 1.0)
    m.link(R)
    return m


def test_ctmc_matches_the_closed_form_mixture():
    """Means must be the absorption-weighted mixture of the two BSCCs."""
    table = SolverCTMC(_model(), cutoff=3).getAvgTable()
    got = {str(row.Station): float(row.QLen) for row in table.itertuples()}
    for place, expected in EXACT.items():
        assert place in got, 'place %s absent from the AvgTable' % place
        assert got[place] == pytest.approx(expected, abs=TOL), (
            '%s: got %.9f, exact %.9f' % (place, got[place], expected))


def test_the_generator_keeps_its_transient_scc():
    """The structure is the point: a transient SCC feeding disjoint BSCCs.
    Collapsing it is what made the means wrong before."""
    solver = SolverCTMC(_model(), cutoff=3, keep=True)
    solver.getAvgTable()
    Q = np.atleast_2d(np.asarray(solver.getInfGen(), dtype=float))
    adj = Q.copy()
    np.fill_diagonal(adj, 0.0)
    # arc by MAGNITUDE, never by sign
    mask = np.abs(adj) > GlobalConstants.ArcTol
    nscc, labels = connected_components(csc_matrix(mask), directed=True,
                                        connection='strong', return_labels=True)
    nbscc = 0
    for c in range(nscc):
        members = np.where(labels == c)[0]
        others = np.setdiff1d(np.arange(Q.shape[0]), members)
        if others.size == 0 or not np.any(mask[np.ix_(members, others)]):
            nbscc += 1
    # counts differ across codebases because enumeration and stochastic
    # complementation retain different numbers of states (MATLAB 9, python 35);
    # the INVARIANT is what matters and holds in both
    assert nscc > 1, 'the chain must stay reducible, got a single SCC'
    assert nbscc >= 2, 'expected at least two recurrent classes, got %d' % nbscc
    assert nscc - nbscc == 1, (
        'expected exactly one transient SCC (the initial marking), got %d; '
        'zero means the declared marking was dropped from the state space'
        % (nscc - nbscc))


def test_mass_splits_by_the_firing_race_not_uniformly():
    """P(A)=1/3 comes from the T1/T2 rates. A uniform 50/50 fallback -- what an
    unseeded decomposition returns -- must not pass."""
    table = SolverCTMC(_model(), cutoff=3).getAvgTable()
    got = {str(row.Station): float(row.QLen) for row in table.itertuples()}
    mass_a = got['P2'] + got['P4']
    mass_b = got['P3'] + got['P5']
    assert mass_a == pytest.approx(1.0, abs=TOL), 'cycle A mass %.9f' % mass_a
    assert mass_b == pytest.approx(2.0, abs=TOL), 'cycle B mass %.9f' % mass_b
