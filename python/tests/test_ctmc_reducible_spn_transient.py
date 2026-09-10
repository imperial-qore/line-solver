"""Transient analysis of the reducible SPN with a transient SCC.

Same model as test_ctmc_reducible_spn: three tokens in P1, T1 (rate 1) drains
them into cycle A, T2 (rate 2) into cycle B, and nothing returns to P1.

The transient has a CLOSED FORM that does not depend on the internal cycle
dynamics at all, which makes it a sharp test of how the transient solver treats
a transient SCC. P1 is left at the total race rate 1 + 2 = 3, so

    E[P1](t)               = 3 exp(-3t)
    E[P2](t) + E[P4](t)    =   (1 - exp(-3t))        (branch weight 1/3)
    E[P3](t) + E[P5](t)    = 2 (1 - exp(-3t))        (branch weight 2/3)

and the three tokens are conserved at every t. As t grows these converge to the
steady-state mixture asserted in test_ctmc_reducible_spn.

The transient SCC is exactly the state this checks: if it were dropped from the
state space (the original defect) the drain curve would not exist, and if the
branch were taken uniformly the 1:2 ratio would come out 1:1.

See _kb/11-conventions-and-gotchas.md.
"""
import numpy as np
import pytest

from line_solver import SolverCTMC

from test_ctmc_reducible_spn import _model

TOL = 1e-6
TSPAN = 4.0
RACE_RATE = 3.0          # T1 (1.0) + T2 (2.0), both enabled by the 3 tokens in P1
POP = 3.0


def _curves():
    """Per-place transient mean curves on a common time grid."""
    solver = SolverCTMC(_model(), cutoff=3, timespan=[0.0, TSPAN])
    QNt, _, _ = solver.getTranAvg()
    by_place = {}
    for i, station in enumerate(solver.model.get_stations()):
        row = QNt[i]
        entry = row[0] if isinstance(row, (list, tuple)) else row
        if entry is None:
            continue
        t = np.asarray(entry.t, dtype=float).ravel()
        m = np.asarray(entry.metric, dtype=float).ravel()
        by_place[station.getName()] = (t, m)
    return by_place


def _at(curve, tq):
    t, m = curve
    return np.interp(tq, t, m)


def test_p1_drains_at_the_race_rate():
    """The transient SCC empties as a pure exponential at rate T1+T2."""
    curves = _curves()
    assert 'P1' in curves, (
        'P1 has no transient curve: the transient SCC is absent from the state space')
    tq = np.linspace(0.0, TSPAN, 21)
    got = _at(curves['P1'], tq)
    exact = POP * np.exp(-RACE_RATE * tq)
    err = float(np.max(np.abs(got - exact)))
    assert err < 1e-4, 'E[P1](t) deviates from 3exp(-3t) by %.3e' % err


def test_drained_mass_splits_one_to_two_at_every_time():
    """The 1:2 split is set by the T1/T2 rates and holds at every t, not just in
    the limit. A uniform branch would give 1:1."""
    curves = _curves()
    tq = np.linspace(0.2, TSPAN, 20)
    mass_a = _at(curves['P2'], tq) + _at(curves['P4'], tq)
    mass_b = _at(curves['P3'], tq) + _at(curves['P5'], tq)
    drained = POP * (1.0 - np.exp(-RACE_RATE * tq))
    assert np.max(np.abs(mass_a - drained / 3.0)) < 1e-4, 'cycle A fill curve'
    assert np.max(np.abs(mass_b - 2.0 * drained / 3.0)) < 1e-4, 'cycle B fill curve'


def test_tokens_are_conserved_at_every_time():
    """Nothing is created or destroyed: the transient states must carry their
    mass rather than leak it."""
    curves = _curves()
    tq = np.linspace(0.0, TSPAN, 21)
    total = sum(_at(curves[p], tq) for p in ('P1', 'P2', 'P3', 'P4', 'P5'))
    assert np.max(np.abs(total - POP)) < 1e-4, (
        'token count drifts from %g: max deviation %.3e'
        % (POP, float(np.max(np.abs(total - POP)))))


def test_the_curves_are_not_constant():
    """Guard against the silent fallback: when the initial state cannot be located
    the solver returns the STATIONARY vector as a flat curve, which satisfies
    conservation and the steady-state limit vacuously. A real transient varies."""
    curves = _curves()
    spans = {p: float(np.ptp(np.asarray(curves[p][1], dtype=float)))
             for p in ('P1', 'P2', 'P3', 'P4', 'P5')}
    assert max(spans.values()) > 1e-3, (
        'every curve is flat (%s): the transient collapsed to the stationary '
        'distribution, so no transient analysis actually ran' % spans)


def test_transient_converges_to_the_steady_state_mixture():
    """The long-run limit of the transient must be the stationary answer."""
    curves = _curves()
    late = TSPAN
    got = {p: float(_at(curves[p], late)) for p in ('P2', 'P3', 'P4', 'P5')}
    exact = {'P2': 108.0 / 175.0, 'P3': 772.0 / 671.0,
             'P4': 67.0 / 175.0, 'P5': 570.0 / 671.0}
    for place, want in exact.items():
        assert got[place] == pytest.approx(want, abs=1e-3), (
            '%s at t=%g: got %.6f, steady state %.6f' % (place, late, got[place], want))
