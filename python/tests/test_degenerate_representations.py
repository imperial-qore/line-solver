"""Degenerate-representation checks for the Markovian distribution family.

A general ME or RAP has no independent reference to test against, so a wrong
accessor can sit undetected for a long time. A DEGENERATE representation does
have one: an ME built from an Erlang must answer exactly like that Erlang, and a
MAP whose D1 restarts the phase process from a fixed vector is a renewal process
that must answer like the corresponding HyperExp. That technique is what exposed
the MAP.getMu/getPhi/getInitProb defect fixed on 2026-07-20, so it is what these
tests use.

Reference values were produced by running MATLAB, which passes every case here.
Compare ProcessType members BY NAME, never by ordinal.
"""

import numpy as np
import pytest

from line_solver import (Network, Source, Sink, Place, Transition, OpenClass,
                         Exp, Erlang, HyperExp, SchedStrategy, SolverCTMC)
from line_solver.api.mam import map_pie
from line_solver.distributions.markovian import MAP, RAP, ME


# A MAP with genuine hidden phase changes: D0 carries off-diagonal mass, so the
# arrival-only quantity D1*e differs from the total exit rate -diag(D0), and the
# time-stationary phase vector differs from the arrival-embedded one.
_D0 = np.array([[-1.1, 0.1], [0.2, -10.2]])
_D1 = np.array([[1.0, 0.0], [0.0, 10.0]])

# HyperExp(p, rates) written as a renewal MAP: after every arrival the phase is
# redrawn from p, independently of the phase that just completed.
_HP = np.array([0.6, 0.4])
_HRATES = np.array([2.0, 0.5])
_HD0 = np.diag(-_HRATES)
_HD1 = np.outer(_HRATES, _HP)


def test_map_accessors_match_analytic_definitions():
    """MAP accessors equal the MATLAB Markovian formulas.

    MATLAB reference for this (D0,D1): getMu = [1.1, 10.2],
    getPhi = [0.90909091, 0.98039216], getInitProb = [0.16666667, 0.83333333].
    """
    m = MAP(_D0, _D1)
    assert np.allclose(m.getMu(), -np.diag(_D0))
    assert np.allclose(m.getMu(), np.array([1.1, 10.2]))
    assert np.allclose(m.getPhi(), -_D1.sum(axis=1) / np.diag(_D0))
    assert np.allclose(m.getPhi(), np.array([0.90909091, 0.98039216]))
    assert np.allclose(m.getInitProb(), np.ravel(map_pie(_D0, _D1)))
    assert np.allclose(m.getInitProb(), np.array([0.16666667, 0.83333333]))


def test_map_init_prob_is_not_the_stationary_vector():
    """Pins the specific confusion that caused the defect.

    getInitProb is the arrival-embedded vector, not the time-stationary phase
    distribution. For a MAP with hidden phase changes the two differ, and
    returning the stationary one silently corrupts any consumer that seeds a
    phase process.
    """
    m = MAP(_D0, _D1)
    stationary = m._pi
    assert not np.allclose(m.getInitProb(), stationary)


@pytest.mark.parametrize('reference,degenerate', [
    (Erlang(2, 2), ME.fromErlang(2, 2)),
    (HyperExp(_HP, _HRATES), ME.fromHyperExp(_HP, _HRATES)),
])
def test_me_degenerate_matches_reference_class(reference, degenerate):
    """An ME that IS a phase-type answers exactly like the phase-type.

    The tolerance absorbs a signed zero in ME.fromErlang's phi (-0.0 against
    0.0) without hiding a sign error: -1.0 against 1.0 would still fail.
    """
    assert np.allclose(degenerate.getMu(), reference.getMu(), rtol=0, atol=1e-12)
    assert np.allclose(degenerate.getPhi(), reference.getPhi(), rtol=0, atol=1e-12)
    assert np.allclose(degenerate.getInitProb(), reference.getInitProb(),
                       rtol=0, atol=1e-12)
    assert degenerate.getMean() == pytest.approx(reference.getMean(), rel=1e-12)
    assert degenerate.getSCV() == pytest.approx(reference.getSCV(), rel=1e-12)


def test_rap_degenerate_agrees_with_map_on_all_accessors():
    """A RAP whose matrices are nonnegative IS a MAP and must answer identically.

    This is the regression guard for a contradiction that lived inside one
    codebase: on this same (D0,D1) the two classes returned different mu and phi
    while agreeing on the initial vector. MATLAB gives mu = [1.1, 10.2],
    phi = [0.90909091, 0.98039216], alpha = [0.16666667, 0.83333333] for both.
    """
    rap, mp = RAP(_D0, _D1), MAP(_D0, _D1)
    assert np.allclose(rap.getMu(), mp.getMu())
    assert np.allclose(rap.getPhi(), mp.getPhi())
    assert np.allclose(rap.getInitProb(), mp.getInitProb())
    assert np.allclose(rap.getMu(), np.array([1.1, 10.2]))
    assert np.allclose(rap.getPhi(), np.array([0.90909091, 0.98039216]))
    assert np.allclose(rap.getInitProb(), np.array([0.16666667, 0.83333333]))


def test_genuine_rap_accessors_match_matlab():
    """A RAP that is NOT a MAP still answers, and answers as MATLAB does.

    H1 carries a negative entry here, so no MAP represents this process. MATLAB
    reference: mu = [5, 6], phi = [0.8, 0.83333333],
    alpha = [1.10204082, -0.10204082]. The initial vector leaves [0,1], which is
    the expected consequence of a RAP having no phase-type representation; it is
    recorded, not repaired.
    """
    rap = RAP(np.array([[-5.0, 1.0], [1.0, -6.0]]),
              np.array([[4.5, -0.5], [4.5, 0.5]]))
    assert np.allclose(rap.getMu(), np.array([5.0, 6.0]))
    assert np.allclose(rap.getPhi(), np.array([0.8, 0.83333333]))
    assert np.allclose(rap.getInitProb(),
                       np.array([1.10204082, -0.10204082]))


def test_renewal_map_reproduces_hyperexp_accessors():
    """A MAP that IS a HyperExp renewal process matches it on all accessors."""
    mp = MAP(_HD0, _HD1)
    ref = HyperExp(_HP, _HRATES)
    assert np.allclose(mp.getMu(), ref.getMu())
    assert np.allclose(mp.getPhi(), ref.getPhi())
    assert np.allclose(mp.getInitProb(), ref.getInitProb())
    assert mp.getMean() == pytest.approx(ref.getMean(), rel=1e-12)
    assert mp.getSCV() == pytest.approx(ref.getSCV(), rel=1e-12)


def _spn_with_firing(firing):
    """M/M/1-shaped SPN whose service completion fires with the given law."""
    model = Network('spn_degenerate')
    source = Source(model, 'source')
    sink = Sink(model, 'sink')
    p_queue = Place(model, 'queue')
    p_service = Place(model, 'service')
    t_begin = Transition(model, 'begin_service')
    t_finish = Transition(model, 'complete_service')
    jc = OpenClass(model, 'jobs')
    source.setArrival(jc, Exp.fit_mean(1 / 0.8))

    mode_begin = t_begin.add_mode('begin')
    t_begin.set_distribution(mode_begin, Exp.fit_mean(1.0))
    t_begin.set_enabling_conditions(mode_begin, jc, p_queue, 1)
    t_begin.set_enabling_conditions(mode_begin, jc, p_service, 0)
    t_begin.set_firing_outcome(mode_begin, jc, p_queue, -1)
    t_begin.set_firing_outcome(mode_begin, jc, p_service, 1)

    mode_finish = t_finish.add_mode('finish')
    t_finish.set_distribution(mode_finish, firing)
    t_finish.set_enabling_conditions(mode_finish, jc, p_service, 1)
    t_finish.set_firing_outcome(mode_finish, jc, p_service, -1)

    R = model.init_routing_matrix()
    R.set(jc, jc, source, p_queue, 1.0)
    R.set(jc, jc, p_queue, t_begin, 1.0)
    R.set(jc, jc, t_begin, p_service, 1.0)
    R.set(jc, jc, p_service, t_finish, 1.0)
    R.set(jc, jc, t_finish, sink, 1.0)
    model.link(R)
    return model


def test_spn_firing_pie_is_the_embedded_vector():
    """sn.nodeparam[...].firingpie seeds the firing process at its start.

    MATLAB reports [0.6, 0.4] for the HyperExp firing on this model; the
    equivalent renewal MAP must report the same. Before the MAP.getInitProb
    fix it reported the time-stationary [0.27272727, 0.72727273].
    """
    for firing in (HyperExp(_HP, _HRATES), MAP(_HD0, _HD1)):
        sn = _spn_with_firing(firing).getStruct()
        idx = [i for i, n in enumerate(sn.nodenames)
               if n == 'complete_service'][0]
        pie = np.asarray(sn.nodeparam[idx].firingpie[0], dtype=float)
        assert np.allclose(pie, _HP)


def test_spn_with_renewal_map_firing_matches_hyperexp_firing():
    """End-to-end: the degenerate MAP must not change the solved SPN.

    This is the check that would have caught the accessor defect at solver
    level; before the fix the throughput differed by roughly 13 percent.
    """
    cols = ['QLen', 'Util', 'Tput']
    ref = SolverCTMC(_spn_with_firing(HyperExp(_HP, _HRATES)),
                     cutoff=3).getAvgTable()[cols].to_numpy()
    got = SolverCTMC(_spn_with_firing(MAP(_HD0, _HD1)),
                     cutoff=3).getAvgTable()[cols].to_numpy()
    assert np.allclose(got, ref)
