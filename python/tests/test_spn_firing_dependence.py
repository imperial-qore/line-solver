"""
Native-Python regression test: marking-dependent transition firing rates
(Transition.setFiringRateDependence).

A closed single-class SPN with N tokens on two places and mass-action rates on
both transitions (rate = k * tokens-at-input-place, single server) makes each
token an independent two-state CTMC, so the stationary marking of P1 is
Binomial(N, k2/(k1+k2)) and E[n1] = N*k2/(k1+k2). SolverCTMC evaluates the
g(marking) multiplier per enumerated state, so it must match this closed form
exactly; the JSON round-trip must preserve the handle; the setter must reject
non-timed / non-exponential modes; and SolverSSA must reject the feature rather
than silently simulate the nominal (unscaled) rate.

The LDES exactness of the same model (QN ~ [3, 2]) is validated separately at
the engine level (see _kb/09-ldes-and-cache.md).
"""
import os

import numpy as np
import pytest

import line_solver
assert os.path.realpath(__file__).rsplit('/python/', 1)[0] in os.path.realpath(
    line_solver.__file__), (
    "test must import the worktree line_solver, got %s" % line_solver.__file__)

from line_solver import (Network, Place, Transition, ClosedClass, Exp, Erlang,
                         SolverCTMC, SolverSSA, SolverOptions, TimingStrategy)
from line_solver.io.linemodel_io import save_model, load_model

N, K1, K2 = 5, 2, 3


def _mass_action_loop(N=N, k1=K1, k2=K2):
    m = Network('firingdep')
    P1 = Place(m, 'P1'); P2 = Place(m, 'P2')
    T1 = Transition(m, 'T1'); T2 = Transition(m, 'T2')
    c = ClosedClass(m, 'C', N, P1, 0)
    i1 = P1.get_index0(); i2 = P2.get_index0()
    a = T1.addMode('m1'); T1.setDistribution(a, Exp(k1))
    T1.setEnablingConditions(a, c, P1, 1); T1.setFiringOutcome(a, c, P2, 1)
    T1.setFiringRateDependence(a, lambda mk: mk[i1, 0])
    b = T2.addMode('m2'); T2.setDistribution(b, Exp(k2))
    T2.setEnablingConditions(b, c, P2, 1); T2.setFiringOutcome(b, c, P1, 1)
    T2.setFiringRateDependence(b, lambda mk: mk[i2, 0])
    rm = m.initRoutingMatrix()
    rm.set(c, c, P1, T1, 1.0); rm.set(c, c, P2, T2, 1.0)
    rm.set(c, c, T1, P2, 1.0); rm.set(c, c, T2, P1, 1.0)
    m.link(rm)
    P1.setState(np.array([N])); P2.setState(np.array([0]))
    return m


def _opt():
    o = SolverOptions(); o.cutoff = N; o.verbose = False
    return o


def test_ctmc_matches_binomial_closed_form():
    avg = SolverCTMC(_mass_action_loop(), _opt()).getAvgTable()
    q = list(avg['QLen'])
    assert q[0] == pytest.approx(N * K2 / (K1 + K2), abs=1e-6)  # E[n1]=3
    assert q[1] == pytest.approx(N * K1 / (K1 + K2), abs=1e-6)  # E[n2]=2


def test_json_round_trip_preserves_dependence(tmp_path):
    fn = str(tmp_path / 'firingdep.json')
    save_model(_mass_action_loop(), fn)
    avg = SolverCTMC(load_model(fn), _opt()).getAvgTable()
    q = list(avg['QLen'])
    assert q[0] == pytest.approx(N * K2 / (K1 + K2), abs=1e-6)
    assert q[1] == pytest.approx(N * K1 / (K1 + K2), abs=1e-6)


def test_setter_rejects_non_exponential():
    m = Network('g'); P = Place(m, 'P'); T = Transition(m, 'T')
    ClosedClass(m, 'C', 3, P, 0)
    mo = T.addMode('m'); T.setDistribution(mo, Erlang(1, 2))
    with pytest.raises(ValueError):
        T.setFiringRateDependence(mo, lambda x: 1.0)


def test_setter_rejects_immediate():
    m = Network('g'); P = Place(m, 'P'); T = Transition(m, 'T')
    ClosedClass(m, 'C', 3, P, 0)
    mo = T.addMode('m'); T.setTimingStrategy(mo, TimingStrategy.IMMEDIATE)
    with pytest.raises(ValueError):
        T.setFiringRateDependence(mo, lambda x: 1.0)


def test_ssa_rejects_dependence():
    o = SolverOptions(); o.seed = 42; o.samples = 3000; o.verbose = False
    with pytest.raises(RuntimeError):
        SolverSSA(_mass_action_loop(), o).getAvgTable()


if __name__ == '__main__':
    test_ctmc_matches_binomial_closed_form()
    test_setter_rejects_non_exponential()
    test_setter_rejects_immediate()
    test_ssa_rejects_dependence()
    print('PASS')
