"""
Native-Python tests for the NRM stochastic-Petri-net path
(api/solvers/ssa/nrm.py::_solver_ssa_nrm_spn). Mirror of the JAR
``SolverSSANrmSpnTest``.

Each test asserts (a) the solver actually ran the NRM method (result.method ==
'nrm', never a silent serial fallback) and (b) the simulated marking means and
transition throughputs agree with the exact SPN CTMC.

Two net classes are covered:

* a closed net with an inhibitor arc (spn_inhibiting), validated directly
  against the SPN CTMC that handles it;
* a closed net with an IMMEDIATE transition (vanishing markings). The SPN CTMC
  path has a known open gap on immediate transitions, so ground truth is an
  EQUIVALENT all-timed net: the immediate transition merely relays a token from
  a vanishing place, so the reduced net (that place and its immediate transition
  removed) has identical steady-state marking on every tangible place. The test
  asserts the vanishing place holds exactly zero tokens and the timed places
  match the reduced-net CTMC.
"""

import os
import numpy as np
import pytest

import line_solver
# The worktree copy of line_solver must be the one under test; a global install
# would silently validate the wrong code.
assert os.path.realpath(__file__).rsplit('/python/', 1)[0] in os.path.realpath(
    line_solver.__file__), (
    "test must import the worktree line_solver, got %s" % line_solver.__file__)

from line_solver import (Network, Place, Transition, ClosedClass, Exp, Immediate,
                         SolverSSA, SolverCTMC, TimingStrategy)

SAMPLES = 300000
SEED = 23000
RTOL = 0.03


def _inhibiting():
    m = Network('spn_inhibiting')
    P1 = Place(m, 'P1'); P2 = Place(m, 'P2'); P3 = Place(m, 'P3')
    T1 = Transition(m, 'T1'); T2 = Transition(m, 'T2'); T3 = Transition(m, 'T3')
    jc = ClosedClass(m, 'Class1', 4, P1, 0)
    a = T1.addMode('Mode1'); T1.setDistribution(a, Exp(2)); T1.setEnablingConditions(a, jc, P1, 2); T1.setFiringOutcome(a, jc, P2, 2)
    b = T1.addMode('Mode2'); T1.setDistribution(b, Exp(1)); T1.setEnablingConditions(b, jc, P1, 1); T1.setFiringOutcome(b, jc, P3, 1)
    c = T2.addMode('Mode3'); T2.setDistribution(c, Exp(4)); T2.setEnablingConditions(c, jc, P2, 1); T2.setFiringOutcome(c, jc, P1, 1)
    d = T3.addMode('Mode4'); T3.setDistribution(d, Exp(1)); T3.setEnablingConditions(d, jc, P3, 3)
    T3.setInhibitingConditions(d, jc, P2, 1); T3.setFiringOutcome(d, jc, P1, 3)
    rm = m.initRoutingMatrix()
    rm.set(jc, jc, P1, T1, 1.0); rm.set(jc, jc, P2, T2, 1.0); rm.set(jc, jc, P2, T3, 1.0); rm.set(jc, jc, P3, T3, 1.0)
    rm.set(jc, jc, T1, P2, 1.0); rm.set(jc, jc, T1, P3, 1.0); rm.set(jc, jc, T2, P1, 1.0); rm.set(jc, jc, T3, P1, 1.0)
    m.link(rm)
    P1.setState(4); P2.setState(0); P3.setState(0)
    return m


def _immediate():
    m = Network('spn_immediate')
    P1 = Place(m, 'P1'); P2 = Place(m, 'P2'); P3 = Place(m, 'P3')
    T1 = Transition(m, 'T1'); T2 = Transition(m, 'T2'); T3 = Transition(m, 'T3')
    jc = ClosedClass(m, 'Class1', 2, P1, 0)
    a = T1.addMode('M1'); T1.setDistribution(a, Exp(3)); T1.setEnablingConditions(a, jc, P1, 1); T1.setFiringOutcome(a, jc, P2, 1)
    b = T2.addMode('M2'); T2.setDistribution(b, Immediate()); T2.setTimingStrategy(b, TimingStrategy.IMMEDIATE)
    T2.setFiringPriorities(b, 1); T2.setFiringWeights(b, 1.0)
    T2.setEnablingConditions(b, jc, P2, 1); T2.setFiringOutcome(b, jc, P3, 1)
    c = T3.addMode('M3'); T3.setDistribution(c, Exp(2)); T3.setEnablingConditions(c, jc, P3, 1); T3.setFiringOutcome(c, jc, P1, 1)
    rm = m.initRoutingMatrix()
    rm.set(jc, jc, P1, T1, 1.0); rm.set(jc, jc, P2, T2, 1.0); rm.set(jc, jc, P3, T3, 1.0)
    rm.set(jc, jc, T1, P2, 1.0); rm.set(jc, jc, T2, P3, 1.0); rm.set(jc, jc, T3, P1, 1.0)
    m.link(rm)
    P1.setState(2); P2.setState(0); P3.setState(0)
    return m


def _immediate_reduced():
    m = Network('spn_immediate_reduced')
    P1 = Place(m, 'P1'); P3 = Place(m, 'P3')
    T1 = Transition(m, 'T1'); T3 = Transition(m, 'T3')
    jc = ClosedClass(m, 'Class1', 2, P1, 0)
    a = T1.addMode('M1'); T1.setDistribution(a, Exp(3)); T1.setEnablingConditions(a, jc, P1, 1); T1.setFiringOutcome(a, jc, P3, 1)
    c = T3.addMode('M3'); T3.setDistribution(c, Exp(2)); T3.setEnablingConditions(c, jc, P3, 1); T3.setFiringOutcome(c, jc, P1, 1)
    rm = m.initRoutingMatrix()
    rm.set(jc, jc, P1, T1, 1.0); rm.set(jc, jc, P3, T3, 1.0)
    rm.set(jc, jc, T1, P3, 1.0); rm.set(jc, jc, T3, P1, 1.0)
    m.link(rm)
    P1.setState(2); P3.setState(0)
    return m


def _nrm(model):
    solver = SolverSSA(model, 'nrm', seed=SEED, samples=SAMPLES, verbose=False)
    q = np.asarray(solver.getAvgQLen()).flatten()
    t = np.asarray(solver.getAvgTput()).flatten()
    ran = getattr(solver._result, 'method', None)
    assert ran == 'nrm', 'SPN model must run the NRM method, got %r' % ran
    return q, t


def test_inhibiting_vs_ctmc():
    ctmc = SolverCTMC(_inhibiting(), cutoff=10, seed=1)
    qc = np.asarray(ctmc.getAvgQLen()).flatten()
    tc = np.asarray(ctmc.getAvgTput()).flatten()
    q, t = _nrm(_inhibiting())
    for i in range(qc.shape[0]):
        assert abs(q[i] - qc[i]) / max(abs(qc[i]), 1e-9) < RTOL, \
            'QLen[%d]: NRM %g vs CTMC %g' % (i, q[i], qc[i])
        assert abs(t[i] - tc[i]) / max(abs(tc[i]), 1e-9) < RTOL, \
            'Tput[%d]: NRM %g vs CTMC %g' % (i, t[i], tc[i])


def test_immediate_vanishing_marking():
    red = SolverCTMC(_immediate_reduced(), cutoff=10, seed=1)
    qr = np.asarray(red.getAvgQLen()).flatten()   # P1(0), P3(1)
    tr = np.asarray(red.getAvgTput()).flatten()
    q, t = _nrm(_immediate())                     # P1(0), P2(1), P3(2)
    # The vanishing place holds exactly zero tokens.
    assert abs(q[1]) < 1e-12, 'vanishing place P2 mean tokens %g' % q[1]
    assert abs(q[0] - qr[0]) / max(abs(qr[0]), 1e-9) < RTOL, 'QLen P1 %g vs %g' % (q[0], qr[0])
    assert abs(q[2] - qr[1]) / max(abs(qr[1]), 1e-9) < RTOL, 'QLen P3 %g vs %g' % (q[2], qr[1])
    assert abs(t[0] - tr[0]) / max(abs(tr[0]), 1e-9) < RTOL, 'Tput T1 %g vs %g' % (t[0], tr[0])
    assert abs(t[2] - tr[0]) / max(abs(tr[0]), 1e-9) < RTOL, 'Tput T3 %g vs %g' % (t[2], tr[0])
