"""
Native-Python tests for the NRM stochastic-Petri-net path on OPEN nets
(api/solvers/ssa/nrm.py::_solver_ssa_nrm_spn). Mirror of the JAR
``SolverSSANrmOpenSpnTest``.

An open SPN has a Source feeding a Place whose tokens drain through a Transition
to a Sink. Before the Source-arrival reaction was added the fed Place stayed
empty and the run raised "Deadlock: no transition is enabled".

Each test asserts (a) the solver actually ran the NRM method (result.method ==
'nrm', never a silent serial fallback) and (b) the simulated marking means and
throughputs match the analytic M/M/1 result: a Source Exp(lambda) feeding a
single-server Transition Exp(mu) is an M/M/1 queue at the Place, with mean
tokens rho/(1-rho) and throughput lambda (rho = lambda/mu < 1).
"""

import os
import numpy as np

import line_solver
# The worktree copy of line_solver must be the one under test; a global install
# would silently validate the wrong code.
assert os.path.realpath(__file__).rsplit('/python/', 1)[0] in os.path.realpath(
    line_solver.__file__), (
    "test must import the worktree line_solver, got %s" % line_solver.__file__)

from line_solver import (Network, Source, Sink, Place, Transition, OpenClass,
                         Exp, SolverSSA)

SAMPLES = 300000
SEED = 23000
RTOL = 0.04


def _mm1(lam, mu):
    """Open M/M/1-as-SPN: Source Exp(lam) -> P1 -> T1 Exp(mu) -> Sink."""
    m = Network('mm1spn')
    source = Source(m, 'Source'); sink = Sink(m, 'Sink')
    P1 = Place(m, 'P1'); T1 = Transition(m, 'T1')
    jc = OpenClass(m, 'Class1', 0)
    source.setArrival(jc, Exp(lam))
    a = T1.addMode('Mode1'); T1.setDistribution(a, Exp(mu))
    T1.setEnablingConditions(a, jc, P1, 1); T1.setFiringOutcome(a, jc, sink, 1)
    rm = m.initRoutingMatrix()
    rm.set(jc, jc, source, P1, 1.0); rm.set(jc, jc, P1, T1, 1.0); rm.set(jc, jc, T1, sink, 1.0)
    m.link(rm)
    return m


def _tandem(lam, mu1, mu2):
    """Open tandem: Source -> P1 -> T1 Exp(mu1) -> P2 -> T2 Exp(mu2) -> Sink."""
    m = Network('tandemspn')
    source = Source(m, 'Source'); sink = Sink(m, 'Sink')
    P1 = Place(m, 'P1'); P2 = Place(m, 'P2')
    T1 = Transition(m, 'T1'); T2 = Transition(m, 'T2')
    jc = OpenClass(m, 'Class1', 0)
    source.setArrival(jc, Exp(lam))
    a = T1.addMode('Mode1'); T1.setDistribution(a, Exp(mu1))
    T1.setEnablingConditions(a, jc, P1, 1); T1.setFiringOutcome(a, jc, P2, 1)
    b = T2.addMode('Mode1'); T2.setDistribution(b, Exp(mu2))
    T2.setEnablingConditions(b, jc, P2, 1); T2.setFiringOutcome(b, jc, sink, 1)
    rm = m.initRoutingMatrix()
    rm.set(jc, jc, source, P1, 1.0); rm.set(jc, jc, P1, T1, 1.0)
    rm.set(jc, jc, T1, P2, 1.0); rm.set(jc, jc, P2, T2, 1.0); rm.set(jc, jc, T2, sink, 1.0)
    m.link(rm)
    return m


def _nrm(model):
    solver = SolverSSA(model, 'nrm', seed=SEED, samples=SAMPLES, verbose=False)
    q = np.asarray(solver.getAvgQLen()).flatten()
    t = np.asarray(solver.getAvgTput()).flatten()
    ran = getattr(solver._result, 'method', None)
    assert ran == 'nrm', 'open SPN must run the NRM method, got %r' % ran
    return q, t


def _close(expected, got, what, rtol=RTOL):
    assert abs(got - expected) / max(abs(expected), 1e-9) < rtol, \
        '%s: NRM %g deviates from %g' % (what, got, expected)


def test_mm1_open_spn_rho05():
    lam, mu = 0.5, 1.0
    q_exact = (lam / mu) / (1.0 - lam / mu)   # = 1.0
    q, t = _nrm(_mm1(lam, mu))                 # stations: Source(0), P1(1)
    _close(q_exact, q[1], 'P1 mean tokens')
    _close(lam, t[0], 'Source throughput')
    _close(lam, t[1], 'P1 (T1) throughput')


def test_mm1_open_spn_rho08():
    lam, mu = 0.8, 1.0
    q_exact = (lam / mu) / (1.0 - lam / mu)   # = 4.0
    q, t = _nrm(_mm1(lam, mu))
    _close(q_exact, q[1], 'P1 mean tokens')
    _close(lam, t[0], 'Source throughput')
    _close(lam, t[1], 'P1 (T1) throughput')


def test_open_tandem_two_places():
    lam, mu1, mu2 = 0.5, 1.0, 2.0
    q1 = (lam / mu1) / (1.0 - lam / mu1)      # = 1.0
    q2 = (lam / mu2) / (1.0 - lam / mu2)      # = 1/3
    q, t = _nrm(_tandem(lam, mu1, mu2))        # stations: Source(0), P1(1), P2(2)
    _close(q1, q[1], 'P1 mean tokens')
    _close(q2, q[2], 'P2 mean tokens')
    _close(lam, t[1], 'P1 (T1) throughput')
    _close(lam, t[2], 'P2 (T2) throughput')
