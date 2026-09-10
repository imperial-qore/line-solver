"""Native-Python regression test: SolverFLD solves stochastic Petri nets on the
'dae' method.

The properties asserted here are the ones that distinguish this route from a
plausible-looking wrong answer:

  * WHERE THE DRIFT IS LINEAR THE FLUID ANSWER IS EXACT, not approximate. A
    single-input mode with an infinite server count makes min() exact on its
    whole support, so an open net Source -> P1 -> T1(inf) -> Sink must return
    lambda/mu to solver tolerance, not near it.
  * CONSERVATION IS AN EQUATION, not a consequence of the drift: a closed net's
    P-invariants hold to the Newton tolerance rather than the integrator's.
  * AN IMMEDIATE TRANSITION IS AN ALGEBRAIC FLOW. Its input place is pinned at
    zero and the marking splits as the exact CTMC says it does.
  * A MULTI-PHASE FIRING TIME carries its running servers as counts latched to
    the enabling degree, so an Erlang(2) net tracks the exact chain.
  * every OTHER fluid method refuses the model rather than integrating an empty
    one and reporting zeros.
"""

import os
import sys

_WORKTREE_PY = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
if sys.path and sys.path[0] != _WORKTREE_PY:
    sys.path.insert(0, _WORKTREE_PY)

import numpy as np
import pytest

from line_solver import (Network, Source, Sink, Place, Transition, OpenClass,
                         ClosedClass, Exp, Erlang, TimingStrategy, FLD, CTMC)

import line_solver as _ls
assert os.path.abspath(_ls.__file__).startswith(_WORKTREE_PY), (
    'wrong line_solver imported: %s' % _ls.__file__)


def _open_net(lam=1.0, mu=4.0):
    m = Network('spn_basic_open')
    src, snk = Source(m, 'Source'), Sink(m, 'Sink')
    P1, T1 = Place(m, 'P1'), Transition(m, 'T1')
    jc = OpenClass(m, 'Class1')
    src.setArrival(jc, Exp(lam))
    a = T1.addMode('Mode1')
    T1.setNumberOfServers(a, float('inf'))
    T1.setDistribution(a, Exp(mu))
    T1.setEnablingConditions(a, jc, P1, 1)
    T1.setFiringOutcome(a, jc, snk, 1)
    R = m.initRoutingMatrix()
    R.set(jc, jc, src, P1, 1.0)
    R.set(jc, jc, P1, T1, 1.0)
    R.set(jc, jc, T1, snk, 1.0)
    m.link(R)
    return m


def _cycle(N, c):
    """P1 -> T1(c servers) -> P2 -> T2(inf) -> P1, closed with N tokens."""
    m = Network('cyc')
    P1, P2 = Place(m, 'P1'), Place(m, 'P2')
    T1, T2 = Transition(m, 'T1'), Transition(m, 'T2')
    jc = ClosedClass(m, 'C', N, P1)
    a = T1.addMode('a')
    T1.setNumberOfServers(a, c)
    T1.setDistribution(a, Exp(1.0))
    T1.setEnablingConditions(a, jc, P1, 1)
    T1.setFiringOutcome(a, jc, P2, 1)
    b = T2.addMode('b')
    T2.setNumberOfServers(b, float('inf'))
    T2.setDistribution(b, Exp(1.0))
    T2.setEnablingConditions(b, jc, P2, 1)
    T2.setFiringOutcome(b, jc, P1, 1)
    R = m.initRoutingMatrix()
    R.set(jc, jc, P1, T1, 1.0)
    R.set(jc, jc, T1, P2, 1.0)
    R.set(jc, jc, P2, T2, 1.0)
    R.set(jc, jc, T2, P1, 1.0)
    m.link(R)
    return m


def _immediate_net(N):
    """P1 -(timed)-> P2 -(IMMEDIATE)-> P3 -(timed)-> P1, closed with N tokens."""
    m = Network('imm')
    P = [Place(m, 'P%d' % (i + 1)) for i in range(3)]
    T = [Transition(m, 'T%d' % (i + 1)) for i in range(3)]
    jc = ClosedClass(m, 'C', N, P[0])
    for i, (rate, timing) in enumerate([(1.0, None), (None, TimingStrategy.IMMEDIATE), (2.0, None)]):
        mode = T[i].addMode('m%d' % i)
        T[i].setNumberOfServers(mode, float('inf'))
        if timing is None:
            T[i].setDistribution(mode, Exp(rate))
        else:
            T[i].setTimingStrategy(mode, timing)
        T[i].setEnablingConditions(mode, jc, P[i], 1)
        T[i].setFiringOutcome(mode, jc, P[(i + 1) % 3], 1)
    R = m.initRoutingMatrix()
    for i in range(3):
        R.set(jc, jc, P[i], T[i], 1.0)
        R.set(jc, jc, T[i], P[(i + 1) % 3], 1.0)
    m.link(R)
    return m


def _erlang_net(N):
    m = Network('erl')
    P1, P2 = Place(m, 'P1'), Place(m, 'P2')
    T1, T2 = Transition(m, 'T1'), Transition(m, 'T2')
    jc = ClosedClass(m, 'C', N, P1)
    a = T1.addMode('a')
    T1.setNumberOfServers(a, 1)
    T1.setDistribution(a, Erlang.fitMeanAndOrder(1.0, 2))
    T1.setEnablingConditions(a, jc, P1, 1)
    T1.setFiringOutcome(a, jc, P2, 1)
    b = T2.addMode('b')
    T2.setNumberOfServers(b, float('inf'))
    T2.setDistribution(b, Exp(1.0))
    T2.setEnablingConditions(b, jc, P2, 1)
    T2.setFiringOutcome(b, jc, P1, 1)
    R = m.initRoutingMatrix()
    R.set(jc, jc, P1, T1, 1.0)
    R.set(jc, jc, T1, P2, 1.0)
    R.set(jc, jc, P2, T2, 1.0)
    R.set(jc, jc, T2, P1, 1.0)
    m.link(R)
    return m


def _qlen(model, method='dae'):
    return np.asarray(FLD(model, method).getAvgQLen()).reshape(-1)


def test_linear_open_net_is_exact():
    """min() exact on its whole support => the fluid mean is the exact mean."""
    q = _qlen(_open_net(lam=1.0, mu=4.0))
    assert q[1] == pytest.approx(0.25, abs=1e-9)


def test_default_method_resolves_to_dae():
    """A Petri net has one fluid route; 'default' must find it."""
    q = np.asarray(FLD(_cycle(4, 2)).getAvgQLen()).reshape(-1)
    assert q.sum() == pytest.approx(4.0, abs=1e-7)


@pytest.mark.parametrize('N,c', [(4, 2), (8, 4), (20, 10)])
def test_conservation_is_an_equation(N, c):
    """The P-invariant holds to the NEWTON tolerance, not the integrator's."""
    q = _qlen(_cycle(N, c))
    assert q.sum() == pytest.approx(float(N), abs=1e-7)


def test_matches_the_matlab_reference():
    """The same closure the MATLAB twin solves, on the same nets.

    Recorded from matlab/src/solvers/FLD/solver_fluid_petri.m; see
    _kb/06-solver-catalog.md for the measurement.
    """
    for (N, c, ref) in [(4, 2, 2.30470), (8, 4, 4.44146), (20, 10, 10.71316)]:
        q = _qlen(_cycle(N, c))
        assert q[0] == pytest.approx(ref, abs=1e-4)


def test_immediate_transition_pins_its_input():
    """An immediate mode's input place holds no mass, and the rest splits as the
    exact chain says."""
    q = _qlen(_immediate_net(4))
    exact = np.asarray(CTMC(_immediate_net(4)).getAvgQLen()).reshape(-1)
    assert q[1] == pytest.approx(0.0, abs=1e-9)
    assert q == pytest.approx(exact, rel=1e-3)
    assert q.sum() == pytest.approx(4.0, abs=1e-7)


@pytest.mark.parametrize('N', [3, 4])
def test_multiphase_firing_tracks_the_exact_chain(N):
    """Erlang(2) firing: the running servers are counts latched to the enabling
    degree, so the marking tracks the exact chain to closure accuracy."""
    q = _qlen(_erlang_net(N))
    exact = np.asarray(CTMC(_erlang_net(N)).getAvgQLen()).reshape(-1)
    assert q == pytest.approx(exact, rel=5e-3)
    assert q.sum() == pytest.approx(float(N), abs=1e-7)


def test_every_other_fluid_method_refuses_a_petri_net():
    """Integrating a Petri net as a queueing network reports zeros with no
    warning, so the gate must refuse it by name."""
    with pytest.raises(Exception) as exc:
        FLD(_cycle(4, 2), 'closing').runAnalyzer()
    msg = str(exc.value).lower()
    assert 'transition' in msg or 'place' in msg or 'petri' in msg
