"""
PNML (ISO/IEC 15909-2) place/transition import and export.

The oracle is not a recorded file. A round trip through save_pnml/load_pnml must
leave every CTMC metric of the model unchanged, which a consistent error on both
sides of the interchange would not survive; and the untimed net read from
another tool is checked against a marginal computed by hand from the chain it
defines.

Twin of jar/src/test/java/jline/io/PnmlIOTest.java and of the MATLAB checks on
pnml_save.m / pnml_load.m.
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

from line_solver import (ClosedClass, Erlang, Exp, Network, OpenClass, Place,
                         Sink, SolverCTMC, Source, TimingStrategy, Transition,
                         load_pnml, save_pnml)


def _two_modes():
    """Two places, one timed mode each, arc weight 2 on one side."""
    m = Network('twomodes')
    p1 = Place(m, 'P1')
    p2 = Place(m, 'P2')
    t1 = Transition(m, 'T1')
    t2 = Transition(m, 'T2')
    jc = ClosedClass(m, 'Class1', 4, p1, 0)

    a = t1.addMode('Mode1')
    t1.setDistribution(a, Exp(2))
    t1.setEnablingConditions(a, jc, p1, 2)
    t1.setFiringOutcome(a, jc, p2, 2)

    b = t2.addMode('Mode2')
    t2.setDistribution(b, Erlang(1.5, 2))
    t2.setEnablingConditions(b, jc, p2, 1)
    t2.setFiringOutcome(b, jc, p1, 1)

    rm = m.initRoutingMatrix()
    rm.set(jc, jc, p1, t1, 1.0)
    rm.set(jc, jc, p2, t2, 1.0)
    rm.set(jc, jc, t1, p2, 1.0)
    rm.set(jc, jc, t2, p1, 1.0)
    m.link(rm)
    p1.setState(4)
    p2.setState(0)
    return m


def _multi_mode():
    """One transition carrying TWO modes, plus an immediate mode with an
    inhibitor arc."""
    m = Network('multimode')
    q1 = Place(m, 'Q1')
    q2 = Place(m, 'Q2')
    u1 = Transition(m, 'U1')
    u2 = Transition(m, 'U2')
    jc = ClosedClass(m, 'Class1', 3, q1, 0)

    fast = u1.addMode('Fast')
    u1.setDistribution(fast, Exp(3))
    u1.setEnablingConditions(fast, jc, q1, 1)
    u1.setFiringOutcome(fast, jc, q2, 1)
    slow = u1.addMode('Slow')
    u1.setDistribution(slow, Exp(1))
    u1.setEnablingConditions(slow, jc, q1, 2)
    u1.setFiringOutcome(slow, jc, q2, 2)

    back = u2.addMode('Back')
    u2.setTimingStrategy(back, TimingStrategy.IMMEDIATE)
    u2.setFiringWeights(back, 2.5)
    u2.setFiringPriorities(back, 1)
    u2.setEnablingConditions(back, jc, q2, 1)
    u2.setInhibitingConditions(back, jc, q1, 3)
    u2.setFiringOutcome(back, jc, q1, 1)

    rm = m.initRoutingMatrix()
    rm.set(jc, jc, q1, u1, 1.0)
    rm.set(jc, jc, u1, q2, 1.0)
    rm.set(jc, jc, q2, u2, 1.0)
    rm.set(jc, jc, q1, u2, 1.0)
    rm.set(jc, jc, u2, q1, 1.0)
    m.link(rm)
    q1.setState(3)
    q2.setState(0)
    return m


def _metrics(model, cutoff):
    solver = SolverCTMC(model, cutoff=cutoff)
    return (np.asarray(solver.getAvgQLen()).ravel(),
            np.asarray(solver.getAvgTput()).ravel())


def test_round_trip_keeps_every_metric(tmp_path):
    model = _two_modes()
    path = str(tmp_path / 'twomodes.pnml')
    save_pnml(model, path)
    back = load_pnml(path)

    q0, x0 = _metrics(model, 4)
    q1, x1 = _metrics(back, 4)
    assert np.allclose(q0, q1, atol=1e-12), (q0, q1)
    assert np.allclose(x0, x1, atol=1e-12), (x0, x1)


def test_round_trip_regroups_modes_and_keeps_inhibitors(tmp_path):
    model = _multi_mode()
    path = str(tmp_path / 'multimode.pnml')
    save_pnml(model, path)

    # Two modes of one transition are written as two PNML transitions and must be
    # regrouped into ONE LINE transition by the reader, not left as two.
    text = open(path, encoding='utf-8').read()
    assert 'id="U1.Fast"' in text
    assert 'id="U1.Slow"' in text
    assert '<type value="inhibitor"/>' in text

    back = load_pnml(path)
    transitions = [nd for nd in back.get_nodes() if isinstance(nd, Transition)]
    assert len(transitions) == 2, 'the two modes of U1 must regroup into one transition'
    u1 = [t for t in transitions if t.get_name() == 'U1'][0]
    u2 = [t for t in transitions if t.get_name() == 'U2'][0]
    assert u1.get_number_of_modes() == 2
    assert u2._timing_strategies[0] == TimingStrategy.IMMEDIATE
    assert u2._firing_weights[0] == pytest.approx(2.5)
    q1 = [nd for nd in back.get_nodes() if isinstance(nd, Place) and nd.get_name() == 'Q1'][0]
    # get_node_index is 1-based; the mode matrices are 0-based numpy arrays.
    assert u2._inhibiting_conditions[0][back.get_node_index(q1) - 1, 0] == pytest.approx(3.0)

    q0, x0 = _metrics(model, 3)
    qb, xb = _metrics(back, 3)
    assert np.allclose(q0, qb, atol=1e-12), (q0, qb)
    assert np.allclose(x0, xb, atol=1e-12), (x0, xb)


@pytest.mark.parametrize('build,cutoff', [(_two_modes, 4), (_multi_mode, 3)])
def test_write_read_write_is_byte_identical(build, cutoff, tmp_path):
    """A field lost in the round trip changes the second file.

    Comparing metrics can only see what the solver reads; comparing the FILE
    sees every attribute the writer emits, so this is the stronger oracle of the
    two and needs no solver at all.
    """
    first = str(tmp_path / 'first.pnml')
    second = str(tmp_path / 'second.pnml')
    save_pnml(build(), first)
    save_pnml(load_pnml(first), second)
    assert open(first, encoding='utf-8').read() == open(second, encoding='utf-8').read()


def test_untimed_net_from_another_tool_reads_and_solves(tmp_path):
    # No toolspecific block and no timing: every transition is read as Exp(1)
    # with one server. Two indistinguishable tokens then cycle between two
    # places, so the number in `busy` is a 3-state birth-death chain with equal
    # rates; the stationary distribution is uniform and each place holds a mean
    # of one token, with the cycle firing at 2/3.
    xml = '''<?xml version="1.0" encoding="UTF-8"?>
<pnml xmlns="http://www.pnml.org/version-2009/grammar/pnml">
  <net id="foreign" type="http://www.pnml.org/version-2009/grammar/ptnet">
    <page id="p0">
      <place id="ready"><initialMarking><text>2</text></initialMarking></place>
      <place id="busy"><initialMarking><text>0</text></initialMarking></place>
      <transition id="start"/>
      <transition id="finish"/>
      <arc id="e1" source="ready" target="start"/>
      <arc id="e2" source="start" target="busy"/>
      <arc id="e3" source="busy" target="finish"/>
      <arc id="e4" source="finish" target="ready"/>
    </page>
  </net>
</pnml>
'''
    path = str(tmp_path / 'foreign.pnml')
    with open(path, 'w', encoding='utf-8') as fh:
        fh.write(xml)

    model = load_pnml(path)
    qlen, tput = _metrics(model, 2)
    assert qlen[0] == pytest.approx(1.0, abs=1e-9)
    assert qlen[1] == pytest.approx(1.0, abs=1e-9)
    assert tput[1] == pytest.approx(2.0 / 3.0, abs=1e-9)


def test_open_net_is_refused(tmp_path):
    model = Network('open')
    Source(model, 'Source')
    Place(model, 'P1')
    Sink(model, 'Sink')
    OpenClass(model, 'Class1', 0)
    with pytest.raises(ValueError) as excinfo:
        save_pnml(model, str(tmp_path / 'bad.pnml'))
    assert 'unbounded token source' in str(excinfo.value)


def test_multiclass_net_is_refused(tmp_path):
    model = _two_modes()
    ClosedClass(model, 'Class2', 1, model.get_nodes()[0], 0)
    with pytest.raises(ValueError) as excinfo:
        save_pnml(model, str(tmp_path / 'bad.pnml'))
    assert 'UNCOLOURED' in str(excinfo.value)
