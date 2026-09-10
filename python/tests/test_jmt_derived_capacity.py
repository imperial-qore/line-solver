"""A DERIVED station capacity is not a buffer, so SolverJMT must not refuse it.

`refresh_capacity` gives EVERY station of a closed model a finite `sn.cap` even
when nobody called `setCapacity`: it derives it as
``min(sum_c chaincap, sum_r classcap)``, and the classcap row is the CHAIN
population repeated once per class served there. A Queue serving two classes of
one 4-job chain therefore carries 8 -- a bound those 4 jobs can never reach.

The writer's reachability test was ``cap == sum(njobs)``, an EQUALITY that only
the SINGLE-class case satisfies, so every multi-class station fell through to
``_jmt_station_cap_assert`` and was refused as a "binding" buffer the model never
declared. It is ``>=`` now. The genuine refusal -- a capacity the population CAN
reach -- is asserted here too, because widening the escape hatch must not open
it for the model BUG-81 is about.
"""

import numpy as np
import pytest

from line_solver import (Network, Delay, Queue, ClosedClass, Exp, SchedStrategy)
from line_solver.api.solvers.jmt.handler import SolverJMTOptions, _write_jsim_file


def _two_class_one_chain(n_each=2, cap=None):
    """Delay+Queue, both classes served at both stations, ONE chain.

    Two INDEPENDENT closed classes would be two chains, and then both chaincap
    columns are set and the two sums agree at the population, so the inflation
    does not arise. The 1/2 class switch at each station is what merges them.
    """
    m = Network('twoclass')
    d = Delay(m, 'Delay')
    q = Queue(m, 'Queue1', SchedStrategy.PS)
    c1 = ClosedClass(m, 'Class1', n_each, d)
    c2 = ClosedClass(m, 'Class2', n_each, d)
    d.setService(c1, Exp(1.0))
    d.setService(c2, Exp(1.0))
    q.setService(c1, Exp(2.0))
    q.setService(c2, Exp(2.0))
    if cap is not None:
        q.setCapacity(cap)
    P = m.initRoutingMatrix()
    # Every (class, node) row sums to 1.
    P.set(c1, c1, d, q, 0.5)
    P.set(c1, c2, d, q, 0.5)
    P.set(c2, c1, d, q, 0.5)
    P.set(c2, c2, d, q, 0.5)
    P.set(c1, c1, q, d, 0.5)
    P.set(c1, c2, q, d, 0.5)
    P.set(c2, c1, q, d, 0.5)
    P.set(c2, c2, q, d, 0.5)
    m.link(P)
    return m


def _sizes(sn, model, tmp_path):
    path = str(tmp_path / 'm.jsimg')
    _write_jsim_file(sn, path, SolverJMTOptions(), model)
    import xml.etree.ElementTree as ET
    tree = ET.parse(path)
    return [e.find('value').text for e in tree.iter('parameter')
            if e.get('name') == 'size']


def test_derived_multiclass_capacity_is_not_a_buffer(tmp_path):
    m = _two_class_one_chain()
    sn = m.getStruct()
    assert int(sn.nchains) == 1, 'the two classes must share one chain'
    assert float(np.sum(sn.njobs)) == 4.0
    # (2 classes served) x (chain population 4), every bit of it derived
    assert float(np.max(np.asarray(sn.cap).flatten())) == 8.0

    sizes = _sizes(sn, m, tmp_path)
    assert sizes, 'the queue section must still carry a size parameter'
    assert all(s == '-1' for s in sizes), \
        'a bound the population cannot reach is JMT unbounded, got %s' % sizes


def test_declared_capacity_above_the_population_is_not_a_buffer(tmp_path):
    """4 jobs cannot fill 100 either: only `cap < total` is a buffer."""
    m = _two_class_one_chain(cap=100)
    sizes = _sizes(m.getStruct(), m, tmp_path)
    assert all(s == '-1' for s in sizes), sizes


def test_binding_closed_capacity_is_still_refused(tmp_path):
    """BUG-81's model must keep its refusal: cap 2 with 6 jobs DOES bind."""
    m = Network('tandem')
    q = [Queue(m, 'Q%d' % (i + 1), SchedStrategy.FCFS) for i in range(3)]
    c = ClosedClass(m, 'C', 6, q[0])
    for qi in q:
        qi.setService(c, Exp(1.0))
    q[1].setCapacity(2)
    m.link(Network.serialRouting(q[0], q[1], q[2]))

    with pytest.raises(Exception) as exc:
        _sizes(m.getStruct(), m, tmp_path)
    msg = str(exc.value)
    assert 'closed class' in msg and 'Q2' in msg, msg
