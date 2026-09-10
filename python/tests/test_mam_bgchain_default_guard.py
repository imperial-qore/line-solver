"""The MAM default sizes the background chain before choosing bgchain.

`bgchain` asked for BY NAME is refused above ``config['bgstates_max']``, which is
right: the model cannot afford the chain the user requested. ``default`` must not
inherit that refusal -- it has to land on a method that answers. Before the guard
(2026-08-14) MATLAB and the JAR both THREW on ``mqn_singleserver_fcfs`` (one
closed chain, N = 100 over 4 stations, ``nchoosek(103,3) = 176851`` states) with
"The background chain of this model has 176851 states, above the limit of 20000".

MATLAB is the reference here and prints ``default/dec.source`` on that model with
the guard in place, and ``default/bgchain`` on the small mixed twin.

Twin of the C++ "the default leaves an oversized background chain to dec.source"
in cpp/tests/test_mam_bgchain.cpp.
"""

import numpy as np

from line_solver import (ClosedClass, Delay, Exp, MAM, Network, OpenClass, Queue,
                         SchedStrategy, Sink, Source)
from line_solver.solvers.solver_mam.algorithms.bgchain import bgchain_states


def _mixed_cycle(npop):
    """Source -> Q -> Sink for the open class, Delay -> Q -> Delay for the closed one."""
    m = Network('bgGuard')
    source = Source(m, 'Source')
    delay = Delay(m, 'Delay')
    q = Queue(m, 'Queue1', SchedStrategy.PS)
    sink = Sink(m, 'Sink')
    openc = OpenClass(m, 'Open')
    closed = ClosedClass(m, 'Closed', npop, delay, 0)
    source.setArrival(openc, Exp(0.4))
    q.setService(openc, Exp(2.0))
    delay.setService(closed, Exp(1.0))
    q.setService(closed, Exp(2.0))
    P = m.initRoutingMatrix()
    P.set(openc, openc, source, q, 1.0)
    P.set(openc, openc, q, sink, 1.0)
    P.set(closed, closed, delay, q, 1.0)
    P.set(closed, closed, q, delay, 1.0)
    m.link(P)
    return m


def test_bgchain_states_counts_the_compositions():
    # 4 closed jobs over the 2 stations the closed class visits: nchoosek(5,1) = 5.
    assert bgchain_states(_mixed_cycle(4).get_struct(), {}) == 5.0
    # 100 jobs over 2 stations: 101, still small; the shape grows as nchoosek(N+M-1,M-1).
    assert bgchain_states(_mixed_cycle(100).get_struct(), {}) == 101.0


def test_default_takes_bgchain_when_the_chain_fits():
    solver = MAM(_mixed_cycle(4), seed=23000)
    solver.avg_table()
    assert solver.result.method == 'bgchain'


def test_default_falls_to_dec_source_when_the_chain_does_not_fit():
    solver = MAM(_mixed_cycle(4), seed=23000)
    solver.options.config['bgstates_max'] = 2   # below the 5 the chain needs
    solver.avg_table()
    assert solver.result.method == 'dec.source'


def test_bgchain_by_name_still_refuses_an_oversized_chain():
    solver = MAM(_mixed_cycle(4), method='bgchain', seed=23000)
    solver.options.config['bgstates_max'] = 2
    try:
        solver.avg_table()
    except Exception as e:
        assert 'above the limit' in str(e)
    else:
        raise AssertionError('expected a refusal by name')


def test_a_large_closed_population_does_not_throw_by_default():
    # The shape that broke: enough closed jobs over enough stations to put the
    # chain past bgstates_max. The default must still return a table.
    m = Network('bgBig')
    source = Source(m, 'Source')
    sink = Sink(m, 'Sink')
    qs = [Queue(m, 'Queue%d' % (i + 1), SchedStrategy.PS) for i in range(4)]
    closed = ClosedClass(m, 'ClosedClass', 100, qs[0], 0)
    openc = OpenClass(m, 'OpenClass', 0)
    for i, q in enumerate(qs):
        q.setService(closed, Exp(i + 1))
        q.setService(openc, Exp(np.sqrt(i + 1)))
    source.setArrival(openc, Exp(1.0 / 3.0))
    P = m.initRoutingMatrix()
    for i in range(4):
        P.set(closed, closed, qs[i], qs[(i + 1) % 4], 1.0)
    P.set(openc, openc, source, qs[0], 1.0)
    P.set(openc, openc, qs[0], qs[1], 1.0)
    P.set(openc, openc, qs[1], qs[2], 1.0)
    P.set(openc, openc, qs[2], sink, 1.0)
    m.link(P)

    assert bgchain_states(m.get_struct(), {}) == 176851.0
    solver = MAM(m, seed=23000)
    table = solver.avg_table()
    assert table is not None
    assert solver.result.method == 'dec.source'
