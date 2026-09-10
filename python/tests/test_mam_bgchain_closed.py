"""The bgchain method on a PURELY CLOSED model.

With no open class the fixed point has nothing to iterate: ``cshare`` never
leaves its initial ``min(e,c)``, which is exactly the capacity the closed jobs
hold when no open work competes for it, so the background chain alone answers
and it answers with the EXACT closed CTMC at chain granularity. The tandem below
is product-form, so exact MVA is an independent oracle and the agreement is
machine precision rather than an approximation tolerance.

The remaining tests pin the DEFAULT chooser: bgchain takes a closed model whose
chain fits in ``bgstates_max`` and whose service laws it represents exactly, and
stands aside otherwise. The background chain is built from the MEAN service time
alone, exact at PS/INF by insensitivity and under any discipline when the law IS
exponential, so a non-exponential law at an FCFS station belongs to mna, which
carries the phase-type.
"""

import numpy as np
import pytest

from line_solver import (ClosedClass, Delay, Erlang, Exp, Network, OpenClass, Queue,
                         SchedStrategy, Sink, Source, SolverMAM, SolverMVA)


def _tandem(mq, n, sched, erlang=False):
    """Delay -> Q1 -> ... -> Qmq -> Delay, one closed class of n jobs."""
    m = Network('closedTandem')
    d = Delay(m, 'Think')
    q = [Queue(m, 'Q%d' % (i + 1), sched) for i in range(mq)]
    c = ClosedClass(m, 'C', n, d, 0)
    d.setService(c, Exp(1.0))
    for i in range(mq):
        rate = 1.0 + 0.3 * (i + 1)
        q[i].setService(c, Erlang.fitMeanAndOrder(1.0 / rate, 3) if erlang else Exp(rate))
    m.link(Network.serialRouting([d] + q))
    return m


def _qlen(solver):
    return np.asarray(solver.getAvgQLen(), dtype=float).ravel()


@pytest.mark.parametrize('mq', [2, 4, 8])
def test_closed_tandem_is_the_exact_closed_ctmc(mq):
    got = _qlen(SolverMAM(_tandem(mq, 5, SchedStrategy.PS), 'bgchain'))
    ref = _qlen(SolverMVA(_tandem(mq, 5, SchedStrategy.PS)))
    assert np.allclose(got, ref, rtol=0, atol=1e-9)
    assert got.sum() == pytest.approx(5.0, abs=1e-9)


def test_closed_fcfs_tandem_is_the_same_chain():
    got = _qlen(SolverMAM(_tandem(3, 4, SchedStrategy.FCFS), 'bgchain'))
    ref = _qlen(SolverMVA(_tandem(3, 4, SchedStrategy.FCFS)))
    assert np.allclose(got, ref, rtol=0, atol=1e-9)
    assert got.sum() == pytest.approx(4.0, abs=1e-9)


def test_default_takes_bgchain_on_a_closed_model():
    solver = SolverMAM(_tandem(4, 5, SchedStrategy.PS))
    solver.getAvgQLen()
    assert solver.result.method == 'bgchain'


def test_default_leaves_a_non_exponential_fcfs_model_to_mna():
    solver = SolverMAM(_tandem(2, 4, SchedStrategy.FCFS, erlang=True))
    solver.getAvgQLen()
    assert solver.result.method != 'bgchain'


def test_a_chain_above_bgstates_max_falls_back_rather_than_throwing():
    solver = SolverMAM(_tandem(6, 40, SchedStrategy.PS))
    QN = _qlen(solver)
    assert solver.result.method != 'bgchain'
    assert QN.size > 0


def test_a_purely_open_model_is_refused():
    m = Network('open')
    src = Source(m, 'Source')
    q = Queue(m, 'Q1', SchedStrategy.FCFS)
    snk = Sink(m, 'Sink')
    o = OpenClass(m, 'Open', 0)
    src.setArrival(o, Exp(0.5))
    q.setService(o, Exp(1.0))
    m.link(Network.serialRouting(src, q, snk))
    with pytest.raises(Exception):
        SolverMAM(m, 'bgchain').getAvgQLen()


def _two_chain_cycle(rate1, rate2, sched):
    """Delay -> Q -> Delay for two closed chains of 2 jobs each."""
    m = Network('twoChain')
    d = Delay(m, 'D')
    q = Queue(m, 'Q1', sched)
    c1 = ClosedClass(m, 'C1', 2, d, 0)
    c2 = ClosedClass(m, 'C2', 2, d, 0)
    d.setService(c1, Exp(1.0))
    d.setService(c2, Exp(1.0))
    q.setService(c1, Exp(rate1))
    q.setService(c2, Exp(rate2))
    P = m.initRoutingMatrix()
    P[c1] = Network.serial_routing([d, q, d])
    P[c2] = Network.serial_routing([d, q, d])
    m.link(P)
    return m


def test_default_leaves_class_dependent_fcfs_rates_alone():
    """The chain splits a station's capacity in proportion to the job COUNTS,
    i.e. random order. Under FCFS that is exact only at one shared rate: with
    Exp(3) against Exp(0.8) it reads 25.2% off SolverCTMC, so the default must
    stand aside."""
    solver = SolverMAM(_two_chain_cycle(3.0, 0.8, SchedStrategy.FCFS))
    solver.getAvgQLen()
    assert solver.result.method != 'bgchain'


def test_default_takes_class_independent_fcfs_rates():
    solver = SolverMAM(_two_chain_cycle(1.5, 1.5, SchedStrategy.FCFS))
    solver.getAvgQLen()
    assert solver.result.method == 'bgchain'


def test_default_takes_class_dependent_rates_under_ps():
    """PS splits by DEMAND, not by count, so class-dependent rates cost nothing
    there: measured 0.0e+00 against SolverCTMC."""
    solver = SolverMAM(_two_chain_cycle(3.0, 0.8, SchedStrategy.PS))
    solver.getAvgQLen()
    assert solver.result.method == 'bgchain'
