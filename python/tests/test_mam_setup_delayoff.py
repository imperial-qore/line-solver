"""
Regression tests: MAM setup/delay-off dispatch.

Queue setup/delay-off reached no MAM code path at all: sn.hassetup was never
populated and no writer emitted the nodeparam entry the reader expected, so an
M/M/1 with setup silently returned the setup-free rho/(1-rho).

These tests cover the OPEN station only, where the delay-off races the aggregate
interarrival Exp(lambda) -- the idle a Poisson stream sees -- so
qbd_setupdelayoff is exact. References are the MATLAB values, which agree with a
4e6-sample JMT simulation (1.2290 and 3.4162).

The CLOSED station is pinned too, since 2026-09: qbd_setupdelayoff_closed solves
the finite level-dependent chain exactly, so the reference is the SIMULATORS
rather than a MATLAB row. LDES reads 0.9009 / 1.1157 / 3.0842 and JMT
0.9002 / 1.1137 / 3.0703 at setup means none / 0.5 / 5.0, and the analysis lands
between them at every point. What stood before was the per-instance cold-start
race p_cold*E[setup] + S, which carried no queueing term and reported the same
number across a tenfold change in the setup mean (BUG-78).
"""
import numpy as np
import pytest

from line_solver import (ClosedClass, Delay, Exp, Network, OpenClass, Queue,
                         SchedStrategy, Sink, SolverMAM, Source)


def _open(setup=None, delayoff=None):
    model = Network('open_setup')
    source = Source(model, 'Source')
    queue = Queue(model, 'Q', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')
    oclass = OpenClass(model, 'C', 0)
    source.setArrival(oclass, Exp(0.5))
    queue.setService(oclass, Exp(1.0))
    if setup is not None:
        queue.setDelayOff(oclass, setup, delayoff)
    model.link(Network.serialRouting(source, queue, sink))
    return model


def _closed(setup=None, delayoff=None):
    model = Network('closed_setup')
    think = Delay(model, 'Think')
    queue = Queue(model, 'Q', SchedStrategy.FCFS)
    cclass = ClosedClass(model, 'C', 3, think)
    think.setService(cclass, Exp(1.0))
    queue.setService(cclass, Exp(2.0))
    if setup is not None:
        queue.setDelayOff(cclass, setup, delayoff)
    model.link(Network.serialRouting(think, queue))
    return model


def test_sn_hassetup_marks_setup_stations():
    sn = _closed(Exp(2.0), Exp(4.0)).getStruct()
    # Station-indexed: the Delay is not a function station, the Queue is.
    assert list(np.asarray(sn.hassetup).flatten()) == [False, True]


@pytest.mark.parametrize('setup_rate,expected', [
    (2.0, 1.227273),  # setup mean 0.5
    (0.2, 3.413793),  # setup mean 5.0
])
def test_open_setup_matches_matlab(setup_rate, expected):
    qlen = float(SolverMAM(_open(Exp(setup_rate), Exp(4.0))).getAvgTable().QLen[1])
    assert np.abs(qlen - expected) / expected < 1e-5


def test_open_setup_is_not_ignored():
    # Without setup the station is a plain M/M/1 holding rho/(1-rho) = 1 job.
    assert np.abs(float(SolverMAM(_open()).getAvgTable().QLen[1]) - 1.0) < 1e-6
    assert float(SolverMAM(_open(Exp(2.0), Exp(4.0))).getAvgTable().QLen[1]) > 1.0


def test_open_setup_monotone():
    qlens = [float(SolverMAM(_open(Exp(1.0 / m), Exp(4.0))).getAvgTable().QLen[1])
             for m in (0.1, 0.5, 2.0, 5.0)]
    assert all(b > a for a, b in zip(qlens, qlens[1:]))


def test_closed_without_setup_is_unaffected():
    # Single-class closed Delay+Queue now routes to the exact ldqbd method
    # (default MAM), which returns the exact 0.9 (matching MVA/NC/CTMC), not the
    # old dec.source approximation 0.7105.
    respt = float(SolverMAM(_closed()).getAvgTable().RespT[1])
    assert np.abs(respt - 0.9) < 1e-5


# ---------------------------------------------------------------------------
# Closed setup/delay-off, BUG-78. The chain is Delay(Z=1) + Queue(FCFS, D=0.5),
# one closed class with N=3 and a delay-off of Exp(4). Two independent
# simulators bracket every row below (LDES / JMT):
#     setup mean 0.5  ->  1.1157 / 1.1137
#     setup mean 5.0  ->  3.0842 / 3.0703
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('setup_mean,lo,hi', [
    (0.5, 1.1137, 1.1157),
    (5.0, 3.0703, 3.0842),
])
def test_closed_setup_lands_between_the_simulators(setup_mean, lo, hi):
    respt = float(SolverMAM(_closed(Exp(1.0 / setup_mean), Exp(4.0))).getAvgTable().RespT[1])
    # 0.5% either side of the simulator bracket: the two simulators are 0.2%
    # apart on these models and the analysis has to sit inside that, not merely
    # near it.
    assert lo * 0.995 <= respt <= hi * 1.005, respt


def test_closed_setup_is_not_ignored():
    # The defect this pins: every analytical row was byte-identical across a
    # tenfold change in the setup mean while the simulators moved 0.90 -> 3.08.
    respts = [float(SolverMAM(_closed(*args)).getAvgTable().RespT[1])
              for args in [(None, None), (Exp(2.0), Exp(4.0)), (Exp(0.2), Exp(4.0))]]
    assert respts[0] < respts[1] < respts[2]
    assert respts[2] / respts[0] > 3.0


def test_closed_setup_throughput_falls_with_the_setup():
    # A slower setup is a slower server, so the closed chain's throughput has to
    # drop with it -- the p_cold formula left it at the setup-free 1.5789.
    tputs = [float(SolverMAM(_closed(*args)).getAvgTable().Tput[1])
             for args in [(None, None), (Exp(2.0), Exp(4.0)), (Exp(0.2), Exp(4.0))]]
    assert tputs[0] > tputs[1] > tputs[2]
    assert abs(tputs[0] - 1.578947) < 1e-5


def test_closed_setup_conserves_the_population():
    # Little's law across the two stations: the jobs are either thinking or at
    # the queue, so the two queue lengths have to sum to N.
    table = SolverMAM(_closed(Exp(0.2), Exp(4.0))).getAvgTable()
    assert abs(float(table.QLen[0]) + float(table.QLen[1]) - 3.0) < 1e-6


def test_closed_setup_degenerates_to_the_plain_queue():
    # An instantaneous setup is a server that is never cold, so the vacation
    # chain has to collapse onto the exact closed queue, R = 0.9.
    respt = float(SolverMAM(_closed(Exp(1e8), Exp(4.0))).getAvgTable().RespT[1])
    assert abs(respt - 0.9) < 1e-5
