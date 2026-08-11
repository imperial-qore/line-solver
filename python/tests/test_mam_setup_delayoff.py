"""
Regression tests: MAM setup/delay-off dispatch.

Queue setup/delay-off reached no MAM code path at all: sn.isfunction was never
populated and no writer emitted the nodeparam entry the reader expected, so an
M/M/1 with setup silently returned the setup-free rho/(1-rho).

These tests cover the OPEN station only, where the delay-off races the aggregate
interarrival Exp(lambda) -- the idle a Poisson stream sees -- so
qbd_setupdelayoff is exact. References are the MATLAB values, which agree with a
4e6-sample JMT simulation (1.2290 and 3.4162).

Closed setup is deliberately not pinned here: MAM's closed analysis is an
approximation (it reads 0.7105 on this model against an exact 0.9 from
MVA/NC/CTMC), and the codebases disagree on whether the setup branch is reached
at all for a plain closed queue. See _kb/log.md.
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


def test_sn_isfunction_marks_setup_stations():
    sn = _closed(Exp(2.0), Exp(4.0)).getStruct()
    # Station-indexed: the Delay is not a function station, the Queue is.
    assert list(np.asarray(sn.isfunction).flatten()) == [False, True]


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
