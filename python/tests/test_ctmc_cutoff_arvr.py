"""
Native-Python test for the CTMC offered-vs-carried ArvR convention at a bound.

A finite bound on an open class is EITHER a physical capacity (setCapacity /
setClassCapacity) OR a state-space cutoff imposed only to keep enumeration
finite. These must report the arrival rate (ArvR) differently, and the
distinguishing signal is the station's drop rule (``State.is_physical_capacity``
/ ``_arrival_is_lost`` in ``api/state/after_event_station.py``, applied by the
CTMC state-space generator ``ctmc_ssg.py``):

  * PHYSICAL cap  -> a job that meets a full buffer is really lost, so the loss
    is a real event: ArvR reports the OFFERED rate (lambda), Tput < lambda.
  * state-space CUTOFF -> the refused job never existed; it is a truncation
    artifact and counts nowhere: ArvR reports the CARRIED rate (== Tput).

This is the USER DECISION of 2026-07-17 (a cutoff is not physical capacity),
documented in ``_kb/06-solver-catalog.md``. Fixture: open M/M/1,
Source(Exp 0.9) -> Queue FCFS(Exp 1.0) -> Sink, so lambda=0.9, mu=1.0.

The BINDING control is that the two configurations share an IDENTICAL stationary
distribution (M/M/1/2 == M/M/1 truncated at 2, both give Tput 0.630996) yet must
report DIFFERENT ArvR: 0.9 (cap) vs 0.630996 (cutoff). A solver that keyed ArvR
off the state space alone -- or folded the cutoff into a physical cap -- would
return the same number for both and fail here.
"""

import numpy as np
import pytest

from line_solver import (Network, Source, Queue, Sink, OpenClass, Exp,
                         SchedStrategy, SolverCTMC)

LAM, MU = 0.9, 1.0
QI, SI = 1, 0  # station rows: 0 = Source, 1 = Queue


def _build(capacity=None):
    model = Network('MM1')
    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')
    oclass = OpenClass(model, 'C1')
    source.setArrival(oclass, Exp(LAM))
    queue.setService(oclass, Exp(MU))
    if capacity is not None:
        queue.setCapacity(capacity)
    model.link(Network.serialRouting(source, queue, sink))
    return model


def _solve(capacity, cutoff):
    solver = SolverCTMC(_build(capacity), cutoff=cutoff, verbose=False)
    AN = np.asarray(solver.getAvgArvR())
    TN = np.asarray(solver.getAvgTput())
    return AN[QI, 0], TN[QI, 0]


@pytest.mark.parametrize("cutoff", [2, 3, 5, 8])
def test_cutoff_reports_carried_arvr(cutoff):
    # infinite capacity + tight cutoff: ArvR must equal the carried throughput,
    # and both must sit strictly below the offered lambda.
    arvr, tput = _solve(capacity=None, cutoff=cutoff)
    assert arvr == pytest.approx(tput, abs=1e-9)
    assert arvr < LAM - 1e-6


def test_physical_cap_reports_offered_arvr():
    # finite physical capacity K with a loose cutoff: ArvR must equal the offered
    # lambda while the carried throughput stays below it (real loss at the cap).
    for K in (2, 3, 5):
        arvr, tput = _solve(capacity=K, cutoff=max(K + 5, 20))
        assert arvr == pytest.approx(LAM, abs=1e-9)
        assert tput < LAM - 1e-6


def test_cutoff_and_cap_share_distribution_but_differ_in_arvr():
    # The decisive control: K=2 (physical) and cutoff=2 (truncation) yield the
    # SAME stationary throughput but must report DIFFERENT ArvR.
    arvr_cut, tput_cut = _solve(capacity=None, cutoff=2)
    arvr_cap, tput_cap = _solve(capacity=2, cutoff=20)
    assert tput_cut == pytest.approx(tput_cap, abs=1e-9)   # identical carried rate
    assert arvr_cut == pytest.approx(tput_cut, abs=1e-9)   # cutoff -> carried
    assert arvr_cap == pytest.approx(LAM, abs=1e-9)        # cap    -> offered
    assert arvr_cap > arvr_cut + 1e-3                      # and they differ
