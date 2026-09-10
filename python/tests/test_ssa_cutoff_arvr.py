"""
Native-Python test for the SSA offered-vs-carried ArvR convention at a bound.

The twin of ``test_ctmc_cutoff_arvr.py`` for the SERIAL SSA engine, which owns
its own copy of the rule (``_recompute_capacity`` returns the ``is_physical_cap``
predicate in ``api/solvers/ssa/serial.py``, and the arrival/departure tally skips
a refusal that predicate calls an artifact). The CTMC test cannot cover it: the
gate there is reached through ``State.afterEventStation``, a different code path.

A finite bound on an open class is EITHER a physical capacity (setCapacity /
setClassCapacity) OR a state-space cutoff imposed only to keep the simulated
buffer finite, and the two must report the arrival rate differently:

  * PHYSICAL cap  -> a job that meets a full buffer is really lost, so the loss
    is a real event: ArvR reports the OFFERED rate (lambda), Tput < lambda.
  * state-space CUTOFF -> the refused job never existed; it is a truncation
    artifact and counts nowhere: ArvR reports the CARRIED rate (== Tput).

This is the USER DECISION of 2026-07-17 (a cutoff is not physical capacity),
documented in ``_kb/06-solver-catalog.md``. Fixture: open M/M/1,
Source(Exp 0.9) -> Queue FCFS(Exp 1.0) -> Sink, so lambda=0.9, mu=1.0.

The BINDING control is that K=2 and cutoff=2 drive an IDENTICAL sample path
(same seed, same refusals, same stationary law, Tput 0.6305) yet must report
DIFFERENT ArvR: 0.9 against 0.6307. An engine that keyed ArvR off the buffer
bound alone -- or folded the cutoff into a physical cap -- returns one number
for both and fails here.

Tolerances are simulation tolerances: ArvR and Tput at a cutoff come from two
accumulators that differ by the jobs in the station when the run ends, so they
converge rather than coincide. The DISCRIMINATOR the test really rests on is
0.27 wide, two orders above that noise.

Note that ``method='serial'`` is explicit: the default SSA method routes to the
NRM engine, which grows its buffers dynamically and does not truncate at all.
"""

import numpy as np
import pytest

from line_solver import (Network, Source, Queue, Sink, OpenClass, Exp,
                         SchedStrategy, SolverSSA)

LAM, MU = 0.9, 1.0
QI = 1              # station rows: 0 = Source, 1 = Queue
SAMPLES, SEED = 50000, 23000

# Exact CTMC values of the same fixture, from test_ctmc_cutoff_arvr.py.
CARRIED = {2: 0.630996, 5: 0.786580}


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
    solver = SolverSSA(_build(capacity), method='serial', cutoff=cutoff,
                       samples=SAMPLES, seed=SEED, verbose=False)
    AN = np.asarray(solver.getAvgArvR())
    TN = np.asarray(solver.getAvgTput())
    return AN[QI, 0], TN[QI, 0]


@pytest.mark.parametrize("cutoff", [2, 5])
def test_cutoff_reports_carried_arvr(cutoff):
    # infinite capacity + tight cutoff: ArvR must track the carried throughput,
    # and both must sit well below the offered lambda.
    arvr, tput = _solve(capacity=None, cutoff=cutoff)
    assert arvr == pytest.approx(tput, rel=2e-2)
    assert arvr == pytest.approx(CARRIED[cutoff], rel=2e-2)
    assert arvr < LAM - 0.05


def test_no_bound_reports_offered_arvr():
    # the unbounded control: with nothing refusing anything, carried == offered,
    # so a solver that always reported the carried rate would still pass the
    # cutoff cases above. This pins the other end.
    arvr, tput = _solve(capacity=None, cutoff=np.inf)
    assert arvr == pytest.approx(LAM, abs=5e-3)
    assert tput == pytest.approx(LAM, abs=2e-2)


def test_physical_cap_reports_offered_arvr():
    # finite physical capacity K with a loose cutoff: ArvR must be the offered
    # lambda while the carried throughput stays well below it (real loss).
    for K in (2, 5):
        arvr, tput = _solve(capacity=K, cutoff=max(K + 5, 20))
        assert arvr == pytest.approx(LAM, abs=5e-3)
        assert tput == pytest.approx(CARRIED[K], rel=2e-2)
        assert tput < LAM - 0.05


def test_cutoff_and_cap_share_path_but_differ_in_arvr():
    # The decisive control: K=2 (physical) and cutoff=2 (truncation) refuse the
    # same jobs and carry the same rate, but must report DIFFERENT ArvR.
    arvr_cut, tput_cut = _solve(capacity=None, cutoff=2)
    arvr_cap, tput_cap = _solve(capacity=2, cutoff=20)
    assert tput_cut == pytest.approx(tput_cap, rel=2e-2)   # same carried rate
    assert arvr_cut == pytest.approx(tput_cut, rel=2e-2)   # cutoff -> carried
    assert arvr_cap == pytest.approx(LAM, abs=5e-3)        # cap    -> offered
    assert arvr_cap - arvr_cut > 0.2                       # and they differ
