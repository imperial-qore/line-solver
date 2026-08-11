"""
Native-Python tests for the NRM's WAITQ-rule finite capacity regions, with the
region on an INTERIOR station. Mirror of the JAR ``SolverSSANrmFcrWaitqTest``.

The WAITQ rule PARKS a refused job in a per-region FIFO and admits it
head-of-line as capacity frees, in contrast to DROP, which destroys it. The two
are observationally identical when the region is fed by a Source (a Source's
population is fictitious), so a genuine WAITQ test must place the region on a
station fed by a real upstream queue. The interior fixture is a closed cycle
Think(INF) -> Q1/FCFS -> Q2/FCFS -> Think with the region over {Q2}; the parked
jobs are conserved (they re-enter Q2 later), so the sum of station queue lengths
is strictly below the closed population N, the difference being the mean FIFO
occupancy. A DROP rule there would destroy jobs and drain the closed network.

Ground truth is the exact CTMC. The exact per-station means below are hard-coded
rather than obtained from a CTMC solve so the assertion cannot drift with the
reference solver, and because Python SolverCTMC currently gates the Region
feature out of its analyzer (its handler solves FCR correctly, so this is a
featset lag, not a numeric gap). The values agree to 5+ significant figures
across the MATLAB, JAR, and Python-native CTMC handlers:

    interior {Q2}, cap 2, N=4:  QLen = [0.966672, 1.221896, 1.379727]
                                Tput = [0.966672, 0.966672, 0.966672]
                                sum(QLen) = 3.5683, parked mean = 0.4317
    source-fed {Q1}, cap 3:     Q1 QLen = 1.175991, Q1 Tput = 0.6

Each case also asserts the result method is 'nrm': the analyzer downgrades an
ineligible model to the serial engine SILENTLY, and the serial engine is correct
here, so a fallback would pass this test while never exercising the WAITQ FIFO it
exists to cover.
"""

import os
import sys

# Ensure the WORKTREE line_solver is imported, not a global checkout whose older
# NRM rejects FCR and reads as a silent serial fallback. The package root is two
# levels up from this test file (python/), so prepend it before any global path.
_WORKTREE_PY = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
if sys.path and sys.path[0] != _WORKTREE_PY:
    sys.path.insert(0, _WORKTREE_PY)

import numpy as np
import pytest

from line_solver import (
    Network, Delay, Queue, Exp, ClosedClass, OpenClass, Source, Sink,
    SchedStrategy, DropStrategy, SolverSSA,
)

# Fail loudly if a stray global package shadowed the worktree one.
import line_solver as _ls
assert os.path.abspath(_ls.__file__).startswith(_WORKTREE_PY), (
    f'wrong line_solver imported: {_ls.__file__} (expected under {_WORKTREE_PY})')

SAMPLES = 300000
SEEDS = [23000, 24000, 25000, 26000]
ATOL = 0.03

# Exact CTMC means (MATLAB == JAR == Python CTMC handler, 5+ sig figs).
INTERIOR_QLEN = np.array([0.966672, 1.221896, 1.379727])
INTERIOR_TPUT = np.array([0.966672, 0.966672, 0.966672])
SRCFED_Q1_QLEN = 1.175991
SRCFED_Q1_TPUT = 0.6


def _interior_model():
    """Closed interior-region WAITQ model: Think -> Q1 -> Q2{region cap} -> Think."""
    model = Network('wq_interior')
    think = Delay(model, 'Think')
    q1 = Queue(model, 'Q1', SchedStrategy.FCFS)
    q2 = Queue(model, 'Q2', SchedStrategy.FCFS)
    c1 = ClosedClass(model, 'C1', 4, think)
    think.setService(c1, Exp(1.0))
    q1.setService(c1, Exp(1.5))
    q2.setService(c1, Exp(1.2))
    model.link(Network.serialRouting(think, q1, q2))
    region = model.addRegion([q2])
    region.setGlobalMaxJobs(2)
    region.setDropRule(c1, DropStrategy.WaitingQueue)
    return model


def _source_fed_model():
    """Open source-fed WAITQ model (non-regression): Source -> Q1{region cap} -> Sink."""
    model = Network('wq_srcfed')
    source = Source(model, 'Source')
    q1 = Queue(model, 'Q1', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')
    c1 = OpenClass(model, 'C1')
    source.setArrival(c1, Exp(0.6))
    q1.setService(c1, Exp(1.0))
    model.link(Network.serialRouting(source, q1, sink))
    region = model.addRegion([q1])
    region.setGlobalMaxJobs(3)
    region.setDropRule(c1, DropStrategy.WaitingQueue)
    return model


def _nrm_mean(build):
    """Mean over SEEDS of (QLen, Tput), asserting the NRM engine actually ran."""
    qs, ts = [], []
    for seed in SEEDS:
        solver = SolverSSA(build(), 'nrm', seed=seed, samples=SAMPLES, verbose=False)
        q = np.atleast_2d(solver.getAvgQLen())
        t = np.atleast_2d(solver.getAvgTput())
        ran = getattr(solver._result, 'method', None)
        assert ran == 'nrm', f'NRM did not run: method was {ran!r}'
        qs.append(q)
        ts.append(t)
    return np.mean(qs, axis=0), np.mean(ts, axis=0)


def test_fcr_waitq_interior_region_matches_ctmc():
    """Interior region {Q2} in a closed cycle: every station queue length and
    throughput must match the exact CTMC, and the parked mass (N minus the sum
    of station queue lengths) must be strictly positive -- the observable WAITQ
    signature that a DROP rule cannot produce."""
    nrm_q, nrm_t = _nrm_mean(_interior_model)
    for i in range(INTERIOR_QLEN.size):
        assert abs(nrm_q[i, 0] - INTERIOR_QLEN[i]) < ATOL, (
            f'interior station {i} QLen: NRM {nrm_q[i, 0]:.5f}, CTMC {INTERIOR_QLEN[i]:.5f}')
        assert abs(nrm_t[i, 0] - INTERIOR_TPUT[i]) < ATOL, (
            f'interior station {i} Tput: NRM {nrm_t[i, 0]:.5f}, CTMC {INTERIOR_TPUT[i]:.5f}')
    # The parked jobs live in the region FIFO, not at any station, so the
    # simulated station queue lengths must sum to strictly below the population.
    sum_nrm = float(np.sum(nrm_q[:INTERIOR_QLEN.size, 0]))
    assert sum_nrm < 3.9, (
        f'WAITQ FIFO must hold jobs off-station: sum of station QLen was '
        f'{sum_nrm:.4f} (population is 4)')


def test_fcr_waitq_source_fed_region_matches_ctmc():
    """Source-fed region {Q1}: the pre-existing configuration, kept as a
    non-regression. A source-fed WAITQ region blocks the source, so the model is
    bounded and its station queue length and throughput match the exact CTMC."""
    nrm_q, nrm_t = _nrm_mean(_source_fed_model)
    # Station row 1 is Q1 (Source is station 0).
    assert abs(nrm_q[1, 0] - SRCFED_Q1_QLEN) < ATOL, (
        f'source-fed Q1 QLen: NRM {nrm_q[1, 0]:.5f}, CTMC {SRCFED_Q1_QLEN:.5f}')
    assert abs(nrm_t[1, 0] - SRCFED_Q1_TPUT) < ATOL, (
        f'source-fed Q1 Tput: NRM {nrm_t[1, 0]:.5f}, CTMC {SRCFED_Q1_TPUT:.5f}')
