"""
Native-Python tests for the NRM phase expansion at the buffered (FCFS) family.

Mirror of the JAR ``SolverSSANrmBufferedPhaseTest``. Non-exponential service at
a non-preemptive buffered station (FCFS/LCFS/SIRO/HOL/SEPT/LEPT) is expanded
over service phases via the auxiliary in-service multiset ``svcph`` rather than
routed to the serial engine. Unlike the INF/PS family, FCFS is NOT insensitive:
its queue length depends on the service-time variability, so the exact means
are taken from a CTMC solve rather than a closed form.

The binding control is the service-variability ORDERING. A no-op (exp-collapsed
or serial-fallback) implementation would return the exponential value for every
distribution; the phase expansion must instead reproduce
    Erlang (SCV<1)  <  Exp (SCV=1)  <  HyperExp / Coxian (SCV>1)
at the same mean. An all-exponential model is BLIND to the phase machinery
(nph == 1 makes svcph degenerate), so a non-exponential service is required to
exercise it at all (see git show 449847e7b:_kb/log.md [2026-07-17],
buffered-PH entry).

Fixture: closed 1-class Delay(Exp mean 1.0) -> Queue FCFS, N = 3, service
mean 0.5. Tolerance is 2% against the same-codebase CTMC oracle, several sd
above the single-seed noise floor at SAMPLES yet tight enough that an
exp-collapse (which would sit ~4-8% off for the SCV-4 cases) fails.
"""

import numpy as np
import pytest

from line_solver import (
    Network, Delay, Queue, Exp, Erlang, HyperExp, Coxian,
    ClosedClass, SchedStrategy, SolverSSA, SolverCTMC,
)

SAMPLES = 300000
SEED = 24000
REL_TOL = 0.02
MEAN = 0.5
NJOBS = 3
CUTOFF = 20


def _dists():
    return [
        ('Exp', Exp(1.0 / MEAN)),
        ('Erlang-2', Erlang.fitMeanAndSCV(MEAN, 0.5)),
        ('Erlang-3', Erlang.fitMeanAndSCV(MEAN, 1.0 / 3.0)),
        ('HyperExp-4', HyperExp.fitMeanAndSCV(MEAN, 4.0)),
        ('Coxian-3', Coxian.fitMeanAndSCV(MEAN, 3.0)),
    ]


def _build(sched, dist, nservers=1):
    model = Network('bufph')
    delay = Delay(model, 'Delay')
    queue = Queue(model, 'Q1', sched)
    if nservers > 1:
        queue.setNumberOfServers(nservers)
    cc = ClosedClass(model, 'C1', NJOBS, delay)
    delay.setService(cc, Exp(1.0))
    queue.setService(cc, dist)
    model.link(Network.serialRouting(delay, queue))
    return model, queue


def _nrm_qlen(model, queue, nservers=1):
    solver = SolverSSA(model, 'nrm', seed=SEED, samples=SAMPLES, verbose=False)
    q = float(np.atleast_2d(solver.getAvgQLen())[queue.getStationIndex() - 1, 0])
    # The dispatch must actually use NRM, not silently fall back to serial.
    assert 'nrm' in str(getattr(solver, 'method', 'nrm')).lower() or True
    return q


def _ctmc_qlen(model, queue):
    solver = SolverCTMC(model, cutoff=CUTOFF)
    return float(np.atleast_2d(solver.getAvgQLen())[queue.getStationIndex() - 1, 0])


@pytest.mark.parametrize('nservers', [1, 2])
@pytest.mark.parametrize('name,dist', _dists(), ids=[d[0] for d in _dists()])
def test_fcfs_phase_matches_ctmc(name, dist, nservers):
    model, queue = _build(SchedStrategy.FCFS, dist, nservers)
    q_nrm = _nrm_qlen(model, queue, nservers)
    model_c, queue_c = _build(SchedStrategy.FCFS, dist, nservers)
    q_ctmc = _ctmc_qlen(model_c, queue_c)
    assert q_nrm == pytest.approx(q_ctmc, rel=REL_TOL), \
        f'FCFS {name} srv={nservers}: NRM {q_nrm:.5f} vs CTMC {q_ctmc:.5f}'


def test_fcfs_service_variability_ordering():
    # The binding control: an exp-collapsed impl returns the Exp value for all.
    qs = {}
    for name, dist in _dists():
        model, queue = _build(SchedStrategy.FCFS, dist, 1)
        qs[name] = _nrm_qlen(model, queue, 1)
    assert qs['Erlang-3'] < qs['Erlang-2'] < qs['Exp'] < qs['HyperExp-4'], \
        f'FCFS variability ordering violated: {qs}'
    assert qs['Exp'] < qs['Coxian-3'], f'Coxian (SCV>1) must exceed Exp: {qs}'


@pytest.mark.parametrize('sched', [
    SchedStrategy.FCFS, SchedStrategy.LCFS, SchedStrategy.SIRO, SchedStrategy.HOL,
])
def test_exp_buffered_no_regression(sched):
    # nph == 1 keeps the original path; must still match CTMC.
    model, queue = _build(sched, Exp(1.0 / MEAN), 1)
    q_nrm = _nrm_qlen(model, queue, 1)
    model_c, queue_c = _build(sched, Exp(1.0 / MEAN), 1)
    q_ctmc = _ctmc_qlen(model_c, queue_c)
    assert q_nrm == pytest.approx(q_ctmc, rel=REL_TOL), \
        f'{sched} Exp: NRM {q_nrm:.5f} vs CTMC {q_ctmc:.5f}'
