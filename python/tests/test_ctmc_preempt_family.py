"""CTMC parity for the eight preempt disciplines.

The family factors as {FCFS,LCFS} x {PR,PI} x {plain,PRIO}. Every reference
value below was produced by MATLAB SolverCTMC on the same model and agrees to
the printed precision; the class-1 queue length under the FCFS PRIO variants is
additionally 13/11 = 1.181818..., which is what the highest-priority class must
see when it is never preempted (it is then an isolated 2-job closed model with
Exp(1) think time and Erlang-2 service, solvable by hand).

The disciplines differ in exactly two ways, and both are pinned here:
  - PR resumes a promoted job in the phase it was preempted at, PI restarts it
    from the entry distribution, so PI != PR on a multi-phase service;
  - within one priority group the base discipline decides whether an arrival
    preempts: LCFS-PR keeps the newest job in service (it does), FCFS-PR does
    not. That is why LCFS*PRIO != FCFS*PRIO even though the priorities are the
    same.
"""

import numpy as np
import pytest

from line_solver import (Network, Delay, Queue, ClosedClass, Exp, Erlang,
                         SchedStrategy, SolverCTMC)

# station-major: [Think-C1, Think-C2, Q1-C1, Q1-C2]
EXPECTED = {
    'FCFSPR':     [0.625000, 0.312500, 1.375000, 0.687500],
    'FCFSPI':     [0.507398, 0.310141, 1.492602, 0.689859],
    'FCFSPRPRIO': [0.818182, 0.117439, 1.181818, 0.882561],
    'FCFSPIPRIO': [0.818182, 0.101112, 1.181818, 0.898888],
    'LCFSPR':     [0.625000, 0.312500, 1.375000, 0.687500],
    'LCFSPI':     [0.495885, 0.320988, 1.504115, 0.679012],
    'LCFSPRPRIO': [0.800000, 0.129146, 1.200000, 0.870854],
    'LCFSPIPRIO': [0.750000, 0.092744, 1.250000, 0.907256],
}


def _build(sched):
    model = Network('preempt')
    d = Delay(model, 'Think')
    q = Queue(model, 'Q1', sched)
    q.setNumberOfServers(1)
    c1 = ClosedClass(model, 'C1', 2, d, 0)   # priority 0 = more urgent
    c2 = ClosedClass(model, 'C2', 1, d, 1)
    d.setService(c1, Exp(1.0))
    d.setService(c2, Exp(2.0))
    # two service phases, so PR and PI cannot coincide
    q.setService(c1, Erlang.fitMeanAndOrder(1.0, 2))
    q.setService(c2, Erlang.fitMeanAndOrder(0.5, 2))
    model.link(Network.serialRouting(d, q))
    return model


@pytest.mark.parametrize('name', sorted(EXPECTED))
def test_preempt_family_qlen(name):
    model = _build(getattr(SchedStrategy, name))
    table = SolverCTMC(model).getAvgTable()
    qlen = np.array([float(x) for x in table['QLen']])
    np.testing.assert_allclose(qlen, EXPECTED[name], rtol=0, atol=5e-6)


def test_prio_isolates_the_urgent_class():
    # class 1 is never preempted under the FCFS PRIO variants, so its queue
    # length is the exact 13/11 of the isolated single-class model
    for name in ('FCFSPRPRIO', 'FCFSPIPRIO'):
        model = _build(getattr(SchedStrategy, name))
        qlen = np.array([float(x) for x in SolverCTMC(model).getAvgTable()['QLen']])
        assert abs(qlen[2] - 13.0 / 11.0) < 5e-6


def test_pi_differs_from_pr():
    # the whole PR-vs-PI difference is the phase a promoted job resumes in
    for pr, pi in (('FCFSPR', 'FCFSPI'), ('LCFSPR', 'LCFSPI'),
                   ('FCFSPRPRIO', 'FCFSPIPRIO'), ('LCFSPRPRIO', 'LCFSPIPRIO')):
        assert EXPECTED[pr] != EXPECTED[pi]
        qs = []
        for name in (pr, pi):
            model = _build(getattr(SchedStrategy, name))
            qs.append(np.array([float(x) for x in SolverCTMC(model).getAvgTable()['QLen']]))
        assert not np.allclose(qs[0], qs[1], atol=1e-6)
