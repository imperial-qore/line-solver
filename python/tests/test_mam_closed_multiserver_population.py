"""
Regression: MAM must conserve the closed population at a multiserver station.

solver_mam_basic converges a per-chain surrogate arrival rate so that the queue
lengths sum to the chain population. Its post-loop second pass then rescaled QN
to that population, floored the response time at one full service time
(RN = max(S, QN/TN)) and restated QN = RN*TN from an untouched TN -- discarding
the rescale. With one server the floor is rarely active and the pass is a near
no-op; with c > 1 the loop's surrogate response time sits well below S, so the
floor roughly doubled RN and both QN and TN came out inflated by the same
per-class factor. A 3-station model with 10 servers and N = 2 returned
sum(QN) = 3.52 instead of 2.

Oracle: SolverCTMC. With 10 servers and 2 jobs no job ever queues, so R = S is
the exact answer and MAM must reproduce CTMC outright, not approximately.
"""

import numpy as np
import pytest

from line_solver import (Network, Delay, Queue, ClosedClass, Exp,
                         SchedStrategy, SolverMAM, SolverCTMC)

RATES = [[96.0, 60.0], [70.0, 84.0], [96.0, 1.0]]


def _cqn(nservers):
    m = Network('cqn')
    st = [Delay(m, 'Delay'),
          Queue(m, 'Queue1', SchedStrategy.PS),
          Queue(m, 'Queue2', SchedStrategy.PS)]
    jc = [ClosedClass(m, 'ClassA', 1, st[0]), ClosedClass(m, 'ClassB', 1, st[0])]
    for i in range(3):
        for r in range(2):
            st[i].setService(jc[r], Exp(RATES[i][r]))
        if i > 0:
            st[i].setNumServers(nservers)
    P = m.initRoutingMatrix()
    P.set(jc[0], jc[0], st[0], st[1], 0.6)
    P.set(jc[0], jc[0], st[0], st[2], 0.4)
    P.set(jc[0], jc[0], st[1], st[0], 1.0)
    P.set(jc[0], jc[0], st[2], st[0], 1.0)
    P.set(jc[1], jc[1], st[0], st[1], 1.0)
    P.set(jc[1], jc[1], st[1], st[0], 1.0)
    P.set(jc[1], jc[1], st[2], st[0], 1.0)
    m.link(P)
    return m


@pytest.mark.parametrize('nservers', [1, 2, 10])
def test_closed_population_is_conserved(nservers):
    QN = SolverMAM(_cqn(nservers), verbose=False).getAvg()[0]
    # one job per class, and each class is its own chain
    assert np.allclose(np.sum(QN, axis=0), [1.0, 1.0], atol=1e-6)


def test_multiserver_matches_ctmc_exactly():
    # 10 servers, 2 jobs: no job ever queues, so MAM's R = S floor is exact
    m = _cqn(10)
    QN, UN, RN, TN = SolverMAM(m, verbose=False).getAvg()[:4]
    Qc, Uc, Rc, Tc = SolverCTMC(_cqn(10), verbose=False).getAvg()[:4]
    assert np.allclose(QN, Qc, atol=1e-9)
    assert np.allclose(RN, Rc, atol=1e-9)
    assert np.allclose(TN, Tc, rtol=1e-9)
