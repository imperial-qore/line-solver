"""SolverMVA on the size-based M/G/1 disciplines, against MATLAB.

MATLAB `solver_mva_qsys_sizebased_analyzer.m` evaluates the Wierman and
Harchol-Balter (SIGMETRICS 2003) response times. Two open classes at a single
Queue: C1 is Exp(1) at rate 0.3, C2 is Erlang(mean 2, SCV 0.5) at rate 0.2, so
rho = 0.7 and the classes differ in both mean size and variability -- which is
what the size-based orderings act on.
"""

import numpy as np
import pytest

from line_solver import (Erlang, Exp, Network, OpenClass, Queue, SchedStrategy,
                         Sink, SolverMVA, Source)


def _build(sched):
    m = Network('sb')
    src = Source(m, 'Source')
    q = Queue(m, 'Q', sched)
    snk = Sink(m, 'Sink')
    c1 = OpenClass(m, 'C1')
    c2 = OpenClass(m, 'C2')
    src.setArrival(c1, Exp(0.3))
    src.setArrival(c2, Exp(0.2))
    q.setService(c1, Exp(1.0))
    q.setService(c2, Erlang.fitMeanAndSCV(2.0, 0.5))
    P = m.initRoutingMatrix()
    P.set(c1, c1, Network.serialRouting(src, q, snk))
    P.set(c2, c2, Network.serialRouting(src, q, snk))
    m.link(P)
    return m


# MATLAB SolverMVA on the same model: (QLen C1, QLen C2, RespT C1, RespT C2)
REFERENCE = {
    'SRPT': (0.4985567765823074, 0.8212546141301067,
             1.6618559219410247, 4.1062730706505333),
    'PSJF': (0.61224489795918369, 3.333333333333333,
             2.0408163265306123, 16.666666666666664),
    'FB': (0.63587346812077949, 2.1712143660569709,
           2.1195782270692649, 10.856071830284854),
    'LRPT': (2.1333333333333333, 0.66666666666666674,
             7.1111111111111107, 3.3333333333333335),
    'SETF': (1.2256856111717471, 2.8758520284392772,
             4.0856187039058236, 14.379260142196385),
}


@pytest.mark.parametrize('name', sorted(REFERENCE.keys()))
def test_sizebased_matches_matlab(name):
    q1, q2, r1, r2 = REFERENCE[name]
    solver = SolverMVA(_build(getattr(SchedStrategy, name)))
    res = solver.getAvg()
    QN, UN, RN, TN = res[0], res[1], res[2], res[3]
    assert QN[1, 0] == pytest.approx(q1, abs=1e-9)
    assert QN[1, 1] == pytest.approx(q2, abs=1e-9)
    assert RN[1, 0] == pytest.approx(r1, abs=1e-9)
    assert RN[1, 1] == pytest.approx(r2, abs=1e-9)
    # discipline-independent: rho_k = lambda_k/mu_k, and the flow is lossless
    assert UN[1, 0] == pytest.approx(0.3, abs=1e-12)
    assert UN[1, 1] == pytest.approx(0.4, abs=1e-12)
    assert TN[1, 0] == pytest.approx(0.3, abs=1e-12)
    assert TN[1, 1] == pytest.approx(0.2, abs=1e-12)
    # Little's law at the station
    assert QN[1, 0] == pytest.approx(TN[1, 0] * RN[1, 0], abs=1e-9)
    assert QN[1, 1] == pytest.approx(TN[1, 1] * RN[1, 1], abs=1e-9)


def test_the_disciplines_are_not_interchangeable():
    # SRPT favours the short class and LRPT the long one: a size-blind analyzer
    # would return the same numbers for both.
    srpt = SolverMVA(_build(SchedStrategy.SRPT)).getAvg()
    lrpt = SolverMVA(_build(SchedStrategy.LRPT)).getAvg()
    assert srpt[2][1, 0] < lrpt[2][1, 0]
    assert srpt[2][1, 1] > lrpt[2][1, 1]


def test_lrpt_leaf_uses_the_general_branch_when_cs_differs():
    # qsys_mg1_lrpt has TWO branches; only the exponential one was ported, so a
    # class with cs != 1 was answered with the exponential formula.
    from line_solver.api.qsys import qsys_mg1_lrpt
    W, _ = qsys_mg1_lrpt([0.3, 0.2], [1.0, 0.5], [1.0, np.sqrt(0.5)])
    assert W[0] == pytest.approx(7.1111111111111107, abs=1e-9)
    assert W[1] == pytest.approx(3.3333333333333335, abs=1e-9)
