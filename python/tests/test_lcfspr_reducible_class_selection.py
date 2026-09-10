"""Which closed communicating class an order-preserving discipline selects.

Under deterministic cyclic routing jobs cannot overtake, so the relative service
order is frozen by the initial condition: the ordered chain is REDUCIBLE and its
stationary distribution is not unique. Kelly's Theorem 1 and BCMP product form
both presume irreducibility, so LCFS-PR is NOT required to equal PS here, and it
does not -- that part is expected and is not what this file guards.

What must hold is that every engine picks the SAME class, namely the one the
declared initial state reaches. The reference is MATLAB/JAR. The failure this
pins is silent in the worst way: each engine returns a valid stationary vector of
the same generator, so nothing errors and every number looks plausible, but they
are stationary vectors of DIFFERENT closed classes. Three distinct answers were
in circulation before this was aligned (MATLAB/JAR, C++/LDES, and a python
equal-weight mixture that corresponded to no initial condition at all).

The second test is the control: on a topology where jobs CAN overtake the chain
is irreducible, and there LCFS-PR must equal PS to machine precision.
"""

import numpy as np
import pytest

from line_solver import (Network, Queue, ClosedClass, Exp, SchedStrategy,
                         SolverCTMC, SolverMVA)

MU_CYCLIC = np.array([[2.0, 3.0], [1.5, 4.0]])

# Kelly weights (product of 1/mu over the jobs) on each closed class, by hand:
# class A = 3/5, class B = 6/11. The declared order C1,C2 must select A.
CLASS_A = 0.6
CLASS_B = 6.0 / 11.0


def _cyclic(order):
    """Q1 -> Q2 -> Q1, one job per class, classes DECLARED in `order`."""
    model = Network('lcfspr_cyclic')
    q1 = Queue(model, 'Q1', SchedStrategy.LCFSPR)
    q2 = Queue(model, 'Q2', SchedStrategy.LCFSPR)
    cls = [None, None]
    for k in range(2):
        r = order[k]
        cls[r] = ClosedClass(model, 'C%d' % (r + 1), 1, q1, 0)
    for r in range(2):
        q1.setService(cls[r], Exp(MU_CYCLIC[0, r]))
        q2.setService(cls[r], Exp(MU_CYCLIC[1, r]))
    model.link(model.serialRouting(q1, q2))
    return model


def _x1(model):
    table = SolverCTMC(model).getAvgTable()
    jobclass = [str(v) for v in table['JobClass']]
    return max(float(v) for v, c in zip(table['Tput'], jobclass) if c == 'C1')


@pytest.mark.parametrize('order,expected', [((0, 1), CLASS_A), ((1, 0), CLASS_B)])
def test_declared_order_selects_the_reachable_class(order, expected):
    # flipping the declaration order flips the initial service order, hence the
    # class; both values are exact Kelly normalizations, not approximations
    assert abs(_x1(_cyclic(order)) - expected) < 1e-9


def test_the_two_classes_are_actually_distinct():
    # guards the test itself: if the model ever stopped being reducible the two
    # parametrizations above would coincide and pass vacuously
    assert abs(CLASS_A - CLASS_B) > 1e-3
    assert abs(_x1(_cyclic((0, 1))) - _x1(_cyclic((1, 0)))) > 1e-3


def test_irreducible_topology_restores_the_kelly_equivalence():
    """Q1 -> {Q2,Q3} w.p. 0.5 -> Q1: jobs can overtake, so the chain is
    irreducible and LCFS-PR must equal PS and exact MVA."""
    mu = np.array([[2.0, 3.0, 1.7], [1.5, 4.0, 2.6], [2.2, 1.9, 3.1]])

    def build(sched):
        model = Network('branch')
        q = [Queue(model, 'Q%d' % (i + 1), sched) for i in range(3)]
        cls = [ClosedClass(model, 'C%d' % (r + 1), 1, q[0], 0) for r in range(3)]
        for r in range(3):
            for i in range(3):
                q[i].setService(cls[r], Exp(mu[i, r]))
        P = model.initRoutingMatrix()
        for r in range(3):
            P.set(cls[r], cls[r], q[0], q[1], 0.5)
            P.set(cls[r], cls[r], q[0], q[2], 0.5)
            P.set(cls[r], cls[r], q[1], q[0], 1.0)
            P.set(cls[r], cls[r], q[2], q[0], 1.0)
        model.link(P)
        return model

    qlen = lambda t: np.array([float(x) for x in t['QLen']])
    ps = qlen(SolverCTMC(build(SchedStrategy.PS)).getAvgTable())
    lpr = qlen(SolverCTMC(build(SchedStrategy.LCFSPR)).getAvgTable())
    mva = qlen(SolverMVA(build(SchedStrategy.PS)).getAvgTable())
    np.testing.assert_allclose(ps, mva, rtol=0, atol=1e-10)
    np.testing.assert_allclose(lpr, ps, rtol=0, atol=1e-10)
