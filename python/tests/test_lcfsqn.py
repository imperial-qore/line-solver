"""Regression test for the specialized LCFS + LCFS-PR product-form solver
(_solver_nc_lcfsqn): class multiplicities N_r > 1 must be handled by
expanding classes into exchangeable single-job copies (returned an empty
all-zero table before the fix). Reference: exact SolverCTMC on the same
model (order-dependent product form, Casale QUESTA 2026)."""

import numpy as np
import pytest

from line_solver import (
    Network, Queue, ClosedClass, Exp, SchedStrategy, SolverCTMC, SolverNC,
)

L1 = [0.20, 0.45]
L2 = [0.30, 0.12]


def _cycle(N):
    model = Network('lcfsqn_cycle')
    q1 = Queue(model, 'Q1', SchedStrategy.LCFS)
    q2 = Queue(model, 'Q2', SchedStrategy.LCFSPR)
    classes = []
    for r in range(2):
        c = ClosedClass(model, 'C%d' % (r + 1), N[r], q1)
        q1.setService(c, Exp(1 / L1[r]))
        q2.setService(c, Exp(1 / L2[r]))
        classes.append(c)
    P = model.initRoutingMatrix()
    for c in classes:
        P.set(c, Network.serialRouting(q1, q2))
    model.link(P)
    return model


@pytest.mark.parametrize('N', [[1, 1], [2, 1], [2, 2]])
def test_lcfsqn_nc_matches_ctmc(N):
    tc = SolverCTMC(_cycle(N)).getAvgTable()
    tn = SolverNC(_cycle(N)).getAvgTable()
    assert len(tn['QLen']) == len(tc['QLen'])
    for col in ('QLen', 'Util', 'Tput'):
        err = np.max(np.abs(np.array(tn[col]) - np.array(tc[col])))
        assert err < 1e-8, 'NC %s deviates from CTMC by %g at N=%s' % (col, err, N)
