"""A capacity set AFTER the struct was materialized must still be honoured.

`Network` caches the compiled `sn` and rebuilds it only when the cache is
invalidated. `sn.cap`/`sn.classcap` are DERIVED from the node's declared
capacity, so a `setCapacity` that does not invalidate is not applied late -- it
is DROPPED, and every sn-reading solver then answers the UNBOUNDED model with no
error and a plausible number.

The model here is the smallest one that makes the difference visible. Two FCFS
queues in a closed cycle, mu = 1 and 0.8, N = 2, with a buffer of 1 at the second
queue. Capped, the chain has two states and is solvable by hand:

    (2,0) --1.0--> (1,1) --0.8--> (2,0)
    pi = [0.8, 1.0]/1.8, X = 0.8*pi(1,1) = 4/9, Q = [13/9, 5/9]

Uncapped it is the product form over D = [1, 1.25], X = 2.25/3.8125 = 0.590164
and Q2 = 1.1475 -- MORE JOBS THAN THE BUFFER HOLDS, which is what the regression
looked like before the setters invalidated.

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

import numpy as np
import pytest

from line_solver import (Network, Queue, ClosedClass, SchedStrategy, Exp,
                         SolverCTMC, SolverAG)

# The capped chain, exactly.
Q_EXACT = np.array([13.0 / 9.0, 5.0 / 9.0])
X_EXACT = 4.0 / 9.0
# The unbounded answer the stale struct used to return.
Q_UNBOUNDED = np.array([3.25 / 3.8125, 4.375 / 3.8125])


def _cycle(cap=None, cap_late=False):
    model = Network('Cycle')
    q1 = Queue(model, 'Queue1', SchedStrategy.FCFS)
    q2 = Queue(model, 'Queue2', SchedStrategy.FCFS)
    cls = ClosedClass(model, 'Class1', 2, q1)
    q1.setService(cls, Exp(1.0))
    q2.setService(cls, Exp(0.8))
    if cap is not None and not cap_late:
        q2.setCapacity(cap)
    model.link(Network.serialRouting(q1, q2))
    if cap is not None and cap_late:
        model.getStruct()      # materialize the cache FIRST
        q2.setCapacity(cap)    # the setter under test
    return model


@pytest.mark.parametrize('cap_late', [False, True])
def test_capacity_reaches_the_struct(cap_late):
    """Declared before or after the first getStruct(), the buffer is the same."""
    sn = _cycle(cap=1, cap_late=cap_late).getStruct()
    assert sn.cap.ravel()[1] == 1.0
    assert sn.classcap.ravel()[1] == 1.0


@pytest.mark.parametrize('cap_late', [False, True])
def test_ctmc_honours_a_late_capacity(cap_late):
    """The exact solver must answer the CAPPED chain either way."""
    solver = SolverCTMC(_cycle(cap=1, cap_late=cap_late))
    QN = np.asarray(solver.getAvgQLen()).ravel()
    TN = np.asarray(solver.getAvgTput()).ravel()
    assert np.allclose(QN, Q_EXACT, atol=1e-9)
    assert np.allclose(TN, X_EXACT, atol=1e-9)
    # and specifically NOT the unbounded answer, which is what a stale struct gave
    assert not np.allclose(QN, Q_UNBOUNDED, atol=1e-4)
    # a station never holds more than its buffer
    assert QN[1] <= 1.0 + 1e-9


def test_number_of_servers_reaches_the_struct():
    """The same invalidation covers the rest of the family."""
    model = _cycle()
    model.getStruct()
    model.getNodeByName('Queue2').setNumberOfServers(3)
    assert model.getStruct().nservers.ravel()[1] == 3.0


def test_ag_refuses_a_binding_capacity():
    """RCAT has no representation of a finite buffer, so it must refuse it.

    Before the gate it returned the same numbers with and without the cap.
    """
    with pytest.raises(Exception) as excinfo:
        SolverAG(_cycle(cap=1)).getAvgTable()
    assert 'apacity' in str(excinfo.value)


def test_ag_still_solves_the_uncapped_model():
    """Only a capacity that can BIND is refused."""
    QN = np.asarray(SolverAG(_cycle()).getAvgQLen()).ravel()
    assert np.isfinite(QN).all() and QN.sum() > 0
