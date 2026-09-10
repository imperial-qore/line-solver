"""A closed job that finds no room must BLOCK, and the NRM engine must say so.

Until 2026-08-19 it did not (BUG-81). ``_capacity_loss`` correctly declined to
DROP a closed job -- a closed network's population is an invariant -- but
nothing then stopped the reaction, so the firing went into the full station
anyway and the engine returned the UNCONSTRAINED product-form answer: on the
model below, QLen [1.96 2.07 1.97] against the exact [3.609 0.971 1.420], a mean
of 2.07 at a station that holds 2. The serial engine, whose producer
(``after_event_station``) already implements the open/closed contract, was right
the whole time, so the two SSA methods disagreed with each other and the default
method (``nrm``) was the wrong one.

Population conservation alone does not catch this: the broken engine conserved
it too. The oracle is the exact CTMC, plus the arithmetic fact that a station
capped at 2 cannot hold 2.07 on average.
"""

import numpy as np
import pytest

from line_solver import (Network, Queue, Source, Sink, ClosedClass, OpenClass,
                         Exp, SchedStrategy, DropStrategy, SolverSSA, SolverCTMC)

SAMPLES = 50000
SEED = 23000


def _capped_tandem(cap=2, n=6):
    """Closed 3-queue tandem, Exp(1) FCFS everywhere, Q2 capped at `cap`."""
    model = Network('tandem')
    q1 = Queue(model, 'Q1', SchedStrategy.FCFS)
    q2 = Queue(model, 'Q2', SchedStrategy.FCFS)
    q3 = Queue(model, 'Q3', SchedStrategy.FCFS)
    c1 = ClosedClass(model, 'C1', n, q1, 0)
    for q in (q1, q2, q3):
        q.setService(c1, Exp(1.0))
    q2.setClassCapacity(c1, cap)
    model.link(Network.serialRouting(q1, q2, q3))
    return model


def _mm1k(K):
    """M/M/1/K, lambda 0.8, mu 1, explicit DROP -- the open half of the contract."""
    model = Network('mm1k')
    src = Source(model, 'Src')
    q = Queue(model, 'Q', SchedStrategy.FCFS)
    snk = Sink(model, 'Snk')
    oc = OpenClass(model, 'O')
    src.setArrival(oc, Exp(0.8))
    q.setService(oc, Exp(1.0))
    q.setClassCapacity(oc, K)
    q.setDropRule(oc, DropStrategy.DROP)
    model.link(Network.serialRouting(src, q, snk))
    return model


def _ssa(model, method):
    s = SolverSSA(model, method=method, samples=SAMPLES, seed=SEED, verbose=False)
    return (np.asarray(s.getAvgQLen()).flatten(),
            np.asarray(s.getAvgTput()).flatten())


def test_nrm_respects_a_binding_closed_class_capacity():
    """The capped station's mean cannot exceed its cap, and did (2.07 > 2)."""
    model = _capped_tandem(cap=2)
    q, _ = _ssa(model, 'nrm')
    assert q[1] <= 2.0, 'Q2 holds at most 2 jobs, so its mean cannot be %g' % q[1]


def test_nrm_matches_the_exact_chain_on_a_capped_closed_model():
    """Against SolverCTMC on the same model: the constrained answer, not the free one."""
    model = _capped_tandem(cap=2)
    exact = np.asarray(SolverCTMC(model).getAvgQLen()).flatten()
    xexact = float(np.asarray(SolverCTMC(model).getAvgTput()).flatten()[0])
    q, t = _ssa(model, 'nrm')
    np.testing.assert_allclose(q, exact, atol=0.15)
    # the departure-rate integral must net out the firings that were blocked;
    # uncorrected it read 0.744 against 0.652 here.
    assert t[0] == pytest.approx(xexact, abs=0.05)


def test_the_two_ssa_engines_agree_on_a_capped_closed_model():
    """They disagreed by a factor of 2 at Q1: the defect was engine-specific."""
    model = _capped_tandem(cap=2)
    q_nrm, _ = _ssa(model, 'nrm')
    q_ser, _ = _ssa(model, 'serial')
    np.testing.assert_allclose(q_nrm, q_ser, atol=0.15)


def test_population_is_conserved():
    """Necessary, not sufficient -- the broken engine passed this too."""
    model = _capped_tandem(cap=2)
    q, _ = _ssa(model, 'nrm')
    assert q.sum() == pytest.approx(6.0, abs=1e-9)


@pytest.mark.parametrize('K,exact', [(1, 4.0 / 9.0), (2, 0.852459016393443)])
def test_open_class_drop_is_unchanged(K, exact):
    """The regression guard: an OPEN arrival that finds no room is still LOST.

    The blocking gate is keyed on the class type, so M/M/1/K must keep the exact
    loss behaviour it had; if it starts blocking instead, the source stops
    offering and these means move.
    """
    q, _ = _ssa(_mm1k(K), 'nrm')
    # station 0 is the Source, which holds no jobs; the Queue is station 1
    assert float(q[1]) == pytest.approx(exact, abs=0.02)
