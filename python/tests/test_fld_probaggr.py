"""SolverFLD getProbAggr under the second-order moment closure.

Every expected value is the answer of the MATLAB reference (`SolverFLD` with the
same method on the same model), which is ground truth for this solver. The
rectangle integral is evaluated by a DETERMINISTIC lattice rule in all four
codebases, so these are exact agreements, not statistical ones.
"""

import numpy as np
import pytest

from line_solver import (ClosedClass, Delay, Exp, GlobalConstants, Network, OpenClass, Queue,
                         SchedStrategy, Sink, SolverFLD, Source, VerboseLevel)

GlobalConstants.setVerbose(VerboseLevel.SILENT)

# MATLAB @SolverFLD/getProbAggr, method='minnormal', on Delay(1) -> Queue(PS, 0.8), N=4.
# RE-RECORDED 2026-09-01 against MATLAB R2026a, to full precision. The previous
# six-digit row [0.016942, 0.100103, 0.281035, 0.351526, 0.250395] was measured
# BEFORE the closure alternation and its inner mean solve were tightened from
# CoarseTol to mom_tol = 1e-6 (solver_fluid_moments.m, minnormal.py,
# fluid_moments.h), and that tightening MOVES the converged answer by ~7e-6 --
# more than the 1e-6 this test asserts at. These are MATLAB's numbers, not
# python's: native python reproduces them to 3.1e-8, which is what licenses the
# re-record. Cf. the same re-pin of MinNormalTest.testOpenMm1 in the JAR.
CELL_EXPECTED = [0.0169393390009142, 0.10009822404675858, 0.28103432766464359,
                 0.35153239268631958, 0.2503958284610216]


def _closed_ps(n_at_queue):
    """Delay(mean 1) <-> Queue(PS, mean 0.8), N=4, with n_at_queue jobs at the queue."""
    N = 4
    model = Network('probaggr')
    delay = Delay(model, 'Delay')
    queue = Queue(model, 'Q1', SchedStrategy.PS)
    cclass = ClosedClass(model, 'C1', N, delay, 0)
    delay.setService(cclass, Exp.fitMean(1.0))
    queue.setService(cclass, Exp.fitMean(0.8))
    model.link(model.serialRouting(delay, queue))
    model.initFromMarginal(np.array([[N - n_at_queue], [n_at_queue]], dtype=float))
    return model


def test_gaussian_cell_matches_matlab():
    """The closure answers from its covariance: the cell of the joint normal."""
    for n, expected in enumerate(CELL_EXPECTED):
        _, prob = SolverFLD(_closed_ps(n), method='minnormal').getProbAggr(2)
        assert prob == pytest.approx(expected, abs=1e-6), 'cell probability at n=%d' % n


def test_folded_cells_tile_the_state_space():
    """The two ends are extended to infinity, so the cells sum to one."""
    total = sum(SolverFLD(_closed_ps(n), method='minnormal').getProbAggr(2)[1]
                for n in range(5))
    assert total == pytest.approx(1.0, abs=1e-6)


def test_open_class_keeps_the_product_form():
    """An open class keeps the EXACT first-order product form, not the normal.

    On M/M/1 at rho = 0.5 the empty queue has probability 0.5; the Gaussian cell
    would return 0.391, so replacing the product form would lose accuracy.
    """
    model = Network('open')
    src = Source(model, 'Source')
    queue = Queue(model, 'Q1', SchedStrategy.PS)
    sink = Sink(model, 'Sink')
    oclass = OpenClass(model, 'C1')
    src.setArrival(oclass, Exp.fitMean(2.0))
    queue.setService(oclass, Exp.fitMean(1.0))
    model.link(model.serialRouting(src, queue, sink))
    for n, expected in enumerate([0.5, 0.25, 0.125, 0.0625]):
        model.initFromMarginal(np.array([[0], [n]], dtype=float))
        # lang is named: the product form is read off the solved rho, and the
        # three fluid ports integrate the ODE with different solvers, so the JAR
        # lands 1.8e-9 from MATLAB here -- a documented cross-codebase deviation
        # that this exactness assertion has no room for.
        _, prob = SolverFLD(model, method='minnormal', lang='python').getProbAggr(2)
        assert prob == pytest.approx(expected, abs=1e-9), 'M/M/1 marginal at n=%d' % n


def test_first_order_method_keeps_the_binomial():
    """Without a second moment there is no joint law: `closing` stays on Schmidt."""
    # lang named for the same reason as the M/M/1 marginal above: Schmidt's
    # binomial raises the mean queue length to the 4th power, so the JAR's
    # 1e-6 ODE deviation from MATLAB arrives here multiplied by four.
    _, prob = SolverFLD(_closed_ps(4), method='closing', lang='python').getProbAggr(2)
    assert prob == pytest.approx(0.223404, abs=1e-6)
