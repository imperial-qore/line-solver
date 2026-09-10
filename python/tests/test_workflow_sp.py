"""
Closed-form checks of the series-parallel composition of Workflow and of the
geometric loop. Mirrors jar/src/test/java/jline/lang/workflow/WorkflowSPTest.java
and the MATLAB Workflow checks.

Copyright (c) 2012-2026, QORE Lab, Imperial College London
All rights reserved.
"""
import numpy as np
import pytest

from line_solver.lang.workflow import Workflow
from line_solver.distributions import Exp, APH

TOL = 1e-9


def moments(alpha, T):
    """Mean and SCV of the PH law (alpha, T)."""
    alpha = np.asarray(alpha, dtype=float).reshape(1, -1)
    T = np.asarray(T, dtype=float)
    e = np.ones((T.shape[0], 1))
    Ti = np.linalg.inv(T)
    m1 = float((-alpha @ Ti @ e).item())
    m2 = float((2 * alpha @ Ti @ Ti @ e).item())
    return m1, (m2 - m1 ** 2) / m1 ** 2


def emax_exp(m1, m2):
    """E[max] of two independent exponentials of means m1, m2."""
    r1, r2 = 1.0 / m1, 1.0 / m2
    return 1 / r1 + 1 / r2 - 1 / (r1 + r2)


def fork_in_loop(mean_c=2.0):
    wf = Workflow('ForkInLoop')
    for nm, v in [('A', 0.5), ('B', 1.0), ('C', mean_c), ('D', 1.5),
                  ('E', 0.25), ('F', 0.75)]:
        wf.addActivity(nm, Exp.fitMean(v))
    wf.addPrecedence(Workflow.Loop('A', ['B', 'F'], 2))
    wf.addPrecedence(Workflow.AndFork('B', ['C', 'D']))
    wf.addPrecedence(Workflow.AndJoin(['C', 'D'], 'E'))
    return wf


def test_serial_of_two_exponentials_is_erlang2():
    wf = Workflow('serial2')
    wf.addActivity('A', Exp.fitMean(1.0))
    wf.addActivity('B', Exp.fitMean(1.0))
    wf.addPrecedence(Workflow.Serial('A', 'B'))

    m, c = moments(*wf.toPH())
    assert m == pytest.approx(2.0, abs=TOL)
    assert c == pytest.approx(0.5, abs=TOL)
    assert wf.getSPTree() is not None


def test_parallel_is_the_maximum_of_the_branches():
    ra, rb = 0.5, 1.0 / 1.5
    alpha, T = Workflow._composeParallel(np.array([[1.0]]), np.array([[-ra]]),
                                         np.array([[1.0]]), np.array([[-rb]]))
    m, _ = moments(alpha, T)
    assert m == pytest.approx(1 / ra + 1 / rb - 1 / (ra + rb), abs=TOL)


def test_geometric_loop_of_an_exponential_is_exponential():
    mean, count = 2.0, 3.0
    alpha, T = Workflow._composeLoopGeometric(np.array([[1.0]]),
                                              np.array([[-1 / mean]]), count)
    m, c = moments(alpha, T)
    assert m == pytest.approx(count * mean, abs=TOL)
    assert c == pytest.approx(1.0, abs=TOL)


def test_geometric_loop_matches_the_compound_moments():
    # Erlang-2 body of mean 2 and SCV 1/2
    alpha_b = np.array([[1.0, 0.0]])
    T_b = np.array([[-1.0, 1.0], [0.0, -1.0]])
    count, scv_body = 3.0, 0.5

    alpha, T = Workflow._composeLoopGeometric(alpha_b, T_b, count)
    m, c = moments(alpha, T)

    assert m == pytest.approx(count * 2.0, abs=TOL)
    assert c == pytest.approx(scv_body / count + 1 - 1 / count, abs=TOL)
    # the order is that of the body, and the generator is cyclic
    assert T.shape[0] == 2
    assert not Workflow.isAcyclicGenerator(T)


def test_fractional_loop_count_runs_the_body_with_that_probability():
    alpha, T = Workflow._composeLoopGeometric(np.array([[1.0]]),
                                              np.array([[-0.5]]), 0.25)
    m, _ = moments(alpha, T)
    assert m == pytest.approx(0.25 * 2.0, abs=1e-6)


def test_loop_workflow_keeps_the_mean_and_the_body_order():
    wf = Workflow('LoopWorkflow')
    wf.addActivity('A', Exp.fitMean(1.0))
    wf.addActivity('B', Exp.fitMean(2.0))
    wf.addActivity('C', Exp.fitMean(0.5))
    wf.addPrecedence(Workflow.Loop('A', ['B', 'C'], 3.0))

    alpha, T = wf.toPH()
    m, c = moments(alpha, T)
    assert m == pytest.approx(1.0 + 3 * 2.0 + 0.5, abs=TOL)
    # A, the geometric loop over B, and C
    assert T.shape[0] == 3
    var_tot = 1.0 + 36.0 + 0.25
    assert c == pytest.approx(var_tot / 7.5 ** 2, abs=TOL)


def test_or_fork_mixes_the_branches():
    wf = Workflow('BranchingWorkflow')
    for nm, v in [('A', 1.0), ('B', 2.0), ('C', 5.0), ('D', 0.5)]:
        wf.addActivity(nm, Exp.fitMean(v))
    wf.addPrecedence(Workflow.OrFork('A', ['B', 'C'], np.array([0.6, 0.4])))
    wf.addPrecedence(Workflow.OrJoin(['B', 'C'], 'D'))

    m, _ = moments(*wf.toPH())
    assert m == pytest.approx(1.0 + 0.6 * 2.0 + 0.4 * 5.0 + 0.5, abs=TOL)
    assert wf.getSPTree() is not None


def test_fork_nested_inside_a_loop_is_reduced_exactly():
    wf = fork_in_loop()
    m, _ = moments(*wf.toPH())
    body = 1.0 + emax_exp(2.0, 1.5) + 0.25
    assert m == pytest.approx(0.5 + 2 * body + 0.75, abs=TOL)
    assert wf.getSPTree() is not None


def test_incremental_refresh_matches_a_rebuilt_workflow():
    wf = fork_in_loop()
    wf.toPH()

    ref = fork_in_loop(mean_c=3.0)
    m_ref, c_ref = moments(*ref.toPH())
    n_ref = ref.toPH()[1].shape[0]

    wf.setActivityDemand('C', Exp.fitMean(3.0))
    alpha, T = wf.refreshPH()
    m_inc, c_inc = moments(alpha, T)

    assert m_inc == pytest.approx(m_ref, abs=TOL)
    assert c_inc == pytest.approx(c_ref, abs=TOL)
    assert T.shape[0] == n_ref


def test_mean_only_rescale_keeps_the_shape():
    wf = Workflow('Rescale')
    wf.addActivity('A', APH.fitMeanAndSCV(2.0, 0.3))
    wf.addActivity('B', Exp.fitMean(1.0))
    wf.addPrecedence(Workflow.Serial('A', 'B'))
    wf.toPH()

    wf.setActivityDemandMean('A', 5.0)
    alpha, T = wf.refreshPH()
    m, c = moments(alpha, T)

    ref = Workflow('RescaleRef')
    ref.addActivity('A', APH.fitMeanAndSCV(5.0, 0.3))
    ref.addActivity('B', Exp.fitMean(1.0))
    ref.addPrecedence(Workflow.Serial('A', 'B'))
    alpha_r, T_r = ref.toPH()
    _, c_ref = moments(alpha_r, T_r)

    assert m == pytest.approx(6.0, abs=TOL)
    assert c == pytest.approx(c_ref, abs=1e-6)
    assert T.shape[0] == T_r.shape[0]


def test_quorum_join_is_refused_by_name():
    wf = Workflow('Quorum')
    for nm in ['A', 'B', 'C', 'D']:
        wf.addActivity(nm, Exp.fitMean(1.0))
    wf.addPrecedence(Workflow.AndFork('A', ['B', 'C']))
    wf.addPrecedence(Workflow.AndJoin(['B', 'C'], 'D', 1))

    with pytest.raises(ValueError, match='quorum'):
        wf.toPH()


def test_full_and_join_is_accepted():
    wf = Workflow('FullJoin')
    for nm in ['A', 'B', 'C', 'D']:
        wf.addActivity(nm, Exp.fitMean(1.0))
    wf.addPrecedence(Workflow.AndFork('A', ['B', 'C']))
    wf.addPrecedence(Workflow.AndJoin(['B', 'C'], 'D'))

    m, _ = moments(*wf.toPH())
    assert m == pytest.approx(1.0 + emax_exp(1.0, 1.0) + 1.0, abs=TOL)


def test_a_graph_that_is_not_series_parallel_falls_back():
    wf = Workflow('NotSP')
    for nm in ['A', 'B', 'C']:
        wf.addActivity(nm, Exp.fitMean(1.0))
    wf.addPrecedence(Workflow.Serial('A', 'B'))
    wf.addPrecedence(Workflow.Serial('A', 'C'))

    m, _ = moments(*wf.toPH())
    assert np.isfinite(m) and m > 0
    assert wf.getSPTree() is None


def test_execution_counts_weight_loop_and_branch():
    wf = Workflow('Execs')
    for nm in ['A', 'B', 'C']:
        wf.addActivity(nm, Exp.fitMean(1.0))
    wf.addPrecedence(Workflow.Loop('A', ['B', 'C'], 3.0))

    tree = wf.getSPTree()
    assert tree is not None
    execs = tree['execs']
    leaf_of = tree['leaf_of']
    assert execs[leaf_of[0]] == pytest.approx(1.0)  # A
    assert execs[leaf_of[1]] == pytest.approx(3.0)  # B, the loop body
    assert execs[leaf_of[2]] == pytest.approx(1.0)  # C
