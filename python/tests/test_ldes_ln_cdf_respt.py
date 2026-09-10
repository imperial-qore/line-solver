"""SolverLDES on a LayeredNetwork: the measured per-entry response time law.

The analytical twin is ``SolverLN.getCdfRespT`` (``tests/test_ln_cdf_respt.py``),
which fits an APH to three moments of each term and convolves. Here the engine
timed every invocation from the instant the request reached the entry to its
reply -- the interval ``RLN`` averages -- so the law's mean reproduces that row
exactly and its tail is observed rather than extrapolated.

The whole chain is exercised, not a stub of it: ``save_model`` writes the
layered ``line-model``, ``jline.cli.LdesCLI`` simulates it, and
``_parse_ln_result_json`` reads the layered ``ldes-result`` back.
"""
import numpy as np
import pytest

from line_solver import (Activity, Entry, Exp, LayeredNetwork, Processor,
                         SchedStrategy, Task)
from line_solver.solvers.wrappers.solver_ldes import SolverLDES


def two_task_lqn():
    """T1 -> T2, one call level, so E2 is driven and E1 waits on it."""
    m = LayeredNetwork('ldescdf')
    p1 = Processor(m, 'P1', 1, SchedStrategy.PS)
    p2 = Processor(m, 'P2', 1, SchedStrategy.PS)
    t1 = Task(m, 'T1', 3, SchedStrategy.REF).on(p1).setThinkTime(Exp(1))
    t2 = Task(m, 'T2', 2, SchedStrategy.FCFS).on(p2)
    e1 = Entry(m, 'E1').on(t1)
    e2 = Entry(m, 'E2').on(t2)
    Activity(m, 'A1', Exp(2)).on(t1).boundTo(e1).synchCall(e2, 1)
    Activity(m, 'A2', Exp(3)).on(t2).boundTo(e2).repliesTo(e2)
    return m


@pytest.fixture(scope='module')
def solved():
    """One simulation, reused: the run is the expensive part, not the getters."""
    model = two_task_lqn()
    solver = SolverLDES(model, samples=50000, seed=23000)
    return model, solver, model.getStruct()


def ecdf_mean(cdf):
    """Mean of an [F(t), t] table, from its jumps."""
    F, t = np.asarray(cdf)[:, 0], np.asarray(cdf)[:, 1]
    return float(np.sum(np.diff(np.concatenate(([0.0], F))) * t))


def test_every_entry_carries_an_empirical_law(solved):
    _, solver, lsn = solved
    RD = solver.getCdfRespTLN()
    assert len(RD) == lsn.nentries
    for e, cdf in enumerate(RD):
        assert cdf is not None, 'entry %d observed nothing' % e
        assert cdf.shape[1] == 2, 'the table is [F(t), t], two columns'
        assert cdf.shape[0] > 1
        F, t = cdf[:, 0], cdf[:, 1]
        assert np.all(np.diff(t) >= 0), 'the time column must be nondecreasing'
        assert np.all(np.diff(F) >= -1e-12), 'the CDF column must be nondecreasing'
        assert F[-1] == pytest.approx(1.0), 'an ecdf reaches 1 at its last observation'
        assert t[0] > 0.0, 'a response time is strictly positive'


def test_the_law_mean_reproduces_the_reported_entry_response_time(solved):
    """The samples must be taken over the interval the table averages.

    This is the check that the two are the same measurement and not merely
    similar ones: a law drawn from a different interval agrees to simulation
    noise, not to machine precision.
    """
    _, solver, lsn = solved
    RD = solver.getCdfRespTLN()
    table = solver.getLNAvgTable()
    for e in range(lsn.nentries):
        reported = float(table.RespT[lsn.eshift + e])
        assert ecdf_mean(RD[e]) == pytest.approx(reported, rel=1e-12)


def test_the_layered_table_carries_the_lqn_mask(solved):
    """A processor has no queue, response time or throughput of its own, and a
    task no response time; that mask is shared with SolverLN and the JAR so the
    three tables compare cell for cell."""
    _, solver, lsn = solved
    table = solver.getLNAvgTable()
    assert len(table) == lsn.nidx
    assert list(table.columns) == ['Node', 'NodeType', 'QLen', 'Util', 'RespT',
                                   'ResidT', 'ArvR', 'Tput']
    for _, row in table.iterrows():
        if row.NodeType == 'Processor':
            assert np.isnan(row.QLen) and np.isnan(row.RespT) and np.isnan(row.Tput)
            assert not np.isnan(row.Util)
        elif row.NodeType in ('Task', 'RefTask'):
            assert np.isnan(row.RespT)
            assert not np.isnan(row.QLen)
        elif row.NodeType == 'Entry':
            assert not np.isnan(row.RespT)
        assert np.isnan(row.ArvR), 'ArvR is not measured on the layered path'


def test_ensemble_avg_spans_the_lqn_index_space(solved):
    _, solver, lsn = solved
    vectors = solver.getEnsembleAvg()
    assert len(vectors) == 5
    for v in vectors:
        assert v.shape == (lsn.nidx,)
    RLN = vectors[2]
    table = solver.getLNAvgTable()
    for e in range(lsn.nentries):
        assert RLN[lsn.eshift + e] == pytest.approx(float(table.RespT[lsn.eshift + e]))


def test_station_indexed_getters_refuse_a_layered_model(solved):
    """They used to die on a missing attribute of LNLDESResult, which named a
    private field instead of the mistake."""
    _, solver, _ = solved
    with pytest.raises(RuntimeError, match='LQN element'):
        solver.getAvgTable()
    with pytest.raises(RuntimeError, match='LQN element'):
        solver.getAvg()
    with pytest.raises(RuntimeError, match='getCdfRespTLN'):
        solver.getTranCdfRespT()


def test_get_cdf_respt_dispatches_a_layered_model(solved):
    """MATLAB @SolverLDES/getCdfRespT.m hands a LayeredNetwork to
    getCdfRespTLN and returns the per-entry laws; so does this port."""
    _, solver, lsn = solved
    RD = solver.getCdfRespT()
    RDLN = solver.getCdfRespTLN()
    assert len(RD) == lsn.nentries
    for e in range(lsn.nentries):
        if RDLN[e] is None:
            assert RD[e] is None
        else:
            assert np.array_equal(RD[e], RDLN[e])


def test_the_entry_getter_refuses_a_flat_model():
    from line_solver import Network, OpenClass, Queue, Sink, Source
    m = Network('flat')
    source = Source(m, 'Source')
    queue = Queue(m, 'Queue', SchedStrategy.FCFS)
    sink = Sink(m, 'Sink')
    cls = OpenClass(m, 'Class1')
    source.setArrival(cls, Exp(1))
    queue.setService(cls, Exp(2))
    m.link(m.serialRouting(source, queue, sink))
    with pytest.raises(RuntimeError, match='LayeredNetwork'):
        SolverLDES(m, samples=1000, seed=23000).getCdfRespTLN()
