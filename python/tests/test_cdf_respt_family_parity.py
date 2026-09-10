"""
The get*Cdf* family serves the same support matrix as the MATLAB reference.

SolverMAM answers getCdfRespT with the exact matrix-analytic passage-time law
(it used to fit an exponential to the mean); SolverSSA refuses by name, as
@SolverSSA/getCdfRespT.m does; SolverCTMC's getTranCdfRespT refuses instead of
returning the steady-state law under a transient name; and the base
NetworkSolver carries MATLAB's exponential fallback for solvers with no
distributional result of their own.
"""

import shutil
import warnings

import numpy as np
import pytest

from line_solver import Exp, Network, OpenClass, Queue, SchedStrategy, Sink, Source


def _mm1(mu=2.0, lam=1.0, sched=SchedStrategy.FCFS):
    m = Network('mm1cdf')
    src = Source(m, 'Source')
    q = Queue(m, 'Queue', sched)
    snk = Sink(m, 'Sink')
    k = OpenClass(m, 'Class1')
    src.setArrival(k, Exp(lam))
    q.setService(k, Exp(mu))
    m.link(m.serialRouting(src, q, snk))
    return m


def test_mam_cdf_respt_is_the_exact_mm1_law():
    from line_solver.solvers.solver_mam import SolverMAM
    s = SolverMAM(_mm1())
    s.runAnalyzer()
    RD = s.getCdfRespT()
    assert len(RD) == 1
    cell = RD[0]
    assert cell['station'] == 2 and cell['class'] == 1
    t = np.asarray(cell['t'])
    p = np.asarray(cell['p'])
    # M/M/1 FCFS sojourn: F(t) = 1 - exp(-(mu - lambda) t)
    exact = 1.0 - np.exp(-(2.0 - 1.0) * t)
    assert np.max(np.abs(p - exact)) < 1e-12


def test_mam_cdf_respt_ps_is_not_the_fcfs_law():
    from line_solver.solvers.solver_mam import SolverMAM
    s = SolverMAM(_mm1(sched=SchedStrategy.PS))
    s.runAnalyzer()
    RD = s.getCdfRespT()
    assert len(RD) == 1
    t = np.asarray(RD[0]['t'])
    p = np.asarray(RD[0]['p'])
    assert np.all(np.diff(p) >= -1e-12)
    # The PS sojourn law of Masuyama-Takine has a heavier tail than the
    # exponential FCFS law with the same mean; it must not coincide with it
    fcfs = 1.0 - np.exp(-(2.0 - 1.0) * t[1:])
    assert np.max(np.abs(p[1:] - fcfs)) > 1e-3


def test_mam_priority_law_reproduces_the_per_class_means():
    """Distinct priorities under HOL: the MMAPPH1NPPR sojourn law, tabulated.

    The mean read off each class's tabulated CDF must reproduce the RN row the
    (independent) dec.source analyzer reports, which is the check that the
    priority flip on the way into BuTools mapped the classes back correctly.
    """
    from line_solver.solvers.solver_mam import SolverMAM
    m = Network('prio')
    src = Source(m, 'Source')
    q = Queue(m, 'Queue', SchedStrategy.HOL)
    snk = Sink(m, 'Sink')
    k1 = OpenClass(m, 'C1', 1)
    k2 = OpenClass(m, 'C2', 2)
    src.setArrival(k1, Exp(0.5))
    src.setArrival(k2, Exp(0.4))
    q.setService(k1, Exp(2.0))
    q.setService(k2, Exp(3.0))
    m.link(m.serialRouting(src, q, snk))
    s = SolverMAM(m)
    s.runAnalyzer()
    RD = s.getCdfRespT()
    assert [e['class'] for e in RD] == [1, 2]
    RN = np.asarray(s.result.RN)
    trap = np.trapezoid if hasattr(np, 'trapezoid') else np.trapz
    for e in RD:
        t = np.asarray(e['t'])
        p = np.asarray(e['p'])
        assert p[0] == 0.0
        assert np.all(np.diff(p) >= -1e-7)
        assert p[-1] > 0.995
        mean_law = trap(1.0 - p, t)
        assert mean_law == pytest.approx(RN[1, e['class'] - 1], rel=0.02)


def test_mam_cdf_respt_warns_and_returns_empty_off_topology():
    from line_solver.solvers.solver_mam import SolverMAM
    m = Network('tandem')
    src = Source(m, 'Source')
    q1 = Queue(m, 'Queue1', SchedStrategy.FCFS)
    q2 = Queue(m, 'Queue2', SchedStrategy.FCFS)
    snk = Sink(m, 'Sink')
    k = OpenClass(m, 'Class1')
    src.setArrival(k, Exp(1.0))
    q1.setService(k, Exp(3.0))
    q2.setService(k, Exp(4.0))
    m.link(m.serialRouting(src, q1, q2, snk))
    s = SolverMAM(m)
    s.runAnalyzer()
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        RD = s.getCdfRespT()
    assert RD == []
    assert any('not supported by SolverMAM' in str(x.message) for x in w)


def test_ssa_refuses_the_response_time_cdf_by_name():
    from line_solver.solvers.solver_ssa import SolverSSA
    s = SolverSSA(_mm1())
    with pytest.raises(RuntimeError, match='does not record per-job response times'):
        s.getCdfRespT()
    with pytest.raises(NotImplementedError, match='getTranCdfRespT'):
        s.getTranCdfRespT()


def test_ctmc_tran_cdf_respt_refuses_instead_of_relabelling_steady_state():
    from line_solver.solvers.solver_ctmc import SolverCTMC
    s = SolverCTMC(_mm1())
    with pytest.raises(NotImplementedError, match='getTranCdfRespT'):
        s.getTranCdfRespT()


def test_base_networksolver_carries_the_exponential_fallback():
    from line_solver.solvers.base import NetworkSolver
    assert hasattr(NetworkSolver, 'getCdfRespT')
    for name in ('getTranCdfRespT', 'getCdfPassT', 'getTranCdfPassT'):
        assert hasattr(NetworkSolver, name)


@pytest.mark.skipif(shutil.which('qnsolver') is None, reason='qnsolver not installed')
def test_qns_inherits_the_exponential_fallback():
    from line_solver.solvers.wrappers.solver_qns import SolverQNS
    s = SolverQNS(_mm1())
    RD = s.getCdfRespT()
    assert isinstance(RD, list) and len(RD) >= 1
    cell = RD[0]
    # Exponential with the right mean: t = -ln(1-F) * R
    t = np.asarray(cell['t'])
    p = np.asarray(cell['p'])
    r = t[0] / (-np.log(1 - p[0]))
    assert np.allclose(t, -np.log(1 - p) * r, rtol=1e-9)
