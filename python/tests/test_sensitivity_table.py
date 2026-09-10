"""Regression tests for getSensitivityTable, the solver-level derivative of the
mean performance measures with respect to the service rates.

Two branches produce the table and the tests cover both plus the dispatch
between them:

  'exact'  analytic differentiation of a product-form recursion (pfqn_sens for
           closed models, closed-form BCMP for open ones), available on
           SolverMVA and SolverNC and only for single-server queues in a
           non-mixed model;
  'fd'     forward or central finite differences on the calling solver's own
           predictions, with the service process rate-scaled by dist_scale_rate,
           which every other solver uses.

The numeric goldens are the MATLAB values of
matlab/src/solvers/@NetworkSolver/getSensitivityTable.m on the same model, so
these are cross-codebase parity checks and not self-consistency checks. The two
branches are additionally cross-validated against each other, which is what
catches a wrong chain rule or a wrong perturbation on either side.

dist_scale_rate (the perturbation primitive, port of
matlab/src/lang/processes/dist_scale_rate.m) is tested directly per family: a
pure time scaling divides the mean by the factor and leaves the SCV, hence the
shape, untouched.

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

import os

import numpy as np
import pytest

from line_solver import (Network, Delay, Queue, Source, Sink, ClosedClass,
                         OpenClass, SchedStrategy, Exp, Erlang, HyperExp,
                         Coxian, Cox2, APH, PH, MAP, MMPP2, Det, Uniform,
                         Gamma, Pareto, Weibull, Lognormal, Zipf,
                         SolverMVA, SolverNC, SolverCTMC, SolverFluid,
                         ClassSwitch, LayeredNetwork, Processor, Task, Entry,
                         Activity, SolverLN)
from line_solver.distributions.scaling import dist_scale_rate

# Absolute tolerance of the MATLAB goldens, which are quoted to five decimals.
GOLDEN_TOL = 5e-5

COLS = ['dTput_dRate', 'dRespT_dRate', 'dQLen_dRate', 'dUtil_dRate']


def _closed_model(nservers=1, service=None):
    """Delay(Exp(1)) -> Queue/PS(Exp(2)), one closed class with N=3."""
    model = Network('SensClosed')
    delay = Delay(model, 'Delay')
    queue = Queue(model, 'Queue', SchedStrategy.PS)
    if nservers > 1:
        queue.setNumberOfServers(nservers)
    jobclass = ClosedClass(model, 'Class1', 3, delay)
    delay.setService(jobclass, Exp(1))
    queue.setService(jobclass, Exp(2) if service is None else service)
    model.link(Network.serialRouting(delay, queue))
    return model


def _open_model():
    """Two-class M/M/1-PS open network."""
    model = Network('SensOpen')
    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', SchedStrategy.PS)
    sink = Sink(model, 'Sink')
    class1 = OpenClass(model, 'Class1')
    class2 = OpenClass(model, 'Class2')
    source.setArrival(class1, Exp(0.5))
    source.setArrival(class2, Exp(0.3))
    queue.setService(class1, Exp(2.0))
    queue.setService(class2, Exp(3.0))
    model.link(Network.serialRouting(source, queue, sink))
    return model


def _mixed_model():
    """One open and one closed class sharing a single-server PS queue."""
    model = Network('SensMixed')
    source = Source(model, 'Source')
    delay = Delay(model, 'Delay')
    queue = Queue(model, 'Queue', SchedStrategy.PS)
    sink = Sink(model, 'Sink')
    openclass = OpenClass(model, 'Open')
    closedclass = ClosedClass(model, 'Closed', 2, delay)
    source.setArrival(openclass, Exp(0.3))
    delay.setService(openclass, Exp(1))
    queue.setService(openclass, Exp(2))
    delay.setService(closedclass, Exp(1))
    queue.setService(closedclass, Exp(2))
    P = model.initRoutingMatrix()
    P.set(openclass, openclass, source, queue, 1.0)
    P.set(openclass, openclass, queue, sink, 1.0)
    P.set(closedclass, closedclass, delay, queue, 1.0)
    P.set(closedclass, closedclass, queue, delay, 1.0)
    model.link(P)
    return model


def _row(table, station, jobclass):
    sel = table[(table['Station'] == station) & (table['JobClass'] == jobclass)]
    assert len(sel) == 1, "expected exactly one (%s, %s) row" % (station,
                                                                jobclass)
    return sel.iloc[0]


# ---------------------------------------------------------------------------
# Dispatch
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('solver', [SolverMVA, SolverNC])
def test_exact_solvers_take_the_analytic_branch(solver):
    table = solver(_closed_model()).getSensitivityTable()
    assert table.attrs['method'] == 'exact'
    assert len(table) == 1


@pytest.mark.parametrize('solver', [SolverCTMC, SolverFluid])
def test_other_solvers_take_the_finite_difference_branch(solver):
    table = solver(_closed_model()).getSensitivityTable()
    assert table.attrs['method'] == 'fd'
    assert len(table) == 1


def test_only_queues_are_reported():
    """The delay carries no row: the table is indexed by queueing station."""
    table = SolverMVA(_closed_model()).getSensitivityTable()
    assert list(table['Station']) == ['Queue']


# ---------------------------------------------------------------------------
# Cross-codebase values (MATLAB ground truth)
# ---------------------------------------------------------------------------

def test_exact_matches_matlab():
    row = _row(SolverMVA(_closed_model()).getSensitivityTable(),
               'Queue', 'Class1')
    expected = [0.4903, -0.59, -0.4903, -0.14958]
    np.testing.assert_allclose([row[c] for c in COLS], expected,
                               atol=GOLDEN_TOL, rtol=0)


def test_ctmc_fd_matches_matlab():
    row = _row(SolverCTMC(_closed_model()).getSensitivityTable(),
               'Queue', 'Class1')
    expected = [0.49028, -0.58993, -0.49028, -0.14958]
    np.testing.assert_allclose([row[c] for c in COLS], expected,
                               atol=GOLDEN_TOL, rtol=0)


# ---------------------------------------------------------------------------
# Agreement between the two branches
# ---------------------------------------------------------------------------

def test_central_fd_reproduces_the_exact_branch():
    """Same solver, same model, both branches: a wrong chain rule on the exact
    side or a wrong perturbation on the fd side breaks this."""
    model = _closed_model()
    exact = SolverMVA(model).getSensitivityTable()
    fd = SolverMVA(model).getSensitivityTable(method='fd', scheme='central')
    assert exact.attrs['method'] == 'exact'
    assert fd.attrs['method'] == 'fd'
    np.testing.assert_allclose(fd[COLS].values.astype(float),
                               exact[COLS].values.astype(float),
                               rtol=1e-4, atol=1e-9)


def test_model_is_restored_after_a_finite_difference_sweep():
    """The sweep perturbs the service processes in place and must put them
    back, otherwise a second call returns different numbers."""
    model = _closed_model()
    first = SolverMVA(model).getSensitivityTable(method='fd')
    second = SolverMVA(model).getSensitivityTable(method='fd')
    np.testing.assert_allclose(second[COLS].values.astype(float),
                               first[COLS].values.astype(float),
                               rtol=0, atol=0)


# ---------------------------------------------------------------------------
# Open networks
# ---------------------------------------------------------------------------

def test_open_exact_has_zero_throughput_derivative():
    """Open throughput is lambda*visits, fixed by the arrival process, so it
    does not move with the service rate."""
    table = SolverMVA(_open_model()).getSensitivityTable()
    assert table.attrs['method'] == 'exact'
    assert len(table) == 2
    assert list(table['dTput_dRate']) == [0.0, 0.0]


def test_open_fd_agrees_with_exact():
    model = _open_model()
    exact = SolverMVA(model).getSensitivityTable()
    fd = SolverMVA(model).getSensitivityTable(method='fd', scheme='central')
    assert fd.attrs['method'] == 'fd'
    np.testing.assert_allclose(fd[COLS].values.astype(float),
                               exact[COLS].values.astype(float),
                               rtol=1e-3, atol=1e-9)


# ---------------------------------------------------------------------------
# Scope of the exact branch
# ---------------------------------------------------------------------------

def test_forced_exact_rejected_on_a_non_product_form_solver():
    with pytest.raises(ValueError, match='SolverMVA and SolverNC'):
        SolverCTMC(_closed_model()).getSensitivityTable(method='exact')


def test_forced_exact_rejected_on_a_multiserver_model():
    with pytest.raises(ValueError, match='single-server'):
        SolverMVA(_closed_model(nservers=2)).getSensitivityTable(
            method='exact')


def test_forced_exact_rejected_on_a_mixed_model():
    with pytest.raises(ValueError, match='mixed'):
        SolverMVA(_mixed_model()).getSensitivityTable(method='exact')


def test_auto_falls_back_to_fd_out_of_scope():
    """Out of the analytic scope, 'auto' must degrade to finite differences on
    a solver that supports the exact branch, not raise."""
    table = SolverMVA(_closed_model(nservers=2)).getSensitivityTable()
    assert table.attrs['method'] == 'fd'
    assert np.all(np.isfinite(table[COLS].values.astype(float)))


# ---------------------------------------------------------------------------
# Option validation
# ---------------------------------------------------------------------------

def test_unknown_method_rejected():
    with pytest.raises(ValueError, match='method'):
        SolverMVA(_closed_model()).getSensitivityTable(method='analytic')


def test_unknown_scheme_rejected():
    with pytest.raises(ValueError, match='scheme'):
        SolverMVA(_closed_model()).getSensitivityTable(scheme='backward')


@pytest.mark.parametrize('step', [0.0, -1e-4, 1.0, 2.5, np.inf])
def test_step_outside_the_unit_interval_rejected(step):
    with pytest.raises(ValueError, match='step'):
        SolverMVA(_closed_model()).getSensitivityTable(method='fd', step=step)


# ---------------------------------------------------------------------------
# dist_scale_rate, the perturbation primitive
# ---------------------------------------------------------------------------

_PH_ALPHA = np.array([0.4, 0.6])
_PH_T = np.array([[-3.0, 1.0], [0.0, -2.0]])
_MAP_D0 = np.array([[-2.0, 0.5], [0.2, -1.5]])
_MAP_D1 = np.array([[1.0, 0.5], [0.3, 1.0]])

_FAMILIES = [
    Exp(2.0),
    Erlang.fitMeanAndOrder(0.5, 3),
    HyperExp(0.3, 4.0, 1.0),
    Coxian.fitMeanAndSCV(0.4, 2.0),
    Cox2(3.0, 1.5, 0.4),
    APH.fitMeanAndSCV(0.5, 3.0),
    PH(_PH_ALPHA, _PH_T),
    MAP(_MAP_D0, _MAP_D1),
    MMPP2(1.0, 3.0, 0.5, 0.7),
    Det(0.4),
    Uniform(0.2, 0.8),
    Gamma(2.0, 0.3),
    Pareto(2.5, 0.4),
    Weibull(1.5, 0.7),
    Lognormal.fitMeanAndSCV(0.5, 1.2),
]


@pytest.mark.parametrize('distrib', _FAMILIES,
                         ids=[type(d).__name__ for d in _FAMILIES])
@pytest.mark.parametrize('factor', [0.5, 2.5])
def test_dist_scale_rate_is_a_pure_time_scaling(distrib, factor):
    scaled = dist_scale_rate(distrib, factor)
    assert type(scaled) is type(distrib)
    np.testing.assert_allclose(scaled.getMean(),
                               distrib.getMean() / factor, rtol=1e-9)
    np.testing.assert_allclose(scaled.getSCV(), distrib.getSCV(), rtol=1e-9)


@pytest.mark.parametrize('factor', [0.5, 2.5])
def test_dist_scale_rate_scales_an_nhpp_schedule(factor):
    # A piecewise-constant intensity is time-scaled by lambda(t) ->
    # factor*lambda(factor*t): the rates go up and the breakpoints compress.
    # getSCV is NaN for an NHPP (it is not a renewal process), so the mean,
    # which is the arrival-stationary interval, is what carries the check.
    from line_solver import NHPP
    nhpp = NHPP([0.0, 1.0, 2.0], [2.0, 4.0], True)
    scaled = dist_scale_rate(nhpp, factor)
    assert type(scaled) is NHPP
    np.testing.assert_allclose(scaled.getMean(), nhpp.getMean() / factor, rtol=1e-9)
    np.testing.assert_allclose(scaled.getRates(), np.asarray(nhpp.getRates()) * factor, rtol=1e-9)
    np.testing.assert_allclose(scaled.getBreakpoints(),
                               np.asarray(nhpp.getBreakpoints()) / factor, rtol=1e-9)
    assert scaled.isCyclic() == nhpp.isCyclic()


@pytest.mark.parametrize('factor', [0.5, 2.5])
def test_dist_scale_rate_scales_a_replayer_trace(factor):
    from line_solver import Replayer
    trace = [1.0, 2.0, 3.0, 6.0]
    replayer = Replayer(trace)
    scaled = dist_scale_rate(replayer, factor)
    assert type(scaled) is Replayer
    np.testing.assert_allclose(scaled.trace, np.asarray(trace) / factor, rtol=1e-12)
    np.testing.assert_allclose(scaled.getMean(), replayer.getMean() / factor, rtol=1e-9)
    np.testing.assert_allclose(scaled.getSCV(), replayer.getSCV(), rtol=1e-9)


def test_dist_scale_rate_rejects_a_file_backed_replayer(tmp_path):
    # Rescaling a user-supplied trace file would mean writing a new file behind
    # the caller's back, so it is refused rather than done silently.
    from line_solver import Replayer
    path = tmp_path / 'trace.txt'
    path.write_text('1.0\n2.0\n3.0\n')
    with pytest.raises(ValueError, match='file-backed'):
        dist_scale_rate(Replayer(str(path)), 2.0)


def test_dist_scale_rate_rejects_an_unsupported_family():
    with pytest.raises(ValueError, match='Zipf'):
        dist_scale_rate(Zipf(1.0, 10), 2.0)


@pytest.mark.parametrize('factor', [0.0, -1.0, np.inf])
def test_dist_scale_rate_rejects_a_nonpositive_factor(factor):
    with pytest.raises(ValueError, match='factor'):
        dist_scale_rate(Exp(1.0), factor)


# ---------------------------------------------------------------------------
# Non-exponential service
# ---------------------------------------------------------------------------

def test_ctmc_fd_with_erlang_service():
    """The fd branch is the only one that applies to a non-product-form
    service process; the rate scaling must keep the Erlang an Erlang."""
    model = _closed_model(service=Erlang.fitMeanAndOrder(0.5, 2))
    table = SolverCTMC(model).getSensitivityTable()
    assert table.attrs['method'] == 'fd'
    values = table[COLS].values.astype(float)
    assert np.all(np.isfinite(values))
    row = _row(table, 'Queue', 'Class1')
    # A faster server drains the queue and raises throughput.
    assert row['dTput_dRate'] > 0
    assert row['dQLen_dRate'] < 0
    assert row['dRespT_dRate'] < 0
    assert row['dUtil_dRate'] < 0

# ---------------------------------------------------------------------------
# Class switching and layered networks


def _class_switch_model():
    """Closed model whose two classes switch into each other at a ClassSwitch
    node, so the chain population sits on the reference class alone."""
    model = Network('cs')
    delay = Delay(model, 'D')
    queue = Queue(model, 'Q', SchedStrategy.PS)
    cs = ClassSwitch(model, 'CS')
    c1 = ClosedClass(model, 'C1', 2, delay)
    c2 = ClosedClass(model, 'C2', 0, delay)
    delay.setService(c1, Exp(1.0))
    delay.setService(c2, Exp(1.0))
    queue.setService(c1, Exp(2.0))
    queue.setService(c2, Exp(3.0))
    csmat = np.zeros((2, 2))
    csmat[0, 1] = 1.0
    csmat[1, 0] = 1.0
    cs.setClassSwitchingMatrix(csmat)
    P = model.initRoutingMatrix()
    P[c1, c2] = Network.serialRouting(delay, queue, cs)
    P[c2, c1] = Network.serialRouting(delay, queue, cs)
    model.link(P)
    return model


def test_class_switching_uses_the_analytic_branch():
    # A chain-based model is differentiated at chain level and disaggregated back
    # to the classes, so the analytic branch applies and must agree with the
    # finite differences of the same solver.
    model = _class_switch_model()
    Te = SolverMVA(model).getSensitivityTable()
    assert Te.attrs['method'] == 'exact'
    assert len(Te) == 2
    Tf = SolverMVA(_class_switch_model()).getSensitivityTable(method='fd',
                                                             scheme='central')
    for col in ('dTput_dRate', 'dRespT_dRate', 'dQLen_dRate', 'dUtil_dRate'):
        np.testing.assert_allclose(Te[col].values, Tf[col].values, rtol=1e-4)


def _nonunit_visit_model():
    """Closed model where the queues are visited twice per cycle."""
    model = Network('visits')
    delay = Delay(model, 'D')
    q1 = Queue(model, 'Q1', SchedStrategy.PS)
    q2 = Queue(model, 'Q2', SchedStrategy.PS)
    c = ClosedClass(model, 'C', 3, delay)
    delay.setService(c, Exp(1.0))
    q1.setService(c, Exp(2.0))
    q2.setService(c, Exp(3.0))
    P = model.initRoutingMatrix()
    P[c, c] = np.array([[0, 1, 0], [0, 0, 1], [0.5, 0.5, 0]])
    model.link(P)
    return model


def test_nonunit_visits_match_finite_differences():
    # With visit ratios other than one, a class throughput at a station is X*v
    # and the reported response time is per visit, not the chain residence time.
    # Reading the chain quantities straight out of pfqn_sens got both wrong by a
    # factor v.
    Te = SolverMVA(_nonunit_visit_model()).getSensitivityTable()
    assert Te.attrs['method'] == 'exact'
    Tf = SolverMVA(_nonunit_visit_model()).getSensitivityTable(method='fd',
                                                              scheme='central')
    for col in ('dTput_dRate', 'dRespT_dRate', 'dQLen_dRate', 'dUtil_dRate'):
        np.testing.assert_allclose(Te[col].values, Tf[col].values, rtol=1e-4)


def _layered_model():
    model = LayeredNetwork('LQN-single')
    p1 = Processor(model, 'P1', 1, SchedStrategy.PS)
    p2 = Processor(model, 'P2', 1, SchedStrategy.PS)
    t1 = Task(model, 'T1', 5, SchedStrategy.REF).on(p1).setThinkTime(Exp(1.0 / 2))
    t2 = Task(model, 'T2', 1, SchedStrategy.FCFS).on(p2).setThinkTime(Exp(1.0 / 3))
    e1 = Entry(model, 'E1').on(t1)
    e2 = Entry(model, 'E2').on(t2)
    Activity(model, 'AS1', Exp(10)).on(t1).bound_to(e1).synch_call(e2, 1)
    Activity(model, 'AS2', Exp(20)).on(t2).bound_to(e2).replies_to(e2)
    return model


@pytest.mark.skipif(
    os.environ.get('LINE_SOLVER_LANG') == 'java',
    reason="asserts SolverLN concatenates its native per-layer solvers; the "
           "java-dispatch path solves the whole LayeredNetwork in ONE jline.jar "
           "subprocess, so no native layer solver runs")
def test_layered_solver_delegates_to_its_layers():
    # SolverLN owns no recursion of its own to differentiate: it concatenates
    # what the layer solvers report, one row block per layer.
    solver = SolverLN(_layered_model(), lambda m: SolverMVA(m), verbose=False)
    T = solver.getSensitivityTable()
    assert list(T.columns) == ['Layer', 'Station', 'JobClass', 'dTput_dRate',
                               'dRespT_dRate', 'dQLen_dRate', 'dUtil_dRate']
    assert len(T) == 3
    assert list(T['Layer']) == ['P1', 'P2', 'T2']
    # A layer submodel is chain-based, which the analytic branch now handles by
    # aggregating to chains, so every layer differentiates analytically.
    assert T.attrs['method'] == 'exact'
    assert T.attrs['layer_methods'] == ['exact', 'exact', 'exact']
    # MATLAB values for the same model, to 3 significant digits. Re-recorded
    # 2026-08-10: T2 declares a think time and is NOT a reference task, so the
    # P2 row moved when that think time stopped being charged as a per-request
    # delay; LDES on the same model puts T2 at Tput 2.3, RespT 0.0496 against
    # the layered fixed point's 2.2962 / 0.05.
    np.testing.assert_allclose(T['dTput_dRate'].values,
                               [0.015212, 0.013198, 0.0031711], rtol=1e-3)
    np.testing.assert_allclose(T['dRespT_dRate'].values,
                               [-0.014408, -0.0025, -0.0030035], rtol=1e-3)
    np.testing.assert_allclose(T['dQLen_dRate'].values,
                               [-0.031260, -0.0050842, -0.0067267], rtol=1e-3)
    np.testing.assert_allclose(T['dUtil_dRate'].values,
                               [-0.021455, -0.0050842, -0.0055855], rtol=1e-3)


def test_layered_sweep_leaves_the_fixed_point_intact():
    # The layer sweep restores every perturbed service process, so the layered
    # averages are the same before and after it.
    solver = SolverLN(_layered_model(), lambda m: SolverMVA(m), verbose=False)
    before = solver.getAvgTable()
    solver.getSensitivityTable()
    after = solver.getAvgTable()
    np.testing.assert_allclose(after['RespT'].values.astype(float),
                               before['RespT'].values.astype(float), atol=1e-3)

# ---------------------------------------------------------------------------
# Simulation solvers: common random numbers and the simulator default step

_SSA_SAMPLES = 20000
_SSA_SEED = 23000


def _ssa_solver(seed=_SSA_SEED):
    from line_solver import SolverSSA
    return SolverSSA(_closed_model(), samples=_SSA_SAMPLES, seed=seed)


def _derivatives(table):
    return table[['dTput_dRate', 'dRespT_dRate',
                  'dQLen_dRate', 'dUtil_dRate']].values.astype(float)


def test_simulation_solver_uses_finite_differences():
    assert _ssa_solver().getSensitivityTable().attrs['method'] == 'fd'


def test_simulation_sweep_is_reproducible_under_common_random_numbers():
    # The load-bearing property of the simulation branch: the base and perturbed
    # runs share a seed, so the whole sweep is a deterministic function of it.
    # This is what fails first if the pairing is ever broken.
    first = _ssa_solver().getSensitivityTable()
    second = _ssa_solver().getSensitivityTable()
    np.testing.assert_allclose(_derivatives(first), _derivatives(second),
                               rtol=0, atol=1e-12)


def test_simulation_pins_an_unset_seed():
    from line_solver import SolverSSA
    solver = SolverSSA(_closed_model(), samples=_SSA_SAMPLES)
    solver.getSensitivityTable()
    assert solver.options.seed == _SSA_SEED


def test_simulation_default_step_is_one_percent():
    # Asserted behaviourally: with the seed fixed, the default-step table is the
    # step=1e-2 table and is not the step=1e-3 one.
    default = _derivatives(_ssa_solver().getSensitivityTable())
    explicit = _derivatives(_ssa_solver().getSensitivityTable(step=1e-2))
    other = _derivatives(_ssa_solver().getSensitivityTable(step=1e-3))
    np.testing.assert_allclose(default, explicit, rtol=0, atol=1e-12)
    assert not np.allclose(default, other, rtol=0, atol=1e-12)


def test_simulation_derivatives_have_the_right_sign_and_magnitude():
    # A faster server raises throughput and shortens the queue. The magnitude
    # band is deliberately wide: at 20k samples the Monte Carlo error is tens of
    # percent, and a tighter band would only make the test flaky. It still
    # catches a missing visit factor or an unpaired seed.
    sim = _derivatives(_ssa_solver().getSensitivityTable())[0]
    exact = _derivatives(SolverMVA(_closed_model()).getSensitivityTable())[0]
    assert sim[0] > 0
    assert sim[1] < 0 and sim[2] < 0 and sim[3] < 0
    for got, want in zip(sim, exact):
        assert 0.5 * abs(want) <= abs(got) <= 2.0 * abs(want)
