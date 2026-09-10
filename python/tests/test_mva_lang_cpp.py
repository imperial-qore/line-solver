"""lang='cpp' regressions for SolverMVA: delegation to the C++ port (line-cli).

The transport is subprocess + JSON, the same one lang='java' uses, so the
assertions here are about the three things that transport has to get right and
nothing else:

* PARITY. Every metric of the station x class table must agree with lang='python'
  to solver tolerance on models both engines accept. These are self-referential
  (native computed in the same test) rather than hardcoded, so the file stays
  correct when a shared method's numbers legitimately change.
* WHAT FALLS BACK. An absent or unrunnable binary degrades to native Python with
  a warning -- that is a platform adaptation, and on an arch with no build there
  is nothing to run. Nothing about the MODEL falls back: a construct the C++
  analyzer refuses must propagate (the refusals themselves are exercised in
  test_solvers_lang_cpp.py and test_nonavg_lang_cpp.py), because silently
  answering with the native engine would report a python number under
  lang='cpp', which is the one outcome this option exists to rule out.
* WHAT CROSSES ON ITS MOMENTS. A distribution family with no JSON form of its
  own (Pareto, Weibull, Lognormal) is written as {mean, scv} and rebuilt from
  those two moments, which is lossless for a moment-based method like MVA. The
  failure mode to rule out is not a refusal but a silent scv=1.
* THE ARITHMETIC KNOB. `arith='exact'` must reach the binary and be reported
  back, so a caller can tell an exact solve from a double one after the fact.

NOTE for anyone extending this file: line-cli is NOT built by `pip install`, so
every test here skips unless the binary is locatable. Point LINE_CLI_BINARY at it,
or build `cpp/build/line-cli`, before assuming a green run means anything.
"""

import os

import numpy as np
import pytest

import line_solver
# The worktree copy of line_solver must be the one under test; a global install
# would silently validate the wrong code.
assert os.path.realpath(__file__).rsplit('/python/', 1)[0] in os.path.realpath(
    line_solver.__file__), (
    "test must import the worktree line_solver, got %s" % line_solver.__file__)

from line_solver import (ClosedClass, Delay, Exp, Network, OpenClass, Pareto, Queue,
                         SchedStrategy, Sink, Source, SolverMVA)

TOL = 1e-8
METRICS = ('QLen', 'Util', 'RespT', 'ResidT', 'ArvR', 'Tput')


def _require_line_mp():
    """Skip unless line-cli can be located: lang='cpp' needs it."""
    try:
        from line_solver.solvers.cpp_dispatch import find_line_cli
        return find_line_cli()
    except Exception as e:
        pytest.skip("line-cli not available: %s" % e)


# --- models both engines accept -------------------------------------------

def _closed_ps():
    m = Network('cqn_ps')
    d = Delay(m, 'Delay')
    q = Queue(m, 'Queue', SchedStrategy.PS)
    c = ClosedClass(m, 'C1', 3, d)
    d.setService(c, Exp(1.0))
    q.setService(c, Exp(2.0))
    P = m.initRoutingMatrix()
    P.set(c, Network.serialRouting(d, q))
    m.link(P)
    return m


def _closed_two_class_fcfs():
    m = Network('cqn_2c')
    d = Delay(m, 'Delay')
    q1 = Queue(m, 'Q1', SchedStrategy.FCFS)
    q2 = Queue(m, 'Q2', SchedStrategy.PS)
    c1 = ClosedClass(m, 'C1', 2, d)
    c2 = ClosedClass(m, 'C2', 1, d)
    for c, (zd, s1, s2) in ((c1, (1.0, 2.0, 3.0)), (c2, (0.5, 1.5, 2.5))):
        d.setService(c, Exp(1.0 / zd))
        q1.setService(c, Exp(s1))
        q2.setService(c, Exp(s2))
    P = m.initRoutingMatrix()
    for c in (c1, c2):
        P.set(c, Network.serialRouting(d, q1, q2))
    m.link(P)
    return m


def _open_mm1():
    m = Network('mm1')
    s = Source(m, 'Source')
    q = Queue(m, 'Queue', SchedStrategy.FCFS)
    k = Sink(m, 'Sink')
    c = OpenClass(m, 'C1')
    s.setArrival(c, Exp(0.5))
    q.setService(c, Exp(1.0))
    P = m.initRoutingMatrix()
    P.set(c, Network.serialRouting(s, q, k))
    m.link(P)
    return m


def _multiserver_closed():
    m = Network('cqn_ms')
    d = Delay(m, 'Delay')
    q = Queue(m, 'Queue', SchedStrategy.PS)
    q.setNumberOfServers(2)
    c = ClosedClass(m, 'C1', 4, d)
    d.setService(c, Exp(1.0))
    q.setService(c, Exp(1.5))
    P = m.initRoutingMatrix()
    P.set(c, Network.serialRouting(d, q))
    m.link(P)
    return m


MODELS = {
    'closed_ps': _closed_ps,
    'closed_two_class_fcfs': _closed_two_class_fcfs,
    'open_mm1': _open_mm1,
    'multiserver_closed': _multiserver_closed,
}


def _table(model, **kw):
    return SolverMVA(model, **kw).getAvgTable()


@pytest.mark.parametrize('name', sorted(MODELS))
def test_cpp_matches_python(name):
    """lang='cpp' reproduces lang='python' on every metric of the avg table."""
    _require_line_mp()
    build = MODELS[name]
    tp = _table(build())
    tc = _table(build(), lang='cpp')
    assert len(tc['QLen']) == len(tp['QLen']), (
        "row count differs: python %d, cpp %d" % (len(tp['QLen']), len(tc['QLen'])))
    assert list(tc['Station']) == list(tp['Station'])
    assert list(tc['JobClass']) == list(tp['JobClass'])
    for col in METRICS:
        np.testing.assert_allclose(
            np.asarray(tc[col], dtype=float), np.asarray(tp[col], dtype=float),
            rtol=TOL, atol=TOL, err_msg="%s: column %s differs" % (name, col))


def test_open_mm1_is_exact_in_both_ports():
    """Both ports take the M/M/1 closed form, so QLen is rho/(1-rho) EXACTLY.

    This began as a divergence test. Native Python excluded the plain M/M/1 from
    the exact qsys route on the ground that "egflin is already exact when
    ca=cs=1"; it is not -- the linearizer stopped 7.6e-7 short of 1. The gate is
    gone and the assertion is now equality with the closed form on BOTH langs,
    which is a stronger statement than agreeing with each other.
    """
    _require_line_mp()
    sc = SolverMVA(_open_mm1(), lang='cpp')
    tc = sc.getAvgTable()
    assert sc.result.method == 'mm1'
    qc = dict(zip(tc['Station'], (float(x) for x in tc['QLen'])))
    tp = _table(_open_mm1())
    qp = dict(zip(tp['Station'], (float(x) for x in tp['QLen'])))
    # rho/(1-rho) with rho = 0.5
    assert abs(qc['Queue'] - 1.0) < 1e-12
    assert abs(qp['Queue'] - 1.0) < 1e-12


def test_cpp_reports_the_engine_that_ran():
    """The container carries the C++ actualmethod and arithmetic, not python's."""
    _require_line_mp()
    s = SolverMVA(_closed_ps(), lang='cpp')
    s.getAvgTable()
    assert s.result.method == 'exact'
    # No --arith given: the port runs in double, which is what is comparable
    # with the other langs.
    assert s.result.arith == 'double'


def test_arith_exact_reaches_the_binary():
    """arith='exact' is forwarded and reported back, and does not change the value.

    An exact solve of a product-form model returns the same doubles: what the
    option buys is that no rounding happened in the middle, which is invisible
    here by construction and would only show on an ill-conditioned model. The
    assertion is therefore that the knob ARRIVES, plus that it did not perturb
    a model where it must not.
    """
    _require_line_mp()
    s = SolverMVA(_closed_ps(), lang='cpp', arith='exact')
    te = s.getAvgTable()
    assert s.result.arith == 'exact'
    tp = _table(_closed_ps())
    for col in METRICS:
        np.testing.assert_allclose(np.asarray(te[col], dtype=float),
                                   np.asarray(tp[col], dtype=float), rtol=TOL, atol=TOL)


def test_absent_binary_falls_back_to_python(monkeypatch):
    """An unrunnable binary degrades to native Python: platform, not model."""
    monkeypatch.setenv('LINE_CLI_BINARY', '/nonexistent/line-cli')
    tc = SolverMVA(_closed_ps(), lang='cpp').getAvgTable()
    tp = _table(_closed_ps())
    for col in METRICS:
        np.testing.assert_allclose(np.asarray(tc[col], dtype=float),
                                   np.asarray(tp[col], dtype=float), rtol=TOL, atol=TOL)


def test_moment_only_family_crosses_on_its_moments_not_on_a_wrong_scv():
    """Pareto has no JSON form of its own and crosses as (mean, scv).

    The writers emit `{type: Pareto, params: {mean, scv}}` for a family the wire
    format cannot represent, and the C++ reader honours those moments (an
    acyclic phase-type matching the two) rather than refusing -- refusing would
    make a model unreadable that the reference itself declares readable. MVA is
    a moment-based method, so this is LOSSLESS for it and the two langs must
    agree exactly.

    The failure this rules out is a silent degradation to scv=1: the second
    assertion pins the answer to Pareto's own scv=1/3, which an exponential of
    the same mean does not reproduce.
    """
    _require_line_mp()

    def build(dist):
        m = Network('par')
        s = Source(m, 'Source')
        q = Queue(m, 'Queue', SchedStrategy.FCFS)
        k = Sink(m, 'Sink')
        c = OpenClass(m, 'C1')
        s.setArrival(c, Exp(0.2))
        q.setService(c, dist)
        P = m.initRoutingMatrix()
        P.set(c, Network.serialRouting(s, q, k))
        m.link(P)
        return m

    par = Pareto(3.0, 2.0)
    assert par.getSCV() == pytest.approx(1.0 / 3.0)
    tp = _table(build(par))
    tc = _table(build(par), lang='cpp')
    for col in METRICS:
        np.testing.assert_allclose(np.asarray(tc[col], dtype=float),
                                   np.asarray(tp[col], dtype=float), rtol=TOL, atol=TOL)

    # M/G/1 at rho = 0.6: QLen is rho + rho^2 (1 + scv) / (2 (1 - rho)), so a
    # dropped scv would land on the exponential answer instead.
    exp_qlen = np.asarray(_table(build(Exp(1.0 / par.getMean())), lang='cpp')['QLen'],
                          dtype=float)
    cpp_qlen = np.asarray(tc['QLen'], dtype=float)
    assert np.max(np.abs(cpp_qlen - exp_qlen)) > 0.1, (
        "the C++ answer does not depend on the Pareto scv, so it was dropped")
