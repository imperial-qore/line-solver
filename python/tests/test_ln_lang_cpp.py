"""lang='cpp' regressions for SolverLN: delegation of a layered solve to line-cli.

The default layered transport is NOT the JSON one the other solvers use:
line-cli reads a LayeredNetwork from `.lqnx` (`-i lqnx`), so the model goes out
through `LayeredNetwork.writeXML`. That import format carries less than the JSON
one, and the three facts this file exists to pin are consequences of it:

* PARITY on a model the interchange carries losslessly. Every metric of the
  layered AvgTable must agree with lang='python' to the outer fixed point's own
  tolerance, and the NaN mask (which metrics are undefined for a Processor, a
  Task, an Entry) must agree exactly. The reference values are computed natively
  in the same test rather than hardcoded, so the file stays correct when a shared
  method's numbers legitimately change.
* A LOSSY MODEL SWITCHES TRANSPORT rather than losing the value. The lqnx schema
  accepts think-time on a reference task only, and on gallery_lqn_basic dropping
  it moves a processor utilization from 0.041 to 0.44. Such a model is therefore
  sent through `save_model` (model.json), whose reader takes think time on any
  task, so it still solves and still agrees with lang='python'.
* REFUSAL on an option no transport carries. A non-MVA/Fluid layer solver has no
  C++ counterpart, and the layered CLI has no --method or relaxation flags. Every
  one of those raises, because answering anyway would report numbers produced by
  another engine.

The single automatic fallback is an absent or unrunnable binary, which is a
platform adaptation and warns.

NOTE: line-cli is not built by `pip install`, so every test skips unless the
binary is locatable. Point LINE_CLI_BINARY at it, or build `cpp/build/line-cli`.
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

from line_solver import (Activity, Entry, Exp, LayeredNetwork, Processor, SchedStrategy,
                         SolverLN, SolverNC, Task)

# The comparison is between two OUTER FIXED POINTS, not two evaluations of one
# formula, so the achievable agreement is set by the LN loop and not by machine
# precision. Both engines default to iter_tol=5e-3 and both are converged in the
# sense that tightening it to 1e-9 moves them by ~1e-5; what remains is a spread
# in how the layer sequence is relaxed, and it is THREE-way rather than specific
# to this transport: on the three_tier model below T1's queue length is 2.513378
# under lang='python', 2.513754 under lang='cpp' and 2.514187 under lang='java',
# i.e. 3e-4 relative across ports. The tolerance is therefore the outer loop's
# own, which still fails by orders of magnitude on a layer that was BUILT
# differently (a missing class, a wrong population, a dropped call), which is the
# failure mode this transport can actually introduce.
TOL = 5e-3
METRICS = ('QLen', 'Util', 'RespT', 'ResidT', 'Tput')


def _require_line_cli():
    """Skip unless line-cli can be located: lang='cpp' needs it."""
    try:
        from line_solver.solvers.cpp_dispatch import find_line_cli
        return find_line_cli()
    except Exception as e:
        pytest.skip("line-cli not available: %s" % e)


# --- models the .lqnx interchange carries losslessly ------------------------

def _two_tier():
    """Client task calling one server task, think time on the REF task only."""
    m = LayeredNetwork('lqn_two_tier')
    p1 = Processor(m, 'P1', 1, SchedStrategy.PS)
    p2 = Processor(m, 'P2', 1, SchedStrategy.PS)
    t1 = Task(m, 'T1', 5, SchedStrategy.REF).on(p1).set_think_time(Exp(1 / 2))
    t2 = Task(m, 'T2', 3, SchedStrategy.FCFS).on(p2)
    e1 = Entry(m, 'E1').on(t1)
    e2 = Entry(m, 'E2').on(t2)
    Activity(m, 'A1', Exp(5)).on(t1).bound_to(e1).synch_call(e2, 2)
    Activity(m, 'A2', Exp(2)).on(t2).bound_to(e2).replies_to(e2)
    return m


def _three_tier():
    """Client -> middle -> database, two processors, multi-server middle tier."""
    m = LayeredNetwork('lqn_three_tier')
    p1 = Processor(m, 'P1', 2, SchedStrategy.PS)
    p2 = Processor(m, 'P2', 1, SchedStrategy.PS)
    t1 = Task(m, 'T1', 10, SchedStrategy.REF).on(p1).set_think_time(Exp(1 / 3))
    t2 = Task(m, 'T2', 4, SchedStrategy.FCFS).on(p1)
    t3 = Task(m, 'T3', 2, SchedStrategy.FCFS).on(p2)
    e1 = Entry(m, 'E1').on(t1)
    e2 = Entry(m, 'E2').on(t2)
    e3 = Entry(m, 'E3').on(t3)
    Activity(m, 'A1', Exp(10)).on(t1).bound_to(e1).synch_call(e2, 1)
    Activity(m, 'A2', Exp(20)).on(t2).bound_to(e2).synch_call(e3, 3).replies_to(e2)
    Activity(m, 'A3', Exp(8)).on(t3).bound_to(e3).replies_to(e3)
    return m


def _three_tier_forwarding():
    """_three_tier with E2 FORWARDING to E3 instead of calling it.

    THE POINT IS THE LAYERING, not the forwarding: 'srvn.ph' -- what the 'srvn'
    default resolves to on every model that admits it -- composes one phase-type
    class per caller and therefore disables the interlocking correction outright
    (both `build_layers_ph` in the C++ and `_build_layers_ph` natively). A
    forwarding call is one of the constructs the composed law cannot represent,
    so this model takes the 'srvn.cs' layering, where interlocking IS the
    difference between two fixed points and a flag that fails to reach the
    binary is therefore visible.
    """
    m = LayeredNetwork('lqn_three_tier_fwd')
    p1 = Processor(m, 'P1', 2, SchedStrategy.PS)
    p2 = Processor(m, 'P2', 1, SchedStrategy.PS)
    t1 = Task(m, 'T1', 10, SchedStrategy.REF).on(p1).set_think_time(Exp(1 / 3))
    t2 = Task(m, 'T2', 4, SchedStrategy.FCFS).on(p1)
    t3 = Task(m, 'T3', 2, SchedStrategy.FCFS).on(p2)
    e1 = Entry(m, 'E1').on(t1)
    e2 = Entry(m, 'E2').on(t2)
    e3 = Entry(m, 'E3').on(t3)
    Activity(m, 'A1', Exp(10)).on(t1).bound_to(e1).synch_call(e2, 1)
    Activity(m, 'A2', Exp(20)).on(t2).bound_to(e2)
    e2.addForwarding(e3, 1.0)
    Activity(m, 'A3', Exp(8)).on(t3).bound_to(e3).replies_to(e3)
    return m


def _multi_entry():
    """One server task with two entries of different demand, called unevenly."""
    m = LayeredNetwork('lqn_multi_entry')
    p1 = Processor(m, 'P1', 1, SchedStrategy.PS)
    p2 = Processor(m, 'P2', 1, SchedStrategy.PS)
    t1 = Task(m, 'T1', 4, SchedStrategy.REF).on(p1).set_think_time(Exp(1 / 1))
    t2 = Task(m, 'T2', 2, SchedStrategy.FCFS).on(p2)
    e1 = Entry(m, 'E1').on(t1)
    e2 = Entry(m, 'E2').on(t2)
    e3 = Entry(m, 'E3').on(t2)
    Activity(m, 'A1', Exp(6)).on(t1).bound_to(e1).synch_call(e2, 1).synch_call(e3, 2)
    Activity(m, 'A2', Exp(4)).on(t2).bound_to(e2).replies_to(e2)
    Activity(m, 'A3', Exp(9)).on(t2).bound_to(e3).replies_to(e3)
    return m


MODELS = {
    'two_tier': _two_tier,
    'three_tier': _three_tier,
    'multi_entry': _multi_entry,
}


def _table(model, **kw):
    # The reference arm of every comparison here is native python, so it is named
    # rather than inherited: LINE_SOLVER_LANG=java would otherwise silently make
    # the "python" side a JAR solve and the file would compare cpp against java.
    kw.setdefault('lang', 'python')
    s = SolverLN(model, **kw)
    s._table_silent = True
    return s.getAvgTable()


@pytest.mark.parametrize('name', sorted(MODELS))
def test_cpp_matches_python(name):
    """lang='cpp' reproduces lang='python' on the whole layered AvgTable."""
    _require_line_cli()
    build = MODELS[name]
    tp = _table(build())
    tc = _table(build(), lang='cpp')
    assert list(tc['Node']) == list(tp['Node']), "element order/naming differs"
    assert list(tc['NodeType']) == list(tp['NodeType'])
    for col in METRICS:
        a = np.asarray(tp[col], dtype=float)
        b = np.asarray(tc[col], dtype=float)
        # An undefined metric must be undefined in both: a NaN read as 0 (or the
        # reverse) is exactly what a wire-format bug looks like in this table.
        assert np.array_equal(np.isnan(a), np.isnan(b)), (
            "%s: column %s differs in which elements it defines" % (name, col))
        np.testing.assert_allclose(b[~np.isnan(b)], a[~np.isnan(a)],
                                   rtol=TOL, atol=TOL,
                                   err_msg="%s: column %s differs" % (name, col))


def test_ensemble_avg_matches_python():
    """getEnsembleAvg (the vector form) agrees index by index, ArvR included."""
    _require_line_cli()
    QNp, UNp, RNp, TNp, ANp, WNp = SolverLN(_three_tier()).getEnsembleAvg()
    QNc, UNc, RNc, TNc, ANc, WNc = SolverLN(_three_tier(), lang='cpp').getEnsembleAvg()
    for nm, a, b in (('QN', QNp, QNc), ('UN', UNp, UNc), ('RN', RNp, RNc),
                     ('TN', TNp, TNc), ('WN', WNp, WNc)):
        assert a.shape == b.shape, "%s: shape differs" % nm
        assert np.array_equal(np.isnan(a), np.isnan(b)), "%s: NaN mask differs" % nm
        np.testing.assert_allclose(b[~np.isnan(b)], a[~np.isnan(a)], rtol=TOL, atol=TOL,
                                   err_msg="%s differs" % nm)
    # AN is unfilled in both: the native layered solve never computes it and the
    # C++ layered result has no arrival-rate vector, so this is agreement by
    # construction rather than a metric quietly lost in transport.
    assert np.all(np.isnan(ANp)) and np.all(np.isnan(ANc))


def test_cpp_reports_convergence_of_the_solve_that_ran():
    """hasconverged comes from line-cli, not from a fixed point that never ran."""
    _require_line_cli()
    s = SolverLN(_two_tier(), lang='cpp')
    assert not s.hasconverged
    s._table_silent = True
    s.getAvgTable()
    assert s.hasconverged, "a converged C++ layered solve must set hasconverged"
    assert getattr(s, 'it', 0) > 0, "the iteration count must be carried back"


def test_arith_extended_precision_reaches_the_binary():
    """A wider arithmetic backend is forwarded and does not perturb the answer.

    `real:64` is used and NOT `exact`: rational arithmetic on an outer fixed point
    grows the coefficients without bound, so an exact layered solve does not
    terminate in practice (measured: still spinning after 7 minutes on the
    two-tier model, against 0.04s for real:64, which agrees with double to every
    printed digit). The fixed-precision backends carry no such blowup.
    """
    _require_line_cli()
    te = _table(_two_tier(), lang='cpp', arith='real:64')
    tp = _table(_two_tier())
    for col in METRICS:
        a = np.asarray(tp[col], dtype=float)
        b = np.asarray(te[col], dtype=float)
        np.testing.assert_allclose(b[~np.isnan(b)], a[~np.isnan(a)], rtol=TOL, atol=TOL,
                                   err_msg="arith='real:64' moved column %s" % col)


def test_arith_is_forwarded_verbatim(monkeypatch):
    """Whatever arith the caller asked for lands on the command line, unmapped.

    Checked at the argv rather than by running it, so the assertion covers
    `exact` too -- which is forwarded like any other value (it is the caller's
    call), and which the test above documents as non-terminating on a layered
    fixed point.
    """
    _require_line_cli()
    from line_solver.solvers import cpp_dispatch

    seen = {}

    def fake_run(binary, cmd, timeout=None):
        seen['cmd'] = cmd
        raise RuntimeError('stop here')

    monkeypatch.setattr(cpp_dispatch, '_run_line_cli', fake_run)
    with pytest.raises(RuntimeError):
        _table(_two_tier(), lang='cpp', arith='exact')
    cmd = seen['cmd']
    assert cmd[cmd.index('--arith') + 1] == 'exact'
    assert cmd[cmd.index('-i') + 1] == 'lqnx', "the layered path must not send -i json"
    assert cmd[cmd.index('--layer-solver') + 1] == 'mva'


def test_absent_binary_falls_back_to_python(monkeypatch):
    """An unrunnable binary degrades to native Python: platform, not model."""
    monkeypatch.setenv('LINE_CLI_BINARY', '/nonexistent/line-cli')
    tc = _table(_two_tier(), lang='cpp')
    tp = _table(_two_tier())
    for col in METRICS:
        a = np.asarray(tp[col], dtype=float)
        b = np.asarray(tc[col], dtype=float)
        np.testing.assert_allclose(b[~np.isnan(b)], a[~np.isnan(a)], rtol=TOL, atol=TOL)


def test_nonref_think_time_switches_transport_instead_of_being_dropped():
    """A think time the .lqnx schema cannot carry routes through model.json.

    LINE gives a non-reference task's think time to its callers; the lqnx schema
    accepts think-time on reference tasks only, and `writeXML` reports the loss
    and writes the file anyway. `_lqnx_lossy_tasks` names such a model, and
    `solve_lqn_via_cpp` then serializes it with `save_model` instead, whose
    reader takes think time on any task. What is asserted is the consequence:
    the model still solves under lang='cpp', and to the SAME numbers as native.
    A silent drop would not fail here by a tolerance -- it caps T2's throughput
    at multiplicity/think and moves the layer by orders of magnitude.
    """
    _require_line_cli()

    def build():
        m = _two_tier()
        for proc in m.processors:
            for task in proc.tasks:
                if task.name == 'T2':
                    task.set_think_time(Exp(1 / 4))
        return m

    from line_solver.solvers.cpp_dispatch import _lqnx_lossy_tasks
    s = SolverLN(build(), lang='cpp')
    assert 'T2' in _lqnx_lossy_tasks(s), "the lossy-task probe must name T2"

    tp = _table(build())
    tc = _table(build(), lang='cpp')
    assert len(tp['QLen']) > 0
    for col in METRICS:
        a = np.asarray(tp[col], dtype=float)
        b = np.asarray(tc[col], dtype=float)
        np.testing.assert_allclose(b[~np.isnan(b)], a[~np.isnan(a)], rtol=TOL, atol=TOL)
        np.testing.assert_array_equal(np.isnan(a), np.isnan(b))


def test_unported_layer_solver_is_refused():
    """A layer solver with no C++ counterpart raises instead of becoming MVA.

    The layer solver is what the outer fixed point is a fixed point OF, so
    substituting MVA for it answers a different question under the caller's lang.

    CTMC is the refused one, NOT NC: `--layer-solver` takes mva, nc, fluid and
    ssa, so an NC factory is CARRIED rather than refused (a stale refusal of it
    is what left lqn_twotasks and lqn_ofbiz with no table at all). A CTMC layer
    exists only as an automatic per-layer SUBSTITUTION the C++ makes for itself,
    never as an ensemble-wide choice, so naming it as the factory has no
    counterpart to run.
    """
    _require_line_cli()
    from line_solver import SolverCTMC
    s = SolverLN(_two_tier(), SolverCTMC, lang='cpp')
    s._table_silent = True
    with pytest.raises(RuntimeError) as excinfo:
        s.getAvgTable()
    msg = str(excinfo.value)
    assert 'CTMC' in msg, "the refusal must name the layer solver, got: %s" % msg


def test_nc_layer_solver_is_carried():
    """An NC layer factory reaches the binary and agrees with lang='python'."""
    _require_line_cli()

    def table(lang):
        s = SolverLN(_two_tier(), SolverNC, lang=lang)
        s._table_silent = True
        return s.getAvgTable()

    tp, tc = table('python'), table('cpp')
    assert list(tc['Node']) == list(tp['Node'])
    for col in METRICS:
        a = np.asarray(tp[col], dtype=float)
        b = np.asarray(tc[col], dtype=float)
        np.testing.assert_allclose(b[~np.isnan(b)], a[~np.isnan(a)], rtol=TOL, atol=TOL)


def test_unforwardable_options_are_refused():
    """Options the layered CLI cannot carry raise rather than being dropped.

    --method does not exist on the layered path, and the relaxation and layering
    settings have no flags at all; each of them changes the iterate, so running
    without them would silently solve under different settings than requested.
    """
    _require_line_cli()
    with pytest.raises(RuntimeError) as excinfo:
        _table(_two_tier(), lang='cpp', method='mwba.upper')
    assert 'method' in str(excinfo.value)

    with pytest.raises(RuntimeError) as excinfo:
        _table(_two_tier(), lang='cpp', config={'layering': 'flat'})
    assert 'layering' in str(excinfo.value)

    with pytest.raises(RuntimeError) as excinfo:
        _table(_two_tier(), lang='cpp', config={'relax': 'adaptive'})
    assert 'relax' in str(excinfo.value)


def test_analyses_the_layered_cli_lacks_are_refused():
    """getTranAvg and getSensitivityTable raise under lang='cpp'.

    The C++ layered path implements -a avg and nothing else, and both of these are
    native computations, so serving them would put a C++ label on python numbers.
    Neither needs the binary to be present: the refusal is about the analysis.
    """
    s = SolverLN(_two_tier(), lang='cpp')
    with pytest.raises(RuntimeError) as excinfo:
        s.getTranAvg()
    assert 'getTranAvg' in str(excinfo.value)

    with pytest.raises(RuntimeError) as excinfo:
        s.getSensitivityTable()
    assert 'getSensitivityTable' in str(excinfo.value)


def test_forwardable_options_reach_the_binary():
    """iter_max/iter_tol and interlocking DO have flags and must be honoured.

    A one-iteration cap is the cheapest observable proof the flag arrived: the
    C++ solve must then report non-convergence, which the default-tolerance solve
    of the same model does not.
    """
    _require_line_cli()
    s = SolverLN(_three_tier(), lang='cpp', iter_max=1)
    s._table_silent = True
    s.getAvgTable()
    assert not s.hasconverged, "iter_max=1 must not report a converged solve"

    # interlocking off is a different fixed point, not a different tolerance, so
    # it must change at least one metric of a model that has interlocked callers.
    # The model is the FORWARDING one: under the 'srvn.ph' layering the default
    # resolves to, the correction is not applied by either port, so no flag can
    # be observed there -- see `_three_tier_forwarding`.
    tp = _table(_three_tier_forwarding(), lang='cpp')
    ti = _table(_three_tier_forwarding(), lang='cpp', config={'interlocking': False})
    a = np.asarray(tp['Util'], dtype=float)
    b = np.asarray(ti['Util'], dtype=float)
    assert not np.allclose(a[~np.isnan(a)], b[~np.isnan(b)], rtol=1e-9, atol=1e-9), (
        "--no-interlocking did not reach the binary: utilizations are identical")
