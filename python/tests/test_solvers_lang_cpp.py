"""lang='cpp' regressions for the Network solvers: NC, CTMC, MAM, FLD, SSA, BA, AUTO.

SolverMVA has its own file (test_mva_lang_cpp.py) and SolverLN another
(test_ln_lang_cpp.py, whose transport is .lqnx rather than JSON). This one covers
the remaining solvers `line-cli` serves on the model-solving path, and the three
properties that are specific to widening the bridge past one solver:

* PARITY, PER SOLVER. Each engine must reproduce its own lang='python' answer,
  at the tolerance that engine can actually be held to (exact for NC/CTMC/BA,
  approximate-but-deterministic for MAM, method-dependent for FLD, sampling
  error for SSA). Reference values are computed natively in the same test rather
  than hardcoded.
* THE KNOBS ARE GATED BY SOLVER, not just by whether the caller set them.
  line-cli REFUSES an option the chosen engine does not have (`--tol` on `-s ba`,
  `--samples` outside `-s ssa`, `--cutoff` outside `-s ctmc`), so a solver whose
  native default merely differs from the generic one must not have that default
  forwarded. A regression here shows up as an exit-2 refusal on a solve that
  passed no options at all.
* WHAT IS NOT REACHABLE. The bridge serves `-a avg`: line-cli prints `-o json`
  for that analysis only, so every other getter (probabilities, response-time
  CDFs, transients, rewards, cache tables, sample paths) must RAISE rather than
  answer natively under a C++ label.

NOTE: line-cli is not built by `pip install`, so the parity tests skip unless the
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

from line_solver import (ClosedClass, Delay, Exp, Network, OpenClass, Queue, SchedStrategy,
                         Sink, Source, SolverAUTO, SolverBA, SolverCTMC, SolverFLD, SolverJMT,
                         SolverAG, SolverMAM, SolverMVA, SolverNC, SolverSSA)

TOL = 1e-8
METRICS = ('QLen', 'Util', 'RespT', 'ResidT', 'ArvR', 'Tput')


def _require_line_cli():
    """Skip unless line-cli can be located: lang='cpp' needs it."""
    try:
        from line_solver.solvers.cpp_dispatch import find_line_cli
        return find_line_cli()
    except Exception as e:
        pytest.skip("line-cli not available: %s" % e)


_JSON_ARMS = {}


def _require_json_capable_arms():
    """
    Skip unless the located binary emits JSON from the SSA and fluid arms.

    Those two arms used to accept `-o json` and print the readable table anyway,
    so a binary built before that fix cannot serve `-s ssa` / `-s fluid` through
    this bridge at all. Probed once and cached: an old binary in the tree is an
    environment fact, and reporting it as a failure of the wiring under test would
    point at the wrong file.
    """
    binary = _require_line_cli()
    if binary not in _JSON_ARMS:
        import subprocess
        import tempfile
        from line_solver.io.linemodel_io import save_model
        with tempfile.NamedTemporaryFile(suffix='.json', delete=False) as fh:
            path = fh.name
        save_model(_closed_ps(), path)
        proc = subprocess.run([binary, '-f', path, '-i', 'json', '-s', 'ssa', '-a', 'avg',
                               '-o', 'json', '--samples', '500'],
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        os.unlink(path)
        _JSON_ARMS[binary] = b'{' in proc.stdout
    if not _JSON_ARMS[binary]:
        pytest.skip("'%s' predates the shared JSON emitter for -s ssa / -s fluid; rebuild it"
                    % binary)
    return binary


# --- models -----------------------------------------------------------------

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


def _closed_two_class():
    m = Network('cqn_2c')
    d = Delay(m, 'Delay')
    q = Queue(m, 'Queue', SchedStrategy.PS)
    c1 = ClosedClass(m, 'C1', 2, d)
    c2 = ClosedClass(m, 'C2', 1, d)
    for c, (z, s) in ((c1, (1.0, 2.0)), (c2, (0.5, 1.5))):
        d.setService(c, Exp(1.0 / z))
        q.setService(c, Exp(s))
    P = m.initRoutingMatrix()
    for c in (c1, c2):
        P.set(c, Network.serialRouting(d, q))
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


def _assert_tables_agree(tp, tc, rtol, atol, label, cols=METRICS):
    assert list(tc['Station']) == list(tp['Station']), "%s: station order differs" % label
    assert list(tc['JobClass']) == list(tp['JobClass']), "%s: class order differs" % label
    for col in cols:
        np.testing.assert_allclose(
            np.asarray(tc[col], dtype=float), np.asarray(tp[col], dtype=float),
            rtol=rtol, atol=atol, err_msg="%s: column %s differs" % (label, col))


# --- parity, per solver -----------------------------------------------------

@pytest.mark.parametrize('ctor,build,name', [
    (SolverNC, _closed_ps, 'NC/closed_ps'),
    (SolverNC, _closed_two_class, 'NC/closed_2c'),
    (SolverCTMC, _closed_ps, 'CTMC/closed_ps'),
    (SolverCTMC, _closed_two_class, 'CTMC/closed_2c'),
    (SolverBA, _closed_ps, 'BA/closed_ps'),
])
def test_exact_solvers_match_python(ctor, build, name):
    """NC, CTMC and BA are closed forms or direct solves: they agree to 1e-8.

    SolverBA has NO lang='java' counterpart -- the JAR CLI's validateSolver has no
    `ba` token -- so lang='cpp' is the only delegation it has, and this is the
    only place it is exercised.
    """
    _require_line_cli()
    tp = ctor(build()).getAvgTable()
    tc = ctor(build(), lang='cpp').getAvgTable()
    _assert_tables_agree(tp, tc, TOL, TOL, name)


def test_mam_matches_python():
    """MAM's decomposition is approximate but deterministic, so it still agrees."""
    _require_line_cli()
    tp = SolverMAM(_open_mm1()).getAvgTable()
    tc = SolverMAM(_open_mm1(), lang='cpp').getAvgTable()
    _assert_tables_agree(tp, tc, 1e-6, 1e-9, 'MAM/mm1')


def test_ag_matches_python():
    """SolverAG's RCAT fixed point is approximate but deterministic.

    `line-cli -s ag` had served this arm all along while MATLAB's CPPLINE
    refused the token, so the C++ AG engine reached lang='cpp' only from here.
    """
    _require_line_cli()
    tp = SolverAG(_open_mm1()).getAvgTable()
    tc = SolverAG(_open_mm1(), lang='cpp').getAvgTable()
    _assert_tables_agree(tp, tc, 1e-6, 1e-9, 'AG/mm1')


def test_ag_knobs_are_gated_by_the_solver(monkeypatch):
    """AG carries tol, iter_max and maxStates -- and no iter_tol, and no backend.

    `AgOptions` has no `iter_tol`, and line-cli refuses the flag on `-s ag` the
    way it does on `-s mam`, so the bridge must not send it. `maxStates` is a
    TRUNCATION LEVEL and therefore part of the answer, so it must; and the
    execution backends have no flag at all, so a non-serial one is refused
    rather than silently serialized.
    """
    _require_line_cli()
    from line_solver.solvers import cpp_dispatch

    seen = {}

    def fake_run(binary, cmd, timeout=None):
        seen['cmd'] = cmd
        raise RuntimeError('stop here')

    monkeypatch.setattr(cpp_dispatch, '_run_line_cli', fake_run)

    with pytest.raises(RuntimeError):
        SolverAG(_open_mm1(), lang='cpp').getAvgTable()
    cmd = seen['cmd']
    assert cmd[cmd.index('-s') + 1] == 'ag'
    assert '--iter_tol' not in cmd, "AgOptions has no iter_tol; line-cli refuses it on -s ag"
    assert '--max-states' not in cmd, "an untouched maxStates must not be sent as if set"

    with pytest.raises(RuntimeError):
        SolverAG(_open_mm1(), lang='cpp', config={'maxStates': 250}).getAvgTable()
    cmd = seen['cmd']
    assert cmd[cmd.index('--max-states') + 1] == '250'

    with pytest.raises(RuntimeError, match='exec'):
        SolverAG(_open_mm1(), lang='cpp', config={'exec': 'parallel'}).getAvgTable()


def test_fluid_matches_python_on_the_same_method():
    """FLD agrees when both ports run the SAME method; the default does not match.

    The reference resolves method='default' to the second-order closure
    'minnormal' wherever it applies, and the C++ has not ported that method: its
    own default resolves to 'matrix'. Under an EXPLICIT shared method the two
    integrators agree to integrator tolerance, which is what this asserts. The
    default case is asserted to REFUSE in the test below, rather than to return
    the matrix answer under a python-default label.
    """
    _require_json_capable_arms()
    tp = SolverFLD(_closed_ps(), method='matrix').getAvgTable()
    tc = SolverFLD(_closed_ps(), method='matrix', lang='cpp').getAvgTable()
    # Different ODE integrators (LSODA both sides, but different drivers), so the
    # bar is the integrator's, not machine precision.
    _assert_tables_agree(tp, tc, 1e-6, 1e-8, 'FLD/matrix')


def test_fluid_default_resolves_to_the_same_method_on_both_sides(monkeypatch):
    """
    `method='default'` must select the SAME fluid method in both ports.

    WHAT IS ASSERTED IS THE AGREEMENT, not who computes it. The two ports once
    resolved 'default' differently -- the reference, the native solver and the JAR
    prefer the second-order closure 'minnormal' wherever it applies while the C++
    resolved to 'matrix' -- and the same closed PS network came back QLen [2, 1]
    under lang='cpp' against [1.575, 1.425] under lang='python', both correct for
    the method that ran and separable only by reading result.method. The C++ now
    carries the same resolution, so nothing is forwarded and the check is that the
    numbers land together; a divergence here means the two resolutions have parted
    again, whichever side moved.
    """
    _require_json_capable_arms()
    assert SolverFLD(_closed_ps()).resolveMethod(
        SolverFLD(_closed_ps()).options) == 'minnormal', (
        "this test assumes the native default resolves to minnormal on this model")

    from line_solver.solvers import cpp_dispatch
    seen = []
    real = cpp_dispatch._run_line_cli
    monkeypatch.setattr(cpp_dispatch, '_run_line_cli',
                        lambda b, c, timeout=None: (seen.append(list(c)),
                                                    real(b, c, timeout=timeout))[1])
    tc = SolverFLD(_closed_ps(), lang='cpp').getAvgTable()
    tp = SolverFLD(_closed_ps()).getAvgTable()
    # 2e-4, and the number is a MEASUREMENT of an open divergence, not a
    # comfortable margin. MATLAB's SolverFLD answers QLen 1.5751190816 on this
    # model. THE GAP IS NOT THE C++ PORT'S, as this comment said while the JAR
    # still declined minnormal here and could not be compared: measured
    # 2026-08-13, native Python returns 1.575219863772 and the C++ 1.575219863697,
    # i.e. the SAME value to ten digits and both 6.4e-5 relative HIGH, while the
    # JAR reproduces the reference to 2.7e-11. So the two ports agree with each
    # other and the pair differs from the reference; the divergence is in the
    # closure arithmetic they share, not in the delegation this test covers. It
    # does not move when iter_tol goes from 5e-3 to 1e-10, so it is not the outer
    # fixed point's residual either. TIGHTEN THIS to 1e-6 once that is chased
    # down, and treat a failure at 2e-4 as the two ports having parted rather
    # than as a flaky bound.
    _assert_tables_agree(tp, tc, 2e-4, 1e-8, 'FLD/default')
    # No --method on the wire: the C++ resolves its own default, and sending a
    # name computed from this side's struct would override that resolution.
    assert not any('--method' in c for c in seen), seen


def test_ssa_is_a_sample_path_and_is_compared_as_one():
    """SSA delegates and lands within sampling error of the exact answer.

    The two engines do NOT share an RNG stream, so equality with lang='python' is
    not a property this can assert at any seed (see the SSA goldens: seeding is
    lang-specific). What IS assertable is that the C++ path ran and produced a
    sample path of the SAME model: its queue lengths must sit near the exact MVA
    values, and the population must be conserved.
    """
    _require_json_capable_arms()
    exact = SolverMVA(_closed_ps()).getAvgTable()
    qe = np.asarray(exact['QLen'], dtype=float)
    tc = SolverSSA(_closed_ps(), samples=200000, seed=23000, lang='cpp').getAvgTable()
    qc = np.asarray(tc['QLen'], dtype=float)
    assert abs(qc.sum() - 3.0) < 1e-6, "closed population not conserved: %s" % qc
    np.testing.assert_allclose(qc, qe, rtol=0.05, atol=0.05,
                               err_msg="SSA is further from exact than sampling error explains")


def test_auto_propagates_lang_and_arith_to_the_chosen_solver():
    """AUTO under lang='cpp' answers with the C++ engine it selected.

    The delegate is what runs, NOT line-cli's own `-s auto`: that would re-run the
    chooser on the C++ side and could pick a different engine than the one AUTO
    reported, leaving the printed selection unattributable.
    """
    _require_line_cli()
    s = SolverAUTO(_closed_ps(), lang='cpp', arith='real:64')
    tc = s.getAvgTable()
    tp = SolverAUTO(_closed_ps()).getAvgTable()
    _assert_tables_agree(tp, tc, 1e-6, 1e-8, 'AUTO/closed_ps')
    inner = getattr(s, 'solver', None) or getattr(s, '_solver', None)
    if inner is not None:
        assert getattr(inner.options, 'lang', None) == 'cpp'
        assert getattr(inner.options, 'arith', None) == 'real:64'


# --- knobs, and the gating that keeps them attributable ---------------------

def test_untouched_options_are_not_forwarded_to_solvers_that_refuse_them():
    """A default the engine has no flag for must not be sent as if the user set it.

    line-cli refuses `--tol` on `-s ba` ("a bound is a closed form, with nothing
    to converge"), `--samples` outside `-s ssa` and `--cutoff` outside `-s ctmc`.
    Each solve below passes NO options, so any exit-2 refusal means the bridge
    invented a flag from a native default.
    """
    _require_json_capable_arms()
    for ctor, build in ((SolverBA, _closed_ps), (SolverCTMC, _closed_ps),
                        (SolverSSA, _closed_ps), (SolverMAM, _open_mm1),
                        (SolverNC, _closed_ps)):
        ctor(build(), lang='cpp').getAvgTable()   # must not raise


def test_ctmc_cutoff_and_ssa_budget_do_reach_the_binary(monkeypatch):
    """The flags each engine DOES have are forwarded, checked at the argv."""
    _require_line_cli()
    from line_solver.solvers import cpp_dispatch

    seen = {}

    def fake_run(binary, cmd, timeout=None):
        seen['cmd'] = cmd
        raise RuntimeError('stop here')

    monkeypatch.setattr(cpp_dispatch, '_run_line_cli', fake_run)

    with pytest.raises(RuntimeError):
        SolverCTMC(_open_mm1(), cutoff=4, lang='cpp').getAvgTable()
    cmd = seen['cmd']
    assert cmd[cmd.index('-s') + 1] == 'ctmc'
    assert float(cmd[cmd.index('--cutoff') + 1]) == 4.0
    assert '--samples' not in cmd

    with pytest.raises(RuntimeError):
        SolverSSA(_closed_ps(), samples=12345, seed=99, lang='cpp').getAvgTable()
    cmd = seen['cmd']
    assert cmd[cmd.index('-s') + 1] == 'ssa'
    assert cmd[cmd.index('--samples') + 1] == '12345'
    assert cmd[cmd.index('--seed') + 1] == '99'
    assert '--tol' not in cmd and '--iter_max' not in cmd, (
        "a sample path is not an iteration; line-cli refuses these for -s ssa")


# --- what is not reachable through the bridge -------------------------------

def _solved(ctor, build, **kw):
    s = ctor(build(), lang='cpp', **kw)
    s.getAvgTable()
    return s


@pytest.mark.parametrize('call,label', [
    (lambda s: s.getProbAggr(0), 'CTMC.getProbAggr'),
    (lambda s: s.getProbSysAggr(), 'CTMC.getProbSysAggr'),
])
def test_state_probability_getters_delegate_under_cpp(call, label):
    """The two `getProb*Aggr` getters DELEGATE: `-a prob` is a wire contract now.

    They were listed here as unreachable, on the reading that line-cli printed
    these analyses in readable form only. That stopped being true when the arm
    landed, and a stale refusal is not a neutral thing to leave standing: it
    sends a caller away from something the port already answers. What must not
    happen instead is a NATIVE answer under lang='cpp', so the delegated value is
    pinned against the native one rather than merely required to exist.
    """
    _require_line_cli()
    delegated = call(_solved(SolverCTMC, _closed_ps))
    native = call(SolverCTMC(_closed_ps()))
    assert float(delegated) == pytest.approx(float(native), rel=1e-9), label


def test_ssa_sample_path_delegates_under_cpp():
    """`SSA.sample()` rides on `-a sample`: the bridge HAS a JSON for it now.

    It used to refuse for want of one. What comes back cannot be pinned against
    the native path -- two engines draw two sample paths -- so what is checked is
    what a trajectory must satisfy whatever the draw: it is as long as it was
    asked to be, its event times do not run backwards, the closed population is
    conserved at every step, and one seed gives one path.
    """
    _require_line_cli()
    s = SolverSSA(_closed_ps(), samples=5000, seed=1, lang='cpp')
    r = s.sample(1, 100)
    t = np.ravel(np.asarray(r.t, dtype=float))
    assert t.size > 0 and np.all(np.diff(t) >= -1e-12)
    state = np.asarray(r.state, dtype=float)
    assert state.shape[0] == 100
    again = SolverSSA(_closed_ps(), samples=5000, seed=1, lang='cpp').sample(1, 100)
    assert np.allclose(t, np.ravel(np.asarray(again.t, dtype=float)))
    assert np.allclose(state, np.asarray(again.state, dtype=float))


def test_fluid_transient_refuses_under_cpp():
    """FLD.getTranAvg() integrates a horizon the bridge cannot ask for."""
    _require_json_capable_arms()
    s = _solved(SolverFLD, _closed_ps, method='matrix')
    with pytest.raises(RuntimeError):
        s.getTranAvg()


def test_unsupported_solver_is_refused_by_name():
    """A solver with no C++ engine says so, instead of failing obscurely.

    LQNS drives binaries (`lqns`/`lqsim`) that `line-cli` genuinely does not
    carry. JMT USED TO STAND HERE and no longer does: `-s jmt` drives JSIM
    through the same JMT jar the wrapper does, so naming it as unsupported would
    be the decay this table's refusals are prone to -- and the token map is the
    authority for what lang='cpp' accepts, so it is asserted directly.
    """
    _require_line_cli()
    from line_solver.solvers.cpp_dispatch import _CPP_SOLVER_TOKENS, _cpp_token

    class _Fake(object):
        def __init__(self, name):
            self._name = name

        def getName(self):
            return self._name

    with pytest.raises(RuntimeError) as excinfo:
        _cpp_token(_Fake('LQNS'))
    msg = str(excinfo.value)
    assert 'LQNS' in msg and 'lang=' in msg
    assert _cpp_token(_Fake('JMT')) == 'jmt'
    assert 'LDES' not in _CPP_SOLVER_TOKENS, \
        "the LDES wrapper is already a subprocess client of the same engine"


def test_absent_binary_falls_back_to_python(monkeypatch):
    """An unrunnable binary degrades to native Python on every wired solver."""
    monkeypatch.setenv('LINE_CLI_BINARY', '/nonexistent/line-cli')
    for ctor in (SolverNC, SolverCTMC, SolverBA):
        tc = ctor(_closed_ps(), lang='cpp').getAvgTable()
        tp = ctor(_closed_ps()).getAvgTable()
        _assert_tables_agree(tp, tc, TOL, TOL, ctor.__name__)
