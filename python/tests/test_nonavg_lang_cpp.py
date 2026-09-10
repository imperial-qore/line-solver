"""
`lang='cpp'` for the analyses that are NOT the average table.

`line-cli` emits `-o json` for every analysis it implements, each under a key
named after its own `-a`, and this suite pins the getters that ride on that:

* the CTMC generator, state space and stationary law (`-a gen`, `-a states`),
  which reach the caller through `_JavaAvgResult`'s lazy attributes and so had
  been answering with jline.jar's numbers under `lang='cpp'`;
* the aggregate state probabilities under MVA and NC (`-a prob`), which the two
  solvers compute differently ON PURPOSE and which are therefore checked against
  each solver's own native answer, never against each other;
* the exact response-time law (`-a cdf`), the declared rewards (`-a reward`), the
  transient occupancy (`-a tranprob`), a marked trajectory (`-a sample`) and the
  fluid ODE export (`-a odes`);
* every getter that is NOT reachable, which must refuse BY ITS OWN REASON rather
  than with one blanket message -- a refusal that says "avg only" when the arm
  exists for another solver would send a caller away from something the port is
  one step from answering.

The binary is resolved through LINE_CLI_BINARY or the usual search path; when it
is absent every test skips, because a missing binary is the one condition
`lang='cpp'` is allowed to degrade from.
"""

import numpy as np
import pytest

from line_solver import (Cache, ClosedClass, Delay, Exp, Network, OpenClass, Queue,
                         ReplacementStrategy, Reward, SchedStrategy, SolverCTMC,
                         SolverFLD, SolverMVA, SolverNC, SolverSSA, Sink, Source, Zipf)
from line_solver.solvers.cpp_dispatch import LineCliNotAvailable, find_line_cli


def _require_line_cli():
    try:
        return find_line_cli()
    except LineCliNotAvailable as e:
        pytest.skip("line-cli is not available: %s" % e)


def _require_json_arm(analysis, solver='ctmc', extra=(), model=None):
    """
    Skip when the resolved binary predates `-o json` on the arm under test.

    A binary from before the shared emitter exits 0 and prints the readable
    table, so the failure would otherwise surface as an unparseable answer rather
    than as the out-of-date build it is. `model` names the probe model, because an
    arm can refuse for a reason of its own -- `-a reward` on a model with no
    declared reward -- and that refusal is not a stale build.
    """
    import os
    import subprocess
    import tempfile

    from line_solver.io.linemodel_io import save_model

    binary = _require_line_cli()
    tmp = tempfile.mkdtemp(prefix='line_cpp_probe_')
    path = os.path.join(tmp, 'model.json')
    save_model(_cqn() if model is None else model(), path)
    cmd = [binary, '-f', path, '-i', 'json', '-s', solver, '-a', analysis, '-o', 'json']
    cmd += list(extra)
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out = proc.stdout.decode('utf-8', errors='replace')
    if proc.returncode != 0 or '{' not in out:
        pytest.skip("'%s' does not serve -s %s -a %s as JSON (stale build?): %s"
                    % (binary, solver, analysis,
                       proc.stderr.decode('utf-8', errors='replace').strip()[:160] or 'no object'))
    return binary


# --- models -----------------------------------------------------------------

def _cqn():
    """Delay -> PS Queue, one closed class of 2 jobs: a 3-state chain."""
    m = Network('cqn')
    d = Delay(m, 'Think')
    q = Queue(m, 'Q1', SchedStrategy.PS)
    c = ClosedClass(m, 'C1', 2, d)
    d.setService(c, Exp(1.0))
    q.setService(c, Exp(2.0))
    m.link(Network.serialRouting(d, q))
    return m


def _asym():
    """Think -> Q1 -> Q2, 3 jobs: the marginals differ station by station."""
    m = Network('asym')
    d = Delay(m, 'Think')
    q1 = Queue(m, 'Q1', SchedStrategy.PS)
    q2 = Queue(m, 'Q2', SchedStrategy.PS)
    c = ClosedClass(m, 'C1', 3, d)
    d.setService(c, Exp(1.0))
    q1.setService(c, Exp(2.0))
    q2.setService(c, Exp(0.7))
    m.link(Network.serialRouting(d, q1, q2))
    return m


def _fcfs():
    """The same shape with an FCFS queue, for the response-time law."""
    m = Network('cqn_fcfs')
    d = Delay(m, 'Think')
    q = Queue(m, 'Q1', SchedStrategy.FCFS)
    c = ClosedClass(m, 'C1', 2, d)
    d.setService(c, Exp(1.0))
    q.setService(c, Exp(2.0))
    m.link(Network.serialRouting(d, q))
    return m


def _open():
    m = Network('mm1')
    s = Source(m, 'S')
    q = Queue(m, 'Q1', SchedStrategy.FCFS)
    k = Sink(m, 'Sink')
    oc = OpenClass(m, 'O1')
    s.setArrival(oc, Exp(0.5))
    q.setService(oc, Exp(1.0))
    m.link(Network.serialRouting(s, q, k))
    return m


def _reward_model():
    m = _cqn()
    m.setReward('qlen_q1', Reward.queue_length(m.getNodeByName('Q1')))
    return m


def _cache_model():
    """Source -> LRU Cache (4 items, 1 list of 2) -> Sink, one read class."""
    m = Network('cache')
    s = Source(m, 'Source')
    c = Cache(m, 'Cache', 4, 2, ReplacementStrategy.LRU)
    k = Sink(m, 'Sink')
    init = OpenClass(m, 'InitClass', 0)
    hit = OpenClass(m, 'HitClass', 0)
    miss = OpenClass(m, 'MissClass', 0)
    s.setArrival(init, Exp(2.0))
    c.setRead(init, Zipf(1.4, 4))
    c.setHitClass(init, hit)
    c.setMissClass(init, miss)
    P = m.initRoutingMatrix()
    P.set(init, init, s, c, 1.0)
    P.set(hit, hit, c, k, 1.0)
    P.set(miss, miss, c, k, 1.0)
    m.link(P)
    return m


# --- the CTMC internals: generator, state space, stationary law --------------

def test_generator_and_state_space_come_from_the_cpp_not_the_jar():
    """
    EXACT agreement, and the point is WHOSE numbers they are.

    `getInfGen`/`getStateSpace`/`getSteadyState` read lazy attributes off the
    shared result container, whose fetch used to go to jline.jar unconditionally:
    a lang='cpp' solver answered with the JAR's generator. Equality with the
    native chain is what a correct delegation looks like here (both enumerate the
    same chain and solve pi Q = 0 exactly), so this pins the numbers; that they
    are the C++'s is pinned by the argv test below.
    """
    _require_json_arm('gen')
    cpp = SolverCTMC(_cqn(), lang='cpp')
    ref = SolverCTMC(_cqn())
    assert np.allclose(cpp.getInfGen(), ref.getInfGen(), rtol=0, atol=1e-12)
    assert np.allclose(cpp.getSteadyState(), ref.getSteadyState(), rtol=0, atol=1e-12)
    assert np.allclose(cpp.getStateSpaceAggr(), ref.getStateSpaceAggr(), rtol=0, atol=1e-12)
    space_cpp, _ = cpp.getStateSpace()
    space_ref, _ = ref.getStateSpace()
    assert np.allclose(space_cpp, space_ref, rtol=0, atol=1e-12)


def test_generator_fetch_runs_the_cpp_binary(monkeypatch):
    """The lazy fetch must dispatch on lang: no `-s ctmc -a gen` run, no claim."""
    _require_json_arm('gen')
    seen = []
    from line_solver.solvers import cpp_dispatch
    real = cpp_dispatch._run_line_cli

    def spy(binary, cmd, timeout=None):
        seen.append(list(cmd))
        return real(binary, cmd, timeout=timeout)

    monkeypatch.setattr(cpp_dispatch, '_run_line_cli', spy)
    SolverCTMC(_cqn(), lang='cpp').getInfGen()
    analyses = [c[c.index('-a') + 1] for c in seen if '-a' in c]
    assert 'gen' in analyses and 'states' in analyses, analyses


def test_event_filters_sum_to_the_off_diagonal_generator():
    """
    `eventFilt` is half of getInfGen, and the halves must add up.

    sum_a filt[a] is the off-diagonal part of Q by construction, so this checks
    the sparse triplets were rebuilt into the right matrices rather than merely
    being present.
    """
    _require_json_arm('gen')
    s = SolverCTMC(_cqn(), lang='cpp')
    Q, filt = s.getGenerator()
    assert len(filt) > 0
    total = np.zeros_like(Q)
    for f in filt:
        total += f
    offdiag = Q - np.diag(np.diag(Q))
    assert np.allclose(total, offdiag, rtol=0, atol=1e-12)


def test_a_generator_that_cannot_be_paired_with_a_state_space_is_refused(monkeypatch):
    """A state count mismatch is an error, never a silently mispaired chain."""
    _require_json_arm('gen')
    from line_solver.solvers import cpp_dispatch
    real = cpp_dispatch.analysis_via_cpp

    def shrink(solver, analysis, timeout=None, flags=None):
        p = real(solver, analysis, timeout=timeout, flags=flags)
        if analysis == 'states':
            p = dict(p)
            p['space'] = p['space'][:-1]
        return p

    monkeypatch.setattr(cpp_dispatch, 'analysis_via_cpp', shrink)
    with pytest.raises(RuntimeError, match='cannot be paired'):
        cpp_dispatch.generator_via_cpp(SolverCTMC(_cqn(), lang='cpp'))


# --- -a prob ----------------------------------------------------------------

def test_nc_prob_aggr_matches_the_native_nc():
    _require_json_arm('prob', solver='nc')
    assert SolverNC(_cqn(), lang='cpp').getProbSysAggr() == pytest.approx(
        SolverNC(_cqn()).getProbSysAggr(), rel=1e-9)
    assert SolverNC(_cqn(), lang='cpp').getProbAggr(0) == pytest.approx(
        float(np.asarray(SolverNC(_cqn()).getProbAggr(0)).reshape(-1)[0]), rel=1e-9)


def test_mva_prob_aggr_matches_the_native_mva():
    """
    Against the native MVA, and 0.36 is the reference's own figure.

    MATLAB `@SolverMVA/getProbSysAggr` on this model returns 0.36 with log
    -1.02165124753; both ports reproduce it. It is NOT SolverNC's 0.4 -- MVA fits
    a binomial to its means while NC takes a ratio of normalizing constants -- so
    the two solvers are compared each against itself and never against the other.
    """
    _require_json_arm('prob', solver='mva')
    _, p_cpp = SolverMVA(_cqn(), lang='cpp').getProbSysAggr()
    _, p_ref = SolverMVA(_cqn()).getProbSysAggr()
    assert p_cpp == pytest.approx(p_ref, rel=1e-9)
    assert p_cpp == pytest.approx(0.36, rel=1e-9)
    _, a_cpp = SolverMVA(_cqn(), lang='cpp').getProbAggr(1)
    _, a_ref = SolverMVA(_cqn()).getProbAggr(1)
    assert a_cpp == pytest.approx(a_ref, rel=1e-9)


def test_mva_and_nc_prob_disagree_by_construction():
    """Two different approximations of one quantity; neither is the other's bug."""
    _require_json_arm('prob', solver='nc')
    _, mva = SolverMVA(_cqn(), lang='cpp').getProbSysAggr()
    nc = SolverNC(_cqn(), lang='cpp').getProbSysAggr()
    assert mva != pytest.approx(nc, rel=1e-6)


def test_a_complete_set_state_is_the_one_answered_for():
    """
    The state trio travels, so the answer is the caller's state, not the default.

    Delay -> PS with 2 jobs is product form: P(Think=2, Q1=0) = 0.4 and
    P(Think=0, Q1=2) = 0.2, both exact in MATLAB. The set state must move the
    answer from the first to the second under BOTH engines, or the wire is
    carrying a state one of them ignores.
    """
    _require_json_arm('prob', solver='nc')

    def placed():
        m = _cqn()
        m.getNodeByName('Q1').setState([2])
        m.getNodeByName('Think').setState([0])
        return m

    for solver in (SolverNC, SolverCTMC):
        assert solver(placed(), lang='cpp').getProbSysAggr() == pytest.approx(0.2, rel=1e-9)
        assert solver(placed()).getProbSysAggr() == pytest.approx(0.2, rel=1e-9)
        assert solver(_cqn(), lang='cpp').getProbSysAggr() == pytest.approx(0.4, rel=1e-9)


def test_a_state_on_some_nodes_only_is_the_default_state_everywhere():
    """
    A PARTIAL state is not an initialization, and every backend must drop it.

    MATLAB reads sn.state through getState, which runs initDefault whenever
    hasInitState is false: with Q1 alone set to 2 of the 2 jobs the answer is the
    DEFAULT marking's 0.4. The three readings used to differ -- python kept Q1=2
    and emptied the reference station (0.2), and the document named only Q1 so
    the C++ rebuilt Think's default beside it and answered 0 for a joint state
    holding 4 of 2 jobs.
    """
    _require_json_arm('prob', solver='nc')

    def partial():
        m = _cqn()
        m.getNodeByName('Q1').setState([2])
        return m

    for solver in (SolverNC, SolverCTMC):
        assert solver(partial()).getProbSysAggr() == pytest.approx(0.4, rel=1e-9)
        assert solver(partial(), lang='cpp').getProbSysAggr() == pytest.approx(0.4, rel=1e-9)


def test_a_state_prior_over_several_rows_is_refused_by_node_name():
    """
    What the wire cannot carry is a MIXTURE, and the refusal names the node.

    The C++ reader takes a one-row state space as the initial state and rebuilds
    the default marking for any other, so a prior over k > 1 rows would be
    answered for one state under a question about a distribution. The reference
    weights the analysis over the prior's support instead, which is an analyzer's
    job and not a bridge's.
    """
    _require_json_arm('prob', solver='nc')
    m = _cqn()
    q = m.getNodeByName('Q1')
    q.setStateSpace(np.array([[0.0], [1.0], [2.0]]))
    q.setStatePrior(np.array([0.5, 0.25, 0.25]))
    m.getNodeByName('Think').setState([2])
    with pytest.raises(RuntimeError, match='Q1'):
        SolverNC(m, lang='cpp').getProbSysAggr()


# --- -a cdf -----------------------------------------------------------------

def test_cdf_respt_is_the_exact_law_and_reaches_one():
    """
    The tagged-chain law, not an exponential fitted to the mean.

    A response-time CDF is non-decreasing, starts at or above 0 and ends at 1;
    the exponential approximation the native getter returns satisfies the same
    three, so what is checked here is the SHAPE plus the fact that the curve is
    finely sampled (the C++ grids are the reference's 100001 points), which the
    approximation's 100 quantiles are not.
    """
    _require_json_arm('cdf')
    rd = SolverCTMC(_fcfs(), lang='cpp').getCdfRespT()
    assert rd, 'no curve was returned'
    for e in rd:
        t, p = e['t'], e['p']
        assert t.size > 1000, 'a 100-point grid is the native approximation, not the exact law'
        assert np.all(np.diff(p) >= -1e-12)
        assert p[-1] == pytest.approx(1.0, abs=1e-3)
        assert e['station'] >= 1 and e['class'] >= 1


def test_cdf_sys_respt_agrees_with_the_native_law_per_chain():
    """
    The `sysrespt` block is now READ, so this getter delegates instead of
    refusing.

    THE QUANTITY IS THE CYCLE TIME, one law per CHAIN, and both implementations
    tag the same chain and integrate the same absorption law -- so the check is
    numeric agreement, not merely shape. The grids are sampled the same way
    (10000 intervals from the same generator), so the curves are compared point
    for point after the shorter one truncates the pair.
    """
    _require_json_arm('cdf')
    cpp = SolverCTMC(_cqn(), lang='cpp').getCdfSysRespT()
    ref = SolverCTMC(_cqn()).getCdfSysRespT()
    assert cpp, 'no system law was returned'
    assert [e['chain'] for e in cpp] == [e['chain'] for e in ref]
    for a, b in zip(cpp, ref):
        n = min(a['t'].size, b['t'].size)
        assert n > 10, 'the law is a curve, not a handful of quantiles'
        assert np.allclose(a['t'][:n], b['t'][:n], rtol=1e-6, atol=1e-9)
        assert np.allclose(a['p'][:n], b['p'][:n], rtol=1e-6, atol=1e-9)


def test_cdf_respt_refuses_a_supplied_mean():
    _require_json_arm('cdf')
    with pytest.raises(RuntimeError, match='no R for it to take'):
        SolverCTMC(_fcfs(), lang='cpp').getCdfRespT(R=np.ones((2, 1)))


# --- -a reward --------------------------------------------------------------

def test_avg_reward_matches_the_native_evaluation():
    """
    A `Reward.*` template survives the wire; its expectation must agree.

    The C++ evaluates the DECLARATION against its own stationary law, so this is
    the same reward against the same distribution and the agreement is exact.
    """
    _require_json_arm('reward', model=_reward_model)
    r_cpp, n_cpp = SolverCTMC(_reward_model(), lang='cpp').getAvgReward()
    r_ref, n_ref = SolverCTMC(_reward_model()).getAvgReward()
    assert list(n_cpp) == list(n_ref)
    assert np.allclose(r_cpp, r_ref, rtol=1e-9, atol=1e-12)


def test_a_model_with_no_declared_reward_is_refused_by_name():
    _require_json_arm('reward', model=_reward_model)
    with pytest.raises(RuntimeError, match='no rewards are defined'):
        SolverCTMC(_cqn(), lang='cpp').getAvgReward()


# --- -a tranprob ------------------------------------------------------------

def test_tran_prob_sys_matches_the_native_matrix_exponential():
    """
    Agreement to the INTEGRATOR's tolerance, not to machine precision.

    The native getter takes a matrix exponential (exact up to rounding) while the
    C++ integrates the forward equation, so the two differ by the ODE tolerance:
    measured 6.5e-5 on this chain at t=2, with MATLAB's own answer
    [0.193915515413, 0.400011459066, 0.406073816757] sitting between them. A
    tolerance tight enough to reject that spread would be asserting that two
    different methods are one method.
    """
    _require_json_arm('tranprob', extra=('--tspan', '0:2'))
    cpp = SolverCTMC(_cqn(), lang='cpp').getTranProbSys(2.0)
    ref = SolverCTMC(_cqn()).getTranProbSys(2.0)
    assert cpp.shape == ref.shape
    assert np.allclose(cpp, ref, rtol=0, atol=1e-3)
    assert cpp.sum() == pytest.approx(1.0, abs=1e-6)
    assert cpp[0] == pytest.approx(0.1939155, abs=1e-3)


def test_tran_prob_refuses_an_unstated_horizon():
    _require_json_arm('tranprob', extra=('--tspan', '0:2'))
    with pytest.raises(RuntimeError, match='positive horizon'):
        SolverCTMC(_cqn(), lang='cpp').getTranProbSys(0.0)


def test_tran_prob_carries_both_label_sets_off_one_integration(monkeypatch):
    """One process, both views: the detailed and the aggregate labels."""
    _require_json_arm('tranprob', extra=('--tspan', '0:2'))
    from line_solver.solvers import cpp_dispatch
    runs = []
    real = cpp_dispatch._run_line_cli
    monkeypatch.setattr(cpp_dispatch, '_run_line_cli',
                        lambda b, c, timeout=None: (runs.append(list(c)),
                                                    real(b, c, timeout=timeout))[1])
    out = cpp_dispatch.tran_prob_via_cpp(SolverCTMC(_cqn(), lang='cpp'), 2.0)
    assert len(runs) == 1, 'two invocations for one horizon'
    assert out['pit'].shape[0] == out['t'].size
    assert out['labels'].shape[0] == out['pit'].shape[1]
    assert out['labelsAggr'].shape[0] == out['pit'].shape[1]


# --- -a sample --------------------------------------------------------------

def test_sample_sys_is_served_although_the_native_network_path_refuses_it():
    """
    MATLAB has @SolverCTMC/sampleSys and native Python does not; the port does.

    So lang='cpp' answers a getter that raises NotImplementedError natively, and
    the trajectory it returns is checked as a trajectory: non-decreasing times,
    every state a row of the chain's space, and the closed population conserved
    at every step.
    """
    _require_json_arm('sample', extra=('--samples', '5'))
    with pytest.raises(NotImplementedError):
        SolverCTMC(_cqn()).sampleSys(10)
    r = SolverCTMC(_cqn(), lang='cpp', seed=4242).sampleSys(40)
    assert r.numEvents == 40
    assert np.all(np.diff(r.t) >= -1e-12)
    assert r.state.shape[0] == 40
    assert np.allclose(r.state.sum(axis=1), 2.0)


def test_sample_sys_is_reproducible_under_one_seed():
    _require_json_arm('sample', extra=('--samples', '5'))
    a = SolverCTMC(_cqn(), lang='cpp', seed=99).sampleSys(30)
    b = SolverCTMC(_cqn(), lang='cpp', seed=99).sampleSys(30)
    assert np.allclose(a.t, b.t) and np.allclose(a.state, b.state)
    c = SolverCTMC(_cqn(), lang='cpp', seed=100).sampleSys(30)
    assert not np.allclose(a.t, c.t), 'two seeds produced one trace'


def test_sample_aggr_narrows_the_same_walk_to_one_node():
    _require_json_arm('sample', extra=('--samples', '5'))
    r = SolverCTMC(_cqn(), lang='cpp', seed=7).sampleAggr(1, 25)
    assert r.isaggregate and r.nodeIndex == 1
    assert r.state.shape[0] == 25


# --- -a odes ----------------------------------------------------------------

def test_export_odes_returns_the_cpp_document_in_both_notations():
    """
    `--notation` reaches the exporter: the matrix document is not the scalar one.

    The method must be named explicitly because `default` resolves to the
    second-order closure `minnormal` on this model, whose symbolic export the C++
    has not ported -- and refuses by name rather than exporting another method's
    drift.
    """
    _require_json_arm('odes', solver='fluid', extra=('--method', 'matrix'))
    scalar = SolverFLD(_cqn(), method='matrix', lang='cpp').exportODEs()
    matrix = SolverFLD(_cqn(), method='matrix', lang='cpp').exportODEs(notation='matrix')
    assert '\\documentclass' in scalar and '\\documentclass' in matrix
    assert 'scalar' in scalar and scalar != matrix


def test_export_odes_refuses_an_unported_method_by_name():
    _require_json_arm('odes', solver='fluid', extra=('--method', 'matrix'))
    with pytest.raises(RuntimeError, match='minnormal'):
        SolverFLD(_cqn(), lang='cpp').exportODEs()


# --- -a cache ---------------------------------------------------------------

def test_cache_table_comes_from_the_cpp_and_matches_the_native_one():
    """
    `-a cache` is a real arm, so the getter DELEGATES rather than refuses.

    The native builder reads the hit and miss probabilities off the Cache NODE
    OBJECTS, which only a native solve writes; under lang='cpp' nothing writes
    them, so a table built here would report a plain cache for any model. The
    numbers are pinned against the native solve of the same chain: both enumerate
    it exactly, so the two tables must agree to machine precision.
    """
    _require_json_arm('cache', model=_cache_model)
    native = SolverCTMC(_cache_model()).getAvgCacheTable()
    delegated = SolverCTMC(_cache_model(), lang='cpp').getAvgCacheTable()
    assert list(delegated['Node']) == list(native['Node'])
    assert list(delegated['JobClass']) == list(native['JobClass'])
    for col in ('HitProb', 'DelayedHitProb', 'MissProb', 'HitRate', 'MissRate', 'ArvR'):
        np.testing.assert_allclose(np.asarray(delegated[col], dtype=float),
                                   np.asarray(native[col], dtype=float), rtol=1e-10)


# --- the refusals, each with its own reason ---------------------------------

@pytest.mark.parametrize('call,fragment', [
    (lambda: SolverCTMC(_cqn(), lang='cpp').getTranAvg(), 'finite horizon'),
    (lambda: SolverCTMC(_cqn(), lang='cpp').getAvgCacheTable(), 'no Cache node'),
    (lambda: SolverSSA(_cqn(), lang='cpp').getCdfRespT(), 'refuses -s ssa -a cdf'),
])
def test_an_unreachable_getter_refuses_with_its_own_reason(call, fragment):
    _require_line_cli()
    with pytest.raises(RuntimeError, match=fragment):
        call()


def test_the_joint_state_probability_delegates_rather_than_refusing():
    """
    `getProb`/`getProbSys` ride on `-a prob` like their aggregate twins.

    They used to refuse with "no SolverCTMC getProb function", which stopped
    being true when the arm landed. The C++ enumerates the same chain, so the
    numbers are pinned against the native ones rather than merely shape-checked;
    a station-asymmetric model is used because the symmetric one answers 0.4 for
    every station and would hide an index that never moved.
    """
    _require_json_arm('prob')
    for i in (0, 1, 2):
        cpp = SolverCTMC(_asym(), lang='cpp').getProb(i)
        assert cpp == pytest.approx(SolverCTMC(_asym()).getProb(i), rel=1e-9)
    assert SolverCTMC(_asym(), lang='cpp').getProbSys() == pytest.approx(
        SolverCTMC(_asym()).getProbSys(), rel=1e-9)


def test_the_fluid_passage_law_rides_on_the_cpp_state_vector():
    """
    `-a statevec` is what makes the passage-time getters reachable under cpp.

    `getCdfPassT` indexes the solved ODE state vector, which is fetched from the
    C++ and no longer refused as having no arm; the law it integrates from there
    must be the native one, curve for curve.
    """
    _require_json_arm('statevec', solver='fluid', extra=('--method', 'matrix'), model=_open)
    cpp = SolverFLD(_open(), method='matrix', lang='cpp').getCdfPassT()
    ref = SolverFLD(_open(), method='matrix').getCdfPassT()
    pairs = [(a, b) for ra, rb in zip(cpp, ref) for a, b in zip(ra, rb)
             if a is not None and b is not None]
    assert pairs, 'no passage-time law was returned'
    for a, b in pairs:
        a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
        assert a.shape == b.shape
        np.testing.assert_allclose(a, b, rtol=1e-6, atol=1e-9)


def test_no_refusal_still_claims_the_bridge_serves_avg_only():
    """
    The blanket reason is gone, and it must not come back.

    Every analysis now emits JSON, so "delegates the average analysis only" would
    be false -- and a false reason is worse than a terse one: it tells a caller to
    stop asking for something that is one arm away. `getTranAvg` is the probe
    because it is one of the few getters still genuinely out of reach, so what it
    says is the shape every remaining refusal has to take.
    """
    _require_line_cli()
    with pytest.raises(RuntimeError) as e:
        SolverCTMC(_cqn(), lang='cpp').getTranAvg()
    assert 'average analysis only' not in str(e.value)


def test_the_transient_marginal_is_bucketed_off_the_cpp_law():
    """
    `getTranProb` DELEGATES: `-a tranprob` carries the full occupancy pi(t).

    It used to refuse, on the reading that the arm's shape (the whole labelled
    state space) was not this getter's (one station's marginal). The bucketing IS
    the getter, so it is applied to the C++'s law over the C++'s own enumeration,
    and the result is pinned against the native matrix exponential -- to the
    integrator's tolerance, the same spread `getTranProbSys` is checked at.
    """
    _require_json_arm('tranprob', extra=('--tspan', '0:1'))
    for node in (0, 1):
        cpp = SolverCTMC(_cqn(), lang='cpp').getTranProb(node, 1.0)
        ref = SolverCTMC(_cqn()).getTranProb(node, 1.0)
        assert cpp.shape == ref.shape
        assert np.allclose(cpp, ref, rtol=0, atol=1e-3)
        assert cpp.sum() == pytest.approx(1.0, abs=1e-6)


def test_the_nc_queue_length_law_comes_from_the_cpp_marginal():
    """
    `SolverNC.getProbMarg` rides on `-a marg`, DELAY STATION INCLUDED.

    The refusal it used to carry read the arm as a per-class quantity, which is
    not what this getter asks for: the total queue-length law is the sum of the
    aggregate marginal over the per-class partitions of n, and that is what both
    engines compute. The delay is checked beside the queue because it is the one
    the native path used to answer with zeros -- procomom has no row for it --
    and a station holding jobs with probability one cannot have a zero law.
    """
    _require_json_arm('marg', solver='nc')
    for ist in (0, 1):
        cpp = SolverNC(_cqn(), lang='cpp').getProbMarg(ist)
        ref = SolverNC(_cqn()).getProbMarg(ist)
        np.testing.assert_allclose(cpp, ref, rtol=1e-9, atol=1e-12)
        assert float(np.sum(cpp)) == pytest.approx(1.0, abs=1e-9)
    # Think(Z=1) against a PS queue(rho=0.5) at N=2: the product form is exact
    np.testing.assert_allclose(SolverNC(_cqn(), lang='cpp').getProbMarg(0),
                               [0.2, 0.4, 0.4], rtol=1e-9)


def test_the_nc_response_time_law_rides_on_the_cpp_cdf_arm():
    """
    `SolverNC.getCdfRespT` DELEGATES: `-s nc -a cdf` is the same `pfqn_stdf` law.

    Its refusal named a SolverCTMC analysis, which was never what this getter
    runs -- NC's is the exact product-form sojourn law at FCFS stations. An FCFS
    model is used because a PS queue has no entry in either codebase, so the
    empty answer would agree without either law being computed.
    """
    _require_json_arm('cdf', solver='nc', model=_fcfs)
    cpp = {(r['station'], r['class']): r for r in SolverNC(_fcfs(), lang='cpp').getCdfRespT()}
    ref = {(r['station'], r['class']): r for r in SolverNC(_fcfs()).getCdfRespT()}
    assert cpp and sorted(cpp) == sorted(ref)
    for k in sorted(cpp):
        np.testing.assert_allclose(np.ravel(np.asarray(cpp[k]['t'], dtype=float)),
                                   np.ravel(np.asarray(ref[k]['t'], dtype=float)),
                                   rtol=1e-9, atol=1e-12)
        np.testing.assert_allclose(np.ravel(np.asarray(cpp[k]['p'], dtype=float)),
                                   np.ravel(np.asarray(ref[k]['p'], dtype=float)),
                                   rtol=1e-9, atol=1e-12)


def test_open_model_prob_matches_the_product_form():
    """
    A mixed/open model must reach the product-form branch, not an OverflowError.

    `schmidt_binomial_prob_aggr` ran `int(N[r])` on an infinite population, so
    every open class raised before any law was applied. M/M/1 at rho=0.5 in the
    empty state has P = 1 - rho exactly, in MATLAB, native Python and the C++.
    """
    _require_json_arm('prob', 'mva', model=_open)

    native = SolverMVA(_open())
    cpp = SolverMVA(_open(), lang='cpp')
    for solver in (native, cpp):
        # station 1 is the Source: EXT carries no population of its own
        assert solver.getProbAggr(1)[1] == pytest.approx(1.0)
        assert solver.getProbAggr(2)[1] == pytest.approx(0.5)
        assert solver.getProbSysAggr()[1] == pytest.approx(0.5)


def test_fluid_prob_aggr_is_the_scalar_state_probability():
    """
    `SolverFLD.getProbAggr` returns P(state), not a distribution.

    It used to return a geometric VECTOR `(1-rho) rho^n` off a mean utilization,
    with a 0-based station argument, where `@SolverFLD/getProbAggr.m` returns the
    scalar probability of the state the model is in and every sibling getter is
    1-based. Same law as MVA's, on the fluid's Q and U.
    """
    logp, p = SolverFLD(_cqn()).getProbAggr(1)
    assert isinstance(p, float)
    assert 0.0 < p < 1.0
    assert logp == pytest.approx(float(np.log(p)))
    # both stations hold the same closed population, so the binomial agrees
    assert SolverFLD(_cqn()).getProbAggr(2)[1] == pytest.approx(p)
    # and the open branch is reachable, at the product-form value
    assert SolverFLD(_open()).getProbAggr(2)[1] == pytest.approx(0.5)
