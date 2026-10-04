"""
Tests for lqn_mol and for config['warmstart']='mol'.

`lqn_mol` is the Method of Layers on the SRVN decomposition, a compact layered
fixed point that solves every submodel with `pfqn_qdamva` instead of building a
Network. Two things are asserted about it here: that it READS THE STRUCT IN THE
INDEX BASE THIS BRANCH USES, which the copy taken from master did not and which
cost it every model in the suite, and that its answer tracks SolverLN's on a
model both serve.

`config['warmstart']='mol'` then writes that answer onto the layered iterate, the
second starting point beside the program's. What must hold is the contract of a
starting point: the run is marked warm, it may stop before the cold `iter_min`
floor, it lands where the cold run lands, and a model or an encoding the seed
cannot serve is refused BY NAME rather than half-seeded and reported as warm.
"""

import numpy as np
import pytest

from line_solver import (Activity, ActivityPrecedence, Entry, Exp, LayeredNetwork,
                         Processor, SchedStrategy, SolverLN, Task)
from line_solver.api.lqn import lqn_mol


def entry_only():
    """Two tiers, entry-only, run at about three quarters of the back end.

    Every entry binds exactly one activity and there is no precedence, which is
    the scope `lqn_mol` serves. Loaded enough that the layered fixed point has
    somewhere to move: on an idle model the seed is already the fixed point, so
    agreement there would say nothing about the iteration having run.
    """
    m = LayeredNetwork('entry_only')
    PC = Processor(m, 'PC', 1, SchedStrategy.PS)
    PS = Processor(m, 'PS', 2, SchedStrategy.PS)
    C = Task(m, 'Client', 6, SchedStrategy.REF).on(PC).setThinkTime(Exp(1.0))
    S = Task(m, 'Server', 6, SchedStrategy.FCFS).on(PS).setThinkTime(Exp(1e6))
    EC = Entry(m, 'EC').on(C)
    ES = Entry(m, 'ES').on(S)
    Activity(m, 'ac', Exp(20)).on(C).boundTo(EC).synchCall(ES, 2)
    Activity(m, 'as', Exp(5.0)).on(S).boundTo(ES).repliesTo(ES)
    return m


def activity_graph():
    """The same two tiers with a two-activity sequence, which lqn_mol refuses."""
    m = LayeredNetwork('graph')
    PC = Processor(m, 'PC', 1, SchedStrategy.PS)
    PS = Processor(m, 'PS', 2, SchedStrategy.PS)
    C = Task(m, 'Client', 6, SchedStrategy.REF).on(PC).setThinkTime(Exp(1.0))
    S = Task(m, 'Server', 6, SchedStrategy.FCFS).on(PS).setThinkTime(Exp(1e6))
    EC = Entry(m, 'EC').on(C)
    ES = Entry(m, 'ES').on(S)
    Activity(m, 'ac', Exp(20)).on(C).boundTo(EC).synchCall(ES, 2)
    a1 = Activity(m, 'as1', Exp(5.0)).on(S).boundTo(ES)
    a2 = Activity(m, 'as2', Exp(5.0)).on(S).repliesTo(ES)
    S.add_precedence(ActivityPrecedence.Serial(a1, a2))
    return m


def solve_layered(model, warm=None, **kw):
    """(solver, table) for the routing encoding, cold or warm started.

    'srvn.cs' rather than the default: the default resolves to 'srvn.ph', whose
    layers carry a composed phase-type law no seed of means can parameterize,
    and which the solver refuses by name.
    """
    if warm is not None:
        kw['config'] = {'warmstart': warm}
    s = SolverLN(model, method='srvn.cs', **kw)
    s._table_silent = True
    return s, s.getAvgTable()


def by_name(df):
    return {str(r['Node']): r for _, r in df.iterrows()}


# --------------------------------------------------------------------------
# lqn_mol itself
# --------------------------------------------------------------------------

def test_lqn_mol_reads_the_one_based_struct():
    """The base mismatch that cost lqn_mol every model in the suite.

    On this branch LayeredNetworkStruct is 1-BASED, slot 0 unused, as MATLAB's;
    master's is 0-based, and the port was written there. Reading a band from
    `shift` instead of `shift+1` walks the TASK band while believing it is the
    entry band, so the scope assertion reports "entry <a task name> binds k
    activities" and refuses a model that is in fact entry-only. The names below
    are a task and an entry, and only the entry may ever be named that way.
    """
    lsn = entry_only().getStruct()
    QN, UN, RN, TN, info = lqn_mol(lsn)

    names = [str(n) for n in lsn.hashnames]
    ec, es = names.index('EC'), names.index('ES')
    cl, sv = names.index('Client'), names.index('Server')

    # the entry band, which carries every reported metric
    for i in (ec, es):
        assert np.isfinite(QN[i]) and np.isfinite(UN[i])
        assert np.isfinite(RN[i]) and np.isfinite(TN[i])
        assert RN[i] > 0 and TN[i] > 0
    # the task band, which carries a throughput and no response time
    for i in (cl, sv):
        assert np.isfinite(TN[i]) and np.isnan(RN[i])
    # a processor has a utilization and nothing else
    for i in (names.index('PC'), names.index('PS')):
        assert np.isfinite(UN[i]) and np.isnan(TN[i]) and np.isnan(QN[i])

    assert info['resid'] < 1e-6
    assert 0 < info['iter'] < 200


def test_lqn_mol_state_is_the_iterate_solver_ln_carries():
    """The four vectors the warm start copies are all present and positive."""
    lsn = entry_only().getStruct()
    _, _, _, _, info = lqn_mol(lsn)
    for key in ('servt', 'residt', 'callservt', 'thinkt'):
        assert key in info
    names = [str(n) for n in lsn.hashnames]
    es = names.index('ES')
    # a called entry's service time is its host residence plus what it calls,
    # and it makes no calls here, so the two coincide and are the demand blown
    # up by contention
    assert info['servt'][es] == pytest.approx(info['residt'][es])
    assert info['servt'][es] >= 1.0 / 5.0
    assert float(info['callservt'][0]) > 0
    assert float(info['thinkt'][names.index('Server')]) >= 0


def test_lqn_mol_tracks_solver_ln():
    """Its throughputs and processor utilizations land on the layered answer.

    Not to machine precision: a submodel here carries one class per client TASK
    with visit-weighted demands where 'srvn.cs' carries one class per activity
    and encodes the call multiplicities as routing, so the two aggregate
    differently. On an entry-only model that difference is small at the
    throughputs and the processors.
    """
    lsn = entry_only().getStruct()
    _, UN, _, TN, _ = lqn_mol(lsn)
    names = [str(n) for n in lsn.hashnames]

    _, df = solve_layered(entry_only(), iter_tol=1e-9, iter_max=400)
    rows = by_name(df)
    for n in ('EC', 'ES'):
        assert TN[names.index(n)] == pytest.approx(float(rows[n]['Tput']), rel=0.05)
    for n in ('PC', 'PS'):
        assert UN[names.index(n)] == pytest.approx(float(rows[n]['Util']), abs=0.05)


def test_lqn_mol_refuses_an_activity_graph_by_name():
    with pytest.raises(ValueError) as exc:
        lqn_mol(activity_graph().getStruct())
    msg = str(exc.value)
    assert 'lqn_mol' in msg
    assert 'ES' in msg or 'as1' in msg or 'as2' in msg
    assert 'SolverLN' in msg


# --------------------------------------------------------------------------
# config['warmstart']='mol'
# --------------------------------------------------------------------------

def test_warmstart_mol_marks_the_run():
    cold, _ = solve_layered(entry_only())
    warm, _ = solve_layered(entry_only(), warm='mol')
    assert cold.warmstarted is False
    assert warm.warmstarted is True


def test_warmstart_mol_stops_before_the_cold_floor():
    """iter_min does not apply to a warm-started run, so it may stop early."""
    cold, _ = solve_layered(entry_only())
    warm, _ = solve_layered(entry_only(), warm='mol')
    floor = max(2 * cold.nlayers, (cold.options.iter_max + 3) // 4)
    assert len(cold.results) > floor
    assert len(warm.results) < len(cold.results)


def test_warmstart_mol_reaches_the_cold_answer():
    """A starting point moves the path, not the destination."""
    cold_s, cold = solve_layered(entry_only())
    _, warm = solve_layered(entry_only(), warm='mol')
    tol = cold_s.options.iter_tol
    c, w = by_name(cold), by_name(warm)
    for n in c:
        for col in ('QLen', 'Util', 'RespT', 'Tput'):
            a, b = float(c[n][col]), float(w[n][col])
            if np.isnan(a) or np.isnan(b):
                assert np.isnan(a) and np.isnan(b)
                continue
            assert abs(a - b) <= tol * max(1.0, abs(a))


def test_warmstart_mol_does_not_move_the_fixed_point():
    """Driven to the same tight tolerance, both arms land on one answer."""
    _, cold = solve_layered(entry_only(), iter_tol=1e-9, iter_max=400)
    _, warm = solve_layered(entry_only(), warm='mol', iter_tol=1e-9, iter_max=400)
    c, w = by_name(cold), by_name(warm)
    for n in c:
        for col in ('QLen', 'Util', 'RespT', 'Tput'):
            a, b = float(c[n][col]), float(w[n][col])
            if np.isnan(a) or np.isnan(b):
                continue
            assert a == pytest.approx(b, rel=1e-6, abs=1e-9)


def test_warmstart_mol_agrees_with_the_program_seed():
    """The two seeds are two routes to one law, so they land together.

    Both price a station by QD-AMVA: the program roots its stationarity
    conditions in one shot, the Method of Layers sweeps to them. Neither moves
    the layered fixed point, so the two warm runs agree to the layered run's own
    stopping tolerance.
    """
    s, mol = solve_layered(entry_only(), warm='mol')
    _, nlp = solve_layered(entry_only(), warm='nlp')
    m, n = by_name(mol), by_name(nlp)
    for name in m:
        for col in ('QLen', 'Util', 'RespT', 'Tput'):
            a, b = float(m[name][col]), float(n[name][col])
            if np.isnan(a) or np.isnan(b):
                continue
            assert abs(a - b) <= s.options.iter_tol * max(1.0, abs(a))


def test_warmstart_mol_refuses_a_phase_type_encoding():
    """Composed phase-type layers have no distribution for a mean to seed."""
    with pytest.raises(ValueError) as exc:
        SolverLN(entry_only(), method='srvn.ph',
                 config={'warmstart': 'mol'}).getAvgTable()
    msg = str(exc.value)
    assert 'PHASE-TYPE' in msg
    assert "'mol'" in msg


def test_warmstart_mol_refuses_method_nlp():
    """method='nlp' runs no iteration, so there is nothing to start."""
    with pytest.raises(ValueError) as exc:
        SolverLN(entry_only(), method='nlp', config={'warmstart': 'mol'})
    assert 'LAYERED' in str(exc.value)


def test_warmstart_mol_refuses_a_delegated_language():
    with pytest.raises(ValueError) as exc:
        SolverLN(entry_only(), method='srvn.cs', lang='java',
                 config={'warmstart': 'mol'})
    assert "lang='python'" in str(exc.value)


def test_warmstart_mol_refuses_a_model_it_cannot_express():
    """An activity graph is outside lqn_mol's scope, and is named, not ignored."""
    with pytest.raises(ValueError) as exc:
        solve_layered(activity_graph(), warm='mol')
    assert 'lqn_mol' in str(exc.value)


def test_warmstart_mol_is_built_to_its_own_budget():
    """The seed answers to mol_iter_max, not to the layered run's iter_max.

    A seed handed over at LN's own 5e-3 stopping rule would still be two digits
    out, so the Method of Layers is driven to its own 1e-6 instead. Cutting its
    budget to a single sweep must leave a DIFFERENT iterate at the end of
    init(), which is what shows the knob is read; it need not cost layered
    sweeps afterwards, and on a small model it does not.
    """
    def seeded(**cfg):
        s = SolverLN(entry_only(), method='srvn.cs',
                     config=dict(warmstart='mol', **cfg))
        s._table_silent = True
        s._construct()
        s.init()
        return s

    tight, loose = seeded(), seeded(mol_iter_max=1)
    assert tight.warmstarted is loose.warmstarted is True
    assert not np.allclose(np.asarray(tight.servt, dtype=float),
                           np.asarray(loose.servt, dtype=float))

    # and both still reach the one fixed point, the budget being a path and not
    # a destination
    _, a = solve_layered(entry_only(), warm='mol', iter_tol=1e-9, iter_max=400)
    s = SolverLN(entry_only(), method='srvn.cs', iter_tol=1e-9, iter_max=400,
                 config={'warmstart': 'mol', 'mol_iter_max': 1})
    s._table_silent = True
    b = s.getAvgTable()
    for n, r in by_name(a).items():
        for col in ('QLen', 'Util', 'RespT', 'Tput'):
            x, y = float(r[col]), float(by_name(b)[n][col])
            if np.isnan(x) or np.isnan(y):
                continue
            assert x == pytest.approx(y, rel=1e-6, abs=1e-9)
