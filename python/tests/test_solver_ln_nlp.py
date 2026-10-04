"""
Tests for SolverLN method='nlp', the QD-AMVA analyzer.

The method states the whole layered model as one nonlinear program and roots it,
so nothing about it is a layered solve: no submodel is built, no fixed point is
driven, and `get_ensemble_avg` is served by `solver_ln_nlp_analyzer` instead.
What must hold anyway is the CONTRACT of a SolverLN method, and that is what is
asserted here: the table has the layered table's shape and its NaN pattern, the
metrics are mutually consistent in the way the layered ones are, the answer
tracks the layered one on a model both can serve, and a feature the algebraic
form cannot express is refused BY NAME rather than answered wrongly.
"""

import numpy as np
import pandas as pd
import pytest

from line_solver import (Activity, ActivityPrecedence, Entry, Exp, LayeredNetwork,
                         Processor, SchedStrategy, SolverLN, Task)
from line_solver.api.lqn import LqnNlpExportError


def client_server():
    """Three tasks in a chain; the processor of the first two runs near 1.0."""
    m = LayeredNetwork('cs')
    P1 = Processor(m, 'P1', 2, SchedStrategy.PS)
    P2 = Processor(m, 'P2', 3, SchedStrategy.PS)
    T1 = Task(m, 'T1', 50, SchedStrategy.REF).on(P1).setThinkTime(Exp(1.0 / 2))
    T2 = Task(m, 'T2', 50, SchedStrategy.FCFS).on(P1).setThinkTime(Exp(1.0 / 3))
    T3 = Task(m, 'T3', 25, SchedStrategy.FCFS).on(P2).setThinkTime(Exp(1.0 / 4))
    E1 = Entry(m, 'E1').on(T1)
    E2 = Entry(m, 'E2').on(T2)
    E3 = Entry(m, 'E3').on(T3)
    Activity(m, 'AS1', Exp(10)).on(T1).boundTo(E1).synchCall(E2, 1)
    Activity(m, 'AS2', Exp(20)).on(T2).boundTo(E2).synchCall(E3, 5).repliesTo(E2)
    Activity(m, 'AS3', Exp(50)).on(T3).boundTo(E3).repliesTo(E3)
    return m


def light_two_tier():
    """A single-server front end calling a four-server back end, lightly loaded."""
    m = LayeredNetwork('light')
    PC = Processor(m, 'PC', 1, SchedStrategy.PS)
    PS = Processor(m, 'PS', 4, SchedStrategy.PS)
    C = Task(m, 'Client', 2, SchedStrategy.REF).on(PC).setThinkTime(Exp(1.0 / 10))
    S = Task(m, 'Server', 6, SchedStrategy.FCFS).on(PS).setThinkTime(Exp(1e6))
    EC = Entry(m, 'EC').on(C)
    ES = Entry(m, 'ES').on(S)
    Activity(m, 'ac', Exp(50)).on(C).boundTo(EC).synchCall(ES, 1)
    Activity(m, 'as', Exp(80)).on(S).boundTo(ES).repliesTo(ES)
    return m


def fork_pair(d1=0.30, d2=0.15, nclient=4, serial=False):
    """The same work either as two concurrent branches or as one chain.

    Same demands, same population, same servers; only the precedence differs.
    The concurrent reading must be the faster of the two, and the sequential one
    is exactly what the additive composition law charges for the fork.
    """
    m = LayeredNetwork('seq' if serial else 'fork')
    PC = Processor(m, 'PC', 100, SchedStrategy.PS)
    PS = Processor(m, 'PS', 2, SchedStrategy.PS)
    C = Task(m, 'Client', nclient, SchedStrategy.REF).on(PC).setThinkTime(Exp(1.0))
    S = Task(m, 'Server', 8, SchedStrategy.FCFS).on(PS).setThinkTime(Exp(1e6))
    EC = Entry(m, 'EC').on(C)
    ES = Entry(m, 'ES').on(S)
    Activity(m, 'ac', Exp(1000.0)).on(C).boundTo(EC).synchCall(ES, 1)
    head = Activity(m, 'head', Exp(1000.0)).on(S).boundTo(ES)
    b1 = Activity(m, 'b1', Exp(1.0 / d1)).on(S)
    b2 = Activity(m, 'b2', Exp(1.0 / d2)).on(S)
    tail = Activity(m, 'tail', Exp(1000.0)).on(S).repliesTo(ES)
    if serial:
        S.add_precedence(ActivityPrecedence.Serial(head, b1, b2, tail))
    else:
        S.add_precedence(ActivityPrecedence.AndFork(head, [b1, b2]))
        S.add_precedence(ActivityPrecedence.AndJoin([b1, b2], tail))
    return m


def solve_nlp(model, **kw):
    s = SolverLN(model, method='nlp', **kw)
    s._table_silent = True
    return s, s.getAvgTable()


def by_name(df):
    return {row.Node: row for row in df.itertuples(index=False)}


def avg_by_name(solver):
    """The six ensemble arrays keyed by element name, before the table snaps them.

    `get_avg_table` sanitizes what it formats (near-tenth entries go onto the
    tenth), so an identity between metrics holds to 1e-6 in the arrays and only
    to the snapping tolerance in the table.
    """
    QN, UN, RN, TN, AN, WN = solver.get_ensemble_avg()
    out = {}
    for idx in range(solver.lqn.nidx):
        out[solver._get_hashname(idx)] = dict(
            QLen=QN[idx], Util=UN[idx], RespT=RN[idx],
            Tput=TN[idx], ArvR=AN[idx], ResidT=WN[idx])
    return out


def test_nlp_builds_no_layers():
    s = SolverLN(client_server(), method='nlp')
    assert s.lnmethod == 'nlp'
    assert s.nlayers == 0
    assert s.ensemble == []


def test_nlp_table_has_the_layered_shape():
    _, nlp = solve_nlp(client_server())
    ref = SolverLN(client_server())
    ref._table_silent = True
    lay = ref.getAvgTable()

    assert list(nlp.columns) == list(lay.columns)
    assert list(nlp['Node']) == list(lay['Node'])
    assert list(nlp['NodeType']) == list(lay['NodeType'])
    # A processor has no response time, a task none of its own, an entry no
    # residence time and nothing has an arrival rate. The method must leave the
    # same cells empty as the layered path, or a caller reading one column finds
    # numbers under a method where the other reports none.
    for col in ('QLen', 'Util', 'RespT', 'ResidT', 'ArvR', 'Tput'):
        assert (list(pd.isna(nlp[col])) == list(pd.isna(lay[col]))), col


def test_nlp_roots_the_qd_amva_laws():
    s, _ = solve_nlp(client_server())
    assert s.nlp_objective < 1e-16
    assert 'QD-AMVA' in s.nlp_source


def test_nlp_metrics_are_mutually_consistent():
    s = SolverLN(client_server(), method='nlp')
    r = avg_by_name(s)
    # a processor holds exactly the shares of the tasks that run on it
    assert r['P1']['Util'] == pytest.approx(r['T1']['Util'] + r['T2']['Util'], rel=1e-9)
    assert r['P2']['Util'] == pytest.approx(r['T3']['Util'], rel=1e-9)
    # Little at an entry, and a task holding what its entries hold
    for e, t in (('E1', 'T1'), ('E2', 'T2'), ('E3', 'T3')):
        assert r[e]['QLen'] == pytest.approx(r[e]['Tput'] * r[e]['RespT'], rel=1e-9)
        assert r[t]['QLen'] == pytest.approx(r[e]['QLen'], rel=1e-9)
        assert r[t]['Tput'] == pytest.approx(r[e]['Tput'], rel=1e-9)
    # flow across the call: AS2 issues five calls per execution
    assert r['E3']['Tput'] == pytest.approx(5.0 * r['E2']['Tput'], rel=1e-9)
    # an entry's response time is its activity's host residence plus its calls
    assert r['E2']['RespT'] == pytest.approx(
        r['T2']['ResidT'] + 5.0 * r['E3']['RespT'], rel=1e-9)


def test_nlp_tracks_the_layered_answer():
    _, nlp = solve_nlp(client_server())
    ref = SolverLN(client_server())
    ref._table_silent = True
    lay = ref.getAvgTable()
    a, b = by_name(nlp), by_name(lay)
    # QD-AMVA is an approximation of the layered fixed point, not a reproduction
    # of it, and this model runs its first processor at 97% -- where the two are
    # furthest apart. Within 10%.
    for node in ('T1', 'T2', 'T3'):
        assert a[node].Tput == pytest.approx(b[node].Tput, rel=0.10), node
        assert a[node].Util == pytest.approx(b[node].Util, rel=0.10), node


def test_nlp_charges_no_queueing_at_an_idle_multiserver():
    _, df = solve_nlp(light_two_tier())
    r = by_name(df)
    # four servers, two clients: the pool can never fill, so a request is served
    # at its bare demand and the entry's response time is exactly that.
    assert r['ES'].RespT == pytest.approx(1.0 / 80.0, rel=1e-6)


def test_nlp_honours_a_declared_lld():
    """A station declared twice as fast halves the residence it charges."""
    base, plain = solve_nlp(client_server())
    _, fast = solve_nlp(client_server(), config={'lld': {'P2': '2.0'}})
    assert by_name(fast)['E3'].RespT == pytest.approx(
        0.5 * by_name(plain)['E3'].RespT, rel=1e-6)


def test_nlp_writes_its_program_when_asked(tmp_path):
    path = tmp_path / 'cs_nlp.py'
    solve_nlp(client_server(), config={'nlp_source': str(path)})
    text = path.read_text()
    assert text.startswith('"""')
    assert 'HOST_RATE = [' in text


def test_nlp_refuses_a_feature_the_form_cannot_express():
    m = LayeredNetwork('async')
    P = Processor(m, 'P', 1, SchedStrategy.PS)
    Q = Processor(m, 'Q', 1, SchedStrategy.PS)
    T1 = Task(m, 'T1', 1, SchedStrategy.REF).on(P).setThinkTime(Exp(1.0))
    T2 = Task(m, 'T2', 1, SchedStrategy.FCFS).on(Q).setThinkTime(Exp(1e6))
    E1 = Entry(m, 'E1').on(T1)
    E2 = Entry(m, 'E2').on(T2)
    Activity(m, 'a1', Exp(2.0)).on(T1).boundTo(E1).asynchCall(E2, 1)
    Activity(m, 'a2', Exp(2.0)).on(T2).boundTo(E2).repliesTo(E2)
    s = SolverLN(m, method='nlp')
    s._table_silent = True
    with pytest.raises(LqnNlpExportError) as exc:
        s.getAvgTable()
    assert 'E2' in str(exc.value) or 'a1' in str(exc.value)


def test_nlp_refuses_a_delegated_language():
    # the JAR and line-cli have no such layered engine, and would answer with
    # their own default method under this name
    with pytest.raises(ValueError) as exc:
        SolverLN(client_server(), method='nlp', lang='java')
    assert "method='nlp'" in str(exc.value)


def test_nlp_has_no_layer_sensitivities():
    s = SolverLN(client_server(), method='nlp')
    with pytest.raises(RuntimeError) as exc:
        s.getSensitivityTable()
    assert 'no layers' in str(exc.value)


def test_nlp_is_reachable_by_its_aliases():
    from line_solver.solvers.solver_ln.solver_ln import ln_requested_method
    for alias in ('nlp', 'NLP', 'qdamva', 'qd-amva'):
        assert ln_requested_method(alias) == 'nlp'
    # and the layered names are untouched
    assert ln_requested_method('srvn') == 'srvn'
    assert ln_requested_method('flat') == 'flat.cs'


def loaded_two_tier():
    """The same two tiers, run at about three quarters of the back end.

    Light enough that both arms of a cold/warm comparison finish in a second,
    loaded enough that the layered fixed point actually has somewhere to move:
    on the idle model the program's answer is already the fixed point, so an
    agreement there says nothing about the iteration having run.
    """
    m = LayeredNetwork('loaded')
    PC = Processor(m, 'PC', 1, SchedStrategy.PS)
    PS = Processor(m, 'PS', 2, SchedStrategy.PS)
    C = Task(m, 'Client', 6, SchedStrategy.REF).on(PC).setThinkTime(Exp(1.0))
    S = Task(m, 'Server', 6, SchedStrategy.FCFS).on(PS).setThinkTime(Exp(1e6))
    EC = Entry(m, 'EC').on(C)
    ES = Entry(m, 'ES').on(S)
    Activity(m, 'ac', Exp(20)).on(C).boundTo(EC).synchCall(ES, 2)
    Activity(m, 'as', Exp(5.0)).on(S).boundTo(ES).repliesTo(ES)
    return m


def solve_layered(model, warm=False, **kw):
    """(solver, table) for the routing encoding, cold or warm started.

    'srvn.cs' rather than the default: the default resolves to 'srvn.ph', whose
    layers carry a composed phase-type law the program's means cannot seed, and
    which the solver refuses by name.
    """
    if warm:
        kw['config'] = {'warmstart': 'nlp'}
    s = SolverLN(model, method='srvn.cs', **kw)
    s._table_silent = True
    return s, s.getAvgTable()


def test_warmstart_marks_the_run():
    cold, _ = solve_layered(loaded_two_tier())
    warm, _ = solve_layered(loaded_two_tier(), warm=True)
    assert cold.warmstarted is False
    assert warm.warmstarted is True


def test_warmstart_stops_before_the_cold_floor():
    """iter_min does not apply to a warm-started run, and the run uses that."""
    cold, _ = solve_layered(loaded_two_tier())
    warm, _ = solve_layered(loaded_two_tier(), warm=True)
    floor = max(2 * cold.nlayers, (cold.options.iter_max + 3) // 4)
    # the cold arm is held AT the floor: its tolerance test passes well before,
    # and the floor is the only reason it keeps sweeping
    assert len(cold.results) > floor
    assert len(warm.results) < floor
    assert len(warm.results) < len(cold.results)


def test_warmstart_reaches_the_cold_answer():
    """At the default stopping rule the two arms agree to that rule's own width.

    Both stop when the iterate moves by less than iter_tol, so the metrics are
    allowed to differ by about that much; what would not be allowed is a
    difference that survives tightening the rule, which the next test rules out.
    """
    cold_s, cold = solve_layered(loaded_two_tier())
    _, warm = solve_layered(loaded_two_tier(), warm=True)
    a, b = by_name(cold), by_name(warm)
    tol = cold_s.options.iter_tol
    for node in ('Client', 'Server', 'EC', 'ES'):
        assert b[node].Tput == pytest.approx(a[node].Tput, rel=tol), node
        assert b[node].QLen == pytest.approx(a[node].QLen, rel=tol), node
    for node in ('EC', 'ES'):
        assert b[node].RespT == pytest.approx(a[node].RespT, rel=tol), node


def test_warmstart_does_not_move_the_fixed_point():
    """Tighten the stopping rule and the two arms land on the same point.

    This is the claim the feature rests on: the seed chooses where the iteration
    STARTS, not which point it ends at. It is worth asserting because a config
    that silently dropped the defaults once made the warm arm iterate a
    different map (interlocking off) and settle somewhere else entirely.
    """
    tight = dict(iter_tol=1e-9, iter_max=400)
    _, cold = solve_layered(loaded_two_tier(), **tight)
    _, warm = solve_layered(loaded_two_tier(), warm=True, **tight)
    a, b = by_name(cold), by_name(warm)
    for node in ('Client', 'Server', 'EC', 'ES'):
        assert b[node].Tput == pytest.approx(a[node].Tput, rel=1e-6), node
        assert b[node].QLen == pytest.approx(a[node].QLen, rel=1e-6), node


def test_warmstart_refuses_a_phase_type_encoding():
    # 'srvn' resolves to 'srvn.ph', whose layers carry a composed phase-type
    # service law; the program solves for means, so there is no distribution in
    # its answer to start those layers from
    with pytest.raises(ValueError) as exc:
        SolverLN(loaded_two_tier(), method='srvn', config={'warmstart': 'nlp'})
    assert 'srvn.ph' in str(exc.value)
    assert "method='srvn.cs'" in str(exc.value)


def test_warmstart_refuses_method_nlp():
    # the program IS the answer under method='nlp'; there is no iteration to seed
    with pytest.raises(ValueError) as exc:
        SolverLN(loaded_two_tier(), method='nlp', config={'warmstart': 'nlp'})
    assert 'warmstart' in str(exc.value)


def test_warmstart_refuses_a_delegated_language():
    with pytest.raises(ValueError) as exc:
        SolverLN(loaded_two_tier(), method='srvn.cs', lang='java',
                 config={'warmstart': 'nlp'})
    assert "lang='python'" in str(exc.value)


def test_warmstart_rejects_an_unknown_value():
    with pytest.raises(ValueError) as exc:
        SolverLN(loaded_two_tier(), method='srvn.cs', config={'warmstart': 'mva'})
    assert "'nlp'" in str(exc.value)


def test_warmstart_refuses_a_model_the_program_cannot_express():
    """An asynchronous call has no algebraic form here, so the seed is refused.

    By name and at the seeding step, rather than by falling back to a cold start
    and reporting the run as warm.
    """
    m = LayeredNetwork('async_warm')
    P = Processor(m, 'P', 1, SchedStrategy.PS)
    Q = Processor(m, 'Q', 1, SchedStrategy.PS)
    T1 = Task(m, 'T1', 1, SchedStrategy.REF).on(P).setThinkTime(Exp(1.0))
    T2 = Task(m, 'T2', 1, SchedStrategy.FCFS).on(Q).setThinkTime(Exp(1e6))
    E1 = Entry(m, 'E1').on(T1)
    E2 = Entry(m, 'E2').on(T2)
    Activity(m, 'a1', Exp(2.0)).on(T1).boundTo(E1).asynchCall(E2, 1)
    Activity(m, 'a2', Exp(2.0)).on(T2).boundTo(E2).repliesTo(E2)
    s = SolverLN(m, method='srvn.cs', config={'warmstart': 'nlp'})
    s._table_silent = True
    with pytest.raises(LqnNlpExportError):
        s.getAvgTable()


def test_a_partial_config_keeps_the_defaults():
    """Naming one config entry must not drop the other sixteen.

    A dataclass default_factory runs only when the field is OMITTED, so
    config={'warmstart': 'nlp'} used to replace the whole table and leave every
    read site on its own fallback: `interlocking` reads with a fallback of False
    against a declared default of True, so naming any entry at all turned
    interlocking off for the rest of the solve.
    """
    s = SolverLN(loaded_two_tier(), method='srvn.cs', config={'warmstart': 'nlp'})
    assert s.options.config['warmstart'] == 'nlp'
    assert s.options.config['interlocking'] is True
    assert s.options.config['layering'] == 'srvn'
    assert s.options.config['relax'] == 'fixed'

    # an explicitly named value still wins over the default it replaces
    s = SolverLN(loaded_two_tier(), method='srvn.cs', config={'interlocking': False})
    assert s.options.config['interlocking'] is False
    assert s.options.config['layering'] == 'srvn'

    # and None means the same as not passing one at all
    s = SolverLN(loaded_two_tier(), method='srvn.cs', config=None)
    assert s.options.config['interlocking'] is True
    assert s.options.config['warmstart'] is None


# ----------------------------------------------------------------------
# config['nlp_fork']: concurrency, off by default
# ----------------------------------------------------------------------

def test_nlp_refuses_a_fork_by_default():
    # The closure is an assumption about branch shape, not a reading of the
    # model, so the default stays the refusal the method shipped with.
    s = SolverLN(fork_pair(), method='nlp')
    s._table_silent = True
    with pytest.raises(LqnNlpExportError) as err:
        s.getAvgTable()
    assert 'AND-fork or AND-join' in str(err.value)


def test_nlp_fork_config_carries_the_model():
    s, df = solve_nlp(fork_pair(), config={'nlp_fork': 'exp'})
    t = by_name(df)
    assert np.isfinite(t['EC'].Tput) and t['EC'].Tput > 0
    assert 0.0 < t['PS'].Util <= 1.0 + 1e-9
    # flow still balances across the call: the server entry runs once per client
    assert t['ES'].Tput == pytest.approx(t['EC'].Tput, rel=1e-6)


def test_nlp_fork_is_faster_than_the_same_work_in_series():
    # The point of the closure. Identical demands and population, branches run
    # together rather than one after the other, so throughput must be higher
    # and the entry's response time lower.
    _, dfork = solve_nlp(fork_pair(), config={'nlp_fork': 'exp'})
    _, dseq = solve_nlp(fork_pair(serial=True))
    f, s = by_name(dfork), by_name(dseq)
    assert f['EC'].Tput > s['EC'].Tput
    assert f['ES'].RespT < s['ES'].RespT


def test_nlp_fork_sits_between_the_sum_and_the_largest_branch():
    # E[max] is bounded by the largest branch below and their sum above, so the
    # concurrent answer must sit strictly between the sequential reading and
    # the one that drops the cheaper branch entirely.
    d1, d2 = 0.30, 0.15
    _, dfork = solve_nlp(fork_pair(d1, d2), config={'nlp_fork': 'exp'})
    _, dseq = solve_nlp(fork_pair(d1, d2, serial=True))
    # the largest branch alone: the cheap branch made free
    _, dmax = solve_nlp(fork_pair(d1, 1e-9, serial=True))
    f, s, mx = by_name(dfork), by_name(dseq), by_name(dmax)
    assert s['EC'].Tput < f['EC'].Tput < mx['EC'].Tput


def test_nlp_fork_rejects_an_unknown_value():
    s = SolverLN(fork_pair(), method='nlp', config={'nlp_fork': 'max'})
    s._table_silent = True
    with pytest.raises(LqnNlpExportError) as err:
        s.getAvgTable()
    assert "'refuse' or 'exp'" in str(err.value)


def test_nlp_fork_default_is_refuse():
    s = SolverLN(fork_pair(), method='nlp', config={'nlp_warm': 3})
    assert s.options.config.get('nlp_fork', 'refuse') == 'refuse'
