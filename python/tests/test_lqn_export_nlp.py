"""
Tests for line_solver.api.lqn.export_nlp, the LQN-to-nonlinear-program exporter.

The contract has three halves, and all three are asserted here: the emitted
script is SELF-CONTAINED (it runs with numpy and scipy alone), it is SQUARE
(free dimensions after the linear block equal the residual count, so a zero
objective is an isolated solution rather than a point on a manifold), and the
point it converges to satisfies the laws the program is built from. Features the
algebraic form cannot express must be refused BY NAME rather than dropped.
"""

import numpy as np
import pytest

from line_solver import (Activity, ActivityPrecedence, Entry, Exp, Immediate,
                         LayeredNetwork, Processor, SchedStrategy, Task)
from line_solver.api.lqn import LqnNlpExportError, export_nlp


def client_server():
    """Three tasks in a chain, the model the balance-equation tests use."""
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


def idle_multiserver():
    """A four-server back end that is nowhere near busy, plus a delay disk."""
    m = LayeredNetwork('idle')
    PC = Processor(m, 'PC', 1, SchedStrategy.PS)
    PS = Processor(m, 'PS', 4, SchedStrategy.PS)
    PD = Processor(m, 'PD', float('inf'), SchedStrategy.INF)
    C = Task(m, 'Client', 2, SchedStrategy.REF).on(PC).setThinkTime(Exp(1.0 / 10))
    S = Task(m, 'Server', 6, SchedStrategy.FCFS).on(PS).setThinkTime(Exp(1e6))
    D = Task(m, 'Disk', float('inf'), SchedStrategy.INF).on(PD).setThinkTime(Exp(1e6))
    EC = Entry(m, 'EC').on(C)
    ES = Entry(m, 'ES').on(S)
    ED = Entry(m, 'ED').on(D)
    Activity(m, 'ac', Exp(50)).on(C).boundTo(EC).synchCall(ES, 1)
    Activity(m, 'as', Exp(80)).on(S).boundTo(ES).synchCall(ED, 1).repliesTo(ES)
    Activity(m, 'ad', Exp(100)).on(D).boundTo(ED).repliesTo(ED)
    return m


def forked(d1=0.30, d2=0.15, nclient=3, host=float('inf'), quorum=None,
           join=True):
    """A server entry that forks two branches and joins them before replying.

    The server's processor is an infinite server by default, so W(a) = D(a)
    exactly and the entry's composition can be checked against the closure with
    no congestion term in the way.
    """
    m = LayeredNetwork('fork')
    PC = Processor(m, 'PC', 100, SchedStrategy.PS)
    PS = Processor(m, 'PS', host,
                   SchedStrategy.INF if host == float('inf') else SchedStrategy.PS)
    C = Task(m, 'Client', nclient, SchedStrategy.REF).on(PC).setThinkTime(Exp(1.0))
    S = Task(m, 'Server', 8, SchedStrategy.FCFS).on(PS).setThinkTime(Exp(1e6))
    EC = Entry(m, 'EC').on(C)
    ES = Entry(m, 'ES').on(S)
    Activity(m, 'ac', Exp(1000.0)).on(C).boundTo(EC).synchCall(ES, 1)
    head = Activity(m, 'head', Exp(1.0 / 0.01)).on(S).boundTo(ES)
    b1 = Activity(m, 'b1', Exp(1.0 / d1) if d1 else Immediate()).on(S)
    b2 = Activity(m, 'b2', Exp(1.0 / d2) if d2 else Immediate()).on(S)
    tail = Activity(m, 'tail', Exp(1.0 / 0.02)).on(S).repliesTo(ES)
    S.add_precedence(ActivityPrecedence.AndFork(head, [b1, b2]))
    if join:
        S.add_precedence(ActivityPrecedence.AndJoin([b1, b2], tail, quorum))
    return m


def run_emitted(model, **kw):
    """Exec the emitted script in a private namespace and solve it there."""
    src = export_nlp(model, **kw)
    ns = {'__name__': 'emitted_nlp'}
    exec(compile(src, '<emitted_nlp>', 'exec'), ns)
    res = ns['solve']()
    return ns, np.maximum(res.x, 0.0)


def block(ns, x, name):
    n = dict(ns['_BLOCKS'])[name]
    off = ns['OFFSET'][name]
    return x[off:off + n]


def test_emitted_script_is_self_contained():
    src = export_nlp(client_server())
    assert 'import line_solver' not in src
    assert 'from line_solver' not in src
    assert 'import numpy as np' in src
    assert 'from scipy.optimize import' in src
    # the model's own parameters must appear as literals, not be baked into a matrix
    assert "ACT_DEMAND = [0.1, 0.05, 0.02]" in src
    assert "CALL_MEAN = [1.0, 5.0]" in src
    assert "TASK_MULT = [50.0, 50.0, 25.0]" in src
    assert "TASK_THINK = [2.0, 0.0, 0.0]" in src
    assert "HOST_MULT = [2.0, 3.0]" in src


def test_program_is_square():
    # Free dimensions after the linear equalities must equal the residual count.
    ns = {'__name__': 'emitted_nlp'}
    exec(compile(export_nlp(client_server()), '<emitted_nlp>', 'exec'), ns)
    rank = np.linalg.matrix_rank(ns['AEQ'])
    assert rank == ns['AEQ'].shape[0]           # no redundant equality
    assert ns['NVAR'] - rank == ns['NRES']


def test_client_server_satisfies_its_own_laws():
    ns, x = run_emitted(client_server())
    assert ns['objective'](x) < 1e-6
    assert np.max(np.abs(ns['AEQ'] @ x - ns['BEQ'])) < 1e-8
    assert np.max(ns['AUB'] @ x - ns['BUB']) < 1e-6

    Xe, Xa, Xc = (block(ns, x, k) for k in ('Xe', 'Xa', 'Xc'))
    W, S, R = (block(ns, x, k) for k in ('W', 'S', 'R'))
    Q, Qt, BR = (block(ns, x, k) for k in ('Q', 'Qt', 'BR'))

    # flow: AS2 issues five calls per execution, so E3 runs five times per E2
    assert Xe[2] == pytest.approx(5.0 * Xe[1], rel=1e-9)
    assert Xc[1] == pytest.approx(5.0 * Xa[1], rel=1e-9)

    # Little at every host, and at every call's task
    assert np.allclose(Q, Xa * W, rtol=1e-6, atol=1e-9)
    assert np.allclose(Qt, Xc * R, rtol=1e-6, atol=1e-9)

    # the reference task's population is split between thinking and its own entry
    assert 2.0 * Xe[0] + BR[0] == pytest.approx(50.0, rel=1e-9)

    # nothing is busier than it has servers
    util_p1 = sum(ns['ACT_DEMAND'][a] * Xa[a] for a in ns['ACTS_OF_HOST'][0]) / 2.0
    assert 0.0 < util_p1 <= 1.0 + 1e-9

    # an entry's service time is its host time plus the calls it dispatches
    assert S[1] == pytest.approx(W[1] + 5.0 * R[1], rel=1e-9)


def test_idle_multiserver_charges_no_waiting():
    # A four-server pool holding well under three jobs queues nothing, so the
    # residence time is the bare demand. QD-AMVA gets this from the rate alone:
    # r(n) = softmin(n, 4) is n while servers are spare, so W = D*n/n = D. No
    # clamp, no slack, no branch.
    ns, x = run_emitted(idle_multiserver())
    assert ns['objective'](x) < 1e-6
    W = block(ns, x, 'W')
    Q = block(ns, x, 'Q')
    srv = ns['ACTS'].index('as')
    assert ns['qbar_host'](Q, 1, srv) < 3.0
    assert W[srv] == pytest.approx(ns['ACT_DEMAND'][srv], rel=1e-6)

    # the pool is serving at the population present, not at its server count
    present = ns['qd_pop_host'](Q, 1)
    assert present < 3.0
    assert block(ns, x, 'GH')[ns['GH_OF'][1]] == pytest.approx(present, rel=1e-6)


def test_call_into_delay_task_never_queues():
    # A call into an INF task has R = S as a LINEAR constraint, not a residual,
    # so the task carries no queue-dependent rate at all.
    ns, x = run_emitted(idle_multiserver())
    disk = ns['ENTRIES'].index('ED')
    c = [i for i in range(ns['NC']) if ns['CALL_ENTRY'][i] == disk][0]
    assert block(ns, x, 'R')[c] == pytest.approx(block(ns, x, 'S')[disk], rel=1e-9)
    assert ns['ENTRY_TASK'][disk] not in ns['GT_OF']


def test_declared_lld_scales_the_station_rate():
    # A processor declared twice as fast serves at twice the rate, so every
    # residence time on it halves. The multiplier reaches the emitted script as
    # source, next to the softmin it multiplies.
    src = export_nlp(idle_multiserver(), lld={'PS': '2.0'})
    assert 'lambda n: softmin(n, 4.0) * (2.0)' in src

    base_ns, base_x = run_emitted(idle_multiserver())
    ns, x = run_emitted(idle_multiserver(), lld={'PS': '2.0'})
    srv = ns['ACTS'].index('as')
    assert ns['objective'](x) < 1e-6
    assert block(ns, x, 'W')[srv] == pytest.approx(
        0.5 * block(base_ns, base_x, 'W')[srv], rel=1e-6)


def test_refuses_an_lld_that_is_not_smooth():
    # The Jacobian differentiates a rate by complex step, so a kink is not a
    # slow path, it is a wrong derivative. Every rejection names its station.
    # a kink that refuses a complex argument outright
    with pytest.raises(LqnNlpExportError) as err:
        export_nlp(idle_multiserver(), lld={'PS': 'max(n, 2.0)'})
    assert 'PS' in str(err.value) and 'complex step' in str(err.value)

    # and the dangerous kind: one that accepts it and quietly returns a real,
    # leaving the complex step at zero where the true slope is one
    with pytest.raises(LqnNlpExportError) as err:
        export_nlp(idle_multiserver(), lld={'PS': '1.0 + np.abs(n - 2.0)'})
    assert 'PS' in str(err.value) and 'not analytic' in str(err.value)

    with pytest.raises(LqnNlpExportError) as err:
        export_nlp(idle_multiserver(), lld={'PS': [1.0, 2.0, 3.0, 4.0]})
    assert 'tabulated' in str(err.value)

    with pytest.raises(LqnNlpExportError) as err:
        export_nlp(idle_multiserver(), lld={'nowhere': '1.0'})
    assert 'no processor or task by that name' in str(err.value)


def test_refuses_asynchronous_call_by_name():
    m = LayeredNetwork('async')
    P = Processor(m, 'P', 1, SchedStrategy.PS)
    T1 = Task(m, 'T1', 5, SchedStrategy.REF).on(P).setThinkTime(Exp(1.0))
    T2 = Task(m, 'T2', 2, SchedStrategy.FCFS).on(P).setThinkTime(Exp(1e6))
    e1 = Entry(m, 'e1').on(T1)
    e2 = Entry(m, 'e2').on(T2)
    Activity(m, 'x', Exp(10)).on(T1).boundTo(e1).asynchCall(e2, 1)
    Activity(m, 'y', Exp(10)).on(T2).boundTo(e2).repliesTo(e2)
    with pytest.raises(LqnNlpExportError) as err:
        export_nlp(m)
    assert 'asynchronous call' in str(err.value)
    assert 'x=>e2' in str(err.value)


def test_refuses_model_with_no_reference_task():
    m = LayeredNetwork('noref')
    P = Processor(m, 'P', 1, SchedStrategy.PS)
    T = Task(m, 'T', 2, SchedStrategy.FCFS).on(P).setThinkTime(Exp(1e6))
    e = Entry(m, 'e').on(T)
    Activity(m, 'a', Exp(10)).on(T).boundTo(e).repliesTo(e)
    with pytest.raises(LqnNlpExportError) as err:
        export_nlp(m)
    assert 'no reference task' in str(err.value)


def test_writes_the_script_to_a_file(tmp_path):
    out = tmp_path / 'cs_nlp.py'
    text = export_nlp(client_server(), str(out))
    assert out.read_text() == text
    assert text.lstrip().startswith('"""')


# ----------------------------------------------------------------------
# AND-fork blocks: concurrency, and the closure that prices it
# ----------------------------------------------------------------------

def test_refuses_an_and_fork_by_default():
    # Concurrency is not a reading of the model the way a loop count is: it
    # needs an assumption about branch SHAPE, so the caller has to ask for it.
    with pytest.raises(LqnNlpExportError) as err:
        export_nlp(forked())
    assert 'AND-fork or AND-join at activity b1' in str(err.value)
    assert "fork='exp'" in str(err.value)


def test_a_model_without_a_fork_carries_no_fork_machinery():
    # The split into elapsed time and thread-seconds exists only where the two
    # differ, so a sequential model must be untouched by this feature.
    ns = {'__name__': 'emitted_nlp'}
    exec(compile(export_nlp(client_server()), '<emitted_nlp>', 'exec'), ns)
    assert ns['NFORK'] == 0
    assert dict(ns['_BLOCKS'])['SA'] == 0
    assert not ns['MNL_K']


def test_fork_program_stays_square():
    # The forking entry gives up its linear equality and gains a residual, so
    # the counting that makes a zero objective an isolated root is preserved.
    ns = {'__name__': 'emitted_nlp'}
    exec(compile(export_nlp(forked(), fork='exp'), '<emitted_nlp>', 'exec'), ns)
    assert ns['NFORK'] == 1
    assert len(ns['MNL_K']) == 1
    rank = np.linalg.matrix_rank(ns['AEQ'])
    assert rank == ns['AEQ'].shape[0]
    assert ns['NVAR'] - rank == ns['NRES']


def test_closure_is_the_exponential_order_statistic():
    # E[max] of k iid exponentials of mean D is D times the k-th harmonic
    # number, and for unequal means it is the inclusion-exclusion sum.
    ns = {'__name__': 'emitted_nlp'}
    exec(compile(export_nlp(forked(), fork='exp'), '<emitted_nlp>', 'exec'), ns)
    fj = ns['fj_max_exp']
    assert fj(np.array([1.0, 1.0])) == pytest.approx(1.5, rel=1e-12)
    assert fj(np.array([1.0] * 3)) == pytest.approx(1 + 1 / 2.0 + 1 / 3.0, rel=1e-12)
    assert fj(np.array([1.0] * 4)) == pytest.approx(1 + 1 / 2.0 + 1 / 3.0 + 0.25,
                                                    rel=1e-12)
    assert fj(np.array([2.0, 1.0])) == pytest.approx(2 + 1 - 1 / (0.5 + 1.0), rel=1e-12)
    # a branch of zero mean is absorbed rather than dividing by zero
    assert fj(np.array([1.0, 0.0])) == pytest.approx(1.0, rel=1e-12)
    # and it is between the max of the means and their sum, never outside
    m = np.array([0.3, 1.7, 0.9])
    assert m.max() < fj(m) < m.sum()


def test_closure_gradient_is_exact():
    # The Newton step differentiates the closure by complex step over a VECTOR
    # argument; a wrong gradient there costs the root-find its convergence and
    # shows up as an objective that stalls rather than as an error.
    ns = {'__name__': 'emitted_nlp'}
    exec(compile(export_nlp(forked(), fork='exp'), '<emitted_nlp>', 'exec'), ns)
    fj, u, h = ns['fj_max_exp'], np.array([0.7, 1.9, 0.3]), 1e-6
    g = ns['mgrad'](fj, u)
    fd = np.array([(fj(u + h * np.eye(3)[j]) - fj(u - h * np.eye(3)[j])) / (2 * h)
                   for j in range(3)])
    assert np.allclose(g, fd, rtol=1e-6, atol=1e-9)


def test_fork_jacobian_matches_a_numerical_one():
    ns = {'__name__': 'emitted_nlp'}
    exec(compile(export_nlp(forked(), fork='exp'), '<emitted_nlp>', 'exec'), ns)
    x = ns['initial_point'](5)
    J, h = ns['jacobian'](x), 1e-7
    Jn = np.zeros_like(J)
    for j in range(ns['NVAR']):
        e = np.zeros(ns['NVAR'])
        e[j] = h
        Jn[:, j] = (ns['residuals'](x + e) - ns['residuals'](x - e)) / (2 * h)
    assert np.max(np.abs(J - Jn)) < 1e-6 * max(1.0, float(np.max(np.abs(J))))


def test_idle_fork_entry_is_held_for_the_max_not_the_sum():
    # On a delay processor W(a) = D(a) exactly, so the entry's service time is
    # the sequential part plus E[max] of the branches, with nothing else in it.
    d1, d2 = 0.30, 0.15
    ns, x = run_emitted(forked(d1, d2), fork='exp')
    assert ns['objective'](x) < 1e-12
    S = block(ns, x, 'S')
    es = ns['ENTRIES'].index('ES')
    concurrent = ns['fj_max_exp'](np.array([d1, d2]))
    assert S[es] == pytest.approx(0.01 + 0.02 + concurrent, rel=1e-9)
    assert S[es] < 0.01 + 0.02 + d1 + d2          # strictly faster than in series
    assert concurrent == pytest.approx(d1 + d2 - 1.0 / (1 / d1 + 1 / d2), rel=1e-12)


def test_fork_entry_occupies_its_pool_additively():
    # The elapsed hold collapses to the max, but the entry runs both branches,
    # so the thread-seconds it takes out of the pool stay the sum. B(c) must be
    # built from the latter or the callee looks less loaded than it is.
    d1, d2 = 0.30, 0.15
    ns, x = run_emitted(forked(d1, d2), fork='exp')
    S, SA, Xc, BC = (block(ns, x, k) for k in ('S', 'SA', 'Xc', 'BC'))
    es = ns['ENTRIES'].index('ES')
    assert SA[es] == pytest.approx(0.01 + 0.02 + d1 + d2, rel=1e-9)
    assert S[es] < SA[es]
    assert BC[0] == pytest.approx(Xc[0] * SA[es], rel=1e-6)


def test_fork_raises_throughput_over_the_same_work_in_series():
    # The whole point of the change: concurrency is faster. Same demands, same
    # population, branches run together instead of one after the other.
    ns_f, xf = run_emitted(forked(0.30, 0.15, nclient=3, host=2), fork='exp')
    m = forked(0.30, 0.15, nclient=3, host=2)
    # the same model read sequentially: strip the fork, chain the activities
    seq = LayeredNetwork('seq')
    PC = Processor(seq, 'PC', 100, SchedStrategy.PS)
    PS = Processor(seq, 'PS', 2, SchedStrategy.PS)
    C = Task(seq, 'Client', 3, SchedStrategy.REF).on(PC).setThinkTime(Exp(1.0))
    S = Task(seq, 'Server', 8, SchedStrategy.FCFS).on(PS).setThinkTime(Exp(1e6))
    EC = Entry(seq, 'EC').on(C)
    ES = Entry(seq, 'ES').on(S)
    Activity(seq, 'ac', Exp(1000.0)).on(C).boundTo(EC).synchCall(ES, 1)
    h = Activity(seq, 'head', Exp(100.0)).on(S).boundTo(ES)
    b1 = Activity(seq, 'b1', Exp(1 / 0.30)).on(S)
    b2 = Activity(seq, 'b2', Exp(1 / 0.15)).on(S)
    t = Activity(seq, 'tail', Exp(50.0)).on(S).repliesTo(ES)
    S.add_precedence(ActivityPrecedence.Serial(h, b1, b2, t))
    ns_s, xs = run_emitted(seq)
    assert ns_f['objective'](xf) < 1e-8 and ns_s['objective'](xs) < 1e-8
    ec = ns_f['ENTRIES'].index('EC')
    assert block(ns_f, xf, 'Xe')[ec] > block(ns_s, xs, 'Xe')[ec]


def test_refuses_a_quorum_join_by_name():
    # A join that fires on k of n branches is the k-th order statistic, which
    # is a different closure; carrying it under this one would be wrong.
    with pytest.raises(LqnNlpExportError) as err:
        export_nlp(forked(quorum=1), fork='exp')
    assert 'quorum join at tail' in str(err.value)
    assert '1 of 2' in str(err.value)


def test_refuses_branches_that_never_rejoin():
    with pytest.raises(LqnNlpExportError) as err:
        export_nlp(forked(join=False), fork='exp')
    assert 'never reconverge' in str(err.value)


def test_refuses_a_branch_with_no_work():
    # A branch of identically zero completion time makes the order statistic
    # degenerate, and its mean is a 0/0 in the closure.
    with pytest.raises(LqnNlpExportError) as err:
        export_nlp(forked(d1=0.0), fork='exp')
    assert 'identically zero' in str(err.value)


def test_refuses_a_fork_nested_in_a_branch():
    m = LayeredNetwork('nested')
    PC = Processor(m, 'PC', 100, SchedStrategy.PS)
    PS = Processor(m, 'PS', 2, SchedStrategy.PS)
    C = Task(m, 'Client', 2, SchedStrategy.REF).on(PC).setThinkTime(Exp(1.0))
    S = Task(m, 'Server', 8, SchedStrategy.FCFS).on(PS).setThinkTime(Exp(1e6))
    EC = Entry(m, 'EC').on(C)
    ES = Entry(m, 'ES').on(S)
    Activity(m, 'ac', Exp(1000.0)).on(C).boundTo(EC).synchCall(ES, 1)
    head = Activity(m, 'head', Exp(100.0)).on(S).boundTo(ES)
    b1 = Activity(m, 'b1', Exp(10.0)).on(S)
    b2 = Activity(m, 'b2', Exp(10.0)).on(S)
    c1 = Activity(m, 'c1', Exp(10.0)).on(S)
    c2 = Activity(m, 'c2', Exp(10.0)).on(S)
    mid = Activity(m, 'mid', Exp(100.0)).on(S)
    tail = Activity(m, 'tail', Exp(100.0)).on(S).repliesTo(ES)
    S.add_precedence(ActivityPrecedence.AndFork(head, [b1, b2]))
    S.add_precedence(ActivityPrecedence.AndFork(b1, [c1, c2]))
    S.add_precedence(ActivityPrecedence.AndJoin([c1, c2], mid))
    S.add_precedence(ActivityPrecedence.AndJoin([mid, b2], tail))
    with pytest.raises(LqnNlpExportError) as err:
        export_nlp(m, fork='exp')
    assert 'nested' in str(err.value) or 'one join target' in str(err.value)


def test_fork_mode_is_validated():
    with pytest.raises(LqnNlpExportError) as err:
        export_nlp(client_server(), fork='max')
    assert "'refuse' or 'exp'" in str(err.value)


def test_fork_header_states_the_assumption():
    src = export_nlp(forked(), fork='exp')
    assert 'THIS MODEL FORKS' in src
    assert 'E[MAX] IS NOT A FUNCTION OF THE BRANCH MEANS' in src
    assert 'FORK_BRANCHES = ' in src
    assert 'THIS MODEL FORKS' not in export_nlp(client_server())
