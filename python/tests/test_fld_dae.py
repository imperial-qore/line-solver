"""The fluid 'dae' method: the min-normal closure solved as one system.

WHAT IS ACTUALLY BEING TESTED. 'dae' is not a new closure -- it is the SAME
closure 'minnormal' computes, taken from the same drift, the same rate factors
and the same Lyapunov equation. So a test that only checked "does it give a
plausible number" would pass on an implementation that quietly reproduced
'minnormal' and threw the DAE away. The properties below are the ones that
distinguish the two:

  * population conservation is an EQUATION here, not a consequence of the drift.
    'minnormal' conserves to integrator tolerance; this conserves to machine
    precision, at the fixed point AND at every point of the trajectory.
  * the transient carries a TIME-VARYING covariance. 'minnormal' evaluates its
    whole transient at the single stationary variance.
  * a finite capacity region is solvable. Every other fluid method refuses one
    outright, because an ODE has nowhere to put a linear inequality.

and the refusals, which have to name the model feature rather than surface from
inside the Newton solve.
"""

import numpy as np
import pytest

from line_solver import (Network, Delay, Queue, ClosedClass, OpenClass, Source,
                         Sink, Exp, SchedStrategy, DropStrategy, SolverFLD)
from line_solver.solvers.solver_fld.methods.dae import DaeSolver
from line_solver.solvers.solver_fld.methods.minnormal import MinNormalSolver
from line_solver.solvers.solver_fld.options import SolverFLDOptions


def _cqn(N=5, mu=2.0, nservers=1):
    model = Network('cqn')
    d = Delay(model, 'Think')
    q = Queue(model, 'Q1', SchedStrategy.PS)
    if nservers > 1:
        q.setNumberOfServers(nservers)
    c = ClosedClass(model, 'C', N, d, 0)
    d.setService(c, Exp(1.0))
    q.setService(c, Exp(mu))
    model.link(Network.serialRouting(d, q))
    return model


def _opts(**kw):
    o = SolverFLDOptions()
    o.method = 'dae'
    o.tol = 1e-8
    for k, v in kw.items():
        setattr(o, k, v)
    return o


def _solve(model, **kw):
    return DaeSolver(model.getStruct(), _opts(**kw)).solve()


# --------------------------------------------------------------------- closure

def test_dae_reproduces_the_minnormal_closure():
    """Same closure, different discharge. They must agree to the accuracy the
    substitution route stops at -- its outer loop stops at a coarse tolerance
    while the Newton runs to options.tol, so the DAE answer is the tighter of
    the two and the gap is the substitution error, not a disagreement."""
    model = _cqn()
    dae = _solve(model)
    o = SolverFLDOptions()
    o.method = 'minnormal'
    mn = MinNormalSolver(_cqn().getStruct(), o).solve()
    np.testing.assert_allclose(dae.QN, mn.QN, rtol=1e-3, atol=1e-4)
    np.testing.assert_allclose(dae.XN, mn.XN, rtol=1e-3, atol=1e-4)


def test_the_newton_actually_converges():
    r = _solve(_cqn())
    assert r.moments['converged'], "the simultaneous solve did not reach options.tol"
    assert r.moments['residual'] < 1e-8


@pytest.mark.parametrize('N', [2, 5, 20])
def test_conservation_is_an_equation_not_a_consequence(N):
    """The whole point of the algebraic row. A drift whose rows sum to zero
    conserves the population only as well as the integrator does; writing it as
    an equation conserves it to the solve's own tolerance, which is machine
    precision here."""
    r = _solve(_cqn(N=N))
    assert abs(float(r.QN.sum()) - N) < 1e-9
    assert r.moments['conservation'] < 1e-9


def test_saturation_is_reached_without_a_kink_probe():
    """At N=20 the queue is saturated and the first-order fixed point sits
    exactly on the sigma2=0 kink of min(n,c) -- the configuration 'minnormal'
    needs its two-sided Jacobian probe for. The DAE seeds the variance POSITIVE
    and never adopts that iterate, so it must simply solve."""
    r = _solve(_cqn(N=20))
    assert r.moments['converged']
    # think time 1 at throughput 2 holds 2 jobs; the rest queue
    np.testing.assert_allclose(r.QN.ravel(), [2.0, 18.0], rtol=1e-6, atol=1e-6)


# ------------------------------------------------------------------- transient

def test_transient_carries_a_time_varying_covariance():
    """The property that separates this from 'minnormal', which evaluates its
    whole transient at the single stationary variance."""
    r = _solve(_cqn(), timespan=(0.0, 8.0))
    S = r.moments['Sigmat']
    assert S is not None, "no transient covariance was produced"
    assert S.shape[0] == S.shape[1]
    assert S.shape[2] == len(r.t)
    # Sigma(0) = 0 is the consistent AND the physically right initialisation:
    # the population at t=0 is a known deterministic state.
    assert np.allclose(S[:, :, 0], 0.0)
    # and it must actually move
    assert np.max(np.abs(S[:, :, -1])) > 1e-3
    qv = r.moments['QVart'][(1, 0)]
    assert qv[0] == 0.0
    assert qv[-1] > 0.0


def test_conservation_holds_along_the_whole_trajectory():
    """The algebraic row is enforced at every reported point, not just at the
    fixed point. An ODE integration of the same drift would drift off it."""
    r = _solve(_cqn(N=5), timespan=(0.0, 8.0))
    tot = np.zeros(len(r.t))
    for i in range(2):
        tot += np.asarray(r.QNt[(i, 0)])
    assert np.max(np.abs(tot - 5.0)) < 1e-8


def test_transient_settles_on_the_steady_state():
    """A long enough horizon must reproduce the table, or the trajectory and the
    fixed point are being read off two different drifts."""
    ss = _solve(_cqn())
    tr = _solve(_cqn(), timespan=(0.0, 60.0))
    last = np.array([tr.QNt[(i, 0)][-1] for i in range(2)])
    np.testing.assert_allclose(last, ss.QN.ravel(), rtol=1e-4, atol=1e-4)


# ----------------------------------------------------- finite capacity regions

def _cqn_region(cap, N=8, rule=None, mu=2.0):
    model = Network('fcr')
    d = Delay(model, 'Think')
    q = Queue(model, 'Q1', SchedStrategy.PS)
    c = ClosedClass(model, 'C', N, d, 0)
    d.setService(c, Exp(1.0))
    q.setService(c, Exp(mu))
    model.link(Network.serialRouting(d, q))
    reg = model.addRegion([q])
    reg.setGlobalMaxJobs(cap)
    reg.setDropRule(c, DropStrategy.WAITQ if rule is None else rule)
    return model


@pytest.mark.parametrize('cap', [3, 2])
def test_a_capacity_region_binds_and_conserves(cap):
    """SolverFLD refuses a finite capacity region for every other method. Here
    the cap is an algebraic equation, so the region sits exactly ON it and the
    mass that does not fit is held in the waiting queue rather than vanishing."""
    r = _solve(_cqn_region(cap))
    capinfo = r.moments['capacity']
    assert capinfo['active'].size == 1, "the cap should bind"
    # the region holds exactly its cap
    assert abs(float(r.QN[1, 0]) - cap) < 1e-6
    # and the population is still all there: at a station, or blocked
    assert abs(float(r.QN.sum()) + float(capinfo['blocked']) - 8.0) < 1e-6
    # a waiting queue drains at a finite positive rate
    assert float(np.min(capinfo['drain'])) > 0.0


def test_a_cap_the_population_cannot_reach_never_binds():
    """A region capped above the whole closed population is not a constraint,
    and must be pruned rather than left in the active-set loop being tested on
    every pass."""
    r = _solve(_cqn_region(50))
    assert r.moments['capacity']['active'].size == 0
    assert float(r.moments['capacity']['blocked']) == 0.0
    # and the answer is the unconstrained one
    free = _solve(_cqn(N=8))
    np.testing.assert_allclose(r.QN, free.QN, rtol=1e-5, atol=1e-6)


def test_the_region_transient_holds_the_cap_and_locates_the_crossing():
    """The transient under a cap is a HYBRID DAE: the system switches every time
    the region fills or drains. Integrating with the binding set frozen would
    report the unconstrained path through a cap the model declares, so what is
    tested is that the path NEVER exceeds the cap, that the crossing is located
    rather than stepped over, and that the trajectory ends where the steady-state
    solve says it should."""
    ss = _solve(_cqn_region(5))
    tr = _solve(_cqn_region(5), timespan=(0.0, 40.0))
    region = np.asarray(tr.QNt[(1, 0)])
    assert float(np.max(region)) <= 5.0 + 1e-6, "the trajectory left the feasible set"
    switches = tr.moments['capacity']['switches']
    assert any(kind == 'activate' for _t, _c, kind in switches), \
        "the cap was reached and must have been activated"
    # and the two formulations describe one model: the drain RATE the steady state
    # solves for and the admitted FLOW the transient carries agree at the fixed point
    np.testing.assert_allclose(region[-1], float(ss.QN[1, 0]), rtol=1e-4, atol=1e-4)


def test_a_region_cap_is_released_when_the_waiting_room_empties():
    """A cap that binds at t=0 and stops binding must be RELEASED, or the
    trajectory would be held on a cap the drift has already left."""
    ss = _solve(_cqn_region(2, mu=8.0))
    tr = _solve(_cqn_region(2, mu=8.0), timespan=(0.0, 40.0))
    switches = tr.moments['capacity']['switches']
    assert any(kind == 'release' for _t, _c, kind in switches)
    region = np.asarray(tr.QNt[(1, 0)])
    assert float(np.max(region)) <= 2.0 + 1e-6
    np.testing.assert_allclose(region[-1], float(ss.QN[1, 0]), rtol=1e-4, atol=1e-4)


# ------------------------------------------------------------ station buffers

def _capped(cap, N=8, mu=2.0, servers=1, cls_caps=None):
    model = Network('capped')
    d = Delay(model, 'Think')
    q = Queue(model, 'Q1', SchedStrategy.FCFS if cls_caps is None else SchedStrategy.PS)
    if servers > 1:
        q.setNumberOfServers(servers)
    if cap is not None:
        q.setCapacity(cap)
    if cls_caps is None:
        c = ClosedClass(model, 'C', N, d, 0)
        d.setService(c, Exp(1.0))
        q.setService(c, Exp(mu))
    else:
        a = ClosedClass(model, 'A', N, d, 0)
        b = ClosedClass(model, 'B', N, d, 0)
        d.setService(a, Exp(1.0)); d.setService(b, Exp(1.0))
        q.setService(a, Exp(mu)); q.setService(b, Exp(mu))
        q.setClassCapacity(a, cls_caps[0]); q.setClassCapacity(b, cls_caps[1])
    model.link(Network.serialRouting(d, q))
    return model


def test_a_station_buffer_binds_and_holds_the_blocked_job_upstream():
    """A station buffer is the one-station case of the same row -- but NOT of the
    same model. LINE refuses to lose a closed job (State.arrivalIsLost) and
    disables the upstream departure instead, so the blocked mass is still AT the
    upstream station and still counted there: the station queues sum to the whole
    population and nothing is staged. That is the opposite of a region, whose
    blocked jobs are reported separately."""
    r = _solve(_capped(2, N=8, mu=2.0))
    capinfo = r.moments['capacity']
    assert capinfo['active'].size == 1
    assert abs(float(r.QN[1, 0]) - 2.0) < 1e-6, "the buffer holds exactly its cap"
    assert abs(float(r.QN.sum()) - 8.0) < 1e-6, "a held job is at the upstream station"
    assert float(capinfo['blocked']) == 0.0, "a station buffer has no waiting room"
    # the multiplier is the fraction of upstream completions the cap admits, so it
    # is a fraction and not a rate: below one where the cap binds
    assert 0.0 < float(capinfo['multiplier'][0]) < 1.0


def test_a_station_buffer_the_population_cannot_reach_is_inert():
    """refreshCapacity derives a finite classcap for every closed model, so a cap
    that cannot bind is the common case, not a corner one: it must be pruned and
    the answer must be the unconstrained one."""
    r = _solve(_capped(50, N=8))
    assert r.moments['capacity']['active'].size == 0
    free = _solve(_capped(None, N=8))
    np.testing.assert_allclose(r.QN, free.QN, rtol=1e-6, atol=1e-8)


def test_the_station_total_implied_by_its_class_buffers_is_pruned():
    """LINE derives the station total from the per-class buffers, so a two-class
    station capped 3 and 3 also declares a total of 6 -- exactly the sum of the two
    class rows. All three would bind together with rank 2, and the multipliers
    would be one arbitrary point of a line of solutions."""
    r = _solve(_capped(None, N=6, mu=2.0, cls_caps=(3, 3)))
    labels = r.moments['capacity']['label']
    assert len(labels) == 2, labels
    assert all('class' in lb for lb in labels), labels


def test_two_class_buffers_bind_at_once():
    """TWO CAPS BINDING TOGETHER used to be refused outright: one region carried a
    single throttle, so a second equality had no control to satisfy it. Each active
    row now carries its own multiplier."""
    r = _solve(_capped(None, N=6, mu=2.0, cls_caps=(2, 3)))
    capinfo = r.moments['capacity']
    assert capinfo['active'].size == 2, capinfo['label']
    np.testing.assert_allclose(capinfo['value'], capinfo['b'], rtol=0, atol=1e-6)
    assert r.moments['converged']


def test_an_arrival_is_lost_at_a_full_buffer_of_an_open_class():
    """The same predicate that holds a closed job LOSES an open one: the external
    stream is memoryless, so a job that finds the buffer full never enters. The
    arrival fires and only the carried flow is admitted, so throughput is the
    service rate exactly and the multiplier is the admitted fraction."""
    model = Network('loss')
    s = Source(model, 'S')
    q = Queue(model, 'Q1', SchedStrategy.FCFS)
    k = Sink(model, 'K')
    q.setCapacity(5)
    c = OpenClass(model, 'C')
    s.setArrival(c, Exp(4.0))
    q.setService(c, Exp(2.0))
    model.link(Network.serialRouting(s, q, k))
    r = _solve(model)
    capinfo = r.moments['capacity']
    assert capinfo['active'].size == 1
    assert abs(float(r.QN[1, 0]) - 5.0) < 1e-6
    # in overload the fluid loss rate is exact: what gets through is the server
    np.testing.assert_allclose(float(r.TN[1, 0]), 2.0, rtol=1e-6, atol=1e-6)
    np.testing.assert_allclose(float(capinfo['multiplier'][0]), 0.5, rtol=1e-4, atol=1e-4)


def test_a_station_buffer_transient_never_leaves_the_feasible_set():
    ss = _solve(_capped(2, N=8, mu=2.0))
    tr = _solve(_capped(2, N=8, mu=2.0), timespan=(0.0, 40.0))
    buf = np.asarray(tr.QNt[(1, 0)])
    assert float(np.max(buf)) <= 2.0 + 1e-6
    np.testing.assert_allclose(buf[-1], float(ss.QN[1, 0]), rtol=1e-4, atol=1e-4)


def test_a_blocking_drop_rule_at_a_station_is_refused_by_name():
    """BAS/BBS/RSRD give the upstream station a blocked-server state and the
    retrial rules add an orbit. Neither is a constraint on this drift."""
    model = _capped(2, N=8)
    q = model.getNodes()[1]
    q.setDropRule(model.getClasses()[0], DropStrategy.BAS)
    with pytest.raises(ValueError, match='BAS'):
        _solve(model)


def test_default_resolves_to_dae_on_a_blocked_model():
    """The refusal below was correct and its advice was useless: it told the
    caller to type the one method the resolution could have picked itself.
    'default' now stands for 'dae' wherever a buffer or a region BINDS and the
    dae route accepts the model, and for nothing else -- an unconstrained twin
    must resolve exactly as it did before."""
    capped = _capped(2, N=4, mu=2.0)
    s = SolverFLD(capped)
    assert s.resolveMethod(s.options) == 'dae'
    assert abs(float(np.asarray(s.getAvgQLen()).ravel()[1]) - 2.0) < 1e-6

    region = _cqn_region(2)
    sr = SolverFLD(region)
    assert sr.resolveMethod(sr.options) == 'dae'

    free = _capped(99, N=4, mu=2.0)  # a cap the population cannot reach
    sf = SolverFLD(free)
    assert sf.resolveMethod(sf.options) != 'dae'


def test_every_other_fluid_method_refuses_a_binding_buffer():
    """Nothing in the fluid tree reads sn.cap, so a capped station was integrated
    as an unbounded one and the table reported more jobs in the buffer than the
    buffer holds. Only the route that can enforce the cap may accept the model."""
    model = _capped(2, N=4, mu=2.0)
    # runAnalyzer rather than a getter: the getters re-raise with their own message
    with pytest.raises(Exception, match="dae"):
        SolverFLD(model, method='minnormal').runAnalyzer()
    q = SolverFLD(model, method='dae').getAvgQLen()
    assert abs(float(np.asarray(q).ravel()[1]) - 2.0) < 1e-6


# ------------------------------------------------------- overlapping limits

def test_two_overlapping_regions_bind_at_once():
    """Two regions sharing a station: an admission into one is an INTERNAL move of
    the other, so the outer cap has to count the mass waiting in the inner one's
    room. Without that the outer cap is met on paper while the region holds more,
    and the Newton system is inconsistent rather than merely inexact."""
    model = Network('two')
    d = Delay(model, 'Think')
    q1 = Queue(model, 'Q1', SchedStrategy.FCFS)
    q2 = Queue(model, 'Q2', SchedStrategy.FCFS)
    q3 = Queue(model, 'Q3', SchedStrategy.FCFS)
    c = ClosedClass(model, 'C', 14, d, 0)
    d.setService(c, Exp(1.0))
    for q in (q1, q2, q3):
        q.setService(c, Exp(4.0))
    model.addRegion([q1, q2]).setGlobalMaxJobs(5)
    model.addRegion([q2, q3]).setGlobalMaxJobs(4)
    model.link(Network.serialRouting(d, q1, q2, q3))
    r = _solve(model)
    capinfo = r.moments['capacity']
    assert capinfo['active'].size == 2
    np.testing.assert_allclose(capinfo['value'], capinfo['b'], rtol=0, atol=1e-6)
    assert r.moments['conservation'] < 1e-8
    assert abs(float(r.QN.sum()) + float(capinfo['blocked']) - 14.0) < 1e-6


def test_a_region_global_and_class_cap_bind_at_once():
    """The case the single per-region throttle had to refuse by name: a global cap
    and a per-class cap of the SAME region, binding together. The rooms are per
    class, so the global row gates both and the class row one, and a room gated
    twice drains at the harmonic composition of the two rates."""
    model = Network('mcreg')
    d = Delay(model, 'Think')
    q = Queue(model, 'Q1', SchedStrategy.PS)
    a = ClosedClass(model, 'A', 8, d, 0)
    b = ClosedClass(model, 'B', 8, d, 0)
    d.setService(a, Exp(1.0)); d.setService(b, Exp(1.0))
    q.setService(a, Exp(3.0)); q.setService(b, Exp(3.0))
    reg = model.addRegion([q])
    reg.setGlobalMaxJobs(6)
    reg.setClassMaxJobs(a, 2)
    model.link(Network.serialRouting(d, q))
    r = _solve(model)
    capinfo = r.moments['capacity']
    assert capinfo['active'].size == 2, capinfo['label']
    np.testing.assert_allclose(capinfo['value'], capinfo['b'], rtol=0, atol=1e-6)
    assert r.moments['conservation'] < 1e-8


def test_a_drop_rule_that_is_not_a_waiting_queue_is_refused_by_name():
    """Only a waiting queue conserves the population and throttles the flow.
    DROP destroys the job, which is a different event set, not a constraint on
    this drift."""
    with pytest.raises(ValueError, match='waiting queue'):
        _solve(_cqn_region(3, rule=DropStrategy.DROP))


# ------------------------------------------------------------------- refusals

def test_the_state_cap_is_refused_by_name():
    with pytest.raises(ValueError, match='dae_maxstate'):
        DaeSolver(_cqn().getStruct(), _opts(config={'dae_maxstate': 1})).solve()


def test_dps_is_refused_by_name():
    """DPS closes on the covariance BETWEEN a station's class coordinates, not
    on the station total, so its closure state is a matrix block rather than the
    scalar this solves for."""
    model = Network('dps')
    d = Delay(model, 'Think')
    q = Queue(model, 'Q1', SchedStrategy.DPS)
    c = ClosedClass(model, 'C', 4, d, 0)
    d.setService(c, Exp(1.0))
    q.setService(c, Exp(2.0))
    model.link(Network.serialRouting(d, q))
    with pytest.raises(ValueError, match='per-station variance only'):
        _solve(model)


def test_dae_is_reachable_through_the_solver_front_door():
    """The dispatch, not the solver: options.method='dae' must resolve, and
    'dae' must appear in listValidMethods or runAnalyzerChecks rejects it."""
    assert 'dae' in SolverFLD.listValidMethods()
    solver = SolverFLD(_cqn(), 'dae')
    QN = np.asarray(solver.getAvgQLen())
    assert abs(float(QN.sum()) - 5.0) < 1e-6


def test_an_open_chain_contributes_no_conservation_row():
    """An open chain has no conserved population, so it must not produce an
    algebraic row -- one would pin the state to a population the model never
    fixed."""
    model = Network('open')
    src = Source(model, 'Source')
    q = Queue(model, 'Q1', SchedStrategy.PS)
    sink = Sink(model, 'Sink')
    c = OpenClass(model, 'C')
    src.setArrival(c, Exp(1.0))
    q.setService(c, Exp(3.0))
    model.link(Network.serialRouting(src, q, sink))
    r = _solve(model)
    # M/M/1 at rho = 1/3: the mean queue length is rho/(1-rho) = 0.5
    assert r.moments['converged']
    assert float(r.QN[1, 0]) > 0.0


# ---------------------------------------------------------------------------
# The non-hyperbolic fallback ladder: minnormal -> dae -> first order
# ---------------------------------------------------------------------------
#
# fluid_lyapunov raises when the drift Jacobian at the converged mean is not
# exponentially stable on range(D). The dominant case is NEUTRAL rather than
# unstable, and it is an artifact of the alternation: MinNormalSolver must start
# at sigma2 = 0, where min(n,c) has no derivative, so a saturated or balanced
# model's first-order fixed point lands on the kink and sits on a continuum of
# equilibria. DaeSolver seeds the variance positive and never adopts sigma2 = 0,
# so the same closure has an isolated, hyperbolic fixed point there.


def _balanced_cycle(N, nservers=1, sched=SchedStrategy.PS):
    """Two identical stations in a closed cycle: the fluid drift is degenerate.

    Every state with both populations at or above the server count is a fluid
    equilibrium, so the first-order methods have no reason to prefer one point of
    that continuum over another and land wherever the integrator stopped.
    """
    model = Network('balanced')
    q1 = Queue(model, 'Q1', sched)
    q2 = Queue(model, 'Q2', sched)
    q1.setNumberOfServers(nservers)
    q2.setNumberOfServers(nservers)
    c = ClosedClass(model, 'C', N, q1, 0)
    q1.setService(c, Exp(1.0))
    q2.setService(c, Exp(1.0))
    model.link(Network.serialRouting(q1, q2))
    return model


@pytest.mark.parametrize('N,nservers,sched', [
    (4, 1, SchedStrategy.PS),
    (6, 1, SchedStrategy.PS),
    (6, 2, SchedStrategy.PS),
    (10, 1, SchedStrategy.PS),
    (10, 2, SchedStrategy.PS),
    (6, 1, SchedStrategy.FCFS),
    (6, 2, SchedStrategy.FCFS),
])
def test_a_declined_minnormal_lands_on_dae_and_on_the_exact_answer(N, nservers, sched):
    """The ladder's first rung, and that it is worth taking.

    By symmetry the exact answer splits the population evenly, which is what the
    closure gives once the degeneracy is broken. The first-order fallback this
    replaces returned [9 1] at N=10 -- a point of the continuum, not the mean.
    """
    solver = SolverFLD(_balanced_cycle(N, nservers, sched), 'minnormal')
    solver.runAnalyzer()
    assert solver._fallback_method == 'dae'
    QN = np.asarray(solver.getAvgQLen()).ravel()
    assert QN == pytest.approx([N / 2.0, N / 2.0], abs=1e-6)


def test_a_genuinely_unstable_fixed_point_walks_past_dae_to_first_order():
    """The ladder's last rung.

    An overloaded open station has no stationary distribution at all, so no
    closure has a stationary covariance there and 'dae' declines on the same
    exception 'minnormal' did. The mean is still reported, by the first-order
    method, which is the whole reason the last rung exists.
    """
    model = Network('overloaded')
    src = Source(model, 'Source')
    q = Queue(model, 'Q1', SchedStrategy.FCFS)
    snk = Sink(model, 'Sink')
    c = OpenClass(model, 'C', 0)
    src.setArrival(c, Exp(1.1))
    q.setService(c, Exp(1.0))
    model.link(Network.serialRouting(src, q, snk))

    solver = SolverFLD(model, 'minnormal')
    solver.runAnalyzer()
    assert solver._fallback_method == 'matrix'


def test_the_dae_rung_is_declined_in_advance_for_the_features_it_cannot_take():
    """fluid_dae_applicable is the STATIC difference set between the closures.

    A rung entered only to be refused a moment later would spend a whole seed
    integration to learn what the model already declares, so each of the three
    conditions is decided before the solve.
    """
    from line_solver.solvers.solver_fld.dae_applicable import fluid_dae_applicable

    # DPS closes on the covariance BETWEEN class coordinates, not the station
    # total, so its closure state is a matrix block rather than the scalar the
    # Newton vector carries.
    dps = Network('dps')
    d = Delay(dps, 'Think')
    q = Queue(dps, 'Q1', SchedStrategy.DPS)
    c = ClosedClass(dps, 'C', 4, d, 0)
    d.setService(c, Exp(1.0))
    q.setService(c, Exp(2.0))
    dps.link(Network.serialRouting(d, q))
    ok, reason = fluid_dae_applicable(dps.getStruct(), _opts())
    assert not ok and 'covariance between' in reason

    # The simultaneous solve is quartic where one Lyapunov solve is cubic, so it
    # carries its own cap, lower than moment_maxstate.
    ok, reason = fluid_dae_applicable(_cqn(N=5).getStruct(), _opts(config={'dae_maxstate': 1}))
    assert not ok and 'dae_maxstate' in reason

    # and the model the ladder does take
    ok, reason = fluid_dae_applicable(_balanced_cycle(6).getStruct(), _opts())
    assert ok and reason == ''
