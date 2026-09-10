"""SolverAUTO ranking tests.

Asserts the selection order: LDES for simulation, MVA over NC over MAM
for analytical solvers (inverted on caches), Fluid where a smooth or high-load
answer is wanted, and per-metric capability routing. Mirrors
line-test.git/test/testsMisc/test_auto_solver_selection.m.
"""

import pytest

from line_solver import (Network, Delay, Queue, Source, Sink, Cache, ClosedClass,
                         OpenClass, SchedStrategy, ReplacementStrategy, Exp, Erlang,
                         MMPP2, DiscreteSampler, SolverAUTO)
from line_solver.solvers.solver_ctmc import SolverCTMC
from line_solver import SolverMVA, SolverNC
import numpy as np


def closed_ps(n=5):
    m = Network('closedPS')
    d = Delay(m, 'd')
    q = Queue(m, 'q', SchedStrategy.PS)
    c = ClosedClass(m, 'c', n, d)
    d.setService(c, Exp(1))
    q.setService(c, Exp(0.5))
    m.link(Network.serialRouting(d, q))
    return m


def open_map():
    m = Network('openMAP')
    s = Source(m, 'src')
    q = Queue(m, 'q', SchedStrategy.FCFS)
    k = Sink(m, 'sink')
    c = OpenClass(m, 'c')
    s.setArrival(c, MMPP2(0.2, 0.6, 0.1, 0.2))
    q.setService(c, Exp(1))
    m.link(Network.serialRouting(s, q, k))
    return m


def cache_model():
    m = Network('cacheModel')
    d = Delay(m, 'd')
    # NC covers RR/FIFO/HLRU caches; LRU is MVA-only, so the inverted
    # NC > MVA order is only observable on a strategy both support.
    cache = Cache(m, 'cache', 4, 2, ReplacementStrategy.FIFO)
    c = ClosedClass(m, 'c', 1, d)
    h = ClosedClass(m, 'hit', 0, d)
    mi = ClosedClass(m, 'miss', 0, d)
    for cls in (c, h, mi):
        d.setService(cls, Exp(1))
    cache.setRead(c, DiscreteSampler(np.ones(4) / 4))
    cache.setHitClass(c, h)
    cache.setMissClass(c, mi)
    P = m.initRoutingMatrix()
    P.set(c, c, d, cache, 1.0)
    P.set(h, c, cache, d, 1.0)
    P.set(mi, c, cache, d, 1.0)
    m.link(P)
    return m


@pytest.mark.parametrize("factory,mode,expected", [
    (lambda: closed_ps(5), 'default', 'MVA'),
    (lambda: closed_ps(60), 'default', 'Fluid'),
    (open_map, 'default', 'MAM'),
    (cache_model, 'default', 'NC'),
    (lambda: closed_ps(5), 'exact', 'NC'),
    (lambda: closed_ps(5), 'sim', 'LDES'),
    (lambda: closed_ps(5), 'fast', 'MVA'),
    (lambda: closed_ps(5), 'accurate', 'Fluid'),
])
def test_mode_ranking(factory, mode, expected):
    solver = SolverAUTO(factory(), mode, verbose=False)
    solver.getAvgTable()
    assert solver.getSelectedSolverName() == expected


@pytest.mark.parametrize("method,expected", [
    ('getCdfRespT', 'Fluid'),
    ('getTranAvg', 'Fluid'),
    ('getProb', 'NC'),
    ('sample', 'SSA'),
    ('sampleAggr', 'LDES'),
    ('getAvgLossTable', 'LDES'),
    ('getMomentTable', 'MVA'),
    ('getSensitivityTable', 'Fluid'),
])
def test_metric_routing(method, expected):
    solver = SolverAUTO(closed_ps(5), verbose=False)
    assert solver._choose_solver_for_method(method) == expected


def test_cache_table_routes_to_nc():
    solver = SolverAUTO(cache_model(), verbose=False)
    assert solver._choose_solver_for_method('getAvgCacheTable') == 'NC'


def test_ldes_leads_ssa_in_sim_mode():
    solver = SolverAUTO(closed_ps(5), 'sim', verbose=False)
    assert solver._select_solver_sim('getCdfRespT') == 'LDES'


def test_jmt_is_never_an_auto_candidate():
    # LDES subsumes SolverJMT's feature set, so automatic selection must never
    # reach the external simulator. The explicit 'jmt' method name is unaffected.
    for mode in ('default', 'sim', 'fast', 'accurate'):
        solver = SolverAUTO(closed_ps(5), mode, verbose=False)
        assert 'JMT' not in solver._candidate_solvers


@pytest.mark.parametrize("mode", ['fast', 'accurate'])
def test_autocorrelated_arrivals_reach_fluid(mode):
    # MAP and MMPP2 are in the Fluid feature set in all three codebases: the ODE
    # takes the stationary rate and conserves flow.
    solver = SolverAUTO(open_map(), mode, verbose=False)
    solver.getAvgTable()
    assert solver.getSelectedSolverName() == 'Fluid'


def test_exact_mode_refuses_rather_than_approximating():
    # An unbounded open MAP model has no exact analytical solver: the request
    # must fail rather than silently return an approximation.
    solver = SolverAUTO(open_map(), 'exact', verbose=False)
    assert solver._select_solver_exact('getMomentTable') in (None, 'CTMC')


def intractable_ctmc():
    """8 PS stations, N=400, Erlang-5 service: logNstates ~ 200."""
    m = Network('intractableCTMC')
    st = [Queue(m, 'q%d' % i, SchedStrategy.PS) for i in range(8)]
    c = ClosedClass(m, 'c', 400, st[0])
    for s in st:
        s.setService(c, Erlang.fitMeanAndOrder(1, 5))
    m.link(Network.serialRouting(*st))
    return m


def test_state_space_gate_accepts_small_and_refuses_large():
    ok_small, _, log_small = SolverCTMC.isStateSpaceTractable(closed_ps(3))
    ok_large, _, log_large = SolverCTMC.isStateSpaceTractable(intractable_ctmc())
    assert ok_small and log_small < 10
    assert not ok_large and log_large > 100


def test_auto_skips_ctmc_when_the_chain_cannot_be_built():
    # getTranProb ranks CTMC alone; an untractable chain must not be proposed.
    assert SolverAUTO(closed_ps(3), verbose=False)._choose_solver_for_method('getTranProb') == 'CTMC'
    assert SolverAUTO(intractable_ctmc(), verbose=False)._choose_solver_for_method('getTranProb') is None


def test_forced_ctmc_token_is_not_screened():
    # An explicit request must reach CTMC and fail with its own diagnostic.
    solver = SolverAUTO(intractable_ctmc(), 'ctmc', verbose=False)
    assert solver._choose_solver_for_method('getAvgTable') in (None, 'CTMC')
    assert 'CTMC' in solver._candidate_solvers


def non_product_form():
    """Multiclass FCFS with unequal per-class rates: no product-form solution."""
    m = Network('nonPF')
    d = Delay(m, 'd')
    q = Queue(m, 'q', SchedStrategy.FCFS)
    c1 = ClosedClass(m, 'c1', 2, d)
    c2 = ClosedClass(m, 'c2', 2, d)
    d.setService(c1, Exp(1))
    d.setService(c2, Exp(1))
    q.setService(c1, Exp(2))
    q.setService(c2, Exp(0.5))
    P = m.initRoutingMatrix()
    P.set(c1, c1, d, q, 1.0)
    P.set(c1, c1, q, d, 1.0)
    P.set(c2, c2, d, q, 1.0)
    P.set(c2, c2, q, d, 1.0)
    m.link(P)
    return m


def test_method_gate_rejects_exact_without_product_form():
    m = non_product_form()
    assert not m.hasProductFormSolution()
    assert not SolverMVA(m).supportsModelMethod('exact')[0]
    assert not SolverNC(m).supportsModelMethod('exact')[0]
    # The other NC methods fall back to comom, so they stay admissible.
    assert SolverNC(m).supportsModelMethod('default')[0]


def test_exact_mode_routes_around_a_non_product_form_model():
    assert SolverAUTO(non_product_form(), 'exact', verbose=False)._select_solver_exact('getAvgTable') == 'CTMC'
    assert SolverAUTO(closed_ps(4), 'exact', verbose=False)._select_solver_exact('getAvgTable') == 'NC'


@pytest.mark.parametrize("n,expected", [(1, 'CTMC'), (2, 'CTMC'), (3, 'MVA')])
def test_small_populations_take_an_exact_solver(n, expected):
    # Total population 2*n: at or below 5 an exact solver is preferred, and on
    # this non-product-form model only CTMC qualifies. Above it the ranking
    # falls back to the approximation.
    m = non_product_form()
    for cls in m.getClasses():
        cls.setNumberOfJobs(n)
    assert SolverAUTO(m, verbose=False)._select_solver_heuristic() == expected


def test_small_product_form_population_still_prefers_mva():
    # MVA leads NC leads CTMC when all three are exact for the model.
    assert SolverAUTO(closed_ps(3), verbose=False)._select_solver_heuristic() == 'MVA'


def test_ctmc_accessors_always_resolve_to_ctmc():
    # State space and generator are CTMC-only concepts: AUTO exposes them under
    # the same names and never routes them through the ranked selection.
    m = closed_ps(2)
    auto = SolverAUTO(m, verbose=False)
    ss = auto.getStateSpace()
    ss = ss[0] if isinstance(ss, tuple) else ss
    ref = SolverCTMC(m).getStateSpace()
    ref = ref[0] if isinstance(ref, tuple) else ref
    assert np.array_equal(np.asarray(ss), np.asarray(ref))

    g = auto.getGenerator()
    g = g[0] if isinstance(g, tuple) else g
    dense = np.asarray(g.todense() if hasattr(g, 'todense') else g)
    assert dense.shape[0] == dense.shape[1] == np.shape(ss)[0]
    assert np.allclose(dense.sum(axis=1), 0)

    # sympy is an optional extra ('symbolic'), so the suite venv may not carry
    # it; the accessor assertions above do not need it and stay unconditional.
    pytest.importorskip('sympy')
    assert auto.getSymbolicGenerator() is not None


def bas_blocking():
    """cqn_bas_blocking: two FCFS queues, one closed class of 2, a finite buffer
    on the second and the BAS drop rule on the first."""
    from line_solver import DropStrategy
    m = Network('cqn_bas_blocking')
    q1 = Queue(m, 'Queue1', SchedStrategy.FCFS)
    q2 = Queue(m, 'Queue2', SchedStrategy.FCFS)
    c = ClosedClass(m, 'Class1', 2, q1, 0)
    q1.setService(c, Exp(1.0))
    q2.setService(c, Exp(0.8))
    q2.setCap(1)
    q1.setDropRule(c, DropStrategy.BAS)
    m.link(Network.serialRouting(q1, q2))
    return m


def test_listvalidmethods_is_gated_on_the_model():
    # listValidMethods USED TO BE the method name universe: it built every family,
    # asked each for its own list and kept everything, so a BAS-blocking model
    # was offered every NC and every LQNS method although both refuse it method
    # by method. It now applies the gate chooseSolverRanked already applies
    # before delegating.
    auto = SolverAUTO(bas_blocking(), verbose=False)
    valid = auto.listValidMethods()
    families = set(t.split('.')[0] for t in valid)

    # SolverNC refuses a binding finite buffer on every method, so neither its
    # methods nor its bare family method name may be offered.
    assert not [t for t in valid if t.startswith('nc.')]
    assert 'nc' not in families
    # SolverLQNS analyzes a LayeredNetwork; its per-layer feature set says yes to
    # a flat closed exponential network only because it was never asked that.
    assert not [t for t in valid if t.startswith('lqns.')]
    assert 'ln' not in families and 'env' not in families
    # The fluid drift ignores sn.cap outside 'dae'/'mol'. SolverFLD declares both
    # the bare and the 'fluid.'-qualified spelling of each name, so qualifying its
    # list a second time yields 'fluid.fluid.dae' as well; that doubling predates
    # the gate and is not what this test is about.
    # 'default' rides with them: on a blocked model SolverFLD resolves it to
    # 'dae' (see _blocked_resolves_to_dae), so it names the same run.
    assert set(t.replace('fluid.', '') for t in valid
               if t.startswith('fluid.')) <= {'dae', 'mol', 'default'}
    assert 'fluid.default' in valid and 'fluid.dae' in valid
    # MVA stays: solver_mva_analyzer routes a BAS model to solver_sqd.
    assert 'mva.sqd' in valid
    # And the simulators, which represent the buffer exactly.
    assert 'ctmc' in families and 'ssa' in families and 'ldes' in families


def test_listallmethods_is_the_unnarrowed_token_universe():
    # The name check must gate on THIS list, so that asking for a method a
    # candidate refuses gets the candidate's own reason rather than a flat
    # "unsupported by this solver".
    auto = SolverAUTO(bas_blocking(), verbose=False)
    every = auto.listAllMethods()
    valid = auto.listValidMethods()
    assert set(valid) <= set(every)
    assert len(every) > len(valid)
    assert 'nc.exact' in every and 'nc.exact' not in valid


def test_an_unconstrained_model_keeps_the_wide_list():
    # The narrowing must be model-sensitive, not a blanket trim: drop the cap
    # and the product-form and fluid families come back.
    valid = SolverAUTO(closed_ps(2), verbose=False).listValidMethods()
    assert 'nc' in set(t.split('.')[0] for t in valid)
    assert len([t for t in valid if t.startswith('fluid.')]) > 2


def test_ba_feature_set_accepts_the_models_the_bounds_are_derived_for():
    # SolverBA carried NO feature set, so it inherited a supports() that accepts
    # everything; MATLAB had the mirror-image defect, a set naming no service
    # distribution, which refused every model including a plain closed
    # exponential network. Both now declare the laws the bounds consume.
    from line_solver.solvers.solver_ba import SolverBA
    assert SolverBA.supports(closed_ps(5))

    # A renewal law is admissible whatever its higher moments: the analyzer
    # reads sn.rates and sn.visits and nothing else.
    m = Network('closedErl')
    q1 = Queue(m, 'q1', SchedStrategy.FCFS)
    q2 = Queue(m, 'q2', SchedStrategy.FCFS)
    c = ClosedClass(m, 'c', 3, q1, 0)
    q1.setService(c, Erlang.fitMeanAndOrder(1, 2))
    q2.setService(c, Exp(0.8))
    m.link(Network.serialRouting(q1, q2))
    assert SolverBA.supports(m)


def _ba_modulated():
    # Single-class closed, one server each, MODULATED service at the first station.
    m = Network('closedMAP')
    q1 = Queue(m, 'q1', SchedStrategy.FCFS)
    q2 = Queue(m, 'q2', SchedStrategy.FCFS)
    c = ClosedClass(m, 'c', 3, q1, 0)
    q1.setService(c, MMPP2(0.2, 0.6, 0.1, 0.2))
    q2.setService(c, Exp(0.8))
    m.link(Network.serialRouting(q1, q2))
    return m


def test_a_modulated_service_process_reaches_mapamva_only():
    # THIS ASSERTED THE OPPOSITE UNTIL 'mapamva' LANDED (2026-09-04), and the
    # premise rather than the rule is what changed. A modulated process was out
    # of the SOLVER envelope because every family here was a function of the
    # demands D = V./rates: the mean rate exists, so the utilization law holds
    # and the formula returns a number, but that number brackets a DIFFERENT
    # system, one whose successive services are independent.
    #
    # MAP-AMVA (Casale-Smirni, DSN 2009) is derived FOR the correlated model,
    # its LP variables being the per-phase QN(i,k) and UN(i,k), so the same
    # reasoning that refuses the others admits it. 'MAP'/'MMPP2' therefore moved
    # INTO the base envelope, and getMethodFeatureSet takes them back from every
    # other family. The direction is forced: a feature set refuses a model for
    # HAVING a construct and never for lacking one, so a law can only be granted
    # to one family by removing it from the rest.
    #
    # Keeping the old assertion would have made the new family unreachable and
    # hidden SolverBA from model.help() on exactly the models it was written
    # for, so the guarantee is re-pinned here at the level where it now lives:
    # the SOLVER accepts the model, and every family but mapamva still refuses
    # it BY NAME. Mirrors SolverBATest.testAModulatedServiceProcessReachesMapamvaOnly
    # in the JAR.
    from line_solver.solvers.solver_ba import SolverBA
    m = _ba_modulated()
    solver = SolverBA(m, verbose=False)
    assert SolverBA.supports(m)

    supported = sorted(name for name in solver.listValidMethods()
                       if solver.supportsModelMethod(name)[0])
    assert supported == ['mapamva.lower', 'mapamva.upper']


def test_a_renewal_family_still_refuses_modulated_service_by_name():
    # The half of the old assertion that must NOT weaken: a demand-parameterized
    # family asked for by name on a modulated model is refused, and the reason
    # names the offending law rather than silently bounding a different system.
    from line_solver.solvers.solver_ba import SolverBA
    solver = SolverBA(_ba_modulated(), verbose=False)
    for name in ('gb.upper', 'gb.lower', 'aba.upper', 'lr.upper'):
        ok, reason = solver.supportsModelMethod(name)
        assert not ok, name + ' no longer refuses a modulated service process'
        assert reason


def test_qns_refuses_a_binding_buffer():
    # Nothing under the QNS tree reads sn.cap: qnsolver's MVA-family algorithms
    # have no representation of a finite buffer, so a capped station was solved
    # as an unbounded one and reported the unconstrained answer.
    from line_solver.solvers.wrappers.solver_qns import SolverQNS
    ok, reason = SolverQNS(bas_blocking()).supportsModelMethod('default')
    assert not ok
    assert 'capacity' in reason and 'Queue2' in reason
    ok_plain, _ = SolverQNS(closed_ps(5)).supportsModelMethod('default')
    assert ok_plain


def test_auto_offers_bounds_but_not_the_solvers_that_cannot_run_the_model():
    # The three gates above, seen through the list a caller enumerates.
    blocked = set(SolverAUTO(bas_blocking(), verbose=False).listValidMethods())
    plain = set(SolverAUTO(closed_ps(5), verbose=False).listValidMethods())

    # SolverBA is offered on the model its bounds are derived for ...
    assert any(t.startswith('ba.') for t in plain)
    # ... and on a blocked model only the qrf.bas family, which carries the
    # blocking tables explicitly, survives.
    assert {t for t in blocked if t.startswith('ba.')}
    assert all('bas' in t or 'rsrd' in t or t == 'ba.default'
               for t in blocked if t.startswith('ba.'))
    # qnsolver and NC have no representation of the buffer at all.
    assert not [t for t in blocked if t.startswith('qns')]
    assert not [t for t in blocked if t.startswith('nc.')]
    # UQ is not a solver for a model with no uncertain parameter, in either.
    assert not [t for t in blocked if t.startswith('uq')]
    assert not [t for t in plain if t.startswith('uq')]
