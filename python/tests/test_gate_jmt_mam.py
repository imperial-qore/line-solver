"""The per-method support gates of SolverJMT and SolverMAM.

WHAT THIS PINS, and why each case is here rather than left to the run.

findSolver / model.help() reports one row per (solver, method) pair, and the
row is produced by asking the solver supportsModelMethod(method) -- the same
gate SolverAUTO.chooseSolverRanked uses before delegating and listValidMethods
projects. A gate that is weaker than the analyzer offers a pair that then
RAISES, or worse returns a table of ZEROS under a method labelled exact. Both
were live:

  * SolverJMT declared ONE envelope for two engines. 'jmva.*' drives the JMVA
    ANALYTICAL engine, whose document carries a station type, a per-chain
    demand, a per-chain visit count, the populations and a reference station --
    no cache, no fork, no region, no discipline. On a three-class LRU cache
    model all eight closed-form jmva methods returned an entirely zero table
    with no error, jmva.mva labelled 'exact' among them.
  * 'dec.mmap' is an open-network departure-process fixed point, and the gate
    knew neither that nor its discipline list, so it was offered on every
    closed model.
  * 'retrial' needs an impatience configuration to analyze, which is a MUST BE
    PRESENT rule no feature set can state.

The assertions below are about the gate, not about numbers, except the one
that runs every jmt row the report offers on a cache model and requires a
non-zero answer -- which is the defect itself, stated directly.
"""

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from line_solver import (Network, Source, Queue, Sink, Delay, Cache, Fork, Join, OpenClass,
                         ClosedClass, Exp, Immediate, Zipf, ReplacementStrategy,
                         SchedStrategy, SolverJMT, SolverMAM)
from line_solver.solvers.wrappers.solver_jmt.solver_jmt import (
    jmt_method_refusal, jmva_is_closed_only)


JMVA_CLOSED_ONLY = ['jmva.amva', 'jmva.recal', 'jmva.comom', 'jmva.chow',
                    'jmva.bs', 'jmva.aql', 'jmva.lin', 'jmva.dmlin']
JMVA_EXACT = ['jmva', 'jmva.mva']


# ---------------------------------------------------------------------------
# models
# ---------------------------------------------------------------------------

def mm1():
    """Source -> FCFS Queue -> Sink: open, single server, product form."""
    m = Network('mm1')
    s, q, k = Source(m, 'Source'), Queue(m, 'Queue', SchedStrategy.FCFS), Sink(m, 'Sink')
    c = OpenClass(m, 'C1')
    s.setArrival(c, Exp(1.0))
    q.setService(c, Exp(2.0))
    m.link(Network.serialRouting(s, q, k))
    return m


def closed_single_server():
    """Delay -> FCFS Queue, N = 3: the shape every jmva algorithm serves."""
    m = Network('repairmen')
    d, q = Delay(m, 'Delay'), Queue(m, 'Queue', SchedStrategy.FCFS)
    c = ClosedClass(m, 'C1', 3, d)
    d.setService(c, Exp(1.0))
    q.setService(c, Exp(2.0))
    m.link(Network.serialRouting(d, q))
    return m


def closed_multiserver():
    """The same with c = 3, which the closed-form jmva algorithms refuse."""
    m = Network('multiserver')
    d, q = Delay(m, 'Delay'), Queue(m, 'Queue', SchedStrategy.FCFS)
    q.setNumberOfServers(3)
    c = ClosedClass(m, 'C1', 4, d)
    d.setService(c, Exp(1.0))
    q.setService(c, Exp(2.0))
    m.link(Network.serialRouting(d, q))
    return m


def open_forkjoin():
    """Source -> Fork -> two FCFS queues -> Join -> Sink, one open class.

    dec.mmap's sweep uses the PLAIN traffic step, which has no synchronization;
    the topology router sends this shape to the FJ solver from 'default' and
    'dec.source' and never to dec.mmap.
    """
    m = Network('fj')
    s, f = Source(m, 'S'), Fork(m, 'F')
    q1 = Queue(m, 'Q1', SchedStrategy.FCFS)
    q2 = Queue(m, 'Q2', SchedStrategy.FCFS)
    j, k = Join(m, 'J'), Sink(m, 'K')
    c = OpenClass(m, 'C')
    s.setArrival(c, Exp(0.5))
    q1.setService(c, Exp(2.0))
    q2.setService(c, Exp(2.0))
    P = m.initRoutingMatrix()
    P.set(c, c, s, f, 1.0)
    P.set(c, c, f, q1, 1.0)
    P.set(c, c, f, q2, 1.0)
    P.set(c, c, q1, j, 1.0)
    P.set(c, c, q2, j, 1.0)
    P.set(c, c, j, k, 1.0)
    m.link(P)
    return m


def mixed_open_closed():
    """Source -> Delay -> PS Queue -> Sink with ONE OPEN and ONE CLOSED class.

    The shape the binding-buffer rule got wrong. Nobody caps anything here, but
    refreshCapacity DERIVES a station capacity from the classes served, and in
    the JAR an unbounded contribution is Integer.MAX_VALUE rather than Inf --
    so this station came out as 2147483647 + 2 = 2147483649 and was read as a
    binding buffer, refusing five MixedExamplesTest models that have no buffer
    at all. This is jline.examples.java.basic.MixedModel.mqn_basic.
    """
    m = Network('mqn_basic')
    d = Delay(m, 'Delay')
    q = Queue(m, 'Queue1', SchedStrategy.PS)
    s, k = Source(m, 'Source'), Sink(m, 'Sink')
    cc = ClosedClass(m, 'ClosedClass', 2, d, 0)
    oc = OpenClass(m, 'OpenClass', 0)
    d.setService(cc, Exp(1.0))
    d.setService(oc, Exp(3.0))
    q.setService(cc, Exp(2.0))
    q.setService(oc, Exp(1.0))
    s.setArrival(oc, Exp(0.1))
    P = m.initRoutingMatrix()
    P.set(cc, cc, d, q, 1.0)
    P.set(cc, cc, q, d, 1.0)
    P.set(oc, oc, s, d, 1.0)
    P.set(oc, oc, d, q, 1.0)
    P.set(oc, oc, q, k, 1.0)
    m.link(P)
    return m


def cqn_two_class():
    """Two closed classes at one PS queue: sn.cap is DERIVED as 2 x N here.

    The multi-class shape is the one that makes "finite cap" and "binding cap"
    different questions, so it is the model the buffer rule must leave alone.
    """
    m = Network('cqn2')
    d, q = Delay(m, 'D'), Queue(m, 'Q', SchedStrategy.PS)
    c1, c2 = ClosedClass(m, 'C1', 2, d), ClosedClass(m, 'C2', 1, d)
    d.setService(c1, Exp(1.0))
    d.setService(c2, Exp(2.0))
    q.setService(c1, Exp(2.0))
    q.setService(c2, Exp(3.0))
    m.link(Network.serialRouting(d, q))
    return m


def open_hol():
    """Two open classes at a head-of-line priority queue.

    JMVA emits NO discipline, so it would answer the priority-blind numbers for
    this model; SolverNC and SolverMVA draw the same line by not declaring HOL.
    """
    m = Network('prio')
    s, q, k = Source(m, 'Source'), Queue(m, 'Queue', SchedStrategy.HOL), Sink(m, 'Sink')
    hi, lo = OpenClass(m, 'Hi', 0), OpenClass(m, 'Lo', 1)
    s.setArrival(hi, Exp(0.4))
    s.setArrival(lo, Exp(0.4))
    q.setService(hi, Exp(2.0))
    q.setService(lo, Exp(2.0))
    m.link(Network.serialRouting(s, q, k))
    return m


def closed_binding_buffer():
    """Delay -> FCFS Queue with cap 2 and N = 4: the buffer BINDS.

    LINE blocks a closed job that finds no room; no JMT drop strategy reproduces
    that, and the JMVA document has no capacity element at all.
    """
    m = Network('blk')
    d, q = Delay(m, 'D'), Queue(m, 'Q', SchedStrategy.FCFS)
    q.setCapacity(2)
    c = ClosedClass(m, 'C', 4, d)
    d.setService(c, Exp(1.0))
    q.setService(c, Exp(2.0))
    m.link(Network.serialRouting(d, q))
    return m


def open_binding_buffer():
    """M/M/1/2: an OPEN binding buffer, which JSIM DOES simulate.

    JMT's queue section carries the drop rule directly, so a refused arrival is
    lost in JMT exactly as it is in LINE. This is the over-tightening guard.
    """
    m = Network('mm1k')
    s, q, k = Source(m, 'S'), Queue(m, 'Q', SchedStrategy.FCFS), Sink(m, 'K')
    q.setCapacity(2)
    c = OpenClass(m, 'C')
    s.setArrival(c, Exp(1.0))
    q.setService(c, Exp(2.0))
    m.link(Network.serialRouting(s, q, k))
    return m


def cachemodel():
    """Client -> LRU Cache -> CacheDelay, hit and miss classes.

    The model on which every jmva method used to return a table of zeros.
    """
    m = Network('cache')
    client = Delay(m, 'Client')
    cache = Cache(m, 'Cache', 4, 2, ReplacementStrategy.LRU)
    hitmiss = Delay(m, 'CacheDelay')
    cc = ClosedClass(m, 'ClientClass', 1, client, 0)
    hc = ClosedClass(m, 'HitClass', 0, client, 0)
    mc = ClosedClass(m, 'MissClass', 0, client, 0)
    client.setService(cc, Immediate())
    hitmiss.setService(hc, Exp.fitMean(0.2))
    hitmiss.setService(mc, Exp.fitMean(1.0))
    cache.setRead(cc, Zipf(1.4, 4))
    cache.setHitClass(cc, hc)
    cache.setMissClass(cc, mc)
    P = m.initRoutingMatrix()
    P.set(cc, cc, client, cache, 1.0)
    P.set(hc, hc, cache, hitmiss, 1.0)
    P.set(mc, mc, cache, hitmiss, 1.0)
    P.set(hc, cc, hitmiss, client, 1.0)
    P.set(mc, cc, hitmiss, client, 1.0)
    m.link(P)
    return m


def offered(model, family):
    """The methods findSolver offers for one solver family, bare of its prefix."""
    table = model.findSolver()
    rows = table[table.Solver == family]
    return [str(m).split('.', 1)[1] for m in rows.Method]


jmt_available = pytest.mark.skipif(not SolverJMT.isAvailable(),
                                   reason='JMT.jar is not present in common/')


# ---------------------------------------------------------------------------
# SolverJMT: the two engines have two envelopes
# ---------------------------------------------------------------------------

class TestJmvaEnvelope:

    def test_jmva_envelope_is_a_strict_subset_of_the_jsim_one(self):
        jsim = SolverJMT.getFeatureSet()
        jmva = SolverJMT.getJMVAFeatureSet()
        assert jmva < jsim
        # The constructs the JMVA document has no element for at all.
        for feat in ('Fork', 'Join', 'Region', 'Place', 'Transition',
                     'Reneging', 'Balking', 'HeteroServers', 'SetupDelayOff',
                     'ServerParallelism'):
            assert feat in jsim
            assert feat not in jmva
        # The writer emits no discipline, so only the BCMP station types survive.
        for feat in ('SchedStrategy_HOL', 'SchedStrategy_DPS', 'SchedStrategy_GPS',
                     'SchedStrategy_POLLING', 'SchedStrategy_SRPT'):
            assert feat not in jmva
        for feat in ('SchedStrategy_INF', 'SchedStrategy_PS', 'SchedStrategy_FCFS',
                     'SchedStrategy_LCFSPR'):
            assert feat in jmva
        # A mean demand is all JMVA reads, so the laws stay.
        for feat in ('Exp', 'Erlang', 'HyperExp', 'Pareto', 'Replayer'):
            assert feat in jmva

    def test_the_closed_only_algorithms_drop_three_more_names(self):
        # 2026-09-05: MultiServer joined the two. SolverJMT.m:255-261 withdraws
        # it from the closed-only algorithms because JMVA runs them
        # single-server, and the sibling test below
        # (test_the_closed_only_algorithms_refuse_a_multiserver_model) is what
        # requires that withdrawal, so leaving it out made the two disagree.
        base = SolverJMT.getJMVAFeatureSet()
        probe = SolverJMT(closed_single_server(), 'jmva')
        assert probe.getMethodFeatureSet('jmva') == base
        for method in JMVA_CLOSED_ONLY:
            assert jmva_is_closed_only(method)
            narrowed = probe.getMethodFeatureSet(method)
            assert narrowed == base - {'OpenClass', 'LoadDependence', 'MultiServer'}
        for method in JMVA_EXACT:
            assert not jmva_is_closed_only(method)

    def test_the_simulation_methods_keep_the_whole_jsim_envelope(self):
        probe = SolverJMT(mm1())
        for method in ('default', 'jsim', 'replication'):
            assert probe.getMethodFeatureSet(method) == SolverJMT.getFeatureSet()


class TestJmtGateRefusals:

    def test_no_jmva_method_is_offered_on_a_cache_model(self):
        model = cachemodel()
        probe = SolverJMT(model)
        for method in JMVA_EXACT + JMVA_CLOSED_ONLY:
            ok, reason = probe.supportsModelMethod(method)
            assert not ok, '%s was offered on a cache model' % method
            assert 'Cache' in reason

    @jmt_available
    def test_no_jmt_row_the_report_offers_returns_an_all_zero_table(self):
        """The defect, stated directly: an offered pair must produce an answer.

        Eight jmva methods used to be offered on the cache model and each
        returned an entirely zero QLen column with no error, jmva.mva labelled
        'exact' among them. The cache half of this sweep is empty in THIS port,
        because its JSIM writer has no cache branch either and the JSIM names
        are withheld too; it is kept because the rule is the contract, and
        because the other three codebases do declare Cache for JSIM and run
        those rows here.
        """
        for model in (cachemodel(), closed_single_server(), mm1()):
            for method in offered(model, 'jmt'):
                solver = SolverJMT(model, method, samples=10000, seed=23000, verbose=0)
                Q = np.asarray(solver.getAvgQLen(), dtype=float)
                assert Q.size, 'jmt.%s produced no queue lengths' % method
                assert not np.allclose(Q, 0.0), 'jmt.%s returned an all-zero table' % method

    def test_closed_only_algorithms_refuse_an_open_model(self):
        probe = SolverJMT(mm1())
        for method in JMVA_CLOSED_ONLY:
            ok, reason = probe.supportsModelMethod(method)
            assert not ok, '%s was offered on an open model' % method
            assert 'OpenClass' in reason
        # and the exact engine still takes it
        for method in JMVA_EXACT:
            assert probe.supportsModelMethod(method)[0]

    def test_closed_only_algorithms_refuse_a_multiserver_model(self):
        model = closed_multiserver()
        probe = SolverJMT(model)
        for method in JMVA_CLOSED_ONLY:
            ok, reason = probe.supportsModelMethod(method)
            assert not ok, '%s was offered on a multi-server model' % method
            assert 'multi-server' in reason
        for method in JMVA_EXACT:
            assert probe.supportsModelMethod(method)[0]

    def test_the_gate_and_the_run_give_the_same_sentence(self):
        """One predicate, two callers: the analyzer must raise what the gate says."""
        model = closed_multiserver()
        _, reason = SolverJMT(model).supportsModelMethod('jmva.amva')
        assert reason == 'jmva.amva does not support multi-server stations.'
        with pytest.raises(Exception) as excinfo:
            SolverJMT(model, 'jmva.amva', verbose=0).getAvgQLen()
        assert 'multi-server' in str(excinfo.value)

    def test_a_priority_model_is_not_offered_the_analytical_engine(self):
        # JMVA writes no discipline, so it would answer the priority-blind
        # numbers under a method the report calls exact.
        probe = SolverJMT(open_hol())
        for method in JMVA_EXACT + JMVA_CLOSED_ONLY:
            ok, reason = probe.supportsModelMethod(method)
            assert not ok
            assert 'SchedStrategy_HOL' in reason
        # the simulator still serves it
        assert probe.supportsModelMethod('jsim')[0]

    def test_replication_needs_a_finite_timespan(self):
        model = mm1()
        ok, reason = SolverJMT(model).supportsModelMethod('replication')
        assert not ok
        assert 'finite timespan' in reason
        assert 'replication' not in offered(model, 'jmt')
        # ... and is admitted once the horizon is stated
        withspan = SolverJMT(model, 'replication', timespan=[0, 10])
        assert withspan.supportsModelMethod('replication')[0]

    def test_replication_raises_the_gate_sentence_when_named_by_hand(self):
        solver = SolverJMT(mm1(), 'replication')
        with pytest.raises(Exception) as excinfo:
            solver.runAnalyzer()
        assert 'finite timespan' in str(excinfo.value)

    def test_a_closed_binding_buffer_is_refused_by_both_engines(self):
        """NEITHER engine carries a binding buffer, for opposite reasons.

        JSIM because no JMT drop strategy reproduces LINE's blocking -- "waiting
        queue" does not enforce the size at all and "BAS blocking" completes the
        service first -- and JMVA because its document has no capacity element,
        so it answered 2.19 jobs at a station that can hold 2 (exact: 1.33).
        """
        model = closed_binding_buffer()
        probe = SolverJMT(model)
        for method in ('default', 'jsim'):
            ok, reason = probe.supportsModelMethod(method)
            assert not ok, 'jmt.%s was offered on a binding buffer' % method
            assert 'binds for the closed class' in reason
        for method in JMVA_EXACT + JMVA_CLOSED_ONLY:
            ok, reason = probe.supportsModelMethod(method)
            assert not ok, 'jmt.%s was offered on a binding buffer' % method
            assert 'no capacity element' in reason
        assert offered(model, 'jmt') == []

    def test_the_gate_and_the_writer_give_the_same_buffer_sentence(self):
        """One predicate, two callers: the JSIM writer's own rule IS the gate."""
        from line_solver.api.solvers.jmt.handler import _jmt_station_cap_assert
        model = closed_binding_buffer()
        _, reason = SolverJMT(model).supportsModelMethod('jsim')
        with pytest.raises(ValueError) as writer:
            _jmt_station_cap_assert(model.getStruct(), 1)
        assert str(writer.value) == reason

    def test_an_open_loss_buffer_keeps_the_simulator(self):
        """The over-tightening guard: JSIM genuinely simulates an M/M/1/K.

        A refused OPEN arrival is lost, which JMT's queue section expresses
        directly, so the row must survive -- and answer the CONSTRAINED value.
        """
        model = open_binding_buffer()
        assert 'jsim' in offered(model, 'jmt')
        assert 'default' in offered(model, 'jmt')
        # the analytical engine still cannot: it has no capacity element at all
        for method in JMVA_EXACT:
            ok, reason = SolverJMT(model).supportsModelMethod(method)
            assert not ok
            assert 'no capacity element' in reason

    @jmt_available
    def test_the_open_loss_buffer_answers_the_constrained_value(self):
        model = open_binding_buffer()
        Q = np.asarray(SolverJMT(model, 'jsim', samples=20000, seed=23000,
                                 verbose=0).getAvgQLen(), dtype=float).ravel()
        # M/M/1/2 at rho = 0.5 holds 0.56 jobs; the UNCONSTRAINED M/M/1 holds 1.
        assert Q[-1] < 0.8

    def test_a_mixed_open_closed_model_has_no_binding_buffer(self):
        """Nobody capped anything, so every jmt method must survive.

        refreshCapacity DERIVES a capacity from the classes served at a station,
        and a station carrying one open and one closed class gets a mixture of
        an unbounded contribution and a population. Reading that as a finite,
        binding buffer refused five MixedExamplesTest models in the JAR, where
        the unbounded contribution is Integer.MAX_VALUE rather than Inf.
        """
        model = mixed_open_closed()
        probe = SolverJMT(model)
        for method in ('default', 'jsim'):
            ok, reason = probe.supportsModelMethod(method)
            assert ok, 'jmt.%s was refused on a model with no buffer: %s' % (method, reason)
        for method in JMVA_EXACT:
            ok, reason = probe.supportsModelMethod(method)
            assert ok, 'jmt.%s was refused on a model with no buffer: %s' % (method, reason)
        methods = offered(model, 'jmt')
        assert 'jsim' in methods and 'jmva' in methods
        # the closed-only algorithms go for the OPEN class, not for a buffer
        for method in JMVA_CLOSED_ONLY:
            ok, reason = probe.supportsModelMethod(method)
            assert not ok
            assert 'OpenClass' in reason

    def test_a_source_is_never_reported_as_carrying_a_buffer(self):
        """A Source has no buffer that can bind: it IS the external world.

        refreshCapacity still writes it a capacity row, and in the JAR that row
        is a sum of unbounded sentinels, which is how a refusal came to name a
        Source. Excluded on node type, so the rule can never speak about one.
        """
        for model in (mixed_open_closed(), open_binding_buffer(), mm1()):
            for method in ('jsim',) + tuple(JMVA_EXACT):
                _, reason = SolverJMT(model).supportsModelMethod(method)
                assert 'Source' not in reason, reason

    def test_an_uncapped_model_is_untouched_by_the_buffer_rule(self):
        """refreshCapacity DERIVES a finite sn.cap for a station nobody capped.

        The rule therefore tests the cap against the population that can REACH
        the station, not against Inf; reading "finite" as "binding" would refuse
        every closed model.
        """
        for model in (closed_single_server(), cqn_two_class(), mixed_open_closed()):
            methods = offered(model, 'jmt')
            assert 'jsim' in methods
            for method in JMVA_EXACT:
                assert method in methods, '%s lost on an uncapped model' % method


class TestJmtIsNotOverTightened:
    """The split must not cost JMVA the models it genuinely solves."""

    def test_a_closed_single_server_network_keeps_every_jmva_method(self):
        methods = offered(closed_single_server(), 'jmt')
        for method in JMVA_EXACT + JMVA_CLOSED_ONLY:
            assert method in methods, '%s was lost on a closed product-form model' % method
        assert 'jsim' in methods

    def test_an_open_product_form_network_keeps_the_exact_engine(self):
        methods = offered(mm1(), 'jmt')
        for method in JMVA_EXACT:
            assert method in methods
        assert 'jsim' in methods and 'default' in methods

    def test_a_multiserver_model_keeps_the_load_dependent_arm(self):
        methods = offered(closed_multiserver(), 'jmt')
        for method in JMVA_EXACT:
            assert method in methods

    @jmt_available
    def test_the_kept_methods_still_answer(self):
        model = closed_single_server()
        for method in JMVA_EXACT + JMVA_CLOSED_ONLY:
            Q = np.asarray(SolverJMT(model, method, verbose=0).getAvgQLen(), dtype=float)
            assert Q.size and not np.allclose(Q, 0.0)


def test_jmt_method_refusal_is_silent_on_an_admissible_pair():
    sn = closed_single_server().getStruct()
    assert jmt_method_refusal(sn, 'jmva.lin', SolverJMT.defaultOptions()) == ''
    assert jmt_method_refusal(sn, 'jsim', SolverJMT.defaultOptions()) == ''


# ---------------------------------------------------------------------------
# SolverMAM: dec.mmap and retrial
# ---------------------------------------------------------------------------

class TestMamDecMmap:

    def test_dec_mmap_is_refused_on_a_closed_model(self):
        model = closed_single_server()
        ok, reason = SolverMAM(model).supportsModelMethod('dec.mmap')
        assert not ok
        assert 'ClosedClass' in reason
        assert 'dec.mmap' not in offered(model, 'mam')

    def test_dec_mmap_is_refused_on_a_delay_station(self):
        # A Delay is an INF station, which the analyzer's opening loop has no
        # branch for; the reason has to name it rather than the class mix alone.
        model = closed_single_server()
        _, reason = SolverMAM(model).supportsModelMethod('dec.mmap')
        assert 'SchedStrategy_INF' in reason

    def test_dec_mmap_raises_rather_than_answering_zeros(self):
        """The analyzer used to warn and return the zero matrices it had built.

        In this port that empty result reached MAMResult, which died with four
        missing positional arguments; in MATLAB and the JAR it was reported as a
        table of zeros. The gate now stops the run first, so the analyzer is
        asked directly here -- it is the backstop for a caller who reaches it
        with the checks disabled, and a silent zero there is still a wrong
        answer.
        """
        from line_solver.api.solvers.mam.handler import (
            solver_mam as handler_solver_mam, SolverMAMOptions as HandlerOptions)
        with pytest.raises(Exception):
            SolverMAM(closed_single_server(), 'dec.mmap').getAvgQLen()

        opts = HandlerOptions()
        opts.method = 'dec.mmap'
        # A CLOSED model with no Delay, so the class-mix refusal is the one that
        # fires: closed_single_server() would be stopped by the INF station
        # first and would not exercise this branch.
        loop = Network('closed_loop')
        q1 = Queue(loop, 'Q1', SchedStrategy.FCFS)
        q2 = Queue(loop, 'Q2', SchedStrategy.FCFS)
        c = ClosedClass(loop, 'C1', 2, q1)
        q1.setService(c, Exp(1.0))
        q2.setService(c, Exp(2.0))
        loop.link(Network.serialRouting(q1, q2))
        with pytest.raises(ValueError) as closed:
            handler_solver_mam(loop.getStruct(), opts)
        assert 'open models only' in str(closed.value)

        # An open model whose station the ladder has no branch for: a Delay in
        # an otherwise open network reaches the same analyzer.
        m = Network('open_with_delay')
        s, d, k = Source(m, 'S'), Delay(m, 'Think'), Sink(m, 'K')
        c = OpenClass(m, 'C1')
        s.setArrival(c, Exp(1.0))
        d.setService(c, Exp(2.0))
        m.link(Network.serialRouting(s, d, k))
        with pytest.raises(ValueError) as sched:
            handler_solver_mam(m.getStruct(), opts)
        assert 'scheduling strategy' in str(sched.value)

    def test_dec_mmap_is_refused_on_a_fork_join_topology(self):
        """The sweep has no synchronization, so a Fork is not a smaller model.

        The topology router sends this shape to the FJ solver from 'default'
        and 'dec.source'; dec.mmap has no such route, and the departure process
        it built here had no recurrent state.
        """
        model = open_forkjoin()
        ok, reason = SolverMAM(model).supportsModelMethod('dec.mmap')
        assert not ok
        assert 'Fork' in reason
        assert 'dec.mmap' not in offered(model, 'mam')
        # the routes that ARE written for it stay
        assert 'default' in offered(model, 'mam')
        assert 'dec.source' in offered(model, 'mam')

    def test_dec_mmap_still_serves_the_open_models_it_is_written_for(self):
        for model in (mm1(), open_hol()):
            ok, reason = SolverMAM(model).supportsModelMethod('dec.mmap')
            assert ok, reason
            assert 'dec.mmap' in offered(model, 'mam')
            Q = np.asarray(SolverMAM(model, 'dec.mmap').getAvgQLen(), dtype=float)
            assert Q.size and not np.allclose(Q, 0.0)


class TestMamRetrial:

    def test_retrial_is_refused_without_an_impatience_configuration(self):
        for model in (mm1(), closed_single_server(), open_hol()):
            ok, reason = SolverMAM(model).supportsModelMethod('retrial')
            assert not ok
            assert 'impatience' in reason
            assert 'retrial' not in offered(model, 'mam')

    def test_the_gate_and_the_run_give_the_same_sentence(self):
        model = mm1()
        _, reason = SolverMAM(model).supportsModelMethod('retrial')
        with pytest.raises(Exception) as excinfo:
            SolverMAM(model, 'retrial').getAvgQLen()
        # runAnalyzerChecks prefixes the reason; the sentence itself is the
        # gate's, because both callers ask mam_retrial_applicable.
        assert reason in str(excinfo.value)
        # and the analyzer arm raises it too, for a caller who reaches it with
        # the checks disabled
        probe = SolverMAM(model, 'retrial')
        with pytest.raises(ValueError) as direct:
            probe._solve_retrial_reneging('retrial')
        assert str(direct.value) == reason

    def test_the_reason_names_the_requirement_that_is_missing(self):
        # qsys_is_retrial reports WHICH requirement failed, and carrying it
        # through is the difference between a bare no and a usable answer.
        _, open_reason = SolverMAM(mm1()).supportsModelMethod('retrial')
        assert 'retrial queue' in open_reason or 'orbit' in open_reason
        _, closed_reason = SolverMAM(closed_single_server()).supportsModelMethod('retrial')
        assert 'open queueing model' in closed_reason


def test_no_mam_row_the_report_offers_raises_or_returns_zeros():
    """The whole family, on the models that used to break it."""
    for model in (mm1(), closed_single_server(), open_hol(), cachemodel(),
                  closed_binding_buffer(), open_forkjoin(), mixed_open_closed()):
        for method in offered(model, 'mam'):
            Q = np.asarray(SolverMAM(model, method).getAvgQLen(), dtype=float)
            assert Q.size, 'mam.%s produced no queue lengths' % method
            assert not np.allclose(Q, 0.0), 'mam.%s returned an all-zero table' % method
