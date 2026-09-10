"""
The per-method support gates of SolverCTMC, SolverFLD and SolverSSA.

WHAT THIS PINS. `model.help()` / `model.findSolver()` reports one row per
(solver, method) pair, and it asks `supportsModelMethod(method)` -- the same
gate `SolverAUTO.chooseSolverRanked` consults before delegating. Until the
rules below reached that gate it was much weaker than what the ANALYZERS
enforce at run time, so the report offered pairs that then raised:

  ctmc.cftp / ctmc.cftp.approx   single-class and closed-only
  ctmc.mdd                       single-class and closed-only, and no fork-join
  ctmc.default / exact / gpu     the state space has to fit memory
  ctmc.* and ssa.*               the fork-join model class sn_fj_validate admits
  fluid.diffusion                closed-only, and no fork-join
  fluid.kp                       no fork-join
  fluid.refined                  closed-only
  fluid.tbi                      closed-only, and no cache
  fluid.dae                      no OPEN fork-join model
  fluid.mol / fluid.mtginf       a finite options.timespan

Each test asserts the refusal AND its converse: a model the method IS derived
for must keep it. Over-tightening a gate hides a method the user could have
run, which is the same defect with the sign flipped -- and two of the cases
here ARE that flipped defect, found and fixed: native python's SolverSSA
declared no Fork/Join at all, and the C++ fluid envelope withheld them from
five methods that do integrate a fork-join model.
"""

import os
import sys
import unittest

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from line_solver import (CTMC, Cache, ClosedClass, Delay, Exp, FLD, Fork, GlobalConstants,
                         Immediate, Join, Network, OpenClass, Queue, ReplacementStrategy,
                         SSA, SchedStrategy, SolverCTMC, SolverFLD, SolverSSA, Sink, Source,
                         VerboseLevel, Zipf)

GlobalConstants.setVerbose(VerboseLevel.SILENT)


# ---------------------------------------------------------------------------
# models
# ---------------------------------------------------------------------------

def mm1():
    """Source -> Queue -> Sink, one open class: the smallest open model."""
    m = Network('mm1')
    s, q, k = Source(m, 'S'), Queue(m, 'Q', SchedStrategy.FCFS), Sink(m, 'K')
    c = OpenClass(m, 'C')
    s.setArrival(c, Exp(1.0))
    q.setService(c, Exp(2.0))
    m.link(Network.serialRouting(s, q, k))
    return m


def repairmen():
    """Delay -> Queue -> Delay, ONE closed class: what cftp and mdd are for."""
    m = Network('rep')
    d, q = Delay(m, 'D'), Queue(m, 'Q', SchedStrategy.FCFS)
    c = ClosedClass(m, 'C', 3, d)
    d.setService(c, Exp(1.0))
    q.setService(c, Exp(2.0))
    m.link(Network.serialRouting(d, q))
    return m


def cqn2():
    """The same shape with TWO closed classes: closed, but not single-class."""
    m = Network('cqn2')
    d, q = Delay(m, 'D'), Queue(m, 'Q', SchedStrategy.PS)
    c1, c2 = ClosedClass(m, 'C1', 2, d), ClosedClass(m, 'C2', 1, d)
    d.setService(c1, Exp(1.0))
    d.setService(c2, Exp(2.0))
    q.setService(c1, Exp(2.0))
    q.setService(c2, Exp(3.0))
    m.link(Network.serialRouting(d, q))
    return m


def prio():
    """Two open HOL classes: the state space the memory gate has to refuse."""
    m = Network('prio')
    s, q, k = Source(m, 'S'), Queue(m, 'Q', SchedStrategy.HOL), Sink(m, 'K')
    c1, c2 = OpenClass(m, 'Hi', 0), OpenClass(m, 'Lo', 1)
    s.setArrival(c1, Exp(0.4))
    s.setArrival(c2, Exp(0.4))
    q.setService(c1, Exp(2.0))
    q.setService(c2, Exp(2.0))
    m.link(Network.serialRouting(s, q, k))
    return m


def cache_model():
    """A closed cache model: three classes, so cftp and mdd must refuse it."""
    m = Network('cache')
    d = Delay(m, 'Client')
    c = Cache(m, 'Cache', 4, 2, ReplacementStrategy.LRU)
    cd = Delay(m, 'CacheDelay')
    cc = ClosedClass(m, 'ClientClass', 1, d, 0)
    hc = ClosedClass(m, 'HitClass', 0, d, 0)
    mc = ClosedClass(m, 'MissClass', 0, d, 0)
    d.setService(cc, Immediate())
    cd.setService(hc, Exp.fitMean(0.2))
    cd.setService(mc, Exp.fitMean(1.0))
    c.setRead(cc, Zipf(1.4, 4))
    c.setHitClass(cc, hc)
    c.setMissClass(cc, mc)
    P = m.initRoutingMatrix()
    P.set(cc, cc, d, c, 1.0)
    P.set(hc, hc, c, cd, 1.0)
    P.set(mc, mc, c, cd, 1.0)
    P.set(hc, cc, cd, d, 1.0)
    P.set(mc, cc, cd, d, 1.0)
    m.link(P)
    return m


def forkjoin(paired=True, closed=True):
    """Fork -> two FCFS queues -> Join.

    `paired` names the Fork on the Join constructor, which is what DECLARES the
    fork-join pairing: sn.fj is read off that declaration and not derived from
    the routing, because a nested model (examples/basic/forkJoin/fj_basic_nesting)
    has two forks and two joins the routing alone does not pair. A Join built
    without it leaves sn.fj empty, which the exact construction refuses.
    """
    m = Network('fj')
    f = Fork(m, 'F')
    q1, q2 = Queue(m, 'Q1', SchedStrategy.FCFS), Queue(m, 'Q2', SchedStrategy.FCFS)
    j = Join(m, 'J', f) if paired else Join(m, 'J')
    if closed:
        d = Delay(m, 'D')
        c = ClosedClass(m, 'C', 2, d)
        d.setService(c, Exp(1.0))
        entry, exit_ = d, d
    else:
        s, k = Source(m, 'S'), Sink(m, 'K')
        c = OpenClass(m, 'C')
        s.setArrival(c, Exp(0.5))
        entry, exit_ = s, k
    q1.setService(c, Exp(2.0))
    q2.setService(c, Exp(2.0))
    P = m.initRoutingMatrix()
    P.set(c, c, entry, f, 1.0)
    P.set(c, c, f, q1, 1.0)
    P.set(c, c, f, q2, 1.0)
    P.set(c, c, q1, j, 1.0)
    P.set(c, c, q2, j, 1.0)
    P.set(c, c, j, exit_, 1.0)
    m.link(P)
    return m


def row(table, method):
    """The findSolver row for a 'family.method' method name, or None."""
    hits = table[table.Method == method]
    return None if len(hits) == 0 else hits.iloc[0]


# ---------------------------------------------------------------------------
# SolverCTMC
# ---------------------------------------------------------------------------

class TestCtmcCftpGate(unittest.TestCase):

    def test_an_open_model_is_refused_in_the_sampler_s_own_words(self):
        ok, reason = SolverCTMC(mm1(), method='cftp').supportsModelMethod('cftp')
        self.assertFalse(ok)
        self.assertIn('closed models only', reason)

    def test_a_multiclass_model_is_refused_and_the_class_count_is_named(self):
        for method in ('cftp', 'cftp.approx'):
            ok, reason = SolverCTMC(cqn2(), method=method).supportsModelMethod(method)
            self.assertFalse(ok, method)
            self.assertIn('single-class models only', reason)
            self.assertIn('2 classes', reason)

    def test_a_single_class_closed_model_keeps_both_cftp_variants(self):
        # The converse: refusing this one would hide a method that runs.
        for method in ('cftp', 'cftp.approx'):
            ok, reason = SolverCTMC(repairmen(), method=method).supportsModelMethod(method)
            self.assertTrue(ok, '%s: %s' % (method, reason))

    def test_the_gate_and_the_analyzer_are_one_predicate(self):
        # Whatever the gate refuses, the run refuses with the SAME sentence.
        model = cqn2()
        _, reason = SolverCTMC(model, method='cftp').supportsModelMethod('cftp')
        with self.assertRaises((ValueError, RuntimeError)) as ctx:
            SolverCTMC(model, method='cftp', samples=100).getAvgQLen()
        self.assertIn(reason, str(ctx.exception))

    def test_the_envelope_drops_what_the_registry_can_name(self):
        feats = SolverCTMC(repairmen(), method='cftp').getMethodFeatureSet('cftp')
        base = SolverCTMC.getFeatureSet()
        self.assertNotIn('OpenClass', feats)
        self.assertNotIn('Source', feats)
        self.assertNotIn('Cache', feats)
        self.assertNotIn('SchedStrategy_HOL', feats)
        self.assertNotIn('Region', feats)
        self.assertNotIn('LoadDependence', feats)
        self.assertNotIn('RoutingStrategy_JSQ', feats)
        # ... and keeps the product-form core it does serve.
        for kept in ('ClosedClass', 'Queue', 'Delay', 'Exp',
                     'SchedStrategy_FCFS', 'SchedStrategy_PS', 'SchedStrategy_INF',
                     'SchedStrategy_SIRO', 'SchedStrategy_LCFSPR'):
            self.assertIn(kept, feats, kept)
            self.assertIn(kept, base, kept)


class TestCtmcMddGate(unittest.TestCase):

    def test_an_open_model_is_refused(self):
        ok, reason = SolverCTMC(mm1(), method='mdd').supportsModelMethod('mdd')
        self.assertFalse(ok)
        self.assertTrue('CLOSED networks' in reason or 'OpenClass' in reason, reason)

    def test_a_multiclass_model_is_refused_and_the_class_count_is_named(self):
        ok, reason = SolverCTMC(cqn2(), method='mdd').supportsModelMethod('mdd')
        self.assertFalse(ok)
        self.assertIn('single-class networks', reason)
        self.assertIn('2 classes', reason)

    def test_a_multiclass_cache_model_is_refused(self):
        ok, reason = SolverCTMC(cache_model(), method='mdd').supportsModelMethod('mdd')
        self.assertFalse(ok)
        self.assertIn('single-class networks', reason)

    def test_a_single_class_closed_model_keeps_mdd(self):
        ok, reason = SolverCTMC(repairmen(), method='mdd').supportsModelMethod('mdd')
        self.assertTrue(ok, reason)

    def test_the_gate_and_the_analyzer_are_one_predicate(self):
        model = cqn2()
        _, reason = SolverCTMC(model, method='mdd').supportsModelMethod('mdd')
        with self.assertRaises(Exception) as ctx:
            SolverCTMC(model, method='mdd').getAvgQLen()
        self.assertIn(reason, str(ctx.exception))


class TestCtmcStateSpaceGate(unittest.TestCase):

    def test_an_intractable_chain_is_refused_by_the_state_space_methods(self):
        solver = SolverCTMC(prio(), method='default')
        tractable, _, _ = SolverCTMC.isStateSpaceTractable(prio(), solver.options)
        self.assertFalse(tractable, 'the fixture must be intractable for this to test anything')
        for method in ('default', 'exact', 'gpu'):
            ok, reason = SolverCTMC(prio(), method=method).supportsModelMethod(method)
            self.assertFalse(ok, method)
            self.assertIn('memory', reason.lower())

    def test_a_small_chain_keeps_the_state_space_methods(self):
        for method in ('default', 'exact', 'gpu'):
            ok, reason = SolverCTMC(repairmen(), method=method).supportsModelMethod(method)
            self.assertTrue(ok, '%s: %s' % (method, reason))

    def test_the_sampler_is_not_gated_on_a_state_space_it_never_builds(self):
        # cftp draws from the balance function, so the memory estimate that
        # stops 'default' says nothing about it. Its own refusal here is the
        # class count, not the state space.
        ok, reason = SolverCTMC(prio(), method='cftp').supportsModelMethod('cftp')
        self.assertFalse(ok)
        self.assertNotIn('memory', reason.lower())


class TestCtmcForkJoinGate(unittest.TestCase):
    """The fork-join model class, which EVERY ctmc method has to clear.

    The tag augmentation runs before the state space, the decision diagram and
    the sampler alike, so a model sn_fj_validate refuses is refused whichever
    method was asked for.
    """

    def test_an_undeclared_pairing_is_refused_for_every_method(self):
        model = forkjoin(paired=False, closed=True)
        for method in ('default', 'exact', 'gpu', 'mdd', 'cftp', 'cftp.approx'):
            ok, reason = SolverCTMC(model, method=method).supportsModelMethod(method)
            self.assertFalse(ok, method)
            self.assertTrue(reason, method)

    def test_the_pairing_refusal_is_the_validator_s_own_sentence(self):
        ok, reason = SolverCTMC(forkjoin(paired=False), method='default').supportsModelMethod('default')
        self.assertFalse(ok)
        self.assertIn('without a matched Join', reason)

    def test_an_open_class_through_a_fork_is_refused(self):
        # The pairing is declared here, so this is the SECOND rule of the same
        # validator: the exact construction is stated for closed chains.
        ok, reason = SolverCTMC(forkjoin(paired=True, closed=False),
                                method='default').supportsModelMethod('default')
        self.assertFalse(ok)
        self.assertIn('Open classes routed through a Fork', reason)

    def test_a_declared_closed_fork_join_keeps_the_state_space_methods(self):
        # The converse: the exact construction IS derived for this model, and
        # refusing it would hide the only exact answer a fork-join model has.
        model = forkjoin(paired=True, closed=True)
        for method in ('default', 'exact', 'gpu'):
            ok, reason = SolverCTMC(model, method=method).supportsModelMethod(method)
            self.assertTrue(ok, '%s: %s' % (method, reason))

    def test_mdd_refuses_a_fork_join_model_it_can_never_be_single_class_for(self):
        # The tag augmentation adds one auxiliary class per branch, so the struct
        # that reaches the analyzer is never single-class however the model was
        # written. Stated as a feature so the reason names the construct.
        model = forkjoin(paired=True, closed=True)
        feats = SolverCTMC(model, method='mdd').getMethodFeatureSet('mdd')
        self.assertNotIn('Fork', feats)
        self.assertNotIn('Join', feats)
        self.assertIn('Fork', SolverCTMC.getFeatureSet())
        ok, reason = SolverCTMC(model, method='mdd').supportsModelMethod('mdd')
        self.assertFalse(ok)
        self.assertIn('Fork', reason)

    def test_the_gate_and_the_analyzer_are_one_predicate(self):
        model = forkjoin(paired=False, closed=True)
        _, reason = SolverCTMC(model, method='default').supportsModelMethod('default')
        with self.assertRaises(Exception) as ctx:
            SolverCTMC(model, method='default').getAvgQLen()
        self.assertIn(reason, str(ctx.exception))


# ---------------------------------------------------------------------------
# SolverFLD
# ---------------------------------------------------------------------------

class TestFluidClosedOnlyMethods(unittest.TestCase):

    def test_an_open_model_loses_diffusion_and_tbi(self):
        for method in ('diffusion', 'tbi'):
            ok, reason = SolverFLD(mm1(), method=method).supportsModelMethod(method)
            self.assertFalse(ok, method)
            self.assertIn('OpenClass', reason)

    def test_a_closed_model_keeps_diffusion_and_tbi(self):
        for model in (repairmen(), cqn2()):
            for method in ('diffusion', 'tbi'):
                ok, reason = SolverFLD(model, method=method).supportsModelMethod(method)
                self.assertTrue(ok, '%s: %s' % (method, reason))

    def test_the_qualified_spelling_carries_the_same_envelope(self):
        # 'fluid.tbi' and 'tbi' are one method; a gate that knew only the bare
        # name would offer the qualified one on a model it cannot run.
        for method in ('fluid.diffusion', 'fluid.tbi'):
            ok, _ = SolverFLD(mm1(), method=method).supportsModelMethod(method)
            self.assertFalse(ok, method)

    def test_tbi_refuses_a_cache_station(self):
        feats = SolverFLD(repairmen(), method='tbi').getMethodFeatureSet('tbi')
        self.assertNotIn('Cache', feats)
        self.assertIn('Cache', SolverFLD.getFeatureSet())


class TestFluidHorizonGate(unittest.TestCase):

    def test_a_time_varying_limit_needs_a_finite_timespan(self):
        for method in ('mol', 'mtginf'):
            ok, reason = SolverFLD(mm1(), method=method).supportsModelMethod(method)
            self.assertFalse(ok, method)
            self.assertIn('finite horizon', reason)
            self.assertIn('timespan', reason)

    def test_a_finite_timespan_restores_them(self):
        # The converse, and the whole point of the rule being on the OPTIONS:
        # the model never changed, only the horizon did.
        for method in ('mol', 'mtginf'):
            solver = SolverFLD(mm1(), method=method, timespan=[0.0, 10.0])
            ok, reason = solver.supportsModelMethod(method)
            self.assertTrue(ok, '%s: %s' % (method, reason))

    def test_a_stationary_single_station_limit_is_not_gated_on_a_horizon(self):
        # 'ggisgi' and 'tga' report a stationary point; they are refused on this
        # model for a different reason (no patience law), never on the horizon.
        for method in ('ggisgi', 'tga'):
            ok, reason = SolverFLD(mm1(), method=method).supportsModelMethod(method)
            self.assertFalse(ok, method)
            self.assertNotIn('finite horizon', reason)

    def test_the_gate_and_the_analyzer_are_one_predicate(self):
        # runAnalyzer rather than getAvgQLen: the fluid getters report "runAnalyzer()
        # must complete before accessing results" for ANY gate refusal, so the
        # sentence only survives one call up.
        model = mm1()
        _, reason = SolverFLD(model, method='mtginf').supportsModelMethod('mtginf')
        with self.assertRaises(Exception) as ctx:
            SolverFLD(model, method='mtginf').runAnalyzer()
        self.assertIn('finite horizon', str(ctx.exception))
        self.assertIn('finite horizon', reason)


class TestSsaForkJoinGate(unittest.TestCase):
    """SSA shares the exact fork-join construction, so it shares the gate.

    The other three codebases have always DECLARED Fork/Join/Forker/Joiner in
    the SSA envelope; this port alone omitted them, so the ssa family never
    appeared on a fork-join model it solves. Declaring them without the
    structural gate would have moved the defect rather than fixed it: the wiring
    rules the names cannot state must still refuse.
    """

    def test_the_four_names_are_declared(self):
        feats = SolverSSA.getFeatureSet()
        for name in ('Fork', 'Join', 'Forker', 'Joiner'):
            self.assertIn(name, feats, name)

    def test_a_declared_closed_fork_join_is_offered_and_agrees_with_the_exact_chain(self):
        model = forkjoin(paired=True, closed=True)
        ok, reason = SolverSSA(model, 'default').supportsModelMethod('default')
        self.assertTrue(ok, reason)
        offered = [r.Method for r in model.findSolver().itertuples()
                   if r.Method.startswith('ssa.')]
        self.assertIn('ssa.default', offered)
        # It is offered because it WORKS: the queues match the exact chain to
        # Monte Carlo error. Only the two queue rows are compared -- they are
        # what the simulation and the chain agree on.
        exact = np.asarray(CTMC(model, 'default').getAvgQLen(), dtype=float).ravel()
        got = np.asarray(SSA(model, 'default', samples=20000, seed=23000).getAvgQLen(),
                         dtype=float).ravel()
        for i in (0, 1):
            self.assertAlmostEqual(exact[i], got[i], delta=0.05)

    def test_an_undeclared_pairing_is_refused_in_the_validator_s_own_words(self):
        ok, reason = SolverSSA(forkjoin(paired=False, closed=True),
                               'default').supportsModelMethod('default')
        self.assertFalse(ok)
        self.assertIn('without a matched Join', reason)

    def test_an_open_class_through_a_fork_is_refused(self):
        ok, reason = SolverSSA(forkjoin(paired=True, closed=False),
                               'default').supportsModelMethod('default')
        self.assertFalse(ok)
        self.assertIn('Open classes routed through a Fork', reason)

    def test_a_model_with_no_fork_is_never_asked_about(self):
        for model in (mm1(), repairmen(), cqn2()):
            ok, reason = SolverSSA(model, 'default').supportsModelMethod('default')
            self.assertTrue(ok, reason)


class TestFluidForkJoinGate(unittest.TestCase):
    """A fork-join model is answered by the MMT fixed point, not by one drift."""

    def test_refined_is_closed_only_wherever_the_model_is_open(self):
        # The reference's runAnalyzer has always refused this by name; the
        # featset never said so, and on an open fork-join model the restriction
        # surfaced as a failure inside the transform instead of a refusal.
        for model in (mm1(), forkjoin(paired=True, closed=False)):
            ok, reason = SolverFLD(model, method='refined').supportsModelMethod('refined')
            self.assertFalse(ok)
            self.assertIn('OpenClass', reason)
        feats = SolverFLD(mm1(), method='refined').getMethodFeatureSet('refined')
        self.assertNotIn('OpenClass', feats)
        self.assertIn('OpenClass', SolverFLD.getFeatureSet())

    def test_refined_keeps_every_closed_model_fork_join_included(self):
        for model in (repairmen(), cqn2(), forkjoin(paired=True, closed=True)):
            ok, reason = SolverFLD(model, method='refined').supportsModelMethod('refined')
            self.assertTrue(ok, reason)

    def test_dae_refuses_an_open_fork_join_model_only(self):
        # The conjunction the feature set cannot state: Fork and OpenClass are
        # both declared names, and it is having BOTH that the DAE form cannot
        # take -- the transform's auxiliary open classes have no unknown in it.
        ok, reason = SolverFLD(forkjoin(paired=True, closed=False),
                               method='dae').supportsModelMethod('dae')
        self.assertFalse(ok)
        self.assertIn('fork-join fixed point', reason)
        self.assertIn('minnormal', reason)
        # ... and each half on its own stays runnable.
        for model in (mm1(), forkjoin(paired=True, closed=True)):
            ok, reason = SolverFLD(model, method='dae').supportsModelMethod('dae')
            self.assertTrue(ok, reason)

    def test_diffusion_and_kp_do_not_integrate_a_fork_join_model_at_all(self):
        # Measured, not assumed. On the SYMMETRIC closed fork-join whose exact
        # chain is Q1 = Q2 = 0.664, 'diffusion' put the whole population on ONE
        # station and zero elsewhere (a different station on a rerun), and 'kp'
        # returned an all-zero table on a symmetric OPEN fork-join fed at rate
        # 0.5. Both answered instead of refusing, which is the reason the names
        # come off the envelope: a refusal is recoverable, a silent wrong answer
        # is not. C++ has always withheld them.
        for method in ('diffusion', 'kp'):
            feats = SolverFLD(mm1(), method=method).getMethodFeatureSet(method)
            for name in ('Fork', 'Join', 'Forker', 'Joiner'):
                self.assertNotIn(name, feats, '%s/%s' % (method, name))
        self.assertIn('Fork', SolverFLD.getFeatureSet())
        ok, reason = SolverFLD(forkjoin(paired=True, closed=True),
                              method='diffusion').supportsModelMethod('diffusion')
        self.assertFalse(ok)
        self.assertIn('Fork', reason)

    def test_the_five_that_do_integrate_it_keep_the_fork(self):
        # The converse, and the half that had to come BACK in C++: these five
        # answer the symmetric closed fork-join symmetrically, which is the one
        # property no approximation of a symmetric model may lose.
        model = forkjoin(paired=True, closed=True)
        for method in ('statedep', 'refined', 'tbi', 'mfq', 'rmf'):
            feats = SolverFLD(model, method=method).getMethodFeatureSet(method)
            self.assertIn('Fork', feats, method)
            ok, reason = SolverFLD(model, method=method).supportsModelMethod(method)
            self.assertTrue(ok, '%s: %s' % (method, reason))
            qlen = np.asarray(FLD(model, method).getAvgQLen(), dtype=float).ravel()
            self.assertTrue(np.all(np.isfinite(qlen)), method)
            # Q1 and Q2 are interchangeable in this model.
            self.assertAlmostEqual(qlen[0], qlen[1], places=6, msg=method)

    def test_minnormal_is_the_advice_and_it_runs(self):
        # The refusal names a replacement; a replacement that does not run would
        # be worse than no advice at all.
        model = forkjoin(paired=True, closed=False)
        ok, reason = SolverFLD(model, method='minnormal').supportsModelMethod('minnormal')
        self.assertTrue(ok, reason)
        qlen = np.asarray(FLD(model, 'minnormal').getAvgQLen(), dtype=float)
        self.assertTrue(np.all(np.isfinite(qlen)))

    def test_the_gate_and_the_analyzer_are_one_predicate(self):
        model = forkjoin(paired=True, closed=False)
        _, reason = SolverFLD(model, method='dae').supportsModelMethod('dae')
        with self.assertRaises(Exception) as ctx:
            SolverFLD(model, method='dae').runAnalyzer()
        self.assertIn('fork-join fixed point', str(ctx.exception))
        self.assertIn('fork-join fixed point', reason)


# ---------------------------------------------------------------------------
# the report the gates feed
# ---------------------------------------------------------------------------

class TestFindSolverAgreesWithTheRun(unittest.TestCase):
    """Every ctmc/fluid row findSolver calls runnable must actually run."""

    FAMILIES = {'ctmc': CTMC, 'fluid': FLD, 'ssa': SSA}

    def _check(self, model):
        table = model.findSolver('', False)
        offered = [r.Method for r in table.itertuples()
                   if r.Method.split('.', 1)[0] in self.FAMILIES]
        for token in offered:
            family, method = token.split('.', 1)
            cls = self.FAMILIES[family]
            try:
                kw = dict(samples=2000, seed=23000) if family == 'ssa' else {}
                qlen = np.asarray(cls(model, method, **kw).getAvgQLen(), dtype=float)
            except Exception as err:      # noqa: BLE001 - the failure IS the assertion
                self.fail('%s was offered but raised: %s: %s'
                          % (token, type(err).__name__, err))
            self.assertTrue(np.all(np.isfinite(qlen)), token)
        return offered

    def test_a_single_class_closed_model_keeps_every_ctmc_method(self):
        offered = self._check(repairmen())
        for token in ('ctmc.default', 'ctmc.exact', 'ctmc.gpu', 'ctmc.mdd',
                      'ctmc.cftp', 'ctmc.cftp.approx',
                      'fluid.diffusion', 'fluid.tbi'):
            self.assertIn(token, offered, token)

    def test_a_multiclass_closed_model_drops_only_the_single_class_methods(self):
        offered = self._check(cqn2())
        for gone in ('ctmc.mdd', 'ctmc.cftp', 'ctmc.cftp.approx'):
            self.assertNotIn(gone, offered, gone)
        for kept in ('ctmc.default', 'fluid.diffusion', 'fluid.tbi'):
            self.assertIn(kept, offered, kept)

    def test_an_open_model_drops_the_closed_only_methods(self):
        offered = self._check(mm1())
        for gone in ('ctmc.mdd', 'ctmc.cftp', 'ctmc.cftp.approx',
                     'fluid.diffusion', 'fluid.tbi', 'fluid.mol', 'fluid.mtginf'):
            self.assertNotIn(gone, offered, gone)
        self.assertIn('ctmc.default', offered)

    def test_an_intractable_chain_drops_the_ctmc_family_entirely(self):
        offered = self._check(prio())
        self.assertEqual([t for t in offered if t.startswith('ctmc.')], [])

    def test_an_undeclared_fork_join_pairing_drops_the_ctmc_family_entirely(self):
        offered = self._check(forkjoin(paired=False, closed=False))
        # Both exact families share the construction, so both withdraw.
        self.assertEqual([t for t in offered if t.startswith('ctmc.')], [])
        self.assertEqual([t for t in offered if t.startswith('ssa.')], [])
        # The fluid family still answers it through the MMT fixed point, minus
        # the methods that cannot run it on an open model.
        self.assertIn('fluid.minnormal', offered)
        for gone in ('fluid.refined', 'fluid.dae', 'fluid.kp'):
            self.assertNotIn(gone, offered, gone)

    def test_a_declared_closed_fork_join_keeps_the_exact_answer(self):
        offered = self._check(forkjoin(paired=True, closed=True))
        for kept in ('ctmc.default', 'ctmc.exact', 'ctmc.gpu', 'ssa.default',
                     'fluid.minnormal', 'fluid.refined', 'fluid.tbi', 'fluid.dae'):
            self.assertIn(kept, offered, kept)
        for gone in ('ctmc.mdd', 'fluid.diffusion'):
            self.assertNotIn(gone, offered, gone)


if __name__ == '__main__':
    unittest.main()
