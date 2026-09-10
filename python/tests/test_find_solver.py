"""
Tests of findSolver: the report of which solvers and solver methods can analyze
a model, why the others cannot, what kind of answer each returns and which
measures it can report.

WHAT IS ASSERTED, and why none of it is a value read back out of the
implementation.

 1. THE PROJECTION IDENTITY. SolverAUTO.listValidMethods is defined as the
    runnable rows of findSolver plus the method names that name no single method.
    That is the point of the refactor -- one gate, two views -- so it is
    asserted directly rather than trusted.

 2. STRUCTURAL FACTS ABOUT THE MODELS. An M/M/1 has a product-form solution and
    is one queueing station fed by a Source, so exact MVA and the QBD are exact
    ON IT; a three-station closed network is still product form but is no
    longer one queue, so the QBD is not. These are properties of the models,
    decided before the code was written.

 3. THE REFUSAL REASONS NAME THE FEATURE. A report whose whole content is the
    explanation is worthless if the explanation is "unsupported".
"""

import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from line_solver import (Network, Source, Queue, Sink, Delay, Cache, OpenClass, ClosedClass,
                         Exp, Pareto, Zipf, Immediate, ReplacementStrategy, SchedStrategy)
from line_solver.solvers.solver_auto.solver_auto import SolverAUTO


def mm1():
    """Source -> FCFS Queue -> Sink, exponential throughout."""
    m = Network('mm1')
    s, q, k = Source(m, 'Source'), Queue(m, 'Queue', SchedStrategy.FCFS), Sink(m, 'Sink')
    c = OpenClass(m, 'C1')
    s.setArrival(c, Exp(1.0))
    q.setService(c, Exp(2.0))
    m.link(Network.serialRouting(s, q, k))
    return m


def mpareto():
    """The same shape with a Pareto service law: SolverNC has no Pareto."""
    m = Network('mpareto')
    s, q, k = Source(m, 'Source'), Queue(m, 'Queue', SchedStrategy.FCFS), Sink(m, 'Sink')
    c = OpenClass(m, 'C1')
    s.setArrival(c, Exp(0.5))
    q.setService(c, Pareto(2.5, 1.0))
    m.link(Network.serialRouting(s, q, k))
    return m


def cqn():
    """Delay -> Queue1 -> Queue2, N = 4: a closed product-form network."""
    m = Network('cqn')
    d = Delay(m, 'Delay')
    q1 = Queue(m, 'Queue1', SchedStrategy.FCFS)
    q2 = Queue(m, 'Queue2', SchedStrategy.FCFS)
    c = ClosedClass(m, 'C1', 4, d)
    d.setService(c, Exp(1.0))
    q1.setService(c, Exp(2.0))
    q2.setService(c, Exp(3.0))
    m.link(Network.serialRouting(d, q1, q2))
    return m


def cachemodel():
    """Client -> LRU Cache -> CacheDelay, hit and miss classes.

    Three closed classes and a cache, which is the shape that separates a
    solver's own envelope from a neighbour's: SolverMVA, SolverCTMC, SolverSSA
    and SolverLDES declare Cache, and SolverBA, SolverMAM and SolverNC do not.
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


def row(table, method):
    hit = table[table.Method == method]
    return None if hit.empty else hit.iloc[0]


class TestFindSolverProjection(unittest.TestCase):

    def test_list_valid_methods_is_the_runnable_rows(self):
        m = mm1()
        auto = SolverAUTO(m, verbose=False)
        all_rows = auto.findSolver('', True)
        valid = auto.listValidMethods()
        self.assertGreater(len(all_rows), 0)

        runnable = 0
        for _, r in all_rows.iterrows():
            if r['Runnable']:
                runnable += 1
                # Every runnable pair is a method name a caller may ask AUTO for.
                self.assertIn(r['Method'], valid)
                self.assertIn(r['Solver'], valid)
                self.assertEqual(r['Reason'], '')
            else:
                # and a refused one is not offered.
                self.assertNotIn(r['Method'], valid)
                self.assertNotEqual(r['Reason'], '')
        self.assertGreater(runnable, 0)

        # The default report is exactly the runnable half.
        only_runnable = auto.findSolver()
        self.assertEqual(len(only_runnable), runnable)
        self.assertTrue(only_runnable.Runnable.all())

        # The selection intents name a ranking rather than an algorithm and are
        # listed whatever the model is.
        for intent in ('default', 'exact', 'bound'):
            self.assertIn(intent, valid)

    def test_a_self_qualified_spelling_is_not_doubled(self):
        # The fluid registry declares both 'dae' and 'fluid.dae' so that its own
        # gate takes either; prefixing the family again would yield
        # 'fluid.fluid.dae', a method name that resolves but names the same method
        # twice and would double every fluid row.
        t = mm1().findSolver('', True)
        self.assertIsNotNone(row(t, 'fluid.dae'))
        self.assertIsNone(row(t, 'fluid.fluid.dae'))
        for _, r in t.iterrows():
            self.assertNotIn('%s.%s.' % (r['Solver'], r['Solver']), r['Method'])

    def test_model_entry_points_and_their_aliases(self):
        m = mm1()
        t = m.findSolver()
        self.assertGreater(len(t), 0)
        self.assertEqual(list(t.columns),
                         ['Solver', 'Method', 'Runnable', 'Class', 'Metrics', 'Reason'])
        # findMethod and help are the same question asked in other words.
        self.assertTrue(t.equals(m.findMethod()))
        self.assertTrue(t.equals(m.help()))
        self.assertTrue(t.equals(m.find_solver()))

    def test_a_report_prints_nothing(self):
        # Asking a solver whether it supports the model warns naming the missing
        # feature, which is right on the solve path and wrong here: the answer
        # IS the table.
        import contextlib
        import io
        m = mpareto()
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
            m.findSolver('', True)
        self.assertEqual(buf.getvalue(), '')


class TestFindSolverClassification(unittest.TestCase):

    def test_exactness_is_claimed_of_the_model_not_the_algorithm(self):
        open_model = mm1()
        self.assertTrue(open_model.has_product_form_solution())
        a = open_model.findSolver('', True)
        # Exact MVA on a product-form model; the QBD on one queueing station
        # fed by a Source.
        self.assertEqual(row(a, 'mva.exact')['Class'], 'exact')
        self.assertEqual(row(a, 'mam.default')['Class'], 'exact')
        # A decomposition of a network into such queues approximates it.
        self.assertEqual(row(a, 'mam.dec.source')['Class'], 'approx')

        # Three stations: still product form, so MVA stays exact, but no longer
        # one queueing station, so the QBD is not.
        closed = cqn()
        b = closed.findSolver('', True)
        self.assertEqual(row(b, 'mva.exact')['Class'], 'exact')
        mam = row(b, 'mam.default')
        if mam is not None:
            self.assertEqual(mam['Class'], 'approx')
        # Every AMVA arm is an approximation whatever the model.
        self.assertEqual(row(b, 'mva.amva')['Class'], 'approx')
        # Bounds are what SolverBA is for, and a simulator is a simulator.
        for _, r in b.iterrows():
            if r['Solver'] == 'ba':
                self.assertEqual(r['Class'], 'bound')
            if r['Solver'] in ('ssa', 'ldes'):
                self.assertEqual(r['Class'], 'simulation')
        # The CTMC generator is solved as written; 'cftp.approx' says otherwise
        # in its own name.
        cftp = row(b, 'ctmc.cftp.approx')
        if cftp is not None:
            self.assertEqual(cftp['Class'], 'approx')
        self.assertEqual(row(b, 'ctmc.exact')['Class'], 'exact')

    def test_a_refusal_names_the_feature_that_caused_it(self):
        t = mpareto().findSolver('', True)
        refused = t[~t.Runnable]
        self.assertGreater(len(refused), 0)
        self.assertTrue((refused.Reason != '').all())
        # SolverNC has no Pareto in its feature set, so it must refuse and say so.
        nc = row(t, 'nc.default')
        self.assertIsNotNone(nc)
        self.assertFalse(nc['Runnable'])
        self.assertIn('Pareto', nc['Reason'])

    def test_a_cache_withdraws_the_exact_claim_from_the_analytic_families(self):
        # sn_has_product_form answers about the QUEUEING network and knows
        # nothing of a cache, so a cache model reads as product form and
        # 'mva.exact' was labelled exact on it. It is not: on this model exact
        # MVA returns 0.2516 at the hit station where the CTMC returns 0.3022
        # and simulation 0.3023.
        m = cachemodel()
        self.assertTrue(m.hasProductFormSolution())
        t = m.findSolver()
        for method in ('mva.exact', 'mva.mva'):
            r = row(t, method)
            if r is not None:
                self.assertEqual('approx', r['Class'],
                                 '%s claims exactness a cache model cannot give' % method)
        # SolverCTMC is NOT conditioned on it: its state space carries the cache
        # contents, so it really is exact there.
        self.assertEqual('exact', row(t, 'ctmc.exact')['Class'])
        # And the condition is about the cache, not about exactness in general:
        # the same claim stands on a cache-free product-form network.
        self.assertEqual('exact', row(cqn().findSolver(), 'mva.exact')['Class'])

    def test_the_single_station_fluid_limits_need_the_shape_they_are_stated_for(self):
        # 'ggingi.tga' is a limit for a queue customers ABANDON; a plain M/M/1
        # has no patience law, and offering the method there invited a caller to
        # ask for an analysis that then raised.
        t = mm1().findSolver('', True)
        tga = row(t, 'fluid.ggingi.tga')
        self.assertIsNotNone(tga)
        self.assertFalse(tga['Runnable'])
        self.assertIn('patience', tga['Reason'])
        # 'mtginf' is stated for the same single-station shape but not for
        # abandonment, so the SHAPE rule must not be what withdraws it: it is
        # withdrawn only for want of a finite horizon, which is an option and
        # not a property of the model. Set options.timespan and it comes back.
        mtginf = row(t, 'fluid.mtginf')
        self.assertFalse(mtginf['Runnable'])
        self.assertNotIn('patience', mtginf['Reason'])
        self.assertIn('timespan', mtginf['Reason'])


class TestFindSolverMetrics(unittest.TestCase):

    def test_a_metric_narrows_the_report_to_the_families_that_answer_it(self):
        m = mm1()
        by_group = m.findSolver('cdf')
        by_accessor = m.findSolver('getCdfRespT')
        self.assertGreater(len(by_group), 0)
        # A group name and the accessor that returns it ask the same question.
        self.assertTrue(by_group.equals(by_accessor))
        # MVA computes no passage-time law and must not appear; the simulators
        # and the transform solvers do.
        self.assertNotIn('mva', set(by_group.Solver))
        self.assertIn('ldes', set(by_group.Solver))
        self.assertTrue(by_group.Metrics.str.contains('cdf').all())
        # Every family answers the mean measures, so 'avg' narrows nothing away.
        self.assertEqual(len(m.findSolver('avg')), len(m.findSolver()))

    def test_an_unknown_measure_is_a_caller_error(self):
        # Not an empty answer, which would read as "nothing can do this".
        with self.assertRaises(ValueError):
            mm1().findSolver('nosuchmeasure')

    def test_the_metric_registry_is_self_consistent(self):
        groups = SolverAUTO.metricGroups()
        # A group name maps to itself.
        for g in groups:
            self.assertEqual(SolverAUTO.metricGroupOf(g), g)
        # Every group a family declares is a registered one; a typo here would
        # silently hide the family from a caller asking for that measure.
        for fam in SolverAUTO.familyNames():
            declared = SolverAUTO.familyMetrics(fam)
            for g in declared:
                self.assertIn(g, groups)
            # Every family answers the mean measures, which is what a solver is for.
            self.assertIn('avg', declared)
        self.assertEqual(SolverAUTO.metricGroupOf(''), '')
        self.assertEqual(SolverAUTO.metricGroupOf('any'), '')
        self.assertEqual(SolverAUTO.metricGroupOf('getAvgTable'), 'avg')
        self.assertEqual(SolverAUTO.metricGroupOf('getTranProbAggr'), 'tranprob')


class TestOneRowPerAlgorithm(unittest.TestCase):
    """A report row is a (family, ALGORITHM) pair, not a (family, spelling) one."""

    def test_an_alias_spelling_is_not_a_second_row(self):
        # SolverMVA advertises every AMVA name twice, plain and 'amva.'-prefixed,
        # and strips the prefix before selecting an algorithm. Both spellings
        # were listed, which was 20 of the 49 mva rows on a cache model and 21 of
        # 52 on a closed one: the same algorithm, offered twice, under two names
        # that cannot give different answers.
        offered = list(cqn().findSolver().query("Solver == 'mva'").Method)
        self.assertTrue(offered)
        for method in offered:
            self.assertFalse(method.startswith('mva.amva.'),
                             '%s respells a method already listed' % method)
        # The plain spelling is the one kept, and it is still there.
        self.assertIn('mva.lin', offered)

    def test_the_alias_rule_needs_the_method_it_aliases(self):
        # An alias only when the thing it aliases is declared beside it, so a
        # genuine method that merely starts with the prefix is not eaten.
        self.assertTrue(SolverAUTO.isMethodAlias('mva', 'amva.lin', ['lin', 'amva.lin']))
        self.assertFalse(SolverAUTO.isMethodAlias('mva', 'amva.lin', ['amva.lin']))
        self.assertFalse(SolverAUTO.isMethodAlias('nc', 'sdr.mva', ['sdr.mva', 'mva']))

    def test_both_spellings_still_reach_the_solver(self):
        # The report drops the duplicate ROW, not the method name: a caller carrying
        # the alias must still be able to pass it.
        from line_solver.solvers.solver_mva import SolverMVA
        declared = SolverMVA(cqn()).listValidMethods()
        self.assertIn('lin', declared)
        self.assertIn('amva.lin', declared)


class TestTheGateIsTheSolversOwn(unittest.TestCase):
    """A family must be judged by ITS OWN envelope, never by another's.

    Both cases below were live defects that this codebase alone had, and both
    showed up as families offered on a cache model that every one of their
    methods refuses at run time. MATLAB, the JAR and C++ reported 78, 76 and 76
    runnable pairs on the model below; native python reported 117.
    """

    def test_ba_is_gated_by_its_own_feature_set_not_by_mvas(self):
        # SolverBA subclasses SolverMVA for its model parsing, and inherited the
        # MVA method gate with it. SolverMVA declares Cache and the replacement
        # strategies; SolverBA deliberately does not, because the bounds are
        # derived for a product-form network and the analyzer has no
        # representation of a cache.
        from line_solver.solvers.solver_ba.solver_ba import SolverBA
        m = cachemodel()
        probe = SolverBA(m)
        # A method's set is drawn from SolverBA's OWN envelope. It is a subset
        # rather than an equality: the envelope is the union over the family, so
        # a method that consumes less than the union (gb.upper, which 'default'
        # resolves to, takes no MAP where the spnlp and MAP-AMVA arms do) sits
        # strictly inside it. That is what per-method gating is for.
        self.assertTrue(set(probe.getMethodFeatureSet('default'))
                        <= set(SolverBA.getFeatureSet()))
        # The actual defect this guards: SolverMVA's names must not leak in
        # through the inherited gate.
        for mva_only in ('Cache', 'CacheClassSwitcher', 'ReplacementStrategy_LRU'):
            self.assertNotIn(mva_only, probe.getMethodFeatureSet('default'))
            self.assertNotIn(mva_only, SolverBA.getFeatureSet())
        for method in probe.listValidMethods():
            ok, reason = probe.supportsModelMethod(method)
            self.assertFalse(ok, 'ba.%s was offered on a cache model' % method)
            self.assertIn('Cache', reason)
        self.assertEqual([], list(m.findSolver().query("Solver == 'ba'").Method))

    def test_mams_structural_predicate_is_not_the_whole_gate(self):
        # supports_network answers about topology and class mix, so a
        # decomposition method sailed through it on a model built from
        # constructs SolverMAM does not declare. The feature-set gate has to run
        # as well, which is how MATLAB ends the same method.
        self.assertEqual([], list(cachemodel().findSolver().query("Solver == 'mam'").Method))

    def test_a_flat_network_family_is_not_offered_for_an_environment(self):
        # familyAcceptsModelClass named the model class each COMPOSITE family
        # takes and said nothing about what the flat-network families take, so
        # every one of them answered yes for an Environment. SolverQNS does not
        # refuse such a model either -- its constructor leaves self.model None
        # rather than raising -- so the gate was asked about nothing and said
        # yes, and seven qns rows were half the report on tut12.
        from line_solver import Environment
        base = Network('base')
        d, q = Delay(base, 'Think'), Queue(base, 'Srv', SchedStrategy.FCFS)
        c = ClosedClass(base, 'Jobs', 3, d)
        d.setService(c, Exp(1.0))
        q.setService(c, Exp(2.0))
        base.link(Network.serialRouting(d, q))
        env = Environment('Modes', 2)
        for stage, (name, kind, rate) in enumerate((('Fast', 'operational', 4.0),
                                                    ('Slow', 'degraded', 1.0))):
            copy = base.copy()
            copy.getNodeByName('Srv').setService(copy.classes[0], Exp(rate))
            env.addStage(stage, name, kind, copy)
        env.addTransition(0, 1, Exp(0.5))
        env.addTransition(1, 0, Exp(1.0))
        offered = env.findSolver()
        self.assertEqual({'env'}, set(offered.Solver),
                         'a family that does not solve an Environment was offered for one')
        self.assertTrue(len(offered) > 0)

    def test_the_bounds_still_answer_the_models_they_are_derived_for(self):
        # The fix must not cost the family the models it exists for: a
        # single-class closed product-form network is one.
        offered = list(cqn().findSolver().query("Solver == 'ba'").Method)
        self.assertIn('ba.aba.lower', offered)
        self.assertIn('ba.aba.upper', offered)


if __name__ == '__main__':
    unittest.main()
