"""The per-method support gate of SolverNC.

WHAT THIS PINS. ``model.help()`` / ``findSolver()`` reports one row per (solver,
method) pair, and it builds the nc rows by asking ``SolverNC.supportsModelMethod``
-- the same gate ``SolverAUTO.chooseSolverRanked`` applies before delegating and
the same one ``listValidMethods`` projects. Until the rules below reached that
gate it was far weaker than what the ANALYZER enforces at run time, so the report
offered pairs that then raised: measured on a plain M/M/1, 15 of the 42 nc.*
rows called runnable raised, and on a single-class closed network 10 did.

TWO ROUTES REACH THE REPORT, and both are asserted here.

  * What the feature registry CAN name rides in ``SolverNC.getMethodFeatureSet``:
    a closed population for 'is' and the six load-dependent evaluators (drop
    OpenClass), no think time for 'divdiff' (drop SchedStrategy_INF). A feature
    set says "I ACCEPT this construct", so it can refuse a model for HAVING one
    and never for lacking one -- which is exactly why the next group cannot live
    there.
  * What it CANNOT -- requires a cache, requires state-dependent routing,
    requires a loss network, requires exactly two stations, requires normal usage
    -- is ``nc_method_refusal``, ONE predicate that ``SolverNC.runAnalyzer``
    raises on and that ``supportsModelMethod`` asks on the way in, so the gate and
    the run are one body of rules rather than two copies that drift.

Each case asserts the refusal AND its converse: a model the method IS derived for
must keep it. Over-tightening hides an answer the user could have had, which is
the same defect with the sign flipped.

The MATLAB twin is matlab/src/solvers/NC/nc_method_refusal.m with
@SolverNC/SolverNC.m; the JAR twin is
jar/src/test/java/jline/solvers/nc/SolverNCGateTest.java and the C++ twin
cpp/tests/test_gate_nc.cpp.
"""

import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np

from line_solver import (Network, Source, Queue, Sink, Delay, OpenClass, ClosedClass,
                         Exp, SchedStrategy, GlobalConstants, VerboseLevel, NC)
from line_solver.solvers.solver_nc.solver_nc import SolverNC

GlobalConstants.setVerbose(VerboseLevel.SILENT)


def mm1():
    """Source -> FCFS Queue -> Sink, one OPEN class."""
    m = Network('ncGateMM1')
    s, q, k = Source(m, 'S'), Queue(m, 'Q', SchedStrategy.FCFS), Sink(m, 'K')
    c = OpenClass(m, 'C')
    s.setArrival(c, Exp(1.0))
    q.setService(c, Exp(2.0))
    m.link(Network.serialRouting(s, q, k))
    return m


def repairmen():
    """Delay -> FCFS Queue, one closed class, N = 3. Saturated: the queue's
    offered load 1.5 exceeds its unit rate, so the model is NOT in normal
    usage and PANACEA's expansion does not apply to it."""
    m = Network('ncGateRepairmen')
    d, q = Delay(m, 'D'), Queue(m, 'Q', SchedStrategy.FCFS)
    c = ClosedClass(m, 'C', 3, d)
    d.setService(c, Exp(1.0))
    q.setService(c, Exp(2.0))
    m.link(Network.serialRouting(d, q))
    return m


def multiserver():
    """The same closed shape with a three-server queue, N = 4: mu(n) = min(n,3)
    at saturation is 3 against an offered load of 2, so this one IS in normal
    usage and 'panald' must survive the gate."""
    m = Network('ncGateMultiserver')
    d, q = Delay(m, 'D'), Queue(m, 'Q', SchedStrategy.FCFS)
    q.setNumberOfServers(3)
    c = ClosedClass(m, 'C', 4, d)
    d.setService(c, Exp(1.0))
    q.setService(c, Exp(2.0))
    m.link(Network.serialRouting(d, q))
    return m


def loaddep():
    """The repairmen shape with a rate LATTICE on the queue, N = 4. It matters
    because pfqn_ncld evaluates 'pana' and 'panald' with the same
    pfqn_panaceald, so on a load-dependent model the load-INDEPENDENT name
    reaches the load-dependent expansion. Saturation rate 1.9 against an offered
    load of 2, so this one is NOT in normal usage and both names must be
    refused; every other nc method on it runs."""
    m = Network('ncGateLoadDep')
    d, q = Delay(m, 'D'), Queue(m, 'Q', SchedStrategy.FCFS)
    c = ClosedClass(m, 'C', 4, d)
    d.setService(c, Exp(1.0))
    q.setService(c, Exp(2.0))
    q.setLoadDependence(np.array([1.0, 1.6, 1.8, 1.9]))
    m.link(Network.serialRouting(d, q))
    return m


def cqn3():
    """Delay -> Q1 -> Q2, one closed class, N = 4: TWO queueing stations, which
    is what the mmint2/gleint/comomld recursions are not stated for."""
    m = Network('ncGateCqn3')
    d = Delay(m, 'D')
    q1, q2 = Queue(m, 'Q1', SchedStrategy.PS), Queue(m, 'Q2', SchedStrategy.PS)
    c = ClosedClass(m, 'C', 4, d)
    d.setService(c, Exp(1.0))
    q1.setService(c, Exp(2.0))
    q2.setService(c, Exp(3.0))
    m.link(Network.serialRouting(d, q1, q2))
    return m


def open_tandem():
    """Source -> Q1 -> Q2 -> Sink: two queueing stations but NO closed
    population, so pfqn_nc answers with the exact open formulas before its
    method switch and the single-station rule must stay inactive."""
    m = Network('ncGateOpenTandem')
    s = Source(m, 'S')
    q1, q2 = Queue(m, 'Q1', SchedStrategy.FCFS), Queue(m, 'Q2', SchedStrategy.FCFS)
    k = Sink(m, 'K')
    c = OpenClass(m, 'C')
    s.setArrival(c, Exp(0.5))
    q1.setService(c, Exp(2.0))
    q2.setService(c, Exp(3.0))
    m.link(Network.serialRouting(s, q1, q2, k))
    return m


def cyclic_delay_free():
    """Queue1 -> Queue2, one closed class, NO delay station: zero think time,
    which is the shape the divided-difference closed form 'divdiff' is derived
    for (Casale, SIGMETRICS 2017, Eqs. 15-16)."""
    m = Network('ncGateCyclic')
    q1, q2 = Queue(m, 'Q1', SchedStrategy.PS), Queue(m, 'Q2', SchedStrategy.PS)
    c = ClosedClass(m, 'C', 2, q1)
    q1.setService(c, Exp(1.0))
    q2.setService(c, Exp(2.0))
    m.link(Network.serialRouting(q1, q2))
    return m


def gate(model, method):
    """(runnable, reason) from the gate the report and the ranking both ask."""
    solver = SolverNC(model, method)
    ok, reason = solver.supportsModelMethod(method)
    return bool(ok), (reason or '')


def run_question(model, method):
    """The predicate asked the RUN's question: '' when the reference performs the
    method by name, even where the report declines to offer it."""
    from line_solver.solvers.solver_nc.nc_method_refusal import nc_method_refusal
    solver = SolverNC(model, method)
    return nc_method_refusal(model.getStruct(), method, solver.options, for_report=False)


def raises(model, method):
    """Did the ANALYZER refuse the same pair? The message when it did, else ''."""
    try:
        NC(model, method).avgTable()
    except Exception as e:      # noqa: BLE001 - the message is the assertion
        return str(e)
    return ''


class TestNCGateStructural(unittest.TestCase):
    """What the feature registry cannot name: nc_method_refusal."""

    def test_the_cache_tokens_name_the_missing_cache(self):
        for method in ('rayint', 'spm'):
            ok, reason = gate(mm1(), method)
            self.assertFalse(ok, method)
            self.assertIn('Cache node', reason)
            self.assertIn('SPM saddle point', reason)

    def test_the_loss_network_tokens_name_the_missing_region(self):
        for method in ('ms', 'erlangfp'):
            ok, reason = gate(repairmen(), method)
            self.assertFalse(ok, method)
            self.assertIn('loss network', reason)

    def test_rec_names_both_routes_it_has(self):
        ok, reason = gate(repairmen(), 'rec')
        self.assertFalse(ok)
        self.assertIn('MDD-rec', reason)
        self.assertIn('Petri net', reason)

    def test_sdr_names_the_routing_the_model_does_not_declare(self):
        for method in ('sdr', 'sdr.mva'):
            ok, reason = gate(mm1(), method)
            self.assertFalse(ok, method)
            self.assertIn('state-dependent routing', reason)

    def test_morrison_names_the_shape_it_is_derived_for(self):
        ok, reason = gate(repairmen(), 'morrison')
        self.assertFalse(ok)
        self.assertIn('DPS', reason)
        self.assertIn('two stations', reason)

    def test_panaceald_names_normal_usage_on_a_saturated_closed_model(self):
        ok, reason = gate(repairmen(), 'panald')
        self.assertFalse(ok)
        self.assertIn('normal usage', reason)

    def test_the_single_station_recursions_name_the_station_count(self):
        # pfqn_nc states 'a model with a delay and a single queueing station' for
        # mmint2/gleint, and pfqn_comomrm_ld raises 'accepts at most a single
        # queueing station'. Neither is a feature name: it is a COUNT.
        model = cqn3()
        for method in ('mmint2', 'gleint', 'comomld'):
            ok, reason = gate(model, method)
            self.assertFalse(ok, method)
            self.assertIn('single queueing station', reason, method)
            self.assertIn('has 2', reason, method)

    def test_mmint2_and_gleint_are_gated_for_the_report_only(self):
        """THE RULING (2026-07-25, reaffirmed when this gate was added): the report
        answers 'should this be offered' and the run answers 'what does the
        reference do'. pfqn_nc answers mmint2/gleint outside their shape with an
        empty constant and a ZERO TABLE (pfqn_nc.m case {'mmint2','gleint'}:
        lG = [] and return, unconditionally), so a caller who names the method
        keeps that answer while model.help() stops offering it.

        comomld is NOT in that bucket: pfqn_comomrm_ld raises natively, so it is
        refused on both paths."""
        model = cqn3()
        for method in ('mmint2', 'gleint'):
            ok, _ = gate(model, method)
            self.assertFalse(ok, method)                       # the report declines
            self.assertEqual('', run_question(model, method),  # the run does not
                             '%s must stay runnable when asked by name' % method)
        ok, _ = gate(model, 'comomld')
        self.assertFalse(ok)
        self.assertNotEqual('', run_question(model, 'comomld'),
                            'comomld raises natively and must be refused on both paths')

    def test_mmint2_and_gleint_return_the_references_zero_table(self):
        """The RUN half of the ruling, asserted NUMERICALLY.

        pfqn_nc.m:287 warns, sets lG = [] and returns on more than one queueing
        station, and getAvg then renders a table of ZEROS. This port used to have
        no such guard: quadrature on a multi-station L returned a PLAUSIBLE WRONG
        NUMBER for 'mmint2' (QLen 2.9037/0.6578/0.4385 against the exact
        1.5548/1.6109/0.8343) and a numpy shape error for 'gleint'. A wrong number
        under the method's own name is the worst of the four codebases' answers,
        which is why it was closed to the reference's zero table rather than to a
        refusal -- the run path is ruled, and the report already declines the pair.
        """
        model = cqn3()
        exact = np.asarray(NC(model, 'exact').getAvgQLen(), dtype=float).flatten()
        # the model itself is not degenerate: the zeros below are the decline
        self.assertTrue(np.all(exact > 0), exact)
        for method in ('mmint2', 'gleint'):
            qn = np.asarray(NC(cqn3(), method).getAvgQLen(), dtype=float).flatten()
            self.assertEqual(exact.shape, qn.shape, method)   # same shape as MATLAB's
            self.assertTrue(np.all(qn == 0), '%s: %s' % (method, qn))

    def test_mmint2_and_gleint_keep_their_answer_on_their_own_shape(self):
        # The converse: one queueing station, or an open model that never reaches
        # pfqn_nc's method switch, must be untouched by the guard.
        one_queue = np.asarray(NC(repairmen(), 'mmint2').getAvgQLen(), dtype=float).flatten()
        self.assertTrue(np.any(one_queue > 0), one_queue)
        open_exact = np.asarray(NC(open_tandem(), 'exact').getAvgQLen(), dtype=float).flatten()
        for method in ('mmint2', 'gleint'):
            got = np.asarray(NC(open_tandem(), method).getAvgQLen(), dtype=float).flatten()
            self.assertTrue(np.allclose(got, open_exact), '%s: %s' % (method, got))

    def test_the_single_station_recursions_survive_on_one_queueing_station(self):
        model = repairmen()
        for method in ('mmint2', 'gleint', 'comomld'):
            ok, reason = gate(model, method)
            self.assertTrue(ok, '%s: %s' % (method, reason))
            self.assertEqual('', raises(repairmen(), method), method)

    def test_the_station_count_rule_is_inactive_without_a_closed_population(self):
        # An open network never reaches pfqn_nc's method switch, so a two-queue
        # OPEN tandem runs these names correctly and must keep them. Gating on the
        # station count alone would have taken three working rows away.
        model = open_tandem()
        for method in ('mmint2', 'gleint'):
            ok, reason = gate(model, method)
            self.assertTrue(ok, '%s: %s' % (method, reason))
            self.assertEqual('', raises(open_tandem(), method), method)

    def test_panacea_is_refused_on_a_load_dependent_model_outside_normal_usage(self):
        # pfqn_ncld's case label is ('pana', 'panald'): on a rate lattice
        # the load-independent NAME is evaluated by the load-dependent kernel, so
        # it raises the panald refusal. The gate has to know that aliasing.
        ok, reason = gate(loaddep(), 'pana')
        self.assertFalse(ok)
        self.assertIn('normal usage', reason)
        self.assertIn("evaluates it as 'panald'", reason)
        self.assertNotEqual('', raises(loaddep(), 'pana'))

    def test_the_load_dependent_model_keeps_every_other_method(self):
        # The converse, and the whole point of not gating 'pana' off the
        # lattice: 31 of the 32 nc rows on this model run, and must go on running.
        model = loaddep()
        for method in ('default', 'exact', 'ca', 'clw', 'comom', 'comomld', 'le',
                       'ble', 'mmint2', 'gleint', 'kt', 'bkt', 'lekt', 'cub', 'rd', 'nrl',
                       'nrp', 'nre', 'propfair', 'ger', 'rgf'):
            ok, reason = gate(model, method)
            self.assertTrue(ok, '%s: %s' % (method, reason))

    def test_panacea_is_left_alone_off_the_load_dependent_route(self):
        # Without a rate lattice 'pana' takes its own pfqn_nc arm, which warns
        # and returns an empty constant rather than raising, so the gate must not
        # refuse it there -- that would be over-tightening on a path this change
        # is not about.
        ok, reason = gate(repairmen(), 'pana')
        self.assertTrue(ok, reason)
        ok, reason = gate(multiserver(), 'pana')
        self.assertTrue(ok, reason)

    def test_panaceald_survives_where_the_expansion_does_apply(self):
        # The converse: a multiserver station whose saturation rate 3 exceeds the
        # offered load 2 IS in normal usage, so the row must stay.
        ok, reason = gate(multiserver(), 'panald')
        self.assertTrue(ok, reason)
        self.assertEqual('', raises(multiserver(), 'panald'))


class TestNCGateFeatureSet(unittest.TestCase):
    """What the registry CAN name: SolverNC.getMethodFeatureSet."""

    def test_explicit_refuses_a_think_time_by_feature_name(self):
        ok, reason = gate(repairmen(), 'divdiff')
        self.assertFalse(ok)
        self.assertIn('SchedStrategy_INF', reason)

    def test_explicit_survives_on_a_delay_free_closed_model(self):
        model = cyclic_delay_free()
        ok, reason = gate(model, 'divdiff')
        self.assertTrue(ok, reason)
        self.assertEqual('', raises(model, 'divdiff'))
        self.assertNotIn('SchedStrategy_INF',
                         SolverNC(model, 'divdiff').getMethodFeatureSet('divdiff'))
        self.assertIn('SchedStrategy_INF',
                      SolverNC(model, 'default').getMethodFeatureSet('default'))

    def test_the_closed_population_methods_refuse_an_open_chain_by_feature_name(self):
        for method in ('is', 'rd', 'nrp', 'nrl', 'nre', 'comomld', 'panald'):
            ok, reason = gate(mm1(), method)
            self.assertFalse(ok, method)
            self.assertIn('OpenClass', reason, method)

    def test_the_closed_population_methods_survive_on_a_closed_model(self):
        model = repairmen()
        for method in ('is', 'rd', 'nrp', 'nrl', 'nre', 'comomld'):
            ok, reason = gate(model, method)
            self.assertTrue(ok, '%s: %s' % (method, reason))
            self.assertEqual('', raises(repairmen(), method), method)


class TestNCGateAgreesWithTheRun(unittest.TestCase):
    """ONE PREDICATE, TWO CALLERS: every pair the report offers must run, and
    every pair it refuses must be one the analyzer also refuses."""

    MODELS = {'mm1': mm1, 'repairmen': repairmen, 'cyclic': cyclic_delay_free,
              'loaddep': loaddep, 'cqn3': cqn3, 'open_tandem': open_tandem}

    def test_every_refused_pair_is_refused_by_the_analyzer_too(self):
        for name, build in self.MODELS.items():
            declared = SolverNC(build(), 'default').listValidMethods()
            for method in declared:
                ok, reason = gate(build(), method)
                if ok:
                    continue
                if not run_question(build(), method):
                    # A REPORT-ONLY refusal, and the ruling says so: the reference
                    # performs this one by name (a warning and a zero table), which
                    # is exactly why the report declines to offer it.
                    continue
                msg = raises(build(), method)
                self.assertNotEqual('', msg, '%s/%s: the gate refused but the run '
                                             'did not' % (name, method))

    def test_the_model_keeps_the_methods_it_can_genuinely_run(self):
        # A closed product-form network must not lose its normalizing-constant
        # methods to this change; these are the ones it has always answered with.
        model = repairmen()
        for method in ('default', 'exact', 'ca', 'comom', 'le', 'ble', 'mmint2',
                       'gleint', 'pana', 'propfair', 'cub', 'kt', 'bkt', 'lekt', 'clw'):
            ok, reason = gate(model, method)
            self.assertTrue(ok, '%s: %s' % (method, reason))


if __name__ == '__main__':
    unittest.main()
