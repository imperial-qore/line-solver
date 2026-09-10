"""The SolverBA method gate: what findSolver offers must be what runs.

model.help() / model.findSolver() reports one row per (solver, method) pair a
model can run, and it builds that report from SolverBA.listValidMethods filtered
by SolverBA.supportsModelMethod. Both were much weaker than the rules the
ANALYZER enforces, so the report offered pairs that raised on contact: measured
on a two-class closed network, 30 of the 36 ba.* rows findSolver called runnable
raised, and on a closed multiserver model 31 of 38 did.

The fix is one predicate, ba_method_refusal, asked by the analyzer, by the gate
and by the list, plus the per-method feature-set deltas for the premises a
feature NAME can state. These tests pin both directions: a refused pair names
the offending thing, and a model the bounds ARE derived for keeps every method.
"""

import unittest

import numpy as np

from line_solver import (BA, ClosedClass, Delay, Erlang, Exp, GlobalConstants,
                         Network, OpenClass, Place, Queue, SchedStrategy, Sink,
                         Source, Transition, VerboseLevel)
from line_solver.solvers.solver_ba.solver_ba import SolverBA as NativeSolverBA
from line_solver.solvers.solver_ba.solver_ba_analyzer import (
    ba_method_degenerate, ba_method_refusal, ba_resolve_method,
    solver_ba_analyzer)

GlobalConstants.setVerbose(VerboseLevel.SILENT)


def cqn_delay_free():
    """Single-class closed, no delay, one server each: the shape every bound
    family in the solver is derived for."""
    model = Network('baGateCyclic')
    q1 = Queue(model, 'Q1', SchedStrategy.FCFS)
    q2 = Queue(model, 'Q2', SchedStrategy.PS)
    c = ClosedClass(model, 'C', 3, q1)
    q1.setService(c, Exp(1.0))
    q2.setService(c, Exp(2.0))
    model.link(Network.serialRouting(q1, q2))
    return model


def cqn_two_class():
    """Two closed classes, and DELAY-FREE so that the class premise is the only
    one this model violates: the single-class families have no demand vector."""
    model = Network('baGate2Class')
    q1 = Queue(model, 'Q1', SchedStrategy.PS)
    q2 = Queue(model, 'Q2', SchedStrategy.PS)
    c1 = ClosedClass(model, 'C1', 2, q1)
    c2 = ClosedClass(model, 'C2', 1, q1)
    q1.setService(c1, Exp(1.0))
    q1.setService(c2, Exp(2.0))
    q2.setService(c1, Exp(2.0))
    q2.setService(c2, Exp(3.0))
    model.link(Network.serialRouting(q1, q2))
    return model


def cqn_multiserver():
    """Single-class closed with a three-server station: only ssd, ldbcmp and the
    auto composite that draws on them survive."""
    model = Network('baGateMultiserver')
    d = Delay(model, 'D')
    q = Queue(model, 'Q', SchedStrategy.FCFS)
    q.setNumberOfServers(3)
    c = ClosedClass(model, 'C', 4, d)
    d.setService(c, Exp(1.0))
    q.setService(c, Exp(2.0))
    model.link(Network.serialRouting(d, q))
    return model


def cqn_with_delay():
    """Single-class closed WITH a delay station: the think-time premise of
    harel/sb/scb/sib is what this one violates, and nothing else."""
    model = Network('baGateDelay')
    d = Delay(model, 'D')
    q = Queue(model, 'Q', SchedStrategy.FCFS)
    c = ClosedClass(model, 'C', 3, d)
    d.setService(c, Exp(1.0))
    q.setService(c, Exp(2.0))
    model.link(Network.serialRouting(d, q))
    return model


def all_delay_closed():
    """Two delay stations and nothing else: no queueing station at all, so the
    ldbcmp bottleneck the open-network occupancy is built on does not exist.
    The general form of the shape a Petri net presents, every Place being an INF
    station -- and the shape that made ba_method_degenerate RAISE from inside
    list_valid_methods before it was guarded."""
    model = Network('baGateAllDelay')
    d1 = Delay(model, 'D1')
    d2 = Delay(model, 'D2')
    c = ClosedClass(model, 'C', 2, d1)
    d1.setService(c, Exp(1.0))
    d2.setService(c, Exp(2.0))
    model.link(Network.serialRouting(d1, d2))
    return model


def forkjoin_spn(ntokens=3):
    """The fork-join stochastic Petri net of jar SpnLpbndTest: four Places and
    three Transitions.

    THE SHAPE THAT SEPARATES THE TWO INDEX SPACES. A Place is a station AND a
    stateful node; a Transition is stateful and NOT a station. So this model has
    4 stations against 7 stateful nodes, and sn.visits -- which is
    stateful-indexed -- cannot be walked against sn.sched or sn.rates, which are
    station-indexed. Reading it wrong is what made ba_method_degenerate raise
    from inside list_valid_methods.
    """
    model = Network('baGateSpn')
    pl = [Place(model, 'P%d' % i) for i in range(4)]
    tf = Transition(model, 'Tf')
    tj = Transition(model, 'Tj')
    tb = Transition(model, 'Tb')
    jc = ClosedClass(model, 'C', ntokens, pl[0])
    mf = tf.addMode('f')
    tf.setDistribution(mf, Exp(1.3))
    tf.setNumberOfServers(mf, 1)
    tf.setEnablingConditions(mf, jc, pl[0], 1)
    tf.setFiringOutcome(mf, jc, pl[1], 1)
    tf.setFiringOutcome(mf, jc, pl[2], 1)
    mj = tj.addMode('j')
    tj.setDistribution(mj, Exp(0.7))
    tj.setNumberOfServers(mj, 1)
    tj.setEnablingConditions(mj, jc, pl[1], 1)
    tj.setEnablingConditions(mj, jc, pl[2], 1)
    tj.setFiringOutcome(mj, jc, pl[3], 1)
    mb = tb.addMode('b')
    tb.setDistribution(mb, Exp(1.9))
    tb.setNumberOfServers(mb, 1)
    tb.setEnablingConditions(mb, jc, pl[3], 1)
    tb.setFiringOutcome(mb, jc, pl[0], 1)
    P = model.initRoutingMatrix()
    P.set(jc, jc, pl[0], tf, 1.0)
    P.set(jc, jc, tf, pl[1], 1.0)
    P.set(jc, jc, tf, pl[2], 1.0)
    P.set(jc, jc, pl[1], tj, 1.0)
    P.set(jc, jc, pl[2], tj, 1.0)
    P.set(jc, jc, tj, pl[3], 1.0)
    P.set(jc, jc, pl[3], tb, 1.0)
    P.set(jc, jc, tb, pl[0], 1.0)
    model.link(P)
    for i in range(4):
        pl[i].setState(np.array([[float(ntokens) if i == 0 else 0.0]]))
    return model


def mm1_open():
    return _open_mm1(Exp(1.0), Exp(2.0), 'baGateMM1')


def mm1_erlang_service():
    """Erlang SERVICE: what all three open families refuse at a queueing
    station."""
    return _open_mm1(Exp(1.0), Erlang.fitMeanAndOrder(0.5, 3), 'baGateErlSvc')


def mm1_erlang_source():
    """Erlang SOURCE with exponential service. THE CONVERSE MODEL: 'snc'
    consumes the arrival law and answers this one, so a gate that dropped the
    Erlang feature outright would hide a bound the user could have had. 'bpt'
    and 'bgt' must still be refused -- they read the mean alone and would bound
    the Poisson system instead."""
    return _open_mm1(Erlang.fitMeanAndOrder(1.0, 3), Exp(2.0), 'baGateErlSrc')


def _open_mm1(arrival, service, name):
    model = Network(name)
    s = Source(model, 'S')
    q = Queue(model, 'Q', SchedStrategy.FCFS)
    k = Sink(model, 'K')
    c = OpenClass(model, 'C')
    s.setArrival(c, arrival)
    q.setService(c, service)
    model.link(Network.serialRouting(s, q, k))
    return model


def cqn_ldbcmp_boundary():
    """Delay + two PS queues, one class, N = 4: exactly the ldbcmp regime
    boundary N == Qhat, where the bound degenerates to the trivial X >= 0."""
    model = Network('baGateLdbcmp')
    d = Delay(model, 'D')
    q1 = Queue(model, 'Q1', SchedStrategy.PS)
    q2 = Queue(model, 'Q2', SchedStrategy.PS)
    c = ClosedClass(model, 'C', 4, d)
    d.setService(c, Exp(1.0))
    q1.setService(c, Exp(2.0))
    q2.setService(c, Exp(3.0))
    model.link(Network.serialRouting(d, q1, q2))
    return model


def ba_rows(model):
    """The (method, runnable) rows findSolver reports for the ba family."""
    t = model.findSolver('', True)
    t = t[t.Solver == 'ba']
    return [(r.Method, bool(r.Runnable), r.Reason) for r in t.itertuples()]


def ba_offered(model):
    """The ba methods the REPORT offers: listValidMethods filtered by the gate,
    which is exactly what findSolver builds its rows from."""
    probe = NativeSolverBA(model)
    return [m for m in probe.listValidMethods()
            if probe.supportsModelMethod(m)[0]]


class TestTheGateNamesTheOffendingThing(unittest.TestCase):

    def test_a_multiclass_model_refuses_the_single_class_families_by_name(self):
        probe = NativeSolverBA(cqn_two_class())
        for method in ('aba.upper', 'gb.lower', 'pbh.upper', 'harel.lower',
                       'ssd.upper', 'ldbcmp.lower', 'auto.upper', 'default'):
            ok, reason = probe.supportsModelMethod(method)
            self.assertFalse(ok, '%s was offered on a two-class model' % method)
            self.assertIn('single-class closed networks only', reason)
            # The reason names the method that will RUN, not the alias asked
            # for: 'default' is the gb.upper it resolves to.
            self.assertIn(ba_resolve_method(method), reason)

    def test_a_multiserver_model_refuses_the_single_server_families_by_name(self):
        probe = NativeSolverBA(cqn_multiserver())
        for method in ('aba.upper', 'bjb.lower', 'gb.upper', 'pbh.upper',
                       'bjbk.lower', 'default'):
            ok, reason = probe.supportsModelMethod(method)
            self.assertFalse(ok, '%s was offered on a multiserver model' % method)
            self.assertIn('multi-server stations', reason)
            # 'ssd' IS the multiserver bound, so the refusal points at it.
            self.assertIn("use 'ssd'", reason)

    def test_the_multiclass_chain_families_refuse_a_multiserver_model_too(self):
        # mwba/cub/mbjb/looping survive a multiclass model but not a multiserver
        # one, and their reason must NOT point at 'ssd': that bound is
        # single-class, so it is no alternative for them.
        probe = NativeSolverBA(cqn_multiserver())
        for method in ('mwba.upper', 'cub.upper', 'mbjb.lower', 'looping.lower'):
            ok, reason = probe.supportsModelMethod(method)
            self.assertFalse(ok)
            self.assertIn('multi-server stations', reason)
            self.assertNotIn("use 'ssd'", reason)

    def test_an_open_model_refuses_the_closed_families_by_name(self):
        probe = NativeSolverBA(mm1_open())
        for method in ('aba.upper', 'gb.upper', 'mwba.upper', 'cub.upper'):
            ok, reason = probe.supportsModelMethod(method)
            self.assertFalse(ok)
            self.assertIn('closed networks only', reason)

    def test_a_delay_station_refuses_the_think_time_families_by_feature_name(self):
        # This half of the gate is the FEATURE SET, not the structural
        # predicate: "does not accept a delay station" is expressible as
        # dropping SchedStrategy_INF, and a feature set can refuse a model for
        # HAVING a construct. The reason therefore names the feature.
        probe = NativeSolverBA(cqn_with_delay())
        for method in ('harel.upper', 'harel.lower', 'sb.upper', 'sb.lower',
                       'scb.upper', 'sib.lower', 'lr.upper', 'lr'):
            ok, reason = probe.supportsModelMethod(method)
            self.assertFalse(ok, '%s was offered on a model with a delay' % method)
            self.assertIn('SchedStrategy_INF', reason)
        # and the ones that DO carry a think time are untouched
        for method in ('gb.upper', 'aba.lower', 'pbh.upper', 'default'):
            ok, _ = probe.supportsModelMethod(method)
            self.assertTrue(ok, '%s was dropped by the delay premise' % method)

    def test_the_open_families_declare_the_mirror_premise(self):
        feats = NativeSolverBA(mm1_open()).getMethodFeatureSet('bpt.lower')
        self.assertNotIn('ClosedClass', feats)
        self.assertNotIn('SchedStrategy_INF', feats)
        self.assertIn('OpenClass', feats)


class TestOnePredicateTwoCallers(unittest.TestCase):

    def test_the_analyzer_raises_the_sentence_the_gate_reports(self):
        # The whole point of factoring the rule out: a caller gets ONE answer
        # whichever gate it meets first. Same predicate, same string.
        for model, method in ((cqn_two_class(), 'aba.upper'),
                              (cqn_multiserver(), 'gb.lower'),
                              (cqn_multiserver(), 'cub.upper'),
                              (mm1_open(), 'mwba.upper')):
            probe = NativeSolverBA(model)
            ok, reason = probe.supportsModelMethod(method)
            self.assertFalse(ok)
            with self.assertRaises(ValueError) as raised:
                solver_ba_analyzer(probe, method, probe.options)
            self.assertEqual(reason, str(raised.exception))

    def test_asking_by_name_is_still_refused_rather_than_answered_with_zeros(self):
        # Never fabricate a result: an inapplicable family must raise, not come
        # back as a table of zeros. This is what MATLAB's silent `if
        # nclasses==1 && nclosedjobs>0` guard used to do.
        with self.assertRaises(Exception):
            BA(cqn_two_class(), 'aba.upper').avgTable()

    def test_the_list_is_a_projection_of_the_same_predicate(self):
        for model in (cqn_two_class(), cqn_multiserver(), cqn_with_delay(),
                      cqn_delay_free(), mm1_open()):
            probe = NativeSolverBA(model)
            sn = model.get_struct()
            for method in probe.listValidMethods():
                self.assertEqual('', ba_method_refusal(sn, method),
                                 '%s is listed and structurally refused' % method)


class TestTheBoundsStillAnswerTheirOwnModels(unittest.TestCase):

    def test_a_delay_free_single_class_closed_network_keeps_every_method(self):
        # Over-tightening is as bad as the leak. The only names this model may
        # lose are the three OPEN families and the four spnlp ones, both of
        # which were already gated before this change.
        model = cqn_delay_free()
        listed = set(ba_offered(model))
        missing = set(NativeSolverBA.list_all_methods()) - listed
        self.assertEqual(
            {'bpt.lower', 'bgt.upper', 'snc.upper',
             'spnlp.upper', 'spnlp.lower', 'spnlp.op.upper', 'spnlp.op.lower'},
            missing)

    def test_every_offered_method_actually_runs(self):
        # The audit condition, in the small: no offered ba.* pair may raise, and
        # none may come back as a table of zeros. OFFERED is the report's own
        # list -- listValidMethods filtered by the gate -- because the two halves
        # of the gate live in different places on purpose: the structural
        # premises narrow the list, and the ones a feature name states are
        # applied by supportsModelMethod, so that a showAll report can still
        # print the offending FEATURE rather than dropping the row silently.
        for model in (cqn_delay_free(), cqn_two_class(), cqn_multiserver(),
                      cqn_with_delay(), mm1_open()):
            offered = ba_offered(model)
            self.assertTrue(offered, 'no bound method offered for %s' % model.getName())
            for method in offered:
                table = BA(model, method).avgTable()
                qlen = table['QLen'].to_numpy().astype(float)
                self.assertFalse(qlen.size and np.allclose(qlen, 0.0),
                                 'ba.%s answered %s with a table of zeros'
                                 % (method, model.getName()))

    def test_the_report_offers_no_ba_row_that_cannot_run(self):
        for model in (cqn_two_class(), cqn_multiserver()):
            for method, runnable, _reason in ba_rows(model):
                if not runnable:
                    continue
                name = method.split('.', 1)[1]
                BA(model, name).avgTable()

    def test_the_multiserver_model_keeps_the_bounds_stated_for_it(self):
        listed = set(ba_offered(cqn_multiserver()))
        for method in ('ssd.upper', 'ssd.lower', 'ldbcmp.lower',
                       'auto.upper', 'auto.lower'):
            self.assertIn(method, listed)

    def test_the_two_class_model_keeps_the_multiclass_bounds(self):
        listed = set(ba_offered(cqn_two_class()))
        for method in ('mwba.upper', 'mwba.lower', 'cub.upper', 'mbjb.lower',
                       'looping.upper', 'looping.lower'):
            self.assertIn(method, listed)


class TestTheExponentialServicePremiseOfTheOpenFamilies(unittest.TestCase):
    """bpt/bgt/snc each need exponential service, and the rule lands in two
    different places because the three do not treat the SOURCE alike."""

    def test_erlang_service_refuses_all_three(self):
        probe = NativeSolverBA(mm1_erlang_service())
        for method in ('bpt.lower', 'bgt.upper', 'snc.upper'):
            ok, reason = probe.supportsModelMethod(method)
            self.assertFalse(ok, '%s was offered on Erlang service' % method)
            self.assertTrue('Erlang' in reason or 'exponential service' in reason, reason)
        self.assertEqual([], ba_offered(mm1_erlang_service()))

    def test_bpt_and_bgt_drop_every_law_but_exp(self):
        # Registry-expressible: both read the mean alone, so a non-exponential
        # law ANYWHERE -- source included -- is silently bounded as if it were
        # Poisson rather than refused. Measured: swapping the Exp(1) source of
        # an M/M/1 for an Erlang of the same mean leaves bgt.upper at QLen
        # 32.6667 and bpt.lower at 1.0, digit for digit.
        probe = NativeSolverBA(mm1_open())
        laws = ('APH', 'Coxian', 'Cox2', 'Erlang', 'HyperExp', 'PH',
                'Det', 'Lognormal', 'Pareto', 'Uniform', 'Weibull')
        for method in ('bpt.lower', 'bgt.upper'):
            feats = probe.getMethodFeatureSet(method)
            self.assertIn('Exp', feats)
            for law in laws:
                self.assertNotIn(law, feats, '%s still declares %s' % (method, law))

    def test_snc_keeps_the_laws_because_it_consumes_the_arrival_one(self):
        probe = NativeSolverBA(mm1_open())
        feats = probe.getMethodFeatureSet('snc.upper')
        for law in ('Erlang', 'Coxian', 'APH', 'PH', 'HyperExp'):
            self.assertIn(law, feats)

    def test_an_erlang_source_keeps_snc_and_drops_bpt_and_bgt(self):
        # The converse of the delta above, and the reason snc's rule is
        # structural: no feature name can say "Erlang at a Queue but not at a
        # Source", so dropping the law would refuse a model snc answers.
        model = mm1_erlang_source()
        probe = NativeSolverBA(model)
        ok, _reason = probe.supportsModelMethod('snc.upper')
        self.assertTrue(ok, 'snc was over-tightened by the Erlang source')
        self.assertEqual(['snc.upper'], ba_offered(model))
        for method in ('bpt.lower', 'bgt.upper'):
            ok, reason = probe.supportsModelMethod(method)
            self.assertFalse(ok)
            self.assertIn('Erlang', reason)
        # ... and it really does run, which is what makes the refusal of the
        # other two a judgement about correctness rather than about coverage.
        BA(model, 'snc.upper').avgTable()

    def test_the_snc_rule_is_structural_and_skips_the_source(self):
        sn_svc = mm1_erlang_service().get_struct()
        sn_src = mm1_erlang_source().get_struct()
        self.assertIn('requires exponential service', ba_method_refusal(sn_svc, 'snc.upper'))
        self.assertEqual('', ba_method_refusal(sn_src, 'snc.upper'))
        # bpt and bgt carry no structural rule at all: theirs is the feature set.
        for method in ('bpt.lower', 'bgt.upper'):
            self.assertEqual('', ba_method_refusal(sn_svc, method))


class TestADegenerateBoundIsNotOffered(unittest.TestCase):
    """A bound whose premises hold but whose formula says nothing is a VALID
    bound and a useless one. It is withheld from the report rather than turned
    into a refusal, because a caller naming it is entitled to the answer."""

    def test_ldbcmp_is_withheld_at_its_regime_boundary(self):
        model = cqn_ldbcmp_boundary()
        sn = model.get_struct()
        why = ba_method_degenerate(sn, 'ldbcmp.lower')
        self.assertIn('Qhat=4.0000', why)
        self.assertIn('N=4', why)
        self.assertIn('trivial bound X >= 0', why)
        # The APPLICABILITY predicate stays silent: the model is inside the
        # method's domain, which is exactly why this is a separate question.
        self.assertEqual('', ba_method_refusal(sn, 'ldbcmp.lower'))
        ok, reason = NativeSolverBA(model).supportsModelMethod('ldbcmp.lower')
        self.assertFalse(ok)
        self.assertEqual(why, reason)
        self.assertNotIn('ldbcmp.lower', ba_offered(model))

    def test_the_analyzer_still_answers_it_when_named(self):
        # Not a refusal: X >= 0 IS a lower bound, so the run publishes it. What
        # changed is that nothing offers it.
        table = BA(cqn_ldbcmp_boundary(), 'ldbcmp.lower').avgTable()
        self.assertTrue(np.allclose(table['QLen'].to_numpy().astype(float), 0.0))

    def test_ldbcmp_survives_where_its_bound_says_something(self):
        for model in (cqn_with_delay(), cqn_multiserver()):
            sn = model.get_struct()
            self.assertEqual('', ba_method_degenerate(sn, 'ldbcmp.lower'))
            self.assertIn('ldbcmp.lower', ba_offered(model))

    def test_a_model_with_no_queueing_station_is_answered_not_raised(self):
        # A PREDICATE MUST NOT RAISE. ba_method_degenerate is asked once per
        # name by list_valid_methods, which runs it before any later sieve has
        # dropped anything, so a model with no queueing station reached
        # pfqn_ldbcmp with an empty demand vector and took the whole listing
        # down with it.
        model = all_delay_closed()
        sn = model.get_struct()
        why = ba_method_degenerate(sn, 'ldbcmp.lower')
        self.assertIn('no queueing station', why)
        # ... and the listing itself completes and withholds the name.
        listed = NativeSolverBA(model).listValidMethods()
        self.assertNotIn('ldbcmp.lower', listed)
        self.assertTrue(listed)

    def test_a_petri_net_is_answered_not_raised(self):
        # THE TWO INDEX SPACES. sn.visits is stateful-indexed while sn.sched and
        # sn.rates are station-indexed, and on this net that is 7 rows against
        # 4. Walking the visit rows against the station space raised
        # "operands could not be broadcast together with shapes (7,) (4,)" here,
        # and IndexOutOfBounds in the JAR twin, from inside list_valid_methods.
        model = forkjoin_spn()
        sn = model.get_struct()
        self.assertNotEqual(int(sn.nstations), int(sn.nstateful))
        why = ba_method_degenerate(sn, 'ldbcmp.lower')
        self.assertIn('no queueing station', why)
        # ... and the listing completes, offering exactly the marking-indexed
        # family and nothing else.
        self.assertEqual(
            {'spnlp.upper', 'spnlp.lower', 'spnlp.op.upper', 'spnlp.op.lower'},
            set(NativeSolverBA(model).listValidMethods()))

    def test_no_other_method_has_a_degenerate_regime(self):
        for model in (cqn_delay_free(), cqn_two_class(), cqn_multiserver(),
                      cqn_with_delay(), mm1_open(), cqn_ldbcmp_boundary(),
                      all_delay_closed(), forkjoin_spn()):
            sn = model.get_struct()
            for method in NativeSolverBA.list_all_methods():
                if ba_resolve_method(method) == 'ldbcmp.lower':
                    continue
                self.assertEqual('', ba_method_degenerate(sn, method), method)


if __name__ == '__main__':
    unittest.main()
