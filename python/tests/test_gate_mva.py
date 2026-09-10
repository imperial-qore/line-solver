"""
Tests of SolverMVA's per-method support gate: every (solver, method) row
findSolver offers for the MVA family must actually run, and must produce the
model's answer rather than a table of zeros.

WHY THIS IS A TEST AND NOT AN INSPECTION. The gate is the same predicate
SolverAUTO.chooseSolverRanked consults before it delegates and that
listValidMethods projects, so a gate weaker than the analyzer is not a cosmetic
defect in a report: it hands a caller a method name that then answers with
zeros. The closed-population AMVA family is where that bit: bs, aql, qsa, sqni,
tay, scat, lcp, chow, pamb, pami, pamt, clust, dmlin, ab, schmidt and
schmidt-ext each recur on a CLOSED population vector N and are handed (L, N, Z)
alone, so on an open model the recursion runs over an empty set of chains and
falls out with every metric at zero -- silently, and with the row still marked
runnable and 'approx'.

WHAT IS ASSERTED:

 1. THE M/M/1 IDENTITY. lambda = 1, mu = 2 gives rho = 1/2 and E[Q] = 1 at the
    queue exactly, for every method that claims to solve it. The number is a
    property of the model, decided before any code was written, so a method that
    returns 0 there is wrong however it was computed.

 2. THE FAMILY IS WITHHELD, NOT SILENTLY WRONG. The sixteen closed-population
    names must be absent from the open model's report, and asking for one by
    name must raise a sentence that says why.

 3. THE CLOSED DIRECTION IS NOT OVER-TIGHTENED. A closed product-form network
    keeps the whole family, and a pure-delay network keeps it too: with no
    queueing station there is no arrival-instant correction to make and every
    one of these algorithms coincides with the exact delay solution.
"""

import os
import sys
import unittest

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from line_solver import (Network, Source, Queue, Sink, Delay, Fork, Join, OpenClass,
                         ClosedClass, Exp, SchedStrategy, MVA, GlobalConstants,
                         VerboseLevel)

GlobalConstants.setVerbose(VerboseLevel.SILENT)

#: The AMVA algorithms whose recursion is over a closed population vector.
CLOSED_POPULATION = ('bs', 'aql', 'qsa', 'sqni', 'tay', 'scat', 'lcp', 'chow',
                     'pamb', 'pami', 'pamt', 'clust', 'dmlin', 'ab',
                     'schmidt', 'schmidt-ext')


def mm1():
    """Source -> FCFS Queue -> Sink with lambda = 1, mu = 2, so E[Q] = 1."""
    m = Network('mm1')
    s, q, k = Source(m, 'Source'), Queue(m, 'Queue', SchedStrategy.FCFS), Sink(m, 'Sink')
    c = OpenClass(m, 'C1')
    s.setArrival(c, Exp(1.0))
    q.setService(c, Exp(2.0))
    m.link(Network.serialRouting(s, q, k))
    return m


def repairmen():
    """Delay -> FCFS Queue, N = 3: a closed product-form network."""
    m = Network('repairmen')
    d, q = Delay(m, 'Delay'), Queue(m, 'Queue', SchedStrategy.FCFS)
    c = ClosedClass(m, 'C1', 3, d)
    d.setService(c, Exp(1.0))
    q.setService(c, Exp(2.0))
    m.link(Network.serialRouting(d, q))
    return m


def two_delays():
    """Delay -> Delay, N = 2: a closed network with NO queueing station."""
    m = Network('twodelays')
    d1, d2 = Delay(m, 'Delay1'), Delay(m, 'Delay2')
    c = ClosedClass(m, 'C1', 2, d1)
    d1.setService(c, Exp(1.0))
    d2.setService(c, Exp(2.0))
    m.link(Network.serialRouting(d1, d2))
    return m


def class_switching():
    """Delay -> PS Queue -> Delay with the class relabelled on each hop.

    C2 is reached only by switching, so its own population is 0 while the CHAIN
    holds 2. It is the shape that exposed the extended Schmidt leak.
    """
    m = Network('cs')
    d, q = Delay(m, 'D'), Queue(m, 'Q', SchedStrategy.PS)
    c1, c2 = ClosedClass(m, 'C1', 2, d), ClosedClass(m, 'C2', 0, d)
    d.setService(c1, Exp(1.0))
    d.setService(c2, Exp(1.0))
    q.setService(c1, Exp(2.0))
    q.setService(c2, Exp(3.0))
    P = m.initRoutingMatrix()
    P.set(c1, c2, d, q, 1.0)
    P.set(c2, c1, q, d, 1.0)
    m.link(P)
    return m


def class_switching_fcfs():
    """The same shape with an FCFS queue: the station the -ext correction is
    formed at, and the empty class it has no customer of to tag."""
    m = Network('csfcfs')
    d, q = Delay(m, 'D'), Queue(m, 'Q', SchedStrategy.FCFS)
    c1, c2 = ClosedClass(m, 'C1', 2, d), ClosedClass(m, 'C2', 0, d)
    d.setService(c1, Exp(1.0))
    d.setService(c2, Exp(1.0))
    q.setService(c1, Exp(2.0))
    q.setService(c2, Exp(3.0))
    P = m.initRoutingMatrix()
    P.set(c1, c2, d, q, 1.0)
    P.set(c2, c1, q, d, 1.0)
    m.link(P)
    return m


def zero_population_chain():
    """Delay -> FCFS Queue with TWO independent chains, one of population 0.

    No class switching here: each class is its own chain, and chain B simply
    holds no customers. It is the shape the extended Schmidt rule is about once
    the family recurs on chains.
    """
    m = Network('zp')
    d, q = Delay(m, 'D'), Queue(m, 'Q', SchedStrategy.FCFS)
    a, b = ClosedClass(m, 'A', 2, d), ClosedClass(m, 'B', 0, d)
    d.setService(a, Exp(1.0))
    d.setService(b, Exp(1.0))
    q.setService(a, Exp(2.0))
    q.setService(b, Exp(3.0))
    P = m.initRoutingMatrix()
    P.set(a, a, d, q, 1.0)
    P.set(a, a, q, d, 1.0)
    P.set(b, b, d, q, 1.0)
    P.set(b, b, q, d, 1.0)
    m.link(P)
    return m


def fork_join():
    """Source -> Fork -> two FCFS queues -> Join -> Sink."""
    m = Network('fj')
    s, f = Source(m, 'Source'), Fork(m, 'Fork')
    q1 = Queue(m, 'Q1', SchedStrategy.FCFS)
    q2 = Queue(m, 'Q2', SchedStrategy.FCFS)
    j, k = Join(m, 'Join'), Sink(m, 'Sink')
    c = OpenClass(m, 'C1')
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


def hol_open():
    """Source -> HOL Queue -> Sink, two open classes at different priorities."""
    m = Network('hol')
    s, q, k = Source(m, 'Source'), Queue(m, 'Queue', SchedStrategy.HOL), Sink(m, 'Sink')
    hi, lo = OpenClass(m, 'Hi', 0), OpenClass(m, 'Lo', 1)
    s.setArrival(hi, Exp(0.4))
    s.setArrival(lo, Exp(0.4))
    q.setService(hi, Exp(2.0))
    q.setService(lo, Exp(2.0))
    m.link(Network.serialRouting(s, q, k))
    return m


def mva_methods(model, runnable_only=True):
    """The MVA rows of the model's findSolver report, as bare method names."""
    t = model.findSolver('', not runnable_only)
    out = []
    for _, r in t.iterrows():
        if r['Solver'] != 'mva':
            continue
        if runnable_only and not bool(r['Runnable']):
            continue
        out.append(str(r['Method']).split('.', 1)[1])
    return out


class TestMvaGateOpenModel(unittest.TestCase):
    """The open direction: a method the report offers must solve the M/M/1."""

    def test_every_offered_method_returns_the_mm1_queue_length(self):
        # rho = 1/2 gives E[Q] = rho/(1-rho) = 1 at the queue. Every method the
        # report offers is offered as a solution of THIS model, so it has to
        # land on that number to its own accuracy; a table of zeros is the
        # failure this gate exists to make impossible.
        model = mm1()
        offered = mva_methods(model)
        self.assertGreater(len(offered), 10)
        for name in offered:
            qlen = MVA(mm1(), name).avgTable()['QLen'].to_numpy().astype(float)
            self.assertTrue(np.all(np.isfinite(qlen)),
                            "method '%s' returned a non-finite queue length: %s" % (name, qlen))
            self.assertFalse(np.allclose(qlen, 0.0),
                             "method '%s' returned an all-zero queue-length table" % name)
            # The queue is the only station holding jobs; the source row is 0.
            self.assertAlmostEqual(float(np.max(qlen)), 1.0, delta=0.8,
                                   msg="method '%s' put %s at the queue, not ~1" % (name, qlen))

    def test_the_closed_population_family_is_not_offered_on_an_open_model(self):
        offered = set(mva_methods(mm1()))
        for name in CLOSED_POPULATION:
            self.assertNotIn(name, offered)
            self.assertNotIn('amva.' + name, offered)

    def test_asking_for_a_closed_population_method_by_name_raises(self):
        # Contract: the gate withholds the row and the analyzer refuses the run.
        # Silence, or zeros, would be worse than either.
        for name in CLOSED_POPULATION:
            with self.assertRaises(Exception) as ctx:
                MVA(mm1(), name).avgTable()
            self.assertTrue(str(ctx.exception),
                            "method '%s' refused without saying why" % name)

    def test_qna_is_withheld_where_its_station_update_has_no_arm(self):
        # solver_qna decomposes each station as an INF, PS or FCFS centre and
        # has no arm for a priority discipline, so it used to leave the HOL
        # station's row of Q, U, R and T at zero and return the table.
        self.assertNotIn('qna', mva_methods(hol_open()))
        self.assertIn('qna', mva_methods(mm1()))

    def test_the_robust_analyzers_are_withheld_on_a_multiclass_model(self):
        # RQNA and RQT build one uncertainty set per flow from the two moments
        # of a single stream.
        offered = mva_methods(hol_open())
        self.assertNotIn('rqna', offered)
        self.assertNotIn('rqt', offered)
        self.assertIn('rqna', mva_methods(mm1()))

    def test_the_summation_method_is_withheld_on_a_priority_station(self):
        # sum/esum pass every station to sum_closed / sum_closing as an INF, PS,
        # LCFS-PR, FCFS or SIRO centre and refuse the rest by name.
        offered = mva_methods(hol_open())
        self.assertNotIn('sum', offered)
        self.assertNotIn('esum', offered)
        self.assertIn('sum', mva_methods(mm1()))


class TestMvaGateClosedModel(unittest.TestCase):
    """The closed direction: nothing that genuinely applies may be lost."""

    def test_a_closed_product_form_network_keeps_the_whole_family(self):
        offered = set(mva_methods(repairmen()))
        for name in CLOSED_POPULATION:
            self.assertIn(name, offered,
                          "'%s' is defined for this closed product-form model" % name)

    def test_the_family_runs_and_conserves_the_population(self):
        # N = 3 jobs are somewhere, whatever approximation is used.
        for name in CLOSED_POPULATION:
            qlen = MVA(repairmen(), name).avgTable()['QLen'].to_numpy().astype(float)
            self.assertAlmostEqual(float(np.sum(qlen)), 3.0, delta=1e-3,
                                   msg="method '%s' lost the closed population" % name)

    def test_a_pure_delay_network_keeps_the_family_and_solves_it_exactly(self):
        # With no queueing station there is no arrival-instant correction to
        # make, so every one of these algorithms coincides with the exact delay
        # solution. Refusing them there would be an over-tightening, and handing
        # them a zero-row demand matrix is what made them raise or report zeros.
        model = two_delays()
        offered = set(mva_methods(model))
        exact = MVA(two_delays(), 'default').avgTable()['QLen'].to_numpy().astype(float)
        for name in CLOSED_POPULATION:
            if name == 'sqni':
                # pfqn_sqni is a closed form for ONE queueing station with a
                # delay, so listValidMethods withholds it on any other shape;
                # that is a shape rule of its own, not the closed-chain rule.
                self.assertNotIn(name, offered)
                continue
            self.assertIn(name, offered)
            qlen = MVA(two_delays(), name).avgTable()['QLen'].to_numpy().astype(float)
            np.testing.assert_allclose(qlen, exact, rtol=1e-9, atol=1e-9,
                                       err_msg="method '%s' on a pure-delay network" % name)

    def test_mvac_is_withheld_where_it_has_no_recursion(self):
        # pfqn_mvac recurs over single-server fixed-rate queues; it refused a
        # multiserver station by name while the report went on offering it.
        m = Network('ms')
        d, q = Delay(m, 'Delay'), Queue(m, 'Queue', SchedStrategy.FCFS)
        q.setNumberOfServers(3)
        c = ClosedClass(m, 'C1', 4, d)
        d.setService(c, Exp(1.0))
        q.setService(c, Exp(2.0))
        m.link(Network.serialRouting(d, q))
        self.assertNotIn('mvac', mva_methods(m))
        self.assertIn('mvac', mva_methods(repairmen()))


class TestMvaGateOnAwkwardShapes(unittest.TestCase):
    """Shapes whose population vector or topology is not what a kernel assumes."""

    def test_schmidt_ext_is_withheld_where_it_has_no_customer_to_tag(self):
        # Schmidt's EXTENSION corrects an FCFS station from the network with one
        # customer of that class TAGGED, i.e. at population N - 1_r. A chain
        # holding no customer has none to tag, so N_r - 1 is negative and the
        # state lattice prod(N+1) collapses to zero. Plain 'schmidt' forms no
        # such sub-problem and must keep running.
        #
        # The vector is the CHAIN one, which is what the arm recurs on: a class
        # that is empty only because its chain's jobs are sitting in a sibling
        # class is not an empty population, and this rule must not fire on it.
        offered = mva_methods(zero_population_chain())
        self.assertNotIn('schmidt-ext', offered)
        self.assertIn('schmidt', offered)
        with self.assertRaises(Exception) as ctx:
            MVA(zero_population_chain(), 'schmidt-ext').avgTable()
        self.assertIn('customer', str(ctx.exception))
        # ... and the class-switching model, whose C2 is empty only as a class,
        # keeps it.
        self.assertIn('schmidt-ext', mva_methods(class_switching_fcfs()))

    def test_a_class_switching_model_offers_no_row_that_raises(self):
        # The reported leak: findSolver offered 'schmidt-ext' on this model and
        # running it raised. Every row the report offers must run.
        model = class_switching()
        offered = mva_methods(model)
        self.assertIn('schmidt-ext', offered)
        for name in offered:
            qlen = MVA(class_switching(), name).avgTable()['QLen'].to_numpy().astype(float)
            self.assertTrue(np.all(np.isfinite(qlen)),
                            "method '%s' returned a non-finite queue length" % name)

    def test_the_robust_analyzers_are_withheld_on_a_fork_join_model(self):
        # A Join is a synchronisation node, not a queue: it carries no service
        # process, so the index-of-dispersion curve RQNA and RQT read off every
        # station does not exist for it. RQNA dereferenced the absent process and
        # RQT reported an infinite queue length at the Join. QNA keeps Fork/Join:
        # its station loop has an explicit Join arm.
        offered = mva_methods(fork_join())
        self.assertNotIn('rqna', offered)
        self.assertNotIn('rqt', offered)
        self.assertIn('qna', offered)
        for name in ('rqna', 'rqt'):
            with self.assertRaises(Exception):
                MVA(fork_join(), name).avgTable()


class TestMvaGateAgreesWithTheRun(unittest.TestCase):
    """The gate and the analyzer must be the one predicate, not two copies."""

    def test_every_offered_row_runs_on_every_shape(self):
        for build in (mm1, repairmen, two_delays, hol_open, class_switching,
                      class_switching_fcfs, zero_population_chain, fork_join):
            model = build()
            for name in mva_methods(model):
                try:
                    table = MVA(build(), name).avgTable()
                except Exception as exc:  # pragma: no cover - the assertion is the report
                    self.fail("%s offers '%s' but the run raised: %s"
                              % (model.getName(), name, exc))
                qlen = table['QLen'].to_numpy().astype(float)
                if qlen.size:
                    self.assertFalse(
                        np.allclose(qlen, 0.0),
                        "%s offers '%s' and it answered with an all-zero table"
                        % (model.getName(), name))


if __name__ == '__main__':
    unittest.main()
