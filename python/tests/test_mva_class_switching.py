"""
The closed-population AMVA family on a CLASS-SWITCHING model.

THE DEFECT THIS PINS. These sixteen algorithms recur on a population vector, and
the recursion presumes each entry of it is CONSERVED: the arrival-instant
estimate E[Q(N - 1_r)] is only meaningful if removing a customer of r leaves a
network of the same shape. Under class switching a job CHANGES CLASS as it
moves, so no per-class population is conserved -- the conserved quantity is the
CHAIN population, and the class populations are not even a partition of the jobs
in flight. MATLAB (`solver_amva.m`, via `sn_get_product_form_chain_params` and
`sn_deaggregate_chain_results`) and the C++ port build their whole product-form
branch on chains for exactly that reason.

Native python handed the kernels the CLASS vector instead, and so solved a
different network. On the model below -- two classes in one chain, C2 reached
only by switching so its own population is 0 -- every one of the sixteen names
returned `[2, 0]`: both jobs parked at the delay, none at the queue. Sixteen
different approximations agreeing BIT FOR BIT is the signature of a structural
failure, not of accuracy, and that is what these tests assert against.

WHAT IS ASSERTED, and why none of it is a number read back out of the
implementation:

 1. THE EXACT ANSWER IS INDEPENDENT. `SolverCTMC` solves the generator as
    written; it does not go near the AMVA kernels and is unaffected by this fix.
    It is the reference every approximation here is measured against.

 2. THE APPROXIMATIONS MUST DISAGREE WITH EACH OTHER. Bard-Schweitzer, Bard's
    LCP, Chow's second approximation and the Hsieh-Lam proportional forms are
    different estimators; on a model with any queueing at all they cannot all
    land on the same number. Identical answers are the bug, so their spread is
    asserted directly.

 3. NOTHING WITHOUT CLASS SWITCHING MOVES. Chain and class populations carry the
    same numbers when every chain holds one class, so the same models must give
    the same answers as the exact recursion does. Asserted here on a closed
    two-class network with no switching.
"""

import os
import sys
import unittest

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from line_solver import (Network, Queue, Delay, ClosedClass, Exp, SchedStrategy,
                         MVA, CTMC, GlobalConstants, VerboseLevel)

GlobalConstants.setVerbose(VerboseLevel.SILENT)

#: The AMVA algorithms whose recursion is over a closed population vector.
CLOSED_POPULATION = ('bs', 'aql', 'qsa', 'tay', 'scat', 'lcp', 'chow',
                     'pamb', 'pami', 'pamt', 'clust', 'dmlin', 'ab',
                     'schmidt', 'schmidt-ext')

#: 'sqni' is left out of the accuracy assertions on purpose. Its closed form
#: reports Q = N - X Z off a square-root estimate of X, and on this model that
#: estimate is the saturation value X = 2, which drives the queue term to
#: exactly zero. That is the algorithm, not the defect: the C++ reference
#: (`common/line-cli -s mva --method sqni`) reports the same 2 and 0. It is
#: covered by the run-and-conserve test instead.
SQNI = 'sqni'


def class_switching(sched=SchedStrategy.PS):
    """Delay -> Queue -> Delay with the class relabelled on each hop.

    One chain of 2 jobs spread over two classes: C1 is served at the delay and
    switches to C2 on the way to the queue, C2 switches back on the way out. C2
    therefore has population 0 of its own while the chain holds 2.
    """
    m = Network('cs')
    d, q = Delay(m, 'D'), Queue(m, 'Q', sched)
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


def no_switching():
    """Delay -> PS Queue, two closed classes that keep their identity."""
    m = Network('ns')
    d, q = Delay(m, 'D'), Queue(m, 'Q', SchedStrategy.PS)
    c1, c2 = ClosedClass(m, 'C1', 2, d), ClosedClass(m, 'C2', 1, d)
    d.setService(c1, Exp(1.0))
    d.setService(c2, Exp(2.0))
    q.setService(c1, Exp(2.0))
    q.setService(c2, Exp(3.0))
    P = m.initRoutingMatrix()
    P.set(c1, c1, d, q, 1.0)
    P.set(c1, c1, q, d, 1.0)
    P.set(c2, c2, d, q, 1.0)
    P.set(c2, c2, q, d, 1.0)
    m.link(P)
    return m


def qlen(model, method):
    return MVA(model, method).avgTable()['QLen'].to_numpy().astype(float)


class TestClassSwitchingAmva(unittest.TestCase):

    def test_the_exact_answer_is_what_ctmc_says(self):
        # The reference, computed by a solver that never touches these kernels.
        exact = CTMC(class_switching()).avgTable()['QLen'].to_numpy().astype(float)
        np.testing.assert_allclose(exact, [1.4117647058823530, 0.5882352941176471],
                                   rtol=1e-9, atol=1e-9)

    def test_the_family_approaches_the_exact_answer(self):
        # Each of these is an approximation, so each is allowed its own error --
        # but none may be off by the whole queue. The old answer, [2, 0], is 0.59
        # out at both stations, which is the entire content of the queue row.
        exact = CTMC(class_switching()).avgTable()['QLen'].to_numpy().astype(float)
        for name in CLOSED_POPULATION:
            q = qlen(class_switching(), name)
            self.assertEqual(q.shape, exact.shape)
            self.assertGreater(q[1], 0.4, "method '%s' left the queue empty: %s" % (name, q))
            np.testing.assert_allclose(
                q, exact, atol=0.25,
                err_msg="method '%s' is not near the exact answer" % name)

    def test_the_approximations_do_not_all_agree(self):
        # Sixteen different estimators returning one number is the signature of
        # the defect. Bard-Schweitzer, Bard's LCP and the Hsieh-Lam proportional
        # form are derived differently and must land differently.
        seen = set()
        for name in CLOSED_POPULATION:
            seen.add(round(float(qlen(class_switching(), name)[1]), 9))
        self.assertGreater(len(seen), 3,
                           "the family collapsed onto %d distinct answers" % len(seen))
        self.assertNotEqual(round(float(qlen(class_switching(), 'bs')[1]), 9),
                            round(float(qlen(class_switching(), 'lcp')[1]), 9))

    def test_every_family_member_runs_and_conserves_the_population(self):
        # N = 2 jobs are somewhere, whatever approximation is used, and the
        # class-level answer has to add up after the chain is split back over
        # its classes.
        for name in CLOSED_POPULATION + (SQNI,):
            q = qlen(class_switching(), name)
            self.assertAlmostEqual(float(np.sum(q)), 2.0, delta=1e-6,
                                   msg="method '%s' lost the closed population" % name)

    def test_an_fcfs_class_switching_model_is_solved_too(self):
        # The same shape with the queue served FCFS, which is where the extended
        # Schmidt correction is formed.
        exact = CTMC(class_switching(SchedStrategy.FCFS)).avgTable()['QLen'].to_numpy().astype(float)
        for name in ('schmidt', 'schmidt-ext', 'ab', 'bs'):
            q = qlen(class_switching(SchedStrategy.FCFS), name)
            np.testing.assert_allclose(
                q, exact, atol=0.25,
                err_msg="method '%s' is not near the exact answer" % name)

    def test_a_model_without_class_switching_still_matches_exact_mva(self):
        # Chain and class populations carry the same numbers when every chain
        # holds one class, so the aggregation must be invisible here. The
        # product-form methods stay exact on this product-form model.
        exact = qlen(no_switching(), 'exact')
        for name in ('schmidt', 'schmidt-ext', 'aql', 'qsa'):
            np.testing.assert_allclose(
                qlen(no_switching(), name), exact, atol=5e-3,
                err_msg="method '%s' moved on a model with no class switching" % name)


if __name__ == '__main__':
    unittest.main()
