"""
Fork-Join with a Variable Forking Level

The Fork emits a number of tasks per link that is NOT one fixed number. Four
modes select which way the degree varies:

    'fixed'   the classic fork, one task on each of two links (the baseline)
    'vector'  three tasks towards Queue2 and one towards Queue1
    'random'  one or three tasks per link, each with probability 1/2
    'prob'    the branch towards Queue2 fires only half the time

Exact under SolverJMT and SolverLDES, which draw the degree at the fork epoch.
SolverMVA's MMT method sees the EXPECTED degree. The exact CTMC/SSA path accepts
'fixed' and 'vector' -- on 'vector' its answer matches JMT to simulation noise --
and refuses 'random' and 'prob' by name, because its tag construction fixes the
sibling count, and the sibling SET, when the state space is built.
"""

from line_solver import *
import numpy as np


def fj_variable_fanout(mode='fixed'):

    model = Network('model')

    delay = Delay(model, 'Delay')
    queue1 = Queue(model, 'Queue1', SchedStrategy.PS)
    queue2 = Queue(model, 'Queue2', SchedStrategy.PS)
    fork = Fork(model, 'Fork')
    join = Join(model, 'Join', fork)

    jobclass1 = ClosedClass(model, 'class1', 4, delay)

    delay.set_service(jobclass1, Exp(1.0))
    queue1.set_service(jobclass1, Exp(2.0))
    queue2.set_service(jobclass1, Exp(2.0))

    if mode == 'fixed':
        pass  # the baseline every other mode is compared against
    elif mode == 'vector':
        fork.setTasksPerLink(3, jobclass1, queue2)
    elif mode == 'random':
        fork.setTasksPerLinkDistribution(jobclass1, DiscreteSampler([0.5, 0.5], [1, 3]))
    elif mode == 'prob':
        # an uncertain branch needs a Join that does not wait for it
        join.setStrategy(jobclass1, JoinStrategy.PARTIAL)
        join.setRequired(jobclass1, 1)
        fork.setBranchProbability(jobclass1, queue2, 0.5)
    else:
        raise ValueError('Unknown mode; use fixed, vector, random or prob.')

    P = model.init_routing_matrix()
    P.set(jobclass1, jobclass1, delay, fork, 1.0)
    P.set(jobclass1, jobclass1, fork, queue1, 1.0)
    P.set(jobclass1, jobclass1, fork, queue2, 1.0)
    P.set(jobclass1, jobclass1, queue1, join, 1.0)
    P.set(jobclass1, jobclass1, queue2, join, 1.0)
    P.set(jobclass1, jobclass1, join, delay, 1.0)
    model.link(P)
    return model


if __name__ == '__main__':
    GlobalConstants.set_verbose(VerboseLevel.STD)
    for mode in ('fixed', 'vector', 'random', 'prob'):
        model = fj_variable_fanout(mode)
        print('\nmode = %s' % mode)
        print(JMT(model, seed=23000).avg_table())
