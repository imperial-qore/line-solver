"""
RENV_GENQN  Two-station closed network used as a random-environment stage.

The stage models of the repairmen examples differ only in their service rates,
so they are all built by this one function: a Delay (think time) and a PS Queue
in a cycle, traversed by a single closed class of N jobs. The random-environment
solver then switches between instances of it, one per stage.
"""

import numpy as np

from line_solver import (ClosedClass, Delay, Exp, GlobalConstants, Network,
                         Queue, SchedStrategy, VerboseLevel)


def renv_genqn(rate, N):
    """Build the two-station stage network.

    Args:
        rate: the two service rates, (Delay, Queue).
        N: the closed-class population.

    Returns:
        The linked Network.
    """
    qnet = Network('qn1')

    node = np.empty(2, dtype=object)
    node[0] = Delay(qnet, 'Queue1')
    node[1] = Queue(qnet, 'Queue2', SchedStrategy.PS)

    jobclass = np.empty(1, dtype=object)
    jobclass[0] = ClosedClass(qnet, 'Class1', N, node[0], 0)

    node[0].set_service(jobclass[0], Exp(rate[0]))
    node[1].set_service(jobclass[0], Exp(rate[1]))

    P = qnet.init_routing_matrix()
    P.set(jobclass[0], jobclass[0], [[0, 1], [1, 0]])
    qnet.link(P)

    return qnet


if __name__ == "__main__":
    from line_solver import MVA

    GlobalConstants.set_verbose(VerboseLevel.STD)

    # one stage on its own: the environment is what varies these rates
    qnet = renv_genqn([1.0, 0.5], 4)
    print(MVA(qnet).get_avg_table())
