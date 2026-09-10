#!/usr/bin/env python3
"""Gallery Example: gallery_fj_quorum"""

from line_solver import *

def gallery_fj_quorum():
    """Closed fork-join network with a 2-of-3 quorum join.

    The join fires on the SECOND of the three sibling tasks; the third is
    discarded when it arrives. Solved exactly by SolverLDES and SolverJMT;
    SolverMVA and SolverNC charge the second order statistic of the branch
    completion times at the join (fj_ordstat_exp).
    """
    model = Network('Fork-Join-Quorum')
    delay = Delay(model, 'Delay')
    queue1 = Queue(model, 'Queue1', SchedStrategy.PS)
    queue2 = Queue(model, 'Queue2', SchedStrategy.PS)
    queue3 = Queue(model, 'Queue3', SchedStrategy.PS)
    fork = Fork(model, 'Fork')
    join = Join(model, 'Join', fork)
    oclass = ClosedClass(model, 'class1', 5, delay)
    delay.setService(oclass, Exp(1.0))
    queue1.setService(oclass, Exp(2.0))
    queue2.setService(oclass, Exp(2.0))
    queue3.setService(oclass, Exp(2.0))
    join.setStrategy(oclass, JoinStrategy.PARTIAL)
    join.setRequired(oclass, 2)
    P = model.init_routing_matrix()
    P.set(oclass, oclass, delay, fork, 1.0)
    P.set(oclass, oclass, fork, queue1, 1.0)
    P.set(oclass, oclass, fork, queue2, 1.0)
    P.set(oclass, oclass, fork, queue3, 1.0)
    P.set(oclass, oclass, queue1, join, 1.0)
    P.set(oclass, oclass, queue2, join, 1.0)
    P.set(oclass, oclass, queue3, join, 1.0)
    P.set(oclass, oclass, join, delay, 1.0)
    model.link(P)
    return model


if __name__ == '__main__':
    model = gallery_fj_quorum()
    print('Model built:', type(model).__name__)
