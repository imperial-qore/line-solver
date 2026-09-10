#!/usr/bin/env python3
"""Gallery Example: gallery_qn_random (randomly generated, reproducible)"""

import random
import numpy as np
from line_solver import *

def gallery_qn_random(seed=23000):
    """Randomly generated mixed queueing network (reproducible).

    Uses NetworkGenerator with a deterministic (cyclic) topology and a fixed
    default seed so repeated calls yield the same model. Closed network:
    3 queues, 1 delay, 2 closed classes (always stable / solvable).
    """
    random.seed(seed)
    np.random.seed(seed)
    gen = NetworkGenerator(
        sched_strat='fcfs', routing_strat='Probabilities', distribution='Exp',
        cclass_job_load='medium', has_varying_service_rates=False,
        has_multi_server_queues=False, has_random_cs_nodes=False,
        has_multi_chain_cs=False, topology_fcn=cyclic_graph)
    return gen.generate(3, 1, 0, 2)



if __name__ == '__main__':
    model = gallery_qn_random()
    print('Model built:', type(model).__name__)
