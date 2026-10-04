#!/usr/bin/env python3
"""Gallery Example: gallery_lqn_random (randomly generated, reproducible)"""

import random
import numpy as np
from line_solver import *

def gallery_lqn_random(seed=23000):
    """Randomly generated layered queueing network (reproducible).

    Uses LayeredNetworkGenerator with a fixed default seed so repeated calls
    yield the same model. 1 client, 2 levels, 4 tasks, 2 processors.
    """
    random.seed(seed)
    np.random.seed(seed)
    gen = LayeredNetworkGenerator()
    return gen.generate(1, 2, 4, 2)

if __name__ == '__main__':
    model = gallery_lqn_random()
    print('Model built:', type(model).__name__)
