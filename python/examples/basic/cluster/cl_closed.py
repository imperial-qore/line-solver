"""
Closed cluster: 8 jobs cycle between Think delay and 3 PS servers.
"""

import numpy as np

from line_solver import (GlobalConstants, MVA, Network, RoutingStrategy,
                         SchedStrategy, VerboseLevel)


if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.STD)

    N = np.array([8])           # population
    Z = np.array([1.0])          # think time
    D = np.full((3, 1), 1.0)    # mean service time at each server
    strategies = [SchedStrategy.PS] * 3

    model = Network.cluster_closed(N, Z, D, strategies, dispatching=RoutingStrategy.RAND)
    print(MVA(model).get_avg_table())
