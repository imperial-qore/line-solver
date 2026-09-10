"""
Basic open cluster.

Source -> Dispatcher (Router) -> 4 PS servers -> Sink, with random dispatching.
The single call to Network.cluster_ps replaces ~20 lines of node wiring.
"""

import numpy as np

from line_solver import GlobalConstants, MVA, Network, RoutingStrategy, VerboseLevel


if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.STD)

    lambda_rates = np.array([0.4])
    D = np.full((4, 1), 1.0)  # mean service time = 1 at each server

    model = Network.cluster_ps(lambda_rates, D, dispatching=RoutingStrategy.RAND)
    print(MVA(model).get_avg_table())
