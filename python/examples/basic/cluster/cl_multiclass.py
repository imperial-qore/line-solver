"""
Two-class open cluster: e.g. interactive vs batch.
"""

import numpy as np

from line_solver import GlobalConstants, MVA, Network, RoutingStrategy, VerboseLevel


if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.STD)

    lambda_rates = np.array([0.3, 0.2])           # two classes
    D = np.array([[1.0, 0.5], [1.0, 0.5]])       # M=2 servers, R=2 classes

    model = Network.cluster_ps(lambda_rates, D, dispatching=RoutingStrategy.RAND)
    print(MVA(model).get_avg_table())
