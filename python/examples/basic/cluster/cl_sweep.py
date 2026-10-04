"""
Arrival-rate sweep on a 2-server PS cluster. Demonstrates how response time grows
as the system approaches saturation.
"""

from line_solver import (GlobalConstants, RoutingStrategy, SchedStrategy,
                         Cluster, MVA, VerboseLevel)


if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.STD)

    cluster = (Cluster().set_num_stations(2).set_arrival_rate(0.1).set_service_rate(1.0)
            .set_scheduling(SchedStrategy.PS)
            .set_dispatching(RoutingStrategy.RAND))

    rates = [0.2, 0.5, 0.9, 1.5]
    results = cluster.sweep_arrival_rate(rates, MVA)
    for lam, table in results.items():
        print(f"\n=== lambda = {lam} ===")
        print(table)
