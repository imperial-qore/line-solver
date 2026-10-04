"""
Compare dispatching strategies on the same cluster.

Uses the Cluster builder to flip between RAND, RROBIN and JSQ and solves
each with SSA (so non-product-form policies work).
"""

from line_solver import GlobalConstants, RoutingStrategy, SchedStrategy, Cluster, SSA, VerboseLevel


if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.STD)

    cluster = (Cluster().set_num_stations(4).set_arrival_rate(1.0).set_service_rate(0.4)
            .set_scheduling(SchedStrategy.PS))

    def solve(model):
        return SSA(model, seed=23000, samples=2000).get_avg_table()

    results = cluster.compare_dispatching(
        solve,
        [RoutingStrategy.RAND, RoutingStrategy.RROBIN])
    for policy, table in results.items():
        print(f"\n=== Dispatching: {policy} ===")
        print(table)
