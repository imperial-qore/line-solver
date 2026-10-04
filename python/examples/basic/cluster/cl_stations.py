"""
CL_STATIONS  Horizontal-scaling sweep: response time against the number of
stations at a fixed total arrival rate.

Cluster.sweep_num_stations replicates the single-station service rate across M
stations and returns a dictionary keyed by M. Adding stations splits the same
arrival stream, so per-station utilization falls and response time drops
towards the bare service time.
"""

from line_solver import (GlobalConstants, RoutingStrategy, SchedStrategy,
                         Cluster, MVA, VerboseLevel)


if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.STD)

    cluster = Cluster().set_num_stations(1).set_arrival_rate(1.6).set_service_rate(1.0)
    cluster.set_scheduling(SchedStrategy.PS).set_dispatching(RoutingStrategy.RAND)

    results = cluster.sweep_num_stations([2, 3, 4, 6], MVA)

    for count in sorted(results):
        table = results[count]
        # the Source row carries no service, so it is excluded from the summary
        rows = table[table['Station'] != 'Source']
        print(f"\n=== M = {count}  (mean RespT {rows['RespT'].mean():.4f}, "
              f"max Util {rows['Util'].max():.4f}) ===")
        print(table)
