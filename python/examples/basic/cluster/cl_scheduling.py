"""
CL_SCHEDULING  Compare scheduling disciplines on the same cluster.

Cluster.compare_scheduling rebuilds the model once per discipline and returns a
dictionary from the discipline to its AvgTable, restoring the original setting
on exit. With exponential service FCFS and PS agree in the mean; the contrast
appears once the service SCV departs from 1.
"""

from line_solver import (GlobalConstants, RoutingStrategy, SchedStrategy,
                         Cluster, SSA, VerboseLevel)


if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.STD)

    cluster = Cluster().set_num_stations(3).set_arrival_rate(0.9).set_service_rate(0.5)
    cluster.set_dispatching(RoutingStrategy.RAND)
    cluster.set_service_scv(4.0)

    def solve(model):
        return SSA(model, seed=23000, samples=20000).get_avg_table()

    disciplines = [SchedStrategy.FCFS, SchedStrategy.PS]
    results = cluster.compare_scheduling(solve, disciplines)

    for discipline, table in results.items():
        print(f"\n=== Scheduling: {getattr(discipline, 'name', discipline)} ===")
        print(table)
