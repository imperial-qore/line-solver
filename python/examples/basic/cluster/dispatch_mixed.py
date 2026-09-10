"""Parity test: mixed cluster built via the Cluster builder.

2 PS server stations shared by one open class (lambda = 0.5) and one closed
class (N = 3 jobs, think time Z = 1.0), mu = 2.0 for the open class and
mu = 1.5 for the closed one, RAND dispatching. Solved by MVA and re-solved
by simulation via LDES for cross-validation.
"""
from line_solver import *
GlobalConstants.set_verbose(VerboseLevel.STD)

cluster = (Cluster()
        .set_num_stations(2)
        .set_mixed(0.5, 3, 1.0)
        .set_service_rates([[2.0, 1.5], [2.0, 1.5]])
        .set_scheduling(SchedStrategy.PS))

print(MVA(cluster.build()).get_avg_table())
print(LDES(cluster.build(), seed=23000).get_avg_table())
