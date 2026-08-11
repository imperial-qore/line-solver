"""Parity test: open cluster built via the Cluster builder.

3 PS server stations, single open class, lambda = 0.4, mu = 1.0, RAND
dispatching. Solved exactly via MVA, then re-solved by stochastic
simulation via SSA for cross-validation.
"""
from line_solver import *
GlobalConstants.set_verbose(VerboseLevel.STD)

cluster = (Cluster()
        .set_num_stations(3)
        .set_arrival_rate(0.4)
        .set_service_rate(1.0)
        .set_scheduling(SchedStrategy.PS))

model = cluster.build()
print(MVA(model).get_avg_table())
print(SSA(model, seed=23000, samples=20000).get_avg_table())
