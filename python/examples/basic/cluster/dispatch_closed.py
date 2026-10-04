"""Parity test: closed cluster built via the Cluster builder.

3 PS server stations, single closed class with N = 10 jobs, think time
Z = 1.0, mu = 1.0, RAND dispatching. Solved exactly via MVA, then
re-solved by stochastic simulation via SSA for cross-validation.
"""
from line_solver import *
GlobalConstants.set_verbose(VerboseLevel.STD)

cluster = (Cluster()
        .set_num_stations(3)
        .set_service_rate(1.0)
        .set_scheduling(SchedStrategy.PS)
        .set_closed(10, 1.0))

# Bound rather than built inline: the parity rows serialize the example's own
# Network, and an unbound `cluster.build()` leaves them nothing to export.
model = cluster.build()
model2 = cluster.build()

print(MVA(model).get_avg_table())
print(SSA(model2, seed=23000, samples=20000).get_avg_table())
