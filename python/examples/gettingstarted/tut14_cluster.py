# %%
# Example 13: Open cluster
#
# Source -> Dispatcher (Router) -> Server[1..M] -> Sink, with a single open
# class. Two equivalent ways to build the model are shown:
#   (a) the one-liner Network.cluster_ps static factory;
#   (b) the chainable Cluster builder, which also exposes helpers to
#       compare dispatching policies on the same cluster.
from line_solver import *
GlobalConstants.set_verbose(VerboseLevel.STD)
# %%
# Block 1: one-liner factory
lam = [0.4]                     # arrival rate of the open class
D = [[1.0], [1.0], [1.0]]       # mean service time = 1 at each of 3 servers
model = Network.cluster_ps(lam, D, dispatching=RoutingStrategy.RAND)

avg_table = MVA(model).get_avg_table()
print(avg_table)
# %%
# Block 2: Cluster builder with non-uniform multi-server queues
cluster = (Cluster().set_num_stations(3).set_arrival_rate(0.4).set_service_rate(1.0)
        .set_scheduling(SchedStrategy.FCFS)
        .set_station_servers([2, 1, 1]))   # Server1 is M/M/2; Server2/3 are M/M/1
avg_table_fcfs = MVA(cluster.build()).get_avg_table()
print(avg_table_fcfs)
# %%
# Block 3: cross-check the same FCFS multi-server model under three simulators.
#  JMT is the Java-based discrete-event simulator (XML-driven); LDES is a
#  SSJ-backed discrete-event simulator (subprocess-invoked); SSA is
#  LINE's native stochastic simulator using the next-reaction method.
#  All three produce statistically equivalent results on this open-class cluster.
print('=== JMT ===')
print(JMT(cluster.build(), seed=23000, samples=20000).get_avg_table())
print('=== LDES ===')
print(LDES(cluster.build(), seed=23000, samples=20000).get_avg_table())
print('=== SSA ===')
print(SSA(cluster.build(), seed=23000, samples=20000).get_avg_table())
# %%
# Block 4: compare dispatching policies via simulation.
#  MVA assumes RAND (product-form). For non-product-form policies such as
#  RROBIN we drop to a simulator with a small sample budget.
cluster2 = (Cluster().set_num_stations(3).set_arrival_rate(0.4).set_service_rate(1.0)
         .set_scheduling(SchedStrategy.PS))
solver_fcn = lambda m: JMT(m, seed=23000, samples=5000).get_avg_table()

for policy in (RoutingStrategy.RAND, RoutingStrategy.RROBIN):
    cluster2.set_dispatching(policy)
    print(f"\n=== Dispatching: {policy} ===")
    print(solver_fcn(cluster2.build()))
