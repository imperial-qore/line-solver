"""
Tutorial 1: Server sizing with differential evolution.

Find the minimum number of servers of an M/M/c queue such that the
utilization does not exceed 50%. Jobs arrive at rate lambda = 3 and each
server works at rate mu = 1, so the utilization with c servers is
rho = 3/c and the optimum is c* = 6 (cost 60).
"""

from line_solver import Network, Queue, Source, Sink, OpenClass, Exp, SchedStrategy
from line_solver import (OptimizationProblem, ServerAllocation,
                      MinimizeCost, UtilizationConstraint)

# LINE model: Source -> Queue -> Sink
model = Network("MMc")
source = Source(model, "Arrivals")
queue = Queue(model, "Server", SchedStrategy.FCFS)
sink = Sink(model, "Departures")
jobs = OpenClass(model, "Jobs")
source.setArrival(jobs, Exp(3.0))   # lambda = 3
queue.setService(jobs, Exp(1.0))    # mu = 1 per server
model.addLink(source, queue)
model.addLink(queue, sink)

# Optimization problem: choose 1..10 servers at 10 cost units each
problem = OptimizationProblem(model)
problem.add_variable(ServerAllocation(queue, bounds=(1, 10)))
problem.set_objective(MinimizeCost(server_cost={queue: 10.0}))
problem.add_constraint(UtilizationConstraint(queue, max_value=0.5))

result = problem.solve(seed=42)

print(f"Optimal servers : {result.variable_values['Server_servers']}")
print(f"Cost            : {result.objective_value:.1f}")
print(f"Feasible        : {result.feasible}")
print(f"LINE evaluations: {result.model_evaluations}")
