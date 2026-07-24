"""
Tutorial 3: Continuous service rate tuning.

Choose the cheapest service rate (a continuous decision) meeting a
response time SLA on an M/M/1 queue. With lambda = 3, the response time
is RT = 1/(mu - 3), so RT <= 0.5 requires mu >= 5; rate is billed at
20 cost units per unit, hence the optimum is mu* = 5 (cost 100).
"""

from line_solver import Network, Queue, Source, Sink, OpenClass, Exp, SchedStrategy
from line_solver import (OptimizationProblem, ServiceRate,
                      MinimizeCost, ResponseTimeConstraint)

model = Network("MM1")
source = Source(model, "Arrivals")
queue = Queue(model, "Server", SchedStrategy.FCFS)
sink = Sink(model, "Departures")
jobs = OpenClass(model, "Jobs")
source.setArrival(jobs, Exp(3.0))
queue.setService(jobs, Exp(4.0))    # initial rate, will be optimized
model.addLink(source, queue)
model.addLink(queue, sink)

problem = OptimizationProblem(model)
problem.add_variable(ServiceRate(queue, jobs, bounds=(3.5, 8.0)))
problem.set_objective(MinimizeCost(rate_cost={queue: 20.0}))
problem.add_constraint(ResponseTimeConstraint(queue, jobs, max_value=0.5))

result = problem.solve(seed=42, max_iterations=60)

rate = result.variable_values['Server_Jobs_rate']
print(f"Optimal rate    : {rate:.3f} (theory: 5.0)")
print(f"Cost            : {result.objective_value:.1f}")
print(f"Feasible        : {result.feasible}")
