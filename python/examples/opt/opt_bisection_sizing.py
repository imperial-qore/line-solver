"""
Tutorial 2: Exact sizing with the bisection solver.

Same sizing question as Tutorial 1, but solved exactly in O(log n) LINE
evaluations with BisectionSolver instead of differential evolution. This
applies whenever a single integer variable has monotone feasibility:
adding servers can only make the constraints easier to satisfy.
"""

from line_solver import Network, Queue, Source, Sink, OpenClass, Exp, SchedStrategy
from line_solver import (OptimizationProblem, ServerAllocation,
                      MinimizeCost, UtilizationConstraint, BisectionSolver)

model = Network("MMc")
source = Source(model, "Arrivals")
queue = Queue(model, "Server", SchedStrategy.FCFS)
sink = Sink(model, "Departures")
jobs = OpenClass(model, "Jobs")
source.setArrival(jobs, Exp(3.0))
queue.setService(jobs, Exp(1.0))
model.addLink(source, queue)
model.addLink(queue, sink)

problem = OptimizationProblem(model)
problem.add_variable(ServerAllocation(queue, bounds=(1, 10)))
problem.set_objective(MinimizeCost(server_cost={queue: 10.0}))
problem.add_constraint(UtilizationConstraint(queue, max_value=0.5))

result = BisectionSolver(problem).solve()

print(f"Optimal servers : {result.variable_values['Server_servers']}")
print(f"Cost            : {result.objective_value:.1f}")
print(f"LINE evaluations: {result.model_evaluations} (vs hundreds for DE)")
