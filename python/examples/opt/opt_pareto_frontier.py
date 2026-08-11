"""
Tutorial 7: Cost-performance tradeoff frontier.

How much does tightening the utilization SLA cost? ParetoSweep re-solves
the sizing problem for a sweep of utilization bounds (epsilon-constraint
method) and reports the non-dominated cost frontier. Each point is
solved exactly with the bisection solver.
"""

from line_solver import Network, Queue, Source, Sink, OpenClass, Exp, SchedStrategy
from line_solver import (OptimizationProblem, ServerAllocation,
                      MinimizeCost, UtilizationConstraint, ParetoSweep)

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
problem.add_variable(ServerAllocation(queue, bounds=(1, 12)))
problem.set_objective(MinimizeCost(server_cost={queue: 10.0}))

sweep = ParetoSweep(
    problem,
    constraint_factory=lambda eps: UtilizationConstraint(queue, max_value=eps),
    epsilons=[0.3, 0.4, 0.5, 0.6, 0.75, 0.9],
    solver='bisection',
)
sweep.solve()

print("Utilization bound -> optimal cost (servers)")
for point in sweep.getFrontier():
    servers = point.result.variable_values['Server_servers']
    print(f"  util <= {point.epsilon:<4} -> {point.objective_value:6.1f}  ({servers})")

# Plot the actual (non-dominated) frontier as a cost staircase. Dominated
# sample points (e.g. util <= 0.9, same 4 servers as util <= 0.75) are shown
# faintly and are not part of the frontier line.
import os
out = os.path.join(os.path.dirname(__file__), 'opt_pareto_frontier.png')
sweep.plot(xlabel='utilization bound', ylabel='optimal cost',
           title='Cost-performance Pareto frontier', save_path=out)
print(f"Frontier plot saved to {out}")
