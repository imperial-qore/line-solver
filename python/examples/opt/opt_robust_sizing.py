"""
Tutorial 6: Robust sizing under workload scenarios.

Size the server pool so that the utilization SLA holds both under the
average load (lambda = 3) and under a peak-load scenario (lambda = 4.5).
Scenario models share the node and class names of the base model; the
solver enforces every constraint on all scenarios. The average load
alone needs 6 servers, but the peak scenario pushes the answer to 9.
"""

from line_solver import Network, Queue, Source, Sink, OpenClass, Exp, SchedStrategy
from line_solver import (OptimizationProblem, ServerAllocation,
                      MinimizeCost, UtilizationConstraint, BisectionSolver)


def build(arrival_rate):
    model = Network("MMc")
    source = Source(model, "Arrivals")
    queue = Queue(model, "Server", SchedStrategy.FCFS)
    sink = Sink(model, "Departures")
    jobs = OpenClass(model, "Jobs")
    source.setArrival(jobs, Exp(arrival_rate))
    queue.setService(jobs, Exp(1.0))
    model.addLink(source, queue)
    model.addLink(queue, sink)
    return model, queue


base_model, queue = build(3.0)        # average load
peak_model, _ = build(4.5)            # peak load scenario

problem = OptimizationProblem(base_model)
problem.add_variable(ServerAllocation(queue, bounds=(1, 12)))
problem.set_objective(MinimizeCost(server_cost={queue: 10.0}))
problem.add_constraint(UtilizationConstraint(queue, max_value=0.5))
problem.add_scenario(peak_model)

result = BisectionSolver(problem).solve()

print(f"Robust servers  : {result.variable_values['Server_servers']} "
      f"(6 for average load only)")
print(f"Cost            : {result.objective_value:.1f}")
print(f"Feasible        : {result.feasible}")
