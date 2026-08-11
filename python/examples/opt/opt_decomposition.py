"""
Tutorial 8: Decomposing a multi-variable problem.

Jointly choosing the number of servers and the service rate couples two
variable types. The decomposition workflow splits them into blocks
(servers first, then rates), solves the blocks in sequence with the
other block's values held fixed, and cycles until the full penalized
objective stabilizes (block coordinate descent). This scales better
than a joint search when many variable types are present, at the price
of possibly stopping in a block-wise optimum.
"""

from line_solver import Network, Queue, Source, Sink, OpenClass, Exp, SchedStrategy
from line_solver import (OptimizationProblem, ServerAllocation, ServiceRate,
                      MinimizeCost, ResponseTimeConstraint,
                      DecompositionWorkflow)

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
problem.add_variable(ServiceRate(queue, jobs, bounds=(1.0, 4.0)))
problem.set_objective(MinimizeCost(server_cost={queue: 10.0},
                                   rate_cost={queue: 20.0}))
problem.add_constraint(ResponseTimeConstraint(queue, jobs, max_value=0.5))

workflow = DecompositionWorkflow(problem)
workflow.auto_decompose()
workflow.set_solver_options(max_iterations=30, popsize=10, seed=42)
result = workflow.solve_sequential(max_cycles=4, tolerance=1e-3)

print(f"Converged       : {result.converged} "
      f"after {result.cycles_completed} cycle(s)")
print(f"Servers         : {result.final_variable_values['Server_servers']}")
print(f"Rate            : {result.final_variable_values['Server_Jobs_rate']:.3f}")
print(f"Objective       : {result.final_objective:.2f}")
