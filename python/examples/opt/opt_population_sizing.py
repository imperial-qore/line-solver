"""
Tutorial 5: Population sizing in a closed network.

An interactive system is modeled as a closed network: a Delay station
(think time, mean 1) and a PS server (rate 2). Find the largest number
of concurrent users the system can sustain while the server response
time stays within an SLA of 2 time units. Feasibility is monotone
non-increasing in the population, so the bisection solver applies with
direction='max_feasible'.
"""

from line_solver import Network, Queue, Delay, ClosedClass, Exp, SchedStrategy
from line_solver import (OptimizationProblem, JobPopulation, MinimizeCost,
                      ResponseTimeConstraint, BisectionSolver)

model = Network("Interactive")
think = Delay(model, "Think")
server = Queue(model, "AppServer", SchedStrategy.PS)
users = ClosedClass(model, "Users", 1, think)
think.setService(users, Exp(1.0))     # think rate 1 (mean think time 1)
server.setService(users, Exp(2.0))    # service rate 2
model.addLink(think, server)
model.addLink(server, think)

problem = OptimizationProblem(model)
problem.add_variable(JobPopulation(users, bounds=(1, 50)))
problem.set_objective(MinimizeCost())   # pure feasibility sizing
problem.add_constraint(ResponseTimeConstraint(server, users, max_value=2.0))

result = BisectionSolver(problem, direction='max_feasible').solve()

print(f"Max users       : {result.variable_values['Users_population']}")
print(f"Feasible        : {result.feasible}")
print(f"LINE evaluations: {result.model_evaluations}")
