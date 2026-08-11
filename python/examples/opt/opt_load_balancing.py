"""
Tutorial 4: Load balancing across heterogeneous servers.

Split an arrival stream (lambda = 2) between a fast queue (mu = 4) and a
slow queue (mu = 2.5) so that the mean end-to-end response time is
minimized. Routing probabilities are the decision variable; the system
response time is computed by Little's law from the chain queue lengths
and throughput. The exact optimum (from M/M/1 formulas) routes about 74%
of the traffic to the fast queue, not 100%: saturating the fast server
is worse than using both.
"""

from line_solver import Network, Queue, Source, Sink, OpenClass, Exp, SchedStrategy
from line_solver import (OptimizationProblem, RoutingProbabilities,
                      MinimizeSystemResponseTime)

model = Network("LoadBalance")
source = Source(model, "S")
fast = Queue(model, "Fast", SchedStrategy.PS)
slow = Queue(model, "Slow", SchedStrategy.PS)
sink = Sink(model, "K")
jobs = OpenClass(model, "Jobs")
source.setArrival(jobs, Exp(2.0))
fast.setService(jobs, Exp(4.0))
slow.setService(jobs, Exp(2.5))
model.addLink(source, fast)
model.addLink(source, slow)
model.addLink(fast, sink)
model.addLink(slow, sink)

problem = OptimizationProblem(model)
problem.add_variable(RoutingProbabilities(jobs, source=source,
                                          targets=[fast, slow]))
problem.set_objective(MinimizeSystemResponseTime(jobs))

result = problem.solve(seed=42, max_iterations=60)

probs = result.variable_values['Jobs_routing_from_S']
print(f"Fraction to Fast: {probs[0]:.3f} (theory: 0.743)")
print(f"Fraction to Slow: {probs[1]:.3f}")
print(f"System RT       : {result.objective_value:.4f} (theory: 0.4264)")
