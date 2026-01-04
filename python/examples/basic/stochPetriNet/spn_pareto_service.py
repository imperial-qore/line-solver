"""
Basic stochastic Petri net with Pareto service time.

This example demonstrates:
- Open network: Source → Place → Transition → Sink
- Single transition T1 with Pareto firing time (shape=3, scale=1)
- Pareto distribution has mean = shape*scale/(shape-1) = 3*1/(3-1) = 1.5
- Non-Markovian (heavy-tailed) firing times in Petri nets
- Single server capacity for transition
"""

from line_solver import *

if __name__ == "__main__":
    GlobalConstants.setVerbose(VerboseLevel.STD)

    model = Network("model")

    # Nodes
    source = Source(model, "Source")
    sink = Sink(model, "Sink")
    P1 = Place(model, "P1")
    T1 = Transition(model, "T1")

    # Source arrival
    jobclass = OpenClass(model, "Class1")
    source.setArrival(jobclass, Exp(0.5))  # arrival rate 0.5

    # T1 with Pareto service time
    # Pareto(shape, scale) - shape must be >= 2
    # With shape=3 and scale=1, mean service time = 3*1/(3-1) = 1.5
    mode = T1.addMode("Mode1")
    T1.setNumberOfServers(mode, 1)
    T1.setDistribution(mode, Pareto(3, 1))  # Pareto with shape=3, scale=1
    T1.setEnablingConditions(mode, jobclass, P1, 1)
    T1.setFiringOutcome(mode, jobclass, sink, 1)

    # Routing
    model.addLink(source, P1)
    model.addLink(P1, T1)
    model.addLink(T1, sink)

    # Solver
    solver = JMT(model, seed=23000, samples=10000)
    table = solver.avg_table()
    print(table)
