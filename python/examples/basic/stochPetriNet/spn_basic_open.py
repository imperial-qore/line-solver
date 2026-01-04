"""
Basic Open Stochastic Petri Net

This example demonstrates:
- Open SPN: Source → P1 → T1 → Sink
- Single transition T1 with exponential firing time
- Infinite servers at transition
"""

from line_solver import *


def build_model() -> Network:
    """Build and return the SPN model."""
    model = Network('model')

    source = Source(model, 'Source')
    sink = Sink(model, 'Sink')
    P1 = Place(model, 'P1')
    T1 = Transition(model, 'T1')

    jobclass = OpenClass(model, 'Class1', 0)
    source.set_arrival(jobclass, Exp(1.0))

    # T1
    mode = T1.add_mode('Mode1')
    T1.set_number_of_servers(mode, GlobalConstants.MaxInt)
    T1.set_distribution(mode, Exp(4.0))
    T1.set_enabling_conditions(mode, jobclass, P1, 1)
    T1.set_firing_outcome(mode, jobclass, sink, 1)

    # Routing
    R = model.init_routing_matrix()
    R.set(jobclass, jobclass, source, P1, 1.0)
    R.set(jobclass, jobclass, P1, T1, 1.0)
    R.set(jobclass, jobclass, T1, sink, 1.0)
    model.link(R)

    return model


if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.STD)

    model = build_model()

    solver = JMT(model, seed=23000)
    avg_table = solver.avg_table()
    print(avg_table)
