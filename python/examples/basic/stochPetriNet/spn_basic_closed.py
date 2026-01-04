"""
Basic Closed Stochastic Petri Net

This example demonstrates:
- Closed SPN: P1 → T1 → P1 (single place loop)
- Single transition T1 with 3 racing modes
- Mode 1: Exponential, Mode 2: Erlang, Mode 3: HyperExp
- 1 token (job) circulates in the system
"""

from line_solver import *


def build_model() -> Network:
    """Build and return the SPN model."""
    model = Network('model')

    P1 = Place(model, 'P1')
    T1 = Transition(model, 'T1')

    jobclass = ClosedClass(model, 'Class1', 1, P1)

    # T1 - Mode 1: Exponential
    mode1 = T1.add_mode('Mode1')
    T1.set_distribution(mode1, Exp.fit_mean(1.0))
    T1.set_enabling_conditions(mode1, jobclass, P1, 1)
    T1.set_firing_outcome(mode1, jobclass, P1, 1)

    # T1 - Mode 2: Erlang
    mode2 = T1.add_mode('Mode2')
    T1.set_distribution(mode2, Erlang.fit_mean_and_order(1.0, 2))
    T1.set_enabling_conditions(mode2, jobclass, P1, 1)
    T1.set_firing_outcome(mode2, jobclass, P1, 1)

    # T1 - Mode 3: HyperExp
    mode3 = T1.add_mode('Mode3')
    T1.set_distribution(mode3, HyperExp.fit_mean_and_scv(1.0, 4.0))
    T1.set_enabling_conditions(mode3, jobclass, P1, 1)
    T1.set_firing_outcome(mode3, jobclass, P1, 1)

    # Routing
    R = model.init_routing_matrix()
    R.set(jobclass, jobclass, P1, T1, 1.0)
    R.set(jobclass, jobclass, T1, P1, 1.0)
    model.link(R)

    # Set initial state
    P1.set_state([jobclass.get_population()])

    return model


if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.STD)

    model = build_model()

    solver = JMT(model, seed=23000)
    avg_table = solver.avg_table()
    print(avg_table)
