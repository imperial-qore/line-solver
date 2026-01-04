"""
Closed SPN with Four Places and Different Distributions

This example demonstrates:
- Closed SPN with 4 places in a serial chain
- Different distributions: Exp, Erlang, HyperExp, Coxian
- Tokens move P1 → P2 → P3 → P4 → P1
"""

from line_solver import *
import numpy as np


def build_model() -> Network:
    """Build and return the SPN model."""
    model = Network('model')

    P1 = Place(model, 'P1')
    P2 = Place(model, 'P2')
    P3 = Place(model, 'P3')
    P4 = Place(model, 'P4')

    T1 = Transition(model, 'T1')
    T2 = Transition(model, 'T2')
    T3 = Transition(model, 'T3')
    T4 = Transition(model, 'T4')

    jobclass = ClosedClass(model, 'Class1', 2, P1, 0)

    # T1 - Exponential distribution
    mode1 = T1.add_mode('Mode1')
    T1.set_distribution(mode1, Exp(2.0))
    T1.set_enabling_conditions(mode1, jobclass, P1, 2)
    T1.set_firing_outcome(mode1, jobclass, P2, 2)

    # T2 - Erlang distribution
    mode2 = T2.add_mode('Mode2')
    T2.set_distribution(mode2, Erlang(3, 4))
    T2.set_enabling_conditions(mode2, jobclass, P2, 2)
    T2.set_firing_outcome(mode2, jobclass, P3, 2)

    # T3 - HyperExp distribution
    mode3 = T3.add_mode('Mode3')
    T3.set_distribution(mode3, HyperExp(0.7, 3, 1.5))
    T3.set_enabling_conditions(mode3, jobclass, P3, 2)
    T3.set_firing_outcome(mode3, jobclass, P4, 2)

    # T4 - Coxian distribution
    mu0 = np.array([[1.0], [2.0]])
    phi0 = np.array([[0.6], [1.0]])
    mode4 = T4.add_mode('Mode4')
    T4.set_distribution(mode4, Coxian(mu0, phi0))
    T4.set_enabling_conditions(mode4, jobclass, P4, 2)
    T4.set_firing_outcome(mode4, jobclass, P1, 2)

    # Routing
    R = model.init_routing_matrix()
    R.set(jobclass, jobclass, P1, T1, 1.0)
    R.set(jobclass, jobclass, P2, T2, 1.0)
    R.set(jobclass, jobclass, P3, T3, 1.0)
    R.set(jobclass, jobclass, P4, T4, 1.0)
    R.set(jobclass, jobclass, T1, P2, 1.0)
    R.set(jobclass, jobclass, T2, P3, 1.0)
    R.set(jobclass, jobclass, T3, P4, 1.0)
    R.set(jobclass, jobclass, T4, P1, 1.0)
    model.link(R)

    # Set initial state
    P1.set_state([jobclass.get_population()])
    P2.set_state([0])
    P3.set_state([0])
    P4.set_state([0])

    return model


if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.STD)

    model = build_model()

    solver = JMT(model, seed=23000)
    avg_table = solver.avg_table()
    print(avg_table)
