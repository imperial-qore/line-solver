"""
Closed SPN with Two Places and Two Classes

This example demonstrates:
- Closed SPN with 2 places (P1, P2) and 3 transitions
- Two job classes with different populations
- Multi-token enabling conditions and firing outcomes
"""

from line_solver import *


def build_model() -> Network:
    """Build and return the SPN model."""
    model = Network('model')

    P1 = Place(model, 'P1')
    P2 = Place(model, 'P2')

    T1 = Transition(model, 'T1')
    T2 = Transition(model, 'T2')
    T3 = Transition(model, 'T3')

    jobclass1 = ClosedClass(model, 'Class1', 10, P1, 0)
    jobclass2 = ClosedClass(model, 'Class2', 7, P1, 0)

    # T1 - Mode 1: Class1, requires 2 tokens, produces 2 tokens
    mode1 = T1.add_mode('Mode1')
    T1.set_distribution(mode1, Exp(2.0))
    T1.set_enabling_conditions(mode1, jobclass1, P1, 2)
    T1.set_firing_outcome(mode1, jobclass1, P2, 2)

    # T1 - Mode 2: Class2, requires 1 token, produces 1 token
    mode2 = T1.add_mode('Mode2')
    T1.set_distribution(mode2, Exp(3.0))
    T1.set_enabling_conditions(mode2, jobclass2, P1, 1)
    T1.set_firing_outcome(mode2, jobclass2, P2, 1)

    # T2 - Mode 3: Class1 with Erlang
    mode3 = T2.add_mode('Mode3')
    T2.set_distribution(mode3, Erlang(1.5, 2))
    T2.set_enabling_conditions(mode3, jobclass1, P2, 1)
    T2.set_firing_outcome(mode3, jobclass1, P1, 1)

    # T3 - Mode 4: Class2, requires 4 tokens, produces 4 tokens
    mode4 = T3.add_mode('Mode4')
    T3.set_distribution(mode4, Exp(0.5))
    T3.set_enabling_conditions(mode4, jobclass2, P2, 4)
    T3.set_firing_outcome(mode4, jobclass2, P1, 4)

    # Routing
    R = model.init_routing_matrix()
    R.set(jobclass1, jobclass1, P1, T1, 1.0)
    R.set(jobclass2, jobclass2, P1, T1, 1.0)
    R.set(jobclass1, jobclass1, P2, T2, 1.0)
    R.set(jobclass2, jobclass2, P2, T3, 1.0)
    R.set(jobclass1, jobclass1, T1, P2, 1.0)
    R.set(jobclass2, jobclass2, T1, P2, 1.0)
    R.set(jobclass1, jobclass1, T2, P1, 1.0)
    R.set(jobclass2, jobclass2, T3, P1, 1.0)
    model.link(R)

    # Set initial state
    P1.set_state([jobclass1.get_population(), jobclass2.get_population()])
    P2.set_state([0, 0])

    return model


if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.STD)

    model = build_model()

    solver = JMT(model, seed=23000)
    avg_table = solver.avg_table()
    print(avg_table)
