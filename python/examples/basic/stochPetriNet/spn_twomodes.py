"""
Closed SPN with Two Transitions (Two Modes)

This example demonstrates:
- Closed SPN with 2 places and 2 transitions
- Tokens move P1 → P2 → P1 in batches
- T1 requires 4 tokens, T2 requires 2 tokens
"""

from line_solver import *


def build_model() -> Network:
    """Build and return the SPN model."""
    model = Network('model')

    P1 = Place(model, 'P1')
    P2 = Place(model, 'P2')
    T1 = Transition(model, 'T1')
    T2 = Transition(model, 'T2')

    jobclass = ClosedClass(model, 'Class1', 10, P1, 0)

    # T1: requires 4 tokens from P1, produces 4 tokens to P2
    mode1 = T1.add_mode('Mode1')
    T1.set_distribution(mode1, Exp(2.0))
    T1.set_enabling_conditions(mode1, jobclass, P1, 4)
    T1.set_firing_outcome(mode1, jobclass, P2, 4)

    # T2: requires 2 tokens from P2, produces 2 tokens to P1
    mode2 = T2.add_mode('Mode2')
    T2.set_distribution(mode2, Exp(3.0))
    T2.set_enabling_conditions(mode2, jobclass, P2, 2)
    T2.set_firing_outcome(mode2, jobclass, P1, 2)

    # Routing
    R = model.init_routing_matrix()
    R.set(jobclass, jobclass, P1, T1, 1.0)
    R.set(jobclass, jobclass, P2, T2, 1.0)
    R.set(jobclass, jobclass, T1, P2, 1.0)
    R.set(jobclass, jobclass, T2, P1, 1.0)
    model.link(R)

    return model


if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.STD)

    model = build_model()

    solver = JMT(model, seed=23000)
    avg_table = solver.avg_table()
    print(avg_table)
