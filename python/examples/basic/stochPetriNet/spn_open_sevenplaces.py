"""
Open SPN with Seven Places and Immediate Transitions

This example demonstrates:
- Open SPN with 7 places and 8 transitions
- Mix of timed and immediate transitions
- Inhibiting conditions on immediate transitions
- Multiple enabling and firing conditions
"""

from line_solver import *


def build_model() -> Network:
    """Build and return the SPN model."""
    model = Network('model')

    source = Source(model, 'Source')
    sink = Sink(model, 'Sink')

    # Create 7 places
    P = [Place(model, f'P{i+1}') for i in range(7)]

    # Create 8 transitions
    T = [Transition(model, f'T{i+1}') for i in range(8)]

    jobclass = OpenClass(model, 'Class1', 0)
    source.set_arrival(jobclass, Exp.fit_mean(1.0))

    # T1 - Timed transition with exponential distribution
    mode1 = T[0].add_mode('Mode1')
    T[0].set_number_of_servers(mode1, GlobalConstants.MaxInt)
    T[0].set_distribution(mode1, Exp(4.0))
    T[0].set_enabling_conditions(mode1, jobclass, P[0], 1)
    T[0].set_firing_outcome(mode1, jobclass, P[1], 1)

    # T2 - Immediate transition
    mode2 = T[1].add_mode('Mode1')
    T[1].set_number_of_servers(mode2, GlobalConstants.MaxInt)
    T[1].set_enabling_conditions(mode2, jobclass, P[1], 1)
    T[1].set_firing_outcome(mode2, jobclass, P[2], 1)
    T[1].setTimingStrategy(mode2, TimingStrategy.IMMEDIATE)
    T[1].setFiringPriorities(mode2, 1)
    T[1].setFiringWeights(mode2, 1.0)

    # T3 - Immediate transition
    mode3 = T[2].add_mode('Mode1')
    T[2].set_number_of_servers(mode3, GlobalConstants.MaxInt)
    T[2].set_enabling_conditions(mode3, jobclass, P[1], 1)
    T[2].set_firing_outcome(mode3, jobclass, P[3], 1)
    T[2].setTimingStrategy(mode3, TimingStrategy.IMMEDIATE)
    T[2].setFiringPriorities(mode3, 1)

    # T4 - Immediate transition with multiple enabling/firing conditions
    mode4 = T[3].add_mode('Mode1')
    T[3].set_number_of_servers(mode4, GlobalConstants.MaxInt)
    T[3].set_enabling_conditions(mode4, jobclass, P[2], 1)
    T[3].set_enabling_conditions(mode4, jobclass, P[4], 1)
    T[3].set_firing_outcome(mode4, jobclass, P[4], 1)
    T[3].set_firing_outcome(mode4, jobclass, P[5], 1)
    T[3].setTimingStrategy(mode4, TimingStrategy.IMMEDIATE)
    T[3].setFiringPriorities(mode4, 1)

    # T5 - Immediate transition with inhibiting conditions
    mode5 = T[4].add_mode('Mode1')
    T[4].set_number_of_servers(mode5, GlobalConstants.MaxInt)
    T[4].set_enabling_conditions(mode5, jobclass, P[3], 1)
    T[4].set_enabling_conditions(mode5, jobclass, P[4], 1)
    T[4].set_firing_outcome(mode5, jobclass, P[6], 1)
    T[4].setInhibitingConditions(mode5, jobclass, P[5], 1)
    T[4].setTimingStrategy(mode5, TimingStrategy.IMMEDIATE)
    T[4].setFiringPriorities(mode5, 1)

    # T6 - Timed transition with Erlang distribution
    mode6 = T[5].add_mode('Mode1')
    T[5].set_number_of_servers(mode6, GlobalConstants.MaxInt)
    T[5].set_distribution(mode6, Erlang(2, 2))
    T[5].set_enabling_conditions(mode6, jobclass, P[5], 1)
    T[5].set_firing_outcome(mode6, jobclass, P[0], 1)

    # T7 - Timed transition with multiple firing outcomes
    mode7 = T[6].add_mode('Mode1')
    T[6].set_number_of_servers(mode7, GlobalConstants.MaxInt)
    T[6].set_distribution(mode7, Exp(2.0))
    T[6].set_enabling_conditions(mode7, jobclass, P[6], 1)
    T[6].set_firing_outcome(mode7, jobclass, P[0], 1)
    T[6].set_firing_outcome(mode7, jobclass, P[4], 1)

    # T8 - Timed transition to sink
    mode8 = T[7].add_mode('Mode1')
    T[7].set_number_of_servers(mode8, GlobalConstants.MaxInt)
    T[7].set_distribution(mode8, Exp(2.0))
    T[7].set_enabling_conditions(mode8, jobclass, P[3], 1)
    T[7].set_firing_outcome(mode8, jobclass, sink, 1)

    # Routing
    R = model.init_routing_matrix()

    # Source routing
    R.set(jobclass, jobclass, source, P[0], 1.0)

    # Place to transition routing
    R.set(jobclass, jobclass, P[0], T[0], 1.0)
    R.set(jobclass, jobclass, P[1], T[1], 1.0)
    R.set(jobclass, jobclass, P[1], T[2], 1.0)
    R.set(jobclass, jobclass, P[2], T[3], 1.0)
    R.set(jobclass, jobclass, P[3], T[4], 1.0)
    R.set(jobclass, jobclass, P[4], T[3], 1.0)
    R.set(jobclass, jobclass, P[4], T[4], 1.0)
    R.set(jobclass, jobclass, P[5], T[4], 1.0)
    R.set(jobclass, jobclass, P[5], T[5], 1.0)
    R.set(jobclass, jobclass, P[6], T[6], 1.0)
    R.set(jobclass, jobclass, P[3], T[7], 1.0)

    # Transition to place/sink routing
    R.set(jobclass, jobclass, T[0], P[1], 1.0)
    R.set(jobclass, jobclass, T[1], P[2], 1.0)
    R.set(jobclass, jobclass, T[2], P[3], 1.0)
    R.set(jobclass, jobclass, T[3], P[4], 1.0)
    R.set(jobclass, jobclass, T[3], P[5], 1.0)
    R.set(jobclass, jobclass, T[4], P[6], 1.0)
    R.set(jobclass, jobclass, T[5], P[0], 1.0)
    R.set(jobclass, jobclass, T[6], sink, 1.0)
    R.set(jobclass, jobclass, T[6], P[0], 1.0)
    R.set(jobclass, jobclass, T[6], P[4], 1.0)
    R.set(jobclass, jobclass, T[7], sink, 1.0)
    model.link(R)

    # Set initial state
    P[0].set_state([2])
    P[1].set_state([0])
    P[2].set_state([0])
    P[3].set_state([0])
    P[4].set_state([1])
    P[5].set_state([0])
    P[6].set_state([0])

    return model


if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.STD)

    model = build_model()

    solver = JMT(model, seed=23000)
    avg_table = solver.avg_table()
    print(avg_table)
