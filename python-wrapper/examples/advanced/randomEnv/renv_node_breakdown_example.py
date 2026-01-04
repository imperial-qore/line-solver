"""
Example demonstrating the node breakdown/repair API for random environments.

This example shows how to easily model a server that can break down and be
repaired using the new Env convenience methods.

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

from line_solver import Network, Environment, Source, Queue, Sink, OpenClass, SchedStrategy
from line_solver import Exp, ENV, FLD, SolverOptions, SolverType, VerboseLevel
import numpy as np


def create_base_model():
    """Create a base queueing network model with a single server.

    Returns:
        tuple: (model, queue) - the network model and the queue node for use with Env API
    """
    model = Network('ServerWithFailures')

    # Define nodes
    source = Source(model, 'Arrivals')
    queue = Queue(model, 'Server', SchedStrategy.FCFS)
    sink = Sink(model, 'Departures')

    # Define job class
    jobclass = OpenClass(model, 'Jobs')

    # Set service and arrival rates (UP state)
    source.set_arrival(jobclass, Exp(0.8))  # Arrival rate
    queue.set_service(jobclass, Exp(2.0))   # Service rate when UP
    queue.set_number_of_servers(1)

    # Set routing
    P = model.init_routing_matrix()
    P.add_connection(jobclass, jobclass, source, queue)
    P.add_connection(jobclass, jobclass, queue, sink)
    model.link(P)

    return model, queue


def example1_basic_failure_repair():
    """Example 1: Using addNodeFailureRepair convenience method with node object (recommended)."""
    print("=" * 70)
    print("Example 1: Using addNodeFailureRepair with node object")
    print("=" * 70)

    model, queue = create_base_model()

    # Create environment with 2 stages (UP and DOWN)
    env = Environment('ServerEnv1', 2)

    # Add failure and repair for the server node in one call
    # Parameters: base_model, node_or_name, breakdown_dist, repair_dist, down_service_dist
    # Note: Can pass either node object (queue) or node name ('Server')
    env.add_node_failure_repair(model, queue, Exp(0.1), Exp(1.0), Exp(0.5))

    # Initialize and print stage table
    env.obj.init()
    print("\nStage table for env1:")
    env.get_stage_table()
    print()


def example2_separate_calls():
    """Example 2: Using separate breakdown and repair calls with node object."""
    print("=" * 70)
    print("Example 2: Using separate breakdown and repair calls with node object")
    print("=" * 70)

    model, queue = create_base_model()

    env = Environment('ServerEnv2', 2)

    # Add breakdown (creates UP and DOWN stages) - using node object
    env.add_node_breakdown(model, queue, Exp(0.1), Exp(0.5))

    # Add repair transition - using node object
    env.add_node_repair(queue, Exp(1.0))

    # Initialize and print stage table
    env.obj.init()
    print("\nStage table for env2:")
    env.get_stage_table()
    print()


def example3_custom_reset_policies():
    """Example 3: With custom reset policies using node object."""
    print("=" * 70)
    print("Example 3: With custom reset policies using node object")
    print("=" * 70)

    model, queue = create_base_model()

    env = Environment('ServerEnv3', 2)

    # Reset policy: clear all queues on breakdown
    def reset_breakdown(q):
        """Clear all queues when server breaks down."""
        return np.zeros_like(q)

    # Reset policy: keep all jobs on repair
    def reset_repair(q):
        """Keep jobs when server is repaired."""
        return q

    # Using node object instead of string name
    env.add_node_failure_repair(model, queue, Exp(0.1), Exp(1.0), Exp(0.5),
                                reset_breakdown, reset_repair)

    # Initialize and print stage table
    env.obj.init()
    print("\nStage table for env3:")
    env.get_stage_table()
    print()


def example4_modify_reset_policies():
    """Example 4: Modifying reset policies after creation using node object."""
    print("=" * 70)
    print("Example 4: Modifying reset policies after creation using node object")
    print("=" * 70)

    model, queue = create_base_model()

    env = Environment('ServerEnv4', 2)

    # Create environment with default reset policies - using node object
    env.add_node_failure_repair(model, queue, Exp(0.1), Exp(1.0), Exp(0.5))

    # Update breakdown reset policy to clear queues - using node object
    env.set_breakdown_reset_policy(queue, lambda q: np.zeros_like(q))

    # Update repair reset policy (keep jobs) - using node object
    env.set_repair_reset_policy(queue, lambda q: q)

    # Initialize and print stage table
    env.obj.init()
    print("\nStage table for env4:")
    env.get_stage_table()
    print()


def example5_solve_environment():
    """Example 5: Solve the environment model and display results."""
    print("=" * 70)
    print("Example 5: Solving environment model with ENV solver")
    print("=" * 70)

    model, queue = create_base_model()

    env = Environment('ServerEnv5', 2)
    # Using node object
    env.add_node_failure_repair(model, queue, Exp(0.1), Exp(1.0), Exp(0.5))
    env.obj.init()

    # Create solver options
    options = SolverOptions(SolverType.ENV)
    options.iter_tol = 0.01
    options.iter_max = 100
    options.timespan = [0, float('inf')]
    options.verbose = VerboseLevel.STD

    # Create fluid solver options for each stage
    fluid_options = SolverOptions(SolverType.FLUID)
    fluid_options.timespan = [0, 1000]
    fluid_options.stiff = False
    fluid_options.set_ode_max_step(0.25)

    # Create solvers for each stage
    num_stages = len(env.get_ensemble())
    solvers = [FLD(env.get_model(e), fluid_options) for e in range(num_stages)]

    # Create and run ENV solver
    solver = ENV(env, solvers, options)

    print("\nAverage Performance Metrics:")
    avg_table = solver.get_avg_table()
    print(avg_table)

    print("\nInterpretation:")
    print("- The system alternates between UP (operational) and DOWN (failed) states")
    print("- UP state: Server processes jobs at rate 2.0")
    print("- DOWN state: Server processes jobs at reduced rate 0.5")
    print("- Breakdown occurs at rate 0.1 (mean time to failure = 10 time units)")
    print("- Repair occurs at rate 1.0 (mean time to repair = 1 time unit)")
    print("- Results show averaged performance across both states")
    print()


def main():
    """Run all examples."""
    try:
        example1_basic_failure_repair()
        example2_separate_calls()
        example3_custom_reset_policies()
        example4_modify_reset_policies()
        example5_solve_environment()

        print("=" * 70)
        print("All examples completed successfully!")
        print("=" * 70)
    except Exception as e:
        print(f"Error running examples: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
