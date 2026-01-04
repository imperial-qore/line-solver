"""
Example: Workflow Scheduling with Hora Scheduler

This example demonstrates how to use the Hora Scheduler module
to schedule DAG-based workflows on heterogeneous machines.
"""

import sys
import os
import numpy as np

# Add parent directory to path for direct execution
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scheduling import (
    Workflow, Task,
    HeterogeneousMachines, IdenticalMachines,
    HEFT
)


def example_programmatic():
    """Example: Programmatic workflow declaration."""
    print("=" * 60)
    print("Example 1: Programmatic Workflow Declaration")
    print("=" * 60)

    # Create a workflow representing a data processing pipeline
    wf = Workflow('DataPipeline')

    # Add tasks with computation weights
    load = Task(wf, 'Load', computation=10.0)
    filter1 = Task(wf, 'Filter1', computation=15.0)
    filter2 = Task(wf, 'Filter2', computation=12.0)
    merge = Task(wf, 'Merge', computation=8.0)
    save = Task(wf, 'Save', computation=5.0)

    # Add dependencies with communication costs
    wf.addDependency(load, filter1, comm_cost=3.0)
    wf.addDependency(load, filter2, comm_cost=3.0)
    wf.addDependency(filter1, merge, comm_cost=2.0)
    wf.addDependency(filter2, merge, comm_cost=2.0)
    wf.addDependency(merge, save, comm_cost=1.0)

    print(f"Workflow: {wf}")
    print(f"Tasks: {[t.name for t in wf.getTasks()]}")
    print(f"Entry tasks: {[t.name for t in wf.getEntryTasks()]}")
    print(f"Exit tasks: {[t.name for t in wf.getExitTasks()]}")

    # Create identical machines
    machines = IdenticalMachines('Cluster', num_machines=3)
    print(f"\nMachines: {machines}")

    # Schedule with HEFT
    solver = HEFT(wf, machines, verbose=True)
    result = solver.schedule()

    print(f"\nSchedule Table:")
    print(result.schedule_table)
    print(f"\nMakespan: {result.makespan:.2f}")
    print(f"Speedup: {result.speedup:.2f}")
    print(f"Efficiency: {result.efficiency:.2f}")

    return result


def example_heterogeneous():
    """Example: Heterogeneous machines with different speeds."""
    print("\n" + "=" * 60)
    print("Example 2: Heterogeneous Machines")
    print("=" * 60)

    # Create workflow
    wf = Workflow('GPUWorkload')
    t1 = Task(wf, 'Preprocess', computation=20.0)
    t2 = Task(wf, 'GPUCompute', computation=100.0)
    t3 = Task(wf, 'Postprocess', computation=10.0)

    wf.addDependency(t1, t2, comm_cost=5.0)
    wf.addDependency(t2, t3, comm_cost=5.0)

    # Create heterogeneous machines
    machines = HeterogeneousMachines('HybridCluster')
    gpu = machines.addMachine('GPU', speed=4.0)    # 4x faster
    cpu1 = machines.addMachine('CPU1', speed=1.0)
    cpu2 = machines.addMachine('CPU2', speed=1.0)

    # Set communication speeds
    machines.setCommSpeed(gpu, cpu1, 2.0)   # Fast GPU-CPU link
    machines.setCommSpeed(gpu, cpu2, 2.0)
    machines.setCommSpeed(cpu1, cpu2, 1.0)  # Slower CPU-CPU

    print(f"Machines: {machines}")
    print(f"Speeds: {machines.getMachineSpeeds()}")

    # Schedule
    result = HEFT(wf, machines, verbose=True).schedule()

    print(f"\nSchedule Table:")
    print(result.schedule_table)
    print(f"\nMakespan: {result.makespan:.2f}")

    return result


def example_dynamic():
    """Example: Dynamic scheduling with release times."""
    print("\n" + "=" * 60)
    print("Example 3: Dynamic Scheduling with Release Times")
    print("=" * 60)

    wf = Workflow('DynamicWorkflow')

    # Tasks arrive at different times
    t1 = Task(wf, 'Task1', computation=10.0, release_time=0.0)
    t2 = Task(wf, 'Task2', computation=8.0, release_time=5.0)   # Arrives at t=5
    t3 = Task(wf, 'Task3', computation=6.0, release_time=12.0)  # Arrives at t=12

    wf.addDependency(t1, t3, comm_cost=1.0)
    wf.addDependency(t2, t3, comm_cost=1.0)

    machines = IdenticalMachines('Cluster', num_machines=2)

    result = HEFT(wf, machines, verbose=True).schedule()

    print(f"\nSchedule Table:")
    print(result.schedule_table)
    print(f"\nNote: Task2 respects release_time=5.0")
    print(f"Note: Task3 respects release_time=12.0")

    return result


def example_gml_import():
    """Example: Loading workflow from GML file."""
    print("\n" + "=" * 60)
    print("Example 4: Loading from GML File (ISO25 format)")
    print("=" * 60)

    gml_path = '/home/gcasale/Dropbox/code/iso25-wflowsched-wenqi.git/static/dag/10/Tau_0.gml'

    if not os.path.exists(gml_path):
        print(f"GML file not found: {gml_path}")
        return None

    # Load workflow from GML
    wf = Workflow.from_gml(gml_path)
    print(f"Loaded workflow: {wf}")
    print(f"Number of tasks: {wf.getNumberOfTasks()}")

    # Create system matching ISO25 configuration
    speeds = np.array([0.5, 1.0, 0.25])
    comm_matrix = np.array([
        [np.inf, 1, 1.5],
        [1, np.inf, 2],
        [1.5, 2, np.inf]
    ])
    machines = HeterogeneousMachines.from_arrays('ISO25System', speeds, comm_matrix)

    print(f"Machines: {machines}")
    print(f"Machine speeds: {machines.getMachineSpeeds()}")

    # Schedule
    result = HEFT(wf, machines, verbose=True).schedule()

    print(f"\nSchedule Table:")
    print(result.schedule_table)
    print(f"\nMakespan: {result.makespan:.2f}")
    print(f"SLR: {result.slr:.2f}")

    return result


if __name__ == '__main__':
    # Run all examples
    example_programmatic()
    example_heterogeneous()
    example_dynamic()
    example_gml_import()

    print("\n" + "=" * 60)
    print("All examples completed successfully!")
    print("=" * 60)
