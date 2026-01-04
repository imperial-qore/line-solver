#!/usr/bin/env python3
"""
Example 11: Random environments and SolverENV

This tutorial illustrates how to model a queueing system operating in
a random environment, where system parameters (e.g., service rates)
change according to an underlying environmental process.

Scenario: A server that alternates between "Fast" and "Slow" modes.
In Fast mode, service rate is 4.0. In Slow mode, service rate is 1.0.
The environment switches from Fast->Slow at rate 0.5 and Slow->Fast at rate 1.0.
"""

import sys
import os

# Add the line_solver package to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from line_solver import Network, ClosedClass, Delay, Queue
from line_solver import Environment, ENV, FLD, MVA
from line_solver.constants import SchedStrategy
from line_solver.distributions import Exp


def main():
    # Block 1: Create base network model
    # First, we define a closed queueing network with a delay and a queue
    base_model = Network('BaseModel')
    delay = Delay(base_model, 'ThinkTime')
    queue = Queue(base_model, 'Fast/Slow Server', SchedStrategy.FCFS)

    # Closed class with 5 jobs
    N = 5
    jobclass = ClosedClass(base_model, 'Jobs', N, delay)
    delay.setService(jobclass, Exp(1.0))  # Think time = 1.0
    queue.setService(jobclass, Exp(2.0))  # Placeholder service rate

    # Connect nodes in a cycle
    base_model.link(Network.serialRouting(delay, queue))

    # Block 2: Create the random environment
    # Define two stages: "Fast" and "Slow" with different service rates
    env = Environment('ServerModes', 2)

    # Stage 0: Fast mode (service rate = 4.0)
    fast_model = base_model.copy()
    fast_queue = fast_model.getNodeByName('Fast/Slow Server')
    fast_queue.setService(fast_model.classes[0], Exp(4.0))
    env.addStage(0, 'Fast', 'operational', fast_model)

    # Stage 1: Slow mode (service rate = 1.0)
    slow_model = base_model.copy()
    slow_queue = slow_model.getNodeByName('Fast/Slow Server')
    slow_queue.setService(slow_model.classes[0], Exp(1.0))
    env.addStage(1, 'Slow', 'degraded', slow_model)

    # Define transitions between stages
    # Fast -> Slow at rate 0.5 (mean time in Fast mode = 2.0)
    env.addTransition(0, 1, Exp(0.5))
    # Slow -> Fast at rate 1.0 (mean time in Slow mode = 1.0)
    env.addTransition(1, 0, Exp(1.0))

    # Block 3: Inspect the environment structure
    print('Environment stages:')
    stage_table = env.getStageTable()
    print(stage_table)

    # Block 4: Solve using SolverENV
    # SolverENV requires a solver factory that creates solvers for each stage
    # We use the Fluid solver (FLD) with transient analysis
    env_solver = ENV(env, lambda m: FLD(m))
    Q, U, T = env_solver.getAvg()

    # Display average results weighted by environment probabilities
    print('\n--- Environment-Averaged Results ---')
    env_avg_table = env_solver.getAvgTable()
    print(env_avg_table)

    # Block 5: Compare with individual stage analysis
    # Analyze each stage network in steady-state using MVA
    print('\n--- Individual Stage Analysis (MVA) ---')
    ensemble = env.getEnsemble()
    for e in range(len(ensemble)):
        stage_model = ensemble[e]
        print(f'\nStage {e}:')
        mva_solver = MVA(stage_model)
        stage_table = mva_solver.getAvgTable()
        print(stage_table)


if __name__ == '__main__':
    main()
