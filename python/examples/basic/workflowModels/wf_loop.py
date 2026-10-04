"""
Loop Workflow Example
A workflow with a repeated activity.

Workflow: A -> [B x 3] -> C
Activity B is executed 3 times.

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""
from line_solver.lang.workflow import Workflow
from line_solver.distributions import Exp


def wf_loop():
    # Define the workflow
    wf = Workflow('LoopWorkflow')

    # Add activities
    A = wf.addActivity('A', Exp.fitMean(1.0))
    B = wf.addActivity('B', Exp.fitMean(2.0))
    C = wf.addActivity('C', Exp.fitMean(0.5))

    # Define loop: A runs once, B runs 3 times, C runs once
    wf.addPrecedence(Workflow.Loop(A, [B], C, count=3))

    return wf


if __name__ == '__main__':
    wf = wf_loop()

    # Convert to phase-type distribution
    alpha, T = wf.toPH()

    # Display results
    print('Loop Workflow: A -> [B x 3] -> C')
    print('Activity means: A=1.0, B=2.0, C=0.5')
    print(f'Expected total mean: 1.0 + 3*2.0 + 0.5 = {1.0 + 3*2.0 + 0.5:.2f}')
    print(f'Computed PH mean: {wf.getMean():.4f}')
    print(f'Number of phases: {T.shape[0]}')
