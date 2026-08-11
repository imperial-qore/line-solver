"""
Serial Workflow Example
A simple workflow with three sequential activities.

Workflow: A -> B -> C

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

import sys
sys.path.insert(0, '/home/gcasale/Dropbox/code/line-dev.git/python')

from line_solver.lang.workflow import Workflow
from line_solver.distributions import Exp


def wf_serial():
    # Define the workflow
    wf = Workflow('SerialWorkflow')

    # Add activities with exponential service times
    A = wf.addActivity('A', Exp.fitMean(1.0))
    B = wf.addActivity('B', Exp.fitMean(2.0))
    C = wf.addActivity('C', Exp.fitMean(1.5))

    # Define serial precedence
    wf.addPrecedence(Workflow.Serial(A, B))
    wf.addPrecedence(Workflow.Serial(B, C))

    return wf


if __name__ == '__main__':
    wf = wf_serial()

    # Convert to phase-type distribution
    alpha, T = wf.toPH()

    # Display results
    print('Serial Workflow: A -> B -> C')
    print('Activity means: A=1.0, B=2.0, C=1.5')
    print(f'Expected total mean: {1.0 + 2.0 + 1.5:.2f}')
    print(f'Computed PH mean: {wf.getMean():.4f}')
    print(f'Computed PH SCV: {wf.getSCV():.4f}')
    print(f'Number of phases: {T.shape[0]}')
