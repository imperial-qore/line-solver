"""
Probabilistic Branching Workflow Example (OR-fork/join)
A workflow with probabilistic choice between branches.

Workflow:
      +-> B (60%) --+
  A --|             |--> D
      +-> C (40%) --+

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""
import numpy as np

from line_solver.lang.workflow import Workflow
from line_solver.distributions import Exp


def wf_branch():
    # Define the workflow
    wf = Workflow('BranchingWorkflow')

    # Add activities
    A = wf.addActivity('A', Exp.fitMean(1.0))
    B = wf.addActivity('B', Exp.fitMean(2.0))
    C = wf.addActivity('C', Exp.fitMean(5.0))
    D = wf.addActivity('D', Exp.fitMean(0.5))

    # Define precedences with probabilities
    wf.addPrecedence(Workflow.OrFork(A, [B, C], np.array([0.6, 0.4])))
    wf.addPrecedence(Workflow.OrJoin([B, C], D))

    return wf


if __name__ == '__main__':
    wf = wf_branch()

    # Convert to phase-type distribution
    alpha, T = wf.toPH()

    # Display results
    print('Branching Workflow: A -> [B(60%) | C(40%)] -> D')
    print('Activity means: A=1.0, B=2.0, C=5.0, D=0.5')
    expected_branch = 0.6 * 2.0 + 0.4 * 5.0
    print(f'Expected branch mean: 0.6*2.0 + 0.4*5.0 = {expected_branch:.2f}')
    print(f'Expected total mean: {1.0 + expected_branch + 0.5:.2f}')
    print(f'Computed PH mean: {wf.getMean():.4f}')
    print(f'Number of phases: {T.shape[0]}')
