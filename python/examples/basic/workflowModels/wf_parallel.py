"""
Parallel Workflow Example (AND-fork/join)
A workflow with parallel activities that synchronize.

Workflow:
      +-> B --+
  A --|       |--> D
      +-> C --+

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""
from line_solver.lang.workflow import Workflow
from line_solver.distributions import Exp


def wf_parallel():
    # Define the workflow
    wf = Workflow('ParallelWorkflow')

    # Add activities
    A = wf.addActivity('A', Exp.fitMean(1.0))
    B = wf.addActivity('B', Exp.fitMean(2.0))
    C = wf.addActivity('C', Exp.fitMean(3.0))
    D = wf.addActivity('D', Exp.fitMean(0.5))

    # Define precedences
    wf.addPrecedence(Workflow.AndFork(A, [B, C]))
    wf.addPrecedence(Workflow.AndJoin([B, C], D))

    return wf


if __name__ == '__main__':
    wf = wf_parallel()

    # Convert to phase-type distribution
    alpha, T = wf.toPH()

    # Display results
    print('Parallel Workflow: A -> [B || C] -> D')
    print('Activity means: A=1.0, B=2.0, C=3.0, D=0.5')

    # Expected max of two exponentials: E[max(X,Y)] = 1/lambda1 + 1/lambda2 - 1/(lambda1+lambda2)
    lambda_B = 0.5   # rate for mean 2.0
    lambda_C = 1/3   # rate for mean 3.0
    expected_max = 2.0 + 3.0 - 1/(lambda_B + lambda_C)
    print(f'Expected max(B,C) mean: {expected_max:.4f}')
    print(f'Expected total mean: {1.0 + expected_max + 0.5:.4f}')
    print(f'Computed PH mean: {wf.getMean():.4f}')
    print(f'Number of phases: {T.shape[0]}')
