"""
Complex Workflow Example
A workflow combining serial, parallel, and branching patterns.

Workflow:
                  +-> C --+
  A -> B -> [AND] |       | [AND] -> F -> G
                  +-> D --+

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""
from line_solver.lang.workflow import Workflow
from line_solver.distributions import Exp


def wf_complex():
    # Define the workflow
    wf = Workflow('ComplexWorkflow')

    # Add activities
    A = wf.addActivity('A', Exp.fitMean(0.5))
    B = wf.addActivity('B', Exp.fitMean(1.0))
    C = wf.addActivity('C', Exp.fitMean(2.0))
    D = wf.addActivity('D', Exp.fitMean(1.5))
    F = wf.addActivity('F', Exp.fitMean(1.0))
    G = wf.addActivity('G', Exp.fitMean(0.5))

    # Define precedences
    wf.addPrecedence(Workflow.Serial(A, B))       # A -> B
    wf.addPrecedence(Workflow.AndFork(B, [C, D])) # B forks to C and D
    wf.addPrecedence(Workflow.AndJoin([C, D], F)) # C and D join to F
    wf.addPrecedence(Workflow.Serial(F, G))       # F -> G

    return wf


if __name__ == '__main__':
    wf = wf_complex()

    # Convert to phase-type distribution
    alpha, T = wf.toPH()

    # Display results
    print('Complex Workflow: A -> B -> [C || D] -> F -> G')
    print('Activity means: A=0.5, B=1.0, C=2.0, D=1.5, F=1.0, G=0.5')

    # Calculate expected max(C, D)
    lambda_C = 0.5    # rate for mean 2.0
    lambda_D = 1/1.5  # rate for mean 1.5
    expected_max = 2.0 + 1.5 - 1/(lambda_C + lambda_D)
    print(f'Expected max(C,D) mean: {expected_max:.4f}')
    print(f'Expected total mean: {0.5 + 1.0 + expected_max + 1.0 + 0.5:.4f}')
    print(f'Computed PH mean: {wf.getMean():.4f}')
    print(f'Number of phases: {T.shape[0]}')

    # Sample from the distribution
    import numpy as np
    np.random.seed(1241)
    print('\nSampling 10000 workflow execution times...')
    samples = wf.sample(10000)
    print(f'Sample mean: {np.mean(samples):.4f}')
    print(f'Sample std: {np.std(samples):.4f}')
    print(f'Sample min: {np.min(samples):.4f}')
    print(f'Sample max: {np.max(samples):.4f}')
