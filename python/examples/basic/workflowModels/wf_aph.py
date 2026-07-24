"""
Workflow with APH Distributions
A workflow using APH (Acyclic Phase-Type) distributions
for flexible variability modeling.

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""
from line_solver.lang.workflow import Workflow
from line_solver.distributions import APH


def wf_aph():
    # Define the workflow
    wf = Workflow('APHWorkflow')

    # Add activities with different SCVs using APH
    A = wf.addActivity('A', APH.fitMeanAndScv(1.0, 0.5))   # Low variability
    B = wf.addActivity('B', APH.fitMeanAndScv(2.0, 2.0))   # High variability
    C = wf.addActivity('C', APH.fitMeanAndScv(1.5, 1.0))   # Standard variability

    # Define serial precedence
    wf.addPrecedence(Workflow.Serial(A, B))
    wf.addPrecedence(Workflow.Serial(B, C))

    return wf


if __name__ == '__main__':
    wf = wf_aph()

    # Convert to phase-type distribution
    alpha, T = wf.toPH()

    # Display results
    print('APH Workflow: A -> B -> C')
    print('Activity A: APH(mean=1.0, SCV=0.5)')
    print('Activity B: APH(mean=2.0, SCV=2.0)')
    print('Activity C: APH(mean=1.5, SCV=1.0)')
    print(f'Expected total mean: {1.0 + 2.0 + 1.5:.2f}')
    print(f'Computed PH mean: {wf.getMean():.4f}')
    print(f'Computed PH SCV: {wf.getSCV():.4f}')
    print(f'Number of phases: {T.shape[0]}')
