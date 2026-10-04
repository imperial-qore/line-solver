"""
Workflow with Erlang Distribution
A workflow using Erlang distributions for lower variability.

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

from line_solver.lang.workflow import Workflow
from line_solver.distributions import Exp, Erlang


def wf_erlang():

    # Define the workflow
    wf = Workflow('ErlangWorkflow')

    # Add activities with Erlang distributions (SCV < 1)
    A = wf.addActivity('A', Erlang.fitMeanAndOrder(2.0, 4))  # 4 phases, SCV=0.25
    B = wf.addActivity('B', Erlang.fitMeanAndOrder(3.0, 2))  # 2 phases, SCV=0.5
    C = wf.addActivity('C', Exp.fitMean(1.0))                # 1 phase, SCV=1.0

    # Define serial precedence
    wf.addPrecedence(Workflow.Serial(A, B))
    wf.addPrecedence(Workflow.Serial(B, C))

    return wf


if __name__ == '__main__':
    wf = wf_erlang()

    # Convert to phase-type distribution
    alpha, T = wf.toPH()

    # Display results
    print('Erlang Workflow: A -> B -> C')
    print('Activity A: Erlang(mean=2.0, phases=4), SCV=0.25')
    print('Activity B: Erlang(mean=3.0, phases=2), SCV=0.50')
    print('Activity C: Exp(mean=1.0), SCV=1.00')
    print(f'Expected total mean: {2.0 + 3.0 + 1.0:.2f}')
    print(f'Computed PH mean: {wf.getMean():.4f}')
    print(f'Computed PH SCV: {wf.getSCV():.4f}')
    print(f'Number of phases: {T.shape[0]}')
