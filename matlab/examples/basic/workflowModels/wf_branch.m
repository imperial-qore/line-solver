%% Probabilistic Branching Workflow Example (OR-fork/join)
% A workflow with probabilistic choice between branches.
%
% Workflow:
%       +-> B (60%) --+
%   A --|             |--> D
%       +-> C (40%) --+
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

clear;
lineStart;

%% Define the workflow
wf = Workflow('BranchingWorkflow');

% Add activities
A = wf.addActivity('A', Exp.fitMean(1.0));
B = wf.addActivity('B', Exp.fitMean(2.0));
C = wf.addActivity('C', Exp.fitMean(5.0));
D = wf.addActivity('D', Exp.fitMean(0.5));

% Define precedences with probabilities
wf.addPrecedence(Workflow.OrFork(A, {B, C}, [0.6, 0.4]));
wf.addPrecedence(Workflow.OrJoin({B, C}, D));

%% Convert to phase-type distribution
ph = wf.toPH();

%% Display results
fprintf('Branching Workflow: A -> [B(60%%) | C(40%%)] -> D\n');
fprintf('Activity means: A=1.0, B=2.0, C=5.0, D=0.5\n');
expected_branch = 0.6 * 2.0 + 0.4 * 5.0;
fprintf('Expected branch mean: 0.6*2.0 + 0.4*5.0 = %.2f\n', expected_branch);
fprintf('Expected total mean: %.2f\n', 1.0 + expected_branch + 0.5);
fprintf('Computed PH mean: %.4f\n', ph.getMean());
fprintf('Number of phases: %d\n', size(ph.getSubgenerator(), 1));
