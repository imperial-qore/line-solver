%% Parallel Workflow Example (AND-fork/join)
% A workflow with parallel activities that synchronize.
%
% Workflow:
%       +-> B --+
%   A --|       |--> D
%       +-> C --+
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

clear;
lineStart;

%% Define the workflow
wf = Workflow('ParallelWorkflow');

% Add activities
A = wf.addActivity('A', Exp.fitMean(1.0));
B = wf.addActivity('B', Exp.fitMean(2.0));
C = wf.addActivity('C', Exp.fitMean(3.0));
D = wf.addActivity('D', Exp.fitMean(0.5));

% Define precedences
wf.addPrecedence(Workflow.AndFork(A, {B, C}));
wf.addPrecedence(Workflow.AndJoin({B, C}, D));

%% Convert to phase-type distribution
ph = wf.toPH();

%% Display results
fprintf('Parallel Workflow: A -> [B || C] -> D\n');
fprintf('Activity means: A=1.0, B=2.0, C=3.0, D=0.5\n');

% Expected max of two exponentials: E[max(X,Y)] = 1/λ1 + 1/λ2 - 1/(λ1+λ2)
lambda_B = 0.5;  % rate for mean 2.0
lambda_C = 1/3;  % rate for mean 3.0
expected_max = 2.0 + 3.0 - 1/(lambda_B + lambda_C);
fprintf('Expected max(B,C) mean: %.4f\n', expected_max);
fprintf('Expected total mean: %.4f\n', 1.0 + expected_max + 0.5);
fprintf('Computed PH mean: %.4f\n', ph.getMean());
fprintf('Number of phases: %d\n', size(ph.getSubgenerator(), 1));
