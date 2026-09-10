%% Loop Workflow Example
% A workflow with a repeated activity.
%
% Workflow: A -> [B x 3] -> C
% Activity B is executed 3 times.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

clear;
lineStart;

%% Define the workflow
wf = Workflow('LoopWorkflow');

% Add activities
A = wf.addActivity('A', Exp.fitMean(1.0));
B = wf.addActivity('B', Exp.fitMean(2.0));
C = wf.addActivity('C', Exp.fitMean(0.5));

% Define loop: A runs once, B runs 3 times, C runs once
wf.addPrecedence(Workflow.Loop(A, {B, C}, 3));

%% Convert to phase-type distribution
ph = wf.toPH();

%% Display results
fprintf('Loop Workflow: A -> [B x 3] -> C\n');
fprintf('Activity means: A=1.0, B=2.0, C=0.5\n');
fprintf('Expected total mean: 1.0 + 3*2.0 + 0.5 = %.2f\n', 1.0 + 3*2.0 + 0.5);
fprintf('Computed PH mean: %.4f\n', ph.getMean());
fprintf('Number of phases: %d\n', size(ph.getSubgenerator(), 1));
