%% Serial Workflow Example
% A simple workflow with three sequential activities.
%
% Workflow: A -> B -> C
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

clear;
lineStart;

%% Define the workflow
wf = Workflow('SerialWorkflow');

% Add activities with exponential service times
A = wf.addActivity('A', Exp.fitMean(1.0));
B = wf.addActivity('B', Exp.fitMean(2.0));
C = wf.addActivity('C', Exp.fitMean(1.5));

% Define serial precedence
wf.addPrecedence(Workflow.Serial(A, B, C));

%% Convert to phase-type distribution
ph = wf.toPH();

%% Display results
fprintf('Serial Workflow: A -> B -> C\n');
fprintf('Activity means: A=1.0, B=2.0, C=1.5\n');
fprintf('Expected total mean: %.2f\n', 1.0 + 2.0 + 1.5);
fprintf('Computed PH mean: %.4f\n', ph.getMean());
fprintf('Computed PH SCV: %.4f\n', ph.getSCV());
fprintf('Number of phases: %d\n', size(ph.getSubgenerator(), 1));
