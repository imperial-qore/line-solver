%% Workflow with APH Distributions
% A workflow using APH (Acyclic Phase-Type) distributions
% for flexible variability modeling.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

clear;
lineStart;

%% Define the workflow
wf = Workflow('APHWorkflow');

% Add activities with different SCVs using APH
A = wf.addActivity('A', APH.fitMeanAndSCV(1.0, 0.5));   % Low variability
B = wf.addActivity('B', APH.fitMeanAndSCV(2.0, 2.0));   % High variability
C = wf.addActivity('C', APH.fitMeanAndSCV(1.5, 1.0));   % Standard variability

% Define serial precedence
wf.addPrecedence(Workflow.Serial(A, B, C));

%% Convert to phase-type distribution
ph = wf.toPH();

%% Display results
fprintf('APH Workflow: A -> B -> C\n');
fprintf('Activity A: APH(mean=1.0, SCV=0.5)\n');
fprintf('Activity B: APH(mean=2.0, SCV=2.0)\n');
fprintf('Activity C: APH(mean=1.5, SCV=1.0)\n');
fprintf('Expected total mean: %.2f\n', 1.0 + 2.0 + 1.5);
fprintf('Computed PH mean: %.4f\n', ph.getMean());
fprintf('Computed PH SCV: %.4f\n', ph.getSCV());
fprintf('Number of phases: %d\n', size(ph.getSubgenerator(), 1));
