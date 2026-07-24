%% Workflow with Erlang Distribution
% A workflow using Erlang distributions for lower variability.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

clear;
lineStart;

%% Define the workflow
wf = Workflow('ErlangWorkflow');

% Add activities with Erlang distributions (SCV < 1)
A = wf.addActivity('A', Erlang.fitMeanAndOrder(2.0, 4));  % 4 phases, SCV=0.25
B = wf.addActivity('B', Erlang.fitMeanAndOrder(3.0, 2));  % 2 phases, SCV=0.5
C = wf.addActivity('C', Exp.fitMean(1.0));                % 1 phase, SCV=1.0

% Define serial precedence
wf.addPrecedence(Workflow.Serial(A, B, C));

%% Convert to phase-type distribution
ph = wf.toPH();

%% Display results
fprintf('Erlang Workflow: A -> B -> C\n');
fprintf('Activity A: Erlang(mean=2.0, phases=4), SCV=0.25\n');
fprintf('Activity B: Erlang(mean=3.0, phases=2), SCV=0.50\n');
fprintf('Activity C: Exp(mean=1.0), SCV=1.00\n');
fprintf('Expected total mean: %.2f\n', 2.0 + 3.0 + 1.0);
fprintf('Computed PH mean: %.4f\n', ph.getMean());
fprintf('Computed PH SCV: %.4f\n', ph.getSCV());
fprintf('Number of phases: %d\n', size(ph.getSubgenerator(), 1));
