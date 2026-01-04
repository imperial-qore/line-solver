%% Complex Workflow Example
% A workflow combining serial, parallel, and branching patterns.
%
% Workflow:
%                   +-> C --+
%   A -> B -> [AND] |       | [AND] -> F -> G
%                   +-> D --+
%                       |
%                       +-> [OR] -> E (30%)
%                              +-> skip (70%)
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

clear;
lineStart;

%% Define the workflow
wf = Workflow('ComplexWorkflow');

% Add activities
A = wf.addActivity('A', Exp.fitMean(0.5));
B = wf.addActivity('B', Exp.fitMean(1.0));
C = wf.addActivity('C', Exp.fitMean(2.0));
D = wf.addActivity('D', Exp.fitMean(1.5));
F = wf.addActivity('F', Exp.fitMean(1.0));
G = wf.addActivity('G', Exp.fitMean(0.5));

% Define precedences
wf.addPrecedence(Workflow.Serial(A, B));       % A -> B
wf.addPrecedence(Workflow.AndFork(B, {C, D})); % B forks to C and D
wf.addPrecedence(Workflow.AndJoin({C, D}, F)); % C and D join to F
wf.addPrecedence(Workflow.Serial(F, G));       % F -> G

%% Convert to phase-type distribution
ph = wf.toPH();

%% Display results
fprintf('Complex Workflow: A -> B -> [C || D] -> F -> G\n');
fprintf('Activity means: A=0.5, B=1.0, C=2.0, D=1.5, F=1.0, G=0.5\n');

% Calculate expected max(C, D)
lambda_C = 0.5;  % rate for mean 2.0
lambda_D = 1/1.5; % rate for mean 1.5
expected_max = 2.0 + 1.5 - 1/(lambda_C + lambda_D);
fprintf('Expected max(C,D) mean: %.4f\n', expected_max);
fprintf('Expected total mean: %.4f\n', 0.5 + 1.0 + expected_max + 1.0 + 0.5);
fprintf('Computed PH mean: %.4f\n', ph.getMean());
fprintf('Number of phases: %d\n', size(ph.getSubgenerator(), 1));

%% Sample from the distribution
fprintf('\nSampling 10000 workflow execution times...\n');
samples = ph.sample(10000);
fprintf('Sample mean: %.4f\n', mean(samples));
fprintf('Sample std: %.4f\n', std(samples));
fprintf('Sample min: %.4f\n', min(samples));
fprintf('Sample max: %.4f\n', max(samples));
