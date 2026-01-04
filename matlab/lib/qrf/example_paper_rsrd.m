%% Example from Paper - Bounding Problem 1 (RS-RD)
% M=3 queues, F1=4, N=10, p=0.90 (balanced demands)
% This is the validation example from Section 7.1 of the TOMACS paper
% All queues use K=2 phases to match the AMPL model structure

clear;
clc;

%% Define parameters
params = struct();

% Network size
params.M = 3;   % Number of queues
params.N = 10;  % Total population (N >= F1+1 for blocking)

% Queue capacities (F1=4 finite, others infinite=N)
params.F = [4; 10; 10];

% Number of phases per queue (K=2 for all to maintain consistency)
params.K = [2; 2; 2];

% Routing probabilities with p=0.90 (balanced demands)
p = 0.90;
params.r = [
    0.10, 0.50, 0.40;   % from queue 1
    p,    0.00, 1-p;    % from queue 2
    0.00, 0.50, 0.50    % from queue 3
];

% Service rates mu{i}(k,h) - completion transition rates
% Queue 1: MAP with D1 matrix (from paper Equation 7)
% D1 = [1.016186165025678  0.000025857082896]
%      [0.001569887597955  0.014132983910493]
params.mu = cell(3, 1);
params.mu{1} = [
    1.016186165025678, 0.000025857082896;
    0.001569887597955, 0.014132983910493
];
% Queues 2,3: exponential with rate 1 (2-phase representation)
params.mu{2} = [
    1.0, 0.0;
    0.0, 1.0
];
params.mu{3} = [
    1.0, 0.0;
    0.0, 1.0
];

% Background transition rates v{i}(k,h)
% For MAP: v comes from D0 off-diagonal elements (=0 in this case)
params.v = cell(3, 1);
params.v{1} = zeros(2, 2);
params.v{2} = zeros(2, 2);
params.v{3} = zeros(2, 2);

%% Solve for utilization bounds at queue 1
fprintf('====================================\n');
fprintf('Paper Example: RS-RD Bounding Problem 1\n');
fprintf('====================================\n');
fprintf('Network: %d queues, %d jobs\n', params.M, params.N);
fprintf('Finite capacity queue: 1 (capacity=%d)\n', params.F(1));
fprintf('Routing parameter p = %.2f (balanced demands)\n', p);
fprintf('Capacities: [%s]\n', num2str(params.F'));
fprintf('Phases: [%s]\n', num2str(params.K'));
fprintf('\n');

fprintf('--- Minimizing Utilization at Queue 1 ---\n');
[result_min, x_min, fval_min, exitflag_min] = qrf_rsrd(params, 1, 'min');

fprintf('\n--- Maximizing Utilization at Queue 1 ---\n');
[result_max, x_max, fval_max, exitflag_max] = qrf_rsrd(params, 1, 'max');

%% Summary
fprintf('\n====================================\n');
fprintf('SUMMARY\n');
fprintf('====================================\n');
fprintf('Utilization bounds for Queue 1:\n');
fprintf('  Lower bound: %.6f\n', fval_min);
fprintf('  Upper bound: %.6f\n', fval_max);
fprintf('\n');
fprintf('Reference from AMPL/GLPK:\n');
fprintf('  Lower bound: 0.733487 (expected)\n');
fprintf('\nAll queue utilizations (at minimum):\n');
for i = 1:params.M
    fprintf('  Queue %d: U = %.6f\n', i, result_min.U(i));
end
