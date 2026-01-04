%% Example RS-RD Network
% MATLAB port of example_rsrd.mod
% RS-RD network with 20 jobs and 5 queues

clear;
clc;

%% Define parameters
params = struct();

% Network size
params.M = 5;   % Number of queues
params.N = 20;  % Total population

% Queue capacities
params.F = [5; 5; 5; 5; 5];

% Number of phases per queue
params.K = [2; 2; 2; 2; 2];

% Routing probabilities
params.r = [
    0.0, 0.5, 0.0, 0.0, 0.5;   % from queue 1
    0.5, 0.0, 0.5, 0.0, 0.0;   % from queue 2
    0.0, 0.5, 0.0, 0.5, 0.0;   % from queue 3
    0.0, 0.0, 0.5, 0.0, 0.5;   % from queue 4
    0.5, 0.0, 0.0, 0.5, 0.0    % from queue 5
];

% Service rates mu{i}(k,h) - completion transition rates
% Each mu{i} is a K(i) x K(i) matrix
params.mu = cell(5, 1);
for i = 1:5
    params.mu{i} = [
        1.016186e+00, 2.585708e-05;
        1.569888e-03, 1.413298e-02
    ];
end

% Background transition rates v{i}(k,h)
% Each v{i} is a K(i) x K(i) matrix
params.v = cell(5, 1);
for i = 1:5
    params.v{i} = zeros(2, 2);
end

% Load-dependent rates alpha{i}(n) - default is 1 for all
% (not specified in example, so using default)

%% Solve for minimum utilization at queue 1
fprintf('====================================\n');
fprintf('QRF RS-RD Example\n');
fprintf('====================================\n');
fprintf('Network: %d queues, %d jobs\n', params.M, params.N);
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
fprintf('\nAll queue utilizations (at minimum):\n');
for i = 1:params.M
    fprintf('  Queue %d: U = %.6f\n', i, result_min.U(i));
end
