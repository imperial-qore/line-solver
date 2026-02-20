%% Example BAS Network (Small)
% Smaller BAS network for faster testing
% 3 queues, 5 jobs

clear;
clc;

%% Define parameters
params = struct();

% Network size
params.M = 3;   % Number of queues
params.N = 5;   % Total population
params.f = 1;   % Finite capacity queue (1-based index)

% Queue capacities
params.F = [3; 5; 5];

% Number of phases per queue
params.K = [2; 2; 2];

% Routing probabilities
params.r = [
    0.0, 0.5, 0.5;   % from queue 1
    0.5, 0.0, 0.5;   % from queue 2
    0.5, 0.5, 0.0    % from queue 3
];

% Service rates mu{i}(k,h) - completion transition rates
params.mu = cell(3, 1);
for i = 1:3
    params.mu{i} = [
        1.0, 0.1;
        0.1, 0.5
    ];
end

% Background transition rates v{i}(k,h)
params.v = cell(3, 1);
for i = 1:3
    params.v{i} = zeros(2, 2);
end

%% Blocking configurations
% MR = number of blocking configurations
params.MR = 3;

% BB(m, i) = 1 if queue i is blocked in configuration m
params.BB = [
    0, 0, 0;   % m=1: no blocking
    0, 1, 0;   % m=2: queue 2 blocked
    0, 0, 1    % m=3: queue 3 blocked
];

% MM(m, order) = queue index at position 'order' in blocking list
params.MM = [
    0, 0;   % m=1
    2, 0;   % m=2: queue 2 first
    3, 0    % m=3: queue 3 first
];

% ZZ(m) = number of blocked queues in configuration m
params.ZZ = [0; 1; 1];

% ZM = maximum blocking depth
params.ZM = 1;

% MM1(m, j) = extended blocking order info
params.MM1 = [
    0, 0, 0;    % m=1
    0, 0, 0;    % m=2
    0, 0, 0     % m=3
];

%% Solve for minimum utilization at queue 1
fprintf('====================================\n');
fprintf('QRF BAS Example (Small)\n');
fprintf('====================================\n');
fprintf('Network: %d queues, %d jobs\n', params.M, params.N);
fprintf('Finite capacity queue: %d (capacity=%d)\n', params.f, params.F(params.f));
fprintf('Capacities: [%s]\n', num2str(params.F'));
fprintf('Phases: [%s]\n', num2str(params.K'));
fprintf('Blocking configurations: %d\n', params.MR);
fprintf('\n');

fprintf('--- Minimizing Utilization at Queue 1 ---\n');
[result_min, x_min, fval_min, exitflag_min] = qrf_bas(params, 1, 'min');

fprintf('\n--- Maximizing Utilization at Queue 1 ---\n');
[result_max, x_max, fval_max, exitflag_max] = qrf_bas(params, 1, 'max');

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
