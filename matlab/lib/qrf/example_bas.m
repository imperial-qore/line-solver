%% Example BAS Network
% MATLAB port of example_bas.mod
% BAS network with 10 jobs and 5 queues

clear;
clc;

%% Define parameters
params = struct();

% Network size
params.M = 5;   % Number of queues
params.N = 10;  % Total population
params.f = 1;   % Finite capacity queue (1-based index)

% Queue capacities
params.F = [5; 10; 10; 10; 10];

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
params.mu = cell(5, 1);
for i = 1:5
    params.mu{i} = [
        1.016186e+00, 2.585708e-05;
        1.569888e-03, 1.413298e-02
    ];
end

% Background transition rates v{i}(k,h)
params.v = cell(5, 1);
for i = 1:5
    params.v{i} = zeros(2, 2);
end

%% Blocking configurations
% MR = number of blocking configurations
params.MR = 5;

% BB(m, i) = 1 if queue i is blocked in configuration m
params.BB = [
    0, 0, 0, 0, 0;   % m=1: no blocking
    0, 1, 0, 0, 0;   % m=2: queue 2 blocked
    0, 1, 0, 0, 1;   % m=3: queues 2,5 blocked
    0, 0, 0, 0, 1;   % m=4: queue 5 blocked
    0, 1, 0, 0, 1    % m=5: queues 2,5 blocked
];

% MM(m, order) = queue index at position 'order' in blocking list
% MM(m,1) = first blocked queue, MM(m,2) = second blocked queue
params.MM = [
    0, 0;   % m=1
    2, 0;   % m=2: queue 2 first
    2, 5;   % m=3: queue 2 first, queue 5 second
    5, 0;   % m=4: queue 5 first
    5, 2    % m=5: queue 5 first, queue 2 second
];

% ZZ(m) = number of blocked queues in configuration m
params.ZZ = [0; 1; 2; 1; 2];

% ZM = maximum blocking depth
params.ZM = 2;

% MM1(m, j) = extended blocking order info
% Used for THM3L constraint
params.MM1 = [
    0, 0, 0, 0, 0;    % m=1
   -1,-1,-1,-1, 3;    % m=2
    0, 0, 0, 0, 0;    % m=3
   -1, 5,-1,-1,-1;    % m=4
    0, 0, 0, 0, 0     % m=5
];

%% Solve for minimum utilization at queue 1
fprintf('====================================\n');
fprintf('QRF BAS Example\n');
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
