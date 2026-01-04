%% Example: Diffusion Method for Fluid Approximation
%
% This example demonstrates the diffusion approximation method for solving
% closed queueing networks. The diffusion method uses Euler-Maruyama SDE
% approximation to model stochastic behavior in fluid systems.
%
% The diffusion method is particularly useful for:
% - Closed queueing networks with moderate to large populations
% - Systems operating under significant load
% - Rapid approximate analysis when detailed simulation is not needed
%
% Supported scheduling disciplines: PS, FCFS, INF, SIRO
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

clear; clc;

%% Create a simple closed queueing network
% This is a two-station network with N jobs that circulate continuously.
% Station 1: Processing queue (FCFS)
% Station 2: Think/Delay node (infinite servers)

fprintf('Creating closed queueing network for diffusion analysis...\n\n');

% Model parameters
N = 10;              % Number of circulating jobs
mu_service = 2.0;    % Service rate at station 1 (mean = 0.5)
mu_delay = 1.0;      % Delay rate at station 2 (mean = 1.0)

% Create network
model = Network('DiffusionNetwork');

% Create nodes
queue = Queue(model, 'ProcessingQueue', SchedStrategy.FCFS);
delay = Delay(model, 'ThinkTime');

% Create job class
jobClass = ClosedClass(model, 'Jobs', N, delay, 0);

% Set service times
queue.setService(jobClass, Exp(1/mu_service));
delay.setService(jobClass, Exp(1/mu_delay));

% Set routing: queue -> delay -> queue
P = model.initRoutingMatrix();
P{jobClass}(queue, delay) = 1.0;
P{jobClass}(delay, queue) = 1.0;
model.link(P);

fprintf('Network created:\n');
fprintf('  - 1 Processing Queue (FCFS)\n');
fprintf('  - 1 Delay node (think time)\n');
fprintf('  - N = %d circulating jobs\n', N);
fprintf('  - Service rate: %.2f (mean = %.2f)\n', mu_service, 1/mu_service);
fprintf('  - Delay rate: %.2f (mean = %.2f)\n\n', mu_delay, 1/mu_delay);

%% Solve using Fluid solver with default method
fprintf('=== Solving with Default Fluid Method ===\n');

try
    solver_default = SolverFluid(model);
    solver_default.runAnalyzer();
    results_default = solver_default.getAvgTable();

    fprintf('Default Fluid solution obtained.\n');
    fprintf('Results:\n');
    disp(results_default);

    % Extract key metrics
    stations = string(results_default.Station);
    queue_idx = find(stations == "ProcessingQueue", 1);

    if ~isempty(queue_idx)
        resp_time_default = results_default.RespT(queue_idx);
        tput_default = results_default.Tput(queue_idx);
        fprintf('\nDefault method - Queue response time: %.4f\n', resp_time_default);
        fprintf('Default method - Throughput: %.4f\n', tput_default);
    end

catch ME
    fprintf('Error with default method: %s\n', ME.message);
    resp_time_default = NaN;
    tput_default = NaN;
end

%% Solve using Fluid solver with diffusion method
fprintf('\n=== Solving with Diffusion Method ===\n');

try
    solver_diffusion = SolverFluid(model);
    % Set the diffusion method
    solver_diffusion.solveropt.method('diffusion');
    solver_diffusion.runAnalyzer();
    results_diffusion = solver_diffusion.getAvgTable();

    fprintf('Diffusion method solution obtained.\n');
    fprintf('Results:\n');
    disp(results_diffusion);

    % Extract key metrics
    stations = string(results_diffusion.Station);
    queue_idx = find(stations == "ProcessingQueue", 1);

    if ~isempty(queue_idx)
        resp_time_diffusion = results_diffusion.RespT(queue_idx);
        tput_diffusion = results_diffusion.Tput(queue_idx);
        fprintf('\nDiffusion method - Queue response time: %.4f\n', resp_time_diffusion);
        fprintf('Diffusion method - Throughput: %.4f\n', tput_diffusion);
    end

catch ME
    fprintf('Error with diffusion method: %s\n', ME.message);
    resp_time_diffusion = NaN;
    tput_diffusion = NaN;
end

%% Compare the two methods
if ~isnan(resp_time_default) && ~isnan(resp_time_diffusion)
    fprintf('\n=== Comparison ===\n');
    fprintf('%-30s %15s %15s\n', 'Metric', 'Default', 'Diffusion');
    fprintf('%s\n', repmat('-', 1, 60));
    fprintf('%-30s %15.6f %15.6f\n', 'Response Time', resp_time_default, resp_time_diffusion);
    fprintf('%-30s %15.6f %15.6f\n', 'Throughput', tput_default, tput_diffusion);

    % Calculate relative difference
    if resp_time_default > 0
        rt_diff = abs(resp_time_diffusion - resp_time_default) / resp_time_default * 100;
        fprintf('%-30s %15s %15.2f%%\n', 'RespT relative diff', '', rt_diff);
    end

    if tput_default > 0
        tp_diff = abs(tput_diffusion - tput_default) / tput_default * 100;
        fprintf('%-30s %15s %15.2f%%\n', 'Tput relative diff', '', tp_diff);
    end
end

%% Example 2: Larger network with multiple service rates
fprintf('\n\n========================================\n');
fprintf('EXAMPLE 2: Multi-Station Network\n');
fprintf('========================================\n\n');

fprintf('Creating a 3-station network for diffusion analysis...\n\n');

% Create network with 3 stations
model2 = Network('MultiStationDiffusion');

% Create nodes
queue1 = Queue(model2, 'Station1', SchedStrategy.FCFS);
queue2 = Queue(model2, 'Station2', SchedStrategy.FCFS);
delay = Delay(model2, 'Delay');

% Create job class
N2 = 15;
jobClass2 = ClosedClass(model2, 'Jobs', N2, delay, 0);

% Set service times (different at each station)
queue1.setService(jobClass2, Exp(0.4));   % Service time = 0.4
queue2.setService(jobClass2, Exp(0.6));   % Service time = 0.6
delay.setService(jobClass2, Exp(0.5));    % Think time = 0.5

% Set routing: 1 -> 2 -> Delay -> 1
P2 = model2.initRoutingMatrix();
P2{jobClass2}(queue1, queue2) = 1.0;
P2{jobClass2}(queue2, delay) = 1.0;
P2{jobClass2}(delay, queue1) = 1.0;
model2.link(P2);

fprintf('Network created:\n');
fprintf('  - 3 stations (2 queues + 1 delay)\n');
fprintf('  - N = %d jobs\n', N2);
fprintf('  - Station 1: FCFS, service time = 0.4\n');
fprintf('  - Station 2: FCFS, service time = 0.6\n');
fprintf('  - Delay: Think time = 0.5\n\n');

%% Solve with diffusion method
fprintf('Solving 3-station network with diffusion method...\n');

try
    solver2 = SolverFluid(model2);
    solver2.solveropt.method('diffusion');
    solver2.runAnalyzer();
    results2 = solver2.getAvgTable();

    fprintf('Diffusion solution for 3-station network:\n');
    disp(results2);

    % Extract utilization at each queue
    stations2 = string(results2.Station);
    st1_idx = find(stations2 == "Station1", 1);
    st2_idx = find(stations2 == "Station2", 1);

    fprintf('\nStation utilizations:\n');
    if ~isempty(st1_idx)
        fprintf('  Station1 utilization: %.4f\n', results2.Util(st1_idx));
    end
    if ~isempty(st2_idx)
        fprintf('  Station2 utilization: %.4f\n', results2.Util(st2_idx));
    end

catch ME
    fprintf('Error: %s\n', ME.message);
end

%% Summary
fprintf('\n========================================\n');
fprintf('DIFFUSION METHOD SUMMARY\n');
fprintf('========================================\n\n');

fprintf('The diffusion method is a fluid-based approximation that:\n');
fprintf('  1. Uses Euler-Maruyama SDE for stochastic modeling\n');
fprintf('  2. Works for closed networks with moderate to large N\n');
fprintf('  3. Supports PS, FCFS, INF, SIRO scheduling\n');
fprintf('  4. Provides rapid analysis compared to simulation\n');
fprintf('  5. Is particularly accurate for heavily loaded systems\n\n');

fprintf('When to use diffusion:\n');
fprintf('  - When you need faster analysis than simulation\n');
fprintf('  - For closed networks with N >= 10 jobs\n');
fprintf('  - For heavily loaded systems (high utilization)\n');
fprintf('  - When analytical methods fail or are slow\n\n');

fprintf('Comparison with other methods:\n');
fprintf('  - Default Fluid: ODE-based, good for general use\n');
fprintf('  - Diffusion: SDE-based, better for stochastic effects\n');
fprintf('  - Closing: For open networks\n');
fprintf('  - StateDependent: For systems with state-dependent arrivals\n');
