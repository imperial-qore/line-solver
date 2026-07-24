%% Example: Reward Aggregation Operations
% This example demonstrates the aggregation capabilities of the RewardState API.
%
% The RewardStateView class supports aggregation operations:
%   - total()  : Sum all values
%   - max()    : Maximum value
%   - min()    : Minimum value
%   - count()  : Count non-zero entries
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

clear; clc;

%% Model Definition
% Create a single queue with multiple job classes

model = Network('RewardAggregationExample');

% Nodes
source = Source(model, 'Source');
queue = Queue(model, 'Queue', SchedStrategy.FCFS);
sink = Sink(model, 'Sink');

% Queue parameters
queue.setNumberOfServers(1);
queue.setCapacity(4);  % Limit state space

% Job classes
class1 = OpenClass(model, 'HighPriority');
class2 = OpenClass(model, 'LowPriority');

% Arrivals
source.setArrival(class1, Exp(1.0));    % High priority arrival rate = 1.0
source.setArrival(class2, Exp(0.8));    % Low priority arrival rate = 0.8

% Service times
queue.setService(class1, Exp(3.0));
queue.setService(class2, Exp(3.0));

% Routing (same for both classes)
P = model.initRoutingMatrix;
P{class1} = Network.serialRouting(source, queue, sink);
P{class2} = Network.serialRouting(source, queue, sink);
model.link(P);

%% Define Rewards Using Aggregation
fprintf('Defining reward metrics using aggregation operations:\n\n');

% Total jobs at Queue (all classes)
model.setReward('TotalJobs', @(state) state.at(queue).total());
fprintf('  TotalJobs = state.at(queue).total()\n');

% Maximum population of any class at Queue
model.setReward('MaxClass', @(state) state.at(queue).max());
fprintf('  MaxClass = state.at(queue).max()\n');

% Count of classes with jobs at Queue
model.setReward('ClassCount', @(state) state.at(queue).count());
fprintf('  ClassCount = state.at(queue).count()\n');

% HighPriority jobs
model.setReward('HP_Jobs', @(state) state.at(queue, class1));
fprintf('  HP_Jobs = state.at(queue, class1)\n');

% LowPriority jobs
model.setReward('LP_Jobs', @(state) state.at(queue, class2));
fprintf('  LP_Jobs = state.at(queue, class2)\n');

% Weighted load (high priority weighted 2x)
model.setReward('WeightedLoad', @(state) ...
    2.0 * state.at(queue, class1) + 1.0 * state.at(queue, class2));
fprintf('  WeightedLoad = 2.0 * HP + 1.0 * LP\n');

%% Solve with CTMC Solver
fprintf('\nSolving with CTMC solver...\n\n');

options = Solver.defaultOptions;
options.verbose = 0;

solver = CTMC(model, options);

%% Get Steady-State Expected Rewards
[R, names] = solver.getAvgReward();

fprintf('=== Steady-State Expected Rewards (Aggregation) ===\n');
fprintf('  %-20s: %10.6f\n', 'TotalJobs', R(1));
fprintf('  %-20s: %10.6f\n', 'MaxClass', R(2));
fprintf('  %-20s: %10.6f\n', 'ClassCount', R(3));
fprintf('  %-20s: %10.6f\n', 'HP_Jobs', R(4));
fprintf('  %-20s: %10.6f\n', 'LP_Jobs', R(5));
fprintf('  %-20s: %10.6f\n', 'WeightedLoad', R(6));

fprintf('\nExample completed successfully.\n');
