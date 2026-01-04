%% Example: Reward Aggregation Operations
% This example demonstrates the aggregation capabilities of the new RewardState API.
%
% The RewardStateView class supports aggregation operations:
%   - total()  : Sum all values
%   - max()    : Maximum value
%   - min()    : Minimum value
%   - count()  : Count non-zero entries
%
% These can be applied to:
%   - All classes at a station: state.at(station).total()
%   - Single class across stations: state.forClass(class).total()
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

clear; clc;

%% Model Definition
% Create a small tandem queueing network with multiple job classes

model = Network('RewardAggregationExample');

% Nodes
source = Source(model, 'Source');
queue1 = Queue(model, 'Queue1', SchedStrategy.FCFS);
queue2 = Queue(model, 'Queue2', SchedStrategy.FCFS);
sink = Sink(model, 'Sink');

% Queue parameters
queue1.setNumberOfServers(1);
queue2.setNumberOfServers(1);

% Job classes
class1 = OpenClass(model, 'HighPriority');
class2 = OpenClass(model, 'LowPriority');

% Arrivals
source.setArrival(class1, Exp(1.0));    % High priority arrival rate = 1.0
source.setArrival(class2, Exp(0.8));    % Low priority arrival rate = 0.8

% Service times
queue1.setService(class1, Exp(1.5));
queue1.setService(class2, Exp(1.5));
queue2.setService(class1, Exp(2.0));
queue2.setService(class2, Exp(2.0));

% Routing: Series queue1 -> queue2 -> sink
model.link(Network.serialRouting(source, queue1, queue2, sink));

%% Define Rewards Using Aggregation
fprintf('Defining reward metrics using aggregation operations:\n\n');

% === Per-Station Aggregations ===

% Total jobs at Queue1 (all classes)
model.setReward('Q1_TotalJobs', @(state) state.at(queue1).total());
fprintf('  ✓ Q1_TotalJobs = state.at(queue1).total()\n');

% Total jobs at Queue2 (all classes)
model.setReward('Q2_TotalJobs', @(state) state.at(queue2).total());
fprintf('  ✓ Q2_TotalJobs = state.at(queue2).total()\n');

% Maximum population of any class at Queue1
model.setReward('Q1_MaxClass', @(state) state.at(queue1).max());
fprintf('  ✓ Q1_MaxClass = state.at(queue1).max()\n');

% Count of classes with jobs at Queue1
model.setReward('Q1_ClassCount', @(state) state.at(queue1).count());
fprintf('  ✓ Q1_ClassCount = state.at(queue1).count()\n');

% === Per-Class Aggregations ===

% Total HighPriority jobs across all stations
model.setReward('HP_TotalJobs', @(state) state.forClass(class1).total());
fprintf('  ✓ HP_TotalJobs = state.forClass(class1).total()\n');

% Total LowPriority jobs across all stations
model.setReward('LP_TotalJobs', @(state) state.forClass(class2).total());
fprintf('  ✓ LP_TotalJobs = state.forClass(class2).total()\n');

% === System-Wide Aggregations ===

% Total jobs in entire system
model.setReward('SystemJobs', @(state) ...
    state.at(queue1).total() + state.at(queue2).total());
fprintf('  ✓ SystemJobs = state.at(queue1).total() + state.at(queue2).total()\n');

% Weighted system load (high priority weighted 2x)
model.setReward('WeightedLoad', @(state) ...
    2.0 * state.forClass(class1).total() + ...
    1.0 * state.forClass(class2).total());
fprintf('  ✓ WeightedLoad = 2.0 * HP + 1.0 * LP\n');

% === Conditional Rewards ===

% Congestion indicator: 1 if Q1 has 5+ jobs, 0 otherwise
model.setReward('Q1_Congested', @(state) (state.at(queue1).total() >= 5));
fprintf('  ✓ Q1_Congested = (state.at(queue1).total() >= 5)\n');

%% Solve with CTMC Solver
fprintf('\nSolving with CTMC solver...\n\n');

options = Solver.defaultOptions;
options.verbose = 0;

solver = CTMC(model, options);

%% Get Steady-State Expected Rewards
[R, names] = solver.getAvgReward();

fprintf('=== Steady-State Expected Rewards (Aggregation) ===\n\n');

% Display grouped by type
fprintf('Per-Station Aggregations:\n');
fprintf('  %-25s: %10.6f\n', 'Q1_TotalJobs', R(1));
fprintf('  %-25s: %10.6f\n', 'Q2_TotalJobs', R(2));
fprintf('  %-25s: %10.6f\n', 'Q1_MaxClass', R(3));
fprintf('  %-25s: %10.6f\n', 'Q1_ClassCount', R(4));

fprintf('\nPer-Class Aggregations:\n');
fprintf('  %-25s: %10.6f\n', 'HP_TotalJobs', R(5));
fprintf('  %-25s: %10.6f\n', 'LP_TotalJobs', R(6));

fprintf('\nSystem-Wide Aggregations:\n');
fprintf('  %-25s: %10.6f\n', 'SystemJobs', R(7));
fprintf('  %-25s: %10.6f\n', 'WeightedLoad', R(8));

fprintf('\nConditional Rewards:\n');
fprintf('  %-25s: %10.6f (fraction of time congested)\n', ...
    'Q1_Congested', R(9));

%% Analysis
fprintf('\n=== Analysis ===\n');

% Expected length of stay in system
exp_system_jobs = R(7);
total_arrival_rate = 1.0 + 0.8;  % lambda1 + lambda2
exp_response_time = exp_system_jobs / total_arrival_rate;

fprintf('Expected jobs in system: %.6f\n', exp_system_jobs);
fprintf('Total arrival rate: %.6f\n', total_arrival_rate);
fprintf('Expected system response time (Little''s Law): %.6f\n', exp_response_time);

fprintf('\nClass composition:\n');
fprintf('  HighPriority jobs in system: %.6f (%.1f%%)\n', ...
    R(5), 100*R(5)/exp_system_jobs);
fprintf('  LowPriority jobs in system: %.6f (%.1f%%)\n', ...
    R(6), 100*R(6)/exp_system_jobs);

fprintf('\n✓ Example completed successfully.\n');
