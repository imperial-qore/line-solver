%% Example: Multi-Class Reward Analysis
% This example demonstrates reward-based CTMC analysis in a multi-class queueing
% network, showing how to define per-class metrics and analyze class-specific behavior.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

clear; clc;

%% Model Definition
% Create an M/M/2 queue with two job classes (prioritized service)

model = Network('MultiClassRewardExample');

% Nodes
source = Source(model, 'Source');
queue = Queue(model, 'Queue', SchedStrategy.PS);  % Processor Sharing
sink = Sink(model, 'Sink');

% Server capacity and buffer limit
queue.setNumberOfServers(2);
queue.setCapacity(6);  % Limit state space

% Job classes with different characteristics
class_interactive = OpenClass(model, 'Interactive');
class_batch = OpenClass(model, 'Batch');

% Different arrival rates
source.setArrival(class_interactive, Exp(2.0));  % Interactive: λ = 2.0
source.setArrival(class_batch, Exp(1.5));        % Batch: λ = 1.5

% Different service time requirements
queue.setService(class_interactive, Exp(0.5));   % Interactive: μ = 2.0 (fast)
queue.setService(class_batch, Exp(1.0));         % Batch: μ = 1.0 (slow)

% Topology (same for both classes)
P = model.initRoutingMatrix;
P{class_interactive} = Network.serialRouting(source, queue, sink);
P{class_batch} = Network.serialRouting(source, queue, sink);
model.link(P);

%% Define Per-Class Rewards
fprintf('Defining per-class reward metrics:\n\n');

% === Interactive Class Metrics ===

% Jobs of Interactive class
model.setReward('Interactive_QLen', @(state) state.at(queue, class_interactive));
fprintf('  ✓ Interactive_QLen = state.at(queue, Interactive)\n');

% Utilization contributed by Interactive class
model.setReward('Interactive_Util', Reward.utilization(queue, class_interactive));
fprintf('  ✓ Interactive_Util = Reward.utilization(queue, Interactive)\n');

% === Batch Class Metrics ===

% Jobs of Batch class
model.setReward('Batch_QLen', @(state) state.at(queue, class_batch));
fprintf('  ✓ Batch_QLen = state.at(queue, Batch)\n');

% Utilization contributed by Batch class
model.setReward('Batch_Util', Reward.utilization(queue, class_batch));
fprintf('  ✓ Batch_Util = Reward.utilization(queue, Batch)\n');

% === Comparative Metrics ===

% Total queue length (both classes)
model.setReward('Total_QLen', @(state) state.at(queue).total());
fprintf('  ✓ Total_QLen = state.at(queue).total()\n');

% Total system utilization
model.setReward('Total_Util', Reward.utilization(queue));
fprintf('  ✓ Total_Util = Reward.utilization(queue)\n');

% Ratio of Interactive to Batch jobs
model.setReward('Interactive_Ratio', @(state) ...
    state.at(queue, class_interactive) / ...
    max(state.at(queue, class_batch), 0.001));  % Avoid division by zero
fprintf('  ✓ Interactive_Ratio = Interactive / Batch\n');

% === Priority-Aware Metrics ===

% Weighted response time cost (Interactive weighted 3x, Batch weighted 1x)
% Approximates weighted system delay
model.setReward('Weighted_Cost', @(state) ...
    3.0 * state.at(queue, class_interactive) + ...
    1.0 * state.at(queue, class_batch));
fprintf('  ✓ Weighted_Cost = 3.0*Interactive + 1.0*Batch\n');

% Service fairness indicator: 1 if both classes present, 0 otherwise
model.setReward('Fairness', @(state) ...
    double(state.at(queue, class_interactive) > 0 && ...
           state.at(queue, class_batch) > 0));
fprintf('  ✓ Fairness = (Interactive > 0) AND (Batch > 0)\n');

%% Solve with CTMC Solver
fprintf('\nSolving with CTMC solver...\n\n');

options = Solver.defaultOptions;
options.verbose = 0;

solver = CTMC(model, options);

%% Get Steady-State Expected Rewards
[R, names] = solver.getAvgReward();

fprintf('=== Steady-State Expected Rewards (Multi-Class) ===\n\n');

fprintf('Interactive Class:\n');
fprintf('  %-25s: %10.6f (jobs)\n', 'Interactive_QLen', R(1));
fprintf('  %-25s: %10.6f (server capacity)\n', 'Interactive_Util', R(2));

fprintf('\nBatch Class:\n');
fprintf('  %-25s: %10.6f (jobs)\n', 'Batch_QLen', R(3));
fprintf('  %-25s: %10.6f (server capacity)\n', 'Batch_Util', R(4));

fprintf('\nComparative Metrics:\n');
fprintf('  %-25s: %10.6f (both classes)\n', 'Total_QLen', R(5));
fprintf('  %-25s: %10.6f (2 servers)\n', 'Total_Util', R(6));
fprintf('  %-25s: %10.6f (ratio)\n', 'Interactive_Ratio', R(7));

fprintf('\nPriority-Aware Metrics:\n');
fprintf('  %-25s: %10.6f (weighted cost)\n', 'Weighted_Cost', R(8));
fprintf('  %-25s: %10.6f (fraction of time both present)\n', ...
    'Fairness', R(9));

%% Analytical Validation
fprintf('\n=== Analysis ===\n');

lambda_int = 2.0;
mu_int = 2.0;
lambda_batch = 1.5;
mu_batch = 1.0;
c = 2;  % Number of servers

rho_int = lambda_int / (c * mu_int);     % 2 / (2*2) = 0.5
rho_batch = lambda_batch / (c * mu_batch);  % 1.5 / (2*1) = 0.75
rho_total = rho_int + rho_batch;          % 1.25

fprintf('System Characteristics:\n');
fprintf('  Interactive: λ=%.1f, μ=%.1f per server\n', lambda_int, mu_int);
fprintf('  Batch:       λ=%.1f, μ=%.1f per server\n', lambda_batch, mu_batch);
fprintf('  Servers: %d\n', c);
fprintf('  Total utilization: %.3f (%.1f%% capacity)\n', ...
    rho_total, 100*rho_total);

% Class composition
total_qlen = R(5);
if total_qlen > 0
    pct_int = 100 * R(1) / total_qlen;
    pct_batch = 100 * R(3) / total_qlen;
    fprintf('\nQueue Composition:\n');
    fprintf('  Interactive jobs: %.2f%% of queue\n', pct_int);
    fprintf('  Batch jobs: %.2f%% of queue\n', pct_batch);
end

% Response time estimation (Little's Law)
total_arrival = lambda_int + lambda_batch;
est_resp_time_int = R(1) / lambda_int;
est_resp_time_batch = R(3) / lambda_batch;

fprintf('\nEstimated Response Times (Little''s Law):\n');
fprintf('  Interactive: %.6f time units\n', est_resp_time_int);
fprintf('  Batch: %.6f time units\n', est_resp_time_batch);

fprintf('\n✓ Example completed successfully.\n');
