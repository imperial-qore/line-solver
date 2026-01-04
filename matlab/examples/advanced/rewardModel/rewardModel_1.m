%% Example: Using setReward for reward-based CTMC analysis
% This example demonstrates the setReward feature for defining custom
% reward functions on a queueing network model and computing steady-state
% expected rewards using the CTMC solver.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

clear; clc;

%% Model Definition
% Create a simple M/M/1/K queue with finite buffer

model = Network('RewardExample');

% Block 1: nodes
source = Source(model, 'Source');
queue = Queue(model, 'Queue', SchedStrategy.FCFS);
sink = Sink(model, 'Sink');

% Set finite buffer capacity
queue.setNumberOfServers(1);
queue.setCapacity(3);  % Maximum jobs in the system

% Block 2: job classes
oclass = OpenClass(model, 'Class1');
source.setArrival(oclass, Exp(2));  % Arrival rate = 2
queue.setService(oclass, Exp(3));   % Service rate = 3 (utilization ~ 0.67)

% Block 3: topology
model.link(Network.serialRouting(source,queue,sink));

%% Define Reward Functions
% setReward(name, rewardFn) where rewardFn uses RewardState for clear state access
% The new API provides intuitive methods: state.at(node, class), state.at(node).total(), etc.

% Reward 1: Queue length (number of jobs in the queue)
% Old way: state(2) - unclear what state(2) means
% New way: state.at(queue, oclass) - clear reference to queue and class
model.setReward('QueueLength', @(state) state.at(queue, oclass));

% Reward 2: Utilization (1 if server busy, 0 if idle)
% Using Reward template for clarity
model.setReward('Utilization', Reward.utilization(queue, oclass));

% Reward 3: Blocking indicator (1 if buffer full, 0 otherwise)
% Using Reward template for common metric
model.setReward('BlockingProb', Reward.blocking(queue));

% Reward 4: Weighted queue cost (quadratic penalty for long queues)
% Custom reward using state accessor
model.setReward('QueueCost', @(state) state.at(queue, oclass)^2);

%% Solve with CTMC Solver
fprintf('Solving with CTMC solver...\n\n');

options = Solver.defaultOptions;
options.verbose = 1;

solver = CTMC(model, options);

%% Get Steady-State Expected Rewards
[R, names] = solver.getAvgReward();

fprintf('\n=== Steady-State Expected Rewards ===\n');
for i = 1:length(names)
    fprintf('%15s: %.6f\n', names{i}, R(i));
end

%% Get Transient Reward Analysis
% Get full reward trajectories over time
[t, V, names, stateSpace] = solver.getReward();

fprintf('\n=== State Space ===\n');
fprintf('Number of states: %d\n', size(stateSpace, 1));
fprintf('State dimensions: %d\n', size(stateSpace, 2));

%% Plot Reward Convergence
figure('Name', 'Reward Convergence', 'Position', [100 100 800 600]);

nRewards = length(names);
for r = 1:nRewards
    subplot(2, 2, r);

    % Plot mean reward rate over time (averaged over all initial states)
    meanReward = mean(V{r}, 2);  % Average over states

    % Compute reward rate (reward per unit time)
    rewardRate = zeros(length(t), 1);
    for k = 2:length(t)
        rewardRate(k) = meanReward(k) / t(k);
    end

    plot(t(2:end), rewardRate(2:end), 'LineWidth', 1.5);
    hold on;
    yline(R(r), 'r--', 'LineWidth', 1.5);
    hold off;

    xlabel('Time');
    ylabel('Expected Reward Rate');
    title(sprintf('%s (Steady-state: %.4f)', names{r}, R(r)));
    legend('Transient', 'Steady-state', 'Location', 'best');
    grid on;
end

sgtitle('Convergence of Reward Functions to Steady-State');

%% Compare with Analytical Results for M/M/1/K
fprintf('\n=== Comparison with M/M/1/K Analytical Results ===\n');

lambda = 2;  % Arrival rate
mu = 3;      % Service rate
rho = lambda / mu;
K = 10;      % Buffer capacity

% For M/M/1/K queue, steady-state probabilities:
% pi(n) = (1-rho) * rho^n / (1 - rho^(K+1))
if rho ~= 1
    pi = zeros(K+1, 1);
    for n = 0:K
        pi(n+1) = (1 - rho) * rho^n / (1 - rho^(K+1));
    end
else
    pi = ones(K+1, 1) / (K+1);
end

% Analytical expected queue length
L_analytical = sum((0:K)' .* pi);

% Analytical utilization (P(server busy) = 1 - pi(0))
U_analytical = 1 - pi(1);

% Analytical blocking probability = pi(K)
B_analytical = pi(end);

% Analytical queue cost (E[N^2])
Cost_analytical = sum(((0:K).^2)' .* pi);

fprintf('%15s: LINE = %.6f, Analytical = %.6f, Error = %.2e\n', ...
    'QueueLength', R(1), L_analytical, abs(R(1) - L_analytical));
fprintf('%15s: LINE = %.6f, Analytical = %.6f, Error = %.2e\n', ...
    'Utilization', R(2), U_analytical, abs(R(2) - U_analytical));
fprintf('%15s: LINE = %.6f, Analytical = %.6f, Error = %.2e\n', ...
    'BlockingProb', R(3), B_analytical, abs(R(3) - B_analytical));
fprintf('%15s: LINE = %.6f, Analytical = %.6f, Error = %.2e\n', ...
    'QueueCost', R(4), Cost_analytical, abs(R(4) - Cost_analytical));

fprintf('\nExample completed successfully.\n');
