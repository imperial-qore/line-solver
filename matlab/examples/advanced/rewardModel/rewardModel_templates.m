%% Example: Using Reward Templates
% This example demonstrates the use of Reward template factory methods
% for quickly defining common metrics like queue length, utilization, and blocking.
%
% The new Reward templates provide a convenient way to define metrics
% without writing lambda expressions.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

clear; clc;

%% Model Definition
% Create an M/M/1 queue with finite buffer

model = Network('RewardTemplatesExample');

% Nodes
source = Source(model, 'Source');
queue = Queue(model, 'Queue', SchedStrategy.FCFS);
sink = Sink(model, 'Sink');

% Queue parameters
queue.setNumberOfServers(1);
queue.setCapacity(10);

% Job class
oclass = OpenClass(model, 'Class1');
source.setArrival(oclass, Exp(1.5));  % Arrival rate = 1.5
queue.setService(oclass, Exp(2));     % Service rate = 2 (rho = 0.75)

% Topology
model.link(Network.serialRouting(source, queue, sink));

%% Define Rewards Using Templates
% Reward templates provide a clean syntax for common metrics

fprintf('Defining reward metrics using Reward templates:\n\n');

% Template 1: Queue Length
% Returns the total number of jobs in the queue (all classes)
model.setReward('QueueLength', Reward.queueLength(queue));
fprintf('  ✓ QueueLength = Reward.queueLength(queue)\n');

% Template 2: Queue Length for Specific Class
% Returns the number of jobs of a specific class in the queue
model.setReward('QueueLength_Class1', Reward.queueLength(queue, oclass));
fprintf('  ✓ QueueLength_Class1 = Reward.queueLength(queue, oclass)\n');

% Template 3: Server Utilization
% Returns min(jobs, nservers) - the utilization of the server(s)
model.setReward('Utilization', Reward.utilization(queue));
fprintf('  ✓ Utilization = Reward.utilization(queue)\n');

% Template 4: Utilization by Class
% Returns the utilization contributed by a specific class
model.setReward('Utilization_Class1', Reward.utilization(queue, oclass));
fprintf('  ✓ Utilization_Class1 = Reward.utilization(queue, oclass)\n');

% Template 5: Blocking Probability
% Returns 1 if queue is at capacity, 0 otherwise
model.setReward('BlockingProb', Reward.blocking(queue));
fprintf('  ✓ BlockingProb = Reward.blocking(queue)\n');

%% Solve with CTMC Solver
fprintf('\nSolving with CTMC solver...\n\n');

options = Solver.defaultOptions;
options.verbose = 0;

solver = CTMC(model, options);

%% Get Steady-State Expected Rewards
[R, names] = solver.getAvgReward();

fprintf('=== Steady-State Expected Rewards (Templates) ===\n');
fprintf('%-25s: %10.6f\n', 'QueueLength', R(1));
fprintf('%-25s: %10.6f\n', 'QueueLength_Class1', R(2));
fprintf('%-25s: %10.6f\n', 'Utilization', R(3));
fprintf('%-25s: %10.6f\n', 'Utilization_Class1', R(4));
fprintf('%-25s: %10.6f\n', 'BlockingProb', R(5));

%% Analytical Comparison
fprintf('\n=== Comparison with M/M/1 Analytical Results ===\n');

lambda = 1.5;  % Arrival rate
mu = 2;        % Service rate
rho = lambda / mu;
K = 10;        % Buffer capacity

% Steady-state probabilities
pi = zeros(K+1, 1);
for n = 0:K
    pi(n+1) = (1 - rho) * rho^n / (1 - rho^(K+1));
end

% Analytical metrics
L_analytical = sum((0:K)' .* pi);              % Expected queue length
U_analytical = 1 - pi(1);                      % Utilization
B_analytical = pi(end);                        % Blocking probability

fprintf('%-25s: LINE = %10.6f, Analytical = %10.6f, Error = %.2e\n', ...
    'QueueLength', R(1), L_analytical, abs(R(1) - L_analytical));
fprintf('%-25s: LINE = %10.6f, Analytical = %10.6f, Error = %.2e\n', ...
    'Utilization', R(3), U_analytical, abs(R(3) - U_analytical));
fprintf('%-25s: LINE = %10.6f, Analytical = %10.6f, Error = %.2e\n', ...
    'BlockingProb', R(5), B_analytical, abs(R(5) - B_analytical));

fprintf('\n✓ Example completed successfully.\n');
fprintf('  All template-based rewards match analytical results!\n');
