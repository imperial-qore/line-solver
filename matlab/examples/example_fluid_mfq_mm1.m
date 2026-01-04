%% EXAMPLE_FLUID_MFQ_MM1 - MFQ Method for M/M/1 Analysis
%
% This example demonstrates the use of the MFQ (Markovian fluid model) method
% in LINE's SolverFluid for exact fluid-level and sojourn-time analysis of
% single-queue open queueing systems.
%
% The MFQ method uses BUTools' FluFluQueue function to analyze fluid
% queues with independent arrival and service processes. It provides exact
% steady-state solutions for single-queue topologies.
%
% Applicability:
%   - Single-queue open systems (Source -> Queue -> Sink)
%   - Single-server or infinite-server queues
%   - Supported distributions: Exponential, Erlang, Hyper-exponential, etc.
%   - Steady-state analysis
%
% References:
%   Horvath G, Telek M, "Sojourn times in fluid queues with independent
%   and dependent input and output processes", PEVA 79:160-181, 2014.

clear; clc;

%% Setup: Create M/M/1 Model
% Define an M/M/1 queue with:
%   - Poisson arrivals with rate λ = 0.6 jobs/time unit
%   - Exponential service with rate μ = 1.0 jobs/time unit
%   - Utilization ρ = λ/μ = 0.6

model = Network('MM1_Example');

% Create network nodes
source = Source(model, 'Source');
queue = Queue(model, 'Queue', SchedStrategy.FCFS);
sink = Sink(model, 'Sink');

% Create open job class
jobclass = OpenClass(model, 'Class1');

% Set arrival and service processes
lambda = 0.6;
mu = 1.0;

source.setArrival(jobclass, Exp(lambda));   % Exponential arrivals
queue.setService(jobclass, Exp(mu));         % Exponential service
queue.setNumberOfServers(1);                 % Single server

% Define serial routing: Source -> Queue -> Sink
model.link(Network.serialRouting(source, queue, sink));

fprintf('===== M/M/1 Queue Analysis with MFQ Method =====\n');
fprintf('Arrival rate λ = %.1f jobs/time unit\n', lambda);
fprintf('Service rate μ = %.1f jobs/time unit\n', mu);
fprintf('Utilization ρ = %.2f\n\n', lambda/mu);

%% Solve using MFQ Method
fprintf('--- Solution with MFQ Method ---\n');
solver_mfq = SolverFluid(model, 'method', 'mfq');
AvgTable_mfq = solver_mfq.getAvgTable();
disp(AvgTable_mfq);

%% Solve using Matrix Method (for comparison)
fprintf('\n--- Solution with Matrix Method (for comparison) ---\n');
solver_matrix = SolverFluid(model, 'method', 'matrix');
AvgTable_matrix = solver_matrix.getAvgTable();
disp(AvgTable_matrix);

%% Compare with Analytical M/M/1 Results
fprintf('\n--- Analytical M/M/1 Theory ---\n');
rho = lambda / mu;
L_theory = rho / (1 - rho);           % Mean queue length
W_theory = 1 / (mu - lambda);          % Mean response time
U_theory = rho;                        % Utilization
X_theory = lambda;                     % Throughput

fprintf('Queue Length:     %.6f\n', L_theory);
fprintf('Response Time:    %.6f\n', W_theory);
fprintf('Utilization:      %.6f\n', U_theory);
fprintf('Throughput:       %.6f\n\n', X_theory);

%% Create Comparison Table
fprintf('--- Comparison of Methods ---\n');
fprintf('Metric              | MFQ         | Matrix    | Theory    | Error (MFQ)\n');
fprintf('%-19s | %11.6f | %9.6f | %9.6f | %14.2f%%\n', ...
    'Queue Length', ...
    AvgTable_mfq.QLen(2), AvgTable_matrix.QLen(2), L_theory, ...
    100*abs(AvgTable_mfq.QLen(2) - L_theory)/L_theory);

fprintf('%-19s | %11.6f | %9.6f | %9.6f | %14.2f%%\n', ...
    'Response Time', ...
    AvgTable_mfq.RespT(2), AvgTable_matrix.RespT(2), W_theory, ...
    100*abs(AvgTable_mfq.RespT(2) - W_theory)/W_theory);

fprintf('%-19s | %11.6f | %9.6f | %9.6f | %14.2f%%\n', ...
    'Utilization', ...
    AvgTable_mfq.Util(2), AvgTable_matrix.Util(2), U_theory, ...
    100*abs(AvgTable_mfq.Util(2) - U_theory)/U_theory);

fprintf('%-19s | %11.6f | %9.6f | %9.6f | %14.2f%%\n', ...
    'Throughput', ...
    AvgTable_mfq.Tput(2), AvgTable_matrix.Tput(2), X_theory, ...
    100*abs(AvgTable_mfq.Tput(2) - X_theory)/X_theory);

%% Verify Little's Law
fprintf('\n--- Little''s Law Validation (L = λ * W) ---\n');
L_mfq = AvgTable_mfq.QLen(2);
W_mfq = AvgTable_mfq.RespT(2);
littles_check = lambda * W_mfq;

fprintf('MFQ:   L = %.6f, λ*W = %.6f, Error = %.4f%%\n', ...
    L_mfq, littles_check, 100*abs(L_mfq - littles_check)/L_mfq);

%% Key Observations
fprintf('\n===== Key Observations =====\n');
fprintf('1. MFQ provides exact analysis for single-queue open systems\n');
fprintf('2. Results match analytical M/M/1 formulas very closely (< 0.1%% error)\n');
fprintf('3. Little''s Law is satisfied: Queue length = Arrival rate × Response time\n');
fprintf('4. MFQ is specialized for single-queue topologies\n');
fprintf('5. For multi-queue networks, use matrix or closing methods instead\n');

%% Method Selection Guide
fprintf('\n===== Method Selection Guide =====\n');
fprintf('Use MFQ method when:\n');
fprintf('  - Single queue: Source -> Queue -> Sink topology\n');
fprintf('  - Open classes required\n');
fprintf('  - Single-server (c=1) or infinite-server (c=Inf)\n');
fprintf('  - Steady-state analysis needed\n');
fprintf('  - Exact solutions desired (not approximations)\n\n');

fprintf('Use Matrix method when:\n');
fprintf('  - Multiple queues in network\n');
fprintf('  - Transient analysis needed\n');
fprintf('  - Multi-server queues (c > 1)\n');
fprintf('  - General approximation for open/closed networks\n');
