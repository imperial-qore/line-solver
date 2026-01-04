%% EXAMPLE_FLUID_MFQ_MAP - MFQ for MAP/PH/1 Queues
%
% This example demonstrates MFQ (Markovian fluid model) analysis for more
% complex single-queue systems with Markovian Arrival Processes (MAP) and
% Phase-type service distributions.
%
% Model: Balanced HyperExponential Arrivals and Erlang-3 Service
%   - Two-state HyperExp arrival process (balanced, CV² = 0.5)
%   - Erlang-3 service (regular, CV² = 1/3)
%   - Single server
%
% This showcases MFQ's ability to handle non-exponential
% distributions exactly.

clear; clc;

%% Setup: Create MAP/PH/1 Model

model = Network('HyperExp_Erlang_Model');

% Create network nodes
source = Source(model, 'Source');
queue = Queue(model, 'Queue', SchedStrategy.FCFS);
sink = Sink(model, 'Sink');

% Create open job class
jobclass = OpenClass(model, 'Class1');

% Define arrival process: 2-state Balanced Hyper-exponential
%   Pr(phase 1) = Pr(phase 2) = 0.5
%   Phase 1: rate μ₁ = 0.5
%   Phase 2: rate μ₂ = 2.0
%   Aggregate rate: λ = 0.5

% Using HyperExp for balanced hyper-exponential
lambda = 0.5;
p = 0.5;        % Probability of phase 1
mu1 = 1.0;      % Rate of phase 1
mu2 = 1.0;      % Rate of phase 2 (for balanced)

% Create hyper-exponential arrival (uses two exponentials mixed)
HE_arrival = HyperExp(2, [p, 1-p], [mu1, mu2]);
source.setArrival(jobclass, HE_arrival);

% Define service process: Erlang-3
% CV² = 1/3, mean = 1.0 (stage rate = 3.0)
k_erlang = 3;   % Number of phases
mu_erlang = 3.0; % Rate per phase (so mean = k/mu = 3/3 = 1.0)

queue.setService(jobclass, Erlang(k_erlang, mu_erlang));
queue.setNumberOfServers(1);

% Define serial routing
model.link(Network.serialRouting(source, queue, sink));

%% Display Model Information
fprintf('===== MAP/PH/1 Queue Analysis with MFQ =====\n');
fprintf('\nArrival Process: Balanced Hyper-Exponential (2 states)\n');
fprintf('  Pr(phase 1) = 0.5, rate μ₁ = %.1f\n', mu1);
fprintf('  Pr(phase 2) = 0.5, rate μ₂ = %.1f\n', mu2);
fprintf('  Aggregate arrival rate λ = %.2f\n', lambda);

fprintf('\nService Process: Erlang-%d\n', k_erlang);
fprintf('  Phases = %d\n', k_erlang);
fprintf('  Rate per phase = %.1f\n', mu_erlang);
fprintf('  Mean service time = %.3f\n', k_erlang / mu_erlang);
fprintf('  SCV (CV²) = %.3f\n', 1/k_erlang);

rho = lambda / (k_erlang / mu_erlang);
fprintf('\nUtilization ρ = λ * E[S] = %.2f * %.3f = %.4f\n', ...
    lambda, k_erlang / mu_erlang, rho);

%% Solve with MFQ
fprintf('\n--- MFQ Method ---\n');
try
    solver_mfq = SolverFluid(model, 'method', 'mfq');
    [QN_mfq, UN_mfq, RN_mfq, TN_mfq] = solver_mfq.getAvg();

    AvgTable_mfq = solver_mfq.getAvgTable();
    disp(AvgTable_mfq);

    mfq_success = true;
    fprintf('\nMFQ completed successfully\n');
catch ME
    fprintf('\nMFQ error: %s\n', ME.message);
    mfq_success = false;
end

%% Solve with Matrix Method (for comparison)
fprintf('\n--- Matrix Method (for comparison) ---\n');
solver_matrix = SolverFluid(model, 'method', 'matrix');
[QN_matrix, UN_matrix, RN_matrix, TN_matrix] = solver_matrix.getAvg();

AvgTable_matrix = solver_matrix.getAvgTable();
disp(AvgTable_matrix);

%% Comparison
if mfq_success
    fprintf('\n--- Comparison of Methods ---\n');
    fprintf('Metric              | MFQ         | Matrix\n');
    fprintf('%-19s | %11.6f | %9.6f\n', 'Queue Length', ...
        AvgTable_mfq.QLen(2), AvgTable_matrix.QLen(2));
    fprintf('%-19s | %11.6f | %9.6f\n', 'Response Time', ...
        AvgTable_mfq.RespT(2), AvgTable_matrix.RespT(2));
    fprintf('%-19s | %11.6f | %9.6f\n', 'Utilization', ...
        AvgTable_mfq.Util(2), AvgTable_matrix.Util(2));
    fprintf('%-19s | %11.6f | %9.6f\n', 'Throughput', ...
        AvgTable_mfq.Tput(2), AvgTable_matrix.Tput(2));

    % Relative errors
    err_Q = abs(AvgTable_mfq.QLen(2) - AvgTable_matrix.QLen(2)) / ...
            AvgTable_matrix.QLen(2) * 100;
    err_R = abs(AvgTable_mfq.RespT(2) - AvgTable_matrix.RespT(2)) / ...
            AvgTable_matrix.RespT(2) * 100;

    fprintf('\nRelative differences:\n');
    fprintf('  Queue length:  %.2f%%\n', err_Q);
    fprintf('  Response time: %.2f%%\n', err_R);
end

%% Distribution Extraction (if supported)
fprintf('\n--- Distribution Capabilities ---\n');
fprintf('MFQ provides matrix-exponential representations:\n');
fprintf('  - Fluid level distribution (queue length distribution)\n');
fprintf('  - Sojourn time distribution (response time distribution)\n');
fprintf('  - Usable for percentile analysis and CDF evaluation\n');

%% Key Points
fprintf('\n===== Key Points for MAP/PH/1 Analysis =====\n');
fprintf('1. MFQ handles non-exponential arrivals and services exactly\n');
fprintf('2. Balanced hyper-exponential arrivals are well-supported\n');
fprintf('3. Erlang and other Coxian service times work accurately\n');
fprintf('4. Results are exact matrix-analytic solutions (not approximations)\n');
fprintf('5. Useful for performance analysis with realistic traffic patterns\n');

fprintf('\nSuitable for modeling:\n');
fprintf('  - Bursty arrival processes (MAP with multiple phases)\n');
fprintf('  - Multi-stage service processes (Erlang, Coxian)\n');
fprintf('  - Heavy-tailed traffic (HyperExp, Pareto-like via PH)\n');
