%% Example: Fork-Join Percentile Analysis using SolverMAM
%
% This example demonstrates automatic Fork-Join percentile analysis
% using SolverMAM with FJ_codes integration. When a valid Fork-Join
% topology is detected, SolverMAM automatically uses FJ_codes to
% compute response time percentiles.
%
% Reference:
% Z. Qiu, J.F. Pérez, and P. Harrison, "Beyond the Mean in Fork-Join Queues:
% Efficient Approximation for Response-Time Tails", IFIP Performance 2015.
% Copyright 2015 Imperial College London

clear;
clc;

%% Create a simple Fork-Join model
% This creates the structure: Source → Fork → K Queues → Join → Sink

model = Network('FJ_Percentile_Example');

% Source with exponential arrivals
source = Source(model, 'Source');

% Fork node
fork = Fork(model, 'Fork');

% K parallel queues with exponential service
K = 10;  % Number of parallel queues
mu = 1.0;  % Service rate (per queue)
queues = cell(1, K);

for i = 1:K
    queues{i} = Queue(model, sprintf('Queue%d', i), SchedStrategy.FCFS);
end

% Join node (paired with Fork)
join = Join(model, 'Join', fork);

% Sink
sink = Sink(model, 'Sink');

% Create an open class (after Sink)
jobClass = OpenClass(model, 'Class1');
lambda = 0.5;  % Arrival rate
source.setArrival(jobClass, Exp(lambda));

for i = 1:K
    queues{i}.setService(jobClass, Exp(mu));
end

% Link the network: Source → Fork → Queues → Join → Sink
% Create routing matrix
P = model.initRoutingMatrix();
P{1}(source, fork) = 1.0;
for i = 1:K
    P{1}(fork, queues{i}) = 1.0;
    P{1}(queues{i}, join) = 1.0;
end
P{1}(join, sink) = 1.0;

model.link(P);

%% Solve with SolverMAM (automatically detects FJ and uses FJ_codes)

fprintf('Creating SolverMAM...\n');
solver = SolverMAM(model);

fprintf('Running analysis (FJ topology will be auto-detected)...\n');
solver.runAnalyzer();

%% Display average metrics

fprintf('\n=== Average Performance Metrics ===\n');
avgTable = solver.getAvgTable();
disp(avgTable);

%% Display percentile metrics

fprintf('\n=== Response Time Percentiles ===\n');
% Request specific percentiles: 50th, 90th, 95th, 99th
percentiles = [50, 90, 95, 99];
[percRT, percTable] = solver.getPerctRespT(percentiles);

disp(percTable);

% Display detailed percentile information
fprintf('\nDetailed Percentile Results:\n');
for i = 1:length(percRT)
    fprintf('  Class: %s\n', percRT(i).class);
    fprintf('  K (parallel queues): %d\n', percRT(i).K);
    fprintf('  Method: %s\n', percRT(i).method);
    for p = 1:length(percRT(i).percentiles)
        fprintf('    %.2f percentile: %.4f\n', ...
            percRT(i).percentiles(p), percRT(i).values(p));
    end
end

%% Compare with theoretical values for exponential case
% For exponential service in Fork-Join, we can compute theoretical bounds

fprintf('\n=== Theoretical Comparison ===\n');
fprintf('System parameters:\n');
fprintf('  Arrival rate (lambda): %.2f\n', lambda);
fprintf('  Service rate per queue (mu): %.2f\n', mu);
fprintf('  Number of parallel queues (K): %d\n', K);
fprintf('  Utilization per queue: %.2f\n', lambda/mu);
if lambda < mu*K
    fprintf('  Stability: Stable\n');
else
    fprintf('  Stability: Unstable\n');
end

% Mean response time approximation for exponential FJ
mean_rt_approx = sum(1./(mu * (1:K)));
fprintf('  Approximate mean response time: %.4f\n', mean_rt_approx);

fprintf('\nNote: FJ_codes provides accurate percentile approximations\n');
fprintf('      for the synchronization delay in Fork-Join systems.\n');
