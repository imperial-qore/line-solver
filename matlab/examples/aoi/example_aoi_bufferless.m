%% Example: Age of Information analysis for bufferless queue (PH/PH/1/1)
%
% This example demonstrates AoI analysis using the Fluid solver with
% the aoi-fluid MFQ method for bufferless (capacity-1) queues.
%
% The model is: Source -> Queue(cap=1) -> Sink
% with Erlang arrivals and Erlang service.
%
% This matches the bufferless example from the aoi-fluid toolbox:
%   - Arrival: Erlang(3,3) with mean 1
%   - Service: Erlang(2,2) with mean 1
%
% Reference:
%   aoi-fluid toolbox by Ozancan Dogan, Nail Akar, Eray Unsal Atay
%   BSD 2-Clause License, 2020

clear;
fprintf('=== Age of Information: Bufferless Queue (PH/PH/1/1) ===\n\n');

%% Model Construction
model = Network('AoI_Bufferless');

% Create nodes
source = Source(model, 'Source');
queue = Queue(model, 'Queue', SchedStrategy.FCFS);
queue.setNumberOfServers(1);
sink = Sink(model, 'Sink');

% Create open class for status updates
jobclass = OpenClass(model, 'Updates');

% Set distributions
% Erlang(3,3) arrival: mean = 3/3 = 1
source.setArrival(jobclass, Erlang(3, 3));

% Erlang(2,2) service: mean = 2/2 = 1
queue.setService(jobclass, Erlang(2, 2));

% Set bufferless capacity (only one job can be in system at a time)
queue.setCapacity(1);

% Link nodes
model.link(Network.serialRouting(source, queue, sink));

%% Solve with MFQ method (triggers AoI analysis)
fprintf('Solving with Fluid solver (method=mfq)...\n');
solver = SolverFLD(model, 'method', 'mfq');

% Get standard performance metrics
AvgTable = solver.getAvgTable();
fprintf('\nPerformance Metrics:\n');
disp(AvgTable);

% Get AoI metrics
[AoI, PAoI, aoiTable] = solver.getAvgAoI();
fprintf('\nAge of Information Results:\n');
disp(aoiTable);

fprintf('Mean AoI:  %.6f\n', AoI.mean);
fprintf('Var AoI:   %.6f\n', AoI.var);
fprintf('Mean PAoI: %.6f\n', PAoI.mean);
fprintf('Var PAoI:  %.6f\n', PAoI.var);

%% Get and plot AoI CDF
[AoI_cdf, PAoI_cdf] = solver.getCdfAoI();

if ~isempty(AoI_cdf)
    figure;
    plot(AoI_cdf(:,2), AoI_cdf(:,1), 'b-', 'LineWidth', 2);
    hold on;
    plot(PAoI_cdf(:,2), PAoI_cdf(:,1), 'r--', 'LineWidth', 2);
    xlabel('Age');
    ylabel('CDF');
    legend('AoI', 'Peak AoI', 'Location', 'southeast');
    title('Age of Information CDF - Bufferless Queue (PH/PH/1/1)');
    grid on;
end

%% Compare with preemptive LCFS (PH/PH/1/1*)
fprintf('\n\n=== Comparing FCFS vs Preemptive LCFS ===\n');

% Create preemptive LCFS model
model_pr = Network('AoI_Bufferless_Preemptive');
source_pr = Source(model_pr, 'Source');
queue_pr = Queue(model_pr, 'Queue', SchedStrategy.LCFSPR);
queue_pr.setNumberOfServers(1);
sink_pr = Sink(model_pr, 'Sink');
jobclass_pr = OpenClass(model_pr, 'Updates');
source_pr.setArrival(jobclass_pr, Erlang(3, 3));
queue_pr.setService(jobclass_pr, Erlang(2, 2));
queue_pr.setCapacity(1);
model_pr.link(Network.serialRouting(source_pr, queue_pr, sink_pr));

solver_pr = SolverFLD(model_pr, 'method', 'mfq');
[AoI_pr, PAoI_pr, ~] = solver_pr.getAvgAoI();

fprintf('\nFCFS (PH/PH/1/1):\n');
fprintf('  Mean AoI:  %.6f\n', AoI.mean);
fprintf('  Mean PAoI: %.6f\n', PAoI.mean);

fprintf('\nPreemptive LCFS (PH/PH/1/1*):\n');
fprintf('  Mean AoI:  %.6f\n', AoI_pr.mean);
fprintf('  Mean PAoI: %.6f\n', PAoI_pr.mean);

fprintf('\nPreemption reduces Mean AoI by %.1f%%\n', ...
    100 * (AoI.mean - AoI_pr.mean) / AoI.mean);
