%% Example: Age of Information analysis for single-buffer queue (M/PH/1/2)
%
% This example demonstrates AoI analysis using the Fluid solver with
% the aoi-fluid MFQ method for single-buffer (capacity-2) queues.
%
% The model is: Source -> Queue(cap=2) -> Sink
% with Poisson arrivals and Erlang service.
%
% This matches the single-buffer example from the aoi-fluid toolbox:
%   - Arrival: Poisson with rate 5
%   - Service: Erlang(2,2) with mean 1
%
% Reference:
%   aoi-fluid toolbox by Ozancan Dogan, Nail Akar, Eray Unsal Atay
%   BSD 2-Clause License, 2020

clear;
fprintf('=== Age of Information: Single-Buffer Queue (M/PH/1/2) ===\n\n');

%% Model Construction
model = Network('AoI_SingleBuffer');

% Create nodes
source = Source(model, 'Source');
queue = Queue(model, 'Queue', SchedStrategy.FCFS);
queue.setNumberOfServers(1);
sink = Sink(model, 'Sink');

% Create open class for status updates
jobclass = OpenClass(model, 'Updates');

% Set distributions
% Poisson arrivals with rate 5
source.setArrival(jobclass, Exp(5));

% Erlang(2,2) service: mean = 2/2 = 1
queue.setService(jobclass, Erlang(2, 2));

% Set single-buffer capacity (1 in service + 1 waiting)
queue.setCapacity(2);

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
    title('Age of Information CDF - Single-Buffer Queue (M/PH/1/2)');
    grid on;
end

%% Compare with LCFS replacement policy (M/PH/1/2*)
fprintf('\n\n=== Comparing FCFS vs LCFS with Replacement ===\n');

% Create LCFS model (replacement policy)
model_lcfs = Network('AoI_SingleBuffer_LCFS');
source_lcfs = Source(model_lcfs, 'Source');
queue_lcfs = Queue(model_lcfs, 'Queue', SchedStrategy.LCFS);
queue_lcfs.setNumberOfServers(1);
sink_lcfs = Sink(model_lcfs, 'Sink');
jobclass_lcfs = OpenClass(model_lcfs, 'Updates');
source_lcfs.setArrival(jobclass_lcfs, Exp(5));
queue_lcfs.setService(jobclass_lcfs, Erlang(2, 2));
queue_lcfs.setCapacity(2);
model_lcfs.link(Network.serialRouting(source_lcfs, queue_lcfs, sink_lcfs));

solver_lcfs = SolverFLD(model_lcfs, 'method', 'mfq');
[AoI_lcfs, PAoI_lcfs, ~] = solver_lcfs.getAvgAoI();

fprintf('\nFCFS (M/PH/1/2):\n');
fprintf('  Mean AoI:  %.6f\n', AoI.mean);
fprintf('  Mean PAoI: %.6f\n', PAoI.mean);

fprintf('\nLCFS with Replacement (M/PH/1/2*):\n');
fprintf('  Mean AoI:  %.6f\n', AoI_lcfs.mean);
fprintf('  Mean PAoI: %.6f\n', PAoI_lcfs.mean);

fprintf('\nLCFS replacement reduces Mean AoI by %.1f%%\n', ...
    100 * (AoI.mean - AoI_lcfs.mean) / AoI.mean);
