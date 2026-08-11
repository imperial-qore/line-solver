function [AoI, PAoI, aoiTable] = getAvgAoI(self)
%GETAVGAOI Get average Age of Information metrics
%
% [AoI, PAoI, aoiTable] = GETAVGAOI()
%
% Returns Age of Information (AoI) and Peak AoI statistics computed by
% the Fluid solver using the aoi-fluid MFQ method.
%
% This method requires the model to have a valid AoI topology:
% - Single open class
% - Source -> Queue -> Sink topology
% - Queue capacity 1 (bufferless) or 2 (single-buffer)
% - Single server
% - FCFS, LCFS, or LCFSPR scheduling
%
% To use this method, first run the solver with method='mfq':
%   solver = SolverFLD(model, 'method', 'mfq');
%   solver.getAvg();
%   [AoI, PAoI, aoiTable] = solver.getAvgAoI();
%
% Parameters:
%   self: SolverFLD instance
%
% Returns:
%   AoI (struct): Age of Information statistics
%       .mean: Mean AoI
%       .var: Variance of AoI
%       .std: Standard deviation of AoI
%   PAoI (struct): Peak Age of Information statistics
%       .mean: Mean Peak AoI
%       .var: Variance of Peak AoI
%       .std: Standard deviation of Peak AoI
%   aoiTable (table): Summary table with all metrics
%
% Example:
%   model = Network('AoI_Example');
%   source = Source(model, 'Source');
%   queue = Queue(model, 'Queue', SchedStrategy.FCFS);
%   queue.setCapacity(1);  % Bufferless
%   sink = Sink(model, 'Sink');
%   jobclass = OpenClass(model, 'Updates');
%   source.setArrival(jobclass, Exp(1));
%   queue.setService(jobclass, Exp(2));
%   model.link(Network.serialRouting(source, queue, sink));
%
%   solver = SolverFLD(model, 'method', 'mfq');
%   [AoI, PAoI, aoiTable] = solver.getAvgAoI();
%   fprintf('Mean AoI: %.4f, Mean PAoI: %.4f\n', AoI.mean, PAoI.mean);
%
% See also: getCdfAoI, solver_mfq_aoi, aoi_is_aoi

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Ensure solver has been run
if isempty(self.result)
    self.getAvg();
end

% Initialize outputs
AoI = struct('mean', NaN, 'var', NaN, 'std', NaN);
PAoI = struct('mean', NaN, 'var', NaN, 'std', NaN);
aoiTable = table();

% Check if AoI results are available
if ~isfield(self.result, 'solverSpecific') || ...
   ~isfield(self.result.solverSpecific, 'aoiResults') || ...
   isempty(self.result.solverSpecific.aoiResults)
    line_warning(mfilename, 'No AoI results available. Ensure model has valid AoI topology and use method=''mfq''.');
    return;
end

aoiResults = self.result.solverSpecific.aoiResults;

% Extract AoI statistics
AoI.mean = aoiResults.AoI_mean;
AoI.var = aoiResults.AoI_var;
AoI.std = sqrt(max(0, AoI.var));

% Extract Peak AoI statistics
PAoI.mean = aoiResults.PAoI_mean;
PAoI.var = aoiResults.PAoI_var;
PAoI.std = sqrt(max(0, PAoI.var));

% Create summary table
Metric = {'AoI'; 'Peak AoI'};
Mean = [AoI.mean; PAoI.mean];
Variance = [AoI.var; PAoI.var];
StdDev = [AoI.std; PAoI.std];
SystemType = {aoiResults.systemType; aoiResults.systemType};
Preemption = [aoiResults.preemption; aoiResults.preemption];

aoiTable = table(Metric, Mean, Variance, StdDev, SystemType, Preemption);

end
