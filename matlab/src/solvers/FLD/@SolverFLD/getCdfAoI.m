function [AoI_cdf, PAoI_cdf] = getCdfAoI(self, t_values)
%GETCDFAOI Get CDF of Age of Information
%
% [AoI_cdf, PAoI_cdf] = GETCDFAOI(t_values)
%
% Returns the cumulative distribution function (CDF) of Age of Information
% and Peak AoI computed using matrix exponential representations.
%
% The CDF is computed using the formula:
%   F(t) = 1 - g * expm(A*t) * h
%
% where g, A, h are the matrix exponential parameters from the aoi-fluid solver.
%
% Parameters:
%   self: SolverFLD instance
%   t_values (optional): Vector of time values at which to evaluate CDF
%                        If not provided, uses automatic range based on mean
%
% Returns:
%   AoI_cdf: [n x 2] matrix with [CDF_values, t_values] for AoI
%   PAoI_cdf: [n x 2] matrix with [CDF_values, t_values] for Peak AoI
%
% Example:
%   model = Network('AoI_Example');
%   % ... model setup ...
%   solver = SolverFLD(model, 'method', 'mfq');
%   solver.getAvg();
%
%   % Get CDF at automatic time points
%   [AoI_cdf, PAoI_cdf] = solver.getCdfAoI();
%
%   % Plot CDFs
%   figure;
%   plot(AoI_cdf(:,2), AoI_cdf(:,1), 'b-', 'LineWidth', 2);
%   hold on;
%   plot(PAoI_cdf(:,2), PAoI_cdf(:,1), 'r--', 'LineWidth', 2);
%   xlabel('Age');
%   ylabel('CDF');
%   legend('AoI', 'Peak AoI');
%
%   % Get CDF at specific time points
%   t = linspace(0, 10, 100);
%   [AoI_cdf, PAoI_cdf] = solver.getCdfAoI(t);
%
% See also: getAvgAoI, solver_mfq_aoi, aoi_is_aoi

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Ensure solver has been run
if isempty(self.result)
    self.getAvg();
end

% Initialize outputs
AoI_cdf = [];
PAoI_cdf = [];

% Check if AoI results are available
if ~isfield(self.result, 'solverSpecific') || ...
   ~isfield(self.result.solverSpecific, 'aoiResults') || ...
   isempty(self.result.solverSpecific.aoiResults)
    line_warning(mfilename, 'No AoI results available. Ensure model has valid AoI topology and use method=''mfq''.');
    return;
end

aoiResults = self.result.solverSpecific.aoiResults;

% Check if matrix exponential parameters are available
if isempty(aoiResults.AoI_A) || isempty(aoiResults.AoI_g) || isempty(aoiResults.AoI_h)
    line_warning(mfilename, 'Matrix exponential parameters not available for AoI CDF computation.');
    return;
end

% Generate time values if not provided
if nargin < 2 || isempty(t_values)
    % Use automatic range based on mean AoI
    mean_aoi = aoiResults.AoI_mean;
    if isnan(mean_aoi) || mean_aoi <= 0
        mean_aoi = 1;  % Default if mean is invalid
    end
    % Range from 0 to 5 times the mean (covers most of the distribution)
    t_max = 5 * mean_aoi;
    t_values = linspace(0, t_max, 200);
end

t_values = t_values(:);  % Ensure column vector
n = length(t_values);

% Extract matrix exponential parameters for AoI
AoI_g = aoiResults.AoI_g;
AoI_A = aoiResults.AoI_A;
AoI_h = aoiResults.AoI_h;

% Compute AoI CDF: F(t) = 1 - g * expm(A*t) * h
AoI_F = zeros(n, 1);
for i = 1:n
    t = t_values(i);
    if t <= 0
        AoI_F(i) = 0;
    else
        % Compute matrix exponential
        expAt = expm(AoI_A * t);
        ccdf = AoI_g * expAt * AoI_h;
        AoI_F(i) = max(0, min(1, 1 - ccdf));
    end
end

AoI_cdf = [AoI_F, t_values];

% Extract matrix exponential parameters for Peak AoI
PAoI_g = aoiResults.PAoI_g;
PAoI_A = aoiResults.PAoI_A;
PAoI_h = aoiResults.PAoI_h;

% Compute Peak AoI CDF: F(t) = 1 - g * expm(A*t) * h
PAoI_F = zeros(n, 1);
for i = 1:n
    t = t_values(i);
    if t <= 0
        PAoI_F(i) = 0;
    else
        % Compute matrix exponential
        expAt = expm(PAoI_A * t);
        ccdf = PAoI_g * expAt * PAoI_h;
        PAoI_F(i) = max(0, min(1, 1 - ccdf));
    end
end

PAoI_cdf = [PAoI_F, t_values];

end
