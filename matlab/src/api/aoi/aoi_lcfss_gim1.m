function [meanAoI, lstAoI, peakAoI] = aoi_lcfss_gim1(Y_lst, mu, E_Y, E_Y2)
%AOI_LCFSS_GIM1 Mean AoI for GI/M/1 non-preemptive LCFS-S queue
%
% [meanAoI, lstAoI, peakAoI] = aoi_lcfss_gim1(Y_lst, mu, E_Y, E_Y2)
%
% Computes the Age of Information metrics for a GI/M/1 queue with
% non-preemptive Last-Come First-Served with Set-aside (LCFS-S) discipline.
%
% In LCFS-S, when a new update arrives while the server is busy:
%   - The new update waits in the queue
%   - When service completes, the most recent update in queue is served next
%   - The older update remains in queue (is "set aside")
%
% Parameters:
%   Y_lst (function_handle): LST of interarrival time, @(s) -> complex
%   mu (double): Service rate (exponential service)
%   E_Y (double): Mean interarrival time (first moment)
%   E_Y2 (double): Second moment of interarrival time
%
% Returns:
%   meanAoI (double): Mean (average) Age of Information
%   lstAoI (function_handle): LST of AoI distribution (empty if not computed)
%   peakAoI (double): Mean Peak Age of Information
%
% Reference:
%   Y. Inoue, H. Masuyama, T. Takine, T. Tanaka, "A General Formula for
%   the Stationary Distribution of the Age of Information and Its
%   Application to Single-Server Queues," IEEE Trans. Information Theory,
%   vol. 65, no. 12, pp. 8305-8324, 2019.
%
% See also: aoi_lcfsd_gim1, aoi_fcfs_gim1, aoi_lcfspr_gim1

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Validate inputs
if mu <= 0
    line_error(mfilename, 'Service rate mu must be positive');
end
if E_Y <= 0
    line_error(mfilename, 'Mean interarrival time E_Y must be positive');
end

% Compute utilization
lambda = 1 / E_Y;
rho = lambda / mu;

% Check stability
if rho >= 1
    line_error(mfilename, 'System unstable: rho = 1/(E_Y*mu) = %.4f >= 1', rho);
end

% Mean service time
E_S = 1 / mu;
E_S2 = 2 / mu^2;  % Second moment of exponential

% Find sigma: probability arriving customer finds server busy
sigma_func = @(sig) Y_lst(mu - mu * sig) - sig;
try
    sigma = fzero(sigma_func, [0.001, 0.999]);
catch
    sigma = rho;  % Fallback
end

% For GI/M/1 LCFS-S, the analysis is similar to M/GI/1 but with
% roles of arrival and service distributions swapped.
%
% Mean system delay
E_D = 1 / (mu * (1 - sigma));

% Mean AoI for LCFS-S (from Section V analysis adapted for GI/M/1)
% E[A] = E[Y] + E[S] + sigma * E[D]
meanAoI = E_Y + E_S + sigma * E_D;

% Mean Peak AoI
peakAoI = E_Y + E_D;

% LST is complex for LCFS-S; return empty handle
lstAoI = [];

end
