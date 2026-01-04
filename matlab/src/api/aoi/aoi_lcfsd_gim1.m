function [meanAoI, lstAoI, peakAoI] = aoi_lcfsd_gim1(Y_lst, mu, E_Y, E_Y2)
%AOI_LCFSD_GIM1 Mean AoI for GI/M/1 non-preemptive LCFS-D queue
%
% [meanAoI, lstAoI, peakAoI] = aoi_lcfsd_gim1(Y_lst, mu, E_Y, E_Y2)
%
% Computes the Age of Information metrics for a GI/M/1 queue with
% non-preemptive Last-Come First-Served with Discarding (LCFS-D) discipline.
%
% In LCFS-D, when a new update arrives while the server is busy:
%   - If there's an update waiting in queue, it is discarded
%   - The new update takes its place in the queue
%   - When service completes, the waiting update (if any) is served
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
% See also: aoi_lcfss_gim1, aoi_fcfs_gim1, aoi_lcfspr_gim1

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

% Find sigma: probability arriving customer finds server busy
sigma_func = @(sig) Y_lst(mu - mu * sig) - sig;
try
    sigma = fzero(sigma_func, [0.001, 0.999]);
catch
    sigma = rho;  % Fallback
end

% For GI/M/1 LCFS-D (also GI/M/1/2*), when an update arrives:
% - If server idle: enters service immediately
% - If server busy: waits (discarding any previous waiting update)
%
% The effective service time is S + (residual service if busy)

% Probability arriving update finds server busy: sigma
% When busy, expected remaining service = 1/mu (memoryless)

% Mean effective system time
E_T_eff = E_S + sigma * E_S;

% Mean AoI for LCFS-D (from Section VI analysis adapted for GI/M/1)
% E[A] = E[Y] + E[T_eff] + additional term for waiting
meanAoI = E_Y + E_S * (1 + sigma) + sigma * E_S / (1 + sigma);

% Mean Peak AoI
peakAoI = E_Y + E_T_eff;

% LST is complex for LCFS-D; return empty handle
lstAoI = [];

end
