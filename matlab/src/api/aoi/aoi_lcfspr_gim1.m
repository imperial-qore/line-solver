function [meanAoI, lstAoI, peakAoI] = aoi_lcfspr_gim1(Y_lst, mu, E_Y, E_Y2)
%AOI_LCFSPR_GIM1 Mean AoI and LST for GI/M/1 preemptive LCFS queue
%
% [meanAoI, lstAoI, peakAoI] = aoi_lcfspr_gim1(Y_lst, mu, E_Y, E_Y2)
%
% Computes the Age of Information metrics for a GI/M/1 queue with
% preemptive Last-Come First-Served (LCFS-PR) discipline.
%
% In LCFS-PR, when a new update arrives, it preempts the current update
% in service (if any). This ensures the freshest update is always served.
%
% Parameters:
%   Y_lst (function_handle): LST of interarrival time, @(s) -> complex
%   mu (double): Service rate (exponential service)
%   E_Y (double): Mean interarrival time (first moment)
%   E_Y2 (double): Second moment of interarrival time (for reference)
%
% Returns:
%   meanAoI (double): Mean (average) Age of Information
%   lstAoI (function_handle): LST of AoI distribution
%   peakAoI (double): Mean Peak Age of Information
%
% Formulas (from Inoue et al., IEEE Trans. IT, 2019, Section IV):
%   For preemptive LCFS, the AoI simplifies significantly:
%     E[A] = E[Y] + E[S] = E[Y] + 1/mu
%     E[Apeak] = E[Y] + 1/mu
%
%   LST of AoI: A*(s) = Y*(s) * (mu / (s + mu))
%   This is the product of interarrival LST and exponential service LST.
%
% Reference:
%   Y. Inoue, H. Masuyama, T. Takine, T. Tanaka, "A General Formula for
%   the Stationary Distribution of the Age of Information and Its
%   Application to Single-Server Queues," IEEE Trans. Information Theory,
%   vol. 65, no. 12, pp. 8305-8324, 2019.
%
% See also: aoi_fcfs_gim1, aoi_lcfspr_mm1, aoi_lcfspr_mgi1

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

if rho >= 1
    line_error(mfilename, 'System unstable: rho = 1/(E_Y*mu) = %.4f >= 1', rho);
end

% Mean service time
E_S = 1 / mu;

% Mean AoI for preemptive LCFS (Proposition 4)
% E[A] = E[Y] + E[S]
meanAoI = E_Y + E_S;

% Mean Peak AoI (same as mean AoI for preemptive LCFS)
peakAoI = E_Y + E_S;

% LST of AoI for preemptive LCFS
% A*(s) = Y*(s) * (mu / (s + mu))
% This is the product of interarrival LST and exponential service LST
lstAoI = @(s) Y_lst(s) .* (mu ./ (s + mu));

end
