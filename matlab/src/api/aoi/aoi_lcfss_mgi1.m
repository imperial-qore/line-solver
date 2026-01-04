function [meanAoI, lstAoI, peakAoI] = aoi_lcfss_mgi1(lambda, H_lst, E_H, E_H2)
%AOI_LCFSS_MGI1 Mean AoI for M/GI/1 non-preemptive LCFS-S queue
%
% [meanAoI, lstAoI, peakAoI] = aoi_lcfss_mgi1(lambda, H_lst, E_H, E_H2)
%
% Computes the Age of Information metrics for an M/GI/1 queue with
% non-preemptive Last-Come First-Served with Set-aside (LCFS-S) discipline.
%
% In LCFS-S, when a new update arrives while the server is busy:
%   - The new update waits in the queue
%   - When service completes, the most recent update in queue is served next
%   - The older update remains in queue (is "set aside")
%
% Parameters:
%   lambda (double): Arrival rate (Poisson arrivals)
%   H_lst (function_handle): LST of service time, @(s) -> complex
%   E_H (double): Mean service time (first moment)
%   E_H2 (double): Second moment of service time
%
% Returns:
%   meanAoI (double): Mean (average) Age of Information
%   lstAoI (function_handle): LST of AoI distribution (empty if not computed)
%   peakAoI (double): Mean Peak Age of Information
%
% Formulas (from Inoue et al., IEEE Trans. IT, 2019, Section V):
%   The analysis is more complex than FCFS or preemptive LCFS.
%   Mean AoI is computed using the results from Theorem 6.
%
% Reference:
%   Y. Inoue, H. Masuyama, T. Takine, T. Tanaka, "A General Formula for
%   the Stationary Distribution of the Age of Information and Its
%   Application to Single-Server Queues," IEEE Trans. Information Theory,
%   vol. 65, no. 12, pp. 8305-8324, 2019.
%
% See also: aoi_lcfsd_mgi1, aoi_fcfs_mgi1, aoi_lcfspr_mgi1

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Validate inputs
if lambda <= 0
    line_error(mfilename, 'Arrival rate lambda must be positive');
end
if E_H <= 0
    line_error(mfilename, 'Mean service time E_H must be positive');
end
if E_H2 < E_H^2
    line_error(mfilename, 'Second moment E_H2 must be >= E_H^2');
end

% Compute utilization
rho = lambda * E_H;

% Check stability
if rho >= 1
    line_error(mfilename, 'System unstable: rho = lambda*E_H = %.4f >= 1', rho);
end

% Mean interarrival time
E_Y = 1 / lambda;

% For LCFS-S, the mean AoI is given by Proposition 5 in the paper.
% The formula involves the busy period distribution of M/G/1.
%
% Key quantities:
% - B: busy period of M/G/1, with LST B*(s) satisfying B*(s) = H*(s + lambda - lambda*B*(s))
% - E[B] = E[H] / (1 - rho)
% - E[B^2] = E[H^2] / (1 - rho)^3

E_B = E_H / (1 - rho);
E_B2 = E_H2 / (1 - rho)^3;

% Mean AoI for LCFS-S (Proposition 5, simplified form)
% E[A] = E[Y] + E[H] + (1 - rho) * lambda * E[B^2] / 2
%      = 1/lambda + E[H] + (1 - rho) * lambda * E[H^2] / (2 * (1 - rho)^3)
%      = 1/lambda + E[H] + lambda * E[H^2] / (2 * (1 - rho)^2)

meanAoI = E_Y + E_H + lambda * E_H2 / (2 * (1 - rho)^2);

% Mean Peak AoI
% For LCFS-S, peak AoI involves the residual busy period
peakAoI = E_Y + E_H + lambda * E_B2 / (2 * E_B);

% LST is complex for LCFS-S; return empty handle
lstAoI = [];

end
