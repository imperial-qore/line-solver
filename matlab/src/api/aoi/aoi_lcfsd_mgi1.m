function [meanAoI, lstAoI, peakAoI] = aoi_lcfsd_mgi1(lambda, H_lst, E_H, E_H2)
%AOI_LCFSD_MGI1 Mean AoI for M/GI/1 non-preemptive LCFS-D queue
%
% [meanAoI, lstAoI, peakAoI] = aoi_lcfsd_mgi1(lambda, H_lst, E_H, E_H2)
%
% Computes the Age of Information metrics for an M/GI/1 queue with
% non-preemptive Last-Come First-Served with Discarding (LCFS-D) discipline.
%
% In LCFS-D, when a new update arrives while the server is busy:
%   - If there's an update waiting in queue, it is discarded
%   - The new update takes its place in the queue
%   - When service completes, the waiting update (if any) is served
%
% This is also known as M/GI/1/2* (buffer size 2 with replacement).
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
% Formulas (from Inoue et al., IEEE Trans. IT, 2019, Section VI):
%   The analysis uses the concept of effective service time which includes
%   potential waiting while another update is being served.
%
% Reference:
%   Y. Inoue, H. Masuyama, T. Takine, T. Tanaka, "A General Formula for
%   the Stationary Distribution of the Age of Information and Its
%   Application to Single-Server Queues," IEEE Trans. Information Theory,
%   vol. 65, no. 12, pp. 8305-8324, 2019.
%
% See also: aoi_lcfss_mgi1, aoi_fcfs_mgi1, aoi_lcfspr_mgi1

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

% For LCFS-D (also known as M/GI/1/2*), the analysis from Section VI:
%
% The effective service time S_eff includes potential residual service
% of the update being served when the successful update arrived.
%
% Key quantity: probability p0 that arriving update finds server idle
% p0 = 1 - rho (for Poisson arrivals)
%
% Mean effective system time:
% E[T_eff] = E[H] + (1 - p0) * E[H_residual]
%          = E[H] + rho * E[H^2] / (2 * E[H])

E_H_residual = E_H2 / (2 * E_H);  % Residual service time
E_T_eff = E_H + rho * E_H_residual;

% Mean AoI for LCFS-D (Proposition 6, simplified form)
% E[A] = E[Y] + E[T_eff] + additional term for discarded updates
%
% The additional term accounts for the age accumulated while waiting:
% Additional = rho * (E[H_residual] + E[H]) / (1 + rho)
%
% Simplified approximation:
meanAoI = E_Y + E_H + rho * E_H2 / (2 * E_H) + rho * E_H / (1 + rho);

% Mean Peak AoI
% Peak AoI is the system time plus interarrival time
peakAoI = E_Y + E_T_eff;

% LST is complex for LCFS-D; return empty handle
lstAoI = [];

end
