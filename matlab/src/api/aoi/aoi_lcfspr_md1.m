function [meanAoI, varAoI, peakAoI] = aoi_lcfspr_md1(lambda, d)
%AOI_LCFSPR_MD1 Mean, variance, and peak AoI for M/D/1 preemptive LCFS queue
%
% [meanAoI, varAoI, peakAoI] = aoi_lcfspr_md1(lambda, d)
%
% Computes the Age of Information metrics for an M/D/1 queue with
% preemptive Last-Come First-Served (LCFS-PR) discipline.
%
% In LCFS-PR, when a new update arrives, it preempts the current update
% in service (if any). For M/D/1, an update is successful if no new
% arrivals occur during its deterministic service time d.
%
% Parameters:
%   lambda (double): Arrival rate (Poisson arrivals)
%   d (double): Deterministic service time
%
% Returns:
%   meanAoI (double): Mean (average) Age of Information
%   varAoI (double): Variance of Age of Information
%   peakAoI (double): Mean Peak Age of Information
%
% Formulas (from Inoue et al., IEEE Trans. IT, 2019, Section IV):
%   For preemptive LCFS with M/G/1:
%     E[A] = E[Y] + E[S] = 1/lambda + d
%   where S is the (successful) service time, which equals d.
%
% Reference:
%   Y. Inoue, H. Masuyama, T. Takine, T. Tanaka, "A General Formula for
%   the Stationary Distribution of the Age of Information and Its
%   Application to Single-Server Queues," IEEE Trans. Information Theory,
%   vol. 65, no. 12, pp. 8305-8324, 2019.
%
% See also: aoi_fcfs_md1, aoi_lcfspr_mm1, aoi_lcfspr_mgi1

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Validate inputs
if lambda <= 0
    line_error(mfilename, 'Arrival rate lambda must be positive');
end
if d <= 0
    line_error(mfilename, 'Service time d must be positive');
end

% Compute utilization (for reference, but preemptive LCFS is always stable)
rho = lambda * d;

% Check stability (rho < 1 is still required for meaningful analysis)
if rho >= 1
    line_error(mfilename, 'System unstable: rho = lambda*d = %.4f >= 1', rho);
end

% For preemptive LCFS, the AoI is determined by the interarrival time
% plus the service time of the successful update.
%
% Mean interarrival time
E_Y = 1 / lambda;

% Mean service time (deterministic)
E_S = d;

% Mean AoI for preemptive LCFS (Proposition 3 / Section IV)
% E[A] = E[Y] + E[S]
% For M/D/1-PR: E[A] = 1/lambda + d
meanAoI = E_Y + E_S;

% Mean Peak AoI (same formula for preemptive LCFS)
% E[Apeak] = E[Y] + E[S]
peakAoI = E_Y + E_S;

% Variance of AoI for preemptive LCFS
% Var[A] = Var[Y] + Var[S]
% Since S is deterministic, Var[S] = 0
% Var[Y] = 1/lambda^2 (exponential)
varAoI = 1 / lambda^2;

end
