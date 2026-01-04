function [meanAoI, varAoI, peakAoI] = aoi_lcfspr_dm1(tau, mu)
%AOI_LCFSPR_DM1 Mean, variance, and peak AoI for D/M/1 preemptive LCFS queue
%
% [meanAoI, varAoI, peakAoI] = aoi_lcfspr_dm1(tau, mu)
%
% Computes the Age of Information metrics for a D/M/1 queue with
% preemptive Last-Come First-Served (LCFS-PR) discipline.
%
% In LCFS-PR, when a new update arrives, it preempts the current update
% in service (if any). For D/M/1, arrivals are deterministic (every tau
% time units) and service is exponential with rate mu.
%
% Parameters:
%   tau (double): Deterministic interarrival time
%   mu (double): Service rate (exponential service)
%
% Returns:
%   meanAoI (double): Mean (average) Age of Information
%   varAoI (double): Variance of Age of Information
%   peakAoI (double): Mean Peak Age of Information
%
% Formulas (from Inoue et al., IEEE Trans. IT, 2019, Section IV):
%   For preemptive LCFS with GI/M/1:
%     E[A] = E[Y] + E[S] = tau + 1/mu
%   where S is the service time (exponential).
%
% Reference:
%   Y. Inoue, H. Masuyama, T. Takine, T. Tanaka, "A General Formula for
%   the Stationary Distribution of the Age of Information and Its
%   Application to Single-Server Queues," IEEE Trans. Information Theory,
%   vol. 65, no. 12, pp. 8305-8324, 2019.
%
% See also: aoi_fcfs_dm1, aoi_lcfspr_mm1, aoi_lcfspr_gim1

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Validate inputs
if tau <= 0
    line_error(mfilename, 'Interarrival time tau must be positive');
end
if mu <= 0
    line_error(mfilename, 'Service rate mu must be positive');
end

% Compute utilization
lambda = 1 / tau;
rho = lambda / mu;

% Check stability (rho < 1 is required)
if rho >= 1
    line_error(mfilename, 'System unstable: rho = 1/(tau*mu) = %.4f >= 1', rho);
end

% For preemptive LCFS, the AoI is determined by the interarrival time
% plus the service time of the successful update.
%
% Mean interarrival time (deterministic)
E_Y = tau;

% Mean service time (exponential)
E_S = 1 / mu;

% Mean AoI for preemptive LCFS (Proposition 4 / Section IV)
% E[A] = E[Y] + E[S]
% For D/M/1-PR: E[A] = tau + 1/mu
meanAoI = E_Y + E_S;

% Mean Peak AoI (same formula for preemptive LCFS)
% E[Apeak] = E[Y] + E[S]
peakAoI = E_Y + E_S;

% Variance of AoI for preemptive LCFS
% Var[A] = Var[Y] + Var[S]
% Since Y is deterministic, Var[Y] = 0
% Var[S] = 1/mu^2 (exponential)
varAoI = 1 / mu^2;

end
