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
%   For preemptive LCFS with GI/M/1 the age is the backward recurrence time
%   of the arrivals, uniform on (0,tau) here, plus an independent Exp(mu):
%     E[A]   = tau/2 + 1/mu
%     Var[A] = tau^2/12 + 1/mu^2
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

% Mean AoI for preemptive LCFS: backward recurrence U(0,tau) plus Exp(mu)
meanAoI = tau / 2 + 1 / mu;

% Mean Peak AoI (exact)
% Success probability per arrival q = P(S < tau) = 1 - exp(-mu*tau);
% E[S|success] = (1/mu - tau*exp(-mu*tau) - exp(-mu*tau)/mu)/q;
% E[Apeak] = E[S|success] + tau/q (validated by DES)
q = 1 - exp(-mu * tau);
ES_succ = (1/mu - tau * exp(-mu * tau) - exp(-mu * tau)/mu) / q;
peakAoI = ES_succ + tau / q;

% Variance of AoI for preemptive LCFS: Var[U(0,tau)] + Var[Exp(mu)]
varAoI = tau^2 / 12 + 1 / mu^2;

end
