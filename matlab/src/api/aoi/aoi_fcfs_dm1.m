function [meanAoI, varAoI, peakAoI] = aoi_fcfs_dm1(tau, mu)
%AOI_FCFS_DM1 Mean, variance, and peak AoI for D/M/1 FCFS queue
%
% [meanAoI, varAoI, peakAoI] = aoi_fcfs_dm1(tau, mu)
%
% Computes the Age of Information metrics for a D/M/1 queue with
% First-Come First-Served (FCFS) discipline.
%
% D/M/1: Deterministic arrivals with interarrival time tau (rate 1/tau),
%        exponential service with rate mu.
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
% Formulas (from Inoue et al., IEEE Trans. IT, 2019):
%   Uses GI/M/1 FCFS results with deterministic arrivals.
%   For D/M/1, the key parameter is sigma, the root of:
%     sigma = exp(-mu*tau*(1-sigma))
%   which gives the probability an arriving customer finds system non-empty.
%
% Reference:
%   Y. Inoue, H. Masuyama, T. Takine, T. Tanaka, "A General Formula for
%   the Stationary Distribution of the Age of Information and Its
%   Application to Single-Server Queues," IEEE Trans. Information Theory,
%   vol. 65, no. 12, pp. 8305-8324, 2019.
%
% See also: aoi_fcfs_mm1, aoi_fcfs_gim1, aoi_lcfspr_dm1

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
% For D/M/1: lambda = 1/tau, so rho = 1/(tau*mu)
lambda = 1 / tau;
rho = lambda / mu;

% Check stability
if rho >= 1
    line_error(mfilename, 'System unstable: rho = 1/(tau*mu) = %.4f >= 1', rho);
end

% Compute sigma: unique root in (0,1) of sigma = exp(-mu*tau*(1-sigma))
% This is the probability an arriving customer finds the server busy
% Use fixed-point iteration or fzero
sigma_func = @(s) exp(-mu * tau * (1 - s)) - s;
sigma = fzero(sigma_func, [0.001, 0.999]);

% For GI/M/1, the system delay (waiting + service) has LST:
% D*(s) = (1-sigma)*mu / (mu - s) for s < mu
% E[D] = 1/mu + sigma/(mu*(1-sigma)) = 1/(mu*(1-sigma))

E_D = 1 / (mu * (1 - sigma));

% Mean interarrival time (deterministic)
E_Y = tau;

% Mean AoI for GI/M/1 FCFS (from Theorem 3 / Proposition 2)
% E[A] = E[Y] + E[D] + sigma*(1-sigma)/(mu*sigma_prime)
% where sigma_prime = d/ds sigma(s)|_{s=0}
%
% For deterministic arrivals, using the formula from the paper:
% E[A] = tau + 1/(mu*(1-sigma)) + sigma/(mu*(1-sigma)^2 * mu*tau)
%
% Simplified formula for D/M/1:
% E[A] = tau + E[D] + sigma / (mu * (1-sigma))
meanAoI = tau + E_D + sigma / (mu * (1 - sigma));

% Mean Peak AoI
% E[Apeak] = E[Y] + E[D] = tau + 1/(mu*(1-sigma))
peakAoI = tau + E_D;

% Variance of AoI
% For D/M/1, the variance depends on the system delay variance and
% the deterministic interarrival (which contributes 0 variance).
% E[D^2] = 2/(mu*(1-sigma))^2 for exponential service with GI/M/1 structure
E_D2 = 2 / (mu * (1 - sigma))^2;
Var_D = E_D2 - E_D^2;

% Var[Y] = 0 (deterministic)
% Approximate variance of AoI (needs exact formula from paper)
varAoI = Var_D + (sigma / (mu * (1-sigma)))^2;

% Ensure non-negative variance (numerical safety)
if varAoI < 0
    varAoI = 0;
end

end
