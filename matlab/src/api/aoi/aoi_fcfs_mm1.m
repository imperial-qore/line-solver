function [meanAoI, varAoI, peakAoI] = aoi_fcfs_mm1(lambda, mu)
%AOI_FCFS_MM1 Mean, variance, and peak AoI for M/M/1 FCFS queue
%
% [meanAoI, varAoI, peakAoI] = aoi_fcfs_mm1(lambda, mu)
%
% Computes the Age of Information metrics for an M/M/1 queue with
% First-Come First-Served (FCFS) discipline.
%
% Parameters:
%   lambda (double): Arrival rate (Poisson arrivals)
%   mu (double): Service rate (exponential service)
%
% Returns:
%   meanAoI (double): Mean (average) Age of Information
%   varAoI (double): Variance of Age of Information
%   peakAoI (double): Mean Peak Age of Information
%
% Formulas (from Inoue et al., IEEE Trans. IT, 2019):
%   Mean AoI: E[A] = (1/mu) * (1 + 1/rho + rho^2/(1-rho))
%   Peak AoI: E[Apeak] = (1/mu) * (1 + 1/rho + rho/(1-rho))
%
% Reference:
%   Y. Inoue, H. Masuyama, T. Takine, T. Tanaka, "A General Formula for
%   the Stationary Distribution of the Age of Information and Its
%   Application to Single-Server Queues," IEEE Trans. Information Theory,
%   vol. 65, no. 12, pp. 8305-8324, 2019.
%
% See also: aoi_lcfspr_mm1, aoi_fcfs_mgi1, aoi_optimal_rate

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Validate inputs
if lambda <= 0
    line_error(mfilename, 'Arrival rate lambda must be positive');
end
if mu <= 0
    line_error(mfilename, 'Service rate mu must be positive');
end

% Compute utilization
rho = lambda / mu;

% Check stability
if rho >= 1
    line_error(mfilename, 'System unstable: rho = lambda/mu = %.4f >= 1', rho);
end

% Mean AoI (Proposition 1 / Corollary from Section III-A)
% E[A] = (1/mu) * (1 + 1/rho + rho^2/(1-rho))
%      = (1/mu) * (1 + mu/lambda + rho^2/(1-rho))
meanAoI = (1/mu) * (1 + 1/rho + rho^2/(1-rho));

% Mean Peak AoI
% E[Apeak] = E[T] + E[Y] where T is system time, Y is interarrival time
% For M/M/1: E[T] = 1/(mu-lambda), E[Y] = 1/lambda
% E[Apeak] = 1/(mu-lambda) + 1/lambda = (1/mu) * (1/rho + 1/(1-rho))
%          = (1/mu) * (1 + 1/rho + rho/(1-rho))
peakAoI = (1/mu) * (1 + 1/rho + rho/(1-rho));

% Variance of AoI
% From the paper, the second moment can be derived from the LST.
% For M/M/1 FCFS, using results from Section III-A:
% E[A^2] can be computed from -d^2/ds^2 A*(s)|_{s=0}
%
% The variance is computed using the known formula:
% Var[A] = E[A^2] - (E[A])^2
%
% For M/M/1 FCFS, using the closed-form results:
% Second moment: E[A^2] = 2/mu^2 * (1/rho^2 + 1/rho + rho^2/(1-rho)^2 + 1/(1-rho))
E_A = meanAoI;
E_A2 = (2/mu^2) * (1/rho^2 + 1/rho + rho^2/(1-rho)^2 + 1/(1-rho) + rho/(1-rho));
varAoI = E_A2 - E_A^2;

% Ensure non-negative variance (numerical safety)
if varAoI < 0
    varAoI = 0;
end

end
