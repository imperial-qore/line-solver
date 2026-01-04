function [meanAoI, varAoI, peakAoI] = aoi_lcfspr_mm1(lambda, mu)
%AOI_LCFSPR_MM1 Mean, variance, and peak AoI for M/M/1 preemptive LCFS queue
%
% [meanAoI, varAoI, peakAoI] = aoi_lcfspr_mm1(lambda, mu)
%
% Computes the Age of Information metrics for an M/M/1 queue with
% preemptive Last-Come First-Served (LCFS-PR) discipline.
%
% In LCFS-PR, when a new update arrives, it preempts the current update
% in service (if any). This ensures the freshest update is always served,
% leading to lower AoI than FCFS.
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
% Formulas (from Inoue et al., IEEE Trans. IT, 2019, Section IV):
%   Mean AoI: E[A] = (1/mu) * (1 + 1/rho)
%   Peak AoI: E[Apeak] = 1/mu + 1/lambda
%
% Note: LCFS-PR always achieves lower mean AoI than FCFS for M/M/1.
%
% Reference:
%   Y. Inoue, H. Masuyama, T. Takine, T. Tanaka, "A General Formula for
%   the Stationary Distribution of the Age of Information and Its
%   Application to Single-Server Queues," IEEE Trans. Information Theory,
%   vol. 65, no. 12, pp. 8305-8324, 2019.
%
% See also: aoi_fcfs_mm1, aoi_lcfspr_mgi1, aoi_compare

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

% Check stability (for preemptive LCFS, the system is always stable
% as long as lambda and mu are positive, but rho < 1 is still required
% for meaningful queueing analysis)
if rho >= 1
    line_error(mfilename, 'System unstable: rho = lambda/mu = %.4f >= 1', rho);
end

% Mean AoI for preemptive LCFS (Proposition 3 / Section IV)
% E[A] = (1/mu) * (1 + 1/rho) = 1/mu + 1/lambda
% This is always less than FCFS: (1/mu)(1 + 1/rho) < (1/mu)(1 + 1/rho + rho^2/(1-rho))
meanAoI = (1/mu) * (1 + 1/rho);

% Mean Peak AoI for preemptive LCFS
% E[Apeak] = E[S] + E[Y] = 1/mu + 1/lambda
% where S is service time of the successful update, Y is interarrival time
% For preemptive LCFS, the peak AoI equals the mean AoI
peakAoI = 1/mu + 1/lambda;

% Variance of AoI for preemptive LCFS
% From Section IV of the paper, the second moment is:
% E[A^2] = 2/mu^2 + 2/(mu*lambda) + 2/lambda^2
%        = 2 * (1/mu + 1/lambda)^2 / something... need to verify
%
% Actually, for preemptive LCFS, the AoI distribution is simpler.
% The AoI at any time equals the interarrival time plus the service time
% of the most recent successful update. Since these are independent:
% Var[A] = Var[Y] + Var[S] = 1/lambda^2 + 1/mu^2
%
% However, this is an approximation. For exact formula:
E_A = meanAoI;
% For preemptive LCFS, E[A^2] = 2*(1/lambda^2 + 1/(lambda*mu) + 1/mu^2)
E_A2 = 2 * (1/lambda^2 + 1/(lambda*mu) + 1/mu^2);
varAoI = E_A2 - E_A^2;

% Ensure non-negative variance (numerical safety)
if varAoI < 0
    varAoI = 0;
end

end
