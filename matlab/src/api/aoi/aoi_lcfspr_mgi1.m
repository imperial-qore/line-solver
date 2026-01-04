function [meanAoI, lstAoI, peakAoI] = aoi_lcfspr_mgi1(lambda, H_lst, E_H, E_H2)
%AOI_LCFSPR_MGI1 Mean AoI and LST for M/GI/1 preemptive LCFS queue
%
% [meanAoI, lstAoI, peakAoI] = aoi_lcfspr_mgi1(lambda, H_lst, E_H, E_H2)
%
% Computes the Age of Information metrics for an M/GI/1 queue with
% preemptive Last-Come First-Served (LCFS-PR) discipline.
%
% In LCFS-PR, when a new update arrives, it preempts the current update
% in service (if any). This ensures the freshest update is always served.
%
% Parameters:
%   lambda (double): Arrival rate (Poisson arrivals)
%   H_lst (function_handle): LST of service time, @(s) -> complex
%   E_H (double): Mean service time (first moment)
%   E_H2 (double): Second moment of service time (for reference)
%
% Returns:
%   meanAoI (double): Mean (average) Age of Information
%   lstAoI (function_handle): LST of AoI distribution
%   peakAoI (double): Mean Peak Age of Information
%
% Formulas (from Inoue et al., IEEE Trans. IT, 2019, Section IV):
%   For preemptive LCFS, the AoI simplifies significantly:
%     E[A] = E[Y] + E[H] = 1/lambda + E[H]
%     E[Apeak] = E[Y] + E[H]
%
%   LST of AoI: A*(s) = (lambda/(s+lambda)) * H*(s)
%   This is the convolution of exponential interarrival and service time.
%
% Reference:
%   Y. Inoue, H. Masuyama, T. Takine, T. Tanaka, "A General Formula for
%   the Stationary Distribution of the Age of Information and Its
%   Application to Single-Server Queues," IEEE Trans. Information Theory,
%   vol. 65, no. 12, pp. 8305-8324, 2019.
%
% See also: aoi_fcfs_mgi1, aoi_lcfspr_mm1, aoi_lcfspr_gim1

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Validate inputs
if lambda <= 0
    line_error(mfilename, 'Arrival rate lambda must be positive');
end
if E_H <= 0
    line_error(mfilename, 'Mean service time E_H must be positive');
end

% Compute utilization (for reference; preemptive LCFS still requires rho < 1)
rho = lambda * E_H;

if rho >= 1
    line_error(mfilename, 'System unstable: rho = lambda*E_H = %.4f >= 1', rho);
end

% Mean interarrival time
E_Y = 1 / lambda;

% Mean AoI for preemptive LCFS (Proposition 3)
% E[A] = E[Y] + E[H]
meanAoI = E_Y + E_H;

% Mean Peak AoI (same as mean AoI for preemptive LCFS)
peakAoI = E_Y + E_H;

% LST of AoI for preemptive LCFS
% A*(s) = (lambda / (s + lambda)) * H*(s)
% This is the product of exponential interarrival LST and service LST
lstAoI = @(s) (lambda ./ (s + lambda)) .* H_lst(s);

end
