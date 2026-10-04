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
%   Looking back from any instant, each earlier update is delivered iff its
%   service ends before the next arrival, so the age is a geometric sum of
%   Exp(lambda) gaps ending at the first delivered update:
%     A*(s)    = lambda*H*(s+lambda) / (s + lambda*H*(s+lambda))
%     E[A]     = 1/(lambda*H*(lambda))
%     E[Apeak] = -H*(lambda)'/H*(lambda) + 1/(lambda*H*(lambda))
%   Both reduce to 1/lambda + 1/mu and the Exp(lambda)*Exp(mu) product only
%   when the service is exponential.
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

% Mean AoI for preemptive LCFS: the delivery rate is lambda*H*(lambda)
meanAoI = 1 / (lambda * H_lst(lambda));

% Mean Peak AoI (exact)
% Success probability per arrival q = P(S < Y) = H*(lambda);
% E[S|success] = -H*'(lambda)/H*(lambda);
% E[Apeak] = E[S|success] + 1/(lambda*H*(lambda)) (validated by DES)
hstep = 1e-6 * max(1, lambda);
dHstar = (H_lst(lambda + hstep) - H_lst(lambda - hstep)) / (2 * hstep);
peakAoI = -dHstar / H_lst(lambda) + 1 / (lambda * H_lst(lambda));

% LST of AoI for preemptive LCFS; A*(0) = 1 identically
lstAoI = @(s) (lambda .* H_lst(s + lambda)) ./ (s + lambda .* H_lst(s + lambda));

end
