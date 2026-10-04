function [meanAoI, lstAoI, peakAoI] = aoi_lcfspr_gim1(Y_lst, mu, E_Y, E_Y2)
%AOI_LCFSPR_GIM1 Mean AoI and LST for GI/M/1 preemptive LCFS queue
%
% [meanAoI, lstAoI, peakAoI] = aoi_lcfspr_gim1(Y_lst, mu, E_Y, E_Y2)
%
% Computes the Age of Information metrics for a GI/M/1 queue with
% preemptive Last-Come First-Served (LCFS-PR) discipline.
%
% In LCFS-PR, when a new update arrives, it preempts the current update
% in service (if any). This ensures the freshest update is always served.
%
% Parameters:
%   Y_lst (function_handle): LST of interarrival time, @(s) -> complex
%   mu (double): Service rate (exponential service)
%   E_Y (double): Mean interarrival time (first moment)
%   E_Y2 (double): Second moment of interarrival time (for reference)
%
% Returns:
%   meanAoI (double): Mean (average) Age of Information
%   lstAoI (function_handle): LST of AoI distribution
%   peakAoI (double): Mean Peak Age of Information
%
% Formulas (from Inoue et al., IEEE Trans. IT, 2019, Section IV):
%   With memoryless service the latest delivered update is the latest arrival
%   before the last Exp(mu) service epoch, so the age is the backward
%   recurrence time of the arrival process plus an independent Exp(mu):
%     A*(s)    = (mu/(s+mu)) * lambda*(1 - Y*(s))/s
%     E[A]     = lambda*E[Y^2]/2 + 1/mu
%     E[Apeak] = E[S|success] + 1/(lambda*(1-Y*(mu)))
%   Both reduce to E[Y] + 1/mu and Y*(s)*mu/(s+mu) only for Poisson arrivals.
%
% Reference:
%   Y. Inoue, H. Masuyama, T. Takine, T. Tanaka, "A General Formula for
%   the Stationary Distribution of the Age of Information and Its
%   Application to Single-Server Queues," IEEE Trans. Information Theory,
%   vol. 65, no. 12, pp. 8305-8324, 2019.
%
% See also: aoi_fcfs_gim1, aoi_lcfspr_mm1, aoi_lcfspr_mgi1

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Validate inputs
if mu <= 0
    line_error(mfilename, 'Service rate mu must be positive');
end
if E_Y <= 0
    line_error(mfilename, 'Mean interarrival time E_Y must be positive');
end
if E_Y2 < E_Y^2
    line_error(mfilename, 'Second moment E_Y2 must be >= E_Y^2');
end

% Compute utilization
lambda = 1 / E_Y;
rho = lambda / mu;

if rho >= 1
    line_error(mfilename, 'System unstable: rho = 1/(E_Y*mu) = %.4f >= 1', rho);
end

% Mean AoI for preemptive LCFS: mean backward recurrence time plus 1/mu
meanAoI = lambda * E_Y2 / 2 + 1 / mu;

% Mean Peak AoI (exact)
% Success probability per arrival q = P(S < Y) = 1 - Y*(mu);
% E[S|success] = (1/mu + Y*'(mu) - Y*(mu)/mu)/q;
% E[Apeak] = E[S|success] + 1/(lambda*q) (validated by DES)
hstep = 1e-6 * max(1, mu);
dYstar = (Y_lst(mu + hstep) - Y_lst(mu - hstep)) / (2 * hstep);
q = 1 - Y_lst(mu);
ES_succ = (1/mu + dYstar - Y_lst(mu)/mu) / q;
peakAoI = ES_succ + 1 / (lambda * q);

% LST of AoI for preemptive LCFS
lstAoI = @(s) aoi_lcfspr_gim1_lst_eval(s, Y_lst, mu, lambda);

end

function val = aoi_lcfspr_gim1_lst_eval(s, Y_lst, mu, lambda)
%AOI_LCFSPR_GIM1_LST_EVAL Evaluate GI/M/1 LCFS-PR AoI LST at point(s) s
    val = zeros(size(s));
    for i = 1:numel(s)
        si = s(i);
        if abs(si) < 1e-12
            val(i) = 1; % removable singularity: lambda*(1-Y*(s))/s -> 1
            continue
        end
        val(i) = (mu / (si + mu)) * lambda * (1 - Y_lst(si)) / si;
    end
end
