function [meanAoI, lstAoI, peakAoI] = aoi_lcfss_gim1(Y_lst, mu, E_Y, E_Y2)
%AOI_LCFSS_GIM1 Mean AoI for GI/M/1 non-preemptive LCFS-S queue
%
% [meanAoI, lstAoI, peakAoI] = aoi_lcfss_gim1(Y_lst, mu, E_Y, E_Y2)
%
% Computes the Age of Information metrics for a GI/M/1 queue with
% non-preemptive Last-Come First-Served with Set-aside (LCFS-S) discipline.
%
% In LCFS-S, when a new update arrives while the server is busy:
%   - The new update waits in the queue
%   - When service completes, the most recent update in queue is served next
%   - Older updates stay in queue ("set aside") and are still served later,
%     although their delivery no longer reduces the age
%
% Parameters:
%   Y_lst (function_handle): LST of interarrival time, @(s) -> complex
%   mu (double): Service rate (exponential service)
%   E_Y (double): Mean interarrival time (first moment)
%   E_Y2 (double): Second moment of interarrival time
%
% Returns:
%   meanAoI (double): Mean (average) Age of Information
%   lstAoI (function_handle): LST of AoI distribution
%   peakAoI (double): Mean Peak Age of Information
%
% Formulas (Inoue et al., IEEE Trans. IT, 2019, Section 3.3, NP-LCFS (D),
% i.e. non-preemptive LCFS without discarding), with G*(s) the interarrival
% LST, rho = 1/(E_Y*mu) and gamma the root of G*(mu - mu*x) = x:
%   E[A]     = 1/mu + E[G^2]/(2E[G]) + rho*(-G*'(mu - mu*gamma))      (eq. 71)
%   E[Apeak] = P0*(E[G] + (1+G*(mu))/mu)
%              + Pw*(1/mu + (E[G] + G*'(mu) + (gamma-G*(mu))/(gamma*mu))/(1-G*(mu))),
%              P0 = (1-gamma)/(1-gamma*G*(mu)), Pw = 1-P0           (eqs. 100-103)
%   A*(s)    = eq. (67)
%
% Reference:
%   Y. Inoue, H. Masuyama, T. Takine, T. Tanaka, "A General Formula for
%   the Stationary Distribution of the Age of Information and Its
%   Application to Single-Server Queues," IEEE Trans. Information Theory,
%   vol. 65, no. 12, pp. 8305-8324, 2019.
%
% See also: aoi_lcfsd_gim1, aoi_fcfs_gim1, aoi_lcfspr_gim1

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

% Check stability: without discarding the backlog of stale updates must drain
if rho >= 1
    line_error(mfilename, 'System unstable: rho = 1/(E_Y*mu) = %.4f >= 1', rho);
end

% gamma: probability an arrival finds the server busy, root of Y*(mu - mu*x) = x
gamma_func = @(x) Y_lst(mu - mu * x) - x;
gam = fzero(gamma_func, [0.001, 0.999]);

gM = Y_lst(mu);
hstep = 1e-6 * max(1, mu);
mdY = -(Y_lst(mu + hstep) - Y_lst(mu - hstep)) / (2 * hstep);
zg = mu - mu * gam;
hz = 1e-6 * max(1, zg);
mdYg = -(Y_lst(zg + hz) - Y_lst(zg - hz)) / (2 * hz);

% Mean AoI, Inoue et al. 2019, Corollary 36(iv), eq. (71)
meanAoI = 1 / mu + E_Y2 / (2 * E_Y) + rho * mdYg;

% Mean peak AoI from eqs. (100), (102), (103): an informative update finds the
% system empty w.p. P0, otherwise it waited behind a residual Exp(mu) service
P0 = (1 - gam) / (1 - gam * gM);
Pw = gam * (1 - gM) / (1 - gam * gM);
E0 = E_Y + (1 + gM) / mu;
Ew = 1 / mu + (E_Y - mdY + (gam - gM) / (gam * mu)) / (1 - gM);
peakAoI = P0 * E0 + Pw * Ew;

% LST of the AoI, Theorem 35(iv), eq. (67)
lstAoI = @(s) aoi_lcfss_gim1_lst_eval(s, Y_lst, mu, E_Y, rho, gam);

end

function val = aoi_lcfss_gim1_lst_eval(s, Y_lst, mu, E_Y, rho, gam)
%AOI_LCFSS_GIM1_LST_EVAL Evaluate eq. (67) at point(s) s
    val = zeros(size(s));
    for i = 1:numel(s)
        si = s(i);
        if abs(si) < 1e-12
            val(i) = 1; % removable singularity of the residual-interarrival LST at s = 0
            continue
        end
        gres = (1 - Y_lst(si)) / (si * E_Y);
        val(i) = (gres + rho * (Y_lst(si + mu - mu * gam) - gam) * mu / (si + mu)) * mu / (si + mu);
    end
end
