function [meanAoI, lstAoI, peakAoI] = aoi_lcfss_mgi1(lambda, H_lst, E_H, E_H2)
%AOI_LCFSS_MGI1 Mean AoI for M/GI/1 non-preemptive LCFS-S queue
%
% [meanAoI, lstAoI, peakAoI] = aoi_lcfss_mgi1(lambda, H_lst, E_H, E_H2)
%
% Computes the Age of Information metrics for an M/GI/1 queue with
% non-preemptive Last-Come First-Served with Set-aside (LCFS-S) discipline.
%
% In LCFS-S, when a new update arrives while the server is busy:
%   - The new update waits in the queue
%   - When service completes, the most recent update in queue is served next
%   - Older updates stay in queue ("set aside") and are still served later,
%     although their delivery no longer reduces the age
%
% Parameters:
%   lambda (double): Arrival rate (Poisson arrivals)
%   H_lst (function_handle): LST of service time, @(s) -> complex
%   E_H (double): Mean service time (first moment)
%   E_H2 (double): Second moment of service time
%
% Returns:
%   meanAoI (double): Mean (average) Age of Information
%   lstAoI (function_handle): LST of AoI distribution
%   peakAoI (double): Mean Peak Age of Information
%
% Formulas (Inoue et al., IEEE Trans. IT, 2019, Section 3.3, NP-LCFS (D),
% i.e. non-preemptive LCFS without discarding), with H*(s) the service LST:
%   E[A]     = lambda*E[H^2]/2 + ((1-rho)^2/(rho*H*(lambda)) + 2)*E[H]     (eq. 70)
%   E[Apeak] = E[H] + E[W] + 1/(lambda*(2 - rho - H*(lambda)))
%   E[W]     = ((1 - H*(lambda))/lambda + H*'(lambda))/(2 - rho - H*(lambda))
%   A*(s)    = eq. (66)
%
% Reference:
%   Y. Inoue, H. Masuyama, T. Takine, T. Tanaka, "A General Formula for
%   the Stationary Distribution of the Age of Information and Its
%   Application to Single-Server Queues," IEEE Trans. Information Theory,
%   vol. 65, no. 12, pp. 8305-8324, 2019.
%
% See also: aoi_lcfsd_mgi1, aoi_fcfs_mgi1, aoi_lcfspr_mgi1

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Validate inputs
if lambda <= 0
    line_error(mfilename, 'Arrival rate lambda must be positive');
end
if E_H <= 0
    line_error(mfilename, 'Mean service time E_H must be positive');
end
if E_H2 < E_H^2
    line_error(mfilename, 'Second moment E_H2 must be >= E_H^2');
end

% Compute utilization
rho = lambda * E_H;

% Check stability: without discarding the backlog of stale updates must drain
if rho >= 1
    line_error(mfilename, 'System unstable: rho = lambda*E_H = %.4f >= 1', rho);
end

% H*(lambda) = P(no arrival during a service), -H*'(lambda) = E[H exp(-lambda H)]
hL = H_lst(lambda);
hstep = 1e-6 * max(1, lambda);
mdH = -(H_lst(lambda + hstep) - H_lst(lambda - hstep)) / (2 * hstep);

% Mean AoI, Inoue et al. 2019, Corollary 36(iii), eq. (70)
meanAoI = lambda * E_H2 / 2 + ((1 - rho)^2 / (rho * hL) + 2) * E_H;

% Mean peak AoI = E[D] + 1/lambda_dagger (eq. 28). An arrival is informative if it
% finds the server idle or no later arrival precedes the end of the current
% service, so lambda_dagger = lambda*(1 - rho + rho*Hres*(lambda)) = lambda*(2 - rho - H*(lambda)),
% and E[W] follows from eq. (98).
den = 2 - rho - hL;
E_W = ((1 - hL) / lambda - mdH) / den;
peakAoI = E_H + E_W + 1 / (lambda * den);

% LST of the AoI, Theorem 35(iii), eq. (66)
lstAoI = @(s) aoi_lcfss_mgi1_lst_eval(s, lambda, H_lst, E_H, rho);

end

function val = aoi_lcfss_mgi1_lst_eval(s, lambda, H_lst, E_H, rho)
%AOI_LCFSS_MGI1_LST_EVAL Evaluate eq. (66) at point(s) s
    val = zeros(size(s));
    for i = 1:numel(s)
        si = s(i);
        if abs(si) < 1e-12
            val(i) = 1; % removable singularity of the residual-service LST at s = 0
            continue
        end
        Hs = H_lst(si);
        HsL = H_lst(si + lambda);
        hres_s = (1 - Hs) / (si * E_H);
        val(i) = lambda / (si + lambda) * Hs * (rho * hres_s ...
            + (1 - rho) * (si + lambda) * (1 - Hs + HsL) / (si + lambda * HsL));
    end
end
