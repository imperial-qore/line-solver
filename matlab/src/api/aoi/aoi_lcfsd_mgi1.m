function [meanAoI, lstAoI, peakAoI] = aoi_lcfsd_mgi1(lambda, H_lst, E_H, E_H2)
%AOI_LCFSD_MGI1 Mean AoI for M/GI/1 non-preemptive LCFS-D queue
%
% [meanAoI, lstAoI, peakAoI] = aoi_lcfsd_mgi1(lambda, H_lst, E_H, E_H2)
%
% Computes the Age of Information metrics for an M/GI/1 queue with
% non-preemptive Last-Come First-Served with Discarding (LCFS-D) discipline.
%
% In LCFS-D, when a new update arrives while the server is busy:
%   - If there's an update waiting in queue, it is discarded
%   - The new update takes its place in the queue
%   - When service completes, the waiting update (if any) is served
%
% This is also known as M/GI/1/2* (buffer size 2 with replacement).
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
% Formulas (Inoue et al., IEEE Trans. IT, 2019, Section 3.3, NP-LCFS (C)),
% with H*(s) the service LST:
%   E[A]     = (lambda*E[H^2]/2 + H*(lambda)/lambda - H*'(lambda))/(rho + H*(lambda))
%              + (1 - H*(lambda))/lambda + H*'(lambda) + E[H]          (eq. 68)
%   E[Apeak] = 1/lambda + H*'(lambda) + 2*E[H]
%   A*(s)    = eq. (64)
% Discarding keeps the system stable for every rho, so rho >= 1 is accepted.
%
% Reference:
%   Y. Inoue, H. Masuyama, T. Takine, T. Tanaka, "A General Formula for
%   the Stationary Distribution of the Age of Information and Its
%   Application to Single-Server Queues," IEEE Trans. Information Theory,
%   vol. 65, no. 12, pp. 8305-8324, 2019.
%
% See also: aoi_lcfss_mgi1, aoi_fcfs_mgi1, aoi_lcfspr_mgi1

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

% Discarding keeps at most one update waiting, so the age process is
% regenerative for every rho (Inoue et al. 2019, Section 3.3): no rho < 1 gate.
rho = lambda * E_H;

% H*(lambda) = P(no arrival during a service), -H*'(lambda) = E[H exp(-lambda H)]
hL = H_lst(lambda);
hstep = 1e-6 * max(1, lambda);
mdH = -(H_lst(lambda + hstep) - H_lst(lambda - hstep)) / (2 * hstep);

% Mean AoI, Inoue et al. 2019, Corollary 36(i), eq. (68)
meanAoI = (lambda * E_H2 / 2 + hL / lambda + mdH) / (rho + hL) ...
    + (1 - hL) / lambda - mdH + E_H;

% Mean peak AoI = E[D] + 1/lambda_dagger, with E[W] = (1-H*(lambda))/lambda - (-H*'(lambda))
% and 1/lambda_dagger = E[H] + H*(lambda)/lambda (eqs. (28), (89))
peakAoI = 1 / lambda - mdH + 2 * E_H;

% LST of the AoI, Theorem 35(i), eq. (64)
lstAoI = @(s) aoi_lcfsd_mgi1_lst_eval(s, lambda, H_lst, E_H, rho, hL);

end

function val = aoi_lcfsd_mgi1_lst_eval(s, lambda, H_lst, E_H, rho, hL)
%AOI_LCFSD_MGI1_LST_EVAL Evaluate eq. (64) at point(s) s
    val = zeros(size(s));
    for i = 1:numel(s)
        si = s(i);
        if abs(si) < 1e-12
            val(i) = 1; % removable singularity of the residual-service LST at s = 0
            continue
        end
        Hs = H_lst(si);
        HsL = H_lst(si + lambda);
        hres_s = (1 - Hs) / (si * E_H);                 % residual service LST at s
        hres_sL = (1 - HsL) / ((si + lambda) * E_H);    % residual service LST at s+lambda
        val(i) = (hL + rho * hres_sL) * Hs ...
            * (rho * hres_s + HsL * lambda / (si + lambda)) / (rho + hL);
    end
end
