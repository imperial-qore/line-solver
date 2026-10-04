function [meanAoI, lstAoI, peakAoI] = aoi_fcfs_mgi1(lambda, H_lst, E_H, E_H2)
%AOI_FCFS_MGI1 Mean AoI and LST for M/GI/1 FCFS queue
%
% [meanAoI, lstAoI, peakAoI] = aoi_fcfs_mgi1(lambda, H_lst, E_H, E_H2)
%
% Computes the Age of Information metrics for an M/GI/1 queue with
% First-Come First-Served (FCFS) discipline.
%
% M/GI/1: Poisson arrivals with rate lambda, general independent service.
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
% Formulas (from Inoue et al., IEEE Trans. IT, 2019, Theorem 2):
%   LST of AoI: A*(s) = (lambda * H*(s)) / (s + lambda - lambda*H*(s)) * W*(s)
%   where W*(s) is the LST of waiting time (Pollaczek-Khinchine).
%
%   Mean AoI (Proposition 1):
%   E[A] = E[H] + E[T] + (1-2*rho)/lambda - d/ds T*(s)|_{s=lambda}  (exact)
%
% Reference:
%   Y. Inoue, H. Masuyama, T. Takine, T. Tanaka, "A General Formula for
%   the Stationary Distribution of the Age of Information and Its
%   Application to Single-Server Queues," IEEE Trans. Information Theory,
%   vol. 65, no. 12, pp. 8305-8324, 2019.
%
% See also: aoi_fcfs_mm1, aoi_fcfs_gim1, aoi_lst_exp, aoi_lst_erlang

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

% Check stability
if rho >= 1
    line_error(mfilename, 'System unstable: rho = lambda*E_H = %.4f >= 1', rho);
end

% Mean interarrival time
E_Y = 1 / lambda;

% Mean waiting time (Pollaczek-Khinchine formula)
% E[W] = lambda * E[H^2] / (2 * (1 - rho))
E_W = lambda * E_H2 / (2 * (1 - rho));

% Mean system time (sojourn time)
E_T = E_W + E_H;

% see _kb/03-api-layer.md (AoI family) for the Inoue et al. 2019 derivation
Tstar = @(s) H_lst(s) .* ((1 - rho) .* s ./ (s - lambda + lambda .* H_lst(s)));
hstep = 1e-6 * max(1, lambda);
dTstar = (Tstar(lambda + hstep) - Tstar(lambda - hstep)) / (2 * hstep);
meanAoI = E_H + E_T + (1 - 2 * rho) / lambda - dTstar;

% Mean Peak AoI
% E[Apeak] = E[T] + E[Y]
peakAoI = E_T + E_Y;

% LST of AoI (Inoue et al. 2019, Theorem 2), via the general age formula
%   A*(s) = (lambda/s) * ( T*(s) - Apeak*(s) )
% which is the cycle average of exp(-s*age) over a departure interval: the age
% starts each cycle at the system time T of the packet just delivered and grows
% linearly to the peak Apeak = T + (next interarrival) at the next delivery.
%
% For M/GI/1 FCFS the peak is max(Y, T) plus a fresh service, because Lindley
% gives W' = max(0, T - Y) and hence W' + Y = max(Y, T). With Y ~ Exp(lambda)
% independent of T,
%   E[exp(-s*max(Y,T))] = T*(s) - (s/(s+lambda)) * T*(s+lambda),
% so   Apeak*(s) = H*(s) * ( T*(s) - (s/(s+lambda)) * T*(s+lambda) ).
%
% THE PREVIOUS FORM WAS NOT AN LST. It read
%   (lambda*H*(s)) / (s + lambda - lambda*H*(s)) * W*(s),
% whose first factor diverges as s -> 0 (the denominator vanishes like
% s*(1+rho)), so A*(0) was +Inf instead of 1 and the value exceeded 1 for small
% s. Checked against simulation on M/E2/1: at s = 0.3 the old form gave 1.3062,
% the form below 0.56935, and the sample path 0.56948.
lstAoI = @(s) aoi_mgi1_lst_eval(s, lambda, H_lst, rho);

end

function T_s = aoi_mgi1_tstar(s, lambda, H_lst, rho)
%AOI_MGI1_TSTAR System-time LST T*(s) = H*(s)*W*(s), with T*(0) = 1 by continuity
    if abs(s) < 1e-12
        T_s = 1;
        return
    end
    H_s = H_lst(s);
    T_s = H_s * (1 - rho) * s / (s - lambda + lambda * H_s);
end

function val = aoi_mgi1_lst_eval(s, lambda, H_lst, rho)
%AOI_MGI1_LST_EVAL Evaluate M/GI/1 FCFS AoI LST at point(s) s

    val = zeros(size(s));
    for i = 1:numel(s)
        si = s(i);
        if abs(si) < 1e-12
            val(i) = 1; % A*(0) = 1 for any proper LST
            continue
        end
        H_s = H_lst(si);
        T_s = aoi_mgi1_tstar(si, lambda, H_lst, rho);
        T_sl = aoi_mgi1_tstar(si + lambda, lambda, H_lst, rho);
        peak_s = H_s * (T_s - (si / (si + lambda)) * T_sl);
        val(i) = (lambda / si) * (T_s - peak_s);
    end
    if isscalar(s)
        val = val(1);
    end
end
