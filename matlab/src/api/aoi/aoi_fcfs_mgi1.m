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
%   E[A] = E[Y] + E[T] + lambda * E[H^2] / (2*(1-rho))
%        = 1/lambda + E[W] + E[H] + lambda * E[H^2] / (2*(1-rho))
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

% Mean AoI (Proposition 1)
% E[A] = E[Y] + E[T] + lambda * E[H^2] / (2*(1-rho))
meanAoI = E_Y + E_T + lambda * E_H2 / (2 * (1 - rho));

% Mean Peak AoI
% E[Apeak] = E[T] + E[Y]
peakAoI = E_T + E_Y;

% LST of AoI (Theorem 2)
% A*(s) = (lambda * H*(s)) / (s + lambda - lambda*H*(s)) * W*(s)
% where W*(s) is the Pollaczek-Khinchine LST:
% W*(s) = (1-rho) * s / (s - lambda + lambda*H*(s))

lstAoI = @(s) aoi_mgi1_lst_eval(s, lambda, H_lst, rho);

end

function val = aoi_mgi1_lst_eval(s, lambda, H_lst, rho)
%AOI_MGI1_LST_EVAL Evaluate M/GI/1 FCFS AoI LST at point(s) s

    if isscalar(s)
        H_s = H_lst(s);
        % Pollaczek-Khinchine LST for waiting time
        W_s = (1 - rho) * s / (s - lambda + lambda * H_s);
        % AoI LST (Theorem 2)
        val = (lambda * H_s) / (s + lambda - lambda * H_s) * W_s;
    else
        % Handle array input
        val = zeros(size(s));
        for i = 1:numel(s)
            H_s = H_lst(s(i));
            W_s = (1 - rho) * s(i) / (s(i) - lambda + lambda * H_s);
            val(i) = (lambda * H_s) / (s(i) + lambda - lambda * H_s) * W_s;
        end
    end
end
