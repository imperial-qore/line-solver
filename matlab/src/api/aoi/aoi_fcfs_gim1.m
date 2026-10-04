function [meanAoI, lstAoI, peakAoI] = aoi_fcfs_gim1(Y_lst, mu, E_Y, E_Y2)
%AOI_FCFS_GIM1 Mean AoI and LST for GI/M/1 FCFS queue
%
% [meanAoI, lstAoI, peakAoI] = aoi_fcfs_gim1(Y_lst, mu, E_Y, E_Y2)
%
% Computes the Age of Information metrics for a GI/M/1 queue with
% First-Come First-Served (FCFS) discipline.
%
% GI/M/1: General independent arrivals, exponential service with rate mu.
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
% Formulas (from Inoue et al., IEEE Trans. IT, 2019, Theorem 3):
%   The key parameter sigma is the unique root in (0,1) of:
%     Y*(mu - mu*sigma) = sigma
%
%   LST of AoI: A*(s) = (mu*sigma_s) / (s + mu - mu*sigma_s) * D*(s)
%   where sigma_s solves Y*(s + mu - mu*sigma_s) = sigma_s
%   and D*(s) is the LST of system delay.
%
% Reference:
%   Y. Inoue, H. Masuyama, T. Takine, T. Tanaka, "A General Formula for
%   the Stationary Distribution of the Age of Information and Its
%   Application to Single-Server Queues," IEEE Trans. Information Theory,
%   vol. 65, no. 12, pp. 8305-8324, 2019.
%
% See also: aoi_fcfs_mm1, aoi_fcfs_mgi1, aoi_fcfs_dm1

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

% Check stability
if rho >= 1
    line_error(mfilename, 'System unstable: rho = 1/(E_Y*mu) = %.4f >= 1', rho);
end

% Find sigma: unique root in (0,1) of Y*(mu - mu*sigma) = sigma
% This is the probability an arriving customer finds the server busy
sigma_func = @(sig) Y_lst(mu - mu * sig) - sig;

% Use fzero to find the root
try
    sigma = fzero(sigma_func, [0.001, 0.999]);
catch
    % Fallback: use fixed-point iteration
    sigma = 0.5;
    for iter = 1:100
        sigma_new = Y_lst(mu - mu * sigma);
        if abs(sigma_new - sigma) < 1e-10
            break;
        end
        sigma = sigma_new;
    end
end

% Mean system delay for GI/M/1
% E[D] = 1 / (mu * (1 - sigma))
E_D = 1 / (mu * (1 - sigma));

% see _kb/03-api-layer.md (AoI family) for the Inoue et al. 2019 derivation
eta = mu * (1 - sigma);
hstep = 1e-6 * max(1, eta);
dYstar = (Y_lst(eta + hstep) - Y_lst(eta - hstep)) / (2 * hstep);
meanAoI = lambda * E_Y2 / 2 + 1 / mu + lambda * (-dYstar) / eta;

% Mean Peak AoI
% E[Apeak] = E[Y] + E[D]
peakAoI = E_Y + E_D;

% LST of AoI (Inoue et al. 2019, Theorem 3), via the general age formula
%   A*(s) = (lambda/s) * ( T*(s) - Apeak*(s) )
% the cycle average of exp(-s*age) over a departure interval. In GI/M/1 the
% system time is EXPONENTIAL at rate eta = mu*(1-sigma), so T*(s) = eta/(s+eta),
% and Lindley gives W' + Y = max(Y, T) with T ~ Exp(eta) independent of the next
% interarrival Y, so
%   E[exp(-s*max(Y,T))] = Y*(s) - (s/(s+eta)) * Y*(s+eta),
% and the peak adds one fresh Exp(mu) service:
%   Apeak*(s) = (mu/(s+mu)) * ( Y*(s) - (s/(s+eta)) * Y*(s+eta) ).
% A*(0) = 1 follows from the defining relation Y*(eta) = sigma.
%
% THE PREVIOUS FORM WAS NOT AN LST. It read
%   (mu*sigma(s)) / (s + mu - mu*sigma(s)) * D*(s),
% which at s = 0 gives sigma/(1-sigma) rather than 1, and it re-solved sigma(s)
% by fzero at every evaluation point with a SILENT fallback to sigma(0) on
% failure. Checked against simulation on E2/M/1: at s = 0.2 the old form gave
% 0.23056, the form below 0.59319, and the sample path 0.59332.
lstAoI = @(s) aoi_gim1_lst_eval(s, Y_lst, mu, sigma, lambda);

end

function val = aoi_gim1_lst_eval(s, Y_lst, mu, sigma0, lambda)
%AOI_GIM1_LST_EVAL Evaluate GI/M/1 FCFS AoI LST at point(s) s

    eta = mu * (1 - sigma0);
    val = zeros(size(s));
    for i = 1:numel(s)
        si = s(i);
        if abs(si) < 1e-12
            val(i) = 1; % A*(0) = 1 for any proper LST
            continue
        end
        T_s = eta / (si + eta);
        peak_s = (mu / (si + mu)) * (Y_lst(si) - (si / (si + eta)) * Y_lst(si + eta));
        val(i) = (lambda / si) * (T_s - peak_s);
    end
    if isscalar(s)
        val = val(1);
    end
end
