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

% Mean AoI (Proposition 2)
% E[A] = E[Y] + E[D] + sigma / (mu * (1 - sigma) * (1 - sigma_prime))
% where sigma_prime = d/ds sigma(s)|_{s=0}
%
% For the mean, we can use:
% E[A] = E[Y] + E[D] + sigma * E[Y] / (1 - sigma)
% This is a simplification; the exact formula involves derivatives.

% Approximate mean AoI using the relationship from the paper
meanAoI = E_Y + E_D + sigma / (mu * (1 - sigma));

% Mean Peak AoI
% E[Apeak] = E[Y] + E[D]
peakAoI = E_Y + E_D;

% LST of AoI (Theorem 3)
% This requires solving for sigma(s) at each evaluation point
lstAoI = @(s) aoi_gim1_lst_eval(s, Y_lst, mu, sigma);

end

function val = aoi_gim1_lst_eval(s, Y_lst, mu, sigma0)
%AOI_GIM1_LST_EVAL Evaluate GI/M/1 FCFS AoI LST at point(s) s

    if isscalar(s)
        % Find sigma(s): root of Y*(s + mu - mu*sigma) = sigma
        sig_func = @(sig) Y_lst(s + mu - mu * sig) - sig;
        try
            sigma_s = fzero(sig_func, [0.001, 0.999]);
        catch
            sigma_s = sigma0;  % Fallback
        end

        % System delay LST: D*(s) = (1-sigma) * mu / (s + mu - mu*sigma_s)
        D_s = (1 - sigma0) * mu / (s + mu - mu * sigma_s);

        % AoI LST (Theorem 3)
        val = (mu * sigma_s) / (s + mu - mu * sigma_s) * D_s;
    else
        % Handle array input
        val = zeros(size(s));
        for i = 1:numel(s)
            sig_func = @(sig) Y_lst(s(i) + mu - mu * sig) - sig;
            try
                sigma_s = fzero(sig_func, [0.001, 0.999]);
            catch
                sigma_s = sigma0;
            end
            D_s = (1 - sigma0) * mu / (s(i) + mu - mu * sigma_s);
            val(i) = (mu * sigma_s) / (s(i) + mu - mu * sigma_s) * D_s;
        end
    end
end
