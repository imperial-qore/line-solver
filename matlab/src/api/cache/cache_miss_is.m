%{ @file cache_miss_is.m
 %  @brief Computes miss rates using importance sampling
 %
 %  @author LINE Development Team
%}

%{
 % @brief Computes cache miss rates using importance sampling
 %
 % @details
 % This function computes global, per-user, and per-item miss rates for
 % cache models using Monte Carlo importance sampling.
 %
 % @par Syntax:
 % @code
 % [M, MU, MI, pi0, lE] = cache_miss_is(gamma, m, lambda)
 % [M, MU, MI, pi0, lE] = cache_miss_is(gamma, m, lambda, samples)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>gamma<td>Item popularity probabilities (n x h matrix)
 % <tr><td>m<td>Cache capacity vector (1 x h)
 % <tr><td>lambda<td>Arrival rates per user per item (u x n x h+1)
 % <tr><td>samples<td>(Optional) Number of Monte Carlo samples (default: 1e5)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>M<td>Global miss rate
 % <tr><td>MU<td>Per-user miss rate (u x 1)
 % <tr><td>MI<td>Per-item miss rate (n x 1)
 % <tr><td>pi0<td>Per-item miss probability (1 x n)
 % <tr><td>lE<td>Log of normalizing constant
 % </table>
%}
function [M, MU, MI, pi0, lE] = cache_miss_is(gamma, m, lambda, samples)
% [M, MU, MI, PI0, LE] = CACHE_MISS_IS(GAMMA, M, LAMBDA, SAMPLES)
%
% Importance sampling estimation of cache miss rates.
%
% Input:
%   gamma   - (n x h) item popularity probabilities
%   m       - (1 x h) cache capacity vector
%   lambda  - (u x n x h+1) arrival rates per user per item per level
%   samples - (optional) number of Monte Carlo samples, default 1e5
%
% Output:
%   M   - global miss rate
%   MU  - per-user miss rate
%   MI  - per-item miss rate
%   pi0 - per-item miss probability
%   lE  - log of normalizing constant

if nargin < 4 || isempty(samples)
    samples = 1e5;
end

[n, h] = size(gamma);

% Compute normalizing constant via importance sampling
[~, lE] = cache_is(gamma, m, samples);

% Compute hit probabilities via importance sampling
pij = cache_prob_is(gamma, m, samples);

% Extract miss probabilities (first column)
pi0 = pij(:, 1)';

% Compute miss rates if lambda provided
if nargin >= 3 && ~isempty(lambda)
    u = size(lambda, 1);  % number of users

    % Per-user miss rate
    MU = zeros(u, 1);
    for v = 1:u
        for k = 1:n
            MU(v) = MU(v) + lambda(v, k, 1) * pi0(k);
        end
    end

    % Per-item miss rate
    MI = zeros(n, 1);
    for k = 1:n
        MI(k) = sum(lambda(:, k, 1)) * pi0(k);
    end

    % Global miss rate (sum of per-item miss rates)
    M = sum(MI);
else
    % Without lambda, return just the mean miss probability
    M = mean(pi0);
    MU = [];
    MI = [];
end

end
