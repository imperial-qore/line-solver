%{ @file cache_is.m
 %  @brief Importance sampling for cache normalizing constant
 %
 %  @author LINE Development Team
%}

%{
 % @brief Computes cache normalizing constant via importance sampling
 %
 % @details
 % This function estimates the normalizing constant for cache models using
 % Monte Carlo importance sampling. It provides an alternative to exact
 % enumeration (cache_erec) and saddle-point approximation (cache_spm)
 % that scales better for large item counts while providing confidence bounds.
 %
 % @par Syntax:
 % @code
 % [E, lE] = cache_is(gamma, m)
 % [E, lE] = cache_is(gamma, m, samples)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>gamma<td>Item popularity probabilities (n x h matrix)
 % <tr><td>m<td>Cache capacity vector (1 x h)
 % <tr><td>samples<td>(Optional) Number of Monte Carlo samples (default: 1e5)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>E<td>Normalizing constant estimate
 % <tr><td>lE<td>Log of normalizing constant
 % </table>
%}
function [E, lE] = cache_is(gamma, m, samples)
% [E, LE] = CACHE_IS(GAMMA, M, SAMPLES)
%
% Importance sampling estimation of cache normalizing constant.
%
% The normalizing constant E is defined as the sum over all valid cache
% configurations of the product of gamma values for items in each level.
%
% Input:
%   gamma   - (n x h) item popularity probabilities at each cache level
%   m       - (1 x h) cache capacity vector
%   samples - (optional) number of Monte Carlo samples, default 1e5
%
% Output:
%   E  - normalizing constant estimate
%   lE - log of normalizing constant

if nargin < 3 || isempty(samples)
    samples = 1e5;
end

% Remove items with zero gamma (no contribution)
gamma = gamma(sum(gamma, 2) > 0, :);

[n, h] = size(gamma);
mt = sum(m);  % total cache capacity

% Edge cases
if n == 0 || mt == 0
    E = 1;
    lE = 0;
    return
end

if n < mt
    line_warning(mfilename, 'Number of items (%d) less than cache capacity (%d).', n, mt);
    E = 0;
    lE = -Inf;
    return
end

if n == mt
    % All items must be in cache - only one valid configuration
    E = cache_erec(gamma, m);
    lE = log(E);
    return
end

% Pre-compute log(gamma)
log_gamma = log(gamma + 1e-300);

% Pre-compute log(m(j)!)
log_m_fact = sum(factln(m));

% Log of number of ways to choose mt items from n
log_combinations = factln(n) - factln(mt) - factln(n - mt);

lZ_samples = zeros(samples, 1);

for s = 1:samples
    % Sample mt items uniformly without replacement
    selected = randperm(n, mt);

    % Assign items to levels uniformly at random
    assignment = assign_items_to_levels(m, selected);

    % Compute log of unnormalized state probability
    log_state_prob = log_m_fact;
    for j = 1:h
        items_in_level = assignment{j};
        for idx = 1:length(items_in_level)
            i = items_in_level(idx);
            log_state_prob = log_state_prob + log_gamma(i, j);
        end
    end

    % Proposal probability is 1/(C(n,mt) * multinomial(mt; m))
    % where multinomial(mt; m) = mt! / prod(m(j)!)
    log_multinomial = factln(mt) - log_m_fact;
    log_proposal = -log_combinations - log_multinomial;

    lZ_samples(s) = log_state_prob - log_proposal;
end

lE = logmeanexp(lZ_samples);
E = exp(lE);

end

%% Helper: Assign items to cache levels uniformly
function assignment = assign_items_to_levels(m, selected)
    h = length(m);
    assignment = cell(1, h);

    % Shuffle and assign
    perm = randperm(length(selected));
    shuffled = selected(perm);

    idx = 1;
    for j = 1:h
        assignment{j} = shuffled(idx:(idx + m(j) - 1));
        idx = idx + m(j);
    end
end
