%{ @file cache_prob_is.m
 %  @brief Computes cache hit probabilities using importance sampling
 %
 %  @author LINE Development Team
%}

%{
 % @brief Computes cache hit probabilities via importance sampling
 %
 % @details
 % This function estimates cache hit probability distribution using
 % Monte Carlo importance sampling. For each item i, it estimates the
 % probability that item i is cached at level j (hit) or not cached (miss).
 %
 % @par Syntax:
 % @code
 % prob = cache_prob_is(gamma, m)
 % prob = cache_prob_is(gamma, m, samples)
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
 % <tr><td>prob<td>Cache hit probability distribution (n x h+1 matrix)
 %         prob(i,1) = miss probability for item i
 %         prob(i,1+j) = hit probability for item i at level j
 % </table>
%}
function prob = cache_prob_is(gamma, m, samples)
% PROB = CACHE_PROB_IS(GAMMA, M, SAMPLES)
%
% Importance sampling estimation of cache hit probabilities.
%
% Input:
%   gamma   - (n x h) item popularity probabilities at each cache level
%   m       - (1 x h) cache capacity vector
%   samples - (optional) number of Monte Carlo samples, default 1e5
%
% Output:
%   prob    - (n x h+1) hit probability matrix
%             prob(i,1) = miss probability for item i
%             prob(i,1+j) = hit probability for item i at cache level j

if nargin < 3 || isempty(samples)
    samples = 1e5;
end

[n, h] = size(gamma);
mt = sum(m);

% Initialize probability matrix
prob = zeros(n, h + 1);

% Edge cases
if n == 0 || mt == 0
    prob(:, 1) = 1;  % all items miss
    return
end

if n < mt
    line_warning(mfilename, 'Number of items (%d) less than cache capacity (%d).', n, mt);
    prob(:, 1) = 1;
    return
end

if n == mt
    % All items must be in cache - use exact method
    prob = cache_prob_erec(gamma, m);
    return
end

log_gamma = log(gamma + 1e-300);
log_m_fact = sum(factln(m));
log_combinations = factln(n) - factln(mt) - factln(n - mt);
log_multinomial = factln(mt) - log_m_fact;

% Accumulators for importance sampling estimates
item_level_weight = zeros(n, h);  % sum of weights when item i is at level j
total_weight = 0;  % sum of all weights

for s = 1:samples
    % Sample mt items uniformly without replacement
    selected = randperm(n, mt);

    % Assign items to levels
    assignment = assign_items_to_levels(m, selected);

    % Compute unnormalized state probability
    log_state_prob = log_m_fact;
    for j = 1:h
        items_in_level = assignment{j};
        for idx = 1:length(items_in_level)
            i = items_in_level(idx);
            log_state_prob = log_state_prob + log_gamma(i, j);
        end
    end

    % Compute proposal probability
    log_proposal = -log_combinations - log_multinomial;

    % Importance weight for this sample
    log_is_weight = log_state_prob - log_proposal;
    is_weight = exp(log_is_weight - 50);  % scale to avoid overflow

    % Update accumulators
    total_weight = total_weight + is_weight;
    for j = 1:h
        items_in_level = assignment{j};
        for idx = 1:length(items_in_level)
            i = items_in_level(idx);
            item_level_weight(i, j) = item_level_weight(i, j) + is_weight;
        end
    end
end

% Compute probabilities from accumulated weights
if total_weight > 0
    for i = 1:n
        for j = 1:h
            prob(i, 1 + j) = item_level_weight(i, j) / total_weight;
        end
        prob(i, 1) = max(0, 1 - sum(prob(i, 2:end)));
    end
else
    % Fallback: all items miss
    prob(:, 1) = 1;
end

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
