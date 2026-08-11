%% Cache RMF Transient Analysis
% Demonstrates refined mean field (RMF) transient and steady-state
% analysis of a multi-list cache with RANDOM(m) replacement.
%
% Uses CacheRMF directly to compute:
%   - Steady-state hit/miss probabilities with 1/N correction
%   - Transient evolution of hit rates via coupled ODE system

clear;

% Cache parameters
n = 10;         % number of items
m = [3, 2];     % list capacities (2-list cache)
alpha = 0.8;    % Zipf exponent

% Zipf popularity distribution
p = (1:n).^(-alpha);
p = p / sum(p);

fprintf('Cache parameters: n=%d, m=[%s], Zipf(%.1f)\n', n, strjoin(arrayfun(@num2str, m, 'UniformOutput', false), ', '), alpha);

% Build DDPP model
model = CacheRMF(p, m);

% Steady-state analysis with 1/N correction
[pi, V, VW] = model.meanFieldExpansionSteadyState(1);
pi_refined = pi(:)' + V(:)' / n;

fprintf('\nSteady-state results (refined mean field):\n');
total_hit = 0;
for k = 1:length(m)
    hr = model.hitRate(pi_refined, k);
    total_hit = total_hit + hr;
    fprintf('  Hit rate (list %d): %.6f\n', k, hr);
end
miss_rate = model.hitRate(pi_refined, 0);
fprintf('  Miss rate:         %.6f\n', miss_rate);
fprintf('  Total hit prob:    %.6f\n', total_hit);
fprintf('  Total miss prob:   %.6f\n', miss_rate);

% Transient analysis
[T, X, Vt, W] = model.meanFieldExpansionTransient(50, 200, 1);

fprintf('\nTransient hit rates (refined, N=%d):\n', n);
time_indices = [1, 21, 51, 101, 200];  % t=0, ~5, ~12.5, ~25, 50
for idx = time_indices
    xt = X(idx, :) + Vt(idx, :) / n;
    hr = 0;
    for k = 1:length(m)
        hr = hr + model.hitRate(xt, k);
    end
    fprintf('  t=%7.3f: hit_rate=%.6f\n', T(idx), hr);
end
