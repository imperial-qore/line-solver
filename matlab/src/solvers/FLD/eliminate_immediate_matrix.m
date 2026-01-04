function [W_red, state_map] = eliminate_immediate_matrix(W, sn, options)
% [W_RED, STATE_MAP] = ELIMINATE_IMMEDIATE_MATRIX(W, SN, OPTIONS)
%
% Eliminate immediate states from infinitesimal generator matrix W
%
% Input:
%   W: [n_states x n_states] infinitesimal generator matrix
%   sn: Model structure
%   options: Solver options
%
% Output:
%   W_red: Reduced generator matrix
%   state_map: Mapping from reduced to original state indices
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Get immediate detection threshold
if isfield(options.config, 'immediate_tol')
    imm_tol = options.config.immediate_tol;
else
    imm_tol = GlobalConstants.Immediate / 10; % Default: 1e7
end

% Identify immediate states (those with very high outgoing rates)
max_rate_per_state = max(abs(W), [], 2);
imm_states = find(max_rate_per_state >= imm_tol);

% If no immediate states, return unchanged
if isempty(imm_states)
    W_red = W;
    state_map = 1:size(W, 1);
    return;
end

% Identify timed states
timed_states = setdiff(1:size(W,1), imm_states);

% Check that we have some timed states remaining
if isempty(timed_states)
    error('LINE:Fluid:AllImmediate', ...
        'All states have immediate transitions - cannot eliminate');
end

try
    % Apply stochastic complement
    [W_red, ~, ~, ~, Q22, ~] = ctmc_stochcomp(W, timed_states);

    % Check conditioning
    if rcond(full(Q22)) < 1e-12
        warning('LINE:Fluid:IllConditioned', ...
            'Immediate transition matrix ill-conditioned (rcond = %g)', rcond(full(Q22)));
    end

    state_map = timed_states;

    % Report reduction
    if isfield(options, 'verbose') && options.verbose > 0
        reduction_pct = 100 * (1 - length(timed_states)/size(W,1));
        line_printf(sprintf('State space reduced from %d to %d states (%.1f%% reduction)\n', ...
            size(W,1), length(timed_states), reduction_pct));
    end

catch ME
    % Elimination failed - fall back to original
    warning('LINE:Fluid:EliminationFailed', ...
        'Immediate transition elimination failed: %s. Using original system.', ME.message);
    W_red = W;
    state_map = 1:size(W, 1);
end

end % eliminate_immediate_matrix
