function [all_jumps_red, rateBase_red, eventIdx_red, state_map] = ...
    ode_eliminate_immediate(all_jumps, rateBase, eventIdx, sn, options)
% [ALL_JUMPS_RED, RATEBASE_RED, EVENTIDX_RED, STATE_MAP] = ...
%     ODE_ELIMINATE_IMMEDIATE(ALL_JUMPS, RATEBASE, EVENTIDX, SN, OPTIONS)
%
% Eliminate immediate transitions from Fluid ODE system using stochastic
% complementation to reduce stiffness
%
% Input:
%   all_jumps: [n_states x n_transitions] matrix of state changes
%   rateBase: [n_transitions x 1] vector of base rates
%   eventIdx: [n_transitions x 1] vector of state indices for rate gating
%   sn: Model structure
%   options: Solver options
%
% Output:
%   all_jumps_red: Reduced jump matrix
%   rateBase_red: Reduced rate vector
%   eventIdx_red: Reduced event index vector
%   state_map: Mapping from reduced to original state indices
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Get immediate detection threshold using GlobalConstants.Immediate with tolerance
% A transition is considered immediate if its rate is within 1% of GlobalConstants.Immediate
imm_ref = GlobalConstants.Immediate;
rel_tol = 0.01;
imm_tol = imm_ref * (1 - rel_tol);

% Identify immediate transitions (rates close to or above GlobalConstants.Immediate)
imm_idx = find(rateBase >= imm_tol);

% If no immediate transitions, return unchanged
if isempty(imm_idx)
    all_jumps_red = all_jumps;
    rateBase_red = rateBase;
    eventIdx_red = eventIdx;
    state_map = 1:size(all_jumps, 1);
    return;
end

% Report detected immediate transitions
if isfield(options, 'verbose') && options.verbose > 0
    line_printf(sprintf('Eliminating %d immediate transitions (out of %d total)\n', ...
        length(imm_idx), length(rateBase)));
end

try
    % Apply stochastic complement elimination
    [all_jumps_red, rateBase_red, eventIdx_red, state_map] = ...
        eliminate_via_stochcomp(all_jumps, rateBase, eventIdx, imm_idx, options);

catch ME
    % Elimination failed - fall back to original system
    warning('LINE:Fluid:EliminationFailed', ...
        'Immediate transition elimination failed: %s. Using original system.', ME.message);
    all_jumps_red = all_jumps;
    rateBase_red = rateBase;
    eventIdx_red = eventIdx;
    state_map = 1:size(all_jumps, 1);
end

end % ode_eliminate_immediate


function [all_jumps_red, rateBase_red, eventIdx_red, state_map] = ...
    eliminate_via_stochcomp(all_jumps, rateBase, eventIdx, imm_idx, options)
% Eliminate immediate transitions using stochastic complementation

n_states = size(all_jumps, 1);
n_transitions = length(rateBase);

% Step 1: Construct infinitesimal generator matrix from jumps and rates
W = sparse(n_states, n_states);

for t = 1:n_transitions
    src_state = eventIdx(t);
    jump = all_jumps(:, t);

    % Find destination state (where jump > 0)
    dst_states = find(jump > 0);

    if ~isempty(dst_states)
        for d = 1:length(dst_states)
            dst_state = dst_states(d);
            % Add transition rate from src to dst
            W(src_state, dst_state) = W(src_state, dst_state) + rateBase(t) * jump(dst_state);
        end
    end
end

% Make W a proper generator (row sums = 0)
W = W - diag(sum(W, 2));

% Step 2: Identify immediate vs timed states
% Immediate states are those that have outgoing immediate transitions
imm_states = unique(eventIdx(imm_idx));
timed_states = setdiff(1:n_states, imm_states);

% Check that we have some timed states remaining
if isempty(timed_states)
    error('LINE:Fluid:AllImmediate', ...
        'All states have immediate transitions - cannot eliminate');
end

% Step 3: Apply stochastic complement
[W_red, ~, ~, ~, Q22, ~] = ctmc_stochcomp(W, timed_states);

% Check conditioning of Q22
if rcond(full(Q22)) < 1e-12
    warning('LINE:Fluid:IllConditioned', ...
        'Immediate transition matrix ill-conditioned (rcond = %g)', rcond(full(Q22)));
end

% Step 4: Convert reduced generator back to jump/rate representation
[all_jumps_red_local, rateBase_red, eventIdx_red_local] = generator_to_jumps(W_red);

% Step 5: Set up state mapping and expand to original state space
state_map = timed_states;

% Expand jump matrix and eventIdx to original state space dimension
% This ensures compatibility with ode_rates_closing which uses original q_indices
all_jumps_red = zeros(n_states, length(rateBase_red));
for t = 1:length(rateBase_red)
    jump_local = all_jumps_red_local(:, t);
    % Map local indices to original indices
    for s = 1:length(state_map)
        if jump_local(s) ~= 0
            all_jumps_red(state_map(s), t) = jump_local(s);
        end
    end
end

% Map eventIdx from reduced to original state space
eventIdx_red = state_map(eventIdx_red_local);

% Report reduction
if isfield(options, 'verbose') && options.verbose > 0
    reduction_pct = 100 * (1 - length(timed_states)/n_states);
    line_printf(sprintf('State space reduced from %d to %d states (%.1f%% reduction)\n', ...
        n_states, length(timed_states), reduction_pct));
end

end % eliminate_via_stochcomp


function x_full = ode_expand_state(x_reduced, state_map, n_original)
% X_FULL = ODE_EXPAND_STATE(X_REDUCED, STATE_MAP, N_ORIGINAL)
%
% Expand reduced state vector back to original dimension
%
% Input:
%   x_reduced: [n_reduced x 1] reduced state vector
%   state_map: [n_reduced x 1] mapping from reduced to original indices
%   n_original: original state space dimension
%
% Output:
%   x_full: [n_original x 1] full state vector
%           Immediate states (not in state_map) are set to 0

x_full = zeros(n_original, 1);
x_full(state_map) = x_reduced;

end % ode_expand_state
