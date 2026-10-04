function [all_jumps, rateBase, eventIdx] = generator_to_jumps(W)
% [ALL_JUMPS, RATEBASE, EVENTIDX] = GENERATOR_TO_JUMPS(W)
%
% Convert infinitesimal generator matrix to jump/rate representation
%
% Input:
%   W: [n_states x n_states] infinitesimal generator matrix
%      W(i,j) = rate of transition from state i to state j (i~=j)
%      W(i,i) = -sum of outgoing rates from state i
%
% Output:
%   all_jumps: [n_states x n_transitions] matrix of state change vectors
%   rateBase: [n_transitions x 1] vector of transition rates
%   eventIdx: [n_transitions x 1] vector of source state indices
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

n_states = size(W, 1);

% Find all non-zero off-diagonal elements
[src_states, dst_states, rates] = find(W);

% Filter out diagonal elements and zero/negative rates
off_diag_mask = (src_states ~= dst_states) & (rates > 0);
src_states = src_states(off_diag_mask);
dst_states = dst_states(off_diag_mask);
rates = rates(off_diag_mask);

n_transitions = length(rates);

% Pre-allocate outputs
all_jumps = zeros(n_states, n_transitions);
rateBase = rates;
eventIdx = src_states;

% Construct jump vectors
for t = 1:n_transitions
    src = src_states(t);
    dst = dst_states(t);

    % Jump: decrease at source, increase at destination
    all_jumps(src, t) = -1;
    all_jumps(dst, t) = +1;
end

end % generator_to_jumps
