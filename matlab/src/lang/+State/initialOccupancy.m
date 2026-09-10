function n = initialOccupancy(sn, ind, r)
% N = INITIALOCCUPANCY(SN, IND, R)
%
% Number of class-R jobs held by node IND in the DECLARED initial state, or 0
% when no initial state is available.
%
% Used to bound the enumerated local state space from below. A station whose
% visit ratio is zero is never re-entered, but it may still hold jobs at time
% zero: an SPN place with no input arc is a TRANSIENT state of the chain, not an
% absent one. Enumerating up to the initial occupancy keeps that state in the
% space; the unreachable-state pruning in solver_ctmc then removes whatever the
% chain cannot actually reach.
%
% See _kb/11-conventions-and-gotchas.md.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

n = 0;
if ~isfield(sn,'state') || isempty(sn.state)
    return
end
% only a Place holds tokens in a class-indexed row; every other node type encodes
% its local state differently, so column r is not an occupancy there
if sn.nodetype(ind) ~= NodeType.Place
    return
end
isf = sn.nodeToStateful(ind);
if isf < 1 || isf > numel(sn.state) || isempty(sn.state{isf})
    return
end
state_i = sn.state{isf}(1,:);
try
    [~, nir] = State.toMarginal(sn, ind, state_i);
catch
    % the marginal decoder needs a state of the enumerated width; before the
    % space exists the declared row can be narrower, in which case the raw
    % per-class counts are the best available reading
    if numel(state_i) >= r
        n = state_i(r);
    end
    return
end
if ~isempty(nir) && numel(nir) >= r
    n = nir(r);
end
end
