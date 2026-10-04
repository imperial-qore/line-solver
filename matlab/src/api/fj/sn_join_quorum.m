function k = sn_join_quorum(sn, joinIdx, r, nbranches)
% K = SN_JOIN_QUORUM(SN, JOINIDX, R, NBRANCHES)
%
% Number of sibling tasks the Join node JOINIDX waits for in class R, given
% that NBRANCHES of them are forked. A standard join, an absent declaration,
% a non-positive quorum and a quorum that is not smaller than the sibling
% count all return NBRANCHES, i.e. the ordinary AND-join: those are the four
% ways a join fires only when every sibling has arrived.
%
% The count is the one the simulation engines apply (SolverLDES fixes it at
% FORK time and discards the stragglers when they reach the join), so an
% analytical solver reading it here charges the same synchronisation event.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

k = nbranches;
if joinIdx < 1 || joinIdx > length(sn.nodeparam) || ~isstruct(sn.nodeparam{joinIdx})
    return
end
np = sn.nodeparam{joinIdx};
if ~isfield(np,'joinStrategy') || length(np.joinStrategy) < r || isempty(np.joinStrategy{r})
    return
end
if np.joinStrategy{r} == JoinStrategy.STD
    return
end
if ~isfield(np,'joinRequired') || length(np.joinRequired) < r || isempty(np.joinRequired{r})
    return
end
q = round(np.joinRequired{r});
if q > 0 && q < nbranches
    k = q;
end
end
