function n = sn_join_siblings(sn, joinIdx, r)
% N = SN_JOIN_SIBLINGS(SN, JOINIDX, R)
%
% Number of sibling tasks forked per parent job of class R on the fork-join
% pair that ends at the Join node JOINIDX. Siblings are counted at the FORK,
% as the simulation engines count them.
%
% THE COUNT IS PER LINK, not the out-degree times a node-wide scalar. A fork
% carries a VARIABLE FORKING LEVEL: `setTasksPerLink(class, n [, dest])` sets
% one link of one class, `setTasksPerLinkDistribution` makes the degree a
% draw, and `setBranchProb` makes a link taken with probability < 1. All three
% land in `sn.nodeparam{f}.fanOutLink` and `.fanOutProb`, both
% (nnodes x nclasses) and indexed by DESTINATION NODE, with the DISTRIBUTION
% case storing its mean. So
%
%   N = sum_d fanOutProb(d,r) * fanOutLink(d,r)
%
% which is the EXPECTED sibling count, and reduces to out-degree times the
% node-wide `fanOut` on a model that sets none of the three (every connected
% link then carries fanOutLink = tasksPerLink and fanOutProb = 1).
%
% Falls back to that older product when `fanOutLink` is absent -- it is built
% by refreshLocalVars, which runs after refreshCapacity -- and to the join's
% in-degree when the matched fork cannot be identified.
%
% N is what a quorum is measured against: the join fires on the k-th of N
% siblings and the remaining N-k are discarded when they arrive.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 3 || isempty(r)
    r = [];
end

n = 0;
if joinIdx < 1 || ~isfield(sn,'connmatrix') || joinIdx > size(sn.connmatrix,1)
    return
end
n = nnz(sn.connmatrix(:,joinIdx));
if ~isfield(sn,'fj') || isempty(sn.fj) || joinIdx > size(sn.fj,2)
    return
end
f = find(sn.fj(:,joinIdx), 1);
if isempty(f)
    return
end

np = [];
if f <= length(sn.nodeparam) && isstruct(sn.nodeparam{f})
    np = sn.nodeparam{f};
end

% The per-link count, when the refresh has built it.
if ~isempty(np) && isfield(np,'fanOutLink') && ~isempty(np.fanOutLink)
    fol = np.fanOutLink;
    if isfield(np,'fanOutProb') && isequal(size(np.fanOutProb), size(fol))
        fop = np.fanOutProb;
    else
        fop = double(fol > 0);
    end
    if isempty(r)
        % No class asked for: the widest fork over the classes, which is the
        % count a class-blind caller must not under-state.
        cols = 1:size(fol,2);
    else
        cols = r(r >= 1 & r <= size(fol,2));
    end
    if ~isempty(cols)
        n = max(sum(fol(:,cols) .* fop(:,cols), 1));
        return
    end
end

% Fallback: out-degree times the node-wide scalar.
w = 1;
if ~isempty(np) && isfield(np,'fanOut') && ~isempty(np.fanOut)
    w = max(1, round(np.fanOut(1)));
end
n = nnz(sn.connmatrix(f,:)) * w;
end
