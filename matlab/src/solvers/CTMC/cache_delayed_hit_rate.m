function delayedRate = cache_delayed_hit_rate(sn, isf, InfGen, StateSpace, probSysState)
% CACHE_DELAYED_HIT_RATE Exact delayed-hit rate of a retrieval-aware cache.
%
% @brief Delayed-hit completion rate per originating class, as a transition
% reward over a solved CTMC.
%
% A fetch of item i completes on exactly the transitions that clear block A bit
% i, and each such transition releases the block-B counts of item i as delayed
% hits. The rate is therefore a TRANSITION reward over the generator, not a
% state reward: the alternative arrival-rate identity lambda_i*phi_i is only
% PASTA-exact. Shared by the CTMC analyzer and by the environment blend, which
% solves the same generators one stage at a time.
%
% @param sn NetworkStruct whose state space STATESPACE enumerates.
% @param isf Stateful index of the cache node.
% @param InfGen Infinitesimal generator of the stage, rows matching STATESPACE.
% @param StateSpace Enumerated state space, one row per state.
% @param probSysState Stationary probability of each row of STATESPACE.
% @return delayedRate (1 x nclasses) delayed-hit rate per originating class;
%         all zeros when the cache has no retrieval sub-system.

delayedRate = zeros(1, sn.nclasses);
ind = sn.statefulToNode(isf);
if sn.nodetype(ind) ~= NodeType.Cache
    return
end
np = sn.nodeparam{ind};
if ~isfield(np,'retrievalSystemCapacity') || np.retrievalSystemCapacity <= 0
    return
end
[~, rcItems, rcOrigClass] = State.cacheRetrievalClassMap(sn, ind);
if isempty(rcItems)
    return
end
nitems = np.nitems;
tcc = np.totalCacheCapacity;

% Column span of this stateful node's local variables inside a system state.
colOff = 0;
for jsf = 1:isf-1
    colOff = colOff + size(sn.space{jsf},2);
end
lvw = size(sn.space{isf},2);
cols = (colOff+1):(colOff+lvw);
lvs = lvw - (tcc + nitems + numel(rcItems)); % per-class server presence width
if lvs < 0
    return
end
aCols = cols(lvs+tcc+(1:nitems));
bCols = cols(lvs+tcc+nitems+(1:numel(rcItems)));

offdiag = InfGen - diag(diag(InfGen));
for j = 1:numel(rcItems)
    i = rcItems(j);
    rows = find(StateSpace(:,aCols(i)) ~= 0 & StateSpace(:,bCols(j)) > 0);
    for rr = rows(:).'
        nz = find(offdiag(rr,:) ~= 0);
        completes = nz(StateSpace(nz, aCols(i)) == 0);
        if isempty(completes), continue; end
        delayedRate(rcOrigClass(j)) = delayedRate(rcOrigClass(j)) ...
            + probSysState(rr) * StateSpace(rr,bCols(j)) * sum(offdiag(rr,completes));
    end
end
end
