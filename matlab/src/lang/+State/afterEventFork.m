function [outspace, outrate, outprob, eventCache] = afterEventFork(sn, ind, event, class, isSimulation, eventCache, space_buf, space_srv, space_var, key) %#ok<INUSL>
% [OUTSPACE, OUTRATE, OUTPROB, EVENTCACHE] = AFTEREVENTFORK(SN, IND, EVENT, CLASS, ISSIMULATION, EVENTCACHE, SPACE_BUF, SPACE_SRV, SPACE_VAR, KEY)
%
% Handle afterEvent logic for stateful Fork nodes (FJ-augmented structs
% only, see ModelAdapter.fjtag). The fork state is a per-class count of
% parent jobs momentarily held before the fork firing. Arrivals are
% buffered here; the atomic multi-branch emission is not a DEP event but
% a fork firing synchronization (sn.fjsync) handled by State.afterFJEvent.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

outspace = [];
outrate = [];
outprob = 1;

switch event
    case EventType.ARV
        space_srv(:,class) = space_srv(:,class) + 1;
        outspace = [space_srv, space_var];
        outrate = -1*ones(size(outspace,1),1); % passive action, rate is unspecified
    case EventType.DEP
        % departures from a Fork are never generated as regular syncs
        % (see refreshSync); they occur only through sn.fjsync firings
end

if isSimulation
    if nargin>=6 && isobject(eventCache)
        eventCache(key) = {outprob, outspace, outrate};
    end
    if size(outspace,1) > 1
        tot_rate = sum(outrate);
        cum_rate = cumsum(outrate) / tot_rate;
        firing_ctr = 1 + max([0,find( rand > cum_rate' )]); % select action
        outspace = outspace(firing_ctr,:);
        outrate = sum(outrate);
        outprob = outprob(firing_ctr,:);
    end
end

end
