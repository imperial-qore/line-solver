function [outspace, outrate, outprob, eventCache] = afterEventRouter(sn, ind, event, class, isSimulation, eventCache, space_buf, space_srv, space_var, key)
% [OUTSPACE, OUTRATE, OUTPROB, EVENTCACHE] = AFTEREVENTROUTER(SN, IND, EVENT, CLASS, ISSIMULATION, EVENTCACHE, SPACE_BUF, SPACE_SRV, SPACE_VAR, KEY)
%
% Handle router afterEvent logic

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

outspace = [];
outrate = [];
outprob = 1;

switch event
    case EventType.ARV
        space_srv(:,class) = space_srv(:,class) + 1;
        outspace = [space_srv, space_var]; % buf is empty
        outrate = -1*ones(size(outspace,1)); % passive action, rate is unspecified
    case EventType.DEP
        if space_srv(class)>0
            space_srv(:,class) = space_srv(:,class) - 1;
            switch sn.routing(ind,class)
                case RoutingStrategy.RROBIN
                    idx = find(space_var(sum(sn.nvars(ind,1:(sn.nclasses+class)))) == sn.nodeparam{ind}{class}.outlinks);
                    if idx < length(sn.nodeparam{ind}{class}.outlinks)
                        space_var(sum(sn.nvars(ind,1:(sn.nclasses+class)))) = sn.nodeparam{ind}{class}.outlinks(idx+1);
                    else
                        space_var(sum(sn.nvars(ind,1:(sn.nclasses+class)))) = sn.nodeparam{ind}{class}.outlinks(1);
                    end
            end
            outspace = [space_srv, space_var]; % buf is empty
            outrate = GlobalConstants.Immediate*ones(size(outspace,1)); % immediate action
        end
end

if isSimulation
    if nargin>=8 && isobject(eventCache)
        eventCache(key) = {outprob, outspace,outrate};
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