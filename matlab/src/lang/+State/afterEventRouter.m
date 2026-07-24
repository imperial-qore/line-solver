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
                    slot = sum(sn.nvars(ind,1:(sn.nclasses+class)));
                    outlinks = sn.nodeparam{ind}{class}.outlinks;
                    idx = find(space_var(slot) == outlinks);
                    if isempty(idx) || idx >= length(outlinks)
                        space_var(slot) = outlinks(1);
                    else
                        space_var(slot) = outlinks(idx+1);
                    end
                case RoutingStrategy.WRROBIN
                    slot = sum(sn.nvars(ind,1:(sn.nclasses+class)));
                    % WRR cycles through weighted_outlinks (each outlink
                    % replicated by its weight). The state slot holds the
                    % POSITION in this list rather than a destination value,
                    % so that repeated outlinks advance correctly.
                    if isfield(sn.nodeparam{ind}{class}, 'weighted_outlinks') ...
                            && ~isempty(sn.nodeparam{ind}{class}.weighted_outlinks)
                        cycle_len = length(sn.nodeparam{ind}{class}.weighted_outlinks);
                    else
                        cycle_len = length(sn.nodeparam{ind}{class}.outlinks);
                    end
                    pos = space_var(slot);
                    if pos < 1 || pos > cycle_len
                        space_var(slot) = 1;
                    elseif pos >= cycle_len
                        space_var(slot) = 1;
                    else
                        space_var(slot) = pos + 1;
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