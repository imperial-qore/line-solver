function space = spaceLocalVars(sn, ind, maxPending)
% SPACE = SPACELOCALVARS(QN, IND, MAXPENDING)
%
% MAXPENDING (optional, default 0) bounds the number of secondary (delayed-hit)
% requests that a cache node may merge onto a single in-flight fetch.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Generate state space for local state variables

%ind: node index
%ist = sn.nodeToStation(ind);
%isf = sn.nodeToStateful(ind);

if nargin < 3 || isempty(maxPending)
    maxPending = 0;
end

space = [];

switch sn.nodetype(ind)
    case NodeType.Cache
        rsCap = 0;
        if isfield(sn.nodeparam{ind}, 'retrievalSystemCapacity')
            rsCap = sn.nodeparam{ind}.retrievalSystemCapacity;
        end
        [~, rcItems] = State.cacheRetrievalClassMap(sn, ind);
        space = State.spaceCache(sn.nodeparam{ind}.nitems, sn.nodeparam{ind}.itemcap, rsCap, maxPending, rcItems);
end

for r=1:sn.nclasses
    switch sn.routing(ind,r)
        case RoutingStrategy.RROBIN
            % RR slot holds the destination node index — enumerate over outlinks.
            space = State.cartesian(space, sn.nodeparam{ind}{r}.outlinks(:));
        case RoutingStrategy.WRROBIN
            % WRR slot holds the POSITION in weighted_outlinks (1..len) —
            % enumerate positions; afterEventRouter advances the position
            % cyclically and sub_wrr maps it back to a destination.
            np = sn.nodeparam{ind}{r};
            if isfield(np, 'weighted_outlinks') && ~isempty(np.weighted_outlinks)
                positions = (1:length(np.weighted_outlinks))';
            else
                positions = (1:length(np.outlinks))';
            end
            space = State.cartesian(space, positions);
    end
end
end
