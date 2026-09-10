function [rcList, rcItems, rcOrigClass] = cacheRetrievalClassMap(sn, ind)
% [RCLIST, RCITEMS, RCORIGCLASS] = CACHERETRIEVALCLASSMAP(SN, IND)
%
% Canonical ordering of the retrieval classes of cache node IND, used to index
% block B of the cache local-variable vector (the per-retrieval-class counts of
% merged secondary requests).
%
%   rcList       retrieval class indices, ascending
%   rcItems      item served by each entry of rcList
%   rcOrigClass  originating (arrival) class of each entry of rcList
%
% Returns empty arrays when the node has no retrieval system.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

rcList = [];
rcItems = [];
rcOrigClass = [];
np = sn.nodeparam{ind};
if ~isfield(np, 'retrievalClasses') || isempty(np.retrievalClasses)
    return
end
rc = np.retrievalClasses; % (nitems x nclasses), entry = retrieval class index or -1
[nitems, nclasses] = size(rc);
for k = 1:nitems
    for c = 1:nclasses
        if rc(k,c) > 0
            rcList(end+1) = rc(k,c); %#ok<AGROW>
            rcItems(end+1) = k; %#ok<AGROW>
            rcOrigClass(end+1) = c; %#ok<AGROW>
        end
    end
end
[rcList, order] = sort(rcList);
rcItems = rcItems(order);
rcOrigClass = rcOrigClass(order);

end
