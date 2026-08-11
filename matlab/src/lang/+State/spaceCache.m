function SS = spaceCache(n, m, retrievalSystemCapacity)
% SS = SPACECACHE(N, M, RETRIEVALSYSTEMCAPACITY)
%
% Generate the local-variable state space of a cache node.
%   n  : total number of distinct items
%   m  : per-level cache capacity (row vector)
%   retrievalSystemCapacity (optional, default 0): number of extra slots tracking
%     items currently in the retrieval system. When > 0 each emitted row has
%     totalCacheCapacity + retrievalSystemCapacity columns; cache slots hold
%     1..n item indices (or 0 when partially filled at low item count k) and
%     retrieval slots hold 0 (empty) or 1..n (item being fetched).

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 3 || isempty(retrievalSystemCapacity)
    retrievalSystemCapacity = 0;
end

totalCacheCapacity = sum(m);
nItems = sum(n);
% The retrieval system is encoded as a per-item occupancy bitmap (one column per
% item) appended after the cache contents; bit i is set iff item i is currently
% being retrieved. With no retrieval system the bitmap is omitted entirely.
if retrievalSystemCapacity > 0
    retrievalWidth = nItems;
else
    retrievalWidth = 0;
end
nVars = totalCacheCapacity + retrievalWidth;

SS = zeros(0, nVars);
try
    if totalCacheCapacity == 0
        % degenerate cache with zero capacity
        SS = zeros(1, max(nVars,1));
    else
        items = 1:nItems;
        % Cache contents: every ordered placement of totalCacheCapacity distinct items.
        cacheCombos = nchoosek(items, totalCacheCapacity);
        for ci = 1:size(cacheCombos,1)
            cacheCombo = cacheCombos(ci,:);
            cachePerms = perms(cacheCombo);
            remaining = setdiff(items, cacheCombo);
            % Retrieval-system occupancy: every subset of the remaining items of
            % size <= retrievalSystemCapacity, as a one-hot bitmap.
            for s = 0:retrievalSystemCapacity
                if s == 0
                    retrCombos = zeros(1,0); % the empty subset
                elseif numel(remaining) < s
                    continue;
                else
                    retrCombos = nchoosek(remaining, s);
                end
                for rci = 1:size(retrCombos,1)
                    bitmap = zeros(1, retrievalWidth);
                    for c = 1:size(retrCombos,2)
                        bitmap(retrCombos(rci,c)) = 1;
                    end
                    for pp = 1:size(cachePerms,1)
                        SS = [SS; [cachePerms(pp,:), bitmap]];
                    end
                end
            end
        end
    end
catch ME
    line_error(mfilename,'The model is too large, LINE cannot generate the state space explicitly.');
end

end
