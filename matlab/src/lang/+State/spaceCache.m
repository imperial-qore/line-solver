function SS = spaceCache(n, m, retrievalSystemCapacity, maxPending, retrievalClassItems)
% SS = SPACECACHE(N, M, RETRIEVALSYSTEMCAPACITY, MAXPENDING, RETRIEVALCLASSITEMS)
%
% Generate the local-variable state space of a cache node.
%   n  : total number of distinct items
%   m  : per-level cache capacity (row vector)
%   retrievalSystemCapacity (optional, default 0): number of items that may be in
%     the retrieval system simultaneously. When > 0 each emitted row carries the
%     retrieval blocks described below after the cache contents.
%   maxPending (optional, default 0): maximum number of secondary (delayed-hit)
%     requests that may be merged onto the in-flight fetches of the cache.
%   retrievalClassItems (optional, default []): item index served by each
%     retrieval class, one entry per retrieval class.
%
% Local-variable layout: [cache contents | block A | block B]
%   cache contents  totalCacheCapacity columns holding 1..n item indices
%   block A         n columns, one per item; 1 iff a fetch of that item is in
%                   flight (the historical occupancy bitmap, semantics unchanged)
%   block B         one column per retrieval class holding the number of
%                   secondary requests of that class merged onto the in-flight
%                   fetch of its item; present only when maxPending > 0
%
% Block B keys the merged requests by retrieval class rather than by item so that
% the originating job class, hence its hit class, is recoverable when the fetch
% completes and the merged requests are released as delayed hits.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 3 || isempty(retrievalSystemCapacity)
    retrievalSystemCapacity = 0;
end
if nargin < 4 || isempty(maxPending)
    maxPending = 0;
end
if nargin < 5
    retrievalClassItems = [];
end

totalCacheCapacity = sum(m);
nItems = sum(n);
if retrievalSystemCapacity > 0
    widthA = nItems;
else
    widthA = 0;
    maxPending = 0;
end
% Block B is part of the layout whenever a retrieval system exists, so that the
% local-variable width matches sn.nvars; maxPending only bounds its counts.
retrievalClassItems = retrievalClassItems(:).';
if widthA > 0
    widthB = numel(retrievalClassItems);
else
    widthB = 0;
end
if widthB == 0
    maxPending = 0;
end
nVars = totalCacheCapacity + widthA + widthB;

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
            % size <= retrievalSystemCapacity.
            for s = 0:retrievalSystemCapacity
                if s == 0
                    retrCombos = zeros(1,0); % the empty subset
                elseif numel(remaining) < s
                    continue;
                else
                    retrCombos = nchoosek(remaining, s);
                end
                for rci = 1:size(retrCombos,1)
                    inflight = retrCombos(rci,:);
                    bitmapA = zeros(1, widthA);
                    for c = 1:numel(inflight)
                        bitmapA(inflight(c)) = 1;
                    end
                    % Only the retrieval classes of items being fetched can carry
                    % merged secondary requests; all other block-B slots are zero.
                    if widthB > 0 && ~isempty(inflight)
                        active = find(ismember(retrievalClassItems, inflight));
                    else
                        active = [];
                    end
                    blocksB = State.spaceCachePendings(numel(active), maxPending);
                    for bi = 1:size(blocksB,1)
                        blockB = zeros(1, widthB);
                        blockB(active) = blocksB(bi,:);
                        for pp = 1:size(cachePerms,1)
                            SS = [SS; [cachePerms(pp,:), bitmapA, blockB]];
                        end
                    end
                end
            end
        end
    end
catch ME
    line_error(mfilename,'The model is too large, LINE cannot generate the state space explicitly.');
end

end
