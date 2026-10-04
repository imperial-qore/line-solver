function varargout = getAvgItemTable(self,varargin)
% [ITEMAVGTABLE] = GETAVGITEMTABLE(SELF)
% Return a table of item-level cache occupancy, one row per Cache node, item
% and cache list (level). Columns:
%
%   Node       cache node name
%   Item       item index (1..nitems)
%   List       cache list (level) index (1..nlists)
%   ListCap    capacity of the list
%   Size       storage cost of the item (NaN unless setItemSizes was called)
%   Prob       steady-state probability that the item resides in that list.
%              SolverCTMC reports the TIME-WEIGHTED (time-stationary) law, a
%              state reward of the exact stationary distribution; the NC/MVA
%              cache algorithms report the EMBEDDED (per-request) law of the
%              cache-content chain seen at request instants. The two coincide
%              only when requests see time averages (PASTA), so they differ once
%              service times distinguish hits from misses.
%   Cost       expected storage cost the item contributes to the list, Size*Prob,
%              so summing Cost over items reproduces the list's ListCost
%   DelayedHitQLen      mean secondary requests waiting on the in-flight fetch of
%                       the item (exact under SolverCTMC; NaN where not computed)
%   DelayedHitQLenFull  as above, including the request that triggered the fetch
%
% The per-item, per-list occupancy is reported only where the solver computes
% a genuine per-item distribution (the NC/MVA cache algorithms, isolated and
% integrated alike, the delayed-hit retrieval algorithms, and SolverCTMC); the
% simulators return an empty table.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
% The result recorder captures the returned table together with the solver
% that produced it, so cross-codebase parity is asserted against the values a
% solver RETURNED rather than the text it printed. Off unless a run asked for
% it (LineResultRecorder.enable), and then it costs one appdata lookup here.
% The wrapper exists so that recording happens on EVERY exit path, including
% the early returns inside the implementation below.
[scope, scopeGuard] = LineResultRecorder.enter(); %#ok<ASGLU>
[varargout{1:max(nargout,1)}] = getAvgItemTable_impl(self,varargin{:});
LineResultRecorder.capture(scope, self, 'item', varargout{1});
end

function [ItemAvgTable] = getAvgItemTable_impl(self)
% GETAVGITEMTABLE_IMPL Implementation of GETAVGITEMTABLE; see the wrapper above.

if GlobalConstants.DummyMode
    ItemAvgTable = IndexedTable(Table());
    return
end

% Ensure the model has been solved so the cache item probabilities are filled.
self.getAvgNode();

sn = self.model.getStruct;
caches = find(sn.nodetype == NodeType.Cache)';
if isempty(caches)
    ItemAvgTable = IndexedTable(Table());
    return
end

Node = {};
[Itemv, Listv, ListCapv, Sizev, Probv, Costv, DHQv, DHQFv] = deal([]);

for ind = caches
    np = sn.nodeparam{ind};
    itemcap = [];
    if isfield(np,'itemcap'), itemcap = np.itemcap(:).'; end
    h = numel(itemcap);
    itemsize = [];
    if isfield(np,'itemsize'), itemsize = np.itemsize(:).'; end
    node = self.model.nodes{ind};
    itemprob = node.getItemProb();   % [nitems x (h+1)], col 1 = miss
    [dhq, dhqf] = node.getDelayedHitQLen();  % [1 x nitems], empty when not computed
    % A solver may compute the per-item occupancy (NC/MVA), the delayed-hit queue
    % length (CTMC), or both; emit rows whenever either is available.
    if h == 0 || (isempty(itemprob) && isempty(dhq))
        continue
    end
    if ~isempty(itemprob)
        n = size(itemprob, 1);
    else
        n = numel(dhq);
    end
    for i = 1:n
        for l = 1:h
            if ~isempty(itemprob) && (l+1) <= size(itemprob, 2)
                p = itemprob(i, l+1);
            else
                p = NaN;
            end
            if numel(itemsize) >= i
                szi = itemsize(i);
            else
                szi = NaN;
            end
            Node{end+1,1} = sn.nodenames{ind}; %#ok<AGROW>
            Itemv(end+1,1) = i; %#ok<AGROW>
            Listv(end+1,1) = l; %#ok<AGROW>
            ListCapv(end+1,1) = itemcap(l); %#ok<AGROW>
            Sizev(end+1,1) = szi; %#ok<AGROW>
            Probv(end+1,1) = p; %#ok<AGROW>
            Costv(end+1,1) = szi * p; %#ok<AGROW>
            if numel(dhq) >= i
                DHQv(end+1,1) = dhq(i); %#ok<AGROW>
                DHQFv(end+1,1) = dhqf(i); %#ok<AGROW>
            else
                DHQv(end+1,1) = NaN; %#ok<AGROW>
                DHQFv(end+1,1) = NaN; %#ok<AGROW>
            end
        end
    end
end

Node = label(Node);
Item = Itemv;
List = Listv;
ListCap = ListCapv;
Size = Sizev;
Prob = Probv;
Cost = Costv;
DelayedHitQLen = DHQv;
DelayedHitQLenFull = DHQFv;
ItemAvgTable = Table(Node, Item, List, ListCap, Size, Prob, Cost, DelayedHitQLen, DelayedHitQLenFull);
ItemAvgTable = IndexedTable(ItemAvgTable);
end
