function [ItemAvgTable] = getAvgItemTable(self)
% [ITEMAVGTABLE] = GETAVGITEMTABLE(SELF)
% Return a table of item-level cache occupancy, one row per Cache node, item
% and cache list (level). Columns:
%
%   Node       cache node name
%   Item       item index (1..nitems)
%   List       cache list (level) index (1..nlists)
%   ListCap    capacity of the list
%   Prob       steady-state probability that the item resides in that list
%
% The per-item, per-list occupancy is reported only where the solver computes
% a genuine per-item distribution (exact NC/MVA cache algorithms and the
% delayed-hit retrieval algorithms); other solvers return an empty table.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

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
[Itemv, Listv, ListCapv, Probv] = deal([]);

for ind = caches
    np = sn.nodeparam{ind};
    itemcap = [];
    if isfield(np,'itemcap'), itemcap = np.itemcap(:).'; end
    h = numel(itemcap);
    node = self.model.nodes{ind};
    itemprob = node.getItemProb();   % [nitems x (h+1)], col 1 = miss
    if isempty(itemprob) || h == 0
        continue
    end
    n = size(itemprob, 1);
    for i = 1:n
        for l = 1:h
            if (l+1) <= size(itemprob, 2)
                p = itemprob(i, l+1);
            else
                p = NaN;
            end
            Node{end+1,1} = sn.nodenames{ind}; %#ok<AGROW>
            Itemv(end+1,1) = i; %#ok<AGROW>
            Listv(end+1,1) = l; %#ok<AGROW>
            ListCapv(end+1,1) = itemcap(l); %#ok<AGROW>
            Probv(end+1,1) = p; %#ok<AGROW>
        end
    end
end

Node = label(Node);
Item = Itemv;
List = Listv;
ListCap = ListCapv;
Prob = Probv;
ItemAvgTable = Table(Node, Item, List, ListCap, Prob);
ItemAvgTable = IndexedTable(ItemAvgTable);
end
