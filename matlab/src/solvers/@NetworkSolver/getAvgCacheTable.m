function [CacheAvgTable] = getAvgCacheTable(self)
% [CACHEAVGTABLE] = GETAVGCACHETABLE(SELF)
% Return a table of detailed per-class performance metrics for every Cache
% node in the model. For each cache node and read (input) class there is a
% total row (List=0) and, where the solver reports per-list (per-level) hit
% probabilities and the cache has more than one list, one extra row per list.
% Columns:
%
%   Node           cache node name
%   JobClass       read (input) class name
%   List           0 = total over all lists; l = cache list (level) l
%   ListCap        capacity of the list (total capacity on the total row)
%   Items          number of items managed by the cache
%   HitProb        (true) hit probability; per list on list rows
%   DelayedHitProb delayed-hit probability (retrieval system; total row only)
%   MissProb       miss probability (total row only)
%   HitRate        hit throughput = ArvR * HitProb
%   DelayedHitRate delayed-hit throughput = ArvR * DelayedHitProb
%   MissRate       miss throughput = ArvR * MissProb
%   ArvR           read-class arrival rate into the cache
%   ResidT         expected retrieval latency / residence time (NaN if not computed)
%
% The hit class throughput reported by getAvgNodeTable aggregates true hits
% and delayed hits; this table separates them.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if GlobalConstants.DummyMode
    CacheAvgTable = IndexedTable(Table());
    return
end

sn = self.model.getStruct;
K = sn.nclasses;

caches = find(sn.nodetype == NodeType.Cache)';
if isempty(caches)
    CacheAvgTable = IndexedTable(Table());
    return
end

% Node throughputs. The cache read-class arrival equals that class's source
% throughput (every read request enters the cache); robust across solvers,
% including simulators where delayed hits are not folded into hit/miss tput.
[~,~,~,TNn] = self.getAvgNode();
srcNode = find(sn.nodetype == NodeType.Source, 1);

Node = {};
JobClass = {};
[Listv, ListCapv, Itemsv, HitPv, DHitPv, MissPv, HitRv, DHitRv, MissRv, ArvRv, Latv] = deal([]);

for ind = caches
    np = sn.nodeparam{ind};
    hitclass = np.hitclass;
    nitems = 0;
    if isfield(np,'nitems'), nitems = np.nitems; end
    itemcap = [];
    if isfield(np,'itemcap'), itemcap = np.itemcap(:).'; end
    h = numel(itemcap);
    totcap = sum(itemcap);
    node = self.model.nodes{ind};
    hitp = node.getHitRatio();
    missp = node.getMissRatio();
    dhitp = []; % delayed-hit retrieval removed; always 0 for plain caches
    hitplist = node.getHitRatioByList();
    lat = node.getResidT();
    for r = 1:K
        % read (input) classes are those with a defined hit class
        if r > length(hitclass) || hitclass(r) <= 0
            continue
        end
        ph = nanGetAt(hitp, r);
        pm = nanGetAt(missp, r);
        pd = nanGetAt(dhitp, r);
        if isnan(ph) && isnan(pm) && isnan(pd)
            continue % no solved cache metrics for this class
        end
        if isnan(ph), ph = 0; end
        if isnan(pm), pm = 0; end
        if isnan(pd), pd = 0; end
        arvr = 0;
        if ~isempty(TNn) && ~isempty(srcNode) && srcNode <= size(TNn,1) && r <= size(TNn,2)
            arvr = TNn(srcNode, r);
        end
        latr = nanGetAt(lat, r);

        % --- total row (List = 0) ---
        % ArvR is the retrieval-system throughput arvr*(missprob+delayedprob),
        % Little-consistent with ResidT. see _kb/09-ldes-and-cache.md for rationale
        arvr_retr = arvr*(pm+pd);
        [Node, JobClass, Listv, ListCapv, Itemsv, HitPv, DHitPv, MissPv, ...
            HitRv, DHitRv, MissRv, ArvRv, Latv] = addrow( ...
            Node, JobClass, Listv, ListCapv, Itemsv, HitPv, DHitPv, MissPv, ...
            HitRv, DHitRv, MissRv, ArvRv, Latv, ...
            sn.nodenames{ind}, sn.classnames{r}, 0, totcap, nitems, ...
            ph, pd, pm, arvr*ph, arvr*pd, arvr*pm, arvr_retr, latr);

        % --- per-list rows (only if a multi-list breakdown is available) ---
        if h > 1 && ~isempty(hitplist) && r <= size(hitplist,1) && ~all(isnan(hitplist(r,:)))
            for l = 1:h
                phl = full(hitplist(r,l));
                if isnan(phl), phl = 0; end
                capl = NaN; if l <= numel(itemcap), capl = itemcap(l); end
                [Node, JobClass, Listv, ListCapv, Itemsv, HitPv, DHitPv, MissPv, ...
                    HitRv, DHitRv, MissRv, ArvRv, Latv] = addrow( ...
                    Node, JobClass, Listv, ListCapv, Itemsv, HitPv, DHitPv, MissPv, ...
                    HitRv, DHitRv, MissRv, ArvRv, Latv, ...
                    sn.nodenames{ind}, sn.classnames{r}, l, capl, nitems, ...
                    phl, NaN, NaN, arvr*phl, NaN, NaN, arvr, NaN);
            end
        end
    end
end

Node = label(Node);
JobClass = label(JobClass);
List = Listv;
ListCap = ListCapv;
Items = Itemsv;
HitProb = HitPv;
DelayedHitProb = DHitPv;
MissProb = MissPv;
HitRate = HitRv;
DelayedHitRate = DHitRv;
MissRate = MissRv;
ArvR = ArvRv;
ResidT = Latv;
CacheAvgTable = Table(Node, JobClass, List, ListCap, Items, HitProb, ...
    DelayedHitProb, MissProb, HitRate, DelayedHitRate, MissRate, ArvR, ResidT);
CacheAvgTable = IndexedTable(CacheAvgTable);
end

function [Node, JobClass, Listv, ListCapv, Itemsv, HitPv, DHitPv, MissPv, ...
    HitRv, DHitRv, MissRv, ArvRv, Latv] = addrow( ...
    Node, JobClass, Listv, ListCapv, Itemsv, HitPv, DHitPv, MissPv, ...
    HitRv, DHitRv, MissRv, ArvRv, Latv, ...
    nodename, classname, listidx, listcap, nitems, ...
    ph, pd, pm, hr, dhr, mr, arvr, lat)
Node{end+1,1} = nodename;
JobClass{end+1,1} = classname;
Listv(end+1,1) = listidx;
ListCapv(end+1,1) = listcap;
Itemsv(end+1,1) = nitems;
HitPv(end+1,1) = ph;
DHitPv(end+1,1) = pd;
MissPv(end+1,1) = pm;
HitRv(end+1,1) = hr;
DHitRv(end+1,1) = dhr;
MissRv(end+1,1) = mr;
ArvRv(end+1,1) = arvr;
Latv(end+1,1) = lat;
end

function v = nanGetAt(vec, r)
% Safe scalar read; returns NaN if out of range or empty.
if isempty(vec) || r > numel(vec)
    v = NaN;
else
    v = full(vec(r));
end
end
