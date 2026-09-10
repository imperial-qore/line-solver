function varargout = getAvgCacheTable(self,varargin)
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
%   Item           item read by a per-item class of a cache network, 0 otherwise
%   Items          number of items managed by the cache
%   HitProb        (true) hit probability; per list on list rows
%   DelayedHitProb delayed-hit probability (retrieval system; total row only)
%   MissProb       miss probability (total row only)
%   HitRate        hit throughput = ArvR * HitProb
%   DelayedHitRate delayed-hit throughput = ArvR * DelayedHitProb
%   MissRate       miss throughput = ArvR * MissProb
%   ArvR           read-class arrival rate into the cache
%   ResidT         expected retrieval latency / residence time (NaN if not computed)
%   ListCost       mean storage cost held by the list (NaN unless item sizes are set)
%
% The hit class throughput reported by getAvgNodeTable aggregates true hits
% and delayed hits; this table separates them.
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
[varargout{1:max(nargout,1)}] = getAvgCacheTable_impl(self,varargin{:});
LineResultRecorder.capture(scope, self, 'cache', varargout{1});
end

function [CacheAvgTable] = getAvgCacheTable_impl(self)
% GETAVGCACHETABLE_IMPL Implementation of GETAVGCACHETABLE; see the wrapper above.

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
[Itemv, Listv, ListCapv, Itemsv, HitPv, DHitPv, MissPv, HitRv, DHitRv, MissRv, ArvRv, Latv, Costv] = deal([]);

for ind = caches
    np = sn.nodeparam{ind};
    hitclass = np.hitclass;
    nitems = 0;
    if isfield(np,'nitems'), nitems = np.nitems; end
    classitem = zeros(1, sn.nclasses);
    if isfield(np,'classitem') && ~isempty(np.classitem)
        nc = min(numel(np.classitem), sn.nclasses);
        classitem(1:nc) = np.classitem(1:nc);
    end
    itemcap = [];
    if isfield(np,'itemcap'), itemcap = np.itemcap(:).'; end
    h = numel(itemcap);
    totcap = sum(itemcap);
    node = self.model.nodes{ind};
    hitp = node.getHitRatio();
    missp = node.getMissRatio();
    dhitp = node.getDelayedHitRatio();
    hitplist = node.getHitRatioByList();
    lat = node.getResidT();
    listcost = [];
    if ismethod(node,'getListCost'), listcost = node.getListCost(); end
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
        totcost = NaN;
        if ~isempty(listcost), totcost = sum(listcost); end
        [Node, JobClass, Itemv, Listv, ListCapv, Itemsv, HitPv, DHitPv, MissPv, ...
            HitRv, DHitRv, MissRv, ArvRv, Latv, Costv] = addrow( ...
            Node, JobClass, Itemv, Listv, ListCapv, Itemsv, HitPv, DHitPv, MissPv, ...
            HitRv, DHitRv, MissRv, ArvRv, Latv, Costv, ...
            sn.nodenames{ind}, sn.classnames{r}, classitem(r), 0, totcap, nitems, ...
            ph, pd, pm, arvr*ph, arvr*pd, arvr*pm, arvr_retr, latr, totcost);

        % --- per-list rows (only if a multi-list breakdown is available) ---
        if h > 1 && ~isempty(hitplist) && r <= size(hitplist,1) && ~all(isnan(hitplist(r,:)))
            for l = 1:h
                phl = full(hitplist(r,l));
                if isnan(phl), phl = 0; end
                capl = NaN; if l <= numel(itemcap), capl = itemcap(l); end
                costl = NaN;
                if ~isempty(listcost) && l <= numel(listcost), costl = full(listcost(l)); end
                [Node, JobClass, Itemv, Listv, ListCapv, Itemsv, HitPv, DHitPv, MissPv, ...
                    HitRv, DHitRv, MissRv, ArvRv, Latv, Costv] = addrow( ...
                    Node, JobClass, Itemv, Listv, ListCapv, Itemsv, HitPv, DHitPv, MissPv, ...
                    HitRv, DHitRv, MissRv, ArvRv, Latv, Costv, ...
                    sn.nodenames{ind}, sn.classnames{r}, classitem(r), l, capl, nitems, ...
                    phl, NaN, NaN, arvr*phl, NaN, NaN, arvr, NaN, costl);
            end
        end
    end
end

Node = label(Node);
JobClass = label(JobClass);
Item = Itemv;
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
ListCost = Costv;
CacheAvgTable = Table(Node, JobClass, Item, List, ListCap, Items, HitProb, ...
    DelayedHitProb, MissProb, HitRate, DelayedHitRate, MissRate, ArvR, ResidT, ListCost);
CacheAvgTable = IndexedTable(CacheAvgTable);
end

function [Node, JobClass, Itemv, Listv, ListCapv, Itemsv, HitPv, DHitPv, MissPv, ...
    HitRv, DHitRv, MissRv, ArvRv, Latv, Costv] = addrow( ...
    Node, JobClass, Itemv, Listv, ListCapv, Itemsv, HitPv, DHitPv, MissPv, ...
    HitRv, DHitRv, MissRv, ArvRv, Latv, Costv, ...
    nodename, classname, itemidx, listidx, listcap, nitems, ...
    ph, pd, pm, hr, dhr, mr, arvr, lat, cost)
Node{end+1,1} = nodename;
JobClass{end+1,1} = classname;
Itemv(end+1,1) = itemidx;
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
Costv(end+1,1) = cost;
end

function v = nanGetAt(vec, r)
% Safe scalar read; returns NaN if out of range or empty.
if isempty(vec) || r > numel(vec)
    v = NaN;
else
    v = full(vec(r));
end
end
