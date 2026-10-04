function space = fromMarginalPool(sn, ind, n)
% SPACE = FROMMARGINALPOOL(SN, IND, N)
%
% Local states with marginal N at a station with heterogeneous server pools,
% as rewritten for the CTMC by solver_ctmc_pools (see
% State.afterEventStationPool for the layout and the semantics). A row is
% [buffer | servers | pool order]: every split of the jobs into in-service
% (per class and pool, within the pool sizes) and waiting ones such that no
% waiting job has a free compatible server, every phase assignment of the jobs
% in service, every order of the waiting jobs (per-class counts under SIRO),
% and, under ALIS/FAIRNESS with a rotating pool order, every order index.
%
% The ordered buffer is W = max(1, min(sum(N), cap)) wide, so the widest
% marginal fixes the width of the station space, as at an FCFS station.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

ist = sn.nodeToStation(ind);
R = sn.nclasses;
pool = sn.nodeparam{ind}.ctmcpool;
K = sn.phasessz(ist,:);
V = sum(sn.nvars(ind,:));
if V ~= double(pool.rotate)
    line_error(mfilename, sprintf('Station ''%s'' carries local variables that a heterogeneous server pool cannot hold.', sn.nodenames{ind}));
end
n = n(:)';
isSiro = sn.sched(ist) == SchedStrategy.SIRO;
cap = sn.cap(ist);
if isSiro
    W = R;
else
    W = max(1, min(sum(n), cap));
end
nsrv = sum(K);
unserved = cellfun(@isempty, pool.pools);
if sum(n) > cap || any(n(unserved) > 0)
    space = zeros(0, W + nsrv + V);
    return
end
% in-service allocations x(r,k): class r on the k-th compatible pool of r
alloc = sub_allocations(pool, n, R);
space = zeros(0, W + nsrv);
for a = 1:numel(alloc)
    x = alloc{a};
    w = n;
    for r = 1:R
        if ~isempty(x{r})
            w(r) = w(r) - sum(x{r});
        end
    end
    % server part: phases of the jobs of each class-pool block
    srvset = zeros(1,0);
    for r = 1:R
        blk = zeros(1,0);
        if isempty(pool.pools{r})
            blk = zeros(1, K(r));
        else
            for k = 1:numel(pool.pools{r})
                blk = State.cartesian(blk, State.spaceClosedSingle(pool.len{r}(k), x{r}(k)));
            end
            if size(blk,2) < K(r)
                blk = [blk, zeros(size(blk,1), K(r)-size(blk,2))]; %#ok<AGROW>
            end
        end
        srvset = State.cartesian(srvset, blk);
    end
    % buffer part
    if isSiro
        bufset = w;
    elseif sum(w) == 0
        bufset = zeros(1, W);
    else
        vi = [];
        for r = 1:R
            vi = [vi, r*ones(1, w(r))]; %#ok<AGROW>
        end
        mi = unique(multiset_perms(vi), 'rows');
        bufset = [zeros(size(mi,1), W - size(mi,2)), mi];
    end
    space = [space; State.cartesian(bufset, srvset)]; %#ok<AGROW>
end
if pool.rotate
    space = State.cartesian(space, (1:size(pool.perms,1))');
end
space = unique(space, 'rows');
space = space(end:-1:1,:);
end

function alloc = sub_allocations(pool, n, R)
% All x{r}(k) >= 0 with sum_k x{r}(k) <= n(r), pool loads within pool.count,
% and every class with waiting jobs finding each compatible pool full.
alloc = {};
x = cell(1,R);
load = zeros(1, pool.ntypes);
alloc = sub_rec(pool, n, R, 1, 1, x, load, alloc);
end

function alloc = sub_rec(pool, n, R, r, k, x, load, alloc)
if r > R
    % invariant check: a class with waiting jobs has no free compatible server
    for s = 1:R
        if isempty(pool.pools{s})
            continue
        end
        if n(s) - sum(x{s}) > 0 && any(load(pool.pools{s}) < pool.count(pool.pools{s}))
            return
        end
    end
    alloc{end+1} = x;
    return
end
if isempty(pool.pools{r})
    alloc = sub_rec(pool, n, R, r+1, 1, x, load, alloc);
    return
end
if k > numel(pool.pools{r})
    alloc = sub_rec(pool, n, R, r+1, 1, x, load, alloc);
    return
end
t = pool.pools{r}(k);
used = 0;
if k > 1
    used = sum(x{r}(1:k-1));
end
for c = 0:min(n(r) - used, pool.count(t) - load(t))
    x{r}(k) = c;
    load(t) = load(t) + c;
    alloc = sub_rec(pool, n, R, r, k+1, x, load, alloc);
    load(t) = load(t) - c;
end
x{r}(k) = 0;
end
