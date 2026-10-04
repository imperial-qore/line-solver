function [outspace, outrate, outprob, eventCache, outstart, outpreempt] = afterEventStationPool(sn, ind, ist, inspace, event, class, isSimulation, eventCache, R, V, key, capacity, classcap)
% [OUTSPACE, OUTRATE, OUTPROB, EVENTCACHE, OUTSTART, OUTPREEMPT] = AFTEREVENTSTATIONPOOL(...)
%
% Event handler for a station with heterogeneous server pools, rewritten by
% solver_ctmc_pools (sn.nodeparam{ind}.ctmcpool). The local state is
% [buffer | servers | local variables] as at any FCFS-family station: the
% buffer holds the WAITING jobs only (class ids right-aligned, the rightmost is
% the oldest; per-class counts under SIRO), and the class-r server block holds
% one sub-block of phases per compatible pool, so the pool of a job in service
% is read from its phase. The semantics are those of LDES:
%
% - ARV: the candidate pools are the compatible pools with a free server. None:
%   the job waits. Otherwise the station's HeteroSchedPolicy picks the pool:
%   ORDER the lowest index, ALFS the least flexible pool (fewest compatible
%   classes, ties by index), FSF the pool with the highest rate for the class
%   (ties by index), RAIS each candidate with equal probability, ALIS and
%   FAIRNESS the first candidate of the rotating pool order, which then moves
%   to the back when there was more than one candidate.
% - DEP: the freed server of pool t takes the first waiting job, in the
%   discipline's service order, that pool t can serve. No rotation happens
%   there: the freed pool is the only one with a free server.
% - PHASE: a transition inside the block of the job's pool.
%
% Invariant kept by every event and by State.fromMarginalPool: no waiting job
% has a free compatible server.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

outspace = [];
outrate = [];
outprob = [];
outstart = [];
outpreempt = [];

pool = sn.nodeparam{ind}.ctmcpool;
K = sn.phasessz(ist,:);
Ks = sn.phaseshift(ist,:);
sched = sn.sched(ist);
isSiro = sched == SchedStrategy.SIRO;
nsrv = sum(K);

for row = 1:size(inspace,1)
    st = inspace(row,:);
    W = numel(st) - nsrv - V;
    buf = st(1:W);
    srv = st((W+1):(W+nsrv));
    var = st((W+nsrv+1):end);
    switch event
        case EventType.ARV
            [os, op, ot] = sub_arrival(sn, ist, pool, buf, srv, var, class, K, Ks, isSiro, capacity, classcap);
            outspace = sub_stack(outspace, os);
            outrate = [outrate; -1*ones(size(os,1),1)]; %#ok<AGROW>
            outprob = [outprob; op]; %#ok<AGROW>
            outstart = [outstart; ot]; %#ok<AGROW>
        case EventType.DEP
            if isempty(pool.pools{class})
                continue
            end
            for k = 1:numel(pool.pools{class})
                t = pool.pools{class}(k);
                for j = 1:pool.len{class}(k)
                    col = Ks(class) + pool.off{class}(k) + j;
                    nj = srv(col);
                    ex = pool.exit{class}{k}(j);
                    if nj <= 0 || ex <= 0
                        continue
                    end
                    srv_d = srv;
                    srv_d(col) = srv_d(col) - 1;
                    [bufs, srvs, ps, cls] = sub_serve_next(sn, ist, pool, buf, srv_d, t, K, Ks, isSiro, sched);
                    for b = 1:size(bufs,1)
                        outspace = sub_stack(outspace, [bufs(b,:), srvs(b,:), var]);
                        outrate = [outrate; ex*nj*ps(b)]; %#ok<AGROW>
                        outprob = [outprob; 1]; %#ok<AGROW>
                        tag = zeros(1,R);
                        if cls(b) > 0
                            tag(cls(b)) = 1;
                        end
                        outstart = [outstart; tag]; %#ok<AGROW>
                    end
                end
            end
        case EventType.PHASE
            if isempty(pool.pools{class})
                continue
            end
            for k = 1:numel(pool.pools{class})
                D0 = pool.D0{class}{k};
                for j = 1:pool.len{class}(k)
                    col = Ks(class) + pool.off{class}(k) + j;
                    nj = srv(col);
                    if nj <= 0
                        continue
                    end
                    for jd = 1:pool.len{class}(k)
                        if jd == j || D0(j,jd) <= 0
                            continue
                        end
                        srv_p = srv;
                        srv_p(col) = srv_p(col) - 1;
                        srv_p(Ks(class) + pool.off{class}(k) + jd) = srv_p(Ks(class) + pool.off{class}(k) + jd) + 1;
                        outspace = sub_stack(outspace, [buf, srv_p, var]);
                        outrate = [outrate; D0(j,jd)*nj]; %#ok<AGROW>
                        outprob = [outprob; 1]; %#ok<AGROW>
                        outstart = [outstart; zeros(1,R)]; %#ok<AGROW>
                    end
                end
            end
    end
end

outstart = State.tagPad(outstart, size(outspace,1), R);
outpreempt = State.tagPad(outpreempt, size(outspace,1), R);

if isSimulation
    if ~isnan(key) && isobject(eventCache)
        eventCache{key} = {outprob, outspace, outrate, outstart, outpreempt};
    end
    if size(outspace,1) > 1
        if event == EventType.ARV
            cum = cumsum(outprob) / sum(outprob);
        else
            cum = cumsum(outrate) / sum(outrate);
        end
        fc = 1 + max([0, find(rand > cum')]);
        outspace = outspace(fc,:);
        if event == EventType.ARV
            outrate = -1;
        else
            outrate = sum(outrate);
        end
        outprob = 1;
        outstart = outstart(fc,:);
        outpreempt = outpreempt(fc,:);
    end
end
end

function busy = sub_busy(pool, srv, Ks)
busy = zeros(1, pool.ntypes);
for r = 1:numel(pool.pools)
    for k = 1:numel(pool.pools{r})
        t = pool.pools{r}(k);
        cols = Ks(r) + pool.off{r}(k) + (1:pool.len{r}(k));
        busy(t) = busy(t) + sum(srv(cols));
    end
end
end

function [os, op, ot] = sub_arrival(sn, ist, pool, buf, srv, var, class, K, Ks, isSiro, capacity, classcap)
R = sn.nclasses;
os = zeros(0, numel(buf)+numel(srv)+numel(var));
op = zeros(0,1);
ot = zeros(0,R);
if isempty(pool.pools{class})
    return % the class is not served here
end
nir = sub_classcounts(buf, srv, K, Ks, isSiro, R);
ni = sum(nir);
if ni >= capacity(ist) || nir(class) >= classcap(ist,class)
    if State.isPhysicalCapacity(sn, ist, class) && State.arrivalIsLost(sn, ist, class)
        % a lost job leaves the state unchanged, so the arrival still fires
        os = [buf, srv, var];
        op = 1;
        ot = zeros(1,R);
    end
    return % blocked, or beyond the state-space cutoff
end
busy = sub_busy(pool, srv, Ks);
cand = pool.pools{class}(busy(pool.pools{class}) < pool.count(pool.pools{class}));
if isempty(cand)
    % no free compatible server: the job waits
    if isSiro
        buf(class) = buf(class) + 1;
    else
        slot = find(buf == 0, 1, 'last');
        if isempty(slot)
            buf = [0, buf];
            slot = 1;
        end
        buf(slot) = class;
    end
    os = [buf, srv, var];
    op = 1;
    ot = zeros(1,R);
    return
end
choice = cand(1);
pchoice = 1;
varc = var;
if numel(cand) > 1
    switch pool.policy
        case HeteroSchedPolicy.ALFS
            for t = pool.alfsorder
                if any(cand == t)
                    choice = t;
                    break
                end
            end
        case HeteroSchedPolicy.FSF
            [~, at] = max(pool.fsfrate(cand, class)); % first maximum wins ties
            choice = cand(at);
        case HeteroSchedPolicy.RAIS
            choice = cand;
            pchoice = ones(1, numel(cand)) / numel(cand);
        case {HeteroSchedPolicy.ALIS, HeteroSchedPolicy.FAIRNESS}
            if pool.rotate
                ord = pool.perms(var(pool.varpos),:);
                for k = 1:numel(ord)
                    if any(cand == ord(k))
                        choice = ord(k);
                        ord = [ord([1:k-1, k+1:end]), ord(k)];
                        break
                    end
                end
                [~, pidx] = ismember(ord, pool.perms, 'rows');
                varc(pool.varpos) = pidx;
            end
    end
end
for c = 1:numel(choice)
    t = choice(c);
    k = find(pool.pools{class} == t, 1);
    a = pool.alpha{class}{k};
    for j = 1:numel(a)
        if a(j) <= 0
            continue
        end
        srv_a = srv;
        col = Ks(class) + pool.off{class}(k) + j;
        srv_a(col) = srv_a(col) + 1;
        os = [os; buf, srv_a, varc]; %#ok<AGROW>
        op = [op; pchoice(c)*a(j)]; %#ok<AGROW>
        tag = zeros(1,R);
        tag(class) = 1;
        ot = [ot; tag]; %#ok<AGROW>
    end
end
end

function [bufs, srvs, ps, cls] = sub_serve_next(sn, ist, pool, buf, srv, t, K, Ks, isSiro, sched)
% The freed server of pool t takes the next compatible waiting job, if any.
R = sn.nclasses;
if isSiro
    elig = find(buf > 0 & pool.compat(t,:));
    if isempty(elig)
        bufs = buf; srvs = srv; ps = 1; cls = 0;
        return
    end
    tot = sum(buf(elig));
    bufs = zeros(0, numel(buf)); srvs = zeros(0, numel(srv)); ps = zeros(0,1); cls = zeros(0,1);
    for s = elig
        b = buf;
        b(s) = b(s) - 1;
        [bb, ss, pp] = sub_start(pool, b, srv, s, t, Ks);
        bufs = [bufs; bb]; srvs = [srvs; ss]; ps = [ps; pp * buf(s) / tot]; cls = [cls; s*ones(numel(pp),1)]; %#ok<AGROW>
    end
    return
end
pos = find(buf > 0);
pos = pos(pool.compat(t, buf(pos)));
if isempty(pos)
    bufs = buf; srvs = srv; ps = 1; cls = 0;
    return
end
switch sched
    case {SchedStrategy.HOL, SchedStrategy.FCFSPRIO, SchedStrategy.LCFSPRIO}
        prio = sn.classprio(buf(pos));
        pos = pos(prio == min(prio)); % a lower value is a higher priority
end
switch sched
    case {SchedStrategy.LCFS, SchedStrategy.LCFSPRIO}
        p = pos(1); % leftmost is the newest
    otherwise
        p = pos(end); % rightmost is the oldest
end
s = buf(p);
b = [0, buf(1:p-1), buf(p+1:end)];
[bufs, srvs, ps] = sub_start(pool, b, srv, s, t, Ks);
cls = s * ones(numel(ps),1);
end

function [bufs, srvs, ps] = sub_start(pool, buf, srv, s, t, Ks)
k = find(pool.pools{s} == t, 1);
a = pool.alpha{s}{k};
bufs = zeros(0, numel(buf)); srvs = zeros(0, numel(srv)); ps = zeros(0,1);
for j = 1:numel(a)
    if a(j) <= 0
        continue
    end
    ss = srv;
    col = Ks(s) + pool.off{s}(k) + j;
    ss(col) = ss(col) + 1;
    bufs = [bufs; buf]; srvs = [srvs; ss]; ps = [ps; a(j)]; %#ok<AGROW>
end
end

function nir = sub_classcounts(buf, srv, K, Ks, isSiro, R)
nir = zeros(1,R);
for r = 1:R
    nir(r) = sum(srv((Ks(r)+1):(Ks(r)+K(r))));
    if isSiro
        nir(r) = nir(r) + buf(r);
    else
        nir(r) = nir(r) + sum(buf == r);
    end
end
end

function out = sub_stack(out, rows)
% Stack rows of possibly different buffer widths, left-padding the narrower ones.
if isempty(rows)
    return
end
if isempty(out)
    out = rows;
elseif size(out,2) > size(rows,2)
    out = [out; zeros(size(rows,1), size(out,2)-size(rows,2)), rows];
elseif size(out,2) < size(rows,2)
    out = [zeros(size(out,1), size(rows,2)-size(out,2)), out; rows];
else
    out = [out; rows];
end
end
