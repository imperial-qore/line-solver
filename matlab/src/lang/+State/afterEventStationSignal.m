function [outspace, outrate, outprob] = afterEventStationSignal(sn, ind, ist, inspace, class, K, Ks, S, pie, space_buf, space_srv, space_var)
% [OUTSPACE, OUTRATE, OUTPROB] = AFTEREVENTSTATIONSIGNAL(SN, IND, IST, INSPACE, CLASS, K, KS, S, PIE, SPACE_BUF, SPACE_SRV, SPACE_VAR)
%
% Passive arrival of a G-network signal class at station IST. A signal never
% joins the station: it acts on the jobs already there and is annihilated.
%
% CATASTROPHE empties the station of every job, ignoring the removal count
% distribution (a catastrophe removes all jobs by definition).
%
% NEGATIVE removes a batch of jobs:
%   - Victims. When the signal declares a target class (setmodel via
%     forJobClass, sn.signaltarget >= 1) only that class is eligible;
%     otherwise every non-signal class is eligible. The untargeted case is
%     the classic Gelenbe negative customer and agrees with SolverMAM and
%     SolverLDES, which are both class-agnostic.
%   - Count. sn.signalremdist gives the batch-size pmf. An oversized batch
%     empties the eligible jobs rather than driving the queue negative, so
%     the pmf tail P(B >= n) lumps onto "remove all n" (the min(B,n)
%     clipping used by LDES and the tail term used by MAM).
%   - Policy. sn.signalrempolicy selects the victim: FCFS removes the oldest
%     waiting job, LCFS the newest, RANDOM draws uniformly over waiting and
%     in-service jobs. FCFS/LCFS only touch an in-service job when no job is
%     waiting (two-tier, as in LDES).
%
% The station state carries no arrival-time order for in-service jobs, so
% when a policy has to reach into the servers the victim is drawn uniformly
% across the occupied phases (exact whenever at most one job of the class is
% in service, which covers every single-server station).
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

[~, nirm, sirm] = State.toMarginal(sn, ind, inspace, K, Ks, space_buf, space_srv, space_var);

% CATASTROPHE: derived from signaltype as well as the iscatastrophe flag.
if State.isCatastropheSignal(sn, class)
    outspace = [zeros(1, size(space_buf, 2)), zeros(1, size(space_srv, 2)), space_var];
    outrate = -1;
    outprob = 1;
    return
end

% Eligible victim classes.
tgt = -1;
if isfield(sn, 'signaltarget') && ~isempty(sn.signaltarget) && numel(sn.signaltarget) >= class
    tgt = sn.signaltarget(class);
end
if tgt >= 1
    tgtclasses = tgt;
else
    tgtclasses = find(~sn.issignal(:)');
end
tgtclasses = tgtclasses(nirm(tgtclasses) > 0);
ntot = sum(nirm(tgtclasses));
if isempty(tgtclasses) || ntot <= 0
    outspace = [space_buf, space_srv, space_var]; % no victim: the signal vanishes
    outrate = -1;
    outprob = 1;
    return
end

% Batch size pmf, clipped at the number of eligible jobs.
[kvals, kprobs] = State.signalBatchPMF(sn, class, ntot);

policy = RemovalPolicy.RANDOM;
if isfield(sn, 'signalrempolicy') && ~isempty(sn.signalrempolicy) && numel(sn.signalrempolicy) >= class
    policy = sn.signalrempolicy(class);
end

outspace = [];
outprob = [];
for ik = 1:numel(kvals)
    if kprobs(ik) <= 0
        continue
    end
    if kvals(ik) <= 0
        outspace = [outspace; space_buf, space_srv, space_var]; %#ok<AGROW>
        outprob = [outprob; kprobs(ik)]; %#ok<AGROW>
        continue
    end
    [sp, pr] = removeBatch(sn, ist, space_buf, space_srv, space_var, ...
        kvals(ik), tgtclasses, policy, K, Ks, S, pie);
    outspace = [outspace; sp]; %#ok<AGROW>
    outprob = [outprob; kprobs(ik) * pr]; %#ok<AGROW>
end

% Merge duplicate destination states so the generator sees one entry each.
[outspace, outprob] = mergeStates(outspace, outprob);
outrate = -1 * ones(size(outspace, 1), 1); % passive action
end

function [outspace, outprob] = removeBatch(sn, ist, buf, srv, var, k, tgtclasses, policy, K, Ks, S, pie)
% Remove k jobs one at a time; sequential uniform draws without replacement
% reproduce a uniform choice of the removed subset.
outspace = [buf, srv, var];
outprob = 1;
for step = 1:k
    nextspace = [];
    nextprob = [];
    for row = 1:size(outspace, 1)
        b = outspace(row, 1:size(buf, 2));
        s = outspace(row, size(buf, 2) + (1:size(srv, 2)));
        v = outspace(row, size(buf, 2) + size(srv, 2) + 1:end);
        [sp, pr] = removeOne(sn, ist, b, s, v, tgtclasses, policy, K, Ks, S, pie);
        if isempty(sp)
            % nothing left to remove: the state is already drained
            nextspace = [nextspace; outspace(row, :)]; %#ok<AGROW>
            nextprob = [nextprob; outprob(row)]; %#ok<AGROW>
        else
            nextspace = [nextspace; sp]; %#ok<AGROW>
            nextprob = [nextprob; outprob(row) * pr]; %#ok<AGROW>
        end
    end
    [outspace, outprob] = mergeStates(nextspace, nextprob);
end
end

function [outspace, outprob] = removeOne(sn, ist, buf, srv, var, tgtclasses, policy, K, Ks, S, pie)
% Enumerate the single-victim outcomes and their probabilities.
outspace = [];
outprob = [];

[waitPos, waitClass, waitWeight, isOrdered, isPairBuf] = waitingVictims(sn, ist, buf, tgtclasses);
[srvClass, srvPhase, srvCount] = inServiceVictims(srv, tgtclasses, K, Ks);
nwait = sum(waitWeight);
nsrv = sum(srvCount);
if nwait == 0 && nsrv == 0
    return
end

% FCFS/LCFS rank the waiting line by age, which only an ordered buffer
% records; a per-class count buffer (SIRO/SEPT/LEPT) carries no age, so an
% age-based policy degenerates to a uniform draw over the waiting jobs.
ageOrdered = isOrdered && (policy == RemovalPolicy.FCFS || policy == RemovalPolicy.LCFS);
if ageOrdered && nwait > 0
    if policy == RemovalPolicy.FCFS
        [~, pick] = max(waitPos); % head of line: the last occupied slot
    else
        [~, pick] = min(waitPos); % most recent arrival: the first occupied slot
    end
    [b2, s2] = dropWaiting(buf, srv, waitPos(pick), isOrdered, isPairBuf, waitClass(pick));
    outspace = [b2, s2, var];
    outprob = 1;
    return
end

if policy == RemovalPolicy.RANDOM
    total = nwait + nsrv; % uniform over waiting and in-service jobs alike
else
    total = nwait; % FCFS/LCFS drain the waiting line before the servers
    if total == 0
        total = nsrv;
    end
end

if nwait > 0
    for w = 1:numel(waitPos)
        [b2, s2] = dropWaiting(buf, srv, waitPos(w), isOrdered, isPairBuf, waitClass(w));
        outspace = [outspace; b2, s2, var]; %#ok<AGROW>
        outprob = [outprob; waitWeight(w) / total]; %#ok<AGROW>
    end
end
if policy == RemovalPolicy.RANDOM || nwait == 0
    for j = 1:numel(srvClass)
        if srvCount(j) <= 0
            continue
        end
        [b2, s2] = dropInService(sn, ist, buf, srv, srvClass(j), srvPhase(j), K, Ks, S, pie, isOrdered, isPairBuf);
        outspace = [outspace; b2, s2, var]; %#ok<AGROW>
        outprob = [outprob; srvCount(j) / total]; %#ok<AGROW>
    end
end
[outspace, outprob] = mergeStates(outspace, outprob);
end

function [pos, cls, weight, isOrdered, isPairBuf] = waitingVictims(sn, ist, buf, tgtclasses)
% Eligible waiting jobs, as (position, class, multiplicity). Three buffer
% layouts are in use: an ordered list of class ids (FCFS family), an ordered
% list of (class, phase) pairs (preemptive family), and a per-class count
% vector (SIRO/SEPT/LEPT). Stations that admit every job into service have
% no waiting line at all. Only the ordered layouts record arrival order.
pos = [];
cls = [];
weight = [];
isOrdered = false;
isPairBuf = false;
if isempty(buf)
    return
end
switch sn.sched(ist)
    case {SchedStrategy.FCFS, SchedStrategy.HOL, SchedStrategy.LCFS, SchedStrategy.LCFSPRIO}
        isOrdered = true;
        for c = 1:numel(buf)
            if buf(c) > 0 && ismember(buf(c), tgtclasses)
                pos(end+1) = c; %#ok<AGROW>
                cls(end+1) = buf(c); %#ok<AGROW>
                weight(end+1) = 1; %#ok<AGROW>
            end
        end
    case {SchedStrategy.FCFSPR, SchedStrategy.FCFSPI, SchedStrategy.FCFSPRPRIO, SchedStrategy.FCFSPIPRIO, ...
            SchedStrategy.LCFSPR, SchedStrategy.LCFSPI, SchedStrategy.LCFSPRPRIO, SchedStrategy.LCFSPIPRIO}
        isOrdered = true;
        isPairBuf = true;
        for c = 1:2:numel(buf)
            if buf(c) > 0 && ismember(buf(c), tgtclasses)
                pos(end+1) = c; %#ok<AGROW>
                cls(end+1) = buf(c); %#ok<AGROW>
                weight(end+1) = 1; %#ok<AGROW>
            end
        end
    case {SchedStrategy.SIRO, SchedStrategy.SEPT, SchedStrategy.LEPT}
        % per-class counts: each eligible class is one victim kind carrying
        % as many interchangeable jobs as its count
        for r = tgtclasses
            if r <= numel(buf) && buf(r) > 0
                pos(end+1) = r; %#ok<AGROW>
                cls(end+1) = r; %#ok<AGROW>
                weight(end+1) = buf(r); %#ok<AGROW>
            end
        end
    otherwise
        % PS/INF/DPS/GPS/LPS and the source hold no waiting jobs
end
end

function [cls, phase, count] = inServiceVictims(srv, tgtclasses, K, Ks)
cls = [];
phase = [];
count = [];
for r = tgtclasses
    for p = 1:K(r)
        n = srv(Ks(r) + p);
        if n > 0
            cls(end+1) = r; %#ok<AGROW>
            phase(end+1) = p; %#ok<AGROW>
            count(end+1) = n; %#ok<AGROW>
        end
    end
end
end

function [buf, srv] = dropWaiting(buf, srv, pos, isOrdered, isPairBuf, cls)
% Remove a waiting job, keeping the layout invariant (empty slots pad the
% left, the head of line stays rightmost).
if isPairBuf
    buf([pos, pos+1]) = [];
    buf = [0, 0, buf];
elseif isOrdered
    buf(pos) = [];
    buf = [0, buf];
else
    buf(cls) = buf(cls) - 1; % per-class count buffer
end
end

function [buf, srv] = dropInService(sn, ist, buf, srv, cls, phase, K, Ks, S, pie, isOrdered, isPairBuf)
% Remove an in-service job and, at a station that keeps a waiting line, pull
% the head of line into the freed server.
srv(Ks(cls) + phase) = srv(Ks(cls) + phase) - 1;
if isempty(buf) || sum(srv) >= S(ist)
    return
end
switch sn.sched(ist)
    case {SchedStrategy.FCFS, SchedStrategy.HOL, SchedStrategy.LCFS, SchedStrategy.LCFSPRIO}
        headpos = find(buf > 0, 1, 'last');
        if ~isempty(headpos)
            promo = buf(headpos);
            % the slot is vacated, not blanked in place: the waiting line
            % stays right-aligned with the empty slots padding the left
            buf(headpos) = [];
            buf = [0, buf];
            srv(Ks(promo) + 1) = srv(Ks(promo) + 1) + 1;
        end
    case {SchedStrategy.FCFSPR, SchedStrategy.FCFSPI, SchedStrategy.FCFSPRPRIO, SchedStrategy.FCFSPIPRIO, ...
            SchedStrategy.LCFSPR, SchedStrategy.LCFSPI, SchedStrategy.LCFSPRPRIO, SchedStrategy.LCFSPIPRIO}
        headpos = [];
        for c = 1:2:numel(buf)
            if buf(c) > 0
                headpos = c;
            end
        end
        if ~isempty(headpos)
            promo = buf(headpos);
            promophase = buf(headpos+1);
            if promophase < 1
                promophase = 1;
            end
            buf([headpos, headpos+1]) = [];
            buf = [0, 0, buf];
            srv(Ks(promo) + promophase) = srv(Ks(promo) + promophase) + 1; % resumes at its stored phase
        end
    case {SchedStrategy.SIRO, SchedStrategy.SEPT, SchedStrategy.LEPT}
        % the freed server takes a waiting job; the count buffer carries no
        % order, so any waiting class is equally eligible. Promote the
        % lowest-indexed waiting class to keep the map single-valued: the
        % SIRO service order is itself resolved by the rate function.
        promo = find(buf > 0, 1, 'first');
        if ~isempty(promo)
            buf(promo) = buf(promo) - 1;
            srv(Ks(promo) + 1) = srv(Ks(promo) + 1) + 1;
        end
    otherwise
        % no waiting line to promote from
end
end

function [space, prob] = mergeStates(space, prob)
if isempty(space)
    return
end
[u, ~, ic] = unique(space, 'rows', 'stable');
p = zeros(size(u, 1), 1);
for i = 1:numel(ic)
    p(ic(i)) = p(ic(i)) + prob(i);
end
space = u;
prob = p;
end
