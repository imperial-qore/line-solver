function [stateSpace,stateSpaceAggr,stateSpaceHashed,Dfilt,sn,basBlockQ,DfiltAux] = solver_ctmc_fcr_waitq(sn, options)
% [SS,SSA,SSH,DFILT,SN,BASBLOCKQ,DFILTAUX]=SOLVER_CTMC_FCR_WAITQ(SN,OPTIONS)
%
% DFILTAUX carries the derived START/PREEMPT filtrations, as in solver_ctmc.
% Here a single transition can start several services: the FIFO release
% cascade admits blocked jobs one after another, and each admission may take a
% server. The tags therefore accumulate along the cascade path and are written
% once, at the settled state, weighted by the same rate as the transition.
% Reachability-based state space and per-action rate filters for models with
% a finite capacity region (FCR) whose drop rule is WAITQ (waiting queue).
%
% JMT WAITQ semantics (reference, mirrored by LDES): a job refused entry to a
% full region leaves the upstream station and waits in a per-region FIFO of
% (class, destination) tokens outside the region; after every transition that
% frees region capacity, tokens are released strictly in FIFO order (head-of-
% line: a stuck head blocks the queue) as long as the admission constraints
% (global cap, per-class caps, memory budget, linear constraints A*x<=b)
% permit; a fresh arrival that satisfies the constraints is admitted even if
% the FIFO is non-empty (it overtakes a head stuck on a different constraint).
% Blocked jobs are counted neither in the region occupancy nor in any station
% state, so station QLen excludes them, matching the JMT report convention.
%
% True-BAS blocking between two stations is orthogonal to the region rule and
% is handled here as it is in the default generator: when a service completion
% at a BAS blocking station finds its destination unable to admit the job, the
% job is HELD AT THE SERVER (blocked marker set) rather than the departure
% being voided. Voiding it instead frees the server to re-serve the same job,
% which for exponential service is NOT equivalent to BAS: on release the held
% job enters the destination immediately, whereas a re-serving server must
% first draw a fresh completion, and throughput is understated. Those
% become-blocked arcs change the chain but are not departures, so they are
% accumulated separately in BASBLOCKQ and never enter DFILT.
%
% The CTMC state is augmented as [h(1:nstateful), buf_1, ..., buf_F] where h
% are the per-node hashed states and buf_f is the token FIFO of region f,
% padded with zeros to its maximum length. Classes whose region rule is DROP
% keep the transition-censoring behavior of the default generator (exact for
% memoryless sources, cross-validated against JMT).
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

nstateful = sn.nstateful;
K = sn.nclasses;
sync = sn.sync;
A = length(sync);
csmask = sn.csmask;
local = sn.nnodes + 1;

%% feature gates: combinations that would need semantics not defined here
if isfield(sn,'gsync') && ~isempty(sn.gsync)
    line_error(mfilename,'WAITQ finite capacity regions are not supported together with stochastic Petri net transitions in SolverCTMC.');
end
if isfield(sn,'fjsync') && ~isempty(sn.fjsync)
    line_error(mfilename,'WAITQ finite capacity regions are not supported together with fork-join in SolverCTMC.');
end
if any(sn.isstatedep(:,3))
    line_error(mfilename,'WAITQ finite capacity regions are not supported together with state-dependent routing in SolverCTMC.');
end

%% region data
F = sn.nregions;
memberMask = false(F, sn.nstations);   % (f,ist) true if station ist in region f
ccap = inf(F, K);                      % per-class caps
gcap = inf(F, 1);                      % global job cap
memcap = inf(F, 1);                    % global memory budget
szrow = ones(F, K);                    % per-class memory footprint
linA = cell(F,1); linb = cell(F,1);
iswaitq = false(F, K);                 % rule per (region, class): true=WAITQ-like, false=DROP
for f = 1:F
    Rmat = sn.region{f};               % M x (K+1)
    % membership: a station is a member if any job-count cap OR the region
    % memory budget is set on its row (a memory-only region has all job-count
    % entries at -1)
    memvec = -ones(sn.nstations,1);
    if isfield(sn,'regionmaxmem') && numel(sn.regionmaxmem) >= f && ~isempty(sn.regionmaxmem{f})
        memvec = sn.regionmaxmem{f}(:);
    end
    members = find(sn_region_members(sn, f, Rmat, memvec));
    memberMask(f, members) = true;
    for r = 1:K
        cv = Rmat(members, r); cv = cv(cv ~= -1);
        if ~isempty(cv); ccap(f,r) = min(cv); end
        iswaitq(f,r) = (sn.regionrule(f,r) ~= DropStrategy.DROP);
    end
    gv = Rmat(members, K+1); gv = gv(gv ~= -1);
    if ~isempty(gv); gcap(f) = min(gv); end
    if isfield(sn,'regionmaxmem') && numel(sn.regionmaxmem) >= f && ~isempty(sn.regionmaxmem{f})
        mv = sn.regionmaxmem{f}(members); mv = mv(mv ~= -1);
        if ~isempty(mv); memcap(f) = min(mv); end
    end
    if isfield(sn,'regionsz') && ~isempty(sn.regionsz)
        szrow(f,:) = sn.regionsz(f,:);
    end
    if isfield(sn,'regionlincon') && size(sn.regionlincon,1) >= f && ~isempty(sn.regionlincon{f,1})
        linA{f} = sn.regionlincon{f,1};
        linb{f} = sn.regionlincon{f,2};
    end
end

% token FIFO length bound per region: at most all closed jobs of WAITQ classes
% plus, per open WAITQ class, the state-space cutoff of that class
if isfield(options,'cutoff') && ~isempty(options.cutoff)
    cutoffMat = options.cutoff;
    if isscalar(cutoffMat)
        cutoffMat = cutoffMat * ones(sn.nstations, K);
    end
else
    cutoffMat = zeros(sn.nstations, K);
end
tokbound = zeros(1,K);
for r = 1:K
    if any(iswaitq(:,r))
        c_ = find(sn.chains(:,r), 1); % chain of class r
        chainpop = sum(sn.njobs(sn.chains(c_,:)));
        if isfinite(chainpop)
            % closed chain: jobs may switch into class r, so bound by the
            % whole chain population rather than njobs(r)
            tokbound(r) = chainpop;
        else
            tokbound(r) = max(cutoffMat(:,r));
        end
    end
end
Lmax = zeros(F,1);
for f = 1:F
    Lmax(f) = sum(tokbound(iswaitq(f,:)));
end
bufoff = nstateful + [0; cumsum(Lmax(1:end-1))]; % column offset of buf_f
width = nstateful + sum(Lmax);

%% initial augmented state (buffers empty)
h0 = zeros(1, nstateful);
for ind = 1:sn.nnodes
    if sn.isstateful(ind)
        isf = sn.nodeToStateful(ind);
        if sn.nodetype(ind) == NodeType.Source
            % canonical Source state (spaceGenerator convention); the stored
            % sn.state row may omit the Inf pool marker or the arrival phase
            hh_ = State.getHash(sn, ind, State.fromMarginal(sn, ind, []));
            h0(isf) = hh_(1);
        else
            h0(isf) = State.getHash(sn, ind, sn.state{isf}(1,:));
        end
        if h0(isf) <= 0
            line_error(mfilename, sprintf('Initial state of node %s not found in its local state space.', sn.nodenames{ind}));
        end
    end
end
row0 = zeros(1, width);
row0(1:nstateful) = h0;
for f = 1:F
    if any(regionAggr(h0, f) > ccap(f,:)) || violates(f, regionAggr(h0, f))
        line_error(mfilename,'The initial state violates the finite capacity region constraints.');
    end
end

%% breadth-first construction of the reachable augmented space
capRows = 1024;
SSH = zeros(capRows, width);
SSH(1,:) = row0;
nrows = 1;
keymap = dictionary(string(rowkey(row0)), 1);
frontier = 1;
% transition triplets per action
capTrip = 4096;
ta = zeros(capTrip,1); ti = zeros(capTrip,1); tj = zeros(capTrip,1); tv = zeros(capTrip,1);
ntrip = 0;
% derived START/PREEMPT triplets, tagged with the station and class they
% belong to: (station, class, from, to, rate)
capAux = 4096;
xk = zeros(capAux,1); xs = zeros(capAux,1); xr = zeros(capAux,1);
xi = zeros(capAux,1); xj = zeros(capAux,1); xv = zeros(capAux,1);
naux = 0;
M = sn.nstations;

while ~isempty(frontier)
    s = frontier(1); frontier(1) = [];
    row = SSH(s,:);
    h = row(1:nstateful);
    bufs = cell(F,1);
    for f = 1:F
        bf = row(bufoff(f)+1:bufoff(f)+Lmax(f));
        bufs{f} = bf(bf > 0);
    end
    % current per-region aggregate populations
    xf = zeros(F, K);
    for f = 1:F
        xf(f,:) = regionAggr(h, f);
    end
    for a = 1:A
        node_a = sync{a}.active{1}.node;
        isf_a = sn.nodeToStateful(node_a);
        class_a = sync{a}.active{1}.class;
        event_a = sync{a}.active{1}.event;
        [new_state_a, rate_a, ~, start_a, preempt_a] = State.afterEventHashed(sn, node_a, h(isf_a), event_a, class_a);
        if isequal(new_state_a, -1)
            continue
        end
        for ia = 1:length(new_state_a)
            if isnan(rate_a(ia)) || rate_a(ia) <= 0 || new_state_a(ia) == -1
                continue
            end
            tagS_a = stationTag(node_a, start_a, ia);
            tagP_a = stationTag(node_a, preempt_a, ia);
            node_p = sync{a}.passive{1}.node;
            if node_p == local
                newh = h;
                newh(isf_a) = new_state_a(ia);
                emit(a, s, newh, bufs, rate_a(ia), [], tagS_a, tagP_a);
            else
                class_p = sync{a}.passive{1}.class;
                event_p = sync{a}.passive{1}.event;
                isf_p = sn.nodeToStateful(node_p);
                % region-entry detection: passive station inside region f,
                % active node outside it, and the passive event is an arrival
                stat_a = 0;
                if node_a <= sn.nnodes && sn.isstation(node_a)
                    stat_a = sn.nodeToStation(node_a);
                end
                stat_p = 0;
                if sn.isstation(node_p)
                    stat_p = sn.nodeToStation(node_p);
                end
                blockedf = 0;
                droppedf = 0;
                if event_p == EventType.ARV && stat_p > 0
                    for f = 1:F
                        if memberMask(f, stat_p) && (stat_a <= 0 || ~memberMask(f, stat_a))
                            xn = xf(f,:);
                            xn(class_p) = xn(class_p) + 1;
                            if violates(f, xn)
                                if ~iswaitq(f, class_p)
                                    droppedf = f; % DROP rule: the job is destroyed
                                else
                                    blockedf = f;
                                end
                                break
                            end
                        end
                    end
                end
                if droppedf > 0
                    % DROP rule (JMT semantics): the refused job is destroyed;
                    % only the active (departing) part of the transition applies
                    newh = h;
                    newh(isf_a) = new_state_a(ia);
                    emit(a, s, newh, bufs, rate_a(ia) * sync{a}.passive{1}.prob, [], tagS_a, tagP_a);
                    continue
                end
                % see _kb/06-solver-catalog.md (CTMC section) for rationale
                switchf = 0;
                if blockedf == 0 && event_p == EventType.ARV && class_p ~= class_a ...
                        && stat_a > 0 && stat_p > 0
                    for f = 1:F
                        if memberMask(f, stat_a) && memberMask(f, stat_p)
                            switchf = f;
                            break
                        end
                    end
                end
                if switchf > 0
                    newh = h;
                    newh(isf_a) = new_state_a(ia);
                    emit(a, s, newh, bufs, rate_a(ia) * sync{a}.passive{1}.prob, ...
                        [switchf, class_p, node_p, iswaitq(switchf, class_p)], tagS_a, tagP_a);
                    continue
                end
                if blockedf > 0
                    if numel(bufs{blockedf}) >= Lmax(blockedf)
                        continue % FIFO truncation boundary (open-class cutoff)
                    end
                    newh = h;
                    newh(isf_a) = new_state_a(ia);
                    newbufs = bufs;
                    newbufs{blockedf}(end+1) = (node_p-1)*K + class_p;
                    % the refused job waits outside the region: it starts nothing
                    emit(a, s, newh, newbufs, rate_a(ia) * sync{a}.passive{1}.prob, [], tagS_a, tagP_a);
                else
                    if node_p == node_a % self-loop
                        [new_state_p, ~, outprob_p, start_p, preempt_p] = State.afterEventHashed(sn, node_p, new_state_a(ia), event_p, class_p);
                    else
                        [new_state_p, ~, outprob_p, start_p, preempt_p] = State.afterEventHashed(sn, node_p, h(isf_p), event_p, class_p);
                    end
                    if isempty(new_state_p) || isequal(new_state_p, -1)
                        % see _kb/06-solver-catalog.md (True BAS blocking) for rationale
                        if event_a == EventType.DEP && ~isempty(sn.isbasblocking) ...
                                && numel(sn.isbasblocking) >= node_a && sn.isbasblocking(node_a) == 1
                            curVecA = sn.space{isf_a}(h(isf_a),:);
                            if curVecA(end) == 0
                                blockedVec = curVecA;
                                blockedVec(end) = 1;
                                blockedIdx = matchrow(sn.space{isf_a}, blockedVec);
                                if blockedIdx > 0
                                    newh = h;
                                    newh(isf_a) = blockedIdx;
                                    % the job is held AT the server: nobody is
                                    % promoted, so this arc carries no tag
                                    emit(0, s, newh, bufs, rate_a(ia) * sync{a}.passive{1}.prob, [], zeros(M,K), zeros(M,K));
                                end
                            end
                        end
                        continue
                    end
                    for ip = 1:size(new_state_p,1)
                        if new_state_p(ip) == -1
                            continue
                        end
                        prob_sync_p = sync{a}.passive{1}.prob * outprob_p(ip);
                        if prob_sync_p <= 0
                            continue
                        end
                        if node_p < local && ~csmask(class_a, class_p) && rate_a(ia) * prob_sync_p > 0 && (sn.nodetype(node_p) ~= NodeType.Source)
                            line_error(mfilename, sprintf('Error: routing at node %d (%s) violates the class switching mask (class %s -> class %s).', node_a, sn.nodenames{node_a}, sn.classnames{class_a}, sn.classnames{class_p}));
                        end
                        newh = h;
                        newh(isf_a) = new_state_a(ia);
                        newh(isf_p) = new_state_p(ip);
                        emit(a, s, newh, bufs, rate_a(ia) * prob_sync_p, [], ...
                            tagS_a + stationTag(node_p, start_p, ip), ...
                            tagP_a + stationTag(node_p, preempt_p, ip));
                    end
                end
            end
        end
    end
end

SSH = SSH(1:nrows,:);

%% assemble outputs
Dfilt = cell(1,A);
for a = 1:A
    sel = (ta(1:ntrip) == a);
    Dfilt{a} = sparse(ti(sel), tj(sel), tv(sel), nrows, nrows);
end
% Sentinel action 0 collects the true-BAS become-blocked arcs: part of the
% generator, but not a departure of any action, so kept out of Dfilt.
selBas = (ta(1:ntrip) == 0);
basBlockQ = sparse(ti(selBas), tj(selBas), tv(selBas), nrows, nrows);
% derived START/PREEMPT filtrations, one matrix per (station, class)
DfiltAux.start = cell(M,K);
DfiltAux.preempt = cell(M,K);
for i_ = 1:M
    for r_ = 1:K
        selS = (xk(1:naux) == 1) & (xs(1:naux) == i_) & (xr(1:naux) == r_);
        DfiltAux.start{i_,r_} = sparse(xi(selS), xj(selS), xv(selS), nrows, nrows);
        selP = (xk(1:naux) == 2) & (xs(1:naux) == i_) & (xr(1:naux) == r_);
        DfiltAux.preempt{i_,r_} = sparse(xi(selP), xj(selP), xv(selP), nrows, nrows);
    end
end
stateSpaceHashed = SSH;
% full state matrix: concatenated per-node states plus the token buffers
cols = 0;
for isf = 1:nstateful
    cols = cols + size(sn.space{isf},2);
end
stateSpace = zeros(nrows, cols + sum(Lmax));
stateSpaceAggr = zeros(nrows, sn.nstations * K);
for s = 1:nrows
    pos = 0;
    for ind = 1:sn.nnodes
        if sn.isstateful(ind)
            isf = sn.nodeToStateful(ind);
            srow = sn.space{isf}(SSH(s,isf),:);
            stateSpace(s, pos+1:pos+length(srow)) = srow;
            pos = pos + size(sn.space{isf},2);
            if sn.isstation(ind)
                ist = sn.nodeToStation(ind);
                [~, nir] = State.toMarginal(sn, ind, srow);
                stateSpaceAggr(s, ((ist-1)*K+1):ist*K) = nir;
            end
        end
    end
    stateSpace(s, cols+1:end) = SSH(s, nstateful+1:end);
end

    function T = stationTag(node, tagrows, irow)
        % T=STATIONTAG(NODE,TAGROWS,IROW) lift one successor row's tag counts
        % to the (station x class) grid; a non-station node contributes none.
        T = zeros(M,K);
        if isempty(tagrows) || irow > size(tagrows,1) || node > sn.nnodes || ~sn.isstation(node)
            return
        end
        T(sn.nodeToStation(node),:) = tagrows(irow,:);
    end

    function tf = violates(f, x)
        % TF=VIOLATES(F,X) true if per-class population vector x breaks any
        % admission constraint of region f
        tf = any(x > ccap(f,:)) || sum(x) > gcap(f) || (x * szrow(f,:)') > memcap(f);
        if ~tf && ~isempty(linA{f})
            tf = any(linA{f} * x(:) > linb{f}(:));
        end
    end

    function x = regionAggr(hvec, f)
        % X=REGIONAGGR(HVEC,F) per-class population of region f under hashed
        % station states hvec
        x = zeros(1, K);
        for ist_ = find(memberMask(f,:))
            ind_ = sn.stationToNode(ist_);
            isf_ = sn.nodeToStateful(ind_);
            [~, nir_] = State.toMarginalAggr(sn, ind_, sn.space{isf_}(hvec(isf_),:));
            x = x + nir_(:)';
        end
    end

    function emit(a, src, newh, newbufs, w, pend, tagS, tagP)
        % EMIT(A,SRC,NEWH,NEWBUFS,W,PEND,TAGS,TAGP) applies the FIFO release
        % cascade to the tentative augmented state and records the transitions
        % of action a. PEND = [f cls destNode] is a pending gated re-entry (a
        % class-switching hop between members of region f): once the cascade
        % settles, the job of class cls is admitted at destNode if region f
        % has capacity, else parked at the tail of the FIFO.
        %
        % TAGS/TAGP are the (station x class) START/PREEMPT counts the two
        % halves of the synchronization have already produced; the cascade adds
        % the tags of every job it releases into a station, since an admitted
        % job may take a server there.
        if w <= 0
            return
        end
        if nargin < 6
            pend = [];
        end
        if nargin < 8
            tagS = zeros(M,K);
            tagP = zeros(M,K);
        end
        % work items: {h, bufs, prob, pend, startTag, preemptTag}
        work = {{newh, newbufs, 1.0, pend, tagS, tagP}};
        while ~isempty(work)
            it = work{1}; work(1) = [];
            hh = it{1}; bb = it{2}; pw = it{3}; pd = it{4};
            ts = it{5}; tp = it{6};
            progressed = false;
            for f_ = 1:F
                if isempty(bb{f_})
                    continue
                end
                x_ = regionAggr(hh, f_);
                tok = bb{f_}(1);
                dest = floor((tok-1)/K) + 1;
                r_ = mod(tok-1, K) + 1;
                xn_ = x_;
                xn_(r_) = xn_(r_) + 1;
                if violates(f_, xn_)
                    continue % head-of-line: this region's FIFO stays blocked
                end
                isf_d = sn.nodeToStateful(dest);
                [hd, ~, opd, sd, pdg] = State.afterEventHashed(sn, dest, hh(isf_d), EventType.ARV, r_);
                if isempty(hd) || isequal(hd, -1)
                    continue
                end
                for id = 1:length(hd)
                    if hd(id) == -1 || opd(id) <= 0
                        continue
                    end
                    hh2 = hh;
                    hh2(isf_d) = hd(id);
                    bb2 = bb;
                    bb2{f_}(1) = [];
                    work{end+1} = {hh2, bb2, pw * opd(id), pd, ...
                        ts + stationTag(dest, sd, id), tp + stationTag(dest, pdg, id)}; %#ok<AGROW>
                end
                progressed = true;
                break
            end
            if ~progressed && ~isempty(pd)
                % cascade settled: resolve the pending gated re-entry
                f_ = pd(1); cls_ = pd(2); dest_ = pd(3);
                x_ = regionAggr(hh, f_);
                xn_ = x_;
                xn_(cls_) = xn_(cls_) + 1;
                if violates(f_, xn_)
                    if numel(pd) >= 4 && ~pd(4)
                        % DROP rule: the switching job is destroyed
                        work{end+1} = {hh, bb, pw, [], ts, tp}; %#ok<AGROW>
                    else
                        % no capacity: park at the tail of the region FIFO
                        bb2 = bb;
                        bb2{f_}(end+1) = (dest_-1)*K + cls_;
                        work{end+1} = {hh, bb2, pw, [], ts, tp}; %#ok<AGROW>
                    end
                else
                    isf_d = sn.nodeToStateful(dest_);
                    [hd, ~, opd, sd, pdg] = State.afterEventHashed(sn, dest_, hh(isf_d), EventType.ARV, cls_);
                    admitted = false;
                    if ~isempty(hd) && ~isequal(hd, -1)
                        for id = 1:length(hd)
                            if hd(id) == -1 || opd(id) <= 0
                                continue
                            end
                            hh2 = hh;
                            hh2(isf_d) = hd(id);
                            work{end+1} = {hh2, bb, pw * opd(id), [], ...
                                ts + stationTag(dest_, sd, id), tp + stationTag(dest_, pdg, id)}; %#ok<AGROW>
                            admitted = true;
                        end
                    end
                    if ~admitted
                        % destination local state missing (e.g. station cap):
                        % park in the FIFO instead
                        bb2 = bb;
                        bb2{f_}(end+1) = (dest_-1)*K + cls_;
                        work{end+1} = {hh, bb2, pw, [], ts, tp}; %#ok<AGROW>
                    end
                end
                continue
            end
            if ~progressed
                % settled: register the augmented state and the transition
                rr = zeros(1, width);
                rr(1:nstateful) = hh;
                for f_ = 1:F
                    rr(bufoff(f_)+1:bufoff(f_)+numel(bb{f_})) = bb{f_};
                end
                kk = rowkey(rr);
                if isKey(keymap, kk)
                    dst = keymap(kk);
                else
                    nrows = nrows + 1;
                    if nrows > size(SSH,1)
                        SSH = [SSH; zeros(size(SSH,1), width)]; %#ok<AGROW>
                    end
                    SSH(nrows,:) = rr;
                    keymap(kk) = nrows;
                    dst = nrows;
                    frontier(end+1) = nrows; %#ok<AGROW>
                end
                ntrip = ntrip + 1;
                if ntrip > numel(ta)
                    ta = [ta; zeros(numel(ta),1)]; %#ok<AGROW>
                    ti = [ti; zeros(numel(ti),1)]; %#ok<AGROW>
                    tj = [tj; zeros(numel(tj),1)]; %#ok<AGROW>
                    tv = [tv; zeros(numel(tv),1)]; %#ok<AGROW>
                end
                ta(ntrip) = a; ti(ntrip) = src; tj(ntrip) = dst; tv(ntrip) = w * pw;
                % derived tags of this settled path, at the same rate
                for i_ = 1:M
                    for r_t = 1:K
                        if ts(i_,r_t) ~= 0
                            [xk,xs,xr,xi,xj,xv,naux] = auxpush(xk,xs,xr,xi,xj,xv,naux, 1, i_, r_t, src, dst, w*pw*ts(i_,r_t));
                        end
                        if tp(i_,r_t) ~= 0
                            [xk,xs,xr,xi,xj,xv,naux] = auxpush(xk,xs,xr,xi,xj,xv,naux, 2, i_, r_t, src, dst, w*pw*tp(i_,r_t));
                        end
                    end
                end
            end
        end
    end

end

function k = rowkey(v)
% K=ROWKEY(V) character key for an augmented state row
k = sprintf('%d,', v);
end

function [xk,xs,xr,xi,xj,xv,n] = auxpush(xk,xs,xr,xi,xj,xv,n, kind, ist, r, src, dst, val)
% Append one derived-filtration triplet; KIND is 1 for START, 2 for PREEMPT.
n = n + 1;
if n > numel(xk)
    xk = [xk; zeros(numel(xk),1)];
    xs = [xs; zeros(numel(xs),1)];
    xr = [xr; zeros(numel(xr),1)];
    xi = [xi; zeros(numel(xi),1)];
    xj = [xj; zeros(numel(xj),1)];
    xv = [xv; zeros(numel(xv),1)];
end
xk(n) = kind; xs(n) = ist; xr(n) = r; xi(n) = src; xj(n) = dst; xv(n) = val;
end
