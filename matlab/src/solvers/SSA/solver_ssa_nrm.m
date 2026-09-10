function [QN, UN, RN, TN, CN, XN, lG, sn, StartN, PreemptN] = solver_ssa_nrm(sn, options)
% SOLVER_SSA_NRM   Steady‑state analysis via the Next‑Reaction Method (SSA)
%
%   [QN, UN, RN, TN, CN, XN, LG, SN] = SOLVER_SSA_NRM(SN, OPTIONS)
%   runs a stochastic simulation of the queueing network described in SN
%   for OPTIONS.samples reaction firings using Gibson & Bruck's
%   Next‑Reaction Method.  During the run it:
%     • computes performance metrics directly during simulation;
%     • returns standard queueing performance measures.
%
%   Outputs
%     QN        – M×K matrix of mean queue lengths
%     UN        – M×K matrix of utilizations
%     RN        – M×K matrix of response times
%     TN        – M×K matrix of throughputs
%     CN        – 1×K vector of cycle times
%     XN        – 1×K vector of system throughputs
%     LG        – Logarithm of normalizing constant (not computed)
%     SN        – (Possibly updated) network structure.
%     STARTN    – M×K rate at which a class-r service starts at station i
%     PREEMPTN  – M×K rate at which a class-r job in service is displaced
%
%   See also SOLVER_SSA_NRM_SPACE, NEXT_REACTION_METHOD_DIRECT.

% ---------------------------------------------------------------------
% Parameters & shorthands
% ---------------------------------------------------------------------
samples = options.samples;
R  = sn.nclasses;
I  = sn.nnodes;
M  = sn.nstations;
K  = sn.nclasses;
state = sn.state;

% ---------------------------------------------------------------------
% Phase slot map --------------------------------------------------------
% ---------------------------------------------------------------------
% The state vector counts jobs per (node, class, PHASE) rather than per
% (node, class), so that phase-type service is represented exactly instead of
% being collapsed onto its mean rate. Phase counts differ per (station, class)
% via sn.phasessz, hence the explicit offset map rather than arithmetic on R.
%
% The layout is chosen so that a single-phase model is bit-for-bit the old one:
% with nph == 1 everywhere, phOff(ind,r) = (ind-1)*R + (r-1) and therefore
% slot(ind,r,1) = (ind-1)*R + r, exactly the flat class index the engine used
% before. Every exponential model must reproduce its previous results, which is
% the self-check for this generalization.
nph = ones(I, R);
for ind = 1:I
    if sn.isstation(ind)
        ist = sn.nodeToStation(ind);
        for r = 1:R
            nph(ind,r) = max(1, sn.phasessz(ist,r));
        end
    end
end
phOff = zeros(I, R);
NS = 0;
for ind = 1:I
    for r = 1:R
        phOff(ind,r) = NS;
        NS = NS + nph(ind,r);
    end
end
isPhaseExpanded = NS > I*R;   % false for a purely exponential model
% Reverse map. Every consumer that used to decode a state index arithmetically
% (floor((slot-1)/R)+1, mod(slot-1,R)+1) must go through this instead: with
% unequal phase counts the flat arithmetic no longer identifies the node.
smap = struct();   % layout descriptor threaded into every slot consumer
slotNode = zeros(NS,1);
slotClass = zeros(NS,1);
slotPhase = zeros(NS,1);
for ind = 1:I
    for r = 1:R
        for kk = 1:nph(ind,r)
            slotNode(phOff(ind,r)+kk) = ind;
            slotClass(phOff(ind,r)+kk) = r;
            slotPhase(phOff(ind,r)+kk) = kk;
        end
    end
end
smap.node = slotNode; smap.class = slotClass; smap.phase = slotPhase;
smap.phOff = phOff; smap.nph = nph; smap.R = R;

% Buffered phase-type service: see _kb/06-solver-catalog.md (SSA/NRM section)
% for the svcph auxiliary-structure rationale.
bufPHSched = [SchedStrategy.FCFS, SchedStrategy.LCFS, SchedStrategy.SIRO, ...
    SchedStrategy.HOL, SchedStrategy.SEPT, SchedStrategy.LEPT];
bufPHClass = false(I, R);
for ind = 1:I
    if sn.isstation(ind) && any(sn.sched(sn.nodeToStation(ind)) == bufPHSched)
        for r = 1:R
            if nph(ind,r) > 1
                bufPHClass(ind,r) = true;
            end
        end
    end
end
bufPHNode = any(bufPHClass, 2);
maxnph = max(nph(:));
smap.bufPHClass = bufPHClass; smap.bufPHNode = bufPHNode;

% Cache nodes. A Cache is an immediate class-switch: a job arrives in a READ
% class, reads an item drawn from pread, and leaves in the hit or miss class
% depending on whether the item is currently cached, after which the replacement
% policy updates the cache contents. The NRM models this as a state-dependent
% class-switch reaction at the cache node (consume [cache,readClass], produce
% [cache,hitClass] or [cache,missClass], chosen at firing by the cache access),
% mirroring State.afterEventCache. The cache contents ride alongside the buffers
% (buffers{cacheNode}) since no rate depends on them; the hit/miss draw and the
% replacement update read and rewrite them at firing. A class r is a READ class
% of cache ind iff its pread entry is a non-empty probability row.
isCacheNode = false(I,1);
isCacheReadClass = false(I,R);
for ind = 1:I
    if sn.nodetype(ind) == NodeType.Cache
        isCacheNode(ind) = true;
        np = sn.nodeparam{ind};
        if isfield(np,'retrievalClassIndices') && ~isempty(np.retrievalClassIndices)
            rci = np.retrievalClassIndices(:)';
        else
            rci = [];
        end
        for r = 1:R
            % A normal read class has a hit class; a retrieval class (created by
            % setRetrievalSystem) reads its own one-hot item to COMPLETE a miss
            % and has hitclass == 0 but a miss class. Both take a cache-access
            % reaction; the outcome class is resolved at firing.
            if r <= numel(np.pread) && ~isempty(np.pread{r}) && all(~isnan(np.pread{r}(:))) ...
                    && ((r <= numel(np.hitclass) && np.hitclass(r) > 0) || any(rci == r))
                isCacheReadClass(ind,r) = true;
            end
        end
    end
end
smap.isCacheNode = isCacheNode;

% Destination slot a BEGUN retrieval routes to (the fetch queue). When a miss
% starts a retrieval the job must go to the retrieval queue, be served (the fetch
% delay), and only THEN return to the cache to complete the miss. If it were left
% at [cache,retrievalClass] the cache-access reaction would fire again and
% complete the miss instantly, collapsing the fetch. So a begin lands the job at
% the retrieval class's routed destination; only queue-returns occupy
% [cache,retrievalClass] and trigger completion. Resolved once from rtnodes.
cacheRetrDest = zeros(I, R);
for ind = 1:I
    if isCacheNode(ind)
        np = sn.nodeparam{ind};
        if isfield(np,'retrievalClassIndices') && ~isempty(np.retrievalClassIndices)
            for rc = np.retrievalClassIndices(:)'
                row = sn.rtnodes((ind-1)*R + rc, :);
                dslot = find(row > 0, 1, 'first');
                if ~isempty(dslot)
                    jnd = floor((dslot-1)/R) + 1; s = mod(dslot-1,R) + 1;
                    cacheRetrDest(ind, rc) = phOff(jnd, s) + 1;
                end
            end
        end
    end
end

% ---------------------------------------------------------------------
% Stochastic Petri net path --------------------------------------------
% ---------------------------------------------------------------------
% A model with Transition nodes is a stochastic Petri net, not a queueing
% network: its dynamics are firings of transition modes over a place marking,
% not job departures routed by the rt matrix. The stoichiometry matrix of the
% reaction network IS the net's incidence matrix, so the NRM is the natural
% simulator, but the generic (node,class) departure grid below does not apply
% (a firing produces to several places deterministically, never a routing
% draw). Route Petri nets to the dedicated builder/runner, which shares the
% Gibson & Bruck clocks but its own firing application and vanishing-marking
% collapse for immediate transitions.
if any(sn.nodetype == NodeType.Transition)
    [QN, UN, RN, TN, CN, XN] = solver_ssa_nrm_spn(sn, options, phOff, nph, NS, smap);
    lG = 0;
    % A Petri net holds tokens in places and fires transitions: there is no
    % service facility to seize or be displaced from, so both derived rates
    % are structurally zero rather than unmeasured.
    StartN = zeros(sn.nstations, sn.nclasses);
    PreemptN = zeros(sn.nstations, sn.nclasses);
    return
end

% ---------------------------------------------------------------------
% Stoichiometry & reaction mapping (self‑loops included) ----------------
% ---------------------------------------------------------------------
S       = zeros(0, NS);   % will transpose at the end
fromIdx = [];
toIdx   = [];
fromIR  = [];

% this is currently M^2*R^2, it can be lowered to M*R decoupling the
% routing
% Departure reactions, one per (node, class, PHASE). A departure is the
% absorption of the phase-type service process, so it fires at mu(k)*phi(k) and
% the job re-enters its destination in an entry phase drawn from pie: the
% destination draw therefore carries the product of the routing probability and
% the entry-phase probability. The existing weighted-destination sampler takes
% that product unchanged.
k = 0;
depPhase = [];       % service phase each departure reaction consumes
isBufSvcRx = false(0,1);   % departure reaction of a buffered-PH class (reads svcph)
isCacheRx = false(0,1);    % cache-access reaction (read -> hit/miss at a cache node)
cacheHitSlot = zeros(0,1);  % nvec slot the hit-class job is produced into
cacheMissSlot = zeros(0,1); % nvec slot the miss-class job is produced into
for ind = 1:I
    for r = 1:R
        for kk = 1:nph(ind,r)
            k = k + 1;
            fromIR(k,:) = [ind, r];
            depPhase(k,1) = kk;
            isCacheRx(k,1) = false;
            cacheHitSlot(k,1) = 0;
            cacheMissSlot(k,1) = 0;
            if isCacheReadClass(ind,r)
                % Cache access: consume the read-class job at the cache; its
                % production (hit or miss class, at the SAME cache node) and the
                % contents update are resolved at firing by cacheAccess. No
                % static routing: rtnodes has no out-edge for the read class.
                fromIdx(k) = phOff(ind,r) + 1;
                np = sn.nodeparam{ind};
                Srow = zeros(1, NS);
                Srow(fromIdx(k)) = -1;
                S(k,:) = Srow;
                probIR{k} = [];
                toIdx{k} = [];
                isCacheRx(k,1) = true;
                isBufSvcRx(k,1) = false;
                % The outcome class (hit/miss/retrieval) is resolved at firing, so
                % these slots are informational only; a retrieval class has
                % hitclass 0, so guard the lookup.
                if r <= numel(np.hitclass) && np.hitclass(r) > 0
                    cacheHitSlot(k,1) = phOff(ind, np.hitclass(r)) + 1;
                end
                if r <= numel(np.missclass) && np.missclass(r) > 0
                    cacheMissSlot(k,1) = phOff(ind, np.missclass(r)) + 1;
                end
                continue
            end
            % At a buffered-PH source only the jobs in service carry a phase and
            % the phase composition lives in svcph, not nvec; nvec holds the whole
            % class population in its first phase slot. A departure therefore
            % removes one job from that total slot regardless of which service
            % phase completed -- the completing phase kk is carried in depPhase
            % and consumed from svcph at firing.
            if bufPHClass(ind,r)
                fromIdx(k) = phOff(ind,r) + 1;
                isBufSvcRx(k,1) = true;
            else
                fromIdx(k) = phOff(ind,r) + kk;
                isBufSvcRx(k,1) = false;
            end
            probIR{k} = [];
            toIdx{k} = [];
            Srow = zeros(1, NS); % build stoichiometry row
            if sn.isslc(r)
                Srow(fromIdx(k))   = -Inf;
            else
                Srow(fromIdx(k)) = -1;
                for jnd = 1:I
                    for s = 1:R
                        p = sn.rtnodes((ind-1)*R+r, (jnd-1)*R+s);
                        if p > 0
                            if bufPHClass(jnd,s)
                                % A job arriving at a buffered-PH destination lands
                                % in the total-population slot; whether it enters
                                % service (and in which entry phase) or waits is
                                % decided at firing from the server occupancy and
                                % pie, not by the routing draw. So the destination
                                % collapses to the single total slot with weight p.
                                dslot = phOff(jnd,s) + 1;
                                toIdx{k}(end+1) = dslot;
                                probIR{k}(end+1) = p;
                                Srow(dslot) = Srow(dslot) + p;
                            else
                                pentry = entryProbs(sn, jnd, s, nph(jnd,s));
                                for ke = 1:nph(jnd,s)
                                    if pentry(ke) <= 0
                                        continue
                                    end
                                    dslot = phOff(jnd,s) + ke;
                                    toIdx{k}(end+1) = dslot;
                                    probIR{k}(end+1) = p * pentry(ke);
                                    Srow(dslot) = Srow(dslot) + p * pentry(ke);
                                end
                            end
                        end
                    end
                end
            end
            S(k,:) = Srow;
        end
    end
end
nDepRx = k;   % departure reactions occupy 1..nDepRx
isBufSvcRx(end+1:k,1) = false;

% Phase-transition reactions, one per (node, class, k -> k'). These move a job
% between the phases of its own service process and so never leave the node;
% D0's off-diagonal carries their rates (State.afterEventStation, EventType.PHASE).
isPhaseRx = false(k,1);
phaseFrom = zeros(k,1);
phaseTo = zeros(k,1);
phaseRate = zeros(k,1);
for ind = 1:I
    if ~sn.isstation(ind)
        continue
    end
    ist = sn.nodeToStation(ind);
    for r = 1:R
        if nph(ind,r) <= 1 || isempty(sn.proc{ist}{r})
            continue
        end
        D0 = sn.proc{ist}{r}{1};
        for ka = 1:nph(ind,r)
            for kb = 1:nph(ind,r)
                if ka == kb || D0(ka,kb) <= 0
                    continue
                end
                k = k + 1;
                fromIR(k,:) = [ind, r];
                probIR{k} = [];
                toIdx{k} = [];
                Srow = zeros(1, NS);
                if bufPHClass(ind,r)
                    % A buffered-PH class keeps its in-service phase counts in
                    % svcph, not in nvec: a phase transition moves a job between
                    % phases of the SAME in-service composition, so it leaves nvec
                    % (the class total) unchanged. The stoichiometry column is
                    % therefore all zeros; the move is applied to svcph at firing
                    % and, like a retry/switchover, its dependency set must be
                    % supplied through a forced refresh (D cannot derive it from S).
                    fromIdx(k) = phOff(ind,r) + 1;
                else
                    fromIdx(k) = phOff(ind,r) + ka;
                    Srow(phOff(ind,r) + ka) = -1;
                    Srow(phOff(ind,r) + kb) = 1;
                end
                S(k,:) = Srow;
                isPhaseRx(k,1) = true;
                phaseFrom(k,1) = ka;
                phaseTo(k,1) = kb;
                phaseRate(k,1) = D0(ka,kb);
            end
        end
    end
end
isPhaseRx(end+1:k,1) = false;
depPhase(end+1:k,1) = 0;

% Reneging: each waiting (queued, not-in-service) class-r job abandons at the
% memoryless rate sn.impatienceMu, so the aggregate rate out of the state is
% (waiting count)*mu and the job leaves the system (the passive half of the
% sync is LOCAL in refreshSync). This is a reaction the (node,class) departure
% grid above cannot express -- it consumes a job without producing one -- so it
% is appended as an extra column whose stoichiometry is a bare -1 at the source
% slot. A renege is not a departure and must not count towards throughput; the
% TN accumulator reads the first reaction with a given source slot, which is
% always the departure, so the appended columns stay out of it.
nDep = k;                       % departure reactions occupy 1..nDep
isRenegeRx = false(nDep,1);
renegeMu = zeros(nDep,1);
if isfield(sn,'impatienceClass') && ~isempty(sn.impatienceClass) ...
        && any(sn.impatienceClass(:) == ImpatienceType.RENEGING)
    for ist = 1:M
        ind = sn.stationToNode(ist);
        for r = 1:R
            if sn.impatienceClass(ist,r) == ImpatienceType.RENEGING && sn.impatienceMu(ist,r) > 0
                k = k + 1;
                fromIR(k,:) = [ind, r];
                fromIdx(k) = (ind-1)*R + r;
                probIR{k} = [];
                toIdx{k} = [];
                Srow = zeros(1, I*R);
                Srow((ind-1)*R + r) = -1;   % job abandons and leaves the system
                S(k,:) = Srow;
                isRenegeRx(k,1) = true;
                renegeMu(k,1) = sn.impatienceMu(ist,r);
            end
        end
    end
end
% Retrial: an orbiting class-r job retries entry at the memoryless rate
% sn.retrialMu and succeeds only when a server is free; otherwise the event is
% a no-op and is simply not generated (State.afterEventStation, EventType.RETRY).
% The orbit needs no new state: orbiting jobs are already counted in the
% station population and held in the buffer, so orbit_r is exactly the buffer
% occupancy the FCFS-family rate law already reads. A retry moves a job from
% the orbit into service WITHOUT changing any population, so its stoichiometry
% column is all zeros -- which is why its dependency set has to be supplied by
% hand below: D is derived from S, and an all-zero column would otherwise leave
% every rate at the node stale after a retry fires.
isRetryRx = false(k,1);
retryMu = zeros(k,1);
retryNode = zeros(k,1);
if isfield(sn,'retrialProc') && ~isempty(sn.retrialProc)
    for ist = 1:M
        if ~any(~cellfun(@isempty, sn.retrialProc(ist,:)))
            continue
        end
        ind = sn.stationToNode(ist);
        for r = 1:R
            if sn.retrialMu(ist,r) > 0
                k = k + 1;
                fromIR(k,:) = [ind, r];
                fromIdx(k) = (ind-1)*R + r;
                probIR{k} = [];
                toIdx{k} = [];
                S(k,:) = zeros(1, I*R);   % a retry moves no job between nodes
                isRetryRx(k,1) = true;
                retryMu(k,1) = sn.retrialMu(ist,r);
                retryNode(k,1) = ind;
            end
        end
    end
end

% Polling switchover reactions. A polling server cycles through the buffers it
% serves, carrying a controller [mode, pos, swphase, ctr] in the auxiliary
% buffer of its node (mode 0 parked, 1 serving pos, 2 switching towards pos).
% A service departure fires only while the controller serves that class (gated
% in the propensity below); when a visit ends the server walks the cyclic order
% (State.pollingNext folds every immediate switchover) and, on meeting a timed
% switchover, dwells in mode 2. That dwell is a genuine timed event with no job
% movement, so it is appended here as one reaction per polling node with a timed
% switchover, exactly as a retry is: an all-zero stoichiometry column whose
% propensity reads the controller and whose firing samples the switchover PH.
% Only exponential service is expanded at a polling station (phaseNrmOK gates PH
% service there); the switchover itself may be phase-type, its phases carried in
% the controller rather than in nvec.
poll = struct('on', false);
poll.isPoll = false(1, I);
poll.pinfo = cell(I, 1);
poll.swRx = zeros(1, I);           % switchover reaction index of each polling node, 0 if none
for ind = 1:I
    if sn.isstation(ind) && sn.sched(sn.nodeToStation(ind)) == SchedStrategy.POLLING
        poll.pinfo{ind} = State.pollingInfo(sn, ind);
        poll.isPoll(ind) = true;
        poll.on = true;
        % The NRM tracks only the controller of a polling station, not the
        % service phase of the single job in service, so phase-type service at a
        % polling station is not expanded here. Reject it rather than spread the
        % class over phases and gate each phase reaction on the same class (which
        % would serve several fictitious phase-jobs at once). Switchover may be
        % phase-type: its phase is carried in the controller.
        for r = 1:R
            if nph(ind,r) > 1
                line_error(mfilename, sprintf('NRM polling supports exponential service only; station %d class %d has phase-type service. Use method=''serial''.', sn.nodeToStation(ind), r));
            end
        end
    end
end
isPollSwRx = false(k, 1);
pollSwNode = zeros(k, 1);
if poll.on
    for ind = 1:I
        pinf = poll.pinfo{ind};
        if isempty(pinf) || ~any(pinf.hasSw)
            continue % no timed switchover: the server never dwells in a walk
        end
        k = k + 1;
        fromIR(k,:) = [ind, 1];      % class field is a sentinel; never read as a class here
        fromIdx(k) = (ind-1)*R + 1;  % unused slot: a switchover consumes no job
        probIR{k} = [];
        toIdx{k} = [];
        S(k,:) = zeros(1, I*R);      % a switchover moves no job between nodes
        isPollSwRx(k,1) = true;
        pollSwNode(k,1) = ind;
        poll.swRx(ind) = k;
    end
end

% Pad every per-reaction marker to the final reaction count, so the reaction
% loops below index them safely regardless of which extra-reaction families
% (renege, retry, switchover) are present.
isRenegeRx(end+1:k,1) = false;
isRetryRx(end+1:k,1) = false;
isPhaseRx(end+1:k,1) = false;
depPhase(end+1:k,1) = 0;
isPollSwRx(end+1:k,1) = false;
pollSwNode(end+1:k,1) = 0;
isBufSvcRx(end+1:k,1) = false;
isCacheRx(end+1:k,1) = false;
cacheHitSlot(end+1:k,1) = 0;
cacheMissSlot(end+1:k,1) = 0;
renegeMu(end+1:k,1) = 0;
retryMu(end+1:k,1) = 0;
retryNode(end+1:k,1) = 0;

S = S.';   % states × reactions

% ---------------------------------------------------------------------
% Initial state vector --------------------------------------------------
% ---------------------------------------------------------------------
nvec0 = zeros(NS,1); % initial state (per node, class and phase)
% Non-preemptive policies that hold waiting jobs in a buffer. They share the
% rate law (a class-r completion fires at mu_r times the class-r jobs actually
% in service) and differ only in which waiting job is promoted on a departure;
% see pickFromBuffer.
bufferedSched = [SchedStrategy.FCFS, SchedStrategy.LCFS, SchedStrategy.SIRO, ...
    SchedStrategy.HOL, SchedStrategy.SEPT, SchedStrategy.LEPT, ...
    SchedStrategy.LCFSPR, SchedStrategy.PAS];
% Preemptive policies: an arrival at a fully busy station takes a server and
% pushes the incumbent it displaced back into the buffer, rather than queueing
% itself (State.afterEventStation, the FCFSPR/LCFSPR arrival group). The rate
% law is unchanged -- it still counts the jobs actually in service -- so no
% extra state is needed: in-service is population minus buffer occupancy, and
% that automatically names the new arrival as the one being served. With
% exponential service preempt-resume needs no stored phase, because a resumed
% job has the same memoryless residual as a fresh one. LCFSPI is deliberately
% absent: SolverSSA.getFeatureSet does not advertise it (nor does SolverCTMC),
% so the NRM must not claim it either.
preemptiveSched = SchedStrategy.LCFSPR;
% Order-independent / pass-and-swap stations keep the FULL ordered job list,
% not just the waiting jobs: there is no server/buffer split at all, and the
% rate is a function mu(c) of the whole list (State.afterEventStationPAS). So
% these carry a different buffer invariant -- numel(buf) == total, rather than
% max(0, total - mi) -- and the list runs OLDEST-FIRST (c(1) is the oldest),
% the reverse of every other buffered policy here.
% sn.sched carries PAS for both PAS and OI stations: OI is canonicalized to
% pass-and-swap with an all-zero swap graph (see MNetwork.refreshLocalVars), so
% reading the graph covers both and OI needs no separate case.
listSched = SchedStrategy.PAS;
% Of those, the policies whose state keeps the buffer as per-class counts
% rather than as an ordered list of class ids (State.fromMarginalAndRunning).
countBufferedSched = [SchedStrategy.SIRO, SchedStrategy.SEPT, SchedStrategy.LEPT];
buffers0 = cell(I,1); % per-node ordered buffer of waiting job classes (FCFS/LCFS)
for ind=1:I
    buffers0{ind} = [];
end
% A cache node has no queueing buffer; its buffers slot instead carries the cache
% CONTENTS (the item held in each of the totalCacheCapacity slots, ordered by
% list as State.afterEventCache lays them out). Any valid ordered placement is a
% correct warm start since the chain is ergodic, so slot i starts holding item i.
for ind=1:I
    if isCacheNode(ind)
        np = sn.nodeparam{ind};
        if isfield(np,'totalCacheCapacity') && ~isempty(np.totalCacheCapacity)
            tcc = np.totalCacheCapacity;
        else
            tcc = sum(np.itemcap);
        end
        % With a retrieval system the contents are followed by block A, a per-item
        % occupancy bitmap (column tcc+i is 1 iff item i is currently being
        % retrieved), and by block B, the per-retrieval-class count of requests
        % merged onto an in-flight fetch (State.spaceCache). Both start empty.
        if isfield(np,'retrievalSystemCapacity') && ~isempty(np.retrievalSystemCapacity) ...
                && any(np.retrievalSystemCapacity > 0)
            rcListInit = State.cacheRetrievalClassMap(sn, ind);
            buffers0{ind} = [1:tcc, zeros(1, np.nitems + numel(rcListInit))];
        else
            buffers0{ind} = 1:tcc;
        end
    end
end
% In-service phase multiset of each buffered-PH node: svcph0{ind}(r,k) counts the
% class-r jobs in service in phase k. Empty for every other node. Populated below
% once the waiting buffer of each buffered-PH node is known (in-service = class
% total minus waiting), so it is filled after the buffer loop.
svcph0 = cell(I,1);
for ind=1:I
    if bufPHNode(ind)
        svcph0{ind} = zeros(R, maxnph);
    else
        svcph0{ind} = [];
    end
end
for ind=1:I
    if sn.isstateful(ind)
        state_i = state{sn.nodeToStateful(ind)};
        [~,nir] = State.toMarginalAggr(sn, ind, state_i);
        for r = 1:R
            if isinf(nir(r))
                if sn.nodetype(ind) == NodeType.Source
                    nir(r) = 1;
                else
                    line_error(mfilename, 'Infinite population error.');
                end
            end
            % Spread the class population across its phases. The marginal the
            % initial state carries is per class, not per phase, so the entry
            % distribution pie is the natural allocation: it is the phase a job
            % starts service in. For a single-phase class this puts everything
            % in slot 1, reproducing the old flat layout exactly. A buffered-PH
            % class keeps its whole population in slot 1 too -- nvec is the class
            % total there and the in-service phase composition lives in svcph0
            % (built below), so the phase slots 2..nph stay empty in nvec.
            if nph(ind,r) <= 1 || bufPHClass(ind,r)
                nvec0(phOff(ind,r) + 1,1) = nir(r);
            else
                pe = entryProbs(sn, ind, r, nph(ind,r));
                left = nir(r);
                for ke = 1:nph(ind,r)
                    if ke == nph(ind,r)
                        take = left;
                    else
                        take = min(left, round(nir(r) * pe(ke)));
                    end
                    nvec0(phOff(ind,r) + ke,1) = take;
                    left = left - take;
                end
            end
        end

        % Populate buffers for buffered nodes from the raw state vector
        % (only stations have buffered scheduling; skip non-station
        % stateful nodes such as RROBIN dispatchers/Routers and Caches)
        ist = sn.nodeToStation(ind);
        if ist >= 1 && any(sn.sched(ist) == bufferedSched)
            sumK = sum(sn.phasessz(ist,:));
            sumNvars = sum(sn.nvars(ind,:));
            bufCols = size(state_i,2) - sumK - sumNvars;
            if any(sn.sched(ist) == listSched)
                % PAS/OI stores the full ordered list left-aligned in the first
                % cap(ist) columns, c(1) oldest, zero-padded on the right --
                % already the order the NRM needs, so it is copied verbatim
                % rather than reversed.
                % Its width is nCols - nvars, NOT the shared bufCols: a PAS
                % station has no server/phase block at all (there is no
                % server/buffer split), yet phasessz still floors to 1 per class
                % as for any other station, so subtracting sumK here would drop
                % the last sum(phasessz) entries of the list. Both PAS
                % authorities, State.afterEventStationPAS and the PAS branch of
                % State.toMarginal, read W = size(inspace,2) - V.
                pasCols = size(state_i,2) - sumNvars;
                for pos = 1:pasCols
                    classId = state_i(1,pos);
                    if classId >= 1 && classId <= R
                        buffers0{ind}(end+1) = classId;
                    end
                end
            elseif any(sn.sched(ist) == countBufferedSched)
                % SIRO/SEPT/LEPT keep an UN-ordered buffer: the first R columns
                % hold the per-class counts of waiting jobs, not class ids (see
                % State.fromMarginalAndRunning). Expand them into the NRM's
                % ordered list; the order within it is immaterial for these
                % disciplines, which select by class and never by position.
                for r = 1:min(R, bufCols)
                    buffers0{ind}(end+1:end+state_i(1,r)) = r;
                end
            else
                % FCFS/HOL/LCFS keep an ordered list of class ids
                for pos = 1:bufCols
                    classId = state_i(1,pos);
                    if classId >= 1 && classId <= R
                        buffers0{ind}(end+1) = classId; % addLast
                    end
                    % classId == 0 means empty position, skip
                end
            end
        end

        % Seed the in-service phase multiset of a buffered-PH node. The jobs in
        % service are the class total minus the ones waiting in the buffer just
        % built; their starting phases are drawn from the entry distribution pie,
        % the same allocation the INF/PS init uses. Only in-service jobs get a
        % phase -- waiting jobs have not started service and carry none.
        if bufPHNode(ind)
            for r = 1:R
                waiting_r = sum(buffers0{ind} == r);
                insvc_r = max(0, nir(r) - waiting_r);
                if nph(ind,r) <= 1
                    svcph0{ind}(r,1) = insvc_r;
                else
                    pe = entryProbs(sn, ind, r, nph(ind,r));
                    left = insvc_r;
                    for ke = 1:nph(ind,r)
                        if ke == nph(ind,r)
                            take = left;
                        else
                            take = min(left, round(insvc_r * pe(ke)));
                        end
                        svcph0{ind}(r,ke) = take;
                        left = left - take;
                    end
                end
            end
        end
    end
end

mi    = zeros(I,1);
rates = zeros(I,R);
for ind=1:I
    if sn.isstation(ind)
        for r=1:R
            ist = sn.nodeToStation(ind);
            muir = sn.rates(ist,r);
            if ~isnan(muir)
                rates(ind,r) = muir;
            end
            mi(ind,1) = sn.nservers(ist);
        end
    else
        for r=1:R
            rates(ind,r) = GlobalConstants.Immediate;
            mi(ind,1) = GlobalConstants.MaxInt;
        end
    end
    mi(isinf(mi)) = GlobalConstants.MaxInt;
end

% Limited load-dependent scaling lld(ist, ntot): a work-conserving factor that
% multiplies the aggregate service rate at total station population ntot (as in
% State.afterEventStation). Default (all ones) for stations without load
% dependence, so it is inert for plain single-/multi-server queues.
if isempty(sn.lldscaling)
    lldMat = []; lldlimit = 0;
else
    lldMat = sn.lldscaling; lldlimit = size(lldMat,2);
end

% Class-dependent scaling cdscaling{ist}: a handle mapping the per-class
% station population vector n to the 1xR vector of rate scalings beta_r(n)
% (as in State.afterEventStation, evaluated per firing on the current state).
% Joint-dependence handles eta_i(n) (sn.jdscaling, non-product-form) enter the
% sample-path rates the same way, so fold them into the effective per-station
% handle eta_i(n).*beta_{i,r}(n), exactly as State.afterEventInit does.
if isempty(sn.cdscaling)
    cdCell = {};
else
    cdCell = sn.cdscaling;
end
if ~isempty(sn.jdscaling)
    M_ = sn.nstations;
    if isempty(cdCell)
        cdCell = cell(M_,1);
    end
    for ist = 1:M_
        if ist <= numel(sn.jdscaling) && ~isempty(sn.jdscaling{ist})
            jdh = sn.jdscaling{ist};
            if ist <= numel(cdCell) && ~isempty(cdCell{ist})
                cdh = cdCell{ist};
                cdCell{ist} = @(ni) cdh(ni) .* jdh(ni);
            else
                cdCell{ist} = jdh;
            end
        end
    end
end

% Scheduling policies whose rate law reads per-class weights from
% sn.schedparam. These are single-server only, as in State.afterEventStation.
weightedSched = [SchedStrategy.DPS, SchedStrategy.GPS, ...
    SchedStrategy.DPSPRIO, SchedStrategy.GPSPRIO];

% Finite capacity regions (DROP rule). A region constrains an aggregate of the
% per-class populations of its member stations, which is a linear function of
% the NRM state vector, so admission is a multiplicative 0/1 gate on the
% routing draw. The DROP rule censors the refused transition, and censoring an
% exponential transition is exactly what zeroing its share of the propensity
% does. WAITQ instead parks refused jobs in a per-region FIFO, which is extra
% state the reaction network does not carry, so those models are routed to the
% serial engine by SOLVER_SSA_ANALYZER and never reach here.
fcr = fcrPrecompute(sn);

% Balking. An arrival that balks is lost: it has left its source but never
% joins the destination, so the departure rate is unchanged and only the
% arrival outcome differs (State.afterEventStation scales the admitted
% branches by 1-balkProb and adds a balked branch of probability balkProb that
% leaves the destination state untouched). Only the QUEUE_LENGTH strategy is a
% pure function of the state vector; EXPECTED_WAIT / COMBINED depend on the
% mean wait and are rejected by the analyzers.
balk = balkPrecompute(sn);

% G-network signals. A signal class never joins the station it reaches: it acts
% on the jobs already there and is annihilated (State.afterEventStationSignal).
% That makes it an arrival-side effect exactly like balking, so the departure
% rate is unchanged and only the arrival outcome differs. The reference
% enumerates every victim subset with its probability because it builds a
% generator; a simulator instead draws the batch size and the victims, which is
% equivalent and avoids the enumeration.
sig = signalPrecompute(sn);

% Round-robin routing. The pointer that RROBIN/WRROBIN walk is a per-(node,
% class) local variable, not a population, and no rate depends on it: it only
% decides where a departure goes. In a generator that makes it a genuine extra
% state dimension, but a simulator can carry it as auxiliary state alongside
% the buffers, which is what happens here. State.afterEventRouter advances the
% pointer on the departure and the routing closure then reads state_AFTER, so
% the destination used is the one the pointer lands on -- advance first, then
% select.
rr = rrPrecompute(sn, state);

% Propensity function ---------------------------------------------------
epstol = GlobalConstants.Zero;
a = {};
classprio = sn.classprio(:)';          % lower value = higher priority in LINE
% Rate of the service-process event each reaction carries: the absorption
% mu(k)*phi(k) for a departure, the off-diagonal D0(k,k') for a phase change.
% For a single-phase class this is just the exponential rate, so an exponential
% model sees exactly the rates it saw before.
rateOf = zeros(length(fromIdx),1);
for j = 1:length(fromIdx)
    ind = fromIR(j,1); r = fromIR(j,2);
    if isPhaseRx(j)
        rateOf(j) = phaseRate(j);
    elseif sn.isstation(ind)
        ist = sn.nodeToStation(ind);
        kk = depPhase(j);
        if nph(ind,r) > 1 && ~isempty(sn.proc{ist}{r})
            rateOf(j) = sn.mu{ist}{r}(kk) * sn.phi{ist}{r}(kk);
        else
            rateOf(j) = rates(ind, r);
        end
    else
        rateOf(j) = rates(ind, r);
    end
end

for j=1:length(fromIdx)
    ind = fromIR(j,1);
    base = (ind-1)*R + 1;              % first per-class state slot of this node
    ldrow = [];                        % load-dependent scaling row of the station
    cdbeta = [];                       % class-dependence handle of the station
    wrow = [];                         % normalized DPS/GPS scheduling weights
    if sn.isstation(ind)
        istj = sn.nodeToStation(ind);
        if istj >= 1 && ~isempty(lldMat), ldrow = lldMat(istj, :); end
        if istj >= 1 && istj <= numel(cdCell) && ~isempty(cdCell{istj})
            cdbeta = cdCell{istj};
        end
        if istj >= 1 && any(sn.sched(istj) == weightedSched)
            wrow = sn.schedparam(istj, 1:R);
            if sum(wrow) <= 0
                line_error(mfilename, sprintf('Station %d has %s scheduling with non-positive total weight.', istj, SchedStrategy.toText(sn.sched(istj))));
            end
            wrow = wrow / sum(wrow);
            % State.afterEventStation rejects multi-server DPS/GPS, so the
            % rate law below is only defined for a single server. Fail here
            % rather than silently simulate a different station.
            if mi(ind) > 1
                line_error(mfilename, sprintf('Multi-server %s stations are not supported yet.', SchedStrategy.toText(sn.sched(istj))));
            end
        end
    end
    % Buffered phase-type service. Only the jobs in service carry a phase, and
    % their per-phase counts live in svc{ind}(r,k), not in nvec. Both the
    % departure (absorption of phase kk) and the internal phase transition
    % (kk -> kb) therefore fire at rate rateOf(j) times the number of class-r
    % jobs currently in service in the source phase kk -- exactly the INF-family
    % law rateOf*kir, but with kir read from the in-service multiset svc rather
    % than from nvec (whose class total also counts the waiting jobs). The
    % load-/class-dependent factors still read the total population, as for the
    % exponential buffered law. A single-phase (exponential) buffered class is
    % NOT bufPHClass and keeps its original rate law below.
    if sn.isstation(ind) && bufPHClass(ind, fromIR(j,2))
        rr_ph = fromIR(j,2);
        if isPhaseRx(j)
            kk_ph = phaseFrom(j);
        else
            kk_ph = depPhase(j);
        end
        a{j} = @(X, bufs, svc) rateOf(j) * svc{ind}(rr_ph, kk_ph) ...
            * lldfac(ldrow, sum(classCounts(X, phOff, nph, ind, R)), lldlimit) ...
            * cdfac(cdbeta, classCounts(X, phOff, nph, ind, R), fromIR(j,2));
        continue
    end
    if sn.isstation(ind)
        switch sn.sched(sn.nodeToStation(ind))
            case SchedStrategy.EXT
                % A Source fires at a constant arrival rate. It has no service
                % phases (nph == 1 there, enforced by phaseNrmOK), so the
                % kir/nir share must NOT be applied: the Source's fictitious
                % token would drive kirFrac to 0 and silence the Source, which
                % deadlocks every open model. rateOf(j) is that constant rate.
                a{j} = @(X, bufs, svc) rateOf(j);
            case SchedStrategy.INF
                a{j} = @(X, bufs, svc) rateOf(j) * kirFrac(X, fromIdx(j), phOff, nph, fromIR(j,1), fromIR(j,2)) * classPop(X, phOff, nph, ind, fromIR(j,2)) ...
                    * cdfac(cdbeta, classCounts(X, phOff, nph, ind, R), fromIR(j,2));
            case {SchedStrategy.PS, SchedStrategy.LPS}
                % LPS shares the PS rate law in State.afterEventStation: the
                % sharing limit is the server count, so min(ni,c) covers both.
                if R == 1 % single class
                    a{j} = @(X, bufs, svc) rateOf(j) * kirFrac(X, fromIdx(j), phOff, nph, fromIR(j,1), fromIR(j,2)) * min( mi(fromIR(j,1)), classPop(X, phOff, nph, ind, fromIR(j,2))) ...
                        * lldfac(ldrow, classPop(X, phOff, nph, ind, fromIR(j,2)), lldlimit) ...
                        * cdfac(cdbeta, classCounts(X, phOff, nph, ind, R), fromIR(j,2));
                else
                    a{j} = @(X, bufs, svc) rateOf(j) * kirFrac(X, fromIdx(j), phOff, nph, fromIR(j,1), fromIR(j,2)) * ( classPop(X, phOff, nph, ind, fromIR(j,2)) ./ ...
                        (epstol+sum( classCounts(X, phOff, nph, ind, R) ) )) * ...
                        min( mi(fromIR(j,1)), (epstol+sum( classCounts(X, phOff, nph, ind, R) )) ) ...
                        * lldfac(ldrow, sum(classCounts(X, phOff, nph, ind, R)), lldlimit) ...
                        * cdfac(cdbeta, classCounts(X, phOff, nph, ind, R), fromIR(j,2));
                end
            case SchedStrategy.DPS
                % Discriminatory PS: class r receives a share w_r*n_r/(w.n) of
                % the single server (State.afterEventStation, case DPS).
                a{j} = @(X, bufs, svc) rateOf(j) * kirFrac(X, fromIdx(j), phOff, nph, fromIR(j,1), fromIR(j,2)) ...
                    * dpsshare(wrow, classCounts(X, phOff, nph, ind, R), fromIR(j,2)) ...
                    * lldfac(ldrow, sum(classCounts(X, phOff, nph, ind, R)), lldlimit) ...
                    * cdfac(cdbeta, classCounts(X, phOff, nph, ind, R), fromIR(j,2));
            case SchedStrategy.GPS
                % Generalized PS: share w_r/(w.c) where c_s = 1{n_s>0}, i.e.
                % weights are split across the *active* classes only.
                a{j} = @(X, bufs, svc) rateOf(j) * kirFrac(X, fromIdx(j), phOff, nph, fromIR(j,1), fromIR(j,2)) ...
                    * gpsshare(wrow, classCounts(X, phOff, nph, ind, R), fromIR(j,2)) ...
                    * lldfac(ldrow, sum(classCounts(X, phOff, nph, ind, R)), lldlimit) ...
                    * cdfac(cdbeta, classCounts(X, phOff, nph, ind, R), fromIR(j,2));
            case SchedStrategy.PSPRIO
                % Below capacity every job is served, so priority is inert;
                % above it, only the most urgent non-empty group shares the
                % servers. lld uses the priority-group population, cd the full
                % one, mirroring State.afterEventStation exactly.
                a{j} = @(X, bufs, svc) rateOf(j) * kirFrac(X, fromIdx(j), phOff, nph, fromIR(j,1), fromIR(j,2)) ...
                    * psprioshare(classCounts(X, phOff, nph, ind, R), fromIR(j,2), mi(fromIR(j,1)), classprio) ...
                    * lldfac(ldrow, prioPop(classCounts(X, phOff, nph, ind, R), fromIR(j,2), mi(fromIR(j,1)), classprio), lldlimit) ...
                    * cdfac(cdbeta, classCounts(X, phOff, nph, ind, R), fromIR(j,2));
            case SchedStrategy.DPSPRIO
                % As DPS, but above capacity restricted to the most urgent
                % non-empty group; cd is evaluated on the priority-restricted
                % population (State.afterEventStation, case DPSPRIO).
                a{j} = @(X, bufs, svc) rateOf(j) * kirFrac(X, fromIdx(j), phOff, nph, fromIR(j,1), fromIR(j,2)) ...
                    * dpsprioshare(wrow, classCounts(X, phOff, nph, ind, R), fromIR(j,2), mi(fromIR(j,1)), classprio) ...
                    * lldfac(ldrow, prioPop(classCounts(X, phOff, nph, ind, R), fromIR(j,2), mi(fromIR(j,1)), classprio), lldlimit) ...
                    * cdfac(cdbeta, prioVec(classCounts(X, phOff, nph, ind, R), fromIR(j,2), mi(fromIR(j,1)), classprio), fromIR(j,2));
            case SchedStrategy.GPSPRIO
                a{j} = @(X, bufs, svc) rateOf(j) * kirFrac(X, fromIdx(j), phOff, nph, fromIR(j,1), fromIR(j,2)) ...
                    * gpsprioshare(wrow, classCounts(X, phOff, nph, ind, R), fromIR(j,2), mi(fromIR(j,1)), classprio) ...
                    * lldfac(ldrow, prioPop(classCounts(X, phOff, nph, ind, R), fromIR(j,2), mi(fromIR(j,1)), classprio), lldlimit) ...
                    * cdfac(cdbeta, prioVec(classCounts(X, phOff, nph, ind, R), fromIR(j,2), mi(fromIR(j,1)), classprio), fromIR(j,2));
            case SchedStrategy.PAS
                % Position p of the ordered list is served at
                % Delta_mu(c1..cp) = mu(c1..cp) - mu(c1..c_{p-1}), and
                % pass-and-swap decides which class that completion ejects. The
                % class-r departure rate is therefore the total Delta_mu over
                % the positions whose pass-and-swap ejects a class-r job, which
                % is exactly what afterEventStationPAS enumerates.
                muFun = sn.nodeparam{ind}.svcRateFun;
                if isempty(muFun)
                    line_error(mfilename, 'PAS/OI station has no service rate function mu(c); set it via setService(@(c) ...).');
                end
                swapG = sn.nodeparam{ind}.swapGraph;
                a{j} = @(X, bufs, svc) oirate(muFun, swapG, bufs{ind}, fromIR(j,2));
            case SchedStrategy.POLLING
                % A polling station has a single server that serves exactly one
                % job, of the class its controller currently attends. The
                % departure of class r therefore fires only while the controller
                % is SERVING class r, at the plain service rate of the one job in
                % service -- never scaled by the class population, since the other
                % class-r jobs wait in the buffer for the server to come back to
                % them. The controller rides in bufs{ind} = [mode, pos, swk, ctr];
                % pollServeGate returns 1 exactly when mode==SERVING and pos==r.
                a{j} = @(X, bufs, svc) rateOf(j) * pollServeGate(bufs{fromIR(j,1)}, fromIR(j,2)) ...
                    * lldfac(ldrow, sum(classCounts(X, phOff, nph, ind, R)), lldlimit) ...
                    * cdfac(cdbeta, classCounts(X, phOff, nph, ind, R), fromIR(j,2));
            case {SchedStrategy.FCFS, SchedStrategy.LCFS, SchedStrategy.SIRO, ...
                    SchedStrategy.HOL, SchedStrategy.SEPT, SchedStrategy.LEPT, ...
                    SchedStrategy.LCFSPR}
                % Invariant: numel(bufs{ind}) == max(0, total - mi(ind)).
                % Rate is proportional to the jobs actually being served, i.e.
                % the class-r population minus the class-r jobs waiting in buffer,
                % scaled by the load-dependent factor at the total population.
                % Every non-preemptive buffered policy shares this law: with
                % exponential service the departure rate depends only on the
                % in-service composition, never on the buffer order, which enters
                % solely through which job is promoted next (pickFromBuffer).
                a{j} = @(X, bufs, svc) rateOf(j) * kirFrac(X, fromIdx(j), phOff, nph, fromIR(j,1), fromIR(j,2)) * ...
                    max(0, classPop(X, phOff, nph, ind, fromIR(j,2)) - sum(bufs{fromIR(j,1)} == fromIR(j,2))) ...
                    * lldfac(ldrow, sum(classCounts(X, phOff, nph, ind, R)), lldlimit) ...
                    * cdfac(cdbeta, classCounts(X, phOff, nph, ind, R), fromIR(j,2));
        end
    else
        a{j} = @(X, bufs, svc) rateOf(j) * kirFrac(X, fromIdx(j), phOff, nph, fromIR(j,1), fromIR(j,2)) * min(1, classPop(X, phOff, nph, ind, fromIR(j,2)));
    end
end

% Reneging propensities ------------------------------------------------
% Only the jobs actually waiting can abandon, so the rate is the class-r
% population minus the class-r jobs in service, exactly the buffer occupancy
% the FCFS-family rate law already relies on.
for j = (nDep+1):length(fromIdx)
    ind = fromIR(j,1);
    a{j} = @(X, bufs, svc) renegeMu(j) * sum(bufs{ind} == fromIR(j,2));
end

% Retrial propensities --------------------------------------------------
% Only jobs actually in orbit retry, and only a free server admits them.
for j = 1:length(fromIdx)
    if ~isRetryRx(j)
        continue
    end
    ind = fromIR(j,1);
    a{j} = @(X, bufs, svc) retryMu(j) * sum(bufs{ind} == fromIR(j,2)) ...
        * double(sum(X(((ind-1)*R+1):(ind*R))) - numel(bufs{ind}) < mi(ind));
end

% Polling switchover propensities ---------------------------------------
% The reneging loop above overwrote these appended columns with a zero-rate
% renege closure; restore the switchover law here. A switchover fires only
% while the controller is walking (mode SWITCHING, bufs{ind}(1)==2), at the
% total leaving rate -D0(swk,swk) of the current phase swk of the switchover
% PH into buffer pos. The competition between advancing to another phase and
% absorbing (arriving at pos) is resolved at firing time, exactly as a routed
% departure resolves its destination after it fires.
for j = 1:length(fromIdx)
    if ~isPollSwRx(j)
        continue
    end
    ind = pollSwNode(j);
    pinf = poll.pinfo{ind};
    a{j} = @(X, bufs, svc) pollSwRate(bufs{ind}, pinf);
end

% Finite capacity regions do NOT gate the propensities ------------------
% Under the DROP rule the refused job is DESTROYED, not held back: the
% departure fires at its full rate and the job simply never reaches the
% destination. Scaling the propensity by the admitted share instead censors
% the transition, which keeps the job at its SOURCE -- a different model, and
% one that diverges as soon as the source is a real queue rather than a Source
% node (an interior region makes the upstream queue grow without bound while
% nothing is ever lost). The two coincide only at a Source, whose population is
% fictitious, which is why every FCR fixture placed a region on a Source-fed
% station and never saw the difference. Refusal is applied at firing time
% instead, on the drawn destination, exactly as a balk is (see balkDraw below):
% the source releases the job and the destination never receives it. This
% matches SOLVER_SSA (which marks the refusal and suppresses only the passive
% application) and the exact CTMC.

% Propensity functions dependencies -----------------------------------
D = cell(1,size(S,2));
for k=1:size(D,2)
    J = find(S(:,k))'; % set of state variables affected by reaction k
    vecd = [];
    for j=1:length(J)
        % Decode through the slot map, never arithmetically: with phase
        % expansion a state index is a (node,class,PHASE) slot, so
        % mod(pos-1,R)+1 names the wrong node as soon as any class has more
        % than one phase. Collect EVERY slot of each affected node, because a
        % rate law reads its node's whole class-count vector (classCounts sums
        % each class over its phases) and the per-phase share reads the sibling
        % phases of its own class.
        ind = slotNode(J(j));
        % NB: not `rr` -- that name holds the round-robin controller in this
        % scope, and shadowing it here would pass an integer to the run loop.
        for rcls = 1:R
            vecd(end+1:end+nph(ind,rcls)) = (phOff(ind,rcls)+1):(phOff(ind,rcls)+nph(ind,rcls));
        end
    end
    % vecd now contains all state variables affected by the firing of
    % reaction k. We now find the propensity functions that depend
    % on those variables
    if isRetryRx(k)
        % A retry has an all-zero stoichiometry column, so the generic
        % derivation below would return an empty dependency set and leave every
        % rate at the node stale. A retry does change the in-service
        % composition, hence every reaction whose source is this node.
        base_k = (fromIR(k,1)-1)*R;
        D{k} = find(ismember(fromIdx, (base_k+1):(base_k+R)));
        continue
    end
    if ~isempty(vecd)
        vecd = unique(vecd);
        vecs = [];
        for j=1:length(vecd)
            % No `fcr.on` widening here: regions no longer gate the
            % propensities (see the FCR note above), so a departure's rate
            % depends only on its own station's populations, as in the
            % unregulated case. Admission is resolved at firing time on the
            % drawn destination and changes no rate.
            vecs = [vecs,find(S(vecd(j),:)<0)];
        end
        D{k} = unique(vecs);
    else
        D{k} = [];
    end
end

% A retry has an all-zero stoichiometry column, so the derivation above -- which
% collects reactions by the sign of their S entries -- can never place it in any
% OTHER reaction's dependency set. It still has to be refreshed whenever the
% node it serves changes, because its rate reads both the orbit occupancy and
% whether a server is free: without this, a retry blocked at a busy server keeps
% its zero rate after the server frees, the orbit never drains and the station
% grows without bound.
for j = find(isRetryRx(:)')
    indj = fromIR(j,1);
    slots = ((indj-1)*R + 1):(indj*R);
    for k = 1:size(S,2)
        if any(S(slots, k) ~= 0) || fromIR(k,1) == indj
            if ~ismember(j, D{k})
                D{k}(end+1) = j;
            end
        end
    end
end

% Having accounted for them in D, we can now remove self-loops markings
S(isinf(S))=0;

% ---------------------------------------------------------------------
% Initialize performance metric matrices
% ---------------------------------------------------------------------
lG = 0; % Not computed in SSA

% ---------------------------------------------------------------------
% Run SSA/NRM with direct metric computation
% ---------------------------------------------------------------------
[QN, UN, RN, TN, CN, XN, cacheProd, StartN, PreemptN] = next_reaction_method_direct(S, D, a, nvec0, buffers0, samples, options, sn, fromIdx, fromIR, mi, fcr, balk, isRenegeRx, sig, rr, isRetryRx, phOff, nph, nDepRx, isPhaseRx, smap, poll, isPollSwRx, pollSwNode, svcph0, isBufSvcRx, bufPHClass, bufPHNode, depPhase, phaseFrom, phaseTo, isCacheRx, cacheHitSlot, cacheMissSlot, isCacheNode, cacheRetrDest);
% Write the measured hit/miss probabilities back into sn so the analyzer can set
% them on each Cache node (State.afterEventCache convention: actualhitprob(r) =
% hit throughput / (hit+miss) throughput at the cache, per read class r).
% The cache hit/miss probability of a read class is the throughput of its hit
% class over hit+miss at the cache -- exactly what cacheProd counts per produced
% class. A retrieval completion produces the miss class, so retrieval misses are
% counted here too, and it releases the requests merged onto that fetch, each
% counted in its own hit class (matching the serial engine). Retrieval
% classes (hitclass == 0) are internal and get no hit/miss probability of their
% own.
for ind = 1:I
    if isCacheNode(ind)
        np = sn.nodeparam{ind};
        % Size to nclasses with NaN defaults, exactly as the serial analyzer does:
        % the arrival-rate reconstruction (sn_get_arvr_from_tput) indexes
        % actual{hit,miss}prob at every origClass whose missclass is set, which
        % includes the internal retrieval classes.
        np.actualhitprob  = NaN(1, R);
        np.actualmissprob = NaN(1, R);
        for r = 1:R
            if isCacheReadClass(ind,r) && r <= numel(np.hitclass) && np.hitclass(r) > 0
                hc = np.hitclass(r); mc = np.missclass(r);
                hcount = cacheProd(ind, hc);
                mcount = cacheProd(ind, mc);
                tot = hcount + mcount;
                if tot > 0
                    np.actualhitprob(r)  = hcount / tot;
                    np.actualmissprob(r) = mcount / tot;
                end
            end
        end
        sn.nodeparam{ind} = np;
    end
end

end  % solver_ssa_nrm

% ======================================================================
% Next-Reaction Method with direct metric computation
% ======================================================================
function [QN, UN, RN, TN, CN, XN, cacheProd, StartN, PreemptN] = next_reaction_method_direct(S, D, a, nvec0, buffers0, samples, options, sn, fromIdx, fromIR, mi, fcr, balk, isRenegeRx, sig, rr, isRetryRx, phOff, nph, nDepRx, isPhaseRx, smap, poll, isPollSwRx, pollSwNode, svcph0, isBufSvcRx, bufPHClass, bufPHNode, depPhase, phaseFrom, phaseTo, isCacheRx, cacheHitSlot, cacheMissSlot, isCacheNode, cacheRetrDest)

numReactions = size(S,2);
R = sn.nclasses;
I = sn.nnodes;
M = sn.nstations;
K = sn.nclasses;
QN = zeros(M, K);
UN = zeros(M, K);
RN = zeros(M, K);
TN = zeros(M, K);
CN = zeros(1, K);
XN = zeros(1, K);

% when a reaction fires, this matrix helps selecting the probability that a
% particular routing or phase is selected as a result ------------------
P = S; P(P<0)=P(P<0)+1';
fromIdxCell = cell(numReactions,1);
toIdxCell = cell(numReactions,1);
cdfVec = cell(numReactions,1);
for r=1:numReactions
    nnzP(r) = nnz(P(:,r));
    if nnzP(r)>1
        fromIdxCell{r} = find(S(:,r)<0);
        toIdxCell{r} = find(P(:,r));
        cdfVec{r} = cumsum(P(toIdxCell{r},r));
    end
end

% JSQ routing: reactions whose source class routes with JSQ select the
% destination node holding the smallest total population at firing time
% (each candidate evaluated on its own queue, never the routing node's;
% ties split uniformly)
isJSQ = false(numReactions,1);
% SQ(d) (shortest queue of d): sample d candidates uniformly WITH
% replacement, join the one holding the smallest total population, ties broken
% by first occurrence in the sampled tuple. This is the sampled form of the
% marginal enumerated by sub_sq in MNetwork.refreshRoutingMatrix and of
% LDES's selectSQDestination; drawing directly is equivalent for
% a simulator and avoids enumerating the ndest^d tuples. Dispatcher memory is
% not supported, so the draw is a pure function of the current populations.
isKCH = false(numReactions,1);
kchK = zeros(numReactions,1);
for r=1:numReactions
    if nnzP(r)>1 && fromIdx(r) > 0
        srcNode = smap.node(fromIdx(r));
        srcClass = smap.class(fromIdx(r));
        if sn.routing(srcNode, srcClass) == RoutingStrategy.JSQ
            isJSQ(r) = true;
        elseif sn.routing(srcNode, srcClass) == RoutingStrategy.SQ
            isKCH(r) = true;
            kk = 2; % sub_sq default when nodeparam carries no d
            if iscell(sn.nodeparam) && srcNode <= numel(sn.nodeparam) ...
                    && iscell(sn.nodeparam{srcNode}) && srcClass <= numel(sn.nodeparam{srcNode})
                np = sn.nodeparam{srcNode}{srcClass};
                if ~isempty(np) && isfield(np,'d') && ~isempty(np.d)
                    kk = np.d;
                end
            end
            kchK(r) = max(1, min(kk, numel(toIdxCell{r})));
        end
    end
end

% initialise Gillespie clocks ------------------------------------------
t   = 0;
buffers = buffers0; % working copy of the per-node ordered buffers
svcph = svcph0;     % working copy of the in-service phase multiset (buffered-PH)
cacheProd = zeros(numel(buffers0), sn.nclasses); % per (cache node, PRODUCED class) count
% Seed each polling controller into the auxiliary buffer of its node. The seed
% is a member of the reachable controller space (State.pollingInit's rule): the
% server walks from a canonical position and settles on the first tangible
% state -- a visit on a class with work, a switchover, or a park -- so the
% initial state carries no controller configuration the dynamics cannot reach.
if poll.on
    for ind = 1:sn.nnodes
        if ~poll.isPoll(ind)
            continue
        end
        pinf = poll.pinfo{ind};
        nbuf = classCounts(nvec0, phOff, nph, ind, R)';   % 1xR per-class populations
        [q0, mode0, budget0] = State.pollingNext(pinf, 1, nbuf, R, true);
        buffers{ind} = pollLandCtrl(pinf, q0, mode0, budget0);
    end
end
% Per-region WAITQ FIFO of parked (dstNode, dstClass) tokens, encoded as
% (dstNode-1)*R + dstClass. Empty and untouched unless a region uses WAITQ.
fcrBuf = {};
if fcr.on
    fcrBuf = repmat({zeros(1,0)}, numel(fcr.classCap), 1);
end
for k=1:size(S,2)
    Ak(k) = a{k}(nvec0, buffers, svcph);
end
nvec   = nvec0;
Pk  = -log(rand(1,numReactions));
Tk  = zeros(1,numReactions);

tau = (Pk - Tk) ./ Ak;

% Performance tracking variables
totalTime = 0;
% Derived START/PREEMPT tallies. The NRM fires one reaction at a time and
% knows exactly which job takes a server and which is displaced, so these are
% COUNTS of events; dividing by the simulated time at the end gives the same
% rate the serial engine estimates from the enabled-transition rates.
startCount = zeros(sn.nnodes, R);
preemptCount = zeros(sn.nnodes, R);
% Departures that were BLOCKED, per (node, class). TN integrates the PROPENSITY,
% which counts a departure the station never makes once its successor is full
% (0.744 against the exact 0.652 on the BUG-81 tandem), so the blocked firings
% are subtracted from that integral before it is normalized. In expectation the
% count IS the integral of the blocked share of the rate, so the difference is
% unbiased -- and unlike recomputing that share it needs no second evaluation of
% a state-dependent dispatcher, whose draw would otherwise have to be replayed.
blockCount = zeros(sn.nnodes, R);
NK = sn.njobs'; % Jobs per class
servers = sn.nservers;
PH = sn.proc; % service-process MAPs/PHs

% Normalized DPS/GPS weights and class priorities, mirroring the propensity
% construction so the utilization accumulators use identical sharing factors.
classprio = sn.classprio(:)';
wnorm = zeros(M, R);
for ist = 1:M
    if any(sn.sched(ist) == [SchedStrategy.DPS, SchedStrategy.GPS, ...
            SchedStrategy.DPSPRIO, SchedStrategy.GPSPRIO])
        wnorm(ist, :) = sn.schedparam(ist, 1:R) / sum(sn.schedparam(ist, 1:R));
    end
end

% Closed-class destination blocking; inert (blk.on false) unless some station
% carries a cap a closed class can actually reach. See capacityBlock.
blk = blockPrecompute(sn, R);

n = 1;
while n <= samples
    [dt, kfire] = min(tau);
    if isinf(dt), line_error(mfilename,'Deadlock. Quitting nrm method.'); end

    totalTime = totalTime + dt;

    % Accumulate state-dependent metrics during this time interval
    for ist = 1:M
        ind = sn.stationToNode(ist);
        for k = 1:K
            % nvec counts jobs per phase now, so the class population is the
            % sum over that class's phases.
            currentPop = classPop(nvec, phOff, nph, ind, k);

            % Accumulate queue length (QN)
            QN(ist, k) = QN(ist, k) + currentPop * dt;

            % Compute throughput contribution from departures
            % Throughput is the total absorption rate of the class: with phase
            % expansion each phase owns its own departure reaction, so they are
            % summed. Phase-change reactions move no job and are excluded, as
            % are the appended renege/retry columns.
            depRate = 0;
            for jd = 1:nDepRx
                if fromIR(jd,1) == ind && fromIR(jd,2) == k && ~isPhaseRx(jd)
                    depRate = depRate + Ak(jd);
                end
            end
            TN(ist, k) = TN(ist, k) + depRate * dt;

            % Compute utilization based on scheduling policy. For the whole PS
            % family the class-k utilization is the share of service capacity
            % it receives divided by the server count, so the same sharing
            % factors that define the propensities are reused here (without the
            % lld/cd rate scalings, which rescale work but not occupancy; the
            % lld and cd stations are overridden with the work-based T*S/peak
            % after the loop, where the two conventions part company).
            switch sn.sched(ist)
                case {SchedStrategy.INF, SchedStrategy.EXT}
                    UN(ist, k) = UN(ist, k) + currentPop * dt;
                case {SchedStrategy.PS, SchedStrategy.LPS}
                    totalPop = sum(classCounts(nvec, phOff, nph, ind, R));
                    if totalPop > 0
                        utilization = (currentPop / totalPop) * min(servers(ist), totalPop) / servers(ist);
                    else
                        utilization = 0;
                    end
                    UN(ist, k) = UN(ist, k) + utilization * dt;
                case SchedStrategy.DPS
                    npop = classCounts(nvec, phOff, nph, ind, R);
                    UN(ist, k) = UN(ist, k) + dpsshare(wnorm(ist,:), npop, k) / servers(ist) * dt;
                case SchedStrategy.GPS
                    npop = classCounts(nvec, phOff, nph, ind, R);
                    UN(ist, k) = UN(ist, k) + gpsshare(wnorm(ist,:), npop, k) / servers(ist) * dt;
                case SchedStrategy.PSPRIO
                    npop = classCounts(nvec, phOff, nph, ind, R);
                    UN(ist, k) = UN(ist, k) + psprioshare(npop, k, servers(ist), classprio) / servers(ist) * dt;
                case SchedStrategy.DPSPRIO
                    npop = classCounts(nvec, phOff, nph, ind, R);
                    UN(ist, k) = UN(ist, k) + dpsprioshare(wnorm(ist,:), npop, k, servers(ist), classprio) / servers(ist) * dt;
                case SchedStrategy.GPSPRIO
                    npop = classCounts(nvec, phOff, nph, ind, R);
                    UN(ist, k) = UN(ist, k) + gpsprioshare(wnorm(ist,:), npop, k, servers(ist), classprio) / servers(ist) * dt;
                case SchedStrategy.PAS
                    % Pass-and-swap / order-independent: utilization is the
                    % time-average number of in-service jobs per class over the
                    % servers, where "in service" means the positions whose
                    % marginal rate increment Delta_mu is positive -- so a job
                    % served by several server types still counts once, not
                    % 1/rate (solver_ctmc_analyzer, case PAS).
                    UN(ist, k) = UN(ist, k) + pasInSvc(sn, ind, buffers{ind}, k) / servers(ist) * dt;
                case {SchedStrategy.FCFS, SchedStrategy.LCFS, SchedStrategy.SIRO, ...
                        SchedStrategy.HOL, SchedStrategy.SEPT, SchedStrategy.LEPT, ...
                        SchedStrategy.LCFSPR}
                    if ~isempty(PH{ist}{k})
                        waiting = sum(buffers{ind} == k);
                        inService = currentPop - waiting;
                        UN(ist, k) = UN(ist, k) + (inService / servers(ist)) * dt;
                    end
                case SchedStrategy.POLLING
                    % The single server is busy on exactly one class-k job while
                    % the controller serves class k, and idle (switching or
                    % parked) otherwise; so class-k utilization is the fraction
                    % of time the controller is SERVING class k.
                    ctrl = buffers{ind};
                    if numel(ctrl) >= 2 && ctrl(1) == 1 && ctrl(2) == k
                        UN(ist, k) = UN(ist, k) + dt / servers(ist);
                    end
            end
        end
    end

    t  = t + dt;

    % update aggregate state
    destPos = [];
    cacheChanged = false;
    % A firing has THREE outcomes, not two: it moves the job, it loses it (the
    % source departs either way), or it is BLOCKED -- cancelled outright, with
    % the source keeping the job and no slot, buffer or controller changing.
    firingBlocked = false;
    if isCacheRx(kfire)
        % Cache access. The read-class job at the cache reads an item drawn from
        % pread, and the cache contents (carried in buffers{cacheNode}) decide a
        % hit or a miss; the replacement policy then rewrites the contents.
        % Mirrors State.afterEventCache (READ, isSimulation). The job leaves in
        % the hit or miss class at the SAME cache node, and the existing
        % immediate forwarding routes it downstream from there.
        cn = fromIR(kfire,1); rdc = fromIR(kfire,2);
        [outClass, newContents, cacheCat, cacheReleased] = cacheAccess(sn, cn, rdc, buffers{cn});
        buffers{cn} = newContents;
        nvec(fromIdx(kfire)) = nvec(fromIdx(kfire)) - 1;   % consume the read-class job
        if outClass > 0
            if cacheCat == 4
                % BEGIN retrieval: the job must travel to the fetch queue and
                % return before the miss completes, so it is placed at the
                % retrieval class's routed destination (the queue), NOT left at
                % the cache where the cache-access reaction would fire again.
                destPos = cacheRetrDest(cn, outClass);
            else
                % Hit or miss/completion: the job leaves in the hit or miss class
                % at the SAME cache node; the existing immediate forwarding routes
                % it downstream. Count the production per produced class so the
                % hit/miss probabilities are the hit/miss-class throughput at the
                % cache (State.afterEventCache convention).
                destPos = phOff(cn, outClass) + 1;
                cacheProd(cn, outClass) = cacheProd(cn, outClass) + 1;
            end
            nvec(destPos) = nvec(destPos) + 1;
        end
        % OUTCLASS == 0 is a request merged onto a pending fetch: it produces
        % nothing now and waits in block B. A completing fetch releases the
        % requests merged onto it, each as a delayed hit in its own hit class.
        for rel = 1:size(cacheReleased,1)
            hc = cacheReleased(rel,1);
            cnt = cacheReleased(rel,2);
            relPos = phOff(cn, hc) + 1;
            nvec(relPos) = nvec(relPos) + cnt;
            cacheProd(cn, hc) = cacheProd(cn, hc) + cnt;
        end
        cacheChanged = true;
    elseif nnzP(kfire)>1
        cand = toIdxCell{kfire};
        % A finite capacity region does NOT filter the routing draw. Routing
        % picks the destination first and the region decides admission at the
        % destination's entry afterwards, dropping the job on refusal; a
        % routing strategy that steered around full regions would be a
        % different (and better-behaved) model than the one SOLVER_SSA and the
        % CTMC implement. The refusal check is applied to the drawn destination
        % below.
        if isJSQ(kfire)
            % JSQ: join the destination node with the smallest total population
            % (ties split uniformly)
            npop = inf(numel(cand),1);
            for x=1:numel(cand)
                jnd = smap.node(cand(x));
                npop(x) = sum(classCounts(nvec, phOff, nph, jnd, R));
            end
            amins = find(npop == min(npop));
            r = amins(1 + floor(rand*length(amins)));
        elseif rr.on && rr.isrr(fromIR(kfire,1), fromIR(kfire,2))
            % Round-robin: advance the pointer, then take the destination it
            % lands on (State.afterEventRouter advances on DEP and the routing
            % closure reads state_after).
            [rr, jnd] = rrNext(rr, fromIR(kfire,1), fromIR(kfire,2));
            % A phase-type destination contributes ONE candidate per entry phase,
            % each weighted by pentry in the routing matrix. The pointer fixes the
            % NODE; the entry PHASE must still be drawn from pentry among that
            % node's candidates. Taking the first match (phase 0) biases the
            % service time -- the RROBIN + phase-type residence bug (RUN-10). Use
            % smap.node (not the flat floor((cand-1)/R) formula, which is wrong
            % once phases expand the state) to find the node's candidates, then
            % sample among them in proportion to their routing weights.
            matches = [];
            for x = 1:numel(cand)
                if smap.node(cand(x)) == jnd
                    matches(end+1) = x; %#ok<AGROW>
                end
            end
            if isempty(matches)
                line_error(mfilename, sprintf('Round-robin selected node %d, which is not a routing destination of node %d.', jnd, fromIR(kfire,1)));
            end
            if numel(matches) == 1
                r = matches(1);
            else
                cd = cdfVec{kfire};
                w = zeros(numel(matches),1);
                for ii = 1:numel(matches)
                    x = matches(ii);
                    if x > 1
                        w(ii) = cd(x) - cd(x-1);
                    else
                        w(ii) = cd(x);
                    end
                end
                wsum = sum(w);
                if wsum <= 0
                    r = matches(1);
                else
                    u = rand * wsum; acc = 0; r = matches(end);
                    for ii = 1:numel(matches)
                        acc = acc + w(ii);
                        if acc > u
                            r = matches(ii);
                            break
                        end
                    end
                end
            end
        elseif isKCH(kfire)
            npop = zeros(numel(cand),1);
            for x=1:numel(cand)
                jnd = smap.node(cand(x));
                npop(x) = sum(classCounts(nvec, phOff, nph, jnd, R));
            end
            % SQ(d): the d candidates are drawn uniformly with replacement
            % and the least loaded wins; the strict comparison retains the
            % first occurrence, which is the tie rule of sub_sq.
            draws = min(kchK(kfire), numel(cand));
            r = 1;
            bestpop = inf;
            for t = 1:draws
                x = 1 + floor(rand*numel(cand));
                if npop(x) < bestpop
                    bestpop = npop(x);
                    r = x;
                end
            end
        else
            % Inverse-CDF sampling: smallest r such that cdfVec(r) > rand. The
            % previous formulation `1+find(rand>=cdfVec,1)` was a misuse of
            % find(...,1) that always returned 2 once rand exceeded cdfVec(1),
            % leaving destinations beyond the second one unreachable (e.g. all
            % traffic skipping Station3 in a 3-way RAND split).
            r = find(cdfVec{kfire} > rand, 1);
            if isempty(r)
                r = length(cdfVec{kfire});
            end
        end
        % Balking is decided on the pre-arrival population, so it is drawn
        % before the state is updated. A balked job is lost: the source still
        % releases it, the destination never receives it.
        % A closed job that finds no room BLOCKS: the firing is cancelled
        % before any gate that would consume it, because a job that cannot
        % leave its station never gets the chance to balk or to be dropped.
        firingBlocked = blk.on && capacityBlock(blk, nvec, toIdxCell{kfire}(r), fromIdx(kfire), R, smap);
        balked = false;
        if ~firingBlocked && balk.on
            balked = balkDraw(balk, nvec, toIdxCell{kfire}(r), R, smap);
        end
        % An open arrival at a full physically-capped destination is lost,
        % exactly as a balked one is: the source releases it, the destination
        % never receives it. Mirrors State.afterEventStation.
        if ~firingBlocked && ~balked && capacityLoss(sn, nvec, toIdxCell{kfire}(r), R, smap)
            balked = true;
        end
        % A region refuses the drawn destination on the same pre-arrival
        % population. Under DROP the refused job is lost, exactly as a balked
        % one is; under WAITQ it is parked in the refusing region's FIFO and
        % admitted later, head-of-line. Either way it does not enter the
        % destination now, so the source still departs and destPos is cleared.
        parkF = 0; parkTok = 0;
        if ~firingBlocked && ~balked && fcr.on
            dstN = smap.node(toIdxCell{kfire}(r));
            dstC = smap.class(toIdxCell{kfire}(r));
            fref = fcrRefusingRegion(fcr, nvec, fromIR(kfire,1), fromIR(kfire,2), dstN, dstC, R, smap);
            if fref ~= 0
                balked = true;
                if fcr.waitq(fref, dstC)
                    parkF = fref; parkTok = (dstN-1)*R + dstC;
                end
            end
        end
        % The source decrement belongs to each outcome separately: a blocked
        % firing is the one case where the job does NOT leave its slot.
        if firingBlocked
            destPos = [];
        elseif balked
            nvec(fromIdxCell{kfire}) = nvec(fromIdxCell{kfire}) - 1;
            destPos = [];
            if parkF > 0
                fcrBuf{parkF}(end+1) = parkTok;
            end
        elseif sig.on && sigIsSignalArrival(sig, toIdxCell{kfire}(r), R, smap)
            % the signal is annihilated on arrival: it never joins the station
            nvec(fromIdxCell{kfire}) = nvec(fromIdxCell{kfire}) - 1;
            [nvec, buffers] = sigApply(sig, nvec, buffers, toIdxCell{kfire}(r), R, mi, smap);
            destPos = [];
        else
            nvec(fromIdxCell{kfire}) = nvec(fromIdxCell{kfire}) - 1;
            nvec(toIdxCell{kfire}(r)) = nvec(toIdxCell{kfire}(r)) + 1;
            destPos = toIdxCell{kfire}(r);
        end
    else
        dpos = find(S(:,kfire) > 0); % deterministic destination (single move)
        % A PHASE change is not an arrival: its destination is another phase slot
        % of the SAME job at the SAME station, so none of the arrival-side gates
        % below may see it. Balking and the capacity gate both LOSE the job when
        % they do; the region gate is inert on a phase change (fcrRefusingRegion
        % discounts the source in the same region, so src==dst cancels and it
        % can never refuse -- measured at 0 refusals in 49998 consultations),
        % and is excluded here to state that invariant rather than to fix a
        % defect. What makes the other two bite HERE and not in the
        % multi-destination branch is the ORDER: there the gates run BEFORE the
        % source decrement and see a true pre-arrival population, while here
        % they run before the update below, so the firing job is still counted
        % at its own station and a phase change at exactly cap reads cap >= cap.
        isPhaseFire = kfire <= numel(isPhaseRx) && isPhaseRx(kfire);
        % Same closed-class block as above; a phase change is not an arrival.
        firingBlocked = blk.on && ~isempty(dpos) && ~isPhaseFire ...
            && ~(kfire <= numel(isRenegeRx) && isRenegeRx(kfire)) ...
            && ~(kfire <= numel(isRetryRx) && isRetryRx(kfire)) ...
            && capacityBlock(blk, nvec, dpos(1), fromIdx(kfire), R, smap);
        balked = false;
        if ~firingBlocked && balk.on && ~isempty(dpos) && ~isPhaseFire
            balked = balkDraw(balk, nvec, dpos(1), R, smap);
        end
        if ~firingBlocked && ~balked && ~isempty(dpos) && ~isPhaseFire ...
                && ~(kfire <= numel(isRenegeRx) && isRenegeRx(kfire)) ...
                && ~(kfire <= numel(isRetryRx) && isRetryRx(kfire)) ...
                && capacityLoss(sn, nvec, dpos(1), R, smap)
            balked = true;
        end
        % Single-destination departures cross region boundaries too, so the
        % region gate applies here exactly as it does to a drawn destination.
        % Renege and retry columns carry no destination and are never gated.
        parkF = 0; parkTok = 0;
        if ~firingBlocked && ~balked && fcr.on && ~isempty(dpos) && ~isPhaseFire ...
                && ~(kfire <= numel(isRenegeRx) && isRenegeRx(kfire)) ...
                && ~(kfire <= numel(isRetryRx) && isRetryRx(kfire))
            dstN = smap.node(dpos(1));
            dstC = smap.class(dpos(1));
            fref = fcrRefusingRegion(fcr, nvec, fromIR(kfire,1), fromIR(kfire,2), dstN, dstC, R, smap);
            if fref ~= 0
                balked = true;
                if fcr.waitq(fref, dstC)
                    parkF = fref; parkTok = (dstN-1)*R + dstC;
                end
            end
        end
        if firingBlocked
            % the departure does not occur at all: nothing moves
            destPos = [];
            dpos = [];
        elseif balked
            % lost or parked on arrival: apply the source departure only
            nvec(fromIdx(kfire)) = nvec(fromIdx(kfire)) - 1;
            destPos = [];
            dpos = [];
            if parkF > 0
                fcrBuf{parkF}(end+1) = parkTok;
            end
        elseif sig.on && ~isempty(dpos) && sigIsSignalArrival(sig, dpos(1), R, smap)
            % the signal is annihilated on arrival: it never joins the station
            nvec(fromIdx(kfire)) = nvec(fromIdx(kfire)) - 1;
            [nvec, buffers] = sigApply(sig, nvec, buffers, dpos(1), R, mi, smap);
            destPos = [];
            dpos = [];
        else
            nvec  = nvec + S(:,kfire);  % zero change for self-loops
        end
        if ~isempty(dpos)
            destPos = dpos(1);
        elseif ~(kfire <= numel(isRenegeRx) && isRenegeRx(kfire)) ...
                && ~(kfire <= numel(isPollSwRx) && isPollSwRx(kfire)) && sn.isslc(fromIR(kfire,2))
            % Self-looping class: the completed job re-enters the same node and
            % class (its stoichiometry is a no-op). At a buffered (FCFS/LCFS)
            % station it must rejoin the buffer so the ordering rotates; point
            % destPos at the source slot so updateBuffers applies the arrival.
            destPos = fromIdx(kfire);
        end
    end

    % maintain the buffers given the source/destination of this firing
    svcChanged = false;
    if firingBlocked
        % nothing moved, so no buffer, controller or service phase may change
        blockCount(fromIR(kfire,1), fromIR(kfire,2)) = blockCount(fromIR(kfire,1), fromIR(kfire,2)) + 1;
    elseif kfire <= numel(isRetryRx) && isRetryRx(kfire)
        % A successful retry moves one orbiting job into the free server. The
        % population is unchanged (it was already counted at the station), so
        % only the orbit shrinks; in-service is read back as population minus
        % orbit occupancy.
        ind = fromIR(kfire,1);
        slot = find(buffers{ind} == fromIR(kfire,2), 1, 'first');
        if ~isempty(slot)
            buffers{ind}(slot) = [];
            % the retrying job seizes the free server: this is where service
            % starts at a retrial station, since its departures never promote
            startCount(ind, fromIR(kfire,2)) = startCount(ind, fromIR(kfire,2)) + 1;
        end
    elseif kfire <= numel(isRenegeRx) && isRenegeRx(kfire)
        % Reneging removes a job that was WAITING, so no server is freed and no
        % queued job is promoted; the abandoning job simply leaves the buffer.
        % State.afterEventStation drops the newest waiting job of the class and
        % notes that for memoryless patience all waiting jobs are exchangeable,
        % so the choice cannot affect the marginal distribution.
        ind = fromIR(kfire,1);
        slot = find(buffers{ind} == fromIR(kfire,2), 1, 'first');
        if ~isempty(slot)
            buffers{ind}(slot) = [];
        end
    elseif kfire <= numel(isPhaseRx) && isPhaseRx(kfire) && bufPHNode(fromIR(kfire,1))
        % A buffered-PH phase transition moves one in-service job between phases
        % of its own service process. It frees no server and adds no arrival, so
        % the buffer is untouched and only svcph changes (INF/PS phase moves are
        % already applied to nvec via the stoichiometry and fall to updateBuffers
        % below as a no-op, as before).
        ind = fromIR(kfire,1); r = fromIR(kfire,2);
        svcph{ind}(r, phaseFrom(kfire)) = svcph{ind}(r, phaseFrom(kfire)) - 1;
        svcph{ind}(r, phaseTo(kfire))   = svcph{ind}(r, phaseTo(kfire))   + 1;
        svcChanged = true;
    elseif isCacheRx(kfire)
        % The cache access already updated the cache contents (buffers{cacheNode})
        % and moved the job to the hit/miss class in the firing block above; there
        % is no job buffer to maintain at a cache node.
    else
        [buffers, svcph, svcChanged, startCount, preemptCount] = updateBuffers(kfire, nvec, buffers, fromIR, destPos, mi, R, sn, smap, svcph, bufPHNode, isBufSvcRx, depPhase, startCount, preemptCount);
    end

    % Polling controller advance. The controller of each polling node lives in
    % its auxiliary buffer as [mode, pos, swk, ctr]; a firing can move it in
    % three ways, mirroring State.afterEventStation exactly (EventType.DEP under
    % SchedStrategy.POLLING, EventType.SWITCH, and the parked-server arrival):
    %   * a service completion at the node ends the visit unless the discipline
    %     still admits another job of the served class, and on ending walks the
    %     cyclic order to the next tangible controller state;
    %   * a switchover reaction advances the switchover PH one phase, or on
    %     absorption arrives at the target buffer and opens a visit or walks on;
    %   * an arrival to a parked server wakes it, and the walk resolves at once
    %     to a visit on the newly present work.
    % Any of these changes a service gate or the switchover rate, so a change is
    % flagged to force a full propensity refresh below (like a WAITQ release).
    pollChanged = false;
    if poll.on && ~firingBlocked
        srcNode = fromIR(kfire,1);
        if kfire <= numel(isPollSwRx) && isPollSwRx(kfire)
            pind = pollSwNode(kfire);
            pinf = poll.pinfo{pind};
            ctrl = buffers{pind};
            posS = ctrl(2); swkS = ctrl(3);
            D0S = pinf.swD0{posS};
            KswS = pinf.Ksw(posS);
            w = zeros(1, KswS + 1);
            for kd = 1:KswS
                if kd ~= swkS && D0S(swkS,kd) > 0
                    w(kd) = D0S(swkS,kd);
                end
            end
            w(KswS + 1) = max(0, -sum(D0S(swkS,:)));   % absorption (D1 row sum)
            pick = drawFromDist(w);
            if pick <= KswS && pick ~= swkS
                ctrl(3) = pick;                        % internal phase advance
                buffers{pind} = ctrl;
            else
                nbufS = classCounts(nvec, phOff, nph, pind, R)';
                [qS, mdS, bgS] = State.pollingNext(pinf, posS, nbufS, R, true);
                buffers{pind} = pollLandCtrl(pinf, qS, mdS, bgS);
            end
            pollChanged = true;
        elseif poll.isPoll(srcNode) && kfire <= nDepRx && ~isPhaseRx(kfire)
            pinf = poll.pinfo{srcNode};
            ctrl = buffers{srcNode};
            posD = ctrl(2); ctrD = ctrl(4);
            nbufD = classCounts(nvec, phOff, nph, srcNode, R)';   % after the departure
            switch pinf.ptype
                case PollingType.EXHAUSTIVE
                    ctrnextD = 0;            goonD = nbufD(posD) > 0;
                case PollingType.GATED
                    ctrnextD = ctrD - 1;     goonD = ctrnextD > 0;
                case PollingType.KLIMITED
                    ctrnextD = ctrD - 1;     goonD = ctrnextD > 0 && nbufD(posD) > 0;
                case PollingType.DECREMENTING
                    ctrnextD = ctrD;         goonD = nbufD(posD) > ctrD;
            end
            if goonD
                buffers{srcNode} = [1, posD, 0, ctrnextD];
            else
                [qD, mdD, bgD] = State.pollingNext(pinf, posD, nbufD, R, false);
                buffers{srcNode} = pollLandCtrl(pinf, qD, mdD, bgD);
            end
            pollChanged = true;
        end
        if ~isempty(destPos) && destPos > 0
            jnd = smap.node(destPos);
            if poll.isPoll(jnd)
                ctrlA = buffers{jnd};
                if ~isempty(ctrlA) && ctrlA(1) == 0
                    pinfA = poll.pinfo{jnd};
                    nbufA = classCounts(nvec, phOff, nph, jnd, R)';   % includes the arrival
                    [qA, mdA, bgA] = State.pollingNext(pinfA, ctrlA(2), nbufA, R, true);
                    buffers{jnd} = pollLandCtrl(pinfA, qA, mdA, bgA);
                    pollChanged = true;
                end
            end
        end
    end

    % WAITQ: admit parked jobs whose regions this firing may have relieved.
    % A release changes populations at arbitrary destination nodes, so when
    % anything is admitted every reaction is refreshed rather than only the
    % dependency set of the fired reaction.
    nReleased = 0;
    if fcr.on && fcr.anyWaitq && ~firingBlocked
        [nvec, buffers, fcrBuf, nReleased, svcph, relChanged, startCount, preemptCount] = ...
            fcrReleaseCascade(fcr, nvec, buffers, fcrBuf, mi, R, sn, smap, svcph, bufPHNode, startCount, preemptCount);
        svcChanged = svcChanged || relChanged;
    end

    Tk = Tk + Ak * dt;

    % update rates for all reactions dependent on the last fired reaction. A
    % polling controller move or a WAITQ release can change rates outside the
    % static dependency set of the fired reaction (a switchover reaction has an
    % all-zero stoichiometry column, and a controller move flips service gates),
    % so either forces a full refresh.
    if nReleased > 0 || pollChanged || svcChanged || cacheChanged
        for k=1:numReactions
            Ak(k) = a{k}(nvec, buffers, svcph);
        end
    else
        for k=D{kfire}
            Ak(k) = a{k}(nvec, buffers, svcph);
        end
    end

    % update clocks
    Pk(kfire)  = Pk(kfire) - log(rand);
    tau        = (Pk - Tk) ./ Ak;
    tau(Ak==0) = inf;

    % do not count immediate events
    n = n + 1;
    print_progress(options, n, t);
end % while
% The counter row is closed here rather than newline-terminated:
% line_printf already ends an open row, so an explicit newline was a
% SECOND one and showed as a blank row before the completion banner.
LineStatus.close();

% Normalize metrics by total time
StartN = zeros(M, K);
PreemptN = zeros(M, K);
if totalTime > 0
    for ist = 1:M
        ind = sn.stationToNode(ist);
        for k = 1:K
            QN(ist, k) = QN(ist, k) / totalTime;
            UN(ist, k) = UN(ist, k) / totalTime;
            % net the blocked firings out of the departure-rate integral
            TN(ist, k) = (TN(ist, k) - blockCount(ind, k)) / totalTime;
            % counts of events over the simulated time: a rate, like TN
            StartN(ist, k) = startCount(ind, k) / totalTime;
            PreemptN(ist, k) = preemptCount(ind, k) / totalTime;
        end
    end
end

% Class-dependent stations report utilization as T*S/peak, where peak is the
% declared per-class peak rate scaling (sn.cdscalingpeak). This matches the
% T*S/c convention of the analytic solvers and serial SSA; the accumulated
% in-service fraction above divides by the server count (1 for a cd station),
% which is not the same quantity. Override those stations here.
if ~isempty(sn.cdscaling) || ~isempty(sn.jdscaling)
    for ist = 1:M
        isCd = ~isempty(sn.cdscaling) && ist <= numel(sn.cdscaling) && ~isempty(sn.cdscaling{ist});
        isJd = ~isempty(sn.jdscaling) && ist <= numel(sn.jdscaling) && ~isempty(sn.jdscaling{ist});
        if isCd || isJd
            for k = 1:K
                % Effective peak = product of the declared cd and jd peaks.
                peak = 1;
                if isCd, peak = peak * sn.cdscalingpeak(ist, k); end
                if isJd, peak = peak * sn.jdscalingpeak(ist, k); end
                if isfinite(sn.rates(ist, k)) && sn.rates(ist, k) > 0 && peak > 0
                    UN(ist, k) = TN(ist, k) / sn.rates(ist, k) / peak;
                else
                    UN(ist, k) = 0;
                end
            end
        end
    end
end

% Load-dependent stations report the same work-based utilization, T*S/peak
% against peak = max(c, max(alpha)). The accumulator above integrates BUSY TIME,
% which is a different quantity once alpha(n) ~= 1: a server running alpha(n)
% times faster does the same work in less time, so busy time reads it as no
% busier than one at its nominal rate. That put NRM at 0.9587 on a 4-job closed
% model with alpha = [1 1.5 2 2.5] where CTMC, MVA, NC and serial SSA all report
% 0.6612, and left the two SSA engines disagreeing with each other. INF keeps
% U = Q, as everywhere else.
if ~isempty(sn.lldscaling)
    for ist = 1:M
        if ist > size(sn.lldscaling, 1) || size(sn.lldscaling, 2) == 0
            continue
        end
        lldrow = sn.lldscaling(ist, :);
        if all(lldrow == 1)
            continue
        end
        if any(sn.sched(ist) == [SchedStrategy.INF, SchedStrategy.EXT])
            continue
        end
        peak = max(servers(ist), max(lldrow));
        for k = 1:K
            if isfinite(sn.rates(ist, k)) && sn.rates(ist, k) > 0 && peak > 0
                UN(ist, k) = TN(ist, k) / sn.rates(ist, k) / peak;
            else
                UN(ist, k) = 0;
            end
        end
    end
end

% Compute derived metrics
for k = 1:K
    % System throughput at reference station
    XN(1, k) = TN(sn.refstat(k), k);

    % Response times
    for ist = 1:M
        if TN(ist, k) > 0
            RN(ist, k) = QN(ist, k) / TN(ist, k);
        else
            RN(ist, k) = 0;
        end
    end

    % Cycle times
    if XN(1, k) > 0
        CN(1, k) = NK(k) / XN(1, k);
    end
end

% Handle NaN values
QN(isnan(QN)) = 0;
UN(isnan(UN)) = 0;
RN(isnan(RN)) = 0;
XN(isnan(XN)) = 0;
TN(isnan(TN)) = 0;
CN(isnan(CN)) = 0;

    function print_progress(opt, samples_collected, tnow)
        if LineConsole.isActive() % the console owns the line; see solver_ssa
            every = max(1,round(opt.samples/20));
            if mod(samples_collected, every) == 0
                LineConsole.iter(samples_collected/every, ...
                    'simulated %d of %g samples (%.0f%%), simulated time %.4g', ...
                    samples_collected, opt.samples, ...
                    100*samples_collected/opt.samples, tnow);
            end
            return
        end
        if ~isfield(opt,'verbose') || ~opt.verbose || batchStartupOptionUsed, return; end
        % ONE REWRITTEN FIELD, not a fixed-width one. LineStatus rewinds by
        % the width it actually wrote, so a counter that only grows needs no
        % padding at all and leaves no trailing blanks -- a fixed %-9d field
        % showed its pad as "SSA samples: 100000   ". It also cannot desync
        % the way a hardcoded run of backspaces does once the count outgrows
        % the field. line_printf closes the row, so the completion banner
        % terminates it without help.
        if opt.verbose == 2 || (samples_collected > 0 && mod(samples_collected,1e3) == 0)
            LineStatus.set('SSA samples: %d', samples_collected);
        end
    end
end  % next_reaction_method_direct

% ======================================================================
% Buffer maintenance for FCFS/LCFS nodes
% ======================================================================
function [buffers, svcph, svcChanged, startCount, preemptCount] = updateBuffers(kfire, nvec, buffers, fromIR, destPos, mi, R, sn, smap, svcph, bufPHNode, isBufSvcRx, depPhase, startCount, preemptCount)
% Maintain the ordered per-node buffers when reaction KFIRE fires. A departure
% frees a server, so the buffered job selected by the station's discipline is
% promoted into service and leaves the buffer; an arrival at a buffered
% destination whose servers are all busy joins the buffer head. At a buffered-PH
% node the same events also move jobs in and out of the in-service phase multiset
% svcph, and SVCCHANGED flags that so the caller forces a full propensity refresh
% (svcph is not part of the stoichiometry, so the static dependency set misses it).
ind = fromIR(kfire,1); % source node of the firing
svcChanged = false;

% Buffered-PH departure: the completing job leaves service, so drop it from the
% in-service phase it occupied (carried in depPhase). The promotion below refills
% the freed server from the buffer at a fresh entry phase.
if bufPHNode(ind) && isBufSvcRx(kfire)
    r = fromIR(kfire,2);
    svcph{ind}(r, depPhase(kfire)) = svcph{ind}(r, depPhase(kfire)) - 1;
    svcChanged = true;
end

% An order-independent station keeps the full ordered list, so a departure is
% not a promotion but a pass-and-swap rewrite: the completing position's chain
% shifts classes along and removes one slot. Which position completed is
% redrawn here in proportion to the Delta_mu of the positions that eject this
% class, which is the same split afterEventStationPAS enumerates.
if isListSched(ind, sn) && ~isempty(buffers{ind})
    oldList = buffers{ind};
    buffers{ind} = oiDepart(sn, ind, buffers{ind}, fromIR(kfire,2));
    % An order-independent station serves every position whose rate increment
    % Delta_mu is positive, so a departure starts whichever positions cross
    % from a zero increment to a positive one -- the same rule
    % State.afterEventStationPAS tags. See its sub_startedHere.
    startCount(ind,:) = startCount(ind,:) + oiStarted(sn, ind, oldList, buffers{ind}, R);
    return
end

% Handle departure from a buffered source node: promote one waiting job. A
% retrial station is the exception -- the freed server is NOT filled from the
% orbit, orbiting jobs re-enter only through RETRY events at the memoryless
% retrial rate (State.afterEventStation suppresses promotion likewise).
if isBuffered(ind, sn) && ~isempty(buffers{ind}) && ~isRetrialStation(ind, sn) ...
        && ~isListSched(ind, sn)
    pos = pickFromBuffer(buffers{ind}, sn, sn.nodeToStation(ind));
    promoted = buffers{ind}(pos);
    buffers{ind}(pos) = [];
    startCount(ind, promoted) = startCount(ind, promoted) + 1; % takes the freed server
    if bufPHNode(ind)
        % The promoted waiting job starts service now, entering a phase drawn
        % from its entry distribution pie (the same allocation the init uses).
        ke = drawEntryPhase(sn, ind, promoted, smap.nph(ind, promoted));
        svcph{ind}(promoted, ke) = svcph{ind}(promoted, ke) + 1;
        svcChanged = true;
    end
end

% Handle arrival at a buffered destination node
if ~isempty(destPos) && destPos > 0
    [buffers, svcph, arrChanged, startCount, preemptCount] = applyArrivalBuffer(smap.node(destPos), smap.class(destPos), ...
        nvec, buffers, mi, R, sn, smap, svcph, bufPHNode, startCount, preemptCount);
    svcChanged = svcChanged || arrChanged;
end
end

function st = oiStarted(sn, ind, cold, cnew, R)
% ST=OISTARTED(SN,IND,COLD,CNEW,R) per-class count of the positions of CNEW
% that are served (Delta_mu > 0) and were not served in COLD. Mirrors
% State.afterEventStationPAS, which tags a PAS start by exactly this rule.
st = zeros(1,R);
muFun = sn.nodeparam{ind}.svcRateFun;
if isempty(muFun)
    return
end
incNew = oiIncrements(muFun, cnew);
incOld = oiIncrements(muFun, cold);
for p = 1:numel(incNew)
    if incNew(p) <= 0
        continue
    end
    if p <= numel(incOld) && incOld(p) > 0
        continue
    end
    st(cnew(p)) = st(cnew(p)) + 1;
end
end

function inc = oiIncrements(muFun, c)
% Per-position increments Delta_mu(c1..cp) = mu(c1..cp) - mu(c1..c_{p-1}).
inc = zeros(1, numel(c));
muPrev = 0;
for p = 1:numel(c)
    muCur = muFun(c(1:p));
    inc(p) = muCur - muPrev;
    muPrev = muCur;
end
end

function [buffers, svcph, svcChanged, startCount, preemptCount] = applyArrivalBuffer(jnd, s, nvec, buffers, mi, R, sn, smap, svcph, bufPHNode, startCount, preemptCount)
% Join a just-arrived class-S job to the ordered buffer of destination node
% JND, if that node is buffered. NVEC already includes the arrival. Shared by
% updateBuffers (routed arrivals) and fcrReleaseCascade (WAITQ releases), so
% the two paths cannot drift. At a buffered-PH destination a job that enters
% service (rather than waiting) is added to the in-service phase multiset svcph
% at a pie-drawn entry phase; SVCCHANGED flags that for a propensity refresh.
    svcChanged = false;
    if isListSched(jnd, sn)
        % PAS/OI: the arrival simply joins the back of the ordered list; there
        % is no server/buffer split, so no capacity test against mi. Capacity
        % is the station's own cap, and an arrival past it is lost.
        if numel(buffers{jnd}) < sn.cap(sn.nodeToStation(jnd))
            oldList = buffers{jnd};
            buffers{jnd}(end+1) = s; % append at the back (newest last)
            startCount(jnd,:) = startCount(jnd,:) + oiStarted(sn, jnd, oldList, buffers{jnd}, R);
        end
    elseif isBuffered(jnd, sn)
        totalAtDest = sum(classCounts(nvec, smap.phOff, smap.nph, jnd, R));
        enteredService = false;
        if isRetrialStation(jnd, sn)
            % A retrial station breaks the buffer invariant the other policies
            % share: because a departure does not promote, the orbit can be
            % occupied while servers sit idle, so "total > mi" no longer means
            % "the servers are busy". An arrival must consult the servers
            % directly and only join the orbit when none is free.
            inSvc = (totalAtDest - 1) - numel(buffers{jnd});
            if inSvc >= mi(jnd)
                buffers{jnd} = [s, buffers{jnd}];
            else
                enteredService = true;
            end
        elseif totalAtDest > mi(jnd)
            if isPreemptive(jnd, sn)
                % Preempt-resume: the arrival seizes a server and the incumbent
                % it displaces is the one that joins the buffer. The victim is
                % drawn in proportion to the class occupancies of the servers,
                % as State.afterEventStation weights its preemption branches by
                % si_preempt/sum(space_srv). Buffering the incumbent rather than
                % the arrival is what leaves the new job in service, since
                % in-service is read back as population minus buffer occupancy.
                c = pickPreempted(nvec, buffers{jnd}, jnd, s, R);
                if c > 0
                    buffers{jnd} = [c, buffers{jnd}]; % addFirst
                    preemptCount(jnd, c) = preemptCount(jnd, c) + 1; % displaced incumbent
                end
                enteredService = true;
            else
                % All servers busy - arriving job joins back of buffer
                buffers{jnd} = [s, buffers{jnd}]; % addFirst
            end
        else
            % A server is free: the job goes straight into service.
            enteredService = true;
        end
        if enteredService
            startCount(jnd, s) = startCount(jnd, s) + 1; % the arrival took a server
        end
        if enteredService && bufPHNode(jnd)
            ke = drawEntryPhase(sn, jnd, s, smap.nph(jnd, s));
            svcph{jnd}(s, ke) = svcph{jnd}(s, ke) + 1;
            svcChanged = true;
        end
    end
end

function ke = drawEntryPhase(sn, jnd, s, nphjs)
% Sample the service phase a class-S job starts in at node JND from its entry
% distribution pie. A single-phase class always enters phase 1.
if nphjs <= 1
    ke = 1;
    return
end
pe = entryProbs(sn, jnd, s, nphjs);
ke = drawFromDist(pe);
end

function pos = pickFromBuffer(buf, sn, ist)
% Index of the waiting job that the discipline at station IST promotes into
% service. BUF is ordered newest-first / oldest-last, matching the convention
% of State.afterEventStation's space_buf (which inserts arrivals at column 1
% and, for HOL, promotes the rightmost job of the urgent priority group).
switch sn.sched(ist)
    case SchedStrategy.FCFS
        pos = numel(buf);                      % oldest
    case {SchedStrategy.LCFS, SchedStrategy.LCFSPR}
        pos = 1;                               % newest / most recently preempted
    case SchedStrategy.SIRO
        % Uniform over the waiting jobs. State.afterEventStation promotes a
        % class-r job with probability (nir(r)-sir(r))/(ni-sum(sir)), i.e. the
        % waiting class-r fraction, which is exactly a uniform draw over buf.
        pos = 1 + floor(rand * numel(buf));
    case SchedStrategy.HOL
        % Highest priority (lowest classprio value); FCFS within the group, so
        % the oldest = the last matching position.
        prio = sn.classprio(buf);
        pos = find(prio == min(prio), 1, 'last');
    case {SchedStrategy.SEPT, SchedStrategy.LEPT}
        % sn.schedparam(ist,r) is the rank of class r's mean service time
        % (ascending for SEPT, descending for LEPT), so the promoted class is
        % the waiting one of least rank. Oldest first within a class.
        ranks = sn.schedparam(ist, buf);
        pos = find(ranks == min(ranks), 1, 'last');
    otherwise
        line_error(mfilename, sprintf('pickFromBuffer: unsupported buffered policy %s.', ...
            SchedStrategy.toText(sn.sched(ist))));
end
end

function n = pasInSvc(sn, ind, c, r)
% Number of class-r jobs in service at a PAS/OI station holding the ordered
% list C: the positions whose marginal rate increment Delta_mu is positive.
% Mirrors the sir the PAS branch of State.toMarginal reports, which is what
% solver_ctmc_analyzer divides by the server count.
n = 0;
if isempty(c)
    return
end
muFun = sn.nodeparam{ind}.svcRateFun;
muPrev = 0;
for p = 1:numel(c)
    muCur = muFun(c(1:p));
    if muCur - muPrev > 0 && c(p) == r
        n = n + 1;
    end
    muPrev = muCur;
end
end

function buf = oiDepart(sn, ind, buf, r)
% Apply the pass-and-swap rewrite for a class-r departure at OI station IND.
% The completing position is drawn among those whose pass-and-swap ejects class
% r, weighted by that position's own service rate Delta_mu.
muFun = sn.nodeparam{ind}.svcRateFun;
G = sn.nodeparam{ind}.swapGraph;
c = buf;
n = numel(c);
pos = [];
w = [];
muPrev = 0;
for p = 1:n
    muCur = muFun(c(1:p));
    ratep = muCur - muPrev;
    muPrev = muCur;
    if ratep <= 0
        continue
    end
    [~, depClass] = State.passAndSwap(c, p, G);
    if depClass == r
        pos(end+1) = p; %#ok<AGROW>
        w(end+1) = ratep; %#ok<AGROW>
    end
end
if isempty(pos)
    return % this class cannot depart from the current list
end
u = rand * sum(w);
acc = 0;
pick = pos(end);
for x = 1:numel(pos)
    acc = acc + w(x);
    if u < acc
        pick = pos(x);
        break
    end
end
buf = State.passAndSwap(c, pick, G);
end

function rt = oirate(muFun, G, c, r)
% Aggregate class-r departure rate of an order-independent / pass-and-swap
% station holding the ordered list C (oldest first). Mirrors the DEP branch of
% State.afterEventStationPAS: every position contributes its own service token
% at Delta_mu, and pass-and-swap decides which class actually leaves.
rt = 0;
n = numel(c);
if n == 0
    return
end
muPrev = 0;   % mu of the empty prefix is 0
for p = 1:n
    muCur = muFun(c(1:p));
    ratep = muCur - muPrev;
    muPrev = muCur;
    if ratep <= 0
        continue   % position p receives no service
    end
    [~, depClass] = State.passAndSwap(c, p, G);
    if depClass == r
        rt = rt + ratep;
    end
end
end

function tf = isListSched(ind, sn)
% True for stations whose buffer holds the FULL ordered job list rather than
% only the waiting jobs.
tf = false;
if sn.isstation(ind)
    tf = (sn.sched(sn.nodeToStation(ind)) == SchedStrategy.PAS);
end
end

function tf = isRetrialStation(ind, sn)
% True for stations with a retrial orbit: their freed servers are not filled by
% promotion, only by a successful RETRY.
tf = false;
if sn.isstation(ind) && isfield(sn,'retrialProc') && ~isempty(sn.retrialProc)
    ist = sn.nodeToStation(ind);
    tf = ist > 0 && any(~cellfun(@isempty, sn.retrialProc(ist,:)));
end
end

function tf = isPreemptive(ind, sn)
% True for the preempt-resume / preempt-independent policies, whose arrivals
% displace an incumbent instead of queueing behind it.
tf = false;
if sn.isstation(ind)
    ist = sn.nodeToStation(ind);
    tf = any(sn.sched(ist) == [SchedStrategy.LCFSPR]);
end
end

function c = pickPreempted(nvec, buf, jnd, arrClass, R)
% Class of the incumbent displaced by an arrival of class ARRCLASS at node JND,
% drawn in proportion to the servers' class occupancies. NVEC already counts
% the arrival, so it is discounted here to recover the pre-arrival in-service
% composition (in-service = population minus buffer occupancy).
base = (jnd-1)*R;
insvc = zeros(1,R);
for r = 1:R
    insvc(r) = nvec(base + r) - sum(buf == r);
    if r == arrClass
        insvc(r) = insvc(r) - 1; % discount the job that just arrived
    end
end
insvc(insvc < 0) = 0;
tot = sum(insvc);
if tot <= 0
    c = 0;
    return
end
u = rand * tot;
acc = 0;
c = find(insvc > 0, 1, 'last');
for r = 1:R
    acc = acc + insvc(r);
    if insvc(r) > 0 && u < acc
        c = r;
        return
    end
end
end

function npop = classCounts(X, phOff, nph, ind, R)
% Per-class populations at node IND, summing each class over its phases. The
% scheduling rate laws are class-level: they are unchanged by phase expansion,
% and only the per-phase share (see kirFrac) is layered on top.
npop = zeros(R,1);
for r = 1:R
    npop(r) = sum(X((phOff(ind,r)+1):(phOff(ind,r)+nph(ind,r))));
end
end

function n = classPop(X, phOff, nph, ind, r)
% Population of class R at node IND, summed over its phases.
n = sum(X((phOff(ind,r)+1):(phOff(ind,r)+nph(ind,r))));
end

function f = kirFrac(X, slot, phOff, nph, ind, r)
% Share of its class that the job population in one phase represents: kir/nir.
% The class-level rate law is split across the class's phases in this ratio,
% which is exactly how State.afterEventStation writes every phase-aware case
% (e.g. DPS uses (kir/nir) * [class share]). For a single-phase class this is
% 1 whenever the class is present, so an exponential model is unaffected.
nir = sum(X((phOff(ind,r)+1):(phOff(ind,r)+nph(ind,r))));
if nir <= 0
    f = 0;
else
    f = X(slot) / nir;
end
end

function pentry = entryProbs(sn, jnd, s, nphjs)
% Entry-phase distribution of a class-s job arriving at node JND: pie of its
% service process there. A non-station node, or a station whose process is
% absent (a disabled class), has a single phase entered with probability 1.
pentry = zeros(1, nphjs);
if nphjs <= 1
    pentry(1) = 1;
    return
end
ist = sn.nodeToStation(jnd);
p = sn.pie{ist}{s};
p = p(:)';
if isempty(p) || all(isnan(p)) || sum(p) <= 0
    % no entry distribution declared: enter the first phase
    pentry(1) = 1;
    return
end
pentry(1:min(nphjs,numel(p))) = p(1:min(nphjs,numel(p)));
pentry = pentry / sum(pentry);
end

function tf = isBuffered(ind, sn)
% True for stations whose waiting jobs are held in an ordered buffer.
tf = false;
if sn.isstation(ind)
    ist = sn.nodeToStation(ind);
    tf = any(sn.sched(ist) == [SchedStrategy.FCFS, SchedStrategy.LCFS, ...
        SchedStrategy.SIRO, SchedStrategy.HOL, SchedStrategy.SEPT, SchedStrategy.LEPT, ...
        SchedStrategy.LCFSPR]);
end
end

function f = lldfac(ldrow, ntot, lldlimit)
% Limited load-dependent scaling factor at total station population NTOT. Returns
% 1 when the station has no load dependence or is empty; otherwise the tabulated
% factor, clamped to the last entry beyond the tabulated limit.
if isempty(ldrow) || ntot < 1
    f = 1;
else
    f = ldrow(min(round(ntot), lldlimit));
end
end

% ======================================================================
% Round-robin routing pointers
% ======================================================================

function rr = rrPrecompute(sn, state)
% Per-(node,class) round-robin pointers, seeded from the initial state. RROBIN
% stores the destination node index in its slot, WRROBIN a POSITION in the
% weighted cycle (each outlink replicated by its weight), matching
% State.fromMarginal and State.afterEventRouter.
rr = struct('on', false);
if ~isfield(sn,'routing') || isempty(sn.routing)
    return
end
if ~any(sn.routing(:) == RoutingStrategy.RROBIN | sn.routing(:) == RoutingStrategy.WRROBIN)
    return
end
R = sn.nclasses;
rr.on = true;
rr.isrr = false(sn.nnodes, R);
rr.iswrr = false(sn.nnodes, R);
rr.cycle = cell(sn.nnodes, R);   % ordered destination list walked per dispatch
rr.pos = ones(sn.nnodes, R);     % current position in that list
for ind = 1:sn.nnodes
    for r = 1:R
        isRR = sn.routing(ind,r) == RoutingStrategy.RROBIN;
        isWRR = sn.routing(ind,r) == RoutingStrategy.WRROBIN;
        if ~isRR && ~isWRR
            continue
        end
        np = sn.nodeparam{ind}{r};
        if isWRR && isfield(np,'weighted_outlinks') && ~isempty(np.weighted_outlinks)
            cyc = np.weighted_outlinks;
        else
            cyc = np.outlinks;
        end
        rr.isrr(ind,r) = true;
        rr.iswrr(ind,r) = isWRR;
        rr.cycle{ind,r} = cyc(:)';
        % seed the pointer from the initial state slot so a warm start is honored
        p0 = 1;
        if sn.isstateful(ind)
            st = state{sn.nodeToStateful(ind)};
            slot = sum(sn.nvars(ind, 1:(R + r)));
            if ~isempty(st) && slot >= 1 && slot <= numel(st)
                v = st(1, slot);
                if isWRR
                    if v >= 1 && v <= numel(cyc), p0 = v; end
                else
                    j = find(cyc == v, 1);
                    if ~isempty(j), p0 = j; end
                end
            end
        end
        rr.pos(ind,r) = p0;
    end
end
end

function [rr, jnd] = rrNext(rr, ind, r)
% Advance the pointer cyclically and return the destination it lands on.
cyc = rr.cycle{ind,r};
p = rr.pos(ind,r);
if p >= numel(cyc)
    p = 1;
else
    p = p + 1;
end
rr.pos(ind,r) = p;
jnd = cyc(p);
end

% ======================================================================
% G-network signals
% ======================================================================

function sig = signalPrecompute(sn)
% Signal classes and their removal parameters. A signal never joins a station:
% it removes jobs there and is annihilated (State.afterEventStationSignal).
sig = struct('on', false);
if ~isfield(sn,'issignal') || isempty(sn.issignal) || ~any(sn.issignal)
    return
end
sig.on = true;
sig.sn = sn;                       % signalBatchPMF and isCatastropheSignal need it
sig.issignal = logical(sn.issignal(:)');
sig.nonsignal = find(~sig.issignal);
end

function tf = sigIsSignalArrival(sig, destPos, R, smap)
% True when the state slot DESTPOS is a signal class at a station.
s = smap.class(destPos);
jnd = smap.node(destPos);
tf = sig.issignal(s) && sig.sn.isstation(jnd);
end

function [nvec, buffers] = sigApply(sig, nvec, buffers, destPos, R, mi, smap)
% Apply the arrival of a signal class at a station: pick the victims and remove
% them. Mirrors State.afterEventStationSignal, sampled instead of enumerated.
sn = sig.sn;
jnd = smap.node(destPos);
cls = smap.class(destPos);
base = smap.phOff(jnd,1);   % first slot of this node
ist = sn.nodeToStation(jnd);

% CATASTROPHE empties the station of every job, ignoring the batch-size
% distribution: a catastrophe removes all jobs by definition.
if State.isCatastropheSignal(sn, cls)
    nvec((base+1):(base+R)) = 0;
    buffers{jnd} = [];
    return
end

% Eligible victim classes. A signal that declares a target (forJobClass,
% sn.signaltarget >= 1) only removes that class; otherwise every non-signal
% class is eligible, which is the classic Gelenbe negative customer and is what
% SolverMAM and SolverLDES both do.
tgt = -1;
if isfield(sn,'signaltarget') && ~isempty(sn.signaltarget) && numel(sn.signaltarget) >= cls
    tgt = sn.signaltarget(cls);
end
if tgt >= 1
    tgtclasses = tgt;
else
    tgtclasses = sig.nonsignal;
end
tgtclasses = tgtclasses(nvec(base + tgtclasses) > 0);
ntot = sum(nvec(base + tgtclasses));
if isempty(tgtclasses) || ntot <= 0
    return % no victim: the signal simply vanishes
end

% Batch size, drawn from the pmf the reference enumerates. It is already
% clipped at the eligible population, so an oversized batch empties it rather
% than driving the queue negative.
[kvals, kprobs] = State.signalBatchPMF(sn, cls, ntot);
k = kvals(find(cumsum(kprobs) >= rand, 1));
if isempty(k)
    k = kvals(end);
end

policy = RemovalPolicy.RANDOM;
if isfield(sn,'signalrempolicy') && ~isempty(sn.signalrempolicy) && numel(sn.signalrempolicy) >= cls
    policy = sn.signalrempolicy(cls);
end

for step = 1:k
    [nvec, buffers, removed] = sigRemoveOne(sn, nvec, buffers, jnd, ist, base, ...
        tgtclasses, policy, R, mi);
    if ~removed
        break % already drained
    end
end
end

function [nvec, buffers, removed] = sigRemoveOne(sn, nvec, buffers, jnd, ist, base, tgtclasses, policy, R, mi)
% Remove one victim under the signal's removal policy. Waiting jobs live in the
% buffer; the rest of each class population is in service.
removed = false;
buf = buffers{jnd};
waitIdx = find(ismember(buf, tgtclasses));      % eligible waiting positions
nwait = numel(waitIdx);
nsrv = 0;
for r = tgtclasses
    nsrv = nsrv + max(0, nvec(base + r) - sum(buf == r));
end
if nwait == 0 && nsrv == 0
    return
end

% FCFS/LCFS rank the waiting line by age, which only an ordered buffer records.
% The NRM buffer is newest-first / oldest-last, so the head of line (oldest) is
% the last eligible position and the most recent arrival the first. A per-class
% count buffer (SIRO/SEPT/LEPT) carries no age, so an age-based policy
% degenerates to a uniform draw there, exactly as in the reference.
isOrdered = any(sn.sched(ist) == [SchedStrategy.FCFS, SchedStrategy.HOL, SchedStrategy.LCFS]);
ageOrdered = isOrdered && (policy == RemovalPolicy.FCFS || policy == RemovalPolicy.LCFS);
if ageOrdered && nwait > 0
    if policy == RemovalPolicy.FCFS
        pick = waitIdx(end);   % head of line: the oldest waiting job
    else
        pick = waitIdx(1);     % the most recent arrival
    end
    victim = buf(pick);
    buffers{jnd}(pick) = [];
    nvec(base + victim) = nvec(base + victim) - 1;
    removed = true;
    return
end

% RANDOM draws uniformly over waiting and in-service alike; FCFS/LCFS drain the
% waiting line before reaching into the servers.
if policy == RemovalPolicy.RANDOM
    total = nwait + nsrv;
else
    total = nwait;
    if total == 0
        total = nsrv;
    end
end
u = rand * total;
if nwait > 0 && (policy ~= RemovalPolicy.RANDOM || u < nwait)
    % a waiting victim, uniform over the eligible positions
    pick = waitIdx(1 + floor(rand * nwait));
    victim = buf(pick);
    buffers{jnd}(pick) = [];
    nvec(base + victim) = nvec(base + victim) - 1;
    removed = true;
    return
end

% an in-service victim, uniform over the eligible in-service jobs
acc = 0;
target = rand * nsrv;
for r = tgtclasses
    cnt = max(0, nvec(base + r) - sum(buf == r));
    acc = acc + cnt;
    if cnt > 0 && target < acc
        nvec(base + r) = nvec(base + r) - 1;
        removed = true;
        % the freed server pulls the head of line in, which in the NRM is just
        % the waiting job leaving the buffer (in-service is derived as
        % population minus buffer occupancy)
        total_new = sum(nvec((base+1):(base+R)));
        if numel(buffers{jnd}) > max(0, total_new - mi(jnd))
            buffers{jnd}(end) = [];   % head of line: the oldest waiting job
        end
        return
    end
end
end

% ======================================================================
% Balking
% ======================================================================

function balk = balkPrecompute(sn)
% Per (station,class) balking threshold table, for the QUEUE_LENGTH strategy.
balk = struct('on', false);
if ~isfield(sn,'balkingStrategy') || isempty(sn.balkingStrategy)
    return
end
if ~any(sn.balkingStrategy(:) == BalkingStrategy.QUEUE_LENGTH)
    return
end
balk.on = true;
balk.strategy = sn.balkingStrategy;
balk.thresholds = sn.balkingThresholds;
% index maps needed by balkDraw, captured so it never needs the whole sn
balk.isstation = sn.isstation;
balk.nodeToStation = sn.nodeToStation;
end

function tf = balkDraw(balk, nvec, destPos, R, smap)
% True if the job routed to state slot DESTPOS balks. The threshold table is
% scanned in order and the FIRST interval containing the pre-arrival total
% station population wins, matching State.afterEventStation.
tf = false;
jnd = smap.node(destPos);
s   = smap.class(destPos);
if ~balk.isstation(jnd)
    return
end
ist = balk.nodeToStation(jnd);
if ist < 1 || balk.strategy(ist, s) ~= BalkingStrategy.QUEUE_LENGTH
    return
end
qlen = sum(classCounts(nvec, smap.phOff, smap.nph, jnd, R)); % pre-arrival total population
th = balk.thresholds{ist, s};
balkProb = 0;
for ti = 1:numel(th)
    t = th{ti};
    if qlen >= t{1} && qlen <= t{2}
        balkProb = t{3};
        break
    end
end
tf = balkProb > 0 && rand < balkProb;
end

function tf = capacityLoss(sn, nvec, destPos, R, smap)
% True if an OPEN-class job routed to state slot DESTPOS is lost because its
% destination station is a physically finite-capacity station that is already
% full. Mirrors the hasRoom gate + State.arrivalIsLost of afterEventStation:
% total occupancy (buffer + in service) is capped at sn.cap, per class at
% sn.classcap (0 = no per-class bound). A refused CLOSED job must block, not
% vanish from the conserved population, so it is NOT dropped here. Inert unless
% the destination declares a physical drop rule. This is the finite-capacity
% loss the NRM reaction network otherwise omits, which let a capped queue
% overflow well past sn.cap under simulation.
tf = false;
jnd = smap.node(destPos);
dstC = smap.class(destPos);
if ~sn.isstation(jnd)
    return
end
ist = sn.nodeToStation(jnd);
if ist < 1
    return
end
if ~State.isPhysicalCapacity(sn, ist, dstC) || ~State.arrivalIsLost(sn, ist, dstC)
    return
end
cc = classCounts(nvec, smap.phOff, smap.nph, jnd, R); % pre-arrival populations
capLimit = sn.cap(ist);
if isfinite(capLimit) && sum(cc) >= capLimit
    tf = true;
    return
end
if ~isempty(sn.classcap) && size(sn.classcap,1) >= ist && size(sn.classcap,2) >= dstC
    classCapLimit = sn.classcap(ist, dstC);
    if classCapLimit > 0 && cc(dstC) >= classCapLimit
        tf = true;
    end
end
end

% ======================================================================
% Destination-side BLOCKING: a refused arrival that may not be dropped
% ======================================================================

function blk = blockPrecompute(sn, R)
% Stations and classes at which a refused arrival must BLOCK rather than be
% lost. CAPACITYLOSS above answers the OPEN half of the same question, and
% returns false for a closed class precisely because a closed network's
% population is an invariant; nothing then stopped the reaction, so the NRM
% fired into the full station anyway and reported the UNCONSTRAINED answer
% (BUG-81). On a closed 3-queue tandem, N=6, Q2 capped at 2 it gave QLen
% [1.96 2.07 1.97] against the exact [3.609 0.971 1.420] -- a mean of 2.07 at a
% station that holds 2 -- while SOLVER_SSA, whose producer already implements
% the contract, gave the exact answer.
%
% Blocking is the THIRD outcome of a firing. The departure does not occur, the
% source is not decremented, no buffer moves; only the reaction's own clock is
% redrawn, which is exact by memorylessness (the residual of an exponential, or
% of the current PH phase, is that same exponential). It is what SOLVER_CTMC
% does when MATCHROW cannot find the over-capacity target and drops the arc,
% which is why the two now agree.
%
% The gate is deliberately NARROWER than State.afterEventStation's: it fires
% only for a CLOSED class, so every open-class drop path -- M/M/1/K included --
% keeps the sample path it had. A cap that cannot bind (the usual cap = N
% default) never enters BLK.CAN, so BLK.ON stays false and the mechanism costs
% nothing on models that do not need it.
M = sn.nstations;
blk.can = false(M, R);
blk.cap = inf(M, 1);
blk.ccap = inf(M, R);
blk.node2st = sn.nodeToStation(:);
if ~isempty(sn.cap)
    capv = sn.cap(:);
    blk.cap(1:min(M,numel(capv))) = capv(1:min(M,numel(capv)));
end
if ~isempty(sn.classcap)
    cc = sn.classcap;
    rows = min(M, size(cc,1));
    cols = min(R, size(cc,2));
    % a non-positive class cap is "unset", not "holds nothing"
    sub = cc(1:rows, 1:cols);
    sub(sub <= 0) = inf;
    blk.ccap(1:rows, 1:cols) = sub;
end
njobs = sn.njobs(:);
% A cap that CANNOT BIND is not a blocking site. Every closed model carries
% cap(ist) = N by default, and a station that can hold the whole population never
% refuses one: the arriving job is itself one of the N, so the pre-arrival count
% is at most N-1. Excluding those is what keeps BLK.ON false -- and the
% per-firing test unpaid -- on ordinary models. An open class present anywhere
% makes a finite station cap binding again, since its jobs are not counted in N.
closedTotal = sum(njobs(isfinite(njobs)));
anyOpen = any(isinf(njobs));
for ist = 1:M
    for r = 1:R
        if r > numel(njobs) || isinf(njobs(r))
            continue % open class: refusal LOSES, handled by capacityLoss
        end
        bindsSt = isfinite(blk.cap(ist)) && (anyOpen || blk.cap(ist) < closedTotal);
        bindsCl = isfinite(blk.ccap(ist, r)) && blk.ccap(ist, r) < njobs(r);
        if bindsSt || bindsCl
            blk.can(ist, r) = true;
        end
    end
end
blk.on = any(blk.can(:));
end

function tf = capacityBlock(blk, nvec, destPos, srcPos, R, smap)
% True when a class-r job routed to state slot DESTPOS cannot be admitted and
% the firing must be cancelled.
%
% The population read is the PRE-arrival one, minus the departing job when it
% currently sits at the destination node: a self-loop or a feedback arc at a
% station already at cap would otherwise block itself forever, while the
% reference producer sees the state AFTER the departure half. FCRREFUSINGREGION
% discounts its source in the same region for the same reason.
tf = false;
if ~blk.on
    return
end
jnd = smap.node(destPos);
if jnd > numel(blk.node2st)
    return
end
ist = blk.node2st(jnd);
if ist < 1
    return
end
r = smap.class(destPos);
if ist > size(blk.can,1) || r > size(blk.can,2) || ~blk.can(ist, r)
    return
end
sameNode = srcPos >= 1 && srcPos <= numel(smap.node) && smap.node(srcPos) == jnd;
cc = classCounts(nvec, smap.phOff, smap.nph, jnd, R);
if isfinite(blk.cap(ist))
    total = sum(cc);
    if sameNode
        total = total - 1;
    end
    if total >= blk.cap(ist)
        tf = true;
        return
    end
end
if isfinite(blk.ccap(ist, r))
    pop = cc(r);
    if sameNode && smap.class(srcPos) == r
        pop = pop - 1;
    end
    if pop >= blk.ccap(ist, r)
        tf = true;
    end
end
end

% ======================================================================
% Finite capacity regions (DROP rule)
% ======================================================================

function fcr = fcrPrecompute(sn)
% Per-region member stations and admission caps, mirroring the FCR precompute
% of SOLVER_SSA (the serial engine) field for field.
fcr = struct('on', false);
if ~isfield(sn,'nregions') || sn.nregions == 0
    return
end
K = sn.nclasses;
F = sn.nregions;
fcr.on = true;
fcr.memberMask = false(F, sn.nstations);
fcr.classCap = cell(F,1);
fcr.globalCap = inf(F,1);
fcr.memCap = inf(F,1);
fcr.sz = cell(F,1);
fcr.A = cell(F,1);
fcr.b = cell(F,1);
% Per-(region,class) admission rule: DROP destroys a refused job, WAITQ parks
% it in the region FIFO and admits it head-of-line as capacity frees. Mirrors
% SOLVER_SSA's fcrRule (regionrule ~= DropStrategy.DROP). Regions with no
% WAITQ class carry no FIFO, so pure-DROP models pay nothing.
fcr.waitq = false(F, K);
if isfield(sn,'regionrule') && ~isempty(sn.regionrule)
    for f = 1:F
        for r = 1:K
            fcr.waitq(f,r) = sn.regionrule(f,r) ~= DropStrategy.DROP;
        end
    end
end
fcr.anyWaitq = any(fcr.waitq(:));
for f = 1:F
    Rmat = sn.region{f};                          % M x (K+1)
    % membership: any job-count cap OR the region memory budget set on the
    % station row (a memory-only region has all job-count entries at -1)
    memvec = -ones(sn.nstations,1);
    if isfield(sn,'regionmaxmem') && numel(sn.regionmaxmem) >= f && ~isempty(sn.regionmaxmem{f})
        memvec = sn.regionmaxmem{f}(:);
    end
    mask = sn_region_members(sn, f, Rmat, memvec);
    fcr.memberMask(f, 1:numel(mask)) = mask;
    members = find(mask);
    ccap = inf(1,K);
    for r = 1:K
        cv = Rmat(members, r); cv = cv(cv ~= -1);
        if ~isempty(cv); ccap(r) = min(cv); end
    end
    fcr.classCap{f} = ccap;
    gv = Rmat(members, K+1); gv = gv(gv ~= -1);
    if ~isempty(gv); fcr.globalCap(f) = min(gv); end
    if isfield(sn,'regionmaxmem') && numel(sn.regionmaxmem) >= f && ~isempty(sn.regionmaxmem{f})
        mv = sn.regionmaxmem{f}(members); mv = mv(mv ~= -1);
        if ~isempty(mv); fcr.memCap(f) = min(mv); end
    end
    fcr.sz{f} = sn.regionsz(f,:);
    if isfield(sn,'regionlincon') && size(sn.regionlincon,1) >= f && ~isempty(sn.regionlincon{f,1})
        fcr.A{f} = sn.regionlincon{f,1};
        fcr.b{f} = sn.regionlincon{f,2};
    end
end
% node-level membership, so the gate can be evaluated straight off the NRM
% state vector without going through station indices on every firing
fcr.memberNode = false(F, sn.nnodes);
for f = 1:F
    for ist = find(fcr.memberMask(f,:))
        fcr.memberNode(f, sn.stationToNode(ist)) = true;
    end
end
end

function tf = fcrViolates(xn, ccap, gcap, memcap, sz, A, b)
% True if per-class population vector XN breaks any admission constraint of
% the region. Mirrors fcr_violates in SOLVER_SSA.
tf = any(xn > ccap) || sum(xn) > gcap || (xn * sz(:) > memcap);
if ~tf && ~isempty(A)
    tf = any(A * xn(:) > b(:));
end
end

function x = fcrRegionPop(nvec, memberNodeRow, R, smap)
% Per-class population of a region, read directly off the NRM state vector.
x = zeros(1, R);
for jnd = find(memberNodeRow)
    x = x + classCounts(nvec, smap.phOff, smap.nph, jnd, R)';
end
end

function tf = fcrAdmits(fcr, nvec, srcNode, srcClass, dstNode, dstClass, R, smap)
% True if a class-DSTCLASS job may enter node DSTNODE, having just left node
% SRCNODE as class SRCCLASS. Only regions containing the destination can
% refuse the move; a move whose source is in the same region frees a slot
% first, so the departure is accounted for before the arrival is tested.
tf = fcrRefusingRegion(fcr, nvec, srcNode, srcClass, dstNode, dstClass, R, smap) == 0;
end

function f = fcrRefusingRegion(fcr, nvec, srcNode, srcClass, dstNode, dstClass, R, smap)
% Index of the FIRST region that refuses a class-DSTCLASS job entering
% DSTNODE, having just left SRCNODE as class SRCCLASS; 0 if every region
% admits it. Same admission test as the DROP path, but it names the refusing
% region so the caller can consult that region's DROP/WAITQ rule. Mirrors the
% first-region `break` of SOLVER_SSA's blockFCR loop.
f = 0;
if ~fcr.on
    return
end
for ff = 1:size(fcr.memberNode,1)
    if ~fcr.memberNode(ff, dstNode)
        continue % this region does not constrain the destination
    end
    x = fcrRegionPop(nvec, fcr.memberNode(ff,:), R, smap);
    % srcNode <= 0 means the mover has no live source in the state (a WAITQ
    % release, whose job already left its source when it was parked), so no
    % source slot is freed.
    if srcNode > 0 && fcr.memberNode(ff, srcNode)
        x(srcClass) = x(srcClass) - 1;
    end
    x(dstClass) = x(dstClass) + 1;
    if fcrViolates(x, fcr.classCap{ff}, fcr.globalCap(ff), fcr.memCap(ff), ...
            fcr.sz{ff}, fcr.A{ff}, fcr.b{ff})
        f = ff;
        return
    end
end
end

function [nvec, buffers, fcrBuf, released, svcph, svcChanged, startCount, preemptCount] = fcrReleaseCascade(fcr, nvec, buffers, fcrBuf, mi, R, sn, smap, svcph, bufPHNode, startCount, preemptCount)
% Strict-FIFO head-of-line release of parked WAITQ tokens: admit each region's
% FIFO head while the admission constraints permit, applying the arrival to
% the destination station (entry-phase slot plus buffer join). Mirrors
% SOLVER_SSA's fcr_release. A token is (dstNode, dstClass); the phase is drawn
% at release, as a routed arrival draws it. Loops until a full pass frees
% nothing, so a release that frees capacity elsewhere cascades.
released = 0;
svcChanged = false;
progress = true;
while progress
    progress = false;
    for f = 1:numel(fcrBuf)
        if isempty(fcrBuf{f})
            continue
        end
        tok = fcrBuf{f}(1);
        dstNode = floor((tok-1)/R) + 1;
        dstClass = mod(tok-1, R) + 1;
        % The parked job already left its source, so admission is tested with
        % the source term absent (srcNode = -1 never matches memberNode).
        if fcrRefusingRegion(fcr, nvec, -1, dstClass, dstNode, dstClass, R, smap) ~= 0
            continue % head-of-line: this FIFO stays blocked
        end
        if bufPHNode(dstNode)
            % Buffered-PH destination: the released job lands in the class total
            % slot; whether it enters service (and its entry phase) is decided in
            % applyArrivalBuffer against the server occupancy, exactly as a routed
            % arrival is.
            nvec(smap.phOff(dstNode,dstClass) + 1) = nvec(smap.phOff(dstNode,dstClass) + 1) + 1;
            [buffers, svcph, arrCh, startCount, preemptCount] = applyArrivalBuffer(dstNode, dstClass, nvec, buffers, mi, R, sn, smap, svcph, bufPHNode, startCount, preemptCount);
            svcChanged = svcChanged || arrCh;
        else
            pentry = entryProbs(sn, dstNode, dstClass, smap.nph(dstNode,dstClass));
            ke = drawFromDist(pentry);
            dslot = smap.phOff(dstNode,dstClass) + ke;
            nvec(dslot) = nvec(dslot) + 1;
            [buffers, svcph, arrCh, startCount, preemptCount] = applyArrivalBuffer(dstNode, dstClass, nvec, buffers, mi, R, sn, smap, svcph, bufPHNode, startCount, preemptCount);
            svcChanged = svcChanged || arrCh;
        end
        fcrBuf{f}(1) = [];
        released = released + 1;
        progress = true;
    end
end
end

function ke = drawFromDist(p)
% Index drawn from the (unnormalized, nonnegative) weight vector P.
tot = sum(p);
if tot <= 0
    ke = 1;
    return
end
c = cumsum(p) / tot;
ke = find(c > rand, 1);
if isempty(ke)
    ke = numel(p);
end
end

function [outClass, var, category, released] = cacheAccess(sn, ind, class, var)
% Simulate one cache READ at cache node IND by a class-CLASS job over the cache
% state VAR (totalCacheCapacity content slots followed, when a retrieval system
% is present, by block A, a per-item retrieval-occupancy bitmap, and by block B,
% the per-retrieval-class count of requests merged onto an in-flight fetch).
% Returns the class the job leaves in -- OUTCLASS = 0 means the request merged
% onto a pending fetch and produces nothing yet -- the rewritten VAR, a CATEGORY
% (1 hit, 2 miss/retrieval-complete, 3 delayed-hit, 4 begin-retrieval), and
% RELEASED, the (hitClass, count) rows freed by a completing fetch. A faithful
% port of State.afterEventCache (READ, isSimulation): non-retrieval hit/miss with
% all replacement policies, plus the retrieval system where a miss for an item
% not yet being fetched begins a retrieval (switch to the item's retrieval class,
% mark block A), a concurrent request for an item already being fetched is held
% in block B as a delayed hit, and a returning retrieval-class read completes the
% miss (clear block A, admit the item, release the merged requests).
np = sn.nodeparam{ind};
m = np.itemcap;
ac = np.accost;
h = length(m);
released = zeros(0,2);
[rcList, rcItems, rcOrigClass] = State.cacheRetrievalClassMap(sn, ind);
replacement_id = np.replacestrat;
if isfield(np,'totalCacheCapacity') && ~isempty(np.totalCacheCapacity)
    totalCacheCapacity = np.totalCacheCapacity;
else
    totalCacheCapacity = sum(m);
end
hitclassArr = np.hitclass;
missclassArr = np.missclass;
if isfield(np,'retrievalClassIndices') && ~isempty(np.retrievalClassIndices)
    rci = np.retrievalClassIndices(:)';
else
    rci = [];
end
isFromRetrieval = any(rci == class);
if isfield(np,'retrievalClasses') && ~isempty(np.retrievalClasses)
    retrClasses = np.retrievalClasses;
else
    retrClasses = [];
end
hasRetrieval = isfield(np,'retrievalSystemCapacity') && ~isempty(np.retrievalSystemCapacity) ...
    && any(np.retrievalSystemCapacity > 0);

p = np.pread{class};
k = drawFromDist(p);                         % requested item
l = drawFromDist(ac{class,k}(1,:));          % target list for a miss (1 => reject)
posk = find(k == var(1:totalCacheCapacity), 1, 'first');
if isFromRetrieval
    posk = [];   % a returning retrieval always COMPLETES its own miss
end

if ~isempty(posk)
    % ===================== CACHE HIT =====================
    outClass = hitclassArr(class);
    category = 1;
    if posk <= sum(m(1:h-1))
        % hit in list i < h: promote toward the last list
        i = find(posk <= cumsum(m), 1);
        j = posk - sum(m(1:i-1));
        accrow = ac{class,k}(1+i, (1+i):end);
        inew = i + drawFromDist(accrow / sum(accrow)) - 1;
        switch replacement_id
            case ReplacementStrategy.FIFO
                if inew ~= i
                    varp = var;
                    varp(cpos(i,j)) = var(cpos(inew,m(inew)));
                    varp(cpos(inew,2):cpos(inew,m(inew))) = var(cpos(inew,1):cpos(inew,m(inew)-1));
                    varp(cpos(inew,1)) = k;
                    var = varp;
                end
            case ReplacementStrategy.RR
                varp = var;
                rpos = randi(m(inew),1,1);
                varp(cpos(i,j)) = var(cpos(inew,rpos));
                varp(cpos(inew,rpos)) = k;
                var = varp;
            case {ReplacementStrategy.LRU, ReplacementStrategy.SFIFO, ...
                    ReplacementStrategy.HLRU, ReplacementStrategy.QLRU}
                varp = var;
                varp(cpos(i,2):cpos(i,j)) = var(cpos(i,1):cpos(i,j-1));
                varp(cpos(i,1)) = var(cpos(inew,m(inew)));
                varp(cpos(inew,2):cpos(inew,m(inew))) = var(cpos(inew,1):cpos(inew,m(inew)-1));
                varp(cpos(inew,1)) = k;
                var = varp;
        end
    else
        % hit in the last list h
        j = posk - sum(m(1:h-1));
        switch replacement_id
            case {ReplacementStrategy.RR, ReplacementStrategy.FIFO, ReplacementStrategy.SFIFO}
                % no reordering
            case {ReplacementStrategy.LRU, ReplacementStrategy.HLRU, ReplacementStrategy.QLRU}
                varp = var;
                varp(cpos(h,2):cpos(h,j)) = var(cpos(h,1):cpos(h,j-1));
                varp(cpos(h,1)) = var(cpos(h,j));
                var = varp;
        end
    end
    return
end

% ===================== CACHE MISS / retrieval =====================
if hasRetrieval && ~isFromRetrieval
    % Consult the retrieval system: an item with a retrieval class is fetched
    % rather than admitted directly on a miss.
    rClass = -1;
    if ~isempty(retrClasses) && k <= size(retrClasses,1) && class <= size(retrClasses,2)
        rClass = retrClasses(k, class);
    end
    if rClass ~= -1
        inRetrieval = (totalCacheCapacity + k <= numel(var)) && var(totalCacheCapacity + k) ~= 0;
        if inRetrieval
            % DELAYED HIT: merges onto the in-flight fetch and is held in block B
            % until it completes, then released in its own hit class.
            bslot = find(rcList == rClass, 1);
            bcol = totalCacheCapacity + np.nitems + bslot;
            if ~isempty(bslot) && bcol <= numel(var)
                var(bcol) = var(bcol) + 1;
            end
            outClass = 0;
            category = 3;
            return
        else
            % BEGIN retrieval: switch to the item's retrieval class and mark the
            % item as being fetched; the job routes to the retrieval queue and
            % returns later to complete the miss.
            var(totalCacheCapacity + k) = 1;
            outClass = rClass;
            category = 4;
            return
        end
    end
end

% COMPLETE the miss: a returning retrieval, or a plain miss with no retrieval
% class. Clear the retrieval bit (if any) and admit item k per the policy.
if isFromRetrieval && (totalCacheCapacity + k <= numel(var))
    var(totalCacheCapacity + k) = 0;
    % Every request merged onto this fetch is released now as a delayed hit, in
    % the hit class of the job class that issued it.
    for bslot = 1:numel(rcList)
        if rcItems(bslot) ~= k
            continue
        end
        bcol = totalCacheCapacity + np.nitems + bslot;
        if bcol <= numel(var) && var(bcol) > 0
            hc = hitclassArr(rcOrigClass(bslot));
            if hc > 0
                released(end+1,:) = [hc, var(bcol)]; %#ok<AGROW>
            end
            var(bcol) = 0;
        end
    end
end
outClass = missclassArr(class);
category = 2;
listidx = l - 1;
switch replacement_id
    case {ReplacementStrategy.FIFO, ReplacementStrategy.LRU, ...
            ReplacementStrategy.SFIFO, ReplacementStrategy.HLRU}
        if listidx > 0
            varp = var;
            varp(cpos(listidx,2):cpos(listidx,m(listidx))) = var(cpos(listidx,1):cpos(listidx,m(listidx)-1));
            varp(cpos(listidx,1)) = k;
            var = varp;
        end
    case ReplacementStrategy.RR
        if listidx > 0
            rpos = randi(m(listidx),1,1);
            var(cpos(listidx,rpos)) = k;
        end
    case ReplacementStrategy.QLRU
        if isfield(np,'qlru') && ~isempty(np.qlru), qadm = np.qlru; else, qadm = 1.0; end
        if listidx > 0 && rand <= qadm
            varp = var;
            varp(cpos(listidx,2):cpos(listidx,m(listidx))) = var(cpos(listidx,1):cpos(listidx,m(listidx)-1));
            varp(cpos(listidx,1)) = k;
            var = varp;
        end
end

    function pos = cpos(ii,jj)
        pos = sum(m(1:ii-1)) + jj;
    end
end

% ======================================================================
% PS-family sharing factors
%
% Each returns the multiplier applied to the class-r service rate, i.e. the
% fraction of total service capacity that class r receives in population
% state NVECPOP. All mirror the corresponding case of
% State.afterEventStation specialized to exponential (single-phase) service,
% where the phase population kir equals the class population nir.
% ======================================================================

function f = dpsshare(w, nvecpop, r)
% DPS: rate_r = mu_r * w_r*n_r / (w.n) on a single server.
den = w(:)' * nvecpop(:);
if den <= 0
    f = 0;
else
    f = w(r) * nvecpop(r) / den;
end
end

function f = gpsshare(w, nvecpop, r)
% GPS: rate_r = mu_r * w_r / (w.c), c_s = 1{n_s>0}, on a single server. The
% weight denominator counts active classes, not jobs, so a class with a
% single job gets the same share as one with many.
if nvecpop(r) <= 0
    f = 0;
    return
end
cir = double(nvecpop(:) > 0);
den = w(:)' * cir;
if den <= 0
    f = 0;
else
    f = w(r) / den;
end
end

function [act, niprio] = prioGroup(nvecpop, r, classprio)
% Population vector restricted to the priority group of class r, and its
% total. Empty classes never define the urgent group.
act = zeros(size(nvecpop));
same = (classprio(:) == classprio(r));
act(same) = nvecpop(same);
niprio = sum(act);
end

function tf = isUrgent(nvecpop, r, classprio)
% True when class r belongs to the most urgent non-empty priority group.
% LINE orders priorities with lower value = more urgent.
occupied = nvecpop(:) > 0;
if ~any(occupied)
    tf = false;
else
    tf = (classprio(r) == min(classprio(occupied)));
end
end

function n = prioPop(nvecpop, r, c, classprio)
% Population that the lld factor is evaluated at: the full station
% population below capacity, the priority-group population above it.
ni = sum(nvecpop);
if ni <= c || ~isUrgent(nvecpop, r, classprio)
    n = ni;
else
    [~, n] = prioGroup(nvecpop, r, classprio);
end
end

function v = prioVec(nvecpop, r, c, classprio)
% Population vector that the cd factor is evaluated at for DPSPRIO/GPSPRIO:
% the priority-restricted vector above capacity, the full one below it.
% Note PSPRIO instead uses the full vector in both branches; that asymmetry
% is inherited from State.afterEventStation and is reproduced here.
ni = sum(nvecpop);
if ni <= c || ~isUrgent(nvecpop, r, classprio)
    v = nvecpop;
else
    v = prioGroup(nvecpop, r, classprio);
end
end

function f = psprioshare(nvecpop, r, c, classprio)
% PSPRIO: PS below capacity; above it only the most urgent non-empty group
% shares the servers and everyone else is frozen.
ni = sum(nvecpop);
if ni <= 0
    f = 0;
elseif ni <= c
    f = (nvecpop(r) / ni) * min(ni, c);
elseif ~isUrgent(nvecpop, r, classprio)
    f = 0;
else
    [~, niprio] = prioGroup(nvecpop, r, classprio);
    f = (nvecpop(r) / niprio) * min(niprio, c);
end
end

function f = dpsprioshare(w, nvecpop, r, c, classprio)
% DPSPRIO: DPS below capacity, DPS restricted to the urgent group above it.
ni = sum(nvecpop);
if ni <= 0
    f = 0;
elseif ni <= c
    f = dpsshare(w, nvecpop, r);
elseif ~isUrgent(nvecpop, r, classprio)
    f = 0;
else
    f = dpsshare(w, prioGroup(nvecpop, r, classprio), r);
end
end

function f = gpsprioshare(w, nvecpop, r, c, classprio)
% GPSPRIO: GPS below capacity, GPS restricted to the urgent group above it.
ni = sum(nvecpop);
if ni <= 0
    f = 0;
elseif ni <= c
    f = gpsshare(w, nvecpop, r);
elseif ~isUrgent(nvecpop, r, classprio)
    f = 0;
else
    f = gpsshare(w, prioGroup(nvecpop, r, classprio), r);
end
end

function f = cdfac(cdbeta, nvecpop, r)
% Class-dependence factor for a class-r completion at a station with per-class
% population vector NVECPOP: the class-r component of the 1xR scaling vector
% returned by the handle CDBETA (see fes_beta_handle and State.cdclassfactor).
% Returns 1 when the station declares no class dependence.
if isempty(cdbeta)
    f = 1;
else
    v = cdbeta(nvecpop(:)');
    f = v(min(r, numel(v)));
end
end


% ======================================================================
% Polling controller helpers
% ======================================================================

function ctrl = pollLandCtrl(pinf, q, mode, budget)
% Controller row [mode, pos, swk, ctr] the server lands in after
% State.pollingNext resolves (q, mode, budget): SERVING q with the visit budget,
% SWITCHING into q with the entry phase drawn from the switchover PH, or PARKED
% at the canonical q. Mirrors State.pollingLand specialized to the single-server
% polling station the NRM carries (exponential service, so no in-service phase).
switch mode
    case 1
        ctrl = [1, q, 0, budget];
    case 2
        swk = drawFromDist(pinf.swpie{q});
        ctrl = [2, q, swk, 0];
    otherwise
        ctrl = [0, q, 0, 0];   % parked
end
end

function g = pollServeGate(ctrl, r)
% 1 when the polling controller CTRL is serving class r, else 0. This is the
% single-server gate that turns a class-r service departure on only while the
% server attends class r.
if numel(ctrl) >= 2 && ctrl(1) == 1 && ctrl(2) == r
    g = 1;
else
    g = 0;
end
end

function rate = pollSwRate(ctrl, pinf)
% Total leaving rate of the switchover phase the controller CTRL currently
% occupies, i.e. -D0(swk,swk) of the switchover PH into buffer pos; 0 unless the
% server is walking (mode SWITCHING). The competition between advancing to
% another phase and absorbing is resolved at firing time by the run loop.
rate = 0;
if numel(ctrl) >= 3 && ctrl(1) == 2
    pos = ctrl(2); swk = ctrl(3);
    D0 = pinf.swD0{pos};
    rate = -D0(swk, swk);
end
end
% ======================================================================
% Stochastic Petri net (Place / Transition) via the Next-Reaction Method
%
% A stochastic Petri net maps onto the reaction network exactly: a Place holds
% a per-class token count (a population slot of the state vector), and a timed
% Transition mode is a reaction whose stoichiometry column is the arc
% incidence -- input (enabling) arcs consume, output (firing) arcs produce.
% Enabling is a propensity gate (all input places at or above their arc weight,
% every inhibitor place strictly below its threshold); a single-server mode
% then fires at its exponential rate, an infinite/k-server mode at that rate
% times its enabling degree. Each firing applies the mode's stoichiometry once
% (consume the input weights, produce the output weights), which is the atomic
% GSPN firing shared by the exact CTMC (single server), JMT and standard GSPN tools.
%
% IMMEDIATE transitions fire in zero time and cannot be an exponential reaction.
% They are resolved by vanishing-marking elimination: after every timed firing
% (and once on the initial marking) every enabled immediate mode is fired,
% highest firing-priority first and, among equal priority, chosen in proportion
% to firing weight, until the marking is tangible (no immediate enabled). The
% timed race only resumes from tangible markings, so the immediate transitions
% never consume simulated time.
%
% Not handled here (rejected upstream by the SSA featset, never reached): a
% Transition whose firing distribution is non-exponential (phase-type or
% general). Representing an in-flight firing's phase needs per-mode phase state
% the reaction network does not carry; the exponential path covers the standard
% GSPN case and every all-exponential validation net (spn_inhibiting,
% spn_twomodes, spn_fourmodes).
% ======================================================================
function [QN, UN, RN, TN, CN, XN] = solver_ssa_nrm_spn(sn, options, phOff, nph, NS, smap)
samples = options.samples;
R = sn.nclasses;
I = sn.nnodes;
M = sn.nstations;
K = sn.nclasses;

% Build the reaction list. Timed modes become reactions (rx); immediate modes
% are collected separately (imm) for the vanishing-marking collapse.
rx = spnEmptyRx();  rx(1) = [];
imm = spnEmptyRx(); imm(1) = [];
% consumers{ind,c}: indices into rx of timed modes that consume from place ind,
% class c. Place throughput is the aggregate firing rate of those modes (once
% per firing, unweighted -- the same depRates the CTMC accumulates from PRE
% events), so this map drives the TN accumulator.
consumers = cell(I, R);
for ind = 1:I
    if sn.nodetype(ind) ~= NodeType.Transition
        continue
    end
    np = sn.nodeparam{ind};
    for m = 1:np.nmodes
        % Marking-dependent firing rates change the propensity with the marking;
        % SolverSSA does not apply the g(marking) multiplier (unlike CTMC and
        % LDES), so reject rather than silently simulate the nominal rate. The
        % sentence is SSA_FIRINGDEP_REFUSAL's, which the gate and the analyzer
        % ask first; this is the enableChecks=false path.
        [fdOk, fdWhy] = ssa_firingdep_refusal(sn);
        if ~fdOk
            line_error(mfilename, fdWhy);
        end
        rec = spnBuildMode(sn, ind, m, phOff, NS);
        if np.timing(m) == TimingStrategy.IMMEDIATE
            imm(end+1) = rec; %#ok<AGROW>
        else
            rx(end+1) = rec; %#ok<AGROW>
            ridx = numel(rx);
            for a = 1:numel(rec.enSlot)
                p = smap.node(rec.enSlot(a));
                c = smap.class(rec.enSlot(a));
                consumers{p, c}(end+1) = ridx;
            end
        end
    end
end

% Source arrivals. A Source is not a Transition, so its Poisson arrival is not
% one of the transition modes above; it needs its own reaction or the fed Place
% stays empty and the net deadlocks. Add one arrival reaction per (Source node,
% open class, routed Place-class edge). Splitting a Poisson stream by the
% independent routing probabilities yields independent Poisson streams, so an
% edge of probability p carries rate lambda*p exactly. The reaction has an EMPTY
% enabling set (always enabled, state-independent propensity = lambda*p) and
% deposits +1 token into the routed Place slot. producers{node,class} indexes
% these so the Source station reports its arrival rate as throughput, which is
% the reference-station throughput of the open class (matching JMT).
producers = cell(I, R);
for ind = 1:I
    if sn.nodetype(ind) ~= NodeType.Source
        continue
    end
    ist = sn.nodeToStation(ind);
    for r = 1:R
        lambda = sn.rates(ist, r);
        if isnan(lambda) || lambda <= 0
            continue
        end
        if sn.procid(ist, r) ~= ProcessType.EXP
            line_error(mfilename, sprintf('Source %s class %d has a non-exponential arrival, which the NRM SPN path does not support; use method=''serial'' or SolverJMT.', sn.nodenames{ind}, r));
        end
        foundPlace = false;
        for jnd = 1:I
            if sn.nodetype(jnd) ~= NodeType.Place
                continue
            end
            for s = 1:R
                p = sn.rtnodes((ind-1)*R + r, (jnd-1)*R + s);
                if p <= 0
                    continue
                end
                foundPlace = true;
                rec = spnEmptyRx();
                rec.node = ind;
                rec.mode = 0;   % arrival, not a transition mode
                Svec = zeros(NS, 1);
                Svec(phOff(jnd, s) + 1) = Svec(phOff(jnd, s) + 1) + 1;
                rec.Svec = Svec;
                rec.enSlot = []; rec.enW = [];
                rec.inhSlot = []; rec.inhThr = [];
                rec.baseRate = lambda * p;   % Poisson thinning by the routing prob
                rec.nservers = 1;            % constant propensity = baseRate
                rec.weight = 1; rec.prio = 1;
                rx(end+1) = rec; %#ok<AGROW>
                producers{ind, r}(end+1) = numel(rx);
            end
        end
        if ~foundPlace
            line_error(mfilename, sprintf('Source %s class %d does not route to any Place; the NRM SPN path needs a Source->Place arc.', sn.nodenames{ind}, r));
        end
    end
end

nR = numel(rx);
if nR == 0
    line_error(mfilename, 'Stochastic Petri net has no timed reaction; nothing to simulate.');
end

% Initial marking: token counts per (place, class), read straight off the
% initial state as the marginal population of each Place.
nvec0 = zeros(NS, 1);
state = sn.state;
for ind = 1:I
    if sn.nodetype(ind) ~= NodeType.Place || ~sn.isstateful(ind)
        continue
    end
    state_i = state{sn.nodeToStateful(ind)};
    [~, nir] = State.toMarginalAggr(sn, ind, state_i);
    for c = 1:R
        if isinf(nir(c))
            line_error(mfilename, 'Infinite marking at a Place is not supported.');
        end
        nvec0(phOff(ind, c) + 1) = nir(c);
    end
end

maxImmSteps = 100000;  % livelock guard for the vanishing-marking collapse

% Finite-capacity Place DROP enforcement. A Place with a finite per-class
% capacity (sn.classcap) or total capacity (sn.cap) loses any arriving token
% that would exceed it (JMT/CTMC loss semantics: an M/M/1/1 Place with cap 1
% holds mean 0.333 at rho=0.5, not the unbounded-M/M/1 value 1.0). Without this
% the deposit nvec+Svec accumulates tokens past capacity. Precompute the per-slot
% per-class caps, the per-place total caps, and each reaction's deposited slots
% so the clamp in the loop touches only what just grew. Mirrors the Python native
% _solver_ssa_nrm_spn.
pcapSlot = inf(NS, 1);            % per-(place,class) slot cap
placeTotalCaps = cell(0, 2);     % {totalCap, slotVec} per capped place
for ind = 1:I
    if sn.nodetype(ind) ~= NodeType.Place || ~sn.isstateful(ind)
        continue
    end
    ist = sn.nodeToStation(ind);
    slotsHere = zeros(1, R);
    for c = 1:R
        slot = phOff(ind, c) + 1;
        slotsHere(c) = slot;
        if ist <= size(sn.classcap, 1)
            cc = sn.classcap(ist, c);
            if isfinite(cc)
                pcapSlot(slot) = cc;
            end
        end
    end
    tcap = Inf;
    if ist <= numel(sn.cap)
        tcap = sn.cap(ist);
    end
    if isfinite(tcap)
        placeTotalCaps(end+1, :) = {tcap, slotsHere}; %#ok<AGROW>
    end
end
hasPlaceCaps = any(isfinite(pcapSlot)) || ~isempty(placeTotalCaps);
depSlots = cell(1, nR);
for k = 1:nR
    depSlots{k} = find(rx(k).Svec > 0);
end

% ---------------------------------------------------------------------
% Next-Reaction Method run loop
% ---------------------------------------------------------------------
nvec = spnCollapse(nvec0, imm, maxImmSteps);
if hasPlaceCaps
    nvec = applyPlaceCaps(nvec, (1:NS)', pcapSlot, placeTotalCaps);
end
Ak = zeros(1, nR);
for k = 1:nR
    Ak(k) = spnProp(nvec, rx(k));
end
Pk = -log(rand(1, nR));
Tk = zeros(1, nR);
tau = (Pk - Tk) ./ Ak;
tau(Ak == 0) = inf;

QN = zeros(M, K); UN = zeros(M, K); RN = zeros(M, K);
TN = zeros(M, K); CN = zeros(1, K); XN = zeros(1, K);
totalTime = 0;
NK = sn.njobs';

n = 1;
while n <= samples
    [dt, kfire] = min(tau);
    if isinf(dt)
        line_error(mfilename, 'Deadlock: no transition is enabled. Quitting nrm method.');
    end
    totalTime = totalTime + dt;

    % Time-average accumulators over the sojourn dt. A Place is an INF station,
    % so its utilization is its mean token count (the SPN convention the CTMC
    % analyzer reports). Its throughput is the summed firing rate of the modes
    % consuming from it.
    for ist = 1:M
        ind = sn.stationToNode(ist);
        for c = 1:K
            tokens = classPop(nvec, phOff, nph, ind, c);
            QN(ist, c) = QN(ist, c) + tokens * dt;
            UN(ist, c) = UN(ist, c) + tokens * dt;
            depr = 0;
            cons = consumers{ind, c};
            for a = 1:numel(cons)
                depr = depr + Ak(cons(a));
            end
            % A Source station has no consuming transition; its throughput is the
            % aggregate arrival rate it injects (producers), so the reference
            % station reports the open-class arrival rate as its throughput.
            prod = producers{ind, c};
            for a = 1:numel(prod)
                depr = depr + Ak(prod(a));
            end
            TN(ist, c) = TN(ist, c) + depr * dt;
        end
    end

    % Fire the selected timed mode (single atomic firing), then collapse any
    % immediate transitions the new marking enabled. A finite-capacity DROP Place
    % loses any token the firing pushed above its capacity, before the immediate
    % cascade sees the new marking.
    nvec = nvec + rx(kfire).Svec;
    if hasPlaceCaps
        nvec = applyPlaceCaps(nvec, depSlots{kfire}, pcapSlot, placeTotalCaps);
    end
    nvec = spnCollapse(nvec, imm, maxImmSteps);

    % Advance the Gibson & Bruck clocks with the pre-firing propensities, then
    % refresh every propensity from the new marking. A firing plus its
    % immediate cascade can change any place, so every reaction is refreshed
    % rather than a dependency subset -- the SPN reaction count is small and
    % this removes any dependency-graph blind spot.
    Tk = Tk + Ak * dt;
    for k = 1:nR
        Ak(k) = spnProp(nvec, rx(k));
    end
    Pk(kfire) = Pk(kfire) - log(rand);
    tau = (Pk - Tk) ./ Ak;
    tau(Ak == 0) = inf;

    n = n + 1;
    if isfield(options, 'verbose') && options.verbose && mod(n, 1e3) == 0 && ~batchStartupOptionUsed
        LineStatus.set('SSA samples: %d', n);
    end
end
LineStatus.close(); % ends the sample-counter row

if totalTime > 0
    QN = QN / totalTime;
    UN = UN / totalTime;
    TN = TN / totalTime;
end
for c = 1:K
    XN(1, c) = TN(sn.refstat(c), c);
    for ist = 1:M
        if TN(ist, c) > 0
            RN(ist, c) = QN(ist, c) / TN(ist, c);
        end
    end
    if XN(1, c) > 0
        CN(1, c) = NK(c) / XN(1, c);
    end
end
QN(isnan(QN)) = 0; UN(isnan(UN)) = 0; RN(isnan(RN)) = 0;
XN(isnan(XN)) = 0; TN(isnan(TN)) = 0; CN(isnan(CN)) = 0;
end

function rec = spnEmptyRx()
% Prototype record for a transition-mode reaction, so struct arrays stay
% homogeneous (MATLAB requires identical fields to concatenate).
rec = struct('node', 0, 'mode', 0, 'Svec', [], 'enSlot', [], 'enW', [], ...
    'inhSlot', [], 'inhThr', [], 'baseRate', 0, 'nservers', 1, ...
    'weight', 1, 'prio', 1);
end

function rec = spnBuildMode(sn, ind, m, phOff, NS)
% Assemble the reaction record of transition IND mode M. Enabling/firing/
% inhibiting are (nnodes x nclasses) matrices; find() gives linear indices
% p+(c-1)*nnodes that decode to the (place, class) whose slot is phOff(p,c)+1.
np = sn.nodeparam{ind};
R = sn.nclasses;
rec = spnEmptyRx();
rec.node = ind;
rec.mode = m;
Svec = zeros(NS, 1);
enSlot = []; enW = [];
en = np.enabling{m};
li = find(en);
for t = 1:numel(li)
    [p, c] = ind2sub([sn.nnodes, R], li(t));
    slot = phOff(p, c) + 1;
    enSlot(end+1) = slot; %#ok<AGROW>
    enW(end+1) = en(li(t)); %#ok<AGROW>
    Svec(slot) = Svec(slot) - en(li(t));
end
fir = np.firing{m};
lf = find(fir);
for t = 1:numel(lf)
    [p, c] = ind2sub([sn.nnodes, R], lf(t));
    slot = phOff(p, c) + 1;
    Svec(slot) = Svec(slot) + fir(lf(t));
end
inhSlot = []; inhThr = [];
inh = np.inhibiting{m};
lh = find(~isinf(inh));
for t = 1:numel(lh)
    [p, c] = ind2sub([sn.nnodes, R], lh(t));
    inhSlot(end+1) = phOff(p, c) + 1; %#ok<AGROW>
    inhThr(end+1) = inh(lh(t)); %#ok<AGROW>
end
rec.Svec = Svec;
rec.enSlot = enSlot; rec.enW = enW;
rec.inhSlot = inhSlot; rec.inhThr = inhThr;
% Exponential firing rate: the single-phase completion rate sum(D1). A
% non-exponential firing distribution is rejected by the featset and must not
% reach here.
if np.timing(m) ~= TimingStrategy.IMMEDIATE
    fK = np.firingphases(m);
    if isnan(fK) || fK ~= 1 || isempty(np.firingproc{m})
        line_error(mfilename, sprintf('Transition %s mode %d has non-exponential firing, which the NRM SPN path does not support.', sn.nodenames{ind}, m));
    end
    D1 = np.firingproc{m}{2};
    rec.baseRate = sum(D1(:));
end
ns = np.nmodeservers(m);
if isinf(ns)
    ns = GlobalConstants.MaxInt();
end
rec.nservers = ns;
rec.weight = np.fireweight(m);
rec.prio = np.firingprio(m);
end

function d = spnEnDegree(nvec, rx)
% Enabling degree of a mode: the number of concurrent firings the marking
% supports, min over input arcs of floor(tokens/weight), zeroed by any active
% inhibitor arc. A mode with no input arc is treated as single-degree.
for i = 1:numel(rx.inhSlot)
    if nvec(rx.inhSlot(i)) >= rx.inhThr(i)
        d = 0;
        return
    end
end
if isempty(rx.enSlot)
    d = 1;
    return
end
d = inf;
for i = 1:numel(rx.enSlot)
    d = min(d, floor(nvec(rx.enSlot(i)) / rx.enW(i)));
end
end

function a = spnProp(nvec, rx)
% Propensity of a timed mode: the exponential rate times the effective number
% of servers, min(enabling degree, mode servers). Single-server modes therefore
% fire at their rate whenever enabled, infinite/k-server modes at the rate
% scaled by the enabling degree.
d = spnEnDegree(nvec, rx);
eff = min(d, rx.nservers);
if eff <= 0
    a = 0;
else
    a = rx.baseRate * eff;
end
end

function nvec = applyPlaceCaps(nvec, deposited, pcapSlot, placeTotalCaps)
% Drop tokens a firing pushed above a Place per-class or total capacity. Only the
% just-deposited slots (Svec > 0) can overflow, so the clamp is local. Mirrors the
% Python native _apply_place_caps.
for a = 1:numel(deposited)
    j = deposited(a);
    if nvec(j) > pcapSlot(j)
        nvec(j) = pcapSlot(j);
    end
end
for p = 1:size(placeTotalCaps, 1)
    tcap = placeTotalCaps{p, 1};
    slots = placeTotalCaps{p, 2};
    excess = sum(nvec(slots)) - tcap;
    if excess > 0
        for a = 1:numel(deposited)
            if excess <= 0
                break
            end
            j = deposited(a);
            if any(slots == j) && nvec(j) > 0
                d = min(excess, nvec(j));
                nvec(j) = nvec(j) - d;
                excess = excess - d;
            end
        end
    end
end
end

function nvec = spnCollapse(nvec, imm, maxsteps)
% Vanishing-marking elimination. Fire enabled immediate transitions until the
% marking is tangible: highest firing priority first, ties resolved in
% proportion to firing weight. Immediate firings take zero time and advance no
% clock, so the timed race only ever samples from tangible markings.
if isempty(imm)
    return
end
steps = 0;
while true
    enabled = [];
    for m = 1:numel(imm)
        if spnEnDegree(nvec, imm(m)) >= 1
            enabled(end+1) = m; %#ok<AGROW>
        end
    end
    if isempty(enabled)
        return
    end
    prios = zeros(1, numel(enabled));
    for i = 1:numel(enabled)
        prios(i) = imm(enabled(i)).prio;
    end
    top = enabled(prios == max(prios));  % larger firing priority = more urgent
    if numel(top) == 1
        pick = top;
    else
        w = zeros(1, numel(top));
        for i = 1:numel(top)
            w(i) = imm(top(i)).weight;
        end
        pick = top(spnWeightedDraw(w));
    end
    nvec = nvec + imm(pick).Svec;
    steps = steps + 1;
    if steps > maxsteps
        line_error(mfilename, 'Immediate-transition livelock: the vanishing-marking collapse did not reach a tangible marking.');
    end
end
end

function idx = spnWeightedDraw(w)
% Index drawn in proportion to the nonnegative weight vector W.
tot = sum(w);
if tot <= 0
    idx = 1;
    return
end
c = cumsum(w) / tot;
idx = find(c > rand, 1);
if isempty(idx)
    idx = numel(w);
end
end

