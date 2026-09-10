function [QN,UN,RN,TN,CN,XN,InfGen,StateSpace,StateSpaceAggr,EventFiltration,runtime,fname,sncopy,AuxFiltration,StartN,PreemptN] = solver_ctmc_analyzer(sn, options)
% [QN,UN,RN,TN,CN,XN,INFGEN,STATESPACE,STATESPACEAGGR,EVENTFILTRATION,RUNTIME,FNAME,sn,AUXFILTRATION,STARTN,PREEMPTN] = SOLVER_CTMC_ANALYZER(sn, OPTIONS)
%
% AUXFILTRATION carries the derived START/PREEMPT filtrations (see
% solver_ctmc); STARTN and PREEMPTN are the corresponding (station x class)
% rates pi*F*e: how often per unit time a class-r service starts at station i,
% and how often a class-r job in service is pushed back into the buffer there.
% For a lossless station with no in-service abandonment they satisfy
% STARTN == TN + PREEMPTN, which is the identity the tags exist to expose.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

%if options.remote
%    sn.rtfun = {};
%    sn.lst = {};
%    qn_json = jsonencode(sn);
%    sn = NetworkStruct.fromJSON(qn_json)
%return
%end

M = sn.nstations;    %number of stations
K = sn.nclasses;    %number of classes
S = sn.nservers;
NK = sn.njobs';  % initial population per class
sched = sn.sched;

Tstart = tic;
PH = sn.proc;

line_debug('CTMC analyzer starting: nstations=%d, nclasses=%d, njobs=%s', M, K, mat2str(NK));

% Note: hide_immediate now selectively preserves Cache immediate transitions
% in solver_ctmc.m, so we no longer need to disable it entirely for Cache nodes

line_debug('Building state space and infinitesimal generator via solver_ctmc');
% Resolve the enumeration default HERE, ahead of its first read: solver_ctmc
% fills it too, but the console line below already names it, so an options
% struct that never passed through SolverOptions('CTMC') failed on the field
% rather than on the analysis.
if ~isfield(options.config, 'state_space_gen')
    options.config.state_space_gen = 'default';
end
LineConsole.step('generating the state space (%s enumeration, cutoff %s)', ...
    options.config.state_space_gen, mat2str(options.cutoff));
[InfGen,StateSpace,StateSpaceAggr,EventFiltration,arvRates,depRates,sn,AuxFiltration] = solver_ctmc(sn, options); % sn is updated with the state space
LineConsole.step(['infinitesimal generator built: %d states, %d transitions, ' ...
    'density %.3g'], size(InfGen,1), nnz(InfGen)-size(InfGen,1), ...
    nnz(InfGen)/max(1,numel(InfGen)));

% if the initial state does not reflect the final size of the state
% vectors, attempt to correct it
for isf=1:sn.nstateful
    if size(sn.state{isf},2) < size(sn.space{isf},2)
        sn.state{isf} = [zeros(1,size(sn.space{isf},2)-size(sn.state{isf},2)),sn.state{isf}];
    end
end
sncopy = sn;

if options.keep
    line_debug('Saving CTMC data to file (options.keep=true)');
    fname = lineTempName;
    save([fname,'.mat'],'InfGen','StateSpace','StateSpaceAggr','EventFiltration')
    if LineConsole.isActive()
        LineConsole.substep('generator and state space saved in %s.mat', fname);
    else
        line_printf('CTMC infinitesimal generator and state space saved in: ');
        line_printf(strrep(sprintf('%s.mat\n',fname),'\','\\'))
    end
else
    fname = '';
end

wset = 1:length(InfGen);

line_debug('State space built: %d states, solving CTMC', length(InfGen));

% Note: solver_ctmc selectively preserves Cache immediate transitions to enable
% hit/miss rate computation while hiding other immediate transitions.
%
% every CTMC solve goes through the block decomposition; the irreducible case is the degenerate one BSCC / no transient states -- see ctmc_stationary
if issym(InfGen)
    % symbolic generators have no numeric SCC decomposition; ctmc_solve carries its own symbolic branch
    LineConsole.step('solving the symbolic generator');
    probSysState = ctmc_solve(InfGen, options);
else
    LineConsole.step('solving for the stationary distribution');
    probSysState = ctmc_stationary(InfGen, StateSpace, sn, options);
end
LineConsole.step('stationary distribution obtained, computing the mean metrics');
probSysState = probSysState(:)';
% clamp removes residues, but an ME stationary vector is genuinely SIGNED, so clamping deletes real mass -- see _kb/11-conventions-and-gotchas.md
if ~(isfield(sn,'isph') && ~isempty(sn.isph) && ~all(sn.isph(:)))
    probSysState(probSysState<GlobalConstants.Zero)=0;
end
probSysState = probSysState/sum(probSysState);

XN = NaN*zeros(1,K);
UN = NaN*zeros(M,K);
QN = NaN*zeros(M,K);
RN = NaN*zeros(M,K);
TN = NaN*zeros(M,K);
CN = NaN*zeros(1,K);

% Column span of each STATION inside a StateSpace row. A row is the
% concatenation of the per-node local states in STATEFUL index order
% (State.spaceGenerator: SS(ctr,:)=cell2mat(u) over isf=1..nstateful), and
% sn.space is keyed the same way. A stateful node need not be a station -- a
% Cache is stateful and is not -- so the running offset must walk every
% stateful node and be read back through sn.stationToStateful. Walking it with
% the station index instead both skipped the non-station widths and took the
% wrong node's width, handing State.toMarginal another node's columns.
sfSpaceShift = zeros(1,sn.nstateful);
for isf=2:sn.nstateful
    sfSpaceShift(isf) = sfSpaceShift(isf-1) + size(sn.space{isf-1},2);
end
istSpaceShift = zeros(1,M);
istSpaceWidth = zeros(1,M);
for ist=1:M
    isf = sn.stationToStateful(ist);
    istSpaceShift(ist) = sfSpaceShift(isf);
    istSpaceWidth(ist) = size(sn.space{isf},2);
end

for k=1:K
    refsf = sn.stationToStateful(sn.refstat(k));
    XN(k) = probSysState*arvRates(wset,refsf,k);
end

% see _kb/06-solver-catalog.md (G-network signals) for rationale
inDropRegion = false(1,M);
if isfield(sn,'nregions') && sn.nregions > 0
    for f=1:sn.nregions
        if sn.regionrule(f) == DropStrategy.DROP
            memb = sn_region_members(sn, f, sn.region{f}, []); % stations constrained by region f
            inDropRegion(1:min(M,numel(memb))) = inDropRegion(1:min(M,numel(memb))) | memb(1:min(M,numel(memb)));
        end
    end
end

for ist=1:M
    isf = sn.stationToStateful(ist);
    ind = sn.stationToNode(ist);
    isSource = sn.nodetype(ind) == NodeType.Source;
    for k=1:K
        TN(ist,k) = probSysState*depRates(wset,isf,k);
        if isSource
            % State.toMarginal encodes an EXT station as nir = Inf, an infinite
            % reservoir, which describes the state space and is not a queue
            % length. Reading it as one gave Q = Inf, and then R = Q/T = Inf.
            QN(ist,k) = 0;
        else
            QN(ist,k) = probSysState*StateSpaceAggr(wset,(ist-1)*K+k);
        end
    end
    if ~isSource
        % see _kb/06-solver-catalog.md (G-network signals) for rationale
        canDropClass = isinf(sn.njobs(:)') & (isfinite(sn.cap(ist)) | isfinite(sn.classcap(ist,:)) | inDropRegion(ist));
        signalLoss = ctmc_signal_lossy(sn, arvRates, probSysState, wset, isf);
        canDropClass = canDropClass | signalLoss;
        switch sched(ist)
            case SchedStrategy.INF
                for k=1:K
                    UN(ist,k) = QN(ist,k);
                end
            case {SchedStrategy.PS, SchedStrategy.DPS, SchedStrategy.GPS, SchedStrategy.LPS}
                if isempty(sn.lldscaling) && isempty(sn.cdscaling) && isempty(sn.jdscaling)
                    for k=1:K
                        if ~isempty(PH{ist}{k})
                            % see _kb/06-solver-catalog.md (Utilization conventions) for rationale
                            UNdep_ik = TN(ist,k)*map_mean(PH{ist}{k})/S(ist); % this is valid because CS in LINE is in a separate node
                            if canDropClass(k)
                                UN(ist,k) = UNdep_ik;
                            else
                                UNarv_ik = probSysState*arvRates(wset,isf,k)*map_mean(PH{ist}{k})/S(ist);
                                UN(ist,k) = max(UNarv_ik,UNdep_ik);
                            end
                        end
                    end
                else % lld/cd/ljd cases
                    % see _kb/06-solver-catalog.md (Utilization conventions) for rationale
                    ind = sn.stationToNode(ist);
                    ceff = S(ist);
                    if ~isempty(sn.lldscaling) && ist <= size(sn.lldscaling,1)
                        ceff = max(ceff, max(sn.lldscaling(ist,:)));
                    end
                    UN(ist,1:K) = 0;
                    for st = wset
                        [ni,nir] = State.toMarginal(sn, ind, StateSpace(st,(istSpaceShift(ist)+1):(istSpaceShift(ist)+istSpaceWidth(ist))));
                        if ni>0
                            lldnow = 1;
                            if ~isempty(sn.lldscaling) && ist <= size(sn.lldscaling,1)
                                lldnow = sn.lldscaling(ist, min(max(sum(ni),1), size(sn.lldscaling,2)));
                            end
                            for k=1:K
                                UN(ist,k) = UN(ist,k) + probSysState(st)*nir(k)*sn.schedparam(ist,k)/(nir*sn.schedparam(ist,:)')*lldnow/ceff;
                            end
                        end
                    end
                end
            case SchedStrategy.PAS
                % see _kb/06-solver-catalog.md (Utilization conventions) for rationale
                ind = sn.stationToNode(ist);
                UN(ist,1:K) = 0;
                for st = wset
                    [~,~,sir] = State.toMarginal(sn, ind, StateSpace(st,(istSpaceShift(ist)+1):(istSpaceShift(ist)+istSpaceWidth(ist))));
                    for k=1:K
                        UN(ist,k) = UN(ist,k) + probSysState(st)*sir(k)/S(ist);
                    end
                end
            otherwise
                if isempty(sn.lldscaling) && isempty(sn.cdscaling) && isempty(sn.jdscaling)
                    for k=1:K
                        if ~isempty(PH{ist}{k})
                            % see _kb/06-solver-catalog.md (Utilization conventions) for rationale
                            UNdep_ik = TN(ist,k)*map_mean(PH{ist}{k})/S(ist); % this is valid because CS in LINE is in a separate node
                            if canDropClass(k)
                                UN(ist,k) = UNdep_ik;
                            else
                                UNarv_ik = probSysState*arvRates(wset,isf,k)*map_mean(PH{ist}{k})/S(ist);
                                UN(ist,k) = max(UNarv_ik,UNdep_ik);
                            end
                        end
                    end
                else % lld/cd/ljd cases
                    % see _kb/06-solver-catalog.md (Utilization conventions) for rationale
                    ind = sn.stationToNode(ist);
                    ceff = S(ist);
                    if ~isempty(sn.lldscaling) && ist <= size(sn.lldscaling,1)
                        ceff = max(ceff, max(sn.lldscaling(ist,:)));
                    end
                    UN(ist,1:K) = 0;
                    for st = wset
                        [ni,~,sir] = State.toMarginal(sn, ind, StateSpace(st,(istSpaceShift(ist)+1):(istSpaceShift(ist)+istSpaceWidth(ist))));
                        if ni>0
                            lldnow = 1;
                            if ~isempty(sn.lldscaling) && ist <= size(sn.lldscaling,1)
                                lldnow = sn.lldscaling(ist, min(max(sum(ni),1), size(sn.lldscaling,2)));
                            end
                            sirtot = sum(sir);
                            for k=1:K
                                if sirtot > 0
                                    UN(ist,k) = UN(ist,k) + probSysState(st)*(sir(k)/sirtot)*lldnow/ceff;
                                end
                            end
                        end
                    end
                end
        end
        % see _kb/06-solver-catalog.md (G-network signals) for rationale
        if any(signalLoss) && sched(ist) ~= SchedStrategy.INF && ...
                isempty(sn.lldscaling) && isempty(sn.cdscaling) && isempty(sn.jdscaling)
            UNb = ctmc_signal_busy(sn, ind, ist, sched(ist), S(ist), StateSpace, istSpaceShift, wset, probSysState);
            UN(ist,signalLoss) = UNb(signalLoss);
        end
    end
end

% Synchronous calls (REPLY signals): a blocked caller still HOLDS its server, so add the held servers to QLen/Util; see _kb/06-solver-catalog.md
QNblocked = zeros(M,K);
if isfield(sn,'replyblock') && ~isempty(sn.replyblock)
    for ist=1:M
        ind = sn.stationToNode(ist);
        if size(sn.replyblock,1) < ind || ~any(sn.replyblock(ind,:) > 0)
            continue
        end
        rinfoU = State.replyBlockInfo(sn, ind);
        cols = (istSpaceShift(ist)+1):(istSpaceShift(ist)+istSpaceWidth(ist));
        bcols = cols(end-rinfoU.width+1:end);
        bmean = probSysState * StateSpace(wset, bcols); % 1 x width
        pos = 0;
        for r = rinfoU.classes
            pos = pos + 1;
            QN(ist,r) = QN(ist,r) + bmean(pos);
            UN(ist,r) = UN(ist,r) + bmean(pos)/S(ist);
            QNblocked(ist,r) = bmean(pos);
        end
    end
end

% Class-, joint- and global-dependence utilization normalization Util=T*S/peak,
% using the declared sn.cdscalingpeak, sn.jdscalingpeak and sn.gdscalingpeak;
% see _kb/06-solver-catalog.md (Utilization conventions) for rationale.
hasgd = isfield(sn,'gdscaling') && ~isempty(sn.gdscaling);
if (~isempty(sn.cdscaling) || ~isempty(sn.jdscaling) || hasgd) && ~any(isinf(sn.njobs))
    for ist=1:M
        hascd = ~isempty(sn.cdscaling) && length(sn.cdscaling) >= ist && ~isempty(sn.cdscaling{ist});
        hasjd = ~isempty(sn.jdscaling) && length(sn.jdscaling) >= ist && ~isempty(sn.jdscaling{ist});
        if ~hascd && ~hasjd && ~hasgd
            continue
        end
        for k=1:K
            % beta_r(n), eta_i(n) and phi(n) scale the SAME rate, so the peaks multiply
            bmax = 1;
            if hascd
                bmax = bmax * sn.cdscalingpeak(ist,k);
            end
            if hasjd
                bmax = bmax * sn.jdscalingpeak(ist,k);
            end
            if hasgd
                bmax = bmax * sn.gdscalingpeak(ist,k);
            end
            if isfinite(sn.rates(ist,k)) && sn.rates(ist,k) > 0 && bmax > 0
                UN(ist,k) = TN(ist,k) / sn.rates(ist,k) / bmax;
            else
                UN(ist,k) = 0;
            end
        end
    end
end


% see _kb/06-solver-catalog.md (True BAS blocking) for rationale
if ~isempty(sn.isbasblocking)
    for ist=1:M
        ind = sn.stationToNode(ist);
        % see _kb/06-solver-catalog.md (True BAS blocking) for rationale
        if numel(sn.isbasblocking) < ind || sn.isbasblocking(ind) ~= 1
            continue % no blocked marker at this station
        end
        % sn.isstation spans nnodes physical nodes plus one virtual entry per
        % finite capacity region; restrict to physical nodes to conform to connmatrix
        isstationNode = sn.isstation(1:sn.nnodes);
        dests = find(sn.connmatrix(ind,:) == 1 & isstationNode(:)' == 1);
        if numel(dests) ~= 1
            continue % ambiguous destination: leave the job where it sits
        end
        jst = sn.nodeToStation(dests(1));
        if isnan(jst) || jst < 1
            continue
        end
        cols = (istSpaceShift(ist)+1):(istSpaceShift(ist)+istSpaceWidth(ist));
        blocked = StateSpace(wset, cols(end)) == 1; % marker is the trailing column
        for k=1:K
            % Only the held job itself moves, not the whole queue at ist: a blocked
            % state holds exactly one completed job, so cap the per-state count at 1.
            shift = probSysState(blocked) * min(StateSpaceAggr(wset(blocked),(ist-1)*K+k), 1);
            if shift > 0
                QN(ist,k) = QN(ist,k) - shift;
                QN(jst,k) = QN(jst,k) + shift;
            end
        end
    end
end

for k=1:K
    for ist=1:M
        if TN(ist,k)>0
            % Response time is time spent AT the station, so the job blocked
            % out at the callee is excluded even though QLen/Util count it
            % (LDES measures the sojourn directly and reports the same).
            RN(ist,k) = (QN(ist,k)-QNblocked(ist,k))./TN(ist,k);
        else
            RN(ist,k)=0;
        end
    end
    CN(k) = NK(k)./XN(k);
end

QN(isnan(QN))=0;
CN(isnan(CN))=0;
RN(isnan(RN))=0;
UN(isnan(UN))=0;
XN(isnan(XN))=0;
TN(isnan(TN))=0;

% Derived service-start and preemption rates: pi*F*e over the same state set
% the metrics above use. They read the aux filtrations only, so no rate,
% probability or state above depends on them.
[StartN, PreemptN] = solver_ctmc_auxrates(AuxFiltration, probSysState, wset, M, K);

runtime = toc(Tstart);

% now update the routing probabilities in nodes with state-dependent routing
TNcache = zeros(sn.nstateful,K);
XNcache = zeros(sn.nstateful,K);
for k=1:K
    for isf=1:sn.nstateful
        ind = sncopy.statefulToNode(isf);
        if sncopy.nodetype(ind) == NodeType.Cache
            TNcache(isf,k) = probSysState*depRates(wset,isf,k);
            XNcache(isf,k) = probSysState*arvRates(wset,isf,k);
        end
    end
end

% Exact delayed-hit queue length of a retrieval-system cache. Block A of the cache
% local-variable vector marks the items being fetched and block B counts, per
% retrieval class, the secondary requests merged onto those fetches, so
%   phi_i    = P(a fetch of item i is in flight)
%   d1_i     = E[secondary requests waiting on the fetch of item i]
%   dfull_i  = d1_i + phi_i  (including the request that triggered the fetch)
% are state rewards of the stationary distribution, hence exact.
stateCols = cell(1, sncopy.nstateful);
colOff = 0;
for isf=1:sncopy.nstateful
    wsf = size(sncopy.space{isf},2);
    stateCols{isf} = (colOff+1):(colOff+wsf);
    colOff = colOff + wsf;
end
% Time-stationary per-item occupancy of each cache list. The cache-contents block
% of the local-variable vector holds the item index resident in each cache
% position, so P(item i is held by list l) is a state reward of the stationary
% distribution. This is the TIME-WEIGHTED occupancy, the CTMC counterpart of the
% EMBEDDED (per-request) occupancy the NC/MVA cache algorithms return; the two
% coincide only when requests see time averages (PASTA).
for isf=1:sncopy.nstateful
    ind = sncopy.statefulToNode(isf);
    if sncopy.nodetype(ind) ~= NodeType.Cache, continue; end
    np = sncopy.nodeparam{ind};
    itemcap = np.itemcap(:).';
    hlists = numel(itemcap);
    nitems = np.nitems;
    if hlists == 0 || nitems == 0, continue; end
    cols = stateCols{isf};
    if isfield(np,'retrievalSystemCapacity') && np.retrievalSystemCapacity > 0
        [~, rcItemsAll] = State.cacheRetrievalClassMap(sncopy, ind);
        lvw = sum(itemcap) + nitems + numel(rcItemsAll); % contents + block A + block B
    else
        lvw = sum(itemcap);
    end
    lvs = numel(cols) - lvw; % per-class server presence width
    if lvs < 0, continue; end
    itemprob = zeros(nitems, hlists+1);
    off = 0;
    for l=1:hlists
        lcols = cols(lvs+off+(1:itemcap(l)));
        off = off + itemcap(l);
        for i=1:nitems
            inlist = any(StateSpace(:,lcols) == i, 2);
            itemprob(i,l+1) = sum(probSysState(inlist));
        end
    end
    itemprob(:,1) = 1 - sum(itemprob(:,2:end),2);
    sncopy.nodeparam{ind}.actualitemprob = itemprob;
end

for isf=1:sncopy.nstateful
    ind = sncopy.statefulToNode(isf);
    if sncopy.nodetype(ind) ~= NodeType.Cache, continue; end
    np = sncopy.nodeparam{ind};
    if ~isfield(np,'retrievalSystemCapacity') || np.retrievalSystemCapacity <= 0, continue; end
    [~, rcItems, rcOrigClass] = State.cacheRetrievalClassMap(sncopy, ind);
    nitems = np.nitems;
    tcc = np.totalCacheCapacity;
    lvw = size(sncopy.space{isf},2);
    cols = stateCols{isf};
    lvs = lvw - (tcc + nitems + numel(rcItems)); % per-class server presence width
    aCols = cols(lvs+tcc+(1:nitems));
    bCols = cols(lvs+tcc+nitems+(1:numel(rcItems)));
    phi = zeros(1,nitems); d1 = zeros(1,nitems);
    for i=1:nitems
        phi(i) = sum(probSysState(StateSpace(:,aCols(i)) ~= 0));
        bsel = bCols(rcItems == i);
        if ~isempty(bsel)
            d1(i) = sum(probSysState(:) .* sum(StateSpace(:,bsel),2));
        end
    end
    sncopy.nodeparam{ind}.delayedhitprobitem = phi;
    sncopy.nodeparam{ind}.delayedhitqlen = d1;
    sncopy.nodeparam{ind}.delayedhitqlenfull = d1 + phi;

    % Exact delayed-hit rate per originating class. A fetch of item i completes on
    % exactly the transitions that clear block A bit i, and each such transition
    % releases the block-B counts of item i as delayed hits. The rate is therefore
    % a TRANSITION reward over the generator, not a state reward: the alternative
    % arrival-rate identity lambda_i*phi_i is only PASTA-exact.
    delayedRate = zeros(1, sncopy.nclasses);
    if ~isempty(rcItems)
        offdiag = InfGen - diag(diag(InfGen));
        for j = 1:numel(rcItems)
            i = rcItems(j);
            rows = find(StateSpace(:,aCols(i)) ~= 0 & StateSpace(:,bCols(j)) > 0);
            for rr = rows(:).'
                nz = find(offdiag(rr,:) ~= 0);
                completes = nz(StateSpace(nz, aCols(i)) == 0);
                if isempty(completes), continue; end
                delayedRate(rcOrigClass(j)) = delayedRate(rcOrigClass(j)) ...
                    + probSysState(rr) * StateSpace(rr,bCols(j)) * sum(offdiag(rr,completes));
            end
        end
    end
    sncopy.nodeparam{ind}.delayedhitrate = delayedRate;
end

% updates cache actual hit and miss data + retrieval-system expected latency.
% Hit / miss for class k are derived from departure rates of the configured
% hitClass / missClass at the cache; for retrieval-aware caches the
% miss-class departure rate equals the true miss rate because
% afterEventCache fires the miss event on a retrieval-complete READ.
retrievalLatencyWarned = false;
for k=1:K
    for isf=1:sncopy.nstateful
        ind = sncopy.statefulToNode(isf);
        if sncopy.nodetype(ind) == NodeType.Cache
            np = sncopy.nodeparam{ind};
            if length(np.hitclass)>=k
                h = np.hitclass(k);
                m = np.missclass(k);
                if h>0 && m>0
                    % Delayed hits depart in the hit class, so the hit-class rate is
                    % (true hits + delayed hits); the exact delayed rate computed above
                    % splits it. hit + delayed + miss = 1, matching the LDES/NC report.
                    denom = sum(TNcache(isf,[h,m]));
                    dRate = 0;
                    if isfield(np,'delayedhitrate') && numel(np.delayedhitrate) >= k
                        dRate = min(np.delayedhitrate(k), TNcache(isf,h));
                    end
                    sncopy.nodeparam{ind}.actualhitprob(k) = (TNcache(isf,h)-dRate)/denom;
                    sncopy.nodeparam{ind}.actualdelayedhitprob(k) = dRate/denom;
                    sncopy.nodeparam{ind}.actualmissprob(k) = TNcache(isf,m)/denom;
                else
                    sncopy.nodeparam{ind}.actualhitprob(k) = NaN;
                    sncopy.nodeparam{ind}.actualdelayedhitprob(k) = NaN;
                    sncopy.nodeparam{ind}.actualmissprob(k) = NaN;
                end

                % The Eq. 8 retrieval-system expected latency is not currently
                % implemented; report NaN whenever a retrieval system is
                % configured for this class.
                expectedLatency = NaN;
                if isfield(np, 'retrievalSystemQueueIndices') ...
                        && isKey(np.retrievalSystemQueueIndices, int32(k-1)) ...
                        && ~isempty(np.retrievalSystemQueueIndices{int32(k-1)})
                    if ~retrievalLatencyWarned
                        line_warning(mfilename, 'Retrieval-system expected latency is not currently implemented; reporting NaN.');
                        retrievalLatencyWarned = true;
                    end
                end
                sncopy.nodeparam{ind}.actualresidt(k) = expectedLatency;
            end
        end
    end
end
end
