function [QN,UN,RN,TN,CN,XN,InfGen,StateSpace,StateSpaceAggr,EventFiltration,runtime,fname,sncopy] = solver_ctmc_analyzer(sn, options)
% [QN,UN,RN,TN,CN,XN,INFGEN,STATESPACE,STATESPACEAGGR,EVENTFILTRATION,RUNTIME,FNAME,sn] = SOLVER_CTMC_ANALYZER(sn, OPTIONS)
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
[InfGen,StateSpace,StateSpaceAggr,EventFiltration,arvRates,depRates,sn] = solver_ctmc(sn, options); % sn is updated with the state space

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
    line_printf('CTMC infinitesimal generator and state space saved in: ');
    line_printf(strrep(sprintf('%s.mat\n',fname),'\','\\'))
else
    fname = '';
end

wset = 1:length(InfGen);

line_debug('State space built: %d states, solving CTMC', length(InfGen));

use_ctmc_solve_stable = true;
if use_ctmc_solve_stable
    % stable version
    % Note: solver_ctmc now selectively preserves Cache immediate transitions
    % to enable hit/miss rate computation while hiding other immediate transitions
    [probSysState, ~, nConnComp, connComp] = ctmc_solve(InfGen, options);

    if nConnComp > 1
        line_debug('CTMC is reducible: %d connected components', nConnComp);
        % the matrix was reducible
        initState = matchrow(StateSpace, cell2mat(sn.state'));
        if initState <= 0
            % Initial state may have been removed by stochcomp (e.g., SPN with
            % immediate ENABLE states). Use the largest connected component.
            compSizes = accumarray(connComp(:), 1);
            [~, largestComp] = max(compSizes);
            wset = find(connComp == largestComp);
        else
            % determine the weakly connected component associated to the initial state
            wset = find(connComp == connComp(initState));
        end
        if initState > 0
            line_debug('Using component %d with %d states (from initial state)', connComp(initState), length(wset));
        else
            line_debug('Using largest component with %d states (initial state removed by stochcomp)', length(wset));
        end
        probSysState = ctmc_solve(InfGen(wset, wset), options);
        InfGen = InfGen(wset, wset);
        % reduce all per-state arrays to the retained component and remap wset
        % to local indices so matrix-form and loop-form estimators stay aligned
        StateSpace = StateSpace(wset,:);
        StateSpaceAggr = StateSpaceAggr(wset,:);
        arvRates = arvRates(wset,:,:);
        depRates = depRates(wset,:,:);
        wset = 1:numel(wset);
    else
        line_debug('CTMC is irreducible, using full state space');
    end
else
    % development version

    % we now find the initial state and then solver the CTMC allowing for the
    % case where it is reducible
    initState = matchrow(StateSpace, cell2mat(sn.state'));
    pi0 = zeros(1,length(InfGen)); pi0(initState) = 1.0;
    [pi,pis,~,scc,~] = ctmc_solve_reducible(InfGen, pi0, options);

    if size(pis,1)==1
        probSysState = pi;
    else
        wset = scc == scc(initState);
        InfGen = InfGen(wset, wset);
        StateSpace = StateSpace(wset,:);
        probSysState = pis(scc(initState),scc == scc(initState));
    end
end
% The clamp removes the tiny negative residues a genuine CTMC solve leaves
% behind. With a matrix-exponential process the stationary vector is a genuinely
% SIGNED measure -- only its aggregates over each phase block are probabilities
% -- so clamping there deletes real mass and every mean measure moves. Mean
% measures are linear in the vector and stay exact without the clamp.
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

istSpaceShift = zeros(1,M);
for ist=1:M
    if ist==1
        istSpaceShift(ist) = 0;
    else
        istSpaceShift(ist) = istSpaceShift(ist-1) + size(sn.space{ist-1},2);
    end
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
            memb = any(sn.region{f} ~= -1, 2); % stations constrained by region f
            memb = memb(:)';
            inDropRegion(1:min(M,numel(memb))) = inDropRegion(1:min(M,numel(memb))) | memb(1:min(M,numel(memb)));
        end
    end
end

for ist=1:M
    isf = sn.stationToStateful(ist);
    ind = sn.stationToNode(ist);
    for k=1:K
        TN(ist,k) = probSysState*depRates(wset,isf,k);
        QN(ist,k) = probSysState*StateSpaceAggr(wset,(ist-1)*K+k);
    end
    if sn.nodetype(ind) ~= NodeType.Source
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
                        [ni,nir] = State.toMarginal(sn, ind, StateSpace(st,(istSpaceShift(ist)+1):(istSpaceShift(ist)+size(sn.space{ist},2))));
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
                    [~,~,sir] = State.toMarginal(sn, ind, StateSpace(st,(istSpaceShift(ist)+1):(istSpaceShift(ist)+size(sn.space{ist},2))));
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
                        [ni,~,sir] = State.toMarginal(sn, ind, StateSpace(st,(istSpaceShift(ist)+1):(istSpaceShift(ist)+size(sn.space{ist},2))));
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
        cols = (istSpaceShift(ist)+1):(istSpaceShift(ist)+size(sn.space{ist},2));
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

% see _kb/06-solver-catalog.md (Utilization conventions) for rationale
if ~isempty(sn.cdscaling) && ~any(isinf(sn.njobs))
    for ist=1:M
        if length(sn.cdscaling) >= ist && ~isempty(sn.cdscaling{ist})
            for k=1:K
                bmax = sn.cdscalingpeak(ist,k);
                if isfinite(sn.rates(ist,k)) && sn.rates(ist,k) > 0 && bmax > 0
                    UN(ist,k) = TN(ist,k) / sn.rates(ist,k) / bmax;
                else
                    UN(ist,k) = 0;
                end
            end
        end
    end
end

% Joint-dependence (non-product-form) utilization normalization: Util=T*S/peak
% using the declared sn.jdscalingpeak, mirroring the class-dependence block.
if ~isempty(sn.jdscaling) && ~any(isinf(sn.njobs))
    for ist=1:M
        if length(sn.jdscaling) >= ist && ~isempty(sn.jdscaling{ist})
            for k=1:K
                bmax = sn.jdscalingpeak(ist,k);
                if isfinite(sn.rates(ist,k)) && sn.rates(ist,k) > 0 && bmax > 0
                    UN(ist,k) = TN(ist,k) / sn.rates(ist,k) / bmax;
                else
                    UN(ist,k) = 0;
                end
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
        cols = (istSpaceShift(ist)+1):(istSpaceShift(ist)+size(sn.space{ist},2));
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
                    sncopy.nodeparam{ind}.actualhitprob(k) = TNcache(isf,h)/sum(TNcache(isf,[h,m]));
                    sncopy.nodeparam{ind}.actualmissprob(k) = TNcache(isf,m)/sum(TNcache(isf,[h,m]));
                else
                    sncopy.nodeparam{ind}.actualhitprob(k) = NaN;
                    sncopy.nodeparam{ind}.actualmissprob(k) = NaN;
                end

                % The Eq. 8 retrieval-system expected latency is not currently
                % implemented; report NaN whenever a retrieval system is
                % configured for this class.
                expectedLatency = NaN;
                if isfield(np, 'retrievalSystemQueueIndices') ...
                        && isKey(np.retrievalSystemQueueIndices, int32(k-1)) ...
                        && ~isempty(np.retrievalSystemQueueIndices(int32(k-1)))
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
