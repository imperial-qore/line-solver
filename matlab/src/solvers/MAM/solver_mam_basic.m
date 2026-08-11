function [QN,UN,RN,TN,CN,XN,totiter] = solver_mam_basic(sn, options)
% [Q,U,R,T,C,X] = SOLVER_MAM_BASIC(QN, PH, OPTIONS)

% This solver uses MAM to solve queues in isolation, but simplifies the
% traffic equations by using visits to rescale the flows into the queue
% inputs.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

config = options.config;
tol = options.tol;

PH = sn.proc;
%% generate local state spaces
I = sn.nnodes;
M = sn.nstations;
K = sn.nclasses;
C = sn.nchains;
N = sn.njobs';
V = cellsum(sn.visits);
S = 1./sn.rates;
Strue = S; % service times as declared, used to report utilization below
Lchain = sn_get_demands_chain(sn);

% SLC interference via inflated service time; see _kb/06-solver-catalog.md for rationale
slcjobs = zeros(M,1);
for k=1:K
    if sn.isslc(k)
        ist_k = sn.refstat(k);
        if isfinite(sn.nservers(ist_k))
            slcjobs(ist_k) = slcjobs(ist_k) + sn.njobs(k);
        end
    end
end
for ist=1:M
    if slcjobs(ist) > 0
        for k=1:K
            if ~sn.isslc(k)
                S(ist,k) = S(ist,k)*(1+slcjobs(ist));
            end
        end
        % The chain demands drive the throughput fixed point and must be
        % consistent with the inflated service times.
        Lchain(ist,:) = Lchain(ist,:)*(1+slcjobs(ist));
    end
end

QN = zeros(M,K);
UN = zeros(M,K);
RN = zeros(M,K);
TN = zeros(M,K);
CN = zeros(1,K);
XN = zeros(1,K);

% Track stations using exact MAP/D/c solver (skip post-processing for these)
mapdcStations = false(M, 1);

pie = {};
D0 = {};

lambda = zeros(1,C);
chainSysArrivals = cell(1,C);
TN_1 = TN+Inf;

it = 0;

% open queueing system (one node is the external world)
% first build the joint arrival process
for ist=1:M
    switch sn.sched(ist)
        case {SchedStrategy.FCFS, SchedStrategy.HOL, SchedStrategy.FCFSPRPRIO, SchedStrategy.PS}
            for k=1:K
                % Skip Det processes - they will be handled by qsys_mapdc
                if sn.procid(ist, k) == ProcessType.DET
                    % For Det, we don't need MAP representation since we use exact solver
                    pie{ist}{k} = [];
                    D0{ist,k} = [];
                    continue;
                end
                % divide service time by nservers, surrogate delay in tandem;
                % ME/RAP rescaled by rate alone, see _kb/06-solver-catalog.md for rationale
                if sn.procid(ist, k) == ProcessType.ME || sn.procid(ist, k) == ProcessType.RAP
                    ratio = map_mean(PH{ist}{k}) / (S(ist,k)/sn.nservers(ist));
                    PH{ist}{k} = {PH{ist}{k}{1}*ratio, PH{ist}{k}{2}*ratio};
                else
                    PH{ist}{k} = map_scale(PH{ist}{k}, S(ist,k)/sn.nservers(ist));
                end
                pie{ist}{k} = map_pie(PH{ist}{k});
                D0{ist,k} = PH{ist}{k}{1};
                if any(isnan(D0{ist,k}))
                    D0{ist,k} = -GlobalConstants.Immediate;
                    pie{ist}{k} = 1;
                    PH{ist}{k} = map_exponential(GlobalConstants.Immediate);
                end
            end
    end
end % i

isOpen = false;
isClosed = false;
if any(isinf(sn.njobs))
    isOpen = true;
end
if any(isfinite(sn.njobs))
    isClosed = true;
end
isMixed = isOpen & isClosed;

% SLC chains excluded from the throughput fixed point; see _kb/06-solver-catalog.md for rationale
isslcchain = false(1,C);
for c=1:C
    isslcchain(c) = all(sn.isslc(sn.inchain{c}));
end

lambdas_inchain = cell(1,C);
for c=1:C
    inchain = sn.inchain{c};
    lambdas_inchain{c} = sn.rates(sn.refstat(inchain(1)),inchain);
    %lambdas_inchain{c} = lambdas_inchain{c}(isfinite(lambdas_inchain{c}));
    lambda(c) = sum(lambdas_inchain{c}(isfinite(lambdas_inchain{c})));
    if isslcchain(c)
        lambda(c) = 0;
        lambdas_inchain{c} = zeros(1,length(inchain));
    end
    ist = sn.refstat(inchain(1)); % identical for all classes in the chain
    if isinf(sum(sn.njobs(inchain))) % if open chain
        % ist here is the source
        % assemble a MMAP for the arrival process from all classes
        for k=1:K
            if isnan(PH{ist}{k}{1})
                PH{ist}{k} = map_exponential(Inf); % no arrivals from this class
            end
        end
        inchain = sn.inchain{c};
        k = inchain(1);
        chainSysArrivals{c} = {PH{ist}{k}{1},PH{ist}{k}{2},PH{ist}{k}{2}};
        for ki=2:length(inchain)
            k = inchain(ki);
            if isnan(PH{ist}{k}{1})
                PH{ist}{k} = map_exponential(Inf); % no arrivals from this class
            end
            chainSysArrivals{c} = mmap_super_safe({chainSysArrivals{c},{PH{ist}{k}{1},PH{ist}{k}{2},PH{ist}{k}{2}}}, config.space_max, 'default');
        end
        TN(ist,inchain') = lambdas_inchain{c};
    end
end

sd = isfinite(sn.nservers);

isclosedchain = false(1,C);
isopenchain = false(1,C);
for c=1:C
    isopenchain(c) = isinf(sum(sn.njobs(sn.inchain{c})));
    isclosedchain(c) = ~isopenchain(c) && ~isslcchain(c);
end
% Mixed-network backoff differs from uniform 1/Umax; see _kb/06-solver-catalog.md for rationale
ismixed = any(isclosedchain) && any(isopenchain);
% Closed chains claim only the capacity open traffic leaves free; stay inside
% the stability region (see _kb/06-solver-catalog.md for rationale)
Ulim = 1 - GlobalConstants.CoarseTol;

% Stations already flagged as carrying a non-phase-type service process, scoped
% to this solve so the fixed-point iteration reports each station once while
% every new model is still warned about.
meWarned = false(1, M);

while max(max(abs(TN-TN_1))) > tol && it <= options.iter_max %#ok<max>
    it = it + 1;
    TN_1 = TN;
    Umax = max(sum(UN(sd,:),2));
    if ismixed || Umax < 1
        for c=1:C
            inchain = sn.inchain{c};
            if isclosedchain(c)
                Nc = sum(sn.njobs(inchain)); % closed population
                QNc = max(tol, sum(sum(QN(:,inchain),2,"omitnan"))); %#ok<NANSUM>
                TNlb = Nc./sum(Lchain(:,c));
                if it == 1
                    lambda(c) = TNlb; % lower bound
                else
                    lambda(c) = lambda(c) * it/options.iter_max + (Nc / QNc)  * lambda(c) * (options.iter_max-it)/options.iter_max; % iteration-averaged regula falsi;
                end
            end
        end
    end
    if ismixed
        % Back closed chains onto the busiest station's residual capacity;
        % see _kb/06-solver-catalog.md for rationale
        Uchain  = Lchain(sd,:) .* lambda;
        Uopen   = sum(Uchain(:,~isclosedchain), 2);
        Uclosed = sum(Uchain(:, isclosedchain), 2);
        binding = Uclosed > tol;
        if any(binding)
            theta = min((Ulim - Uopen(binding)) ./ Uclosed(binding));
            if theta < 1
                lambda(isclosedchain) = lambda(isclosedchain) * max(0, theta);
            end
        end
    elseif Umax >= 1
        lambda = lambda * 1/Umax;
    end

    for c=1:C
        inchain = sn.inchain{c};
        if isslcchain(c)
            % SLC vanishing-rate arrival surrogate; see _kb/06-solver-catalog.md for rationale
            chainSysArrivals{c} = mmap_exponential(GlobalConstants.Zero);
        elseif ~isinf(sum(sn.njobs(inchain)))
            % Closed chain: Poisson surrogate at current throughput iterate;
            % see _kb/06-solver-catalog.md for rationale
            chainSysArrivals{c} = mmap_exponential(lambda(c));
        end
        for m=1:M
            TN(m,inchain) = V(m,inchain) .* lambda(c);
        end
    end

    for ind=1:I
        if sn.isstation(ind)
            ist = sn.nodeToStation(ind);
            switch sn.nodetype(ind)
                case NodeType.Join
                    for c=1:C
                        inchain = sn.inchain{c};
                        for k=inchain
                            fanin = nnz(sn.rtnodes(:, (ind-1)*K+k));
                            TN(ist,k) = lambda(c)*V(ist,k)/fanin;
                            UN(ist,k) = 0;
                            QN(ist,k) = 0;
                            RN(ist,k) = 0;
                        end
                    end
                otherwise
                    switch sn.sched(ist)
                        case SchedStrategy.INF
                            for c=1:C
                                inchain = sn.inchain{c};
                                for k=inchain
                                    if V(ist,k) == 0
                                        % Non-visiting class NaN guard; see _kb/06-solver-catalog.md for rationale
                                        TN(ist,k) = 0;
                                        UN(ist,k) = 0;
                                        QN(ist,k) = 0;
                                        RN(ist,k) = 0;
                                        continue;
                                    end
                                    TN(ist,k) = lambda(c)*V(ist,k);
                                    % INF station: U = QLen = TN*S (no /c), single-V Little;
                                    % see _kb/06-solver-catalog.md for rationale
                                    UN(ist,k) = S(ist,k)*TN(ist,k);
                                    QN(ist,k) = TN(ist,k).*S(ist,k);
                                    RN(ist,k) = QN(ist,k)/TN(ist,k);
                                end
                            end
                        case SchedStrategy.PS
                            for c=1:C
                                inchain = sn.inchain{c};
                                for k=inchain
                                    if V(ist,k) == 0
                                        % Non-visiting class NaN guard; see _kb/06-solver-catalog.md for rationale
                                        TN(ist,k) = 0;
                                        UN(ist,k) = 0;
                                        continue;
                                    end
                                    TN(ist,k) = lambda(c)*V(ist,k);
                                    % Utilization Law: a c-server station holds TN*S/c of its capacity.
                                    UN(ist,k) = S(ist,k)*TN(ist,k)/sn.nservers(ist);
                                end
                                %Nc = sum(sn.njobs(inchain)); % closed population
                                Uden = min([1-GlobalConstants.FineTol,sum(UN(ist,:))]);
                                for k=inchain
                                    if V(ist,k) == 0
                                        QN(ist,k) = 0;
                                        RN(ist,k) = 0;
                                        continue;
                                    end
                                    %QN(ist,k) = (UN(ist,k)-UN(ist,k)^(Nc+1))/(1-Uden); % geometric bound type approximation
                                    QN(ist,k) = UN(ist,k)/(1-Uden);
                                    RN(ist,k) = QN(ist,k)/TN(ist,k);
                                end
                            end
                        case {SchedStrategy.FCFS, SchedStrategy.HOL, SchedStrategy.FCFSPRPRIO}
                            chainArrivalAtNode = cell(1,C);
                            rates = cell(M,C);
                            mapdcUsed = false;  % Flag for exact MAP/D/c solver
                            finiteCapUsed = false;  % Flag for finite-capacity truncate-and-renormalize
                            finiteCapLossPerClass = [];  % per-class loss, exact MMAP/G/1/K branch only
                            for c=1:C %for each chain
                                rates{ist,c} = V(ist,:) .* lambda(c); % visits of classes within the chain
                                inchain = find(sn.chains(c,:))';
                                markProb = rates{ist,c}(inchain) / sum(rates{ist,c}(inchain));
                                markProb(isnan(markProb)) = 0;
                                chainArrivalAtNode{c} = mmap_mark(chainSysArrivals{c}, markProb);
                                % Normalize only Markovian arrivals, not ME/RAP;
                                % see _kb/06-solver-catalog.md for rationale
                                if mam_chain_arrival_is_markovian(sn, c)
                                    chainArrivalAtNode{c} = mmap_normalize(chainArrivalAtNode{c});
                                end
                                % mmap_scale targets marks 1:C; pass chain's own visit
                                % rates, skip zero-rate SLC; see _kb/06-solver-catalog.md for rationale
                                tgtrates = rates{ist,c}(inchain);
                                if any(tgtrates > 0)
                                    chainArrivalAtNode{c} = mmap_scale(chainArrivalAtNode{c}, 1./tgtrates, 0); % non-iterative approximation
                                end
                                chainIsMarkovian = mam_chain_arrival_is_markovian(sn, c);
                                if c == 1
                                    if chainIsMarkovian
                                        aggrArrivalAtNode = mmap_super_safe({chainArrivalAtNode{c}, mmap_exponential(0,1)}, config.space_max, 'default');
                                        aggrArrivalAtNode = {aggrArrivalAtNode{1} aggrArrivalAtNode{2} aggrArrivalAtNode{2}};
                                    else
                                        % ME/RAP: build {D0,D1,D1} shape directly (avoid
                                        % mmap_super normalize); see _kb/06-solver-catalog.md for rationale
                                        aggrArrivalAtNode = {chainArrivalAtNode{c}{1}, chainArrivalAtNode{c}{2}, chainArrivalAtNode{c}{2}};
                                    end
                                    lc = map_lambda(chainArrivalAtNode{c});
                                    if lc>0
                                        aggrArrivalAtNode = mmap_scale(aggrArrivalAtNode, 1/lc, 0); % non-iterative approximation
                                    end
                                else
                                    if ~chainIsMarkovian
                                        % ME/RAP correlation only approximated under superposition;
                                        % see _kb/06-solver-catalog.md for rationale
                                        line_warning_always(mfilename, ...
                                            'Chain %d has a matrix-exponential or rational arrival process, which the superposition of several chains at station %s cannot represent exactly. Its autocorrelation is approximated by the closest Markovian arrival process.', ...
                                            c, sn.nodenames{sn.stationToNode(ist)});
                                    end
                                    aggrArrivalAtNode = mmap_super_safe({aggrArrivalAtNode, chainArrivalAtNode{c}}, config.space_max, 'default');
                                end
                            end
                            Qret = cell(1,K);
                            if (sn.sched(ist)==SchedStrategy.HOL && any(sn.classprio ~= sn.classprio(1))) % if priorities are not identical; sched holds numeric ids so == must be used (strcmp on numerics is always false)
                                [uK,iK] = unique(sn.classprio);
                                % BUTools convention: D1=lowest priority, DK=highest priority
                                % LINE convention: lower value = higher priority
                                % unique() returns ascending order, so we need to reverse for BUTools
                                iK = flipud(iK(:));
                                if length(uK) == length(sn.classprio) % if all priorities are different
                                    [Qret{iK'}] = MMAPPH1NPPR({aggrArrivalAtNode{[1;2+iK]}}, {pie{ist}{iK}}, {D0{ist,iK}}, 'ncMoms', 1);
                                else
                                    line_error(mfilename,'Solver MAM requires either identical priorities or all distinct priorities');
                                end
                            elseif (sn.sched(ist)==SchedStrategy.FCFSPRPRIO && any(sn.classprio ~= sn.classprio(1))) % FCFS preemptive resume priority
                                [uK,iK] = unique(sn.classprio);
                                % BUTools convention: D1=lowest priority, DK=highest priority
                                % LINE convention: lower value = higher priority
                                % unique() returns ascending order, so we need to reverse for BUTools
                                iK = flipud(iK(:));
                                if length(uK) == length(sn.classprio) % if all priorities are different
                                    [Qret{iK'}] = MMAPPH1PRPR({aggrArrivalAtNode{[1;2+iK]}}, {pie{ist}{iK}}, {D0{ist,iK}}, 'ncMoms', 1);
                                else
                                    line_error(mfilename,'Solver MAM requires either identical priorities or all distinct priorities');
                                end
                            else
                                aggrUtil = sum(mmap_lambda(aggrArrivalAtNode)./(GlobalConstants.FineTol+sn.rates(ist,1:K)*sn.nservers(ist)), 'omitnan');   %% to debug
                                aggrLambda = mmap_lambda(aggrArrivalAtNode);
                                if aggrUtil < 1-GlobalConstants.FineTol
                                    if any(isinf(N))
                                        % Check for MAP/D/c: single-class open model with deterministic service
                                        isMapDc = (K == 1) && (sn.procid(ist, 1) == ProcessType.DET);
                                        % Check for D/M/c: single-class open, exp service at this queue,
                                        % deterministic source elsewhere.
                                        isDMc = false;
                                        dmcSourceIdx = -1;
                                        if K == 1 && ~isMapDc && sn.procid(ist, 1) == ProcessType.EXP
                                            for jst = 1:sn.nstations
                                                if jst ~= ist && sn.procid(jst, 1) == ProcessType.DET
                                                    isDMc = true;
                                                    dmcSourceIdx = jst;
                                                    break;
                                                end
                                            end
                                        end
                                        % Check for PH/M/c: single-class open, exp service at this queue,
                                        % source has PH inter-arrivals (Erlang/Coxian/MAP/PH/...).
                                        % Covers both c=1 (PH/M/1) and c>1 (PH/M/c via matrix-geometric).
                                        isPhM1 = false;
                                        phM1SourceIdx = -1;
                                        if K == 1 && ~isMapDc && ~isDMc && sn.procid(ist, 1) == ProcessType.EXP ...
                                                && isfinite(sn.nservers(ist)) && sn.nservers(ist) >= 1
                                            for jst = 1:sn.nstations
                                                if jst == ist
                                                    continue;
                                                end
                                                srcProc = sn.procid(jst, 1);
                                                if srcProc ~= ProcessType.EXP && srcProc ~= ProcessType.DET ...
                                                        && srcProc ~= ProcessType.IMMEDIATE && srcProc ~= ProcessType.DISABLED
                                                    % PH/M/c gate is RENEWAL, not non-exponential;
                                                    % see _kb/06-solver-catalog.md for rationale
                                                    if ~mam_srcproc_is_renewal(sn, jst)
                                                        continue;
                                                    end
                                                    isPhM1 = true;
                                                    phM1SourceIdx = jst;
                                                    break;
                                                elseif srcProc == ProcessType.EXP && sn.nservers(ist) > 1 ...
                                                        && sn.nstations == 2 ...
                                                        && sn.nodetype(sn.stationToNode(jst)) == NodeType.Source
                                                    % M/M/c: the single-fast-server surrogate used by the
                                                    % generic path is inexact for c>1; dispatch to the exact
                                                    % PH/M/c matrix-geometric solver (Erlang-C result)
                                                    isPhM1 = true;
                                                    phM1SourceIdx = jst;
                                                    break;
                                                end
                                            end
                                        end
                                        % Finite buffer takes precedence over infinite-buffer
                                        % closed forms; see _kb/06-solver-catalog.md for rationale
                                        isFiniteCap = isfinite(sn.cap(ist));
                                        if isFiniteCap
                                            isPhM1 = false;
                                            isDMc = false;
                                            isMapDc = false;
                                        end
                                        if isPhM1
                                            muQ = 1.0 / S(ist, 1);
                                            try
                                                phPair = sn.proc{phM1SourceIdx}{1};
                                                D0_ph = phPair{1};
                                                D1_ph = phPair{2};
                                                pie_src = map_pie({D0_ph, D1_ph});
                                                phm1Result = qsys_phmc(pie_src, D0_ph, muQ, sn.nservers(ist));
                                                Qret{1} = phm1Result.meanQueueLength;
                                                mapdcUsed = true;
                                                mapdcStations(ist) = true;
                                                line_debug(options, 'Using exact PH/M/%d solver: Q=%.4f, W=%.4f', ...
                                                    sn.nservers(ist), phm1Result.meanQueueLength, phm1Result.meanWaitingTime);
                                            catch
                                                isPhM1 = false;
                                            end
                                        end
                                        if isDMc
                                            muQ = 1.0 / S(ist, 1);
                                            lamD = sn.rates(dmcSourceIdx, 1);
                                            try
                                                dmcResult = qsys_dmc(lamD, muQ, sn.nservers(ist));
                                                Qret{1} = dmcResult.meanQueueLength;
                                                mapdcUsed = true;  % skip surrogate delay correction
                                                mapdcStations(ist) = true;
                                                line_debug(options, 'Using exact D/M/%d solver: Q=%.4f, W=%.4f', ...
                                                    sn.nservers(ist), dmcResult.meanQueueLength, dmcResult.meanWaitingTime);
                                            catch
                                                isDMc = false;
                                            end
                                        end
                                        if isPhM1 || isDMc
                                            % handled above; skip the MAP/D/c, finite-cap, and MMAPPH1FCFS branches
                                        elseif isMapDc
                                            % Use exact MAP/D/c solver from Q-MAM
                                            D0_arr = aggrArrivalAtNode{1};
                                            D1_arr = aggrArrivalAtNode{2};
                                            detServiceTime = S(ist, 1);  % 1/rate = service time
                                            numServers = sn.nservers(ist);

                                            mapdcResult = qsys_mapdc(D0_arr, D1_arr, detServiceTime, numServers);
                                            Qret{1} = mapdcResult.meanQueueLength;
                                            % Store result for later use to skip surrogate delay adjustment
                                            mapdcUsed = true;
                                            mapdcStations(ist) = true;  % Mark for skipping post-processing
                                            line_debug('Using exact MAP/D/%d solver: Q=%.4f, W=%.4f', ...
                                                numServers, mapdcResult.meanQueueLength, mapdcResult.meanWaitingTime);
                                        elseif isFiniteCap
                                            % Finite-buffer FCFS. Detect exact M/M/c/K
                                            % (Poisson arrivals + same exp service across
                                            % classes); else use MMAP[K]/PH[K]/1/FCFS
                                            % truncate-and-renormalize approximation.
                                            mapdcUsed = false;
                                            capK = sn.cap(ist);
                                            finiteCapLossPerClass = [];
                                            [isMmck, muMmck] = mam_detect_mmck(sn, ist, K, aggrArrivalAtNode);
                                            if isMmck
                                                aggrLambdaTotal = sum(aggrLambda, 'omitnan');
                                                exactRes = qsys_mmck(aggrLambdaTotal, muMmck, sn.nservers(ist), capK);
                                                finiteCapMeanQ = exactRes.meanQueueLength;
                                                finiteCapLossProb = exactRes.lossProbability;
                                                finiteCapP0 = exactRes.queueLengthDist(1);
                                                line_debug('Using exact M/M/%d/%d: Q=%.4f, ploss=%.4f', ...
                                                    sn.nservers(ist), capK, finiteCapMeanQ, finiteCapLossProb);
                                            elseif sn.nservers(ist) == 1
                                                % Exact MMAP[K]/G/1/K with per-class loss ratio;
                                                % see _kb/06-solver-catalog.md for rationale
                                                svcMix = mam_svc_mixture({aggrArrivalAtNode{[1,3:end]}}, ...
                                                    {pie{ist}{:}}, {D0{ist,:}});
                                                exRes = qsys_mmapg1k(aggrArrivalAtNode{1}, ...
                                                    aggrArrivalAtNode(3:end), svcMix, capK);
                                                finiteCapMeanQ = exRes.meanQueueLength;
                                                finiteCapLossProb = exRes.lossAggregate;
                                                finiteCapLossPerClass = exRes.lossRatio;
                                                finiteCapP0 = exRes.p0;
                                                line_debug('Using exact MMAP/G/1/%d: Q=%.4f, ploss=%.4f, per-class ploss=%s', ...
                                                    capK, finiteCapMeanQ, finiteCapLossProb, mat2str(finiteCapLossPerClass, 4));
                                            else
                                                [meanQ_fc, lossProb_fc, p_norm_fc] = mam_truncate_renorm( ...
                                                    {aggrArrivalAtNode{[1,3:end]}}, {pie{ist}{:}}, {D0{ist,:}}, capK);
                                                finiteCapMeanQ = meanQ_fc;
                                                finiteCapLossProb = lossProb_fc;
                                                finiteCapP0 = p_norm_fc(1);
                                                line_debug('Using MMAP/PH/%d/%d truncate-renorm: Q=%.4f, ploss=%.4f', ...
                                                    sn.nservers(ist), capK, meanQ_fc, lossProb_fc);
                                            end
                                            finiteCapUsed = true;
                                            for k=1:K
                                                Qret{k} = NaN;  % filled per class in result loop
                                            end
                                            mapdcStations(ist) = true;  % skip surrogate delay 2nd pass
                                        elseif sn.isfunction(ist)
                                            % Open setup/delay-off via qbd_setupdelayoff; nodeparam
                                            % is NODE-indexed; see _kb/06-solver-catalog.md for rationale
                                            mapdcUsed = false;
                                            alpharate = map_lambda(sn.nodeparam{sn.stationToNode(ist)}{end}.setupTime);
                                            alphascv = map_scv(sn.nodeparam{sn.stationToNode(ist)}{end}.setupTime);
                                            betarate = map_lambda(sn.nodeparam{sn.stationToNode(ist)}{end}.delayoffTime);
                                            betascv = map_scv(sn.nodeparam{sn.stationToNode(ist)}{end}.delayoffTime);
                                            mu_k = zeros(1, K);
                                            lambda_k = zeros(1, K);
                                            active_k = false(1, K);
                                            for k=1:K
                                                pie_k = pie{ist}{k};
                                                if ~isempty(pie_k) && ~isnan(pie_k(1))
                                                    mu_k(k) = 1 / (pie_k * inv(-D0{ist,k}) * ones(size(pie_k))');
                                                    c = find(sn.chains(:,k), 1);
                                                    lambda_k(k) = rates{ist,c}(k);
                                                    active_k(k) = true;
                                                end
                                            end
                                            rho_k = lambda_k ./ mu_k;
                                            rho_k(~active_k) = 0;
                                            rho_total = sum(rho_k);
                                            if rho_total > 0
                                                % Aggregate service rate, split queue by load;
                                                % see _kb/06-solver-catalog.md for rationale
                                                aggrLambdaTotal = sum(aggrLambda, 'omitnan');
                                                aggrRate = aggrLambdaTotal / rho_total;
                                                Q_total = qbd_setupdelayoff(aggrLambdaTotal, aggrRate, alpharate, alphascv, betarate, betascv);
                                                for k=1:K
                                                    if active_k(k)
                                                        Qret{k} = Q_total * rho_k(k) / rho_total;
                                                    else
                                                        Qret{k} = 0;
                                                    end
                                                end
                                            else
                                                for k=1:K
                                                    Qret{k} = 0;
                                                end
                                            end
                                        else
                                            mapdcUsed = false;
                                            % MMAPPH1FCFS renewal marginal vs exact MAP/MAP/1 and
                                            % RAP/RAP/1; see _kb/06-solver-catalog.md for rationale
                                            isMEorRAPsvc = any(sn.procid(ist,:) == ProcessType.ME | ...
                                                sn.procid(ist,:) == ProcessType.RAP);
                                            useRapRap1 = false;
                                            if isMEorRAPsvc
                                                if (K == 1) && (sn.nservers(ist) == 1)
                                                    useRapRap1 = true;
                                                else
                                                    % Multi-class/server ME falls back to MMAPPH1FCFS
                                                    % with a warning; see _kb/06-solver-catalog.md for rationale
                                                    if ~meWarned(ist)
                                                        meWarned(ist) = true;
                                                        line_warning_always(mfilename, ...
                                                            'Station %s has a matrix-exponential or rational service process, which the RAP/RAP/1 analysis supports only with a single class at a single server (here %d classes, %g servers). Falling back to the phase-type approximation MMAPPH1FCFS, which is not exact for this service process.', ...
                                                            sn.nodenames{sn.stationToNode(ist)}, K, sn.nservers(ist));
                                                    end
                                                end
                                            end
                                            useMapMap1 = ~useRapRap1 && (K == 1) && (sn.nservers(ist) == 1) && ...
                                                (abs(map_acf(PH{ist}{1}, 1)) > GlobalConstants.CoarseTol);
                                            if useRapRap1
                                                % qbd_raprap1 gives QLen only, RN via Little;
                                                % see _kb/06-solver-catalog.md for rationale
                                                [~, QNrap] = qbd_raprap1({aggrArrivalAtNode{1}, aggrArrivalAtNode{3}}, ...
                                                    {PH{ist}{1}{1}, PH{ist}{1}{2}});
                                                Qret{1} = QNrap;
                                            elseif useMapMap1
                                                ql = Q_CT_MAP_MAP_1(aggrArrivalAtNode{1}, aggrArrivalAtNode{3}, ...
                                                    PH{ist}{1}{1}, PH{ist}{1}{2}, 'MaxNumComp', 100000);
                                                ql = ql(:);
                                                Qret{1} = sum((0:numel(ql)-1)' .* ql);
                                            else
                                                [Qret{1:K}] = MMAPPH1FCFS({aggrArrivalAtNode{[1,3:end]}}, {pie{ist}{:}}, {D0{ist,:}}, 'ncMoms', 1);
                                            end
                                        end
                                    else % all closed classes
                                        maxLevel = sum(N(isfinite(N)))+1;
                                        D = {aggrArrivalAtNode{[1,3:end]}};
                                        pdistr_k = cell(1,K);
                                        if map_lambda(D)< GlobalConstants.FineTol
                                            for k=1:K
                                                pdistr_k = [1-GlobalConstants.FineTol, GlobalConstants.FineTol];
                                                Qret{k} = GlobalConstants.FineTol / sn.rates(ist);
                                            end
                                        else
                                            if sn.isfunction(ist)
                                                alpharate = map_lambda(sn.nodeparam{sn.stationToNode(ist)}{end}.setupTime);
                                                betarate = map_lambda(sn.nodeparam{sn.stationToNode(ist)}{end}.delayoffTime);
                                                betascv = map_scv(sn.nodeparam{sn.stationToNode(ist)}{end}.delayoffTime);
                                                lambda_k = zeros(1, K);
                                                active_k = false(1, K);
                                                for k=1:K
                                                    pie_k = pie{ist}{k};
                                                    if ~isnan(pie_k(1))
                                                        c = find(sn.chains(:,k), 1);
                                                        lambda_k(k) = rates{ist,c}(k);
                                                        active_k(k) = true;
                                                    end
                                                end
                                                if any(active_k)
                                                    % Closed setup/delay-off: per-instance cold-start race,
                                                    % p_cold = LST_delayoff(1/ZT); see _kb/06-solver-catalog.md for rationale
                                                    setupMean = 1/alpharate;
                                                    if betascv == 1.0
                                                        % rate taken directly, as in qbd_setupdelayoff
                                                        betaProc = {-betarate, betarate};
                                                    else
                                                        betaProc = APH.fitMeanAndSCV(1/betarate, betascv).getProcess;
                                                    end
                                                    Tb = betaProc{1};
                                                    pie_b = map_pie(betaProc);
                                                    infstat = isinf(sn.nservers);
                                                    for k=1:K
                                                        if active_k(k) && lambda_k(k) > 0 && isfinite(S(ist,k))
                                                            c = find(sn.chains(:,k), 1);
                                                            inchain_k = find(sn.chains(c,:));
                                                            % ZT = chain think demand per visit;
                                                            % see _kb/06-solver-catalog.md for rationale
                                                            Vtot = sum(V(ist, inchain_k));
                                                            ZT = sum(Lchain(infstat, c), 'omitnan') / max(Vtot, GlobalConstants.FineTol);
                                                            nu = 1 / max(ZT, GlobalConstants.FineTol);
                                                            pcold = pie_b * ((nu*eye(size(Tb,1)) - Tb) \ (-Tb*ones(size(Tb,1),1)));
                                                            % service part S/nservers, delay correction adds rest;
                                                            % see _kb/06-solver-catalog.md for rationale
                                                            Qret{k} = lambda_k(k) * (pcold*setupMean + S(ist,k)/sn.nservers(ist));
                                                        elseif active_k(k)
                                                            % Zero-load class holds no jobs (NaN guard);
                                                            % see _kb/06-solver-catalog.md for rationale
                                                            Qret{k} = 0;
                                                        else
                                                            % Inactive class holds no jobs (NaN guard);
                                                            % see _kb/06-solver-catalog.md for rationale
                                                            Qret{k} = 0;
                                                        end
                                                    end
                                                else
                                                    for k=1:K
                                                        Qret{k} = NaN;
                                                    end
                                                end
                                            else
                                                % Capture all K per-class distributions;
                                                % see _kb/06-solver-catalog.md for rationale
                                                pdistr_all = cell(1,K);
                                                [pdistr_all{1:K}] = MMAPPH1FCFS(D, {pie{ist}{:}}, {D0{ist,:}}, 'ncDistr', maxLevel);
                                                for k=1:K
                                                    pdistr = pdistr_all{k};
                                                    pdistr_k = abs(pdistr(1:(N(k)+1)));
                                                    % Truncate at N(k), complement over truncated vector;
                                                    % see _kb/06-solver-catalog.md for rationale
                                                    pdistr_k(end) = abs(1-sum(pdistr_k(1:end-1)));
                                                    pdistr_k = pdistr_k / sum(pdistr_k(1:(N(k)+1)));
                                                    Qret{k} = max(0,min(N(k),(0:N(k))*pdistr_k(1:(N(k)+1))'));
                                                end
                                            end
                                        end
                                    end
                                else
                                    for k=1:K
                                        Qret{k} = sn.njobs(k);
                                    end
                                end
                            end
                            if finiteCapUsed
                                % Finite-cap per-class decomposition R_k = W_q + S_k;
                                % see _kb/06-solver-catalog.md for rationale
                                lambdaInflow = zeros(1, K);
                                for k=1:K
                                    cidx = find(sn.chains(:,k),1);
                                    lambdaInflow(k) = rates{ist,cidx}(k);
                                end
                                lambdaInflow(isnan(lambdaInflow)) = 0;
                                if ~isempty(finiteCapLossPerClass)
                                    % Exact MMAP/G/1/K branch: loss differs by class
                                    TN_eff = lambdaInflow .* (1 - finiteCapLossPerClass(1:K));
                                else
                                    TN_eff = lambdaInflow * (1 - finiteCapLossProb);
                                end
                                sumTN = sum(TN_eff);
                                if sumTN > 0
                                    Savg_eff = sum(TN_eff .* S(ist,1:K), 'omitnan') / sumTN;
                                    Wq = max(0, finiteCapMeanQ / sumTN - Savg_eff);
                                else
                                    Wq = 0;
                                end
                                for k=1:K
                                    TN(ist,k) = TN_eff(k);
                                    UN(ist,k) = TN(ist,k) * S(ist,k) / sn.nservers(ist);
                                    if TN(ist,k) > 0
                                        RN(ist,k) = Wq + S(ist,k);
                                        QN(ist,k) = TN(ist,k) * RN(ist,k);
                                    else
                                        RN(ist,k) = 0;
                                        QN(ist,k) = 0;
                                    end
                                end
                            else
                                QN(ist,:) = cell2mat(Qret);
                                for k=1:K
                                    c = find(sn.chains(:,k),1);
                                    TN(ist,k) = rates{ist,c}(k);
                                    UN(ist,k) = TN(ist,k) * S(ist,k) / sn.nservers(ist);
                                    if sn.isfunction(ist) && ~isfinite(UN(ist,k))
                                        % Unfreeze NaN-service util at setup station only;
                                        % see _kb/06-solver-catalog.md for rationale
                                        UN(ist,k) = 0;
                                    end
                                    QN(ist,k) = Qret{k};
                                    if mapdcUsed
                                        % For MAP/D/c, D/M/c, or PH/M/1, use exact results
                                        % (no surrogate delay adjustment).
                                        if isPhM1
                                            RN(ist,k) = phm1Result.meanSojournTime;
                                        elseif isDMc
                                            RN(ist,k) = dmcResult.meanSojournTime;
                                        else
                                            RN(ist,k) = mapdcResult.meanSojournTime;
                                        end
                                    else
                                        % Add surrogate-delay jobs, NaN-service guard;
                                        % see _kb/06-solver-catalog.md for rationale
                                        if isfinite(S(ist,k))
                                            QN(ist,k) = QN(ist,k) + TN(ist,k)*S(ist,k) * (sn.nservers(ist)-1)/sn.nservers(ist);
                                        end
                                        RN(ist,k) = QN(ist,k) ./ TN(ist,k);
                                    end
                                end
                            end
                    end
            end
        else % not a station
            switch sn.nodetype(ind)
                case NodeType.Fork
                    %                    line_error(mfilename,'Fork nodes not supported yet by MAM solver.');
            end
        end
    end
    %it
    %max(max(abs(TN-TN_1)))
end
totiter = it + 2;
CN = sum(RN,1);
QN = abs(QN);
for it=1:2 % second pass to rescale again QN based on RN correction
    for c=1:C
        inchain = sn.inchain{c};
        Nc = sum(sn.njobs(inchain));
        if isfinite(Nc)
            QNc = sum(sum(QN(:,inchain)));
            QN(:,inchain) = QN(:,inchain) * (Nc / QNc);
        end
        for ind=1:I
            for k=inchain
                if sn.isstation(ind)
                    ist = sn.nodeToStation(ind);
                    % Skip stations using exact MAP/D/c solver (already have exact values)
                    if mapdcStations(ist)
                        continue;
                    end
                    if V(ist,k)>0
                        if isinf(sn.nservers(ist))
                            RN(ist,k) = S(ist,k);
                        else
                            RN(ist,k) = max([S(ist,k), QN(ist,k) ./ TN(ist,k)]);
                        end
                    else
                        RN(ist,k) = 0;
                    end
                    QN(ist,k) = RN(ist,k) .* TN(ist,k);
                end
            end
        end
        Nc = sum(sn.njobs(inchain)); % closed population
        if Nc == 0 % if closed chain
            QN(:,c)=0;
            UN(:,c)=0;
            RN(:,c)=0;
            TN(:,c)=0;
            CN(c)=0;
            XN(c)=0;
        end
    end
end

% SLC clamp: closed-form leftover-capacity U_slc = 1 - sum_j U_j, applied last;
% see _kb/06-solver-catalog.md for rationale
slcAll = find(sn.isslc);
if ~isempty(slcAll)
    % Non-SLC util from declared (uninflated) service time;
    % see _kb/06-solver-catalog.md for rationale
    for ist=1:M
        if slcjobs(ist) > 0
            for k=1:K
                if ~sn.isslc(k)
                    UN(ist,k) = Strue(ist,k)*TN(ist,k);
                end
            end
        end
    end
    QN(:,slcAll) = 0;
    UN(:,slcAll) = 0;
    RN(:,slcAll) = 0;
    TN(:,slcAll) = 0;
    for ist=1:M
        slck = slcAll(sn.refstat(slcAll) == ist);
        if isempty(slck)
            continue;
        end
        if isinf(sn.nservers(ist))
            % Delay station: no contention, every customer is always in
            % service, so the class completes at its full aggregate rate.
            for k=slck(:)'
                QN(ist,k) = sn.njobs(k);
                TN(ist,k) = sn.njobs(k)*sn.rates(ist,k);
                RN(ist,k) = Strue(ist,k);
                UN(ist,k) = Strue(ist,k)*TN(ist,k);
            end
        else
            nsrv = sn.nservers(ist);
            Uleft = max(0, 1 - sum(UN(ist,~sn.isslc)));
            % Several self-looping classes at one station share the free
            % capacity in proportion to the service rate they offer.
            w = sn.njobs(slck(:)').*sn.rates(ist,slck(:)');
            if sum(w) <= 0
                continue;
            end
            for i=1:numel(slck)
                k = slck(i);
                UN(ist,k) = min(Uleft*w(i)/sum(w), min(sn.njobs(k),nsrv)/nsrv);
                TN(ist,k) = sn.rates(ist,k)*UN(ist,k)*nsrv;
                QN(ist,k) = sn.njobs(k);
                if TN(ist,k) > 0
                    RN(ist,k) = QN(ist,k) ./ TN(ist,k);
                else
                    RN(ist,k) = 0;
                end
            end
        end
    end
    CN = sum(RN,1);
end
end
