function [QN,UN,RN,TN,CN,XN,totiter] = solver_mam_basic_mmap_inner(sn, options, lambda)
% [QN,UN,RN,TN,CN,XN,TOTITER] = SOLVER_MAM_BASIC_MMAP_INNER(SN, OPTIONS, LAMBDA)
%
% MAM/MMAP fork-join decomposition algorithm parameterised by per-class arrival
% rates LAMBDA. Performs the departure-process refinement loop only;
% population enforcement (closed networks) is the wrapper's responsibility.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

config = options.config;
if ~isfield(config, 'fj_sync_q_len')
    config.fj_sync_q_len = 2;
end
if ~isfield(config, 'etaqa_trunc')
    config.etaqa_trunc = 8;
end

PH = sn.proc;
I = sn.nnodes;
M = sn.nstations;
K = sn.nclasses;
V = cellsum(sn.visits);
S = 1./sn.rates;

QN = zeros(M,K);
UN = zeros(M,K);
RN = zeros(M,K);
TN = zeros(M,K);
CN = zeros(1,K);
XN = zeros(1,K);

% Build FJ synchronization map
fjSyncMap = sn_build_fj_sync_map(sn);

% Prepare PH service distributions
pie = {};
D0 = {};
for ist=1:M
    switch sn.sched(ist)
        case SchedStrategy.EXT
            TN(ist,:) = sn.rates(ist,:);
            TN(ist,isnan(TN(ist,:))) = 0;
        case {SchedStrategy.FCFS, SchedStrategy.HOL, SchedStrategy.FCFSPRPRIO}
            for k=1:K
                PH{ist}{k} = map_scale(PH{ist}{k}, map_mean(PH{ist}{k})/sn.nservers(ist));
                pie{ist}{k} = map_pie(PH{ist}{k});
                D0{ist,k} = PH{ist}{k}{1};
                if any(isnan(D0{ist,k}))
                    D0{ist,k} = -GlobalConstants.Immediate;
                    pie{ist}{k} = 1;
                    PH{ist}{k} = map_exponential(GlobalConstants.Immediate);
                end
            end
        case SchedStrategy.INF
            for k=1:K
                pie{ist}{k} = map_pie(PH{ist}{k});
                D0{ist,k} = PH{ist}{k}{1};
                if any(isnan(D0{ist,k}))
                    D0{ist,k} = -GlobalConstants.Immediate;
                    pie{ist}{k} = 1;
                    PH{ist}{k} = map_exponential(GlobalConstants.Immediate);
                end
            end
        case SchedStrategy.PS
            for k=1:K
                PH{ist}{k} = map_scale(PH{ist}{k}, map_mean(PH{ist}{k})/sn.nservers(ist));
                pie{ist}{k} = map_pie(PH{ist}{k});
                D0{ist,k} = PH{ist}{k}{1};
                if any(isnan(D0{ist,k}))
                    D0{ist,k} = -GlobalConstants.Immediate;
                    pie{ist}{k} = 1;
                    PH{ist}{k} = map_exponential(GlobalConstants.Immediate);
                end
            end
    end
end

% departure-process fixed point (FJ parametric decomposition), driven on
% the station queue lengths by the generic DA driver
DEP = {};
fpopts = options;
fpopts.config.da_miniter = 3; % legacy loop tested convergence only from the third sweep
fpopts.config.da_norm = @(xn,xr) max(abs(xn(:)-xr(:))./(xr(:)+GlobalConstants.FineTol)); % relative difference
[~, totiter] = da_fpi(@mmap_dec_sweep, QN, fpopts);
if options.verbose
    line_printf('\nMAM FJ parametric decomposition completed in %d iterations.', totiter);
end

% Join: derive QN/RN from parallel branch means
for joinIdx = find(sn.nodetype == NodeType.Join)'
    joinStat = sn.nodeToStation(joinIdx);
    if isnan(joinStat)
        continue;
    end
    syncGroups = unique(fjSyncMap.nodeSync(joinIdx, :));
    syncGroups = syncGroups(syncGroups > 0);
    for r = 1:K
        if TN(joinStat, r) <= 0
            continue;
        end
        syncDelay = 0;
        joinArrivalRate = 0;
        for gid = syncGroups
            branchNodes = find(fjSyncMap.nodeSync(joinIdx, :) == gid);
            branchRt = zeros(1, numel(branchNodes));
            branchTput = zeros(1, numel(branchNodes));
            used = 0;
            for b = 1:numel(branchNodes)
                branchStat = sn.nodeToStation(branchNodes(b));
                if isnan(branchStat) || RN(branchStat, r) <= 0
                    continue;
                end
                used = used + 1;
                branchRt(used) = RN(branchStat, r);
                branchTput(used) = TN(branchStat, r);
            end
            branchRt = branchRt(1:used);
            branchTput = branchTput(1:used);
            if numel(branchRt) < 2
                continue;
            end
            lambdai = 1 ./ branchRt;
            maxBranchRt = 0;
            for pow = 0:(numel(branchRt) - 1)
                maxBranchRt = maxBranchRt + (-1)^pow * sum(1 ./ sum(nchoosek(lambdai, pow + 1), 2));
            end
            syncDelay = syncDelay + max(maxBranchRt - mean(branchRt), 0);
            joinArrivalRate = joinArrivalRate + sum(branchTput);
        end
        RN(joinStat, r) = syncDelay;
        QN(joinStat, r) = joinArrivalRate * syncDelay;
        UN(joinStat, r) = 0;
    end
end

CN = sum(RN,1);
QN(isnan(QN)) = 0;
RN(isnan(RN)) = 0;
UN(isnan(UN)) = 0;
TN(isnan(TN)) = 0;

    function [xnew, xref] = mmap_dec_sweep(~, itnum)
    % Initialize departure processes (node-indexed: DEP{ind,r})
    if itnum == 1
        DEP = cell(I,K);
        for ind=1:I
            isForkJoin = (sn.nodetype(ind) == NodeType.Fork || sn.nodetype(ind) == NodeType.Join);
            if sn.isstation(ind) && ~isForkJoin
                ist = sn.nodeToStation(ind);
                for r=1:K
                    if V(ist,r) > 0 && lambda(r) > 0
                        DEP{ind,r} = map_scale(PH{ist}{r}, 1 / (lambda(r) * V(ist,r)));
                    elseif sn.isslc(r)
                        % SLC vanishing-rate departure process; see _kb/06-solver-catalog.md for rationale
                        DEP{ind,r} = map_exponential(1/GlobalConstants.Zero);
                    else
                        DEP{ind,r} = PH{ist}{r};
                    end
                end
            else
                for r=1:K
                    if lambda(r) > 0
                        DEP{ind,r} = map_exponential(1/lambda(r));
                    else
                        DEP{ind,r} = map_exponential(1/GlobalConstants.Immediate);
                    end
                end
            end
        end
    end

    % Compute arrival processes with FJ synchronization
    ARV = solver_mam_traffic_mmap(sn, DEP, config, fjSyncMap);

    xref = QN;
    for ist=1:M
        ind = sn.stationToNode(ist);
        switch sn.nodetype(ind)
            case NodeType.Join
                for k=1:K
                    TN(ist,k) = lambda(k);
                    UN(ist,k) = 0;
                    QN(ist,k) = 0;
                    RN(ist,k) = 0;
                end
            case NodeType.Queue
                if ~isempty(ARV{ind}) && iscell(ARV{ind})
                    if length(ARV{ind}{1}) > config.space_max
                        if options.verbose
                            line_printf('\nArrival process at node %d is now at %d states. Compressing.', ind, length(ARV{ind}{1}));
                        end
                        ARV{ind} = mmap_compress(ARV{ind}, config);
                    end

                    finiteCapUsed = false;
                    switch sn.sched(ist)
                        case {SchedStrategy.FCFS, SchedStrategy.HOL, SchedStrategy.FCFSPRPRIO}
                            isFiniteCap = isfinite(sn.cap(ist));
                            if isFiniteCap
                                capK = sn.cap(ist);
                                [isMmck, muMmck] = mam_detect_mmck(sn, ist, K, ARV{ind});
                                if isMmck
                                    aggrLambda_ist = sum(mmap_lambda(ARV{ind}), 'omitnan');
                                    exactRes = qsys_mmck(aggrLambda_ist, muMmck, sn.nservers(ist), capK);
                                    meanQ_fc = exactRes.meanQueueLength;
                                    lossProb_fc = exactRes.lossProbability;
                                else
                                    [meanQ_fc, lossProb_fc, ~] = mam_truncate_renorm( ...
                                        {ARV{ind}{[1,3:end]}}, {pie{ist}{:}}, {D0{ist,:}}, capK);
                                end
                                lambdaInflow = mmap_lambda(ARV{ind});
                                lambdaInflow(isnan(lambdaInflow)) = 0;
                                TN_eff = lambdaInflow * (1 - lossProb_fc);
                                sumTN = sum(TN_eff);
                                S_actual = zeros(1, K);
                                for k=1:K
                                    S_actual(k) = map_mean(PH{ist}{k}) * sn.nservers(ist);
                                end
                                if sumTN > 0
                                    Savg_eff = sum(TN_eff .* S_actual, 'omitnan') / sumTN;
                                    Wq = max(0, meanQ_fc / sumTN - Savg_eff);
                                else
                                    Wq = 0;
                                end
                                for k=1:K
                                    TN(ist,k) = TN_eff(k);
                                    UN(ist,k) = TN(ist,k) * map_mean(PH{ist}{k});
                                    if TN(ist,k) > 0
                                        RN(ist,k) = Wq + S_actual(k);
                                        QN(ist,k) = TN(ist,k) * RN(ist,k);
                                    else
                                        RN(ist,k) = 0;
                                        QN(ist,k) = 0;
                                    end
                                end
                                finiteCapUsed = true;
                            else
                                rho_ist_classes = mmap_lambda(ARV{ind}) .* arrayfun(@(k) map_mean(PH{ist}{k}), 1:K);
                                rho_ist_classes(isnan(rho_ist_classes)) = 0;
                                % Exclude SLC from FCFS saturation test; see _kb/06-solver-catalog.md for rationale
                                rho_ist = sum(rho_ist_classes(~sn.isslc));
                                if rho_ist < 1 - GlobalConstants.FineTol
                                    % Exact MAP/MAP/1 for correlated single-class service;
                                    % see _kb/06-solver-catalog.md for rationale
                                    useMapMap1 = (K == 1) && (sn.nservers(ist) == 1) && ...
                                        (abs(map_acf(PH{ist}{1}, 1)) > GlobalConstants.CoarseTol);
                                    if useMapMap1
                                        Carv0 = ARV{ind}{1};
                                        Carv1 = ARV{ind}{3};
                                        ql = Q_CT_MAP_MAP_1(Carv0, Carv1, PH{ist}{1}{1}, PH{ist}{1}{2}, 'MaxNumComp', 100000);
                                        ql = ql(:);
                                        QN(ist,1) = sum((0:numel(ql)-1)' .* ql);
                                    else
                                        [Qret{1:K}, ~] = MMAPPH1FCFS({ARV{ind}{[1,3:end]}}, {pie{ist}{:}}, {D0{ist,:}}, 'ncMoms', 1, 'ncDistr', 2);
                                        for k=1:K
                                            QN(ist,k) = sum(Qret{k});
                                        end
                                    end
                                else
                                    % Bound queue lengths under overload (no NaN);
                                    % see _kb/06-solver-catalog.md for rationale
                                    for k=1:K
                                        if isfinite(sn.njobs(k))
                                            QN(ist,k) = sn.njobs(k);
                                        else
                                            QN(ist,k) = 1/GlobalConstants.FineTol;
                                        end
                                    end
                                end
                                TN(ist,:) = mmap_lambda(ARV{ind});
                            end
                        case SchedStrategy.PS
                            TN(ist,:) = mmap_lambda(ARV{ind});
                            for k=1:K
                                UN(ist,k) = TN(ist,k) * S(ist,k);
                            end
                            % Exclude SLC from PS sharing denominator; see _kb/06-solver-catalog.md for rationale
                            Uden = min([1-GlobalConstants.FineTol, sum(UN(ist,~sn.isslc))]);
                            for k=1:K
                                QN(ist,k) = UN(ist,k)/(1-Uden);
                            end
                    end

                    if ~finiteCapUsed
                        for k=1:K
                            UN(ist,k) = TN(ist,k) * map_mean(PH{ist}{k});
                            QN(ist,k) = QN(ist,k) + TN(ist,k)*(map_mean(PH{ist}{k})*sn.nservers(ist)) * (sn.nservers(ist)-1)/sn.nservers(ist);
                            RN(ist,k) = QN(ist,k) ./ TN(ist,k);
                        end
                    end
                end
            otherwise
                switch sn.sched(ist)
                    case SchedStrategy.INF
                        if ~isempty(ARV{ind}) && iscell(ARV{ind})
                            TN(ist,:) = mmap_lambda(ARV{ind});
                        end
                        for k=1:K
                            if TN(ist,k) > 0
                                UN(ist,k) = S(ist,k)*TN(ist,k);
                                QN(ist,k) = TN(ist,k)*S(ist,k);
                                RN(ist,k) = S(ist,k);
                            end
                        end
                    case SchedStrategy.EXT
                        % Source: TN already set above
                end
        end
    end

    % Update departure processes
    for ist=1:M
        ind = sn.stationToNode(ist);
        switch sn.nodetype(ind)
            case NodeType.Queue
                if ~isempty(ARV{ind}) && iscell(ARV{ind})
                    switch sn.sched(ist)
                        case {SchedStrategy.FCFS, SchedStrategy.HOL, SchedStrategy.FCFSPRPRIO}
                            for r=1:K
                                A = mmap_hide(ARV{ind}, setdiff(1:K,r));
                                Srv = PH{ist}{r};
                                na = length(A{1});
                                ns = length(Srv{1});
                                etaqa_n = config.etaqa_trunc;
                                etaqa_sz = (etaqa_n+1)*na*ns;
                                rho = sum(UN(ist,:));
                                if etaqa_sz <= config.space_max && rho < 1-GlobalConstants.FineTol
                                    try
                                        DEP{ind,r} = qbd_depproc_etaqa(A, Srv, etaqa_n);
                                        DEP{ind,r} = map_normalize(DEP{ind,r});
                                    catch
                                        DEP{ind,r} = Srv;
                                    end
                                else
                                    DEP{ind,r} = Srv;
                                end
                                if V(ist,r) > 0 && lambda(r) > 0
                                    DEP{ind,r} = map_scale(DEP{ind,r}, 1 / (lambda(r) * V(ist,r)));
                                end
                            end
                        case SchedStrategy.PS
                            for r=1:K
                                A = mmap_hide(ARV{ind}, setdiff(1:K,r));
                                Srv = PH{ist}{r};
                                na = length(A{1});
                                ns = length(Srv{1});
                                etaqa_n = config.etaqa_trunc;
                                etaqa_sz = (etaqa_n+1)*na*ns;
                                rho = sum(UN(ist,:));
                                if V(ist,r) > 0 && lambda(r) > 0
                                    if etaqa_sz <= config.space_max && rho < 1-GlobalConstants.FineTol
                                        try
                                            DEP{ind,r} = qbd_depproc_etaqa_ps(A, Srv, etaqa_n);
                                            DEP{ind,r} = map_normalize(DEP{ind,r});
                                        catch
                                            DEP{ind,r} = Srv;
                                        end
                                    else
                                        DEP{ind,r} = Srv;
                                    end
                                    DEP{ind,r} = map_scale(DEP{ind,r}, 1 / (lambda(r) * V(ist,r)));
                                end
                            end
                    end
                end
            case NodeType.Join
                for r=1:K
                    if TN(ist,r) > 0
                        DEP{ind,r} = map_exponential(1/TN(ist,r));
                    end
                end
        end
    end
    xnew = QN;
    end
end
