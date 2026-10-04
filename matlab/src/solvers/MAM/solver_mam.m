function [QN,UN,RN,TN,CN,XN,totiter,method,runtime] = solver_mam(sn, options)
%[Q,U,R,T,C,X,totiter] = SOLVER_MAM(QN, PH, OPTIONS)

%Copyright (c) 2012-2026, Imperial College London
%All rights reserved.

method = options.method;
config = options.config;
totiter = NaN;
PH = sn.proc;
I = sn.nnodes;
M = sn.nstations;
K = sn.nclasses;
C = sn.nchains;
N = sn.njobs';
V = cellsum(sn.visits);
Tstart=tic;
QN = zeros(M,K);
UN = zeros(M,K);
RN = zeros(M,K);
TN = zeros(M,K);
CN = zeros(1,K);
XN = zeros(1,K);

lambda = zeros(1,K);
for c=1:C
    inchain = sn.inchain{c};
    lambdas_inchain = sn.rates(sn.refstat(inchain(1)),inchain);
    lambdas_inchain = lambdas_inchain(isfinite(lambdas_inchain));
    lambda(inchain) = sum(lambdas_inchain);
end

chain = zeros(1,K);
for k=1:K
    chain(k) = find(sn.chains(:,k));
end

% THE TWO RESTRICTIONS OF THIS ANALYZER, RAISED RATHER THAN RETURNED. Both used
% to end in a warning and a result: an unsupported discipline returned empty
% matrices and a closed model fell through the branch below with QN..XN still
% at their zero initialization, so SolverMAM reported an entirely zero table as
% if it were the answer. SolverMAM.getMethodFeatureSet states the same two rules
% declaratively -- 'dec.mmap' drops SchedStrategy_INF, ClosedClass and
% SelfLoopingClass -- so findSolver no longer offers the pair; this is what a
% caller who names the method by hand meets.
for ist=1:sn.nstations
    switch sn.sched(ist)
        case SchedStrategy.EXT
            % no-op
        case {SchedStrategy.FCFS, SchedStrategy.HOL, SchedStrategy.FCFSPRPRIO, SchedStrategy.PS}
            % no-op
        otherwise
            line_error(mfilename, sprintf(['The dec.mmap method does not support the %s ' ...
                'scheduling strategy at station %d: the departure-process fixed point is ' ...
                'built for EXT, FCFS, HOL, FCFSPRPRIO and PS stations only. Use the ' ...
                'dec.source method.'], SchedStrategy.toText(sn.sched(ist)), ist));
    end
end

if all(isinf(sn.njobs)) % is open
    %    open queueing system (one node is the external world)
    pie = {};
    D0 = {};
    for ist=1:M
        switch sn.sched(ist)
            case SchedStrategy.EXT
                TN(ist,:) = sn.rates(ist,:);
                TN(ist,isnan(TN(ist,:)))=0;
            case {SchedStrategy.FCFS, SchedStrategy.HOL, SchedStrategy.FCFSPRPRIO, SchedStrategy.PS}
                for k=1:K
                    %                    divide service time by number of servers and put
                    %                    later a surrogate delay server in tandem to compensate
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

    % departure-process fixed point (parametric decomposition), driven on
    % the station queue lengths by the generic DA driver
    DEP = {};
    fpopts = options;
    fpopts.config.da_miniter = 3; % legacy loop tested convergence only from the third sweep
    fpopts.config.da_norm = @(xn,xr) max(abs(xn(:)-xr(:))./xr(:)); % relative difference
    [~, totiter] = da_fpi(@mam_dec_sweep, QN, fpopts);
    if options.verbose
        line_printf('\nMAM parametric decomposition completed in %d iterations.',totiter);
    end
else
    line_error(mfilename, ['The dec.mmap method supports open models only: the ' ...
        'departure-process fixed point iterates on arrival streams that a closed ' ...
        'population does not have. Use the dec.source method, or method ''default'', ' ...
        'which routes a closed model to an analyzer that solves it.']);
end
runtime = toc(Tstart);

    function [xnew, xref] = mam_dec_sweep(~, itnum)
        %it
        %        now estimate arrival processes
        if itnum == 1
            % initially form departure processes using scaled service; DEP/PH/V
            % are STATION-indexed; see _kb/06-solver-catalog.md for rationale
            DEP = PH;
            for ist=1:M
                for r=1:K
                    DEP{ist,r} = map_scale(PH{ist}{r}, 1 / (lambda(r) * V(ist,r)) );
                end
            end
        end

        ARV = solver_mam_traffic(sn, DEP, config);

        xref = QN;
        for ist=1:M
            ind = sn.stationToNode(ist);
            finiteCapUsed = false;
            switch sn.nodetype(ind)
                case NodeType.Queue
                    if length(ARV{ind}{1}) > config.space_max
                        line_printf('\nArrival process at node %d is now at %d states. Compressing.',ind,length(ARV{ind}{1}));
                        ARV{ind} = mmap_compress(ARV{ind});
                    end
                    TN(ist,:) = mmap_lambda(ARV{ind});
                    switch sn.sched(ist)
                        case {SchedStrategy.FCFS, SchedStrategy.HOL, SchedStrategy.FCFSPRPRIO}
                            isFiniteCap = isfinite(sn.cap(ist));
                            if isFiniteCap
                                capK = sn.cap(ist);
                                [isMmck, muMmck] = mam_detect_mmck(sn, ist, K, ARV{ind});
                                if isMmck
                                    aggrLambda_ist = sum(TN(ist,:), 'omitnan');
                                    exactRes = qsys_mmck(aggrLambda_ist, muMmck, sn.nservers(ist), capK);
                                    meanQ_fc = exactRes.meanQueueLength;
                                    lossProb_fc = exactRes.lossProbability;
                                else
                                    [meanQ_fc, lossProb_fc, ~] = mam_truncate_renorm( ...
                                        {ARV{ind}{[1,3:end]}}, {pie{ist}{:}}, {D0{ist,:}}, capK);
                                end
                                lambdaInflow = TN(ist,:);  % per-class rates from mmap_lambda
                                lambdaInflow(isnan(lambdaInflow)) = 0;
                                TN_eff = lambdaInflow * (1 - lossProb_fc);
                                sumTN = sum(TN_eff);
                                % Actual per-class service mean (PH was scaled by 1/c)
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
                                [Qret{1:K}, ~] = MMAPPH1FCFS({ARV{ind}{[1,3:end]}}, {pie{ist}{:}}, {D0{ist,:}}, 'ncMoms', 1, 'ncDistr',2);
                                for k=1:K
                                    QN(ist,k) = sum(Qret{k});
                                end
                            end
                        case SchedStrategy.PS
                            for k=1:K
                                UN(ist,k) = TN(ist,k) * map_mean(PH{ist}{k});
                            end
                            Uden = min([1-GlobalConstants.FineTol, sum(UN(ist,:))]);
                            for k=1:K
                                QN(ist,k) = UN(ist,k)/(1-Uden);
                            end
                    end
            end
            if ~finiteCapUsed
                for k=1:K
                    UN(ist,k) = TN(ist,k) * map_mean(PH{ist}{k});
                    %add number of jobs at the surrogate delay server
                    QN(ist,k) = QN(ist,k) + TN(ist,k)*(map_mean(PH{ist}{k})*sn.nservers(ist)) * (sn.nservers(ist)-1)/sn.nservers(ist);
                    RN(ist,k) = QN(ist,k) ./ TN(ist,k);
                end
            end
        end

        for ist=1:M
            ind = sn.stationToNode(ist);
            switch sn.nodetype(ind)
                case NodeType.Queue
                    for r=1:K
                        % extract class-r arrival MAP
                        A = mmap_hide(ARV{ind},setdiff(1:K,r));
                        S = PH{ist}{r};
                        na = length(A{1});
                        ns = length(S{1});
                        etaqa_n = config.etaqa_trunc;
                        etaqa_sz = (etaqa_n+1)*na*ns;
                        rho = sum(UN(ist,:));
                        % use ETAQA if state space is manageable and queue is stable
                        if etaqa_sz <= config.space_max && rho < 1-GlobalConstants.FineTol
                            try
                                switch sn.sched(ist)
                                    case {SchedStrategy.FCFS, SchedStrategy.HOL, SchedStrategy.FCFSPRPRIO}
                                        DEP{ist,r} = qbd_depproc_etaqa(A, S, etaqa_n);
                                    case SchedStrategy.PS
                                        DEP{ist,r} = qbd_depproc_etaqa_ps(A, S, etaqa_n);
                                end
                                DEP{ist,r} = map_normalize(DEP{ist,r});
                            catch
                                % fall back to scaled service on ETAQA failure
                                DEP{ist,r} = PH{ist}{r};
                            end
                        else
                            DEP{ist,r} = PH{ist}{r};
                        end
                        DEP{ist,r} = map_scale(DEP{ist,r}, 1 / (lambda(r) * V(ist,r)) );
                        SCVd(ist,r) = map_scv(DEP{ist,r});
                        IDCd(ist,r) = map_idc(DEP{ist,r});
                    end
            end
        end
    xnew = QN;
    end
end
