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

for ist=1:sn.nstations
    switch sn.sched(ist)
        case SchedStrategy.EXT
            % no-op
        case {SchedStrategy.FCFS, SchedStrategy.HOL, SchedStrategy.FCFSPRPRIO}
            % no-op
        otherwise
            if options.verbose
                line_warning(mfilename,'The dec.mmap method does not support non-FCFS queues.\n');
            end
            [QN,UN,RN,TN,CN,XN] = deal([],[],[],[],[],[]);
            totiter = 0;
            method = '';
            runtime = toc(Tstart);
            return
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
            case {SchedStrategy.FCFS, SchedStrategy.HOL, SchedStrategy.FCFSPRPRIO}
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

    it_max = options.iter_max;
    for it=1:it_max
        %it
        %        now estimate arrival processes
        if it == 1
            %            initially form departure processes using scaled service
            DEP = PH;
            for ind=1:M
                for r=1:K
                    ist = sn.nodeToStation(ind);
                    DEP{ind,r} = map_scale(PH{ist}{r}, 1 / (lambda(r) * V(ind,r)) );
                end
            end
        end

        ARV = solver_mam_traffic(sn, DEP, config);

        QN_1 = QN;
        for ist=1:M
            ind = sn.stationToNode(ist);
            switch sn.nodetype(ind)
                case NodeType.Queue
                    if length(ARV{ind}{1}) > config.space_max
                        line_printf('\nArrival process at node %d is now at %d states. Compressing.',ind,length(ARV{ind}{1}));
                        ARV{ind} = mmap_compress(ARV{ind});
                    end
                    [Qret{1:K}, ~] = MMAPPH1FCFS({ARV{ind}{[1,3:end]}}, {pie{ist}{:}}, {D0{ist,:}}, 'ncMoms', 1, 'ncDistr',2);
                    for k=1:K
                        QN(ist,k) = sum(Qret{k});
                    end
                    TN(ist,:) = mmap_lambda(ARV{ind});
            end
            for k=1:K
                UN(ist,k) = TN(ist,k) * map_mean(PH{ist}{k});
                %add number of jobs at the surrogate delay server
                QN(ist,k) = QN(ist,k) + TN(ist,k)*(map_mean(PH{ist}{k})*sn.nservers(ist)) * (sn.nservers(ist)-1)/sn.nservers(ist);
                RN(ist,k) = QN(ist,k) ./ TN(ist,k);
            end
        end

        if it >=3 && max(abs(QN(:)-QN_1(:))./QN_1(:)) < options.iter_tol
            break;
        end

        for ist=1:M
            ind = sn.stationToNode(ist);
            switch sn.nodetype(ind)
                case NodeType.Queue
                    [Ret{1:2*K}] = MMAPPH1FCFS({ARV{ind}{[1,3:end]}}, {pie{ist}{:}}, {D0{ist,:}}, 'stDistrPH');
                    for r=1:K
                        %obtain response time distribution for class r
                        %alpha = Ret{(r-1)*2+1}; T0 = Ret{(r-1)*2+2};
                        %RD = {T0,(-T0)*ones(length(alpha),1)*alpha(:)'}; %PH to MAP
                        %tRD = sum(RD{2},2);
                        %pieRD = map_pie(RD);

                        %define a ph for the arrival process of class r
                        A = mmap_hide(ARV{ind},setdiff(1:K,r));
                        tA = sum(A{2},2);
                        pieA = map_pie(A);

                        %define a ph for the service process of class r
                        S = PH{ist}{r};
                        pieS = map_pie(S);
                        tS = sum(S{2},2);

                        %with probability rho the output is the
                        %inter-arrival time+service time, with probability 1-rho is the
                        %response time
                        rho = sum(UN(ist,:));
                        AQ = sum(QN(ist,:)); % queue seen upon arrival with PASTA
                        Afull = AQ/rho; % from AQ = (1-rho)*0 + rho*AQfull
                        pfullonarrival = 1 - (Afull)/(1+Afull);

                        A = mmap_scale(A,map_mean(A)-map_mean(S));

                        zAS = 0*tA*pieS;
                        zSA = 0*tS*pieA;
                        zA = 0*A{2};

                        % this is only for single-class
                        DEP0ir = [S{1}, sum(rho)*tS*pieA;
                            zAS, A{1}];

                        DEP1ir = [(1-sum(rho))*S{2}, zSA;
                            tA*pieS, zA];

                        DEP{ind,r} = map_normalize({DEP0ir,DEP1ir});

                        DEP{ind,r} = map_scale(DEP{ind,r}, 1 / (lambda(r) * V(ind,r)) );
                        SCVd(ind,r) = map_scv(DEP{ind,r});
                        IDCd(ind,r) = map_idc(DEP{ind,r});
                    end
            end
        end
    end
    totiter = it;
    if options.verbose
        line_printf('\nMAM parametric decomposition completed in %d iterations.',it);
    end
else
    if options.verbose
        line_warning(mfilename,'This model is not supported by SolverMAM yet. Returning with no result.\n');
    end
end
runtime = toc(Tstart);
end
