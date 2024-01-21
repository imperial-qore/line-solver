function [QN,UN,RN,TN,CN,XN] = solver_mam_basic(sn, options)
% [Q,U,R,T,C,X] = SOLVER_MAM_BASIC(QN, PH, OPTIONS)

% This solver uses MAM to solve queues in isolation, but simplifies the
% traffic equations by using visits to rescale the flows into the queue
% inputs.
%
% Copyright (c) 2012-2024, Imperial College London
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
Lchain = snGetDemandsChain(sn);

QN = zeros(M,K);
UN = zeros(M,K);
RN = zeros(M,K);
TN = zeros(M,K);
CN = zeros(1,K);
XN = zeros(1,K);

pie = {};
D0 = {};

lambda = zeros(1,C);
chainSysArrivals = cell(1,C);
TN_1 = TN+Inf;

it = 0;

% open queueing system (one node is the external world)
% first build the joint arrival process
for ist=1:M
    switch sn.schedid(ist)
        case {SchedStrategy.ID_FCFS, SchedStrategy.ID_HOL, SchedStrategy.ID_PS}
            for k=1:K
                % divide service time by number of servers and put
                % later a surrogate delay server in tandem to compensate
                PH{ist}{k} = map_scale(PH{ist}{k}, S(ist,k)/sn.nservers(ist));
                pie{ist}{k} = map_pie(PH{ist}{k});
                D0{ist,k} = PH{ist}{k}{1};
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
if isMixed
    % treat open as having higher priority than closed
    %sn.classprio(~isfinite(sn.njobs)) = 1 + max(sn.classprio(isfinite(sn.njobs)));
end

lambdas_inchain = cell(1,C);
for c=1:C
    inchain = sn.inchain{c};
    lambdas_inchain{c} = sn.rates(sn.refstat(inchain(1)),inchain);
    %lambdas_inchain{c} = lambdas_inchain{c}(isfinite(lambdas_inchain{c}));
    lambda(c) = sum(lambdas_inchain{c}(isfinite(lambdas_inchain{c})));
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

while max(max(abs(TN-TN_1))) > tol && it <= options.iter_max %#ok<max>
    it = it + 1;
    TN_1 = TN;
    Umax = max(sum(UN(sd,:),2));
    if Umax > 1
        lambda = lambda * 1/Umax;
    else
        for c=1:C
            inchain = sn.inchain{c};
            if ~isinf(sum(sn.njobs(inchain))) % if closed chain
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

    for c=1:C
        inchain = sn.inchain{c};
        chainSysArrivals{c} = mmap_exponential(sum(lambda(c)));
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
                    switch sn.schedid(ist)
                        case SchedStrategy.ID_INF
                            for c=1:C
                                inchain = sn.inchain{c};
                                for k=inchain
                                    TN(ist,k) = lambda(c)*V(ist,k);
                                    UN(ist,k) = S(ist,k)*TN(ist,k);
                                    QN(ist,k) = TN(ist,k).*S(ist,k)*V(ist,k);
                                    RN(ist,k) = QN(ist,k)/TN(ist,k);
                                end
                            end
                        case SchedStrategy.ID_PS
                            for c=1:C
                                inchain = sn.inchain{c};
                                for k=inchain
                                    TN(ist,k) = lambda(c)*V(ist,k);
                                    UN(ist,k) = S(ist,k)*TN(ist,k);
                                end
                                %Nc = sum(sn.njobs(inchain)); % closed population
                                Uden = min([1-GlobalConstants.FineTol,sum(UN(ist,:))]);
                                for k=inchain
                                    %QN(ist,k) = (UN(ist,k)-UN(ist,k)^(Nc+1))/(1-Uden); % geometric bound type approximation
                                    QN(ist,k) = UN(ist,k)/(1-Uden);
                                    RN(ist,k) = QN(ist,k)/TN(ist,k);
                                end
                            end
                        case {SchedStrategy.ID_FCFS, SchedStrategy.ID_HOL}
                            chainArrivalAtNode = cell(1,C);
                            rates = cell(M,C);
                            for c=1:C %for each chain
                                rates{ist,c} = V(ist,:) .* lambda(c); % visits of classes within the chain
                                inchain = find(sn.chains(c,:))';
                                markProb = rates{ist,c}(inchain) / sum(rates{ist,c}(inchain));
                                markProb(isnan(markProb)) = 0;
                                chainArrivalAtNode{c} = mmap_mark(chainSysArrivals{c}, markProb);
                                chainArrivalAtNode{c} = mmap_normalize(chainArrivalAtNode{c});
                                chainArrivalAtNode{c} = mmap_scale(chainArrivalAtNode{c}, 1./rates{ist,c}, 0); % non-iterative approximation
                                if c == 1
                                    aggrArrivalAtNode = mmap_super_safe({chainArrivalAtNode{c}, mmap_exponential(0,1)}, config.space_max, 'default');
                                    aggrArrivalAtNode = {aggrArrivalAtNode{1} aggrArrivalAtNode{2} aggrArrivalAtNode{2}};
                                    lc = map_lambda(chainArrivalAtNode{c});
                                    if lc>0
                                        aggrArrivalAtNode = mmap_scale(aggrArrivalAtNode, 1/lc, 0); % non-iterative approximation
                                    end
                                else
                                    aggrArrivalAtNode = mmap_super_safe({aggrArrivalAtNode, chainArrivalAtNode{c}}, config.space_max, 'default');
                                end
                            end
                            Qret = cell(1,K);
                            if (strcmp(sn.schedid(ist),SchedStrategy.ID_HOL) && any(sn.classprio ~= sn.classprio(1))) % if priorities are not identical
                                [uK,iK] = unique(sn.classprio);
                                if length(uK) == length(sn.classprio) % if all priorities are different
                                    [Qret{iK'}] = MMAPPH1NPPR({aggrArrivalAtNode{[1;2+iK]}}, {pie{ist}{iK}}, {D0{ist,iK}}, 'ncMoms', 1);
                                else
                                    line_error(mfilename,'Solver MAM requires either identical priorities or all distinct priorities');
                                end
                            else
                                aggrUtil = nansum(mmap_lambda(aggrArrivalAtNode)./(GlobalConstants.FineTol+sn.rates(ist,1:K)));   %% to debug
                                aggrLambda = mmap_lambda(aggrArrivalAtNode);
                                if aggrUtil < 1-GlobalConstants.FineTol
                                    if any(isinf(N))
                                        [Qret{1:K}] = MMAPPH1FCFS({aggrArrivalAtNode{[1,3:end]}}, {pie{ist}{:}}, {D0{ist,:}}, 'ncMoms', 1);
                                    else % all closed classes
                                        maxLevel = max(N(isfinite(N)))+1;
                                        D = {aggrArrivalAtNode{[1,3:end]}};
                                        pdistr = cell(1,K);
                                        if map_lambda(D)< GlobalConstants.FineTol
                                            for k=1:K
                                                pdistr{k} = [1-GlobalConstants.FineTol, GlobalConstants.FineTol];
                                                Qret{k} = GlobalConstants.FineTol / sn.rates(ist);
                                            end
                                        else
                                            if sn.isfunction(ist)
                                                [Qret{1:K}] = qbd_setupdelayoff(ist, sn, aggrLambda);                                               
                                            else
                                                [pdistr{1:K}] = MMAPPH1FCFS(D, {pie{ist}{:}}, {D0{ist,:}}, 'ncDistr', maxLevel);
                                            
                                            for k=1:K
                                                pdistr{k} = abs(pdistr{k}(1:(N(k)+1)));
                                                pdistr{k}(end) = abs(1-sum(pdistr{k}(1:end-1)));
                                                pdistr{k} = pdistr{k} / sum(pdistr{k});
                                                Qret{k} = max(0,min(N(k),(0:N(k))*pdistr{k}(1:(N(k)+1))'));
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
                            QN(ist,:) = cell2mat(Qret);
                            for k=1:K
                                c = find(sn.chains(:,k),1);
                                TN(ist,k) = rates{ist,c}(k);
                                UN(ist,k) = TN(ist,k) * S(ist,k);
                                % add number of jobs at the surrogate delay server
                                QN(ist,k) = Qret{k};
                                QN(ist,k) = QN(ist,k) + TN(ist,k)*(S(ist,k)*sn.nservers(ist)) * (sn.nservers(ist)-1)/sn.nservers(ist);
                                RN(ist,k) = QN(ist,k) ./ TN(ist,k);
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
end

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
end
