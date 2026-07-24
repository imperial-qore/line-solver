function [Q,U,R,T,C,X,lG,totiter] = solver_mna_open(sn, options)

config = options.config;
config.space_max = 32;
if ~isfield(config, 'dep_scv')
    config.dep_scv = 'qna'; % 'qna' for Whitt formula, 'etaqa' for QBD joint moments
end

K = sn.nclasses;
rt = sn.rt;
S = 1./sn.rates;
scv = sn.scv; scv(isnan(scv))=0;

PH = sn.proc;
I = sn.nnodes;
M = sn.nstations;
C = sn.nchains;
V = cellsum(sn.visits);
Q = zeros(M,K);
%QN_1 = Q+Inf;

U = zeros(M,K);
R = zeros(M,K);
T = zeros(M,K);
X = zeros(1,K);

lambda = zeros(1,C);

it = 0;
pie = {};
D0 = {};

% get service process 
for ist=1:M
    switch sn.sched(ist)
        case {SchedStrategy.FCFS, SchedStrategy.INF,SchedStrategy.PS}
            for k=1:K
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

a1 = zeros(M,K);
a2 = zeros(M,K);
mubar = [];
c2 = [];
d2 = zeros(M,1);
f2 = zeros(M*K,M*K); 
for ist=1:M
    for jst=1:M
        if sn.nodetype(sn.stationToNode(jst)) ~= NodeType.Source
            for r=1:K
                for s=1:K
                    if rt((ist-1)*K+r, (jst-1)*K+s)>0
                            f2((ist-1)*K+r, (jst-1)*K+s) = 1; % C^2ij,r
                    end
                 end
             end
         end
     end
end   
lambdas_inchain = cell(1,C);
scvs_inchain = cell(1,C);
d2c = [];
for c=1:C
    inchain = sn.inchain{c};
    sourceIdx = sn.refstat(inchain(1));
    lambdas_inchain{c} = sn.rates(sourceIdx,inchain);
    scvs_inchain{c} = scv(sourceIdx,inchain);
    lambda(c) = sum(lambdas_inchain{c}(isfinite(lambdas_inchain{c})));
    d2c(c) = da_traffic_superpos(lambdas_inchain{c},scvs_inchain{c});
    T(sourceIdx,inchain') = lambdas_inchain{c};
   
end
d2(sourceIdx)=d2c(sourceIdx,:)*lambda'/sum(lambda);

%% main iteration
% flow fixed point on the per-station arrival rates a1 and SCVs a2, driven
% by the generic DA successive-substitution driver
fpopts = options;
fpopts.iter_max = options.iter_max + 1; % legacy while-loop executed one extra sweep at the cap
fpopts.config.da_nanstop = true; % legacy while-loop exited on NaN convergence measure
[~, it] = da_fpi(@mna_sweep, [a1(:); a2(:)], fpopts);

    function [xnew, xref] = mna_sweep(~, itnum)
    xref = [a1(:); a2(:)];
    % update throughputs at all stations
    if itnum==1
        for c=1:C
            inchain = sn.inchain{c};
            for m=1:M
                T(m,inchain) = V(m,inchain) .* lambda(c);
            end
        end
    end

    % superposition
    for ist=1:M
        a1(ist,:) = 0;
        a2(ist,:) = 0;
        lambda_i = sum(T(ist,:));        
        for jst=1:M
            for r=1:K
                for s=1:K
                    a1(ist,r) = a1(ist,r) + T(jst,s)*rt((jst-1)*K+s, (ist-1)*K+r);
                    a2(ist,r) = a2(ist,r) + (1/lambda_i) * f2((jst-1)*K+s, (ist-1)*K+r)*T(jst,s)*rt((jst-1)*K+s, (ist-1)*K+r);
                end
            end
        end
    end

    % update flow trhough queueing station
    for ind=1:I
        if sn.isstation(ind)
            ist = sn.nodeToStation(ind);
            
                    switch sn.sched(ist)
                        case SchedStrategy.INF
                                for r=1:K
                                    for s=1:K
                                        d2(ist,s) = a2(ist,s);
                                    end
                            end
                            for c=1:C
                                inchain = sn.inchain{c};
                                for k=inchain
                                    T(ist,k) = a1(ist,k);
                                    U(ist,k) = S(ist,k)*T(ist,k);
                                    Q(ist,k) = T(ist,k).*S(ist,k)*V(ist,k);
                                    R(ist,k) = Q(ist,k)/T(ist,k);
                                end
                            end
                        case SchedStrategy.PS
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
                        case {SchedStrategy.FCFS}
                            mu_ist = sn.rates(ist,1:K);
                            mu_ist(isnan(mu_ist))=0;                            
                            rho_ist_class = a1(ist,1:K)./(GlobalConstants.FineTol+sn.rates(ist,1:K));
                            rho_ist_class(isnan(rho_ist_class))=0;
                            lambda_ist = sum(a1(ist,:));
                            mi = sn.nservers(ist);
                            rho_ist = sum(rho_ist_class) / mi;
                            if rho_ist < 1-options.tol
                                if strcmp(config.dep_scv, 'etaqa') && mi == 1
                                    % ETAQA-based departure SCV via QBD joint moments
                                    try
                                        % fit aggregate arrival MAP from parametric info
                                        arri_agg = APH.fitMeanAndSCV(1/lambda_ist, sum(a2(ist,:))).getProcess;
                                        serv_agg = APH.fitMeanAndSCV(1/sum(mu_ist(mu_ist>0).*a1(ist,mu_ist>0))/lambda_ist, ...
                                            sum(a1(ist,:).*scv(ist,:))/lambda_ist).getProcess;
                                        JM = qbd_depproc_jointmom(arri_agg, serv_agg, [1,0; 2,0]);
                                        E1 = JM(1); E2 = JM(2);
                                        d2(ist) = (E2 - E1^2) / E1^2;
                                    catch
                                        % fall back to QNA formula on failure
                                        for k=1:K
                                            mubar(ist) = lambda_ist ./ rho_ist;
                                            c2(ist) = -1;
                                            for r=1:K
                                                if mu_ist(r)>0
                                                    c2(ist) = c2(ist) + a1(ist,r)/lambda_ist * (mubar(ist)/mi/mu_ist(r))^2 * (scv(ist,r)+1 );
                                                end
                                            end
                                        end
                                        d2(ist) = 1 + rho_ist^2*(c2(ist)-1)/sqrt(mi) + (1 - rho_ist^2) *(sum(a2(ist,:))-1);
                                    end
                                else
                                    % QNA Whitt formula
                                    for k=1:K
                                        mubar(ist) = lambda_ist ./ rho_ist;
                                        c2(ist) = -1;
                                        for r=1:K
                                            if mu_ist(r)>0
                                                c2(ist) = c2(ist) + a1(ist,r)/lambda_ist * (mubar(ist)/mi/mu_ist(r))^2 * (scv(ist,r)+1 );
                                            end
                                        end
                                    end
                                    d2(ist) = 1 + rho_ist^2*(c2(ist)-1)/sqrt(mi) + (1 - rho_ist^2) *(sum(a2(ist,:))-1);
                                end
                            else
                                for k=1:K
                                    Q(ist,k) = sn.njobs(k);
                                end
                                d2(ist) = 1;
                            end
                            for k=1:K
                                T(ist,k) = a1(ist,k);
                                U(ist,k) = T(ist,k) * S(ist,k) /sn.nservers(ist);
                               
                            end
                    end
            
        else % not a station
            switch sn.nodetype(ind)
                case NodeType.Fork
                    line_error(mfilename,'Fork nodes not supported yet by QNA solver.');
            end
        end
    end    


    % splitting - update flow scvs
    for ist=1:M
        for jst=1:M
            if sn.nodetype(sn.stationToNode(jst)) ~= NodeType.Source
                for r=1:K
                    for s=1:K
                        if rt((ist-1)*K+r, (jst-1)*K+s)>0
                            f2((ist-1)*K+r, (jst-1)*K+s) = 1 + rt((ist-1)*K+r, (jst-1)*K+s) * (d2(ist)-1); 
                        end
                    end
                end
            end
        end
    end   
    xnew = [a1(:); a2(:)];
    end
                

for ind=1:I
        if sn.isstation(ind)
           ist = sn.nodeToStation(ind);
           switch sn.sched(ist)
            case {SchedStrategy.FCFS}
                mu_ist = sn.rates(ist,1:K);
                mu_ist(isnan(mu_ist))=0;                            
                rho_ist_class = a1(ist,1:K)./(GlobalConstants.FineTol+sn.rates(ist,1:K));
                rho_ist_class(isnan(rho_ist_class))=0;
                lambda_ist = sum(a1(ist,:));
                mi = sn.nservers(ist);
                rho_ist = sum(rho_ist_class) / mi;
                if rho_ist < 1-options.tol
                    for k=1:K
                        if a1(ist,k)==0
                            arri_class = map_exponential(Inf);
                         else
                            arri_class = APH.fitMeanAndSCV(1/a1(ist,k),a2(ist,k)).getProcess; % MMAP repres of arrival process for class k at node ist
                            %arri_class = Erlang.fit(1/a1(ist,k),a2(ist,k)).getProcess;
                            %arri_class = map_exponential(1/a1(ist,k));
                            arri_class = {arri_class{1},arri_class{2},arri_class{2}};
                         end
                         if k==1
                            arri_node = arri_class;
                         else
                            arri_node = mmap_super(arri_node,arri_class,  'default');

                            %arri_node = mmap_super_safe({arri_node,arri_class}, config.space_max, 'default'); % combine arrival process from different class
                         end
                    end
                isFiniteCap = isfinite(sn.cap(ist));
                if isFiniteCap
                    capK = sn.cap(ist);
                    [isMmck, muMmck] = mam_detect_mmck(sn, ist, K, arri_node);
                    if isMmck
                        aggrLambda_ist = sum(a1(ist,1:K), 'omitnan');
                        exactRes = qsys_mmck(aggrLambda_ist, muMmck, sn.nservers(ist), capK);
                        meanQ_fc = exactRes.meanQueueLength;
                        lossProb_fc = exactRes.lossProbability;
                    else
                        [meanQ_fc, lossProb_fc, ~] = mam_truncate_renorm( ...
                            {arri_node{[1,3:end]}}, {pie{ist}{:}}, {D0{ist,:}}, capK);
                    end
                    lambdaInflow = a1(ist,1:K);
                    lambdaInflow(isnan(lambdaInflow)) = 0;
                    T_eff = lambdaInflow * (1 - lossProb_fc);
                    sumT = sum(T_eff);
                    if sumT > 0
                        Savg_eff = sum(T_eff .* S(ist,1:K), 'omitnan') / sumT;
                        Wq = max(0, meanQ_fc / sumT - Savg_eff);
                    else
                        Wq = 0;
                    end
                    for k=1:K
                        T(ist,k) = T_eff(k);
                        U(ist,k) = T(ist,k) * S(ist,k) / sn.nservers(ist);
                        if T(ist,k) > 0
                            R(ist,k) = Wq + S(ist,k);
                            Q(ist,k) = T(ist,k) * R(ist,k);
                        else
                            R(ist,k) = 0;
                            Q(ist,k) = 0;
                        end
                    end
                else
                    Qret = cell(1,K);
                    [Qret{1:K}] = MMAPPH1FCFS({arri_node{[1,3:end]}}, {pie{ist}{:}}, {D0{ist,:}}, 'ncMoms', 1);
                    Q(ist,:) = cell2mat(Qret);
                    for k=1:K
                        R(ist,k) = Q(ist,k) ./ T(ist,k);
                    end
                end
                else
                    for k=1:K
                        Q(ist,k) = sn.njobs(k);
                        R(ist,k) = Q(ist,k) ./ T(ist,k);
                    end
                end
            end
            
        end
end

C = sum(R,1);
Q = abs(Q);
Q(isnan(Q))=0;
U(isnan(U))=0;
R(isnan(R))=0;
C(isnan(C))=0;
X(isnan(X))=0;
lG = 0;
totiter = it;
end
