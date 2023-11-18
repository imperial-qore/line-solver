function [Q,U,R,T,C,X,lG,totiter] = solver_qna_mam_closed(sn, options)

config = options.config;
config.space_max = 16;

K = sn.nclasses;
rt = sn.rt;
S = 1./sn.rates;
scv = sn.scv; scv(isnan(scv))=0;
N = sn.njobs';
PH = sn.proc;
I = sn.nnodes;
M = sn.nstations;
C = sn.nchains;
V = cellsum(sn.visits);
Q = zeros(M,K);


U = zeros(M,K);
R = zeros(M,K);
T = zeros(M,K);
X = zeros(1,K);

lambda = zeros(1,C);


pie = {};
D0 = {};
% get service process 
for ist=1:M
    switch sn.schedid(ist)
        case {SchedStrategy.ID_FCFS}
            for k=1:K
                pie{ist}{k} = map_pie(PH{ist}{k});
                D0{ist,k} = PH{ist}{k}{1};
            end
    end
end 

QNc = sn.njobs;
lambda_lb = zeros(1,K);
lambda_ub = zeros(1,K);
for k=1:k
    lambda_ub(k) = min(sn.rates(find(sn.nservers<Inf),k));
end
QN = zeros(1,K);
it_out = 0;

%% main iteration
while max(abs(QN-QNc))>options.iter_tol && it_out < 200
    it_out = it_out+1;
    

    if it_out ~=1
        [~,idx] = max(abs(QN-QNc).*(lambda_ub-lambda_lb)./QNc);
        
            if QN(idx)<QNc(idx)
                lambda_lb(idx)=lambda(idx);
                
            elseif QN(idx)>QNc(idx)
                lambda_ub(idx)=lambda(idx);
               
            end
            lambda=(lambda_ub+lambda_lb)/2;
    else
        lambda = lambda_ub;
    end 
    it =0;
    Q = zeros(M,K);
    

    U = zeros(M,K);
    R = zeros(M,K);
    T = zeros(M,K);
    X = zeros(1,K);
  
    a1 = zeros(M,K);
    a2 = zeros(M,K);
    a1_1 = a1+Inf;
    a2_1 = a2+Inf;
    
    d2 = zeros(M,1);
    f2 = zeros(M*K,M*K); 
    for i=1:M
        for j=1:M
            if sn.nodetype(sn.stationToNode(j)) ~= NodeType.Source
                for r=1:K
                    for s=1:K
                        if rt((i-1)*K+r, (j-1)*K+s)>0
                                f2((i-1)*K+r, (j-1)*K+s) = 1; % C^2ij,r
                        end
                     end
                 end
             end
         end
    end   

    while (max(max(abs(a1_1-a1))) > options.iter_tol || max(max(abs(a2_1-a2))) > options.iter_tol )&& it <= options.iter_max %#ok<max>
        it = it + 1;
        a1_1 = a1;
        a2_1 = a2;
        
    
        % update throughputs at all stations
        if it==1
            for c=1:C
                inchain = sn.inchain{c};
                for m=1:M
                    T(m,inchain) = V(m,inchain) .* lambda(c);
                end
            end
        end
    
        % superposition
        for i=1:M
            a1(i,:) = 0;
            a2(i,:) = 0;
            lambda_i = sum(T(i,:));        
            for j=1:M
                for r=1:K
                    for s=1:K
                        a1(i,r) = a1(i,r) + T(j,s)*rt((j-1)*K+s, (i-1)*K+r);
                        a2(i,r) = a2(i,r) + (1/lambda_i) * f2((j-1)*K+s, (i-1)*K+r)*T(j,s)*rt((j-1)*K+s, (i-1)*K+r);
                    end
                end
            end
        end
    
        % update flow trhough queueing station
        for ind=1:I
            if sn.isstation(ind)
                ist = sn.nodeToStation(ind);
                switch sn.nodetype(ind)
                    case NodeType.Join
                        % no-op
                  
                    otherwise
                        switch sn.schedid(ist)
                            case SchedStrategy.ID_INF
                                for i=1:M
                                    for r=1:K
                                        for s=1:K
                                            d2(ist,s) = a2(ist,s);
                                        end
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
                        
                            case {SchedStrategy.ID_FCFS}
                             
                                
                                mu_ist = sn.rates(ist,1:K);
                                mu_ist(isnan(mu_ist))=0;                            
                                rho_ist_class = a1(ist,1:K)./(GlobalConstants.FineTol+sn.rates(ist,1:K));
                                rho_ist_class(isnan(rho_ist_class))=0;
                                lambda_ist = sum(a1(ist,:));
                                mi = sn.nservers(ist);
                                rho_ist = sum(rho_ist_class) / mi;
                                if rho_ist < 1-options.tol
                                    
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
                                else
                                    for k=1:K
                                        Q(ist,k) = sn.njobs(k);
                                    end
                                    d2(ist) = 1;
                                end
                                for k=1:K
                                    T(ist,k) = a1(ist,k);
                                    U(ist,k) = T(ist,k) * S(ist,k) /sn.nservers(ist);
                                    R(ist,k) = Q(ist,k) ./ T(ist,k);
                                end

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
        for i=1:M
            for j=1:M
                if sn.nodetype(sn.stationToNode(j)) ~= NodeType.Source
                    for r=1:K
                        for s=1:K
                            if rt((i-1)*K+r, (j-1)*K+s)>0
                                f2((i-1)*K+r, (j-1)*K+s) = 1 + rt((i-1)*K+r, (j-1)*K+s) * (d2(i)-1); 
                            end
                        end
                    end
                end
            end
        end   
    end

    for ind=1:I
        if sn.isstation(ind)
           ist = sn.nodeToStation(ind);
           switch sn.schedid(ist)
           case {SchedStrategy.ID_FCFS}
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
                                arri_class = APH.fitMeanAndSCV(1/a1(ist,k),a2(ist,k)).getRepres; % MMAP repres of arrival process for class k at node ist
                                %arri_class = map_exponential(1/a1(ist,k));
                                arri_class = {arri_class{1},arri_class{2},arri_class{2}};
                            end
                            if k==1
                                arri_node = arri_class; 
                            else
                        %arri_node = mmap_super(arri_node,arri_class,  'default');
                       
                            arri_node = mmap_super_safe({arri_node,arri_class}, config.space_max, 'default'); % combine arrival process from different class 
                            end
                   
                        end
                        Qret = cell(1,K);
                        maxLevel = max(N(isfinite(N)))+1;
                        D = {arri_node{[1,3:end]}};
                        pdistr = cell(1,K);
                        if map_lambda(D)< GlobalConstants.FineTol
                            for k=1:K
                                pdistr{k} = [1-GlobalConstants.FineTol, GlobalConstants.FineTol];
                                Qret{k} = GlobalConstants.FineTol / sn.rates(ist);
                            end
                        else
                            [pdistr{1:K}] = MMAPPH1FCFS(D, {pie{ist}{:}}, {D0{ist,:}}, 'ncDistr', maxLevel);
                            for k=1:K
                               pdistr{k} = abs(pdistr{k}(1:(N(k)+1)));
                               pdistr{k}(end) = abs(1-sum(pdistr{k}(1:end-1)));
                               pdistr{k} = pdistr{k} / sum(pdistr{k});
                               Qret{k} = max(0,min(N(k),(0:N(k))*pdistr{k}(1:(N(k)+1))'));
                            end
                        end
                        Q(ist,:) = cell2mat(Qret);
                    else
                        for k=1:K
                            Q(ist,k) = sn.njobs(k);
                        end
                    end
        
                    for k=1:K                         
                        R(ist,k) = Q(ist,k) ./ T(ist,k);
                    end
            end
        end
    end
    QN = sum(Q,1);
end
%it_out
C = sum(R,1);
Q = abs(Q);
Q(isnan(Q))=0;
U(isnan(U))=0;
R(isnan(R))=0;
C(isnan(C))=0;
X(isnan(X))=0;
lG = 0;
totiter = it_out;
end