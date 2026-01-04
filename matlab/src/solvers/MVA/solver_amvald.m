function [Q,U,R,T,C,X,lG,totiter] = solver_amvald(sn, Lchain,STchain,Vchain,alpha,Nchain,SCVchain,refstatchain, options)
totiter = 0;

M = sn.nstations;
K = sn.nchains;
Nt = sum(Nchain(isfinite(Nchain)));
delta  = (Nt - 1) / Nt;
deltaclass = (Nchain - 1) ./ Nchain;
deltaclass(isinf(Nchain)) = 1;
tol = options.iter_tol;
nservers = sn.nservers;
schedparam = sn.schedparam;
lldscaling = sn.lldscaling;
cdscaling = sn.cdscaling;
ljdscaling = sn.ljdscaling;
ljdcutoffs = sn.ljdcutoffs;
ljcdscaling = sn.ljcdscaling;
ljcdcutoffs = sn.ljcdcutoffs;
sched = sn.sched;
classprio = sn.classprio;

Uchain = zeros(M,K);
Tchain = zeros(M,K);
Cchain_s = zeros(1,K);

%% initialize Q,X, U
Qchain = options.init_sol;
if isempty(Qchain)
    % balanced initialization
    Qchain = ones(M,K);
    Qchain = Qchain ./ repmat(sum(Qchain,1),size(Qchain,1),1) .* repmat(Nchain,size(Qchain,1),1);
    Qchain(isinf(Qchain))=0; % open classes
    for r=find(isinf(Nchain)) % open classes
        Qchain(refstatchain(r),r)=0;
    end
end

nnzclasses = find(Nchain>0);
Xchain = 1./sum(STchain,1);
for r=find(isinf(Nchain)) % open classes
    Xchain(r) = 1 ./ STchain(refstatchain(r),r);
end

for k=1:M
    for r=nnzclasses
        if isinf(nservers(k)) % infinite server
            Uchain(k,r) = Vchain(k,r)*STchain(k,r)*Xchain(r);
        else
            Uchain(k,r) = Vchain(k,r)*STchain(k,r)*Xchain(r)/nservers(k);
        end
    end
end

switch options.method
    case {'lin','qdlin'}
        gamma = zeros(K,M,K); % class-based customer fraction corrections
        tau = zeros(K,K); % throughput difference
    otherwise
        gamma = zeros(K,M); % total customer fraction corrections
        tau = zeros(K,K); % throughput difference
end

%% main loop
omicron = 0.5; % under-relaxation parameter
outer_iter = 0;
while (outer_iter < 2 || max(max(abs(Qchain-QchainOuter_1))) > tol) && outer_iter < sqrt(options.iter_max)
    outer_iter = outer_iter + 1;
    QchainOuter_1 = Qchain;
    XchainOuter_1 = Xchain;
    UchainOuter_1 = Uchain;

    if isfinite(Nt) && Nt>0
        switch options.method
            %case {'aql','qdaql'}
            %line_error(mfilename,'AQL is currently disabled in SolverMVA, please use the SolverJMT implementation (method jmva.aql).');
            case {'lin','qdlin'}
                % iteration at population N-1_s
                for s=1:K
                    if isfinite(Nchain(s)) % don't recur on open classes
                        iter_s = 0;
                        Nchain_s = oner(Nchain,s);
                        Qchain_s = Qchain * (Nt-1)/Nt;
                        Xchain_s = Xchain * (Nt-1)/Nt;
                        Uchain_s = Uchain * (Nt-1)/Nt;
                        while (iter_s < 2 || max(max(abs(Qchain_s-Qchain_s_1))) > tol) && iter_s <= sqrt(options.iter_max)
                            iter_s = iter_s + 1;

                            Qchain_s_1 = Qchain_s;
                            Xchain_s_1 = Xchain_s;
                            Uchain_s_1 = Uchain_s;

                            [Wchain_s, STeff_s] = solver_amvald_forward(M, K, nservers, schedparam, lldscaling, cdscaling, ljdscaling, ljdcutoffs, ljcdscaling, ljcdcutoffs, sched, classprio, gamma, tau, Qchain_s_1, Xchain_s_1, Uchain_s_1, STchain, Vchain, Nchain_s, SCVchain, options);
                            totiter = totiter + 1;
                            if totiter > options.iter_max
                                break
                            end

                            %% update other metrics
                            for r=nnzclasses
                                if sum(Wchain_s(:,r)) == 0
                                    Xchain_s(r) = 0;
                                else
                                    if isinf(Nchain_s(r))
                                        Cchain_s(r) = Vchain(:,r)' * Wchain_s(:,r);
                                        % X(r) remains constant
                                    elseif Nchain(r)==0
                                        Xchain_s(r) = 0;
                                        Cchain_s(r) = 0;
                                    else
                                        Cchain_s(r) = Vchain(:,r)' * Wchain_s(:,r);
                                        Xchain_s(r) = omicron * Nchain_s(r) / Cchain_s(r) + (1-omicron) * Xchain_s_1(r);
                                    end
                                end
                                for k=1:M
                                    Rchain_s(k,r) = Vchain(k,r) * Wchain_s(k,r);
                                    Qchain_s(k,r) = omicron * Xchain_s(r) * Vchain(k,r) * Wchain_s(k,r) + (1-omicron) * Qchain_s_1(k,r);
                                    Tchain_s(k,r) = Xchain_s(r) * Vchain(k,r);
                                    Uchain_s(k,r) = omicron * Vchain(k,r) * STeff_s(k,r) * Xchain_s(r) + (1-omicron) * Uchain_s_1(k,r);
                                end
                            end
                        end

                        switch options.method
                            case {'lin'}
                                for k=1:M
                                    for r=nnzclasses
                                        if ~isinf(Nchain(r)) && Nchain_s(r)>0
                                            gamma(s,k,r) = Qchain_s_1(k,r)./Nchain_s(r) - QchainOuter_1(k,r)./Nchain(r);
                                        end
                                    end
                                end
                            otherwise
                                for k=1:M
                                    gamma(s,k) = sum(Qchain_s_1(k,:),2)/(Nt-1) - sum(QchainOuter_1(k,:),2)/Nt;
                                end
                        end

                        for r=nnzclasses
                            tau(s,r) = Xchain_s_1(r) - XchainOuter_1(r); % save throughput for priority AMVA
                        end
                    end
                end
        end
    end

    iter = 0;
    % iteration at population N
    while (iter < 2 || max(max(abs(Qchain-Qchain_1))) > tol) && iter <= sqrt(options.iter_max)
        iter = iter + 1;

        Qchain_1 = Qchain;
        Xchain_1 = Xchain;
        Uchain_1 = Uchain;

        [Wchain, STeff] = solver_amvald_forward(M, K, nservers, schedparam, lldscaling, cdscaling, ljdscaling, ljdcutoffs, ljcdscaling, ljcdcutoffs, sched, classprio, gamma, tau, Qchain_1, Xchain_1, Uchain_1, STchain, Vchain, Nchain, SCVchain, options);
        totiter = totiter + 1;
        if totiter > options.iter_max
            break
        end


        %% update other metrics
        for r=nnzclasses
            if sum(Wchain(:,r)) == 0
                Xchain(r) = 0;
            else
                if isinf(Nchain(r))
                    Cchain_s(r) = Vchain(:,r)'*Wchain(:,r);
                    % X(r) remains constant
                elseif Nchain(r)==0
                    Xchain(r) = 0;
                    Cchain_s(r) = 0;
                else
                    Cchain_s(r) = Vchain(:,r)'*Wchain(:,r);
                    Xchain(r) = omicron * Nchain(r) / Cchain_s(r) + (1-omicron) *Xchain_1(r);
                end
            end
            for k=1:M
                Rchain(k,r) = Vchain(k,r) * Wchain(k,r);
                Qchain(k,r) = omicron * Xchain(r) * Vchain(k,r) * Wchain(k,r) + (1-omicron) * Qchain_1(k,r);
                Tchain(k,r) = Xchain(r) * Vchain(k,r);
                Uchain(k,r) = omicron * Vchain(k,r) * STeff(k,r) * Xchain(r) + (1-omicron) * Uchain_1(k,r);
            end
        end
    end
end


% the next block is a coarse approximation for LD and CD, would need
% cdterm and qterm in it but these are hidden within the iteration calls
for k=1:M
    for r=1:K
        if Vchain(k,r) * STeff(k,r) >0
            switch sn.sched(k)
                case {SchedStrategy.FCFS, SchedStrategy.SIRO, SchedStrategy.PS, SchedStrategy.LCFSPR, SchedStrategy.DPS, SchedStrategy.HOL}
                    if sum(Uchain(k,:))>1
                        Uchain(k,r) = min(1,sum(Uchain(k,:))) * Vchain(k,r) * STeff(k,r) * Xchain(r) / ((Vchain(k,:) .* STeff(k,:)) * Xchain(:));
                    end
            end
        end
    end
end

Rchain = Qchain./Tchain;
Xchain(~isfinite(Xchain))=0;
Uchain(~isfinite(Uchain))=0;
%Qchain(~isfinite(Qchain))=0;
Rchain(~isfinite(Rchain))=0;

Xchain(Nchain==0)=0;
Uchain(:,Nchain==0)=0;
%Qchain(:,Nchain==0)=0;
Rchain(:,Nchain==0)=0;
Tchain(:,Nchain==0)=0;

if isempty(sn.lldscaling) && isempty(sn.cdscaling) && ~sn_has_joint_dependence(sn)
    [Q,U,R,T,C,X] = sn_deaggregate_chain_results(sn, Lchain, [], STchain, Vchain, alpha, [], [], Rchain, Tchain, [], Xchain);
else
    [Q,U,R,T,C,X] = sn_deaggregate_chain_results(sn, Lchain, [], STchain, Vchain, alpha, [], Uchain, Rchain, Tchain, [], Xchain);
end

% estimate normalizing constant for closed classes
ccl = isfinite(Nchain);
Nclosed = Nchain(ccl);
Xclosed = Xchain(ccl);
lG = - Nclosed(Xclosed>options.tol) * log(Xclosed(Xclosed>options.tol))'; % asymptotic approximation
end