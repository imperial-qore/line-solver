function [Q,U,R,T,C,X,lG] = solver_mva(sn,options)
% [Q,U,R,T,C,X,LG] = SOLVER_MVA(SN,OPTIONS)

% Copyright (c) 2012-2024, Imperial College London
% All rights reserved.

% aggregate chains

if nargin < 2
    options = SolverMVA.defaultOptions;
end

[Lchain,STchain,Vchain,alpha,Nchain,~,refstatchain] = snGetDemandsChain(sn);

nservers = sn.nservers;
schedid = sn.schedid;
M = sn.nstations;
K = sn.nchains;

infSET =[]; % set of infinite server stations
qSET =[]; % set of other product-form stations
for i=1:M
    switch schedid(i)
        case SchedStrategy.ID_EXT
            % no-op
        case SchedStrategy.ID_INF
            infSET(1,end+1) = i;
        case {SchedStrategy.ID_PS, SchedStrategy.ID_LCFSPR}
            qSET(1,end+1) = i; 
        case {SchedStrategy.ID_FCFS, SchedStrategy.ID_SIRO}
            if range(sn.rates(i,:))==0
                qSET(end+1) = i;
            end
        otherwise
            line_error(mfilename, sprintf('Unsupported exact MVA analysis for %s scheduling.',SchedStrategy.toFeature(SchedStrategy.fromId(schedid(i)))));
    end
end

Uchain = zeros(M,K); Tchain = zeros(M,K); C = zeros(1,K); Wchain = zeros(M,K); Qchain = zeros(M,K);

lambda = zeros(1,K);
ocl = find(isinf(Nchain));
if any(isinf(Nchain))
    for r=ocl % open classes
        lambda(r) = 1 ./ STchain(refstatchain(r),r);
        Qchain(refstatchain(r),r) = Inf;
    end
end
rset = setdiff(1:K,find(Nchain==0));

[Xchain,Qpf,Uchain,~,lG] = pfqn_mvams(lambda,STchain(qSET,:).*Vchain(qSET,:),Nchain,STchain(infSET,:).*Vchain(infSET,:),ones(length(qSET),1),nservers(qSET));
Qchain(qSET,:) = Qpf;
Qchain(infSET,:) = repmat(Xchain,numel(infSET),1) .* STchain(infSET,:) .* Vchain(infSET,:);

ccl = find(isfinite(Nchain));
for r=rset
    for k=infSET(:)'
        Wchain(k,r) = STchain(k,r);
    end
    for k=qSET(:)'
        if isinf(nservers(k)) % infinite server
            Wchain(k,r) = STchain(k,r);
        else
            if Vchain(k,r) == 0 || Xchain(r) == 0
                Wchain(k,r) = 0;
            else
                Wchain(k,r) = Qchain(k,r) / (Xchain(r) * Vchain(k,r));
            end
        end
    end
end

for r=rset
    if sum(Wchain(:,r)) == 0
        Xchain(r) = 0;
    else
        if isinf(Nchain(r))
            C(r) = Vchain(:,r)'*Wchain(:,r);
            % X(r) remains constant
        elseif Nchain(r)==0
            Xchain(r) = 0;
            C(r) = 0;
        else
            C(r) = Vchain(:,r)'*Wchain(:,r);
            Xchain(r) = Nchain(r) / C(r);
        end
    end
    
    for k=1:M
        Qchain(k,r) = Xchain(r) * Vchain(k,r) * Wchain(k,r);
        Tchain(k,r) = Xchain(r) * Vchain(k,r);
    end
end

for k=1:M
    for r=rset
        if isinf(nservers(k)) % infinite server
            Uchain(k,r) = Vchain(k,r)*STchain(k,r)*Xchain(r);
        else
            Uchain(k,r) = Vchain(k,r)*STchain(k,r)*Xchain(r)/nservers(k);
        end
    end
end

for k=1:M
    for r=1:K
        if Vchain(k,r)*STchain(k,r) > options.tol
            switch schedid(k)
                case {SchedStrategy.ID_FCFS,SchedStrategy.ID_PS}
                    if sum(Uchain(k,:))>1+options.tol
                        Uchain(k,r) = min(1,sum(Uchain(k,:))) * Vchain(k,r)*STchain(k,r)*Xchain(r) / ((Vchain(k,:).*STchain(k,:))*Xchain(:));
                    end
            end
        end
    end
end

Vsink = cellsum(sn.nodevisits);
Vsink = Vsink(find(sn.nodetype==NodeType.ID_SINK),:);
for r=find(isinf(Nchain)) % open classes
    Xchain(r) = Vsink(r) ./ STchain(refstatchain(r),r);
end

Rchain = Qchain./Tchain;
Xchain(~isfinite(Xchain))=0;
Uchain(~isfinite(Uchain))=0;
Qchain(~isfinite(Qchain))=0;
Rchain(~isfinite(Rchain))=0;

Xchain(Nchain==0)=0;
Uchain(:,Nchain==0)=0;
Qchain(:,Nchain==0)=0;
Rchain(:,Nchain==0)=0;
Tchain(:,Nchain==0)=0;
Wchain(:,Nchain==0)=0;

[Q,U,R,T,C,X] = snDeaggregateChainResults(sn, Lchain, [], STchain, Vchain, alpha, [], [], Rchain, Tchain, [], Xchain);
end
