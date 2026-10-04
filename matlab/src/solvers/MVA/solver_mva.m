function [Q,U,R,T,C,X,lG] = solver_mva(sn,options)
% [Q,U,R,T,C,X,LG] = SOLVER_MVA(SN,OPTIONS)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 2
    options = SolverMVA.defaultOptions;
end

[Lchain,STchain,Vchain,alpha,Nchain,~,refstatchain] = sn_get_demands_chain(sn);

nservers = sn.nservers;
sched = sn.sched;
M = sn.nstations;
K = sn.nchains;

% METHOD 'mva' IS THE DELIBERATE APPROXIMATION: SolverMVA.mvaDispatch already
% warns that the exact recursion is being run outside its hypotheses and
% promises an answer, so erroring here would contradict its own message. Only
% an implicit or 'exact' request is refused.
if ~sn_has_product_form(sn) && ~(isfield(options,'method') && strcmp(options.method,'mva'))
    line_error(mfilename, 'Unsupported exact MVA analysis, the model does not have a product form');
end

% Check for special LCFS + LCFS-PR 2-station network. The pairing, the station
% and server counts, the closed population and the self-loop rule are judged by
% SolverMVA.supportsLcfs, the predicate the report gates on, so a pair the
% report calls runnable is one that runs here.
lcfsStat = find(sched == SchedStrategy.LCFS);
if ~isempty(lcfsStat)
    [lcfsOk, lcfsReason] = SolverMVA.supportsLcfs(sn, 'mva');
    if ~lcfsOk
        line_error(mfilename, lcfsReason);
    end
    lcfsprStat = find(sched == SchedStrategy.LCFSPR);
    % Call specialized LCFS MVA solver
    [Q,U,R,T,C,X,lG] = solver_mva_lcfsqn(sn, options, lcfsStat, lcfsprStat);
    return;
end

infSET =[]; % set of infinite server stations
qSET =[]; % set of other product-form stations
for ist=1:M
    switch sched(ist)
        case SchedStrategy.EXT
            % no-op
        case SchedStrategy.INF
            infSET(1,end+1) = ist;
        case {SchedStrategy.PS, SchedStrategy.LCFSPR}
            qSET(1,end+1) = ist;
        case {SchedStrategy.FCFS, SchedStrategy.SIRO}
            qSET(end+1) = ist;
        otherwise
            line_error(mfilename, sprintf('Unsupported exact MVA analysis for %s scheduling.',SchedStrategy.toText(sched(ist))));
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

% Interlocked flow (Franks 1999, Eq. 4.7): a request cannot queue behind work
% that its own submission caused, so the arrival-instant queue drops the
% interlocked share of the other chains. SolverLN supplies the matrix.
IL = [];
if isfield(options,'config') && isfield(options.config,'interlock') && ~isempty(options.config.interlock)
    IL = sn_interlock_chain(sn, options.config.interlock);
end

% the interlocked recursion is a separate entry point: pfqn_mvams and the
% pfqn_mva family it dispatches to carry the standard arrival theorem only
if isempty(IL)
    [Xchain,Qpf,Uchain,~,lG] = pfqn_mvams(lambda,STchain(qSET,:).*Vchain(qSET,:),Nchain,STchain(infSET,:).*Vchain(infSET,:),ones(length(qSET),1),nservers(qSET));
else
    [Xchain,Qpf,Uchain,~,lG] = pfqn_mvams_ilock(lambda,STchain(qSET,:).*Vchain(qSET,:),Nchain,STchain(infSET,:).*Vchain(infSET,:),ones(length(qSET),1),nservers(qSET),IL);
end
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
            switch sched(k)
                case {SchedStrategy.FCFS,SchedStrategy.PS}
                    if sum(Uchain(k,:))>1+options.tol
                        Uchain(k,r) = min(1,sum(Uchain(k,:))) * Vchain(k,r)*STchain(k,r)*Xchain(r) / ((Vchain(k,:).*STchain(k,:))*Xchain(:));
                    end
            end
        end
    end
end

%Vsink = cellsum(sn.nodevisits);
%Vsink = Vsink(find(sn.nodetype==NodeType.Sink),:);
%for r=find(isinf(Nchain)) % open classes
    %Xchain(r) = Vsink(r) ./ STchain(refstatchain(r),r);
%end

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

[Q,U,R,T,C,X] = sn_deaggregate_chain_results(sn, Lchain, [], STchain, Vchain, alpha, [], [], Rchain, Tchain, [], Xchain);
end
