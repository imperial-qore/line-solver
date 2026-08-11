function [Q,U,R,T,C,X,lG,iter] = solver_mva_sum(sn, options)
% [Q,U,R,T,C,X,LG,ITER] = SOLVER_MVA_SUM(SN, OPTIONS)
%
% Summation method (SUM/ESUM) analyzer. Closed models are solved with
% sum_closed; open and mixed models with sum_closing (closing method,
% Kclosed=5000). FCFS and SIRO stations use their service SCVs (ESUM
% corrections for scv~=1); PS, LCFS-PR and infinite-server stations are
% insensitive and are passed scv=1. See api/sum and Bolch et al., Secs.
% 9.2, 10.1.4.4 and 10.1.5.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

[Lchain,STchain,Vchain,alpha,Nchain,SCVchain,refstatchain] = sn_get_demands_chain(sn); %#ok<ASGLU>

M = sn.nstations;
K = sn.nchains;
sched = sn.sched;
nservers = sn.nservers;

% station rows passed to the summation method (all but the source)
rows = [];
mi = [];
for ist=1:M
    switch sched(ist)
        case SchedStrategy.EXT
            % no-op, external world handled by lambda/sum_closing
        case SchedStrategy.INF
            rows(1,end+1) = ist;
            mi(end+1,1) = Inf;
        case {SchedStrategy.PS, SchedStrategy.LCFSPR, SchedStrategy.FCFS, SchedStrategy.SIRO}
            rows(1,end+1) = ist;
            mi(end+1,1) = nservers(ist);
        otherwise
            line_error(mfilename, sprintf('The summation method does not support %s scheduling.', SchedStrategy.toText(sched(ist))));
    end
end

L = STchain(rows,:).*Vchain(rows,:);
scv = ones(length(rows),K);
for j=1:length(rows)
    ist = rows(j);
    if any(sched(ist) == [SchedStrategy.FCFS, SchedStrategy.SIRO])
        for r=1:K
            if isfinite(SCVchain(ist,r)) && SCVchain(ist,r) > 0
                scv(j,r) = SCVchain(ist,r);
            end
        end
    end
end

Z = zeros(1,K);
ocl = find(isinf(Nchain));
if isempty(ocl)
    [Xchain,Qrows,Urows,~,iter] = sum_closed(L,Nchain,Z,mi,scv,options.iter_tol,options.iter_max);
else
    lambda = zeros(1,K);
    scva = ones(1,K);
    for r=ocl(:)'
        lambda(r) = 1 ./ STchain(refstatchain(r),r);
        if isfinite(SCVchain(refstatchain(r),r)) && SCVchain(refstatchain(r),r) > 0
            scva(r) = SCVchain(refstatchain(r),r); % interarrival SCV at the source
        end
    end
    [Xchain,Qrows,Urows,~,~,iter] = sum_closing(lambda,scva,L,mi,scv,Nchain,Z,5000,options.iter_tol,options.iter_max);
end

Qchain = zeros(M,K);
Uchain = zeros(M,K);
Tchain = zeros(M,K);
Qchain(rows,:) = Qrows;
Uchain(rows,:) = Urows;
for k=1:M
    for r=1:K
        Tchain(k,r) = Xchain(r) * Vchain(k,r);
    end
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

lG = NaN;
[Q,U,R,T,C,X] = sn_deaggregate_chain_results(sn, Lchain, [], STchain, Vchain, alpha, [], [], Rchain, Tchain, [], Xchain);
end
