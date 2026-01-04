function [Q,U,R,T,C,X,lG,iter] = solver_mvald(sn,options)
% [Q,U,R,T,C,X,LG,iter] = SOLVER_MVALD(SN, OPTIONS)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% aggregate chains

if nargin < 2
    options = SolverMVA.defaultOptions;
end

[Lchain,STchain,Vchain,alpha,Nchain,~,~] = sn_get_demands_chain(sn);
ST = 1 ./ sn.rates;
ST(isnan(ST))=0;

M = size(STchain,1);
C = sn.nchains;
S = sn.nservers;
N = sn.njobs;

mu_chain = ones(M,sum(Nchain(isfinite(Nchain))));
for ist=1:M
    if isinf(S(ist)) % infinite server
        mu_chain(ist,1:sum(Nchain)) = 1:sum(Nchain(isfinite(Nchain)));
    elseif ~isempty(sn.lldscaling)
        mu_chain(ist,1:sum(Nchain)) = sn.lldscaling(ist,1:sum(Nchain(isfinite(Nchain))));
    else
        mu_chain(ist,1:sum(Nchain)) = ones(1,sum(Nchain(isfinite(Nchain))));
    end
end

lambda = zeros(1,C);
for c=1:sn.nchains
    ocl = find(isinf(N) & sn.chains(c,:));  % open classes in this chain
    if ~isempty(ocl)
        Nchain(c) = Inf;
    end
    for r=ocl
        lambda(c) = lambda(c) + 1 ./ ST(refstat(r),r);
    end
end
[Xchain,Qchain,Uchain] = pfqn_mvaldmx(lambda,Lchain,Nchain,0*Nchain,mu_chain,S);
Tchain = repmat(Xchain,M,1) .* Vchain;
Rchain = Qchain ./ repmat(Xchain,M,1);
lG = NaN;

%% This is likely wrong as it uses Little's law for the utilization computation
[Q,U,R,T,C,X] = sn_deaggregate_chain_results(sn, Lchain, [], STchain, Vchain, alpha, [], Uchain, Rchain, Tchain, [], Xchain);
iter = 1;
end
