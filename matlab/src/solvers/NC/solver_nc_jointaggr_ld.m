function [Pr,G,lG,runtime] = solver_nc_jointaggr_ld(sn, options)
% [PR,G,LG,RUNTIME] = SOLVER_NC_JOINTAGGR(QN, OPTIONS)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

%% initialization
state = sn.state;
S = sn.nservers;
rates = sn.rates;
% determine service times
ST  = 1./rates;
ST(isnan(ST))=0;
[Lchain,STchain,~,alpha,Nchain] = sn_get_demands_chain(sn);

Tstart = tic;

[M,K]=size(STchain);

mu_chain = ones(M,sum(Nchain));
for ist=1:M
    if isinf(S(ist)) % infinite server
        mu_chain(ist,1:sum(Nchain)) = 1:sum(Nchain);
    else
        mu_chain(ist,1:sum(Nchain)) = min(1:sum(Nchain), S(ist)*ones(1,sum(Nchain)));
    end
end

lG = pfqn_ncld(Lchain, Nchain, 0*Nchain, mu_chain);
lPr = 0;
for ist=1:M
    isf = sn.stationToStateful(ist);
    [~,nivec] = State.toMarginal(sn, ist, state{isf});
    nivec_chain = nivec * sn.chains';
    lF_i = pfqn_ncld(Lchain(ist,:), nivec_chain, 0*nivec_chain, mu_chain(ist,:), options);
    lg0_i = pfqn_ncld(ST(ist,:).*alpha(ist,:), nivec, 0*nivec, mu_chain(ist,:), options);
    lG0_i = pfqn_ncld(STchain(ist,:),nivec_chain, 0*nivec_chain, mu_chain(ist,:), options);
    lPr = lPr + lF_i + (lg0_i - lG0_i);
end
Pr = exp(lPr - lG);
G=exp(lG);
runtime = toc(Tstart);

%if options.verbose
%    line_printf('\nNormalizing constant (NC) analysis completed. Runtime: %f seconds.\n',runtime);
%end
end
