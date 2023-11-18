function [Pr,G,lG,runtime] = solver_nc_joint(sn, options)
% [PR,G,LG,RUNTIME] = SOLVER_NC_JOINT(QN, OPTIONS)

% Copyright (c) 2012-2023, Imperial College London
% All rights reserved.

%% initialization
state = sn.state;
S = sn.nservers;
rates = sn.rates;
% determine service times
ST  = 1./rates;
ST(isnan(ST))=0;

[~,STchain,Vchain,alpha,Nchain] = snGetDemandsChain(sn);

Tstart = tic;

[M,K]=size(STchain);

Lchain = zeros(M,K);
mu_chain = ones(M,sum(Nchain));
for i=1:M
    Lchain(i,:) = STchain(i,:) .* Vchain(i,:);
    if isinf(S(i)) % infinite server
        mu_chain(i,1:sum(Nchain)) = 1:sum(Nchain);
    else
        mu_chain(i,1:sum(Nchain)) = min(1:sum(Nchain), S(i)*ones(1,sum(Nchain)));
    end
end

lG = pfqn_ncld(Lchain, Nchain, 0*Nchain, mu_chain, options);
lPr = 0;
for i=1:M
    isf = sn.stationToStateful(i);
    [~,nivec] = State.toMarginal(sn, i, state{isf});
    nivec_chain = nivec * sn.chains';
    lF_i = pfqn_ncld(Lchain(i,:), nivec_chain, 0*nivec_chain,mu_chain(i,:), options);
    lg0_i = pfqn_ncld(ST(i,:).*alpha(i,:), nivec, 0*nivec, mu_chain(i,:), options);
    lG0_i = pfqn_ncld(STchain(i,:), nivec_chain, 0*nivec_chain, mu_chain(i,:), options);
    lPr = lPr + lF_i + (lg0_i - lG0_i);
end
Pr = exp(lPr - lG);

runtime = toc(Tstart);

G = exp(lG);
%if options.verbose
%    line_printf('\nNormalizing constant (NC) analysis completed. Runtime: %f seconds.\n',runtime);
%end
end