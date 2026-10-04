function [Pr,G,lG,runtime] = solver_nc_margaggr(sn, options, lG)
% [PR,G,LG,RUNTIME] = SOLVER_NC_MARGAGGR(QN, OPTIONS)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

M = sn.nstations;    %number of stations
K = sn.nclasses;    %number of classes
state = sn.state;
S = sn.nservers;
NK = sn.njobs';  % initial population per class
C = sn.nchains;

PH = sn.proc;

%% initialization

% determine service times
ST = zeros(M,K);
for k = 1:K
    for ist=1:M
        ST(ist,k) = 1 ./ map_lambda(PH{ist}{k});
    end
end
ST(isnan(ST))=0;

[Lchain,STchain,~,~,Nchain] = sn_get_demands_chain(sn);

V = zeros(sn.nstations,sn.nclasses);
for c=1:sn.nchains
    inchain = sn.inchain{c};
    for ist=1:sn.nstations
        for k=inchain
            V(ist,k) = sn.visits{c}(ist,k);
        end
    end
end

Tstart = tic;

[M,~]=size(STchain);

mu = ones(M,sum(Nchain));
for ist=1:M
    if isinf(S(ist)) % infinite server
        mu(ist,1:sum(Nchain)) = 1:sum(Nchain);
    else
        mu(ist,1:sum(Nchain)) = min(1:sum(Nchain), S(ist)*ones(1,sum(Nchain)));
    end
end

if nargin == 2
    lG = pfqn_ncld(Lchain, Nchain, 0*Nchain, mu, options);
end
G = exp(lG);

lPr = zeros(sn.nstations,1);
for ist=1:sn.nstations
    ind = sn.stationToNode(ist);
    isf = sn.stationToStateful(ist);
    [~,nivec] = State.toMarginal(sn, ind, state{isf});
    if min(nivec) < 0 % user flags that state of i should be ignored
        lPr(ist) = NaN;
    else
        set_ist = setdiff(1:sn.nstations,ist);
        nivec_chain = nivec * sn.chains';
        lG_minus_i = pfqn_ncld(Lchain(set_ist,:), Nchain-nivec_chain, 0*Nchain, mu(set_ist,:), options);
        lF_i = pfqn_ncld(ST(ist,:).*V(ist,:), nivec, 0*nivec, mu(ist,:), options);
        lPr(ist) =  lF_i + lG_minus_i - lG;
    end
end
Pr = exp(lPr);
Pr(isnan(Pr))=0;
lG = log(G);
runtime = toc(Tstart);
%if options.verbose
%    line_printf('Normalizing constant (NC) analysis completed. Runtime: %f seconds.\n',runtime);
%end
end