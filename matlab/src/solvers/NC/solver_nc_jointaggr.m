function [Pr,G,runtime] = solver_nc_jointaggr(sn, options)
% [PR,G,RUNTIME] = SOLVER_NC_JOINTAGGR(QN, OPTIONS)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

Tstart = tic;

V = cellsum(sn.visits);
M = sn.nstations;    %number of stations
C = sn.nchains;
Nchain = zeros(1,C);

for c=1:C
    inchain = sn.inchain{c};
    Nchain(c) = sum(sn.njobs(inchain));
end

% transform everything into a LD model
nservers = sn.nservers;
mu = ones(M,sum(Nchain));
for ist=1:M
    if isinf(nservers(ist)) % infinite server
        mu(ist,1:sum(Nchain)) = 1:sum(Nchain);
    else
        mu(ist,1:sum(Nchain)) = min(1:sum(Nchain), nservers(ist)*ones(1,sum(Nchain)));
    end
end
state = sn.state;

switch options.method
    case 'exact'
        Tstart = tic;
        [Lchain,~,~,~,Nchain] = sn_get_demands_chain(sn);
        lG = pfqn_ncld(Lchain, Nchain, 0*Nchain, mu);
        ST = 1 ./ sn.rates;
        ST(isnan(ST))=0;
    otherwise
        % unclear if this is correct as it doesn't consider the
        % transformation to ld model
        [~,~,~,~,~,~,lG,ST] = solver_nc(sn, options);
end
G = exp(lG);

lPr = 0;
for ist=1:M
    isf = sn.stationToStateful(ist);
    [~,nivec] = State.toMarginal(sn, ist, state{isf});
    nivec = unique(nivec,'rows');
    nivec_chain = nivec * sn.chains';
    if any(nivec_chain>0)
        lF_i = pfqn_ncld(ST(ist,:).*V(ist,:), nivec, 0*nivec, mu(ist,:), options);
        lPr = lPr + lF_i ;
    end
end
lPr = lPr - lG;
Pr = exp(lPr);

runtime = toc(Tstart);
return
end

