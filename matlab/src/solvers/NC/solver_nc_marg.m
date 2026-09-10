function [lPr,G,runtime] = solver_nc_marg(sn, options, lG)
% [PR,G,RUNTIME] = SOLVER_NC_MARG(QN, OPTIONS)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.


%% initialization
M = sn.nstations;    %number of stations
K = sn.nclasses;    %number of classes
state = sn.state;
S = sn.nservers;
V = cellsum(sn.visits);
rates = sn.rates;
ST  = 1./rates;
ST(isnan(ST))=0;


[Lchain,STchain,~,~,Nchain] = sn_get_demands_chain(sn);

Tstart = tic;

[M,K]=size(STchain);

mu = ones(M,sum(Nchain));
for ist=1:M
    if isinf(S(ist)) % infinite server
        mu(ist,1:sum(Nchain)) = 1:sum(Nchain);
    else
        mu(ist,1:sum(Nchain)) = min(1:sum(Nchain), S(ist)*ones(1,sum(Nchain)));
    end
end

if nargin < 3
    lG = pfqn_ncld(Lchain, Nchain, 0*Nchain, mu, options);
end
G = exp(lG);

lPr = zeros(sn.nstations,1);
for ist=1:sn.nstations
    ind = sn.stationToNode(ist);
    isf = sn.stationToStateful(ist);
    [~,nirvec,sivec,kirvec] = State.toMarginal(sn, ind, state{isf});
    if min(nirvec) < 0 % user flags that state of i should be ignored
        lPr(ist) = NaN;
    else
        set_ist = setdiff(1:sn.nstations,ist);
        nivec_chain = nirvec * sn.chains';
        lG_minus_i = pfqn_ncld(Lchain(set_ist,:), Nchain-nivec_chain, 0*Nchain, mu(set_ist,:), options);
        lF_i = 0;
        switch sn.sched(ist)
            case SchedStrategy.FCFS
                for r=1:K
                    PHr = sn.proc{ist}{r};
                    if ~isempty(PHr)
                        kir = kirvec(1,r,:); kir=kir(:)';
                        if length(kir)>1
                            line_error(mfilename,'Cannot return state probability because the product-form solution requires exponential service times at FCFS nodes.');
                        end
                        if ST(ist,r)~=max(ST(ist,:))
                            line_error(mfilename,'Cannot return state probability because the product-form solution requires identical service times at FCFS nodes.');
                        end
                    end
                end
                ci = find(sivec, 1);
                if ~isempty(ci)
                    lF_i = lF_i + sum((nirvec(1,:).*log(V(ist,r)))) - sum(log(mu(ist,1:sum(kirvec(:)))));
                else
                    lF_i = 0;
                end
            case SchedStrategy.SIRO
                for r=1:K
                    PHr = sn.proc{ist}{r};
                    if ~isempty(PHr)
                        kir = kirvec(1,r,:); kir=kir(:)';
                        if length(kir)>1
                            line_error(mfilename,'Cannot return state probability because the product-form solution requires exponential service times at RAND nodes.');
                        end
                        if ST(ist,r)~=max(ST(ist,:))
                            line_error(mfilename,'Cannot return state probability because the product-form solution requires identical service times at RAND nodes.');
                        end
                    end
                end
                ci = find(sivec);
                if ~isempty(ci)
                    lF_i = lF_i + log(nirvec(ci)) - log(sum(nirvec)) + pfqn_ncld(ST(ist,:).*V(ist,:), nirvec, 0*nirvec, mu(ist,:), options);
                else
                    lF_i = 0;
                end
            case SchedStrategy.PS
                for r=1:K
                    PHr = sn.proc{ist}{r};
                    if ~isempty(PHr)
                        kir = kirvec(1,r,:); kir=kir(:)';
                        Ar = -PHr{1}\ map_pie(PHr);
                        lF_i = lF_i + sum(kir.*log(V(ist,r)*Ar)) - sum(factln(kir));
                    end
                end
                lF_i = lF_i + factln(sum(kirvec(:))) - sum(log(mu(ist,1:sum(kirvec(:)))));
            case SchedStrategy.INF
                for r=1:K
                    PHr = sn.proc{ist}{r};
                    if ~isempty(PHr)
                        kir = kirvec(1,r,:); kir=kir(:)';
                        Ar = -PHr{1}\map_pie(PHr);
                        lF_i = lF_i + sum(kir.*log(V(ist,r)*Ar)) - sum(factln(kir));
                    end
                end
        end
        lPr(ist) =  lF_i + lG_minus_i - lG;
    end
end

runtime = toc(Tstart);
lPr(isnan(lPr))=0;
%if options.verbose
%    line_printf('\nNormalizing constant (NC) analysis completed. Runtime: %f seconds.\n',runtime);
%end
end
