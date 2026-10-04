function [Q,U,R,T,C,X,lG,iter] = solver_mvald(sn,options)
% [Q,U,R,T,C,X,LG,iter] = SOLVER_MVALD(SN, OPTIONS)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% aggregate chains

if nargin < 2
    options = SolverMVA.defaultOptions;
end

[Lchain,STchain,Vchain,alpha,Nchain,~,refstatchain] = sn_get_demands_chain(sn);
ST = 1 ./ sn.rates;
ST(isnan(ST))=0;

M = size(STchain,1);
C = sn.nchains;
S = sn.nservers;

% Nt is the CLOSED population. On a mixed model sum(Nchain) is Inf, so a column
% range built from it reads 1:Inf and every open load-dependent model under
% 'exact'/'mva' died here before pfqn_mvaldmx, which takes the open chains
% through lambda, ever ran. It is also the width pfqn_mvaldmx validates the rate
% matrix against.
Nt = sum(Nchain(isfinite(Nchain)));

% Chain arrival rates. The reference station of an open chain is its Source and
% STchain holds one over the SUM of the class arrival rates there, so reading
% the rate at chain level also covers a chain whose classes arrive at several
% rates, or one carrying a class that is only ever reached by a switch: that
% class has no arrival process, its own service time at the source is zero and
% the per-class form 1/ST made lambda infinite.
lambda = zeros(1,C);
openChains = find(isinf(Nchain));
for c=openChains(:)'
    if STchain(refstatchain(c),c) > 0
        lambda(c) = 1 / STchain(refstatchain(c),c);
    end
end

if isempty(openChains)
    % PURELY CLOSED. Every station enters the recursion, an infinite server as
    % the load-dependent rate mu(n)=n, which is exact because n cannot then
    % exceed the closed population.
    mu_chain = ones(M,Nt);
    for ist=1:M
        if isinf(S(ist)) % infinite server
            mu_chain(ist,1:Nt) = 1:Nt;
        elseif ~isempty(sn.lldscaling)
            mu_chain(ist,1:Nt) = sn.lldscaling(ist,1:Nt);
        end
    end
    [Xchain,Qchain,Uchain] = pfqn_mvaldmx(lambda,Lchain,Nchain,0*Nchain,mu_chain,S);
else
    % MIXED OR PURELY OPEN. Three kinds of row are not the same thing to
    % pfqn_mvaldmx and have to be separated before it is called. This is the
    % partition solver_ncld makes for the same recursion.
    %  - THE SOURCE IS NOT A STATION. Its chain demand is the interarrival time
    %    1/lambda, so it carries offered load Lo=1 exactly; pfqn_ldmx_ec then
    %    forms 1/(1-Lo/mu)=Inf and the 0*Inf in the residence-time sum turned
    %    every chain into NaN, which the deaggregation reported as zeros.
    %  - A DELAY IS AN INFINITE SERVER FOR THE OPEN CHAINS TOO. mu(n)=n cut at
    %    the closed population declares it saturated at Nt jobs, i.e. an
    %    Nt-server queue. It enters as chain think time instead and its queue
    %    length is X*L, which is exact.
    %  - A QUEUEING STATION KEEPS ITS WHOLE RATE ROW. pfqn_ldmx_ec reads the
    %    limited-load-dependence level b off the row itself, so a row cut at the
    %    closed population is read as a slower station, and with no closed class
    %    at all it collapses to mu(1), a single fixed-rate server.
    sourceStations = reshape(unique(refstatchain(openChains)),1,[]);
    delayStations = setdiff(reshape(find(isinf(S)),1,[]), sourceStations);
    queueStations = setdiff(1:M, [delayStations, sourceStations]);
    nq = numel(queueStations);
    Zchain = zeros(1,C);
    for c=1:C
        Zchain(c) = sum(Lchain(delayStations,c));
    end
    lldWidth = size(sn.lldscaling,2);
    ncol = max(1,Nt);
    for qi=1:nq
        % first column of the trailing constant run, the level b of pfqn_ldmx_ec
        b = lldWidth;
        while b > 1 && sn.lldscaling(queueStations(qi),b-1) == sn.lldscaling(queueStations(qi),b)
            b = b - 1;
        end
        ncol = max(ncol, b);
    end
    mu_chain = ones(nq,ncol);
    for qi=1:nq
        if lldWidth > 0
            avail = min(ncol, lldWidth);
            mu_chain(qi,1:avail) = sn.lldscaling(queueStations(qi),1:avail);
            mu_chain(qi,(avail+1):ncol) = sn.lldscaling(queueStations(qi),lldWidth); % saturated tail
        end
    end
    [Xchain,Qqueue,Uqueue] = pfqn_mvaldmx(lambda,Lchain(queueStations,:),Nchain,Zchain,mu_chain,S(queueStations));
    Qchain = zeros(M,C);
    Uchain = zeros(M,C);
    Qchain(queueStations,:) = Qqueue;
    Uchain(queueStations,:) = Uqueue;
    for di=1:numel(delayStations)
        % infinite server: X*L for a closed chain, lambda*L for an open one
        Qchain(delayStations(di),:) = Lchain(delayStations(di),:) .* Xchain;
    end
end
Tchain = repmat(Xchain,M,1) .* Vchain;
% Rchain is the PER-VISIT response time, Qchain./Tchain, not the residence time
% Qchain./Xchain. sn_deaggregate_chain_results rebuilds Q(i,k) as
% Rchain*Xchain*Vchain(i,c)/Vchain(refstat,c), i.e. it multiplies the visit
% ratio back in; dividing by Xchain here counted that ratio twice and the class
% queue lengths then did not sum to the population. Invisible whenever the LD
% station carries the reference station's chain visits (every shipped LD
% example), decisive when it does not (the closed delayed-hit retrieval cache,
% whose fetch station is visited once per miss). solver_mva.m:153 and the Python
% exact LD path already divide by Tchain.
Rchain = Qchain ./ Tchain;
Rchain(~isfinite(Rchain)) = 0;
lG = NaN;

%% This is likely wrong as it uses Little's law for the utilization computation
[Q,U,R,T,C,X] = sn_deaggregate_chain_results(sn, Lchain, [], STchain, Vchain, alpha, [], Uchain, Rchain, Tchain, [], Xchain);

% Busy-server FRACTION under load-dependent scaling (NC convention): report
% the carried load over the effective capacity ceff = max(nservers, peak
% lldscaling) instead of the mvaldmx P(busy)-style estimator.
if ~isempty(sn.lldscaling)
    for ist=1:M
        if isfinite(S(ist))
            ceff = max(S(ist), max(sn.lldscaling(ist,:)));
            for r=1:sn.nclasses
                if ST(ist,r) > 0
                    U(ist,r) = T(ist,r) * ST(ist,r) / ceff;
                end
            end
        end
    end
end
iter = 1;
end
