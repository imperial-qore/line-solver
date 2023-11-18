function [Q,U,R,T,C,X,lG,STeff,it,method] = solver_nc(sn, options)
M = sn.nstations;    %number of stations
nservers = sn.nservers;
NK = sn.njobs';  % initial population per class
schedid = sn.schedid;
C = sn.nchains;
SCV = sn.scv;
V = cellsum(sn.visits);
ST = 1 ./ sn.rates;
ST(isnan(ST)) = 0;
ST0 = ST;

Nchain = zeros(1,C);
for c=1:C
    inchain = sn.inchain{c};
    Nchain(c) = sum(NK(inchain)); %#ok<FNDSB>
end
openChains = find(isinf(Nchain));
closedChains = find(~isinf(Nchain));

gamma = zeros(1,M);
eta_1 = zeros(1,M);
eta = ones(1,M);
C = sn.nchains;

if all(schedid~=SchedStrategy.ID_FCFS) options.iter_max=1; end
it = 0;
while max(abs(1-eta./eta_1)) > options.iter_tol & it < options.iter_max
    it = it + 1;
    eta_1 = eta;

    if it==1
        lambda = zeros(1,C);
        [Lchain,STchain,Vchain,alpha,Nchain] = snGetDemandsChain(sn);
        for c=1:C
            inchain = sn.inchain{c};
            isOpenChain = any(isinf(sn.njobs(inchain)));
            for i=1:M
                % we assume that the visits in L(i,inchain) are equal to 1
                if isOpenChain && i == sn.refstat(inchain(1)) % if this is a source ST = 1 / arrival rates
                    lambda(c) = 1 ./ STchain(i,c);
                end
            end
        end
    else
        for c=1:C
            inchain = sn.inchain{c};
            for i=1:M
                % we assume that the visits in L(i,inchain) are equal to 1
                STchain(i,c) = ST(i,inchain) * alpha(i,inchain)';
                Lchain(i,c) = Vchain(i,c) * STchain(i,c);
            end
        end
    end

    STchain(~isfinite(STchain))=0;
    Lchain(~isfinite(Lchain))=0;

    Lms = zeros(M,C);
    Z = zeros(M,C);
    Zms = zeros(M,C);
    infServers = [];
    for i=1:M
        if isinf(nservers(i)) % infinite server
            infServers(end+1) = i;
            Lms(i,:) = 0;
            Z(i,:) = Lchain(i,:);
            Zms(i,:) = 0;
        else
            if strcmpi(options.method,'exact') && nservers(i)>1
                %options.method = 'default';
                line_warning(mfilename,sprintf('%s does not support exact multiserver yet. Switching to approximate method.\n', 'SolverNC'));
            end
            % Seidmann's approximation
            Lms(i,:) = Lchain(i,:) / nservers(i);
            Z(i,:) = 0;
            Zms(i,:) = Lchain(i,:) * (nservers(i)-1)/nservers(i);
        end
    end

    % step 1
    [lG, Xchain, Qchain, method] = pfqn_nc(lambda,Lms,Nchain,sum(Z,1)+sum(Zms,1), options);

    if sum(Zms,1) > GlobalConstants.FineTol
        % in this case, we need to use the iterative approximation below
        Xchain=[];
        Qchain=[];
    end

    if isempty(Xchain)
        Xchain=lambda;
        Qchain = zeros(M,C);
        for r=closedChains
            lGr(r) = pfqn_nc(lambda,Lms,oner(Nchain,r),sum(Z,1)+sum(Zms,1), options);
            Xchain(r) = exp(lGr(r) - lG);
            for i=1:M
                if Lchain(i,r)>0
                    if isinf(nservers(i)) % infinite server
                        Qchain(i,r) = Lchain(i,r) * Xchain(r);
                    else
                        % add replica of station i and move job of class r
                        % in separate class
                        lGar(i,r) = pfqn_nc([lambda,0],[Lms(setdiff(1:size(Lms,1),i),:),zeros(size(Lms,1)-1,1); Lms(i,:),1], [oner(Nchain,r),1], [sum(Z,1)+sum(Zms,1),0], options);
                        Qchain(i,r) = Zms(i,r) * Xchain(r) + Lms(i,r) * exp(lGar(i,r) - lG);
                    end
                end
            end
            Qchain(isnan(Qchain))=0;
        end
        for r=openChains
            for i=1:M
                Qchain(i,r) = lambda(r)*Lchain(i,r)/(1-lambda(openChains)*Lchain(i,openChains)'/nservers(i))*(1+sum(Qchain(i,closedChains)));
            end
        end
    else
        % just fill the delay servers
        for r=1:C
            for i=1:M
                if Lchain(i,r)>0
                    if isinf(nservers(i)) % infinite server
                        Qchain(i,r) = Lchain(i,r) * Xchain(r);
                    end
                end
            end
        end
    end

    if isnan(Xchain)
        line_warning(mfilename,'Normalizing constant computations produced a floating-point range exception. Model is likely too large.\n');
    end

    Z = sum(Z(1:M,:),1);

    Rchain = Qchain ./ repmat(Xchain,M,1) ./ Vchain;
    Rchain(infServers,:) = Lchain(infServers,:) ./ Vchain(infServers,:);
    Tchain = repmat(Xchain,M,1) .* Vchain;
    
    %Uchain = Tchain .* Lchain;
    %Cchain = Nchain ./ Xchain - Z;

    [Q,U,R,T,~,X] = snDeaggregateChainResults(sn, Lchain, ST, STchain, Vchain, alpha, [], [], Rchain, Tchain, [], Xchain);
    STeff = ST;% effective service time at the last iteration
    [ST,gamma,~,~,~,~,eta] = npfqn_nonexp_approx(options.config.highvar,sn,ST0,V,SCV,T,U,gamma,nservers);    
end

Q=abs(Q); R=abs(R); X=abs(X); U=abs(U);
X(~isfinite(X))=0; U(~isfinite(U))=0; Q(~isfinite(Q))=0; R(~isfinite(R))=0;
end
