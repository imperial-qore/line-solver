function [Q,U,R,T,C,X,lG,STeff,it,method] = solver_nc(sn, options)

M = sn.nstations;    %number of stations
nservers = sn.nservers;
NK = sn.njobs';  % initial population per class
sched = sn.sched;
C = sn.nchains;
SCV = sn.scv;
V = cellsum(sn.visits);
ST = 1 ./ sn.rates;
ST(isnan(ST)) = 0;
ST0 = ST;

% Check for special LCFS + LCFS-PR 2-station network
lcfsStat = find(sched == SchedStrategy.LCFS);
lcfsprStat = find(sched == SchedStrategy.LCFSPR);
if ~isempty(lcfsStat) && ~isempty(lcfsprStat)
    % Validate LCFS network topology
    if length(lcfsStat) ~= 1 || length(lcfsprStat) ~= 1
        line_error(mfilename, 'LCFS NC requires exactly one LCFS and one LCFS-PR station.');
    end
    Nchain = zeros(1,C);
    for c=1:C
        inchain = sn.inchain{c};
        Nchain(c) = sum(NK(inchain)); %#ok<FNDSB>
    end
    if any(isinf(Nchain))
        line_error(mfilename, 'LCFS NC requires a closed queueing network.');
    end
    % Check for self-loops in routing matrix
    rt = sn.rt;
    nclasses = sn.nclasses;
    for ist = [lcfsStat, lcfsprStat]
        for r = 1:nclasses
            if rt((ist-1)*nclasses+r, (ist-1)*nclasses+r) > 0
                line_error(mfilename, 'LCFS NC does not support self-loops at stations.');
            end
        end
    end
    % Call specialized LCFS NC solver
    [Q,U,R,T,C,X,lG] = solver_nc_lcfsqn(sn, options, lcfsStat, lcfsprStat);
    STeff = ST;
    it = 1;
    method = 'lcfsqn.ca';
    return;
elseif ~isempty(lcfsStat)
    % LCFS without LCFS-PR is not supported
    line_error(mfilename, 'LCFS scheduling requires a paired LCFS-PR station.');
end

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

if all(sched~=SchedStrategy.FCFS) options.iter_max=1; end
it = 0;
while max(abs(1-eta./eta_1)) > options.iter_tol & it < options.iter_max
    it = it + 1;
    eta_1 = eta;

    if it==1
        lambda = zeros(1,C);
        [Lchain,STchain,Vchain,alpha,Nchain] = sn_get_demands_chain(sn);
        for c=1:C
            inchain = sn.inchain{c};
            isOpenChain = any(isinf(sn.njobs(inchain)));
            for ist=1:M
                % we assume that the visits in L(i,inchain) are equal to 1
                if isOpenChain && ist == sn.refstat(inchain(1)) % if this is a source ST = 1 / arrival rates
                    lambda(c) = 1 ./ STchain(ist,c);
                end
            end
        end
    else
        for c=1:C
            inchain = sn.inchain{c};
            for ist=1:M
                % we assume that the visits in L(i,inchain) are equal to 1
                STchain(ist,c) = ST(ist,inchain) * alpha(ist,inchain)';
                Lchain(ist,c) = Vchain(ist,c) * STchain(ist,c);
            end
        end
    end

    STchain(~isfinite(STchain))=0;
    Lchain(~isfinite(Lchain))=0;

    Lms = zeros(M,C);
    Z = zeros(M,C);
    Zms = zeros(M,C);
    infServers = [];
    for ist=1:M
        if isinf(nservers(ist)) % infinite server
            infServers(end+1) = ist;
            Lms(ist,:) = 0;
            Z(ist,:) = Lchain(ist,:);
            Zms(ist,:) = 0;
        else
            if strcmpi(options.method,'exact') && nservers(ist)>1
                %options.method = 'default';
                if options.verbose
                    line_warning(mfilename,sprintf('%s does not support exact multiserver yet. Switching to approximate method.\n', 'SolverNC'));
                end
            end
            % Seidmann's approximation
            Lms(ist,:) = Lchain(ist,:) / nservers(ist);
            Z(ist,:) = 0;
            Zms(ist,:) = Lchain(ist,:) * (nservers(ist)-1)/nservers(ist);
        end
    end

    % step 1
    [lG, Xchain, Qchain, method] = pfqn_nc(lambda,Lms,Nchain,sum(Z,1)+sum(Zms,1), options);

    if sum(Zms,1) > GlobalConstants.FineTol
        % in this case, we need to use the iterative approximation below
        Xchain=[];
        Qchain=[];
    end

    if isempty(lG)
        [Q,U,R,T,C,X,lG,STeff,it,method] = deal([],[],[],[],[],[],[],[],1,options.method);
        return
    end

    if isempty(Xchain)
        Xchain=lambda;
        Qchain = zeros(M,C);
        for r=closedChains
            lGr(r) = pfqn_nc(lambda,Lms,oner(Nchain,r),sum(Z,1)+sum(Zms,1), options);
            Xchain(r) = exp(lGr(r) - lG);
            for ist=1:M
                if Lchain(ist,r)>0
                    if isinf(nservers(ist)) % infinite server
                        Qchain(ist,r) = Lchain(ist,r) * Xchain(r);
                    else
                        % add replica of station i and move job of class r
                        % in separate class
                        % this is the most costly operation so we record
                        % method here
                        [lGar(ist,r),~,~,method] = pfqn_nc([lambda,0],[Lms(setdiff(1:size(Lms,1),ist),:),zeros(size(Lms,1)-1,1); Lms(ist,:),1], [oner(Nchain,r),1], [sum(Z,1)+sum(Zms,1),0], options);
                        Qchain(ist,r) = Zms(ist,r) * Xchain(r) + Lms(ist,r) * exp(lGar(ist,r) - lG);
                    end
                end
            end
            Qchain(isnan(Qchain))=0;
        end
        for r=openChains
            for ist=1:M
                Qchain(ist,r) = lambda(r)*Lchain(ist,r)/(1-lambda(openChains)*Lchain(ist,openChains)'/nservers(ist))*(1+sum(Qchain(ist,closedChains)));
            end
        end
    else
        % just fill the delay servers
        for r=1:C
            for ist=1:M
                if Lchain(ist,r)>0
                    if isinf(nservers(ist)) % infinite server
                        Qchain(ist,r) = Lchain(ist,r) * Xchain(r);
                    end
                end
            end
        end
    end

    if isnan(Xchain)
        if options.verbose
            line_warning(mfilename,'Normalizing constant computations produced a floating-point range exception. Model is likely too large.\n');
        end
    end

    Z = sum(Z(1:M,:),1);

    Rchain = Qchain ./ repmat(Xchain,M,1) ./ Vchain;
    Rchain(infServers,:) = Lchain(infServers,:) ./ Vchain(infServers,:);
    Tchain = repmat(Xchain,M,1) .* Vchain;

    %Uchain = Tchain .* Lchain;
    %Cchain = Nchain ./ Xchain - Z;

    [Q,U,R,T,~,X] = sn_deaggregate_chain_results(sn, Lchain, ST, STchain, Vchain, alpha, [], [], Rchain, Tchain, [], Xchain);
    STeff = ST;% effective service time at the last iteration
    [ST,gamma,~,~,~,~,eta] = npfqn_nonexp_approx(options.config.highvar,sn,ST0,V,SCV,T,U,gamma,nservers);
end

Q=abs(Q); R=abs(R); X=abs(X); U=abs(U);
X(~isfinite(X))=0; U(~isfinite(U))=0; Q(~isfinite(Q))=0; R(~isfinite(R))=0;

% renormalize qlen and tput to correct for unforeseen population constraint deviations
for c=1:sn.nchains
    inchain = sn.inchain{c};
    Nchain(c) = sum(NK(inchain)); %#ok<FNDSB>
    if isfinite(Nchain(c))
        q_den = sum(sum(Q(:,inchain)));
        if  q_den > 0
            ratio = Nchain(c)/ q_den;
        else
            ratio = 0;
        end
        Q(:,inchain) = ratio * Q(:,inchain);
        X(inchain) = ratio * X(inchain);
        T(:,inchain) = ratio * T(:,inchain);
        U(:,inchain) = ratio * U(:,inchain);
        R(:,inchain) = Q(:,inchain) ./ T(:,inchain);
    end
end
end
