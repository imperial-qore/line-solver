function [Q,U,R,T,C,X,lG,runtime,iter,method] = solver_ncld(sn, options)
% [Q,U,R,T,C,X,LG,RUNTIME,ITER.METHOD] = SOLVER_NCLD(QN, OPTIONS)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
M = sn.nstations;    %number of stations
K = sn.nclasses;
nservers = sn.nservers;
if nservers(isfinite(nservers))>1
    if isempty(sn.lldscaling) && M==2 && all(isfinite(sn.njobs))
        for ist=1:M
            Nt = sum(sn.njobs);
            sn.lldscaling(ist,1:Nt) = min(1:Nt,sn.nservers(ist));
        end
    else
        line_error(mfilename,'The load-dependent solver does not support multi-server stations yet. Specify multi-server stations via limited load-dependence.');
    end
end

if ~isempty(sn.cdscaling) && strcmpi(options.method, 'exact')
    line_error(mfilename,'Exact class-dependent solver not yet available in NC.');
end

if sn_has_joint_dependence(sn) && strcmpi(options.method, 'exact')
    line_error(mfilename,'Exact joint-dependent solver not yet available in NC.');
end

NK = sn.njobs';  % initial population per class

if isinf(max(NK))
    line_error(mfilename,'The load-dependent solver does not support open classes yet.');
end

sched = sn.sched;
%chains = sn.chains;
C = sn.nchains;
SCV = sn.scv;
gamma = zeros(M,1);
V = cellsum(sn.visits);
ST = 1 ./ sn.rates;
ST(isnan(ST))=0;
ST0=ST;
lldscaling = sn.lldscaling;
Nt = sum(NK(isfinite(NK)));
if isempty(lldscaling)
    lldscaling = ones(M,ceil(Nt));
end

[~,~,Vchain,alpha] = sn_get_demands_chain(sn);

eta_1 = zeros(1,M);
eta = ones(1,M);
if all(sched~=SchedStrategy.FCFS) options.iter_max=1; end
iter = 0;
while max(abs(1-eta./eta_1)) > options.iter_tol & iter < options.iter_max
    iter = iter + 1;
    eta_1 = eta;
    M = sn.nstations;    %number of stations
    K = sn.nclasses;    %number of classes
    C = sn.nchains;
    Lchain = zeros(M,C);
    STchain = zeros(M,C);

    SCVchain = zeros(M,C);
    Nchain = zeros(1,C);
    refstatchain = zeros(C,1);
    for c=1:C
        inchain = sn.inchain{c};
        isOpenChain = any(isinf(sn.njobs(inchain)));
        for ist=1:M
            % we assume that the visits in L(i,inchain) are equal to 1
            Lchain(ist,c) = Vchain(ist,c) * ST(ist,inchain) * alpha(ist,inchain)';
            STchain(ist,c) = ST(ist,inchain) * alpha(ist,inchain)';
            if isOpenChain && ist == sn.refstat(inchain(1)) % if this is a source ST = 1 / arrival rates
                STchain(ist,c) = sumfinite(ST(ist,inchain)); % ignore degenerate classes with zero arrival rates
            else
                STchain(ist,c) = ST(ist,inchain) * alpha(ist,inchain)';
            end
            SCVchain(ist,c) = SCV(ist,inchain) * alpha(ist,inchain)';
        end
        Nchain(c) = sum(NK(inchain));
        refstatchain(c) = sn.refstat(inchain(1));
        if any((sn.refstat(inchain(1))-refstatchain(c))~=0)
            line_error(mfilename,sprintf('Classes in chain %d have different reference station.',c));
        end
    end
    STchain(~isfinite(STchain))=0;
    Lchain(~isfinite(Lchain))=0;
    Tstart = tic;
    Nt = sum(Nchain(isfinite(Nchain)));

    L = zeros(M,C);
    mu = zeros(M,ceil(Nt));
    infServers = [];
    Z = zeros(M,C);
    for ist=1:M
        if isinf(nservers(ist)) % infinite server
            %mu_chain(i,1:sum(Nchain)) = 1:sum(Nchain);
            infServers(end+1) = ist;
            L(ist,:) = Lchain(ist,:);
            Z(ist,:) = Lchain(ist,:);
            mu(ist,1:Nt) = 1:Nt;
        else
            if strcmpi(options.method,'exact') && nservers(ist)>1
                %options.method = 'default';
                line_warning(mfilename,sprintf('%s does not support exact multiserver yet. Switching to approximate method.\n', 'SolverNC'));
            end
            L(ist,:) = Lchain(ist,:);
            mu(ist,1:Nt) = lldscaling(ist,1:Nt);
        end
    end
    Qchain = zeros(M,C);
    % Solve original system
    [lG,~,method] = pfqn_ncld(L, Nchain, 0*Nchain, mu, options);
    lG = real(lG);
    Xchain=[];

    % Solve systems with a job less
    if isempty(Xchain)
        for r=1:C
            Nchain_r =oner(Nchain,r);
            lGr(r) = pfqn_ncld(L,Nchain_r,0*Nchain,mu,options);
            lGr = real(lGr);
            Xchain(r) = exp(lGr(r) - lG);
            for ist=1:M
                Qchain(ist,r)=0;
            end
            CQchain_r = zeros(M,1);

            if M==2 && any(isinf(sn.nservers)) % repairmen model
                firstDelay = find(isinf(sn.nservers),1);
                Qchain(firstDelay,r) = real(Lchain(firstDelay,r) * Xchain(r));
                Qchain(setdiff(1:M,firstDelay),r) = Nchain(r) - real(Lchain(firstDelay,r) * Xchain(r));
            else
                % Add queue replicas for queue-length
                for ist=1:M
                    Lms_i = L; Lms_i(ist,:) = [];
                    mu_i = mu; mu_i(ist,:) = [];
                    muhati = mu; muhati = pfqn_mushift(mu,ist); %#ok<NASGU>
                    [muhati_f,c] = pfqn_fnc(muhati(ist,:));
                    if Lchain(ist,r)>0
                        if isinf(nservers(ist)) % infinite server
                            Qchain(ist,r) = real(Lchain(ist,r) * Xchain(r));
                        else
                            if  ist==M && sum(isfinite(nservers))==1 % normalize queue-lengths to Nchain(r)
                                Qchain(ist,r) = max(0,real(Nchain(r) - sum(Lchain(isinf(nservers),r)) * Xchain(r)) - sum(Qchain(setdiff(1:(M-1),find(isinf(nservers))),r)));
                            else
                                [lGhat_fnci(r)] = pfqn_ncld([L;L(ist,:)],Nchain_r, 0*Nchain, [muhati;muhati_f], options);
                                [lGhatir(r)] = pfqn_ncld(L,Nchain_r, 0*Nchain, muhati, options);
                                [lGr_i(r)] = pfqn_ncld(Lms_i,Nchain_r, 0*Nchain, mu_i, options);
                                [lGhati(r)] = pfqn_ncld(L,Nchain_r, 0*Nchain, muhati, options);
                                dlGa = real(lGhat_fnci(r)) - real(lGhatir(r));
                                dlG_i = real(lGr_i(r)) - real(lGhatir(r));
                                CQchain(ist) = (exp(dlGa) - 1) + c*(exp(dlG_i)-1); % conditional qlen
                                ldDemand(ist,r) = log(L(ist,r)) + real(lGhati(r)) - log(mu(ist,1)) - real(lGr(r));
                                Qchain(ist,r) = exp(ldDemand(ist,r)) * Xchain(r) * (1+CQchain(ist)); % conditional MVA formula
                            end
                        end
                    end
                end
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
        line_warning(mfilename,'Normalizing constant computations produced a floating-point range exception. Model is likely too large.\n');
    end

    Z = sum(Z(1:M,:),1);

    Rchain = Qchain ./ repmat(Xchain,M,1) ./ Vchain;
    Rchain(infServers,:) = Lchain(infServers,:) ./ Vchain(infServers,:);
    Tchain = repmat(Xchain,M,1) .* Vchain;
    Uchain = Tchain .* Lchain;
    Cchain = Nchain ./ Xchain - Z;

    Xchain=real(Xchain);
    Uchain=real(Uchain);
    Qchain=real(Qchain);
    Rchain=real(Rchain);

    Xchain(~isfinite(Xchain))=0;
    Uchain(~isfinite(Uchain))=0;
    Qchain(~isfinite(Qchain))=0;
    Rchain(~isfinite(Rchain))=0;

    Xchain(Nchain==0)=0;
    Uchain(:,Nchain==0)=0;
    Qchain(:,Nchain==0)=0;
    Rchain(:,Nchain==0)=0;
    Tchain(:,Nchain==0)=0;

    [Q,U,R,T,C,X] = sn_deaggregate_chain_results(sn, Lchain, ST, STchain, Vchain, alpha, [], [], Rchain, Tchain, [], Xchain);

    [ST,gamma,~,~,~,~,eta] = npfqn_nonexp_approx(options.config.highvar,sn,ST0,V,SCV,T,U,gamma,nservers);
end


[lambda,L]= sn_get_product_form_params(sn);
runtime = toc(Tstart);
Q=abs(Q); R=abs(R); X=abs(X); U=abs(U);
for ist=1:M
    if sn.nservers(ist)>1 && sn.nservers(ist)<Inf
        openClasses = find(isinf(NK));
        closedClasses = setdiff(1:K, openClasses);
        for r=closedClasses
            c = find(sn.chains(:,r));
            if X(r) > 0
                U(ist,r) = X(r) * sn.visits{c}(ist,r) / sn.visits{c}(sn.refstat(r),r) * ST(ist,r)/sn.nservers(ist);
            end
        end
        for r=openClasses
            c = find(sn.chains(:,r));
            if lambda(r)>0
                U(ist,r) = lambda(r) * sn.visits{c}(ist,r) / sn.visits{c}(sn.refstat(r),r) * ST(ist,r)/sn.nservers(ist);
            end
        end
    elseif isinf(sn.nservers(ist))
        openClasses = find(isinf(NK));
        closedClasses = setdiff(1:K, openClasses);
        for r=closedClasses
            if X(r) > 0
                c = find(sn.chains(:,r));
                U(ist,r) = X(r) * sn.visits{c}(ist,r) / sn.visits{c}(sn.refstat(r),r) * ST(ist,r);
            end
        end
        for r=openClasses
            if lambda(r)>0
                c = find(sn.chains(:,r));
                U(ist,r) = lambda(r) * sn.visits{c}(ist,r) / sn.visits{c}(sn.refstat(r),r) * ST(ist,r);
            end
        end
    else
        U(ist,:) = U(ist,:) / max(lldscaling(ist,:));
        if sum(U(ist,:)) > 1
            U(ist,:) = U(ist,:) / sum(U(ist,:),"omitnan");
        end
    end
end

X(~isfinite(X))=0; U(~isfinite(U))=0; Q(~isfinite(Q))=0; R(~isfinite(R))=0;

% renormalize qlen and tput to correct for unforeseen population constraint deviations
for c=1:sn.nchains
    inchain = sn.inchain{c};
    Nchain(c) = sum(NK(inchain)); %#ok<FNDSB>
    if isfinite(Nchain(c))
        q_den = sum(sum(Q(:,inchain)));
        if q_den > 0
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