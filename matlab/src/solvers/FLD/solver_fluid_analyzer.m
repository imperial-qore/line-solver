function [QN, UN, RN, TN, CN, XN, t, QNt, UNt, TNt, xvec, iter, aoiResults] = solver_fluid_analyzer(sn, options)
% [QN, UN, RN, TN, CN, XN, T, QNT, UNT, TNT, XVEC, iter, aoiResults] = SOLVER_FLUID_ANALYZER(QN, OPTIONS)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
%global GlobalConstants.Immediate
%global GlobalConstants.FineTol

% Initialize aoiResults (will be populated if AoI solver is used)
aoiResults = [];

M = sn.nstations;
K = sn.nclasses;
S = sn.nservers;
SCV = sn.scv;
V = cellsum(sn.visits);
gamma = zeros(M,1);
sched = sn.sched;

line_debug('Fluid analyzer starting: method=%s, nstations=%d, nclasses=%d', options.method, M, K);
phases = sn.phases;
phases_last = sn.phases;
rates0 = sn.rates;

if isempty(options.init_sol)
    options.init_sol = solver_fluid_initsol(sn, options);
end

outer_iters = 1;
outer_runtime = tic;
switch options.method
    case {'matrix','fluid.matrix','default','pnorm','fluid.pnorm'}
        % pnorm uses matrix method with pstar smoothing parameter
        if strcmpi(options.method, 'default')
            line_printf('Default method: using matrix/pnorm fluid method\n');
        end
        line_debug('Using matrix/pnorm method, calling solver_fluid_matrix');
        [QN, UN, RN, TN, xvec_iter, QNt, UNt, TNt, ~, t] = solver_fluid_matrix(sn, options);
    case {'closing','statedep','softmin','fluid.closing','fluid.statedep','fluid.softmin'}
        line_debug('Using closing/statedep method, calling solver_fluid_closing');
        [QN, UN, RN, TN, xvec_iter, QNt, UNt, TNt, ~, t] = solver_fluid_closing(sn, options);
    case {'diffusion','fluid.diffusion'}
        % Diffusion approximation using Euler-Maruyama SDE solver (closed networks only)
        line_debug('Using diffusion method, calling solver_fluid_diffusion');
        [QN, UN, RN, TN, xvec_iter, QNt, UNt, TNt, ~, t] = solver_fluid_diffusion(sn, options);
    case {'mfq','fluid.mfq'}
        % Markovian fluid queue method using BUTools for exact single-queue analysis
        % Also supports Age of Information (AoI) analysis for valid topologies

        % Check for AoI topology first (more specific than general single-queue)
        [isAoI, aoiInfo] = aoi_is_aoi(sn);

        if isAoI
            % Route to AoI solver for Age of Information analysis
            line_debug('AoI topology detected, calling solver_mfq_aoi');
            [QN, UN, RN, TN, xvec_iter, QNt, UNt, TNt, ~, t, ~, ~, aoiResults] = solver_mfq_aoi(sn, options);
        else
            % Check for general single-queue topology
            [isSingleQueue, fluidInfo] = fluid_is_single_queue(sn);

            if isSingleQueue
                line_debug('MFQ topology check passed, calling solver_mfq');
                [QN, UN, RN, TN, xvec_iter, QNt, UNt, TNt, ~, t] = solver_mfq(sn, options);
            else
                % Fallback to matrix method if topology is not suitable
                line_warning(mfilename, 'MFQ not applicable: %s. Falling back to matrix method.', fluidInfo.errorMsg);
                options.method = 'matrix';
                [QN, UN, RN, TN, xvec_iter, QNt, UNt, TNt, ~, t] = solver_fluid_matrix(sn, options);
            end
        end
    otherwise
        line_error(mfilename,sprintf('The ''%s'' method is unsupported by this solver.',options.method));
end
outer_runtime = toc(outer_runtime);


switch options.method
    case {'matrix','closing'}
        % approximate FCFS nodes as state-independent stations
        if any(sched==SchedStrategy.FCFS)
            line_debug('FCFS nodes detected, starting iterative approximation');
            iter = 0;
            eta_1 = zeros(1,M);
            eta = Inf*ones(1,M);
            tol = GlobalConstants.CoarseTol;

            while max(abs(1-eta./eta_1)) > tol & iter <= options.iter_max %#ok<AND2>
                iter = iter + 1;
                eta_1 = eta;
                for ist=1:M
                    sd = rates0(ist,:)>0;
                    UN(ist,sd) = TN(ist,sd) ./ rates0(ist,sd);
                end
                ST0 = 1./rates0;
                ST0(isinf(ST0)) = GlobalConstants.Immediate;
                ST0(isnan(ST0)) = GlobalConstants.FineTol;

                XN = zeros(1,K);
                for k=1:K
                    if sn.refstat(k)>0 % ignore artificial classes
                        XN(k) = TN(sn.refstat(k),k);
                    end
                end
                [ST,gamma,~,~,~,~,eta] = npfqn_nonexp_approx(options.config.highvar,sn,ST0,V,SCV,TN,UN,gamma,S);

                rates = 1./ST;
                rates(isinf(rates)) = GlobalConstants.Immediate;
                rates(isnan(rates)) = GlobalConstants.FineTol; %#ok<NASGU>

                for ist=1:M
                    switch sn.sched(ist)
                        case SchedStrategy.FCFS
                            for k=1:K
                                if rates(ist,k)>0 && SCV(ist,k)>0
                                    [cx,muik,phiik] = Coxian.fitMeanAndSCV(1/rates(ist,k), SCV(ist,k));
                                    % we now handle the case that due to either numerical issues
                                    % or different relationship between scv and mean the size of
                                    % the phase-type representation has changed
                                    phases(ist,k) = length(muik);
                                    if phases(ist,k) ~= phases_last(ist,k) % if number of phases changed
                                        % before we update sn we adjust the initial state
                                        isf = sn.stationToStateful(ist);
                                        [~, nir, sir] = State.toMarginal(sn, ist, sn.state{isf});
                                    end
                                    sn.proc{ist}{k} = cx.getProcess;
                                    sn.mu{ist}{k} = muik;
                                    sn.phi{ist}{k} = phiik;
                                    % For Coxian, jobs always start in phase 1 (entry probability 1 on first phase)
                                    sn.pie{ist}{k} = [1, zeros(1, length(muik)-1)];
                                    sn.phases = phases;
                                    sn.phasessz = max(sn.phases,ones(size(sn.phases)));
                                    sn.phaseshift = [zeros(size(phases,1),1),cumsum(sn.phasessz,2)];                                    
                                    if phases(ist,k) ~= phases_last(ist,k)
                                        isf = sn.stationToStateful(ist);
                                        % we now initialize the new service process
                                        sn.state{isf} = State.fromMarginalAndStarted(sn, ist, nir, sir, options);
                                        sn.state{isf} = sn.state{isf}(1,:); % pick one as the marginals won't change
                                    end
                                end
                            end
                    end

                    options.init_sol = xvec_iter{end}(:);
                    if any(phases_last-phases~=0) % If there is a change of phases reset
                        options.init_sol = solver_fluid_initsol(sn);
                    end
                end
                sn.phases = phases;
                switch options.method
                    case {'matrix'}
                        [~, UN, ~, TN, xvec_iter, ~, ~, ~, ~, ~, inner_iters, inner_runtime] = solver_fluid_matrix(sn, options);
                    case {'closing','statedep'}
                        [~, UN, ~, TN, xvec_iter, ~, ~, ~, ~, ~, inner_iters, inner_runtime] = solver_fluid_closing(sn, options);
                end
                phases_last = phases;
                outer_iters = outer_iters + inner_iters;
                outer_runtime = outer_runtime + inner_runtime;
            end % FCFS iteration ends here
            % The FCFS iteration reinitializes at the solution of the last
            % iterative step. We now have converged in the substitution of the
            % model parameters and we rerun everything from the true initial point
            % so that we get the correct transient.
            options.init_sol = solver_fluid_initsol(sn, options);
            switch options.method
                case {'matrix'}
                    [QN, UN, RN, TN, xvec_iter, QNt, UNt, TNt, ~, t] = solver_fluid_matrix(sn, options);
                case {'closing','statedep'}
                    [QN, UN, RN, TN, xvec_iter, QNt, UNt, TNt, ~, t] = solver_fluid_closing(sn, options);
            end
        end
    case 'statedep'
        % do nothing, a single iteration is sufficient
end

if t(1) == 0
    t(1) = GlobalConstants.FineTol;
end

for ist=1:M
    for k=1:K
        %Qfull_t{i,k} = cumsum(Qfull_t{i,k}.*[0;diff(t)])./t;
        %Ufull_t{i,k} = cumsum(Ufull_t{i,k}.*[0;diff(t)])./t;
    end
end

Ufull0 = UN;
for ist=1:M
    sd = find(QN(ist,:)>0);
    UN(ist,QN(ist,:)==0)=0;
    switch sn.sched(ist)
        case SchedStrategy.INF
            for k=sd
                UN(ist,k) = QN(ist,k);
                UNt{ist,k} = QNt{ist,k};
                TNt{ist,k}  = UNt{ist,k}*sn.rates(ist,k);
            end
        case SchedStrategy.DPS
            %w = sn.schedparam(i,:);
            %wcorr = w(:)*QN(i,:)/(w(sd)*QN(i,sd)');
            for k=sd
                % correct for the real rates, instead of the diffusion
                % approximation rates
                UN(ist,k) = min([1,QN(ist,k)/S(ist),sum(Ufull0(ist,sd)) * (TN(ist,k)./(rates0(ist,k)))/sum(TN(ist,sd)./(rates0(ist,sd)))]);
                TNt{ist,k}  = UNt{ist,k}*sn.rates(ist,k)*sn.nservers(ist); % not sure if this is needed
            end
        otherwise
            for k=sd
                % correct for the real rates, instead of the diffusion
                % approximation rates
                UN(ist,k) = min([1,QN(ist,k)/S(ist),sum(Ufull0(ist,sd)) * (TN(ist,k)./rates0(ist,k))/sum(TN(ist,sd)./rates0(ist,sd))]);
                TNt{ist,k}  = UNt{ist,k}*sn.rates(ist,k)*sn.nservers(ist);
            end
    end
end
UN(isnan(UN))=0;

%switch options.method
%case {'closing','statedep'}
%         for i=1:M
%             if sn.nservers(i) > 0 % not INF
%                 for k = 1:K
%                     UNt{i,k} = min(QNt{i,k} / S(i), QNt{i,k} ./ cellsum({QNt{i,:}}) ); % if not an infinite server then this is a number between 0 and 1
%                     UNt{i,k}(isnan(UNt{i,k})) = 0; % fix cases where qlen is 0
%                 end
%             else % infinite server
%                 for k = 1:K
%                     UNt{i,k} = QNt{i,k};
%                 end
%             end
%         end

for ist=1:M
    sd = find(QN(ist,:)>0);
    RN(ist,QN(ist,:)==0)=0;
    for k=sd
        switch sn.sched(ist)
            case SchedStrategy.INF
                % no-op
            otherwise
                RN(ist,k) = QN(ist,k) / TN(ist,k);
        end
    end
end
RN(isnan(RN))=0;
%end

XN = zeros(1,K);
CN = zeros(1,K);
for k=1:K
    if sn.refstat(k)>0 % ignore artificial classes
        XN(k) = TN(sn.refstat(k),k);
        CN(k) = sn.njobs(k) ./ XN(k);
    end
end
xvec.odeStateVec = xvec_iter{end};
xvec.sn = sn;
iter = outer_iters;
end
  