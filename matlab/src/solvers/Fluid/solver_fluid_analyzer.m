function [QN, UN, RN, TN, CN, XN, t, QNt, UNt, TNt, xvec, iter] = solver_fluid_analyzer(sn, options)
% [QN, UN, RN, TN, CN, XN, T, QNT, UNT, TNT, XVEC, iter] = SOLVER_FLUID_ANALYZER(QN, OPTIONS)

% Copyright (c) 2012-2024, Imperial College London
% All rights reserved.
%global GlobalConstants.Immediate
%global GlobalConstants.FineTol

M = sn.nstations;
K = sn.nclasses;
S = sn.nservers;
SCV = sn.scv;
V = cellsum(sn.visits);
gamma = zeros(M,1);
schedid = sn.schedid;
phases = sn.phases;
phases_last = sn.phases;
rates0 = sn.rates;

if isempty(options.init_sol)
    options.init_sol = solver_fluid_initsol(sn, options);
end

outer_iters = 1;
outer_runtime = tic;
switch options.method
    case {'matrix','fluid.matrix','default'}
        [QN, UN, RN, TN, xvec_iter, QNt, UNt, TNt, ~, t] = solver_fluid_matrix(sn, options);
    case {'closing','statedep','softmin','fluid.closing','fluid.statedep','fluid.softmin'}        
        [QN, UN, RN, TN, xvec_iter, QNt, UNt, TNt, ~, t] = solver_fluid_closing(sn, options);
    otherwise
        line_error(mfilename,sprintf('The ''%s'' method is unsupported by this solver.',options.method));
end
outer_runtime = toc(outer_runtime);


switch options.method
    case {'matrix','closing'}
        % approximate FCFS nodes as state-independent stations
        if any(schedid==SchedStrategy.ID_FCFS)
            iter = 0;
            eta_1 = zeros(1,M);
            eta = Inf*ones(1,M);
            tol = GlobalConstants.CoarseTol;

            while max(abs(1-eta./eta_1)) > tol & iter <= options.iter_max %#ok<AND2>
                iter = iter + 1;
                eta_1 = eta;
                for i=1:M
                    sd = rates0(i,:)>0;
                    UN(i,sd) = TN(i,sd) ./ rates0(i,sd);
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

                for i=1:M
                    switch sn.schedid(i)
                        case SchedStrategy.ID_FCFS
                            for k=1:K
                                if rates(i,k)>0 && SCV(i,k)>0
                                    [cx,muik,phiik] = Coxian.fitMeanAndSCV(1/rates(i,k), SCV(i,k));
                                    % we now handle the case that due to either numerical issues
                                    % or different relationship between scv and mean if the size of
                                    % the phase-type representation has changed
                                    phases(i,k) = length(muik);
                                    if phases(i,k) ~= phases_last(i,k) % if number of phases changed
                                        % before we update sn we adjust the initial state
                                        isf = sn.stationToStateful(i);
                                        [~, nir, sir] = State.toMarginal(sn, i, sn.state{isf}, options);
                                    end
                                    sn.proc{i}{k} = cx.getRepres;
                                    sn.mu{i}{k} = muik;
                                    sn.phi{i}{k} = phiik;
                                    sn.phases = phases;
                                    sn.phasessz = max(sn.phases,ones(size(sn.phases)));
                                    sn.phaseshift = [zeros(size(phases,1),1),cumsum(sn.phasessz,2)];                                    
                                    if phases(i,k) ~= phases_last(i,k)
                                        isf = sn.stationToStateful(i);
                                        % we now initialize the new service process
                                        sn.state{isf} = State.fromMarginalAndStarted(sn, i, nir, sir, options);
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

for i=1:M
    for k=1:K
        %Qfull_t{i,k} = cumsum(Qfull_t{i,k}.*[0;diff(t)])./t;
        %Ufull_t{i,k} = cumsum(Ufull_t{i,k}.*[0;diff(t)])./t;
    end
end

Ufull0 = UN;
for i=1:M
    sd = find(QN(i,:)>0);
    UN(i,QN(i,:)==0)=0;
    switch sn.schedid(i)
        case SchedStrategy.ID_INF
            for k=sd
                UN(i,k) = QN(i,k);
                UNt{i,k} = QNt{i,k};
                TNt{i,k}  = UNt{i,k}*sn.rates(i,k);
            end
        case SchedStrategy.ID_DPS
            %w = sn.schedparam(i,:);
            %wcorr = w(:)*QN(i,:)/(w(sd)*QN(i,sd)');
            for k=sd
                % correct for the real rates, instead of the diffusion
                % approximation rates
                UN(i,k) = min([1,QN(i,k)/S(i),sum(Ufull0(i,sd)) * (TN(i,k)./(rates0(i,k)))/sum(TN(i,sd)./(rates0(i,sd)))]);
                TNt{i,k}  = UNt{i,k}*sn.rates(i,k)*sn.nservers(i); % not sure if this is needed
            end
        otherwise
            for k=sd
                % correct for the real rates, instead of the diffusion
                % approximation rates
                UN(i,k) = min([1,QN(i,k)/S(i),sum(Ufull0(i,sd)) * (TN(i,k)./rates0(i,k))/sum(TN(i,sd)./rates0(i,sd))]);
                TNt{i,k}  = UNt{i,k}*sn.rates(i,k)*sn.nservers(i);
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

for i=1:M
    sd = find(QN(i,:)>0);
    RN(i,QN(i,:)==0)=0;
    for k=sd
        switch sn.schedid(i)
            case SchedStrategy.ID_INF
                % no-op
            otherwise
                RN(i,k) = QN(i,k) / TN(i,k);
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
