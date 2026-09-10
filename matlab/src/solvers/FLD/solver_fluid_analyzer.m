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

% A STOCHASTIC PETRI NET TAKES ITS OWN ANALYZER, and returns from here.
%
% The route is the one SOLVER_SSA_ANALYZER takes for the same models: a model
% with any Transition node is a different formalism, so it goes to a dedicated
% runner rather than through the queueing dispatch. It returns early because the
% post-processing below is queueing-specific -- it rewrites UN from sn.rates
% (NaN at a Place) and RN from the station scheduling -- and would overwrite the
% Petri conventions SOLVER_FLUID_PETRI reports (a Place is an INF station whose
% utilization is its token count, and whose throughput is the firing rate of the
% modes consuming from it).
if any(sn.nodetype == NodeType.Transition)
    line_debug('Fluid analyzer: stochastic Petri net, calling solver_fluid_petri');
    [QN, UN, RN, TN, xvec_iter, QNt, UNt, TNt, xvec_t, t, iter, ~, momentResults] = solver_fluid_petri(sn, options);
    XN = zeros(1,K);
    CN = zeros(1,K);
    for k = 1:K
        if sn.refstat(k) > 0
            XN(k) = TN(sn.refstat(k),k);
            if XN(k) > 0
                CN(k) = sn.njobs(k) ./ XN(k);
            end
        end
    end
    if t(1) == 0
        t(1) = GlobalConstants.FineTol;
    end
    xvec = struct();
    if ~isempty(xvec_iter)
        xvec.odeStateVec = xvec_iter{end};
    end
    xvec.sn = sn;
    xvec.moments = momentResults;
    xvec.petriTrajectory = xvec_t;
    return
end

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
else
    % A warm start handed in by the caller (SolverLN carries the previous
    % outer iteration's terminal ODE state) indexes the phase expansion of the
    % model it came from: rebuild it whenever that expansion has moved since.
    expected_init_sol = solver_fluid_initsol(sn, options);
    if numel(options.init_sol) ~= numel(expected_init_sol)
        options.init_sol = expected_init_sol;
    end
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
    case {'closing','statedep','softmin','tbi','fluid.closing','fluid.statedep','fluid.softmin','fluid.tbi'}
        line_debug('Using closing/statedep/tbi method, calling solver_fluid_closing');
        [QN, UN, RN, TN, xvec_iter, QNt, UNt, TNt, ~, t] = solver_fluid_closing(sn, options);
    case {'minnormal','fluid.minnormal','refined','fluid.refined'}
        % Second-order moment closures: linear noise approximation, Gaussian
        % closure of the min() term, and the O(1/N) refined mean field
        if any(sn.nodetype == NodeType.Cache)
            % A cache model is a DECOMPOSITION, not one ODE: the caches are
            % solved in isolation and the network with them relabeled as class
            % switches. The moment closure applies to the queueing layer of
            % that alternation, so the route is the same cacheqn analyzer that
            % 'rmf' takes, with the closure inside NETSOLVE.
            line_debug('Using moment closure on a cache model, calling solver_fld_cacheqn_analyzer');
            [QN, UN, RN, TN, CN, XN, t, QNt, UNt, TNt, xvec_iter, cacheHitProb, cacheMissProb, ~, ~, momentResults] = solver_fld_cacheqn_analyzer(sn, options);
        else
            line_debug('Using moment-closure method %s, calling solver_fluid_moments', options.method);
            [QN, UN, RN, TN, xvec_iter, QNt, UNt, TNt, ~, t, ~, ~, momentResults] = solver_fluid_moments(sn, options);
        end
    case {'dae','fluid.dae'}
        % Same min-normal closure as 'minnormal', stated and solved as one
        % differential-algebraic system instead of by successive substitution:
        % population conservation becomes an equation rather than a consequence
        % of the drift, and the transient carries a time-varying covariance.
        line_debug('Using dae method, calling solver_fluid_dae');
        [QN, UN, RN, TN, xvec_iter, QNt, UNt, TNt, ~, t, ~, ~, momentResults] = solver_fluid_dae(sn, options);
    case {'kp','fluid.kp'}
        % Ko-Pender fluid and diffusion limits of the open (MAP_t/Ph_t/inf)^N
        % network: the mean and the covariance are integrated jointly, so this
        % returns a second moment, as the moment-closure methods above also do
        line_debug('Using kp method, calling solver_fluid_kp');
        [QN, UN, RN, TN, xvec_iter, QNt, UNt, TNt, ~, t] = solver_fluid_kp(sn, options);
    case {'diffusion','fluid.diffusion'}
        % Diffusion approximation using Euler-Maruyama SDE solver (closed networks only)
        line_debug('Using diffusion method, calling solver_fluid_diffusion');
        [QN, UN, RN, TN, xvec_iter, QNt, UNt, TNt, ~, t] = solver_fluid_diffusion(sn, options);
    case {'mfq','fluid.mfq'}
        % Markovian fluid queue method using BUTools for exact single-queue analysis
        % Also supports Age of Information (AoI) analysis for valid topologies

        % THE DOCUMENTED FALLBACK: off the single-queue shape this arm warns and
        % re-enters the matrix method, a tested contract (test_fluid_mfq_mm1's
        % reject_* cases and the sanity goldens). FLUID_MFQ_ADMITS is the ONE
        % predicate deciding the shape: SolverFLD.resolveMethod maps 'mfq' onto
        % 'matrix' on exactly the shapes it refuses, so the gate and the result
        % label say 'matrix' where this arm answers as matrix. The AoI shape is
        % tried first, as before.
        [mfqOk, mfqWhy, isAoI] = fluid_mfq_admits(sn);
        if ~mfqOk
            line_warning(mfilename, 'MFQ not applicable: %s Falling back to matrix method.', mfqWhy);
            options.method = 'matrix';
            [QN, UN, RN, TN, xvec_iter, QNt, UNt, TNt, ~, t] = solver_fluid_matrix(sn, options);
        elseif isAoI
            % Route to AoI solver for Age of Information analysis
            line_debug('AoI topology detected, calling solver_mfq_aoi');
            [QN, UN, RN, TN, xvec_iter, QNt, UNt, TNt, ~, t, ~, ~, aoiResults] = solver_mfq_aoi(sn, options);
        elseif numel(unique(sn.classprio)) > 1
            % Priority classes: use the fluid priority queue
            line_debug('MFQ priority topology, calling solver_mfq_prio');
            [QN, UN, RN, TN, xvec_iter, QNt, UNt, TNt, ~, t] = solver_mfq_prio(sn, options);
        else
            line_debug('MFQ topology check passed, calling solver_mfq');
            [QN, UN, RN, TN, xvec_iter, QNt, UNt, TNt, ~, t] = solver_mfq(sn, options);
        end
    case {'rmf','fluid.rmf'}
        % Refined mean field method for cache analysis
        % Use cacheqn analyzer for integrated cache+queueing network
        % With no Cache node the decomposition has nothing to alternate over and
        % its network step, the matrix method, is the whole answer (a tested
        % contract, the sanity goldens carry it). SolverFLD.resolveMethod maps
        % 'rmf' onto 'matrix' on that shape, so the gate and the label say so.
        if ~any(sn.nodetype == NodeType.Cache)
            line_warning(mfilename, 'RMF not applicable: the model has no Cache node. Falling back to matrix method.');
        end
        line_debug('Using refined mean field method, calling solver_fld_cacheqn_analyzer');
        [QN, UN, RN, TN, CN, XN, t, QNt, UNt, TNt, xvec_iter, cacheHitProb, cacheMissProb, ~, ~] = solver_fld_cacheqn_analyzer(sn, options);
    otherwise
        line_error(mfilename,sprintf('The ''%s'' method is unsupported by this solver.',options.method));
end
outer_runtime = toc(outer_runtime);


switch options.method
    case {'matrix','closing','tbi','minnormal','fluid.minnormal','refined','fluid.refined','dae','fluid.dae'}
        % approximate FCFS nodes as state-independent stations.
        %
        % A CACHE MODEL IS EXCLUDED, as it already is under 'rmf': the solve
        % was a decomposition over a mutated struct (caches relabeled as class
        % switches, routing and visits refreshed), and this loop would re-enter
        % the plain network solver on the ORIGINAL sn, which still holds the
        % cache nodes. The non-exponential refit is therefore not available to
        % cache models under any fluid method, which is a pre-existing limit of
        % the decomposition and not of the closure.
        if any(sched==SchedStrategy.FCFS) && ~any(sn.nodetype == NodeType.Cache)
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
                                    % see _kb/06-solver-catalog.md for rationale
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

                    % see _kb/06-solver-catalog.md for rationale
                    expected_init_sol = solver_fluid_initsol(sn);
                    if ~isempty(xvec_iter) && numel(xvec_iter{end}) == numel(expected_init_sol)
                        options.init_sol = xvec_iter{end}(:);
                    else
                        options.init_sol = expected_init_sol;
                    end
                    if any(phases_last-phases~=0) % If there is a change of phases reset
                        options.init_sol = solver_fluid_initsol(sn);
                    end
                end
                sn.phases = phases;
                switch options.method
                    case {'matrix'}
                        [~, UN, ~, TN, xvec_iter, ~, ~, ~, ~, ~, inner_iters, inner_runtime] = solver_fluid_matrix(sn, options);
                    case {'closing','statedep','tbi'}
                        [~, UN, ~, TN, xvec_iter, ~, ~, ~, ~, ~, inner_iters, inner_runtime] = solver_fluid_closing(sn, options);
                    case {'minnormal','fluid.minnormal','refined','fluid.refined'}
                        [~, UN, ~, TN, xvec_iter, ~, ~, ~, ~, ~, inner_iters, inner_runtime] = solver_fluid_moments(sn, options);
                    % 'dae' IS ADMITTED BY THE OUTER CASE ABOVE, so it has to be
                    % re-solved here too. Without this arm the loop runs the
                    % refit, falls through the switch without assigning
                    % INNER_ITERS, and dies on the next line with 'Unrecognized
                    % function or variable' -- so every model with an FCFS
                    % station was unreachable under options.method='dae'. The
                    % C++, JAR and native Python ports all pass their dae solve
                    % into the same loop (fluid_runner.h fluid_fcfs_nonexp_refit).
                    case {'dae','fluid.dae'}
                        [~, UN, ~, TN, xvec_iter, ~, ~, ~, ~, ~, inner_iters, inner_runtime] = solver_fluid_dae(sn, options);
                end
                phases_last = phases;
                outer_iters = outer_iters + inner_iters;
                outer_runtime = outer_runtime + inner_runtime;
            end % FCFS iteration ends here
            % see _kb/06-solver-catalog.md for rationale
            options.init_sol = solver_fluid_initsol(sn, options);
            switch options.method
                case {'matrix'}
                    [QN, UN, RN, TN, xvec_iter, QNt, UNt, TNt, ~, t] = solver_fluid_matrix(sn, options);
                case {'closing','statedep','tbi'}
                    [QN, UN, RN, TN, xvec_iter, QNt, UNt, TNt, ~, t] = solver_fluid_closing(sn, options);
                case {'minnormal','fluid.minnormal','refined','fluid.refined'}
                    [QN, UN, RN, TN, xvec_iter, QNt, UNt, TNt, ~, t, ~, ~, momentResults] = solver_fluid_moments(sn, options);
                case {'dae','fluid.dae'}
                    % the final solve on the refitted struct, so a dae answer is
                    % read off the SAME phases the loop converged on
                    [QN, UN, RN, TN, xvec_iter, QNt, UNt, TNt, ~, t, ~, ~, momentResults] = solver_fluid_dae(sn, options);
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

% A load-dependent station clears alpha(n) times the nominal work, so its
% utilization normalises by the PEAK scaling, not by the server count:
% Seff = max(c_i, max_n alpha_i(n)), the same T*S/peak convention
% solver_ctmc_avg_from_pi applies (its `ceff`). Without load dependence
% Seff == S and every expression below is unchanged.
Seff = S;
if isfield(sn,'lldscaling') && ~isempty(sn.lldscaling)
    for ist=1:min(M,size(sn.lldscaling,1))
        Seff(ist) = max(S(ist), max(sn.lldscaling(ist,:)));
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
                UN(ist,k) = min([1,QN(ist,k)/Seff(ist),sum(Ufull0(ist,sd)) * (TN(ist,k)./(rates0(ist,k)))/sum(TN(ist,sd)./(rates0(ist,sd)))]);
                TNt{ist,k}  = UNt{ist,k}*sn.rates(ist,k)*sn.nservers(ist); % not sure if this is needed
            end
        otherwise
            for k=sd
                % correct for the real rates, instead of the diffusion
                % approximation rates
                UN(ist,k) = min([1,QN(ist,k)/Seff(ist),sum(Ufull0(ist,sd)) * (TN(ist,k)./rates0(ist,k))/sum(TN(ist,sd)./rates0(ist,sd))]);
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

% A CLASS THE MODEL NEVER ROUTES HERE HAS NO RESPONSE TIME, and QN alone does
% not say so: it holds a remnant of the initial state that the integrator was
% still draining when it stopped, 1.3e-12 on test_CQN_Cox_CS_7. THIS LOOP
% OVERWRITES what the method returned, so the guard the methods carry never
% reached the answer and that test read RN(3,1) = 10.0000086 -- the station's
% own service time -- on picard04 and picard05 against 0 on picard09. See
% FLUID_VISITED_PAIRS for why the visit ratios decide it and a threshold on
% QN or TN cannot.
visited = fluid_visited_pairs(sn, M, K);
for ist=1:M
    served = QN(ist,:)>0 & visited(ist,:);
    sd = find(served);
    RN(ist,~served)=0;
    for k=sd
        switch sn.sched(ist)
            case SchedStrategy.INF
                % no-op
            otherwise
                if TN(ist,k) > GlobalConstants.Zero
                    RN(ist,k) = QN(ist,k) / TN(ist,k);
                else
                    % NO DEPARTURES MEANS NO RESIDENCE TIME TO READ, and a
                    % QN the integrator has not finished draining must not be
                    % divided by a TN of the same origin: the ratio is O(1) and
                    % FILTERMETRIC then rounds the QN away, leaving a queue
                    % length of 0 beside a response time of 10. The methods
                    % guard their own division the same way.
                    RN(ist,k) = 0;
                end
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
if ~isempty(xvec_iter)
    xvec.odeStateVec = xvec_iter{end};
else
    xvec = [];  % Signal failure to caller (runAnalyzer checks isempty)
end
xvec.sn = sn;
if exist('momentResults', 'var')
    xvec.moments = momentResults;
end
if exist('cacheHitProb', 'var')
    xvec.cacheHitProb = cacheHitProb;
    xvec.cacheMissProb = cacheMissProb;
end
iter = outer_iters;
end
