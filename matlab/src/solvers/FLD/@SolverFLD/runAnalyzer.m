function runtime = runAnalyzer(self, options)
% RUNTIME = RUN()
% Run the solver

T0=tic;
iter = NaN;
if nargin<2
    options = self.getOptions;
end
% Wall-clock time-budget launch marker (see options.timeout / lineTimeoutExceeded)
options.timeout_tic = T0;

[pyHandled, options, pyRuntime] = self.runAnalyzerPreamble(options, 'FLD');
if pyHandled
    runtime = pyRuntime;
    return
end

sn = getStruct(self); % this gets modified later on so pass by copy

% Convert non-Markovian distributions to PH
% SSA draws a sample path and the fluid ODEs read mu*phi as a flow, so their
% surrogate must be a genuine phase-type: a matrix exponential has neither.
options.config.phfit = 'ph';
sn = sn_nonmarkov_toph(sn, options);

orig_method = options.method;

% Resolve 'default' BEFORE the feature gate. NetworkSolver.runAnalyzerChecks
% validates getMethodFeatureSet(options.method), and 'default' is not
% 'minnormal', so deferring this to the dispatch switch below would let the
% gate reject a GPS model that the resolved method supports. Only lang='matlab'
% is resolved here: the java path maps 'default' onto JLINE.SolverFluid and
% would warn-and-reset anything else.
if strcmp(options.method,'default') && strcmp(options.lang,'matlab')
    [options.method, mnWhy] = fluid_resolve_default_method(sn, options, self.model);
    if strcmp(options.method,'minnormal')
        line_debug(options, 'FLD default method resolved to: minnormal');
    else
        line_debug(options, 'FLD default resolved to %s, minnormal declined (%s)', ...
            options.method, mnWhy);
    end
end
verboseGuard = self.runAnalyzerChecks(options); %#ok<NASGU> restores the caller verbosity on return

% Finite Capacity Region: the fluid ODEs do not enforce the aggregate
% per-region job limit and would silently return the unconstrained answer.
%
% 'dae' is the exception, and the only one. A region cap is a linear inequality
% on the state and blocking is a throttle on the admission flow that keeps it
% satisfied, so the DAE form has somewhere to put it -- an algebraic equation
% beside the drift -- where an ODE has not. FLUID_CAPACITY_CONSTRAINTS refuses
% the forms that are NOT constraints on this drift (BAS/BBS/RSRD, retrial,
% per-class admission weights) by name. Every other method keeps the blanket
% refusal, because for them it is still true.
if isfield(sn,'nregions') && sn.nregions > 0 && ~any(strcmp(options.method,{'dae','fluid.dae'}))
    line_error(mfilename,'This model uses a Finite Capacity Region (addRegion), which is not supported by SolverFLD. Use options.method=''dae'', or SolverCTMC, SolverJMT, SolverSSA or SolverLDES.');
end

% THE RULES THE REGISTRY CANNOT NAME -- the single-queue shape of 'mfq', the
% single-station shape and horizon of the fluid limits, a Cache node for 'rmf',
% the moment-closure applicability, the server count of 'diffusion', the
% fork-join route and the binding buffer nothing in the FLD tree but 'dae'
% reads (a closed Delay->Queue(cap 2) model once returned 2.19 jobs in a buffer
% of 2) -- live in FLUID_METHOD_REFUSAL, which SUPPORTSMODELMETHOD asks too.
% One predicate, two callers: what the report offers, this run accepts, and
% with enableChecks off this is the only place the rules are still applied.
[okMethod, whyMethod] = fluid_method_refusal(sn, options.method, options, self.model);
if ~okMethod
    line_error(mfilename, whyMethod);
end

% A STOCHASTIC PETRI NET HAS ONE FLUID ROUTE, and it is 'dae'. Every other
% method builds its drift from the station/class/phase encoding, where an
% ordinary Place declares no service process and therefore contributes NO
% coordinate at all: the net would be integrated as an empty model and the table
% would report zeros with no warning. GETMETHODFEATURESET states the same limit
% so the gate refuses it one step earlier, and this names the alternative.
if any(sn.nodetype == NodeType.Transition) && ~any(strcmp(options.method,{'dae','fluid.dae','default'}))
    line_error(mfilename, sprintf(['This model is a stochastic Petri net, which the ''%s'' method has no ' ...
        'drift for. Use options.method=''dae'' (the default for a Petri net), or SolverCTMC, SolverJMT, ' ...
        'SolverSSA or SolverLDES.'], options.method));
end

Solver.resetRandomGeneratorSeed(options.seed);


%options.lang = 'java';

hasOpenClasses = sn_has_open_classes(sn);
switch options.lang
    case {'java'}
        line_debug(options, 'FLD: using lang=java, delegating to JLINE');
        jmodel = LINE2JLINE(self.model);
        %M = jmodel.getNumberOfStatefulNodes;
        M = jmodel.getNumberOfStations;
        R = jmodel.getNumberOfClasses;
        switch options.method
            % 'kp' is forwarded because the JAR implements it (KoPenderAnalyzer);
            % without this arm it fell back to 'default' and the java row silently
            % solved a DIFFERENT model than the matlab one.
            case {'default', 'closing', 'matrix', 'rmf', 'kp'}
                jsolver = JLINE.SolverFluid(jmodel, options);
            otherwise
                line_warning(mfilename,'This solver does not support the specified method. Setting to default.\n');
                options.method = 'default';
                jsolver = JLINE.SolverFluid(jmodel, options);
        end
        % getAvgTable(true) is the UNFILTERED grid, and the reshape below needs it:
        % the no-argument getter DROPS every (station,class) cell whose six metrics
        % are all zero, so on a model with a disabled pair it returns fewer than M*R
        % entries and reshape(...,R,M) errors out. MATLAB applies its own filter when
        % the table is PRINTED, so the bridge must carry the whole grid, zeros included.
        [QN,UN,RN,WN,AN,TN] = JLINE.arrayListToResults(jsolver.getAvgTable(true));
        runtime = toc(T0);
        CN = [];
        XN = [];
        QN = reshape(QN',R,M)';
        UN = reshape(UN',R,M)';
        RN = reshape(RN',R,M)';
        WN = reshape(WN',R,M)';
        AN = reshape(AN',R,M)';
        TN = reshape(TN',R,M)';
        lG = NaN;
        lastiter = NaN;
        self.setAvgResults(QN,UN,RN,TN,AN,WN,CN,XN,runtime,options.method,lastiter);
        self.result.Prob.logNormConstAggr = lG;
        self.result.solverSpecific.sn = JLINE.from_jline_struct(jmodel);
        self.result.solverSpecific.odeStateVec = JLINE.from_jline_matrix(jsolver.result.odeStateVec);
        return
    case 'matlab'
        line_debug(options, 'FLD: using lang=matlab');
        switch options.method
            case {'matrix','fluid.matrix','pnorm','fluid.pnorm'}
                options.method = erase(options.method,'fluid.');
                line_debug(options, 'FLD method: matrix/pnorm requested');
                % Matrix method now supports mixed/open models per Ruuskanen et al., PEVA 151 (2021)
                if sn_has_dps(sn)
                    line_error(mfilename,'The matrix solver does not support DPS scheduling. Use options.method=''closing'' instead.');
                end
            case {'closing','fluid.closing'}
                options.method = 'closing';
                line_debug(options, 'FLD method: closing');
            case {'default'}
                % unreachable for lang='matlab': the resolution above runs
                % before the feature gate. Kept as the single fallback so the
                % two can never disagree.
                options.method = fluid_resolve_default_method(sn, options, self.model);
                line_debug(options, 'FLD default method resolved to: %s', options.method);
            case {'statedep','fluid.statedep','softmin','fluid.softmin'}
                options.method = erase(options.method,'fluid.');
                line_debug(options, 'FLD method: %s', options.method);
                % do nothing
            case {'tbi','fluid.tbi'}
                options.method = 'tbi';
                line_debug(options, 'FLD method: trajectory-based iteration');
                % Trajectory-based iteration (Sheldon, Tuncer, Casale, IEEE T-ITS):
                % cell-decomposed waveform relaxation of the closing ODEs
                if hasOpenClasses
                    line_error(mfilename,'The tbi method supports closed models only.');
                end
                if any(sn.nodetype == NodeType.Cache)
                    line_error(mfilename,'The tbi method does not support caching stations. Use options.method=''rmf'' instead.');
                end
            case {'minnormal','fluid.minnormal','refined','fluid.refined'}
                line_debug(options, 'FLD method: moment closure (%s)', options.method);
                % Open and mixed models are supported by projecting the EXT
                % source pool out of the covariance (see FLUID_MOMENT_TERMS);
                % a multi-phase source is refused there, since those
                % coordinates are the phase of a single arrival process rather
                % than a population. That projection covers 'minnormal' only:
                % FLUID_REFINE_MEANFIELD still solves its correction on
                % orth(D) over the FULL state and would add a perturbation to
                % the source pool mass, a normalisation constant rather than a
                % population. Only minnormal was validated open, so refined
                % keeps the closed-model restriction instead of being declared
                % on an untested path.
                if hasOpenClasses && any(strcmp(options.method,{'refined','fluid.refined'}))
                    line_error(mfilename,'The refined method supports closed models only: its 1/N correction is solved over the full state, including the source pool. Use options.method=''minnormal'' for open or mixed models.');
                end
                % A cache model is answered by the decomposition analyzer with
                % the closure inside its network step (SOLVER_FLUID_ANALYZER),
                % so 'minnormal' is accepted here. 'refined' is not: its 1/N
                % correction is solved over the full state of ONE network, and
                % the cache route has no such state -- the caches are solved in
                % isolation and the network with them relabeled.
                if any(sn.nodetype == NodeType.Cache) && any(strcmp(options.method,{'refined','fluid.refined'}))
                    line_error(mfilename,'The refined method does not support caching stations. Use options.method=''minnormal'' for the moment closure on the queueing layer, or ''rmf''.');
                end
            case {'dae','fluid.dae'}
                line_debug(options, 'FLD method: min-normal closure as a DAE');
                % The DAE route solves the closure as one algebraic system and
                % integrates the transient with a singular mass matrix. It
                % inherits the moment-closure envelope, with two extra limits
                % that are refused here so the message names the model feature
                % rather than surfacing from inside the Newton solve.
                %
                % A CACHE MODEL is a decomposition, not one system: the caches
                % are solved in isolation and the network with them relabeled,
                % so there is no single drift for the constraint to be attached
                % to. SOLVER_FLD_CACHEQN_ANALYZER carries the closure inside its
                % network step for 'minnormal' and 'rmf'; no such route exists
                % for the DAE form.
                if any(sn.nodetype == NodeType.Cache)
                    line_error(mfilename,'The dae method does not support caching stations: a cache model is solved by decomposition, so it has no single drift to constrain. Use options.method=''minnormal'' for the same closure, or ''rmf''.');
                end
                % DPS and GPS close on the covariance BETWEEN class coordinates
                % rather than on the station total, so their closure state is a
                % matrix block. Carrying those blocks simultaneously restores
                % the quartic cost the DAE form exists to avoid; carrying them
                % by substitution restores the alternation it exists to remove.
                if sn_has_dps(sn)
                    line_error(mfilename,'The dae method closes on the per-station variance only, and DPS closes on the covariance between its class coordinates. Use options.method=''minnormal''.');
                end
                % The fork-join route (FLUID_FORKJOIN_ADMITS) is asked above,
                % through FLUID_METHOD_REFUSAL, for every method at once.
            case {'diffusion','fluid.diffusion'}
                line_debug(options, 'FLD method: diffusion approximation');
                % Diffusion approximation - validated in solver_fluid_diffusion
                % do nothing here, validation happens in the solver itself
            case {'mfq','fluid.mfq','butools','aoi','fluid.aoi'}
                % 'butools' names the backend the MFQ branch calls and 'aoi' its
                % age-of-information reading (getAvgAoI/getCdfAoI need
                % method='mfq'); both are aliases rather than separate methods,
                % which is how native python spells them.
                options.method = 'mfq';
                line_debug(options, 'FLD method: Markovian fluid queue');
                % Markovian fluid queue method for single-queue analysis
                % Validation and fallback happens in solver_fluid_analyzer
            case {'kp','fluid.kp'}
                line_debug(options, 'FLD method: Ko-Pender fluid and diffusion limits');
                % Open (MAP_t/Ph_t/inf)^N network; validation in solver_fluid_kp
            case {'rmf','fluid.rmf'}
                line_debug(options, 'FLD method: refined mean field (cache analysis)');
                % Refined mean field method for multi-list cache analysis
                % Uses DDPP framework with 1/N correction
            otherwise
                line_error(mfilename,sprintf('The ''%s'' method is unsupported by this solver.',options.method));
        end
end

% 'mfq' AND 'rmf' ARE LABELLED WITH THE METHOD THAT ANSWERS. Off the shape
% FLUID_MFQ_ADMITS decides, or with no Cache node, their analyzer arms warn and
% re-enter the matrix method (the documented fallback, a tested contract);
% RESOLVEMETHOD tells the gate the same, and the result carries the same name,
% so what the report offered as "runs as 'matrix'" is what the result says.
reportedMethod = '';
if any(strcmp(SolverFLD.canonicalMethod(orig_method), {'mfq','rmf'}))
    reportedMethod = self.resolveMethod(options);
end

% Limited load dependence is honoured only by the closing family (the drift
% multiplies the scheduling share by alpha(n_i) in ODE_RATES_CLOSING_FACTORS).
% Network.getUsedLangFeatures does NOT emit a 'LoadDependence' flag, so no
% feature set can gate this and every other FLD method would silently return
% the answer for alpha == 1. Refuse explicitly instead.
if isfield(sn,'lldscaling') && ~isempty(sn.lldscaling) && any(abs(sn.lldscaling(:)-1) > GlobalConstants.Zero) ...
        && ~any(strcmp(options.method, {'closing','fluid.closing', ...
            'minnormal','fluid.minnormal','refined','fluid.refined'}))
    line_error(mfilename,sprintf(['This model uses load dependence (setLoadDependence), which the ''%s'' method ' ...
        'does not evaluate. Use options.method=''closing'' for the mean-field answer, or ''minnormal''/''refined'' ' ...
        'for the moment-closure correction.'], options.method));
end

if isinf(options.timespan(1))
    if options.verbose  == 2
        line_warning(mfilename,'%s requires options.timespan(1) to be finite. Setting it to 0.\n',mfilename);
    end
    options.timespan(1) = 0;
end

if options.timespan(1) == options.timespan(2)
    line_warning(mfilename,'%s: timespan is a single point, unsupported. Setting options.timespace(1) to 0.\n',mfilename);
    options.timespan(1) = 0;
end

% THE COARSE GATE, and it is method-agnostic: self.supports reads the STATIC
% SolverFLD.getFeatureSet, not getMethodFeatureSet. RUNANALYZERCHECKS above has
% already applied the method-aware gate, which is strictly finer -- it starts
% from the same static set and then narrows it per method. The two can only
% disagree where a method WIDENS the set, and 'dae' is the one that does:
% it declares Region, which the static set must keep false because every other
% method has to go on rejecting a finite capacity region. Deferring to the finer
% verdict for that method is what lets the declaration stand; this check still
% runs, unchanged, for every other method and whenever enableChecks is off
% upstream.
% The single-station fluid limits widen the set too: three of them declare
% Reneging, which the static set must keep false because the network drift
% carries no abandonment flow.
if self.enableChecks && ~any(strcmp(options.method,{'dae','fluid.dae', ...
        'ggisgi.fluid','fluid.ggisgi','ggisgi','ggingi.tga','fluid.tga','tga','tvms','fluid.tvms', ...
        'mtginf','fluid.mtginf','mol','fluid.mol'})) ...
        && ~self.supports(self.model)
    line_error(mfilename,'This model contains features not supported by the solver.');
end

% The resolution of 'default' is PER SOLVE and must not be persisted: the
% non-hyperbolic fallback below fires only for a RESOLVED 'default', so a second
% solve on the same object (getCdfRespT clears the result and re-runs getAvg)
% would request the resolved 'minnormal' EXPLICITLY and die at the Lyapunov step.
persistedOptions = options;
persistedOptions.method = orig_method;
self.setOptions(persistedOptions);

% Fork-join: the same solver-agnostic fixed point MVA and NC drive
% (@NetworkSolver/fjFixedPoint.m), with fldDispatch as the inner solve. The MMT
% transformation emits only Source, Delay, Queue, Router and ClassSwitch, all of
% which the fluid drift already carries, so no fork-join code is added here.
% Intercepted after the feature gate so an unsupported feature is still named by
% its own message rather than by a failure inside the transformed model.
if self.model.hasFork
    fjres = self.fjFixedPoint(options, @(sn_, opt_) self.fldDispatch(sn_, opt_));
    QN = fjres.QN; UN = fjres.UN; RN = fjres.RN; TN = fjres.TN;
    CN = fjres.CN; XN = fjres.XN;
    runtime = fjres.runtime; iter = fjres.iter;
    sn = self.model.getStruct();
    AN = sn_get_arvr_from_tput(sn, TN, self.getAvgTputHandles());
    WN = sn_get_residt_from_respt(sn, RN, self.getAvgResidTHandles());
    fjLabel = fjres.method;
    if ~isempty(reportedMethod)
        fjLabel = reportedMethod;
    end
    self.setAvgResults(QN,UN,RN,TN,AN,WN,CN,XN,runtime,fjLabel,iter);
    if lineTimeoutExceeded(options)
        self.result.Avg.timedOut = true;
        line_warning(mfilename,'Solver stopped after the wall-clock time budget (options.timeout=%gs) was exceeded; returning the interim solution.\n', options.timeout);
    end
    return
end

% The single-station fluid limits are closed forms, not integrations of the
% network drift: they take the whole model in one call and have no initial
% state to average over, so they return here rather than entering the loop
% below. SOLVER_FLUID_QSYS_ANALYZER refuses any model that is not the
% Source -> Queue -> Sink shape they are stated for.
if any(strcmp(options.method, {'ggisgi.fluid','fluid.ggisgi','ggisgi','ggingi.tga','fluid.tga','tga', ...
        'tvms','fluid.tvms','mtginf','fluid.mtginf','mol','fluid.mol'}))
    [QN,UN,RN,TN,CN,XN,Qt,Ut,Tt,runtime,resolved] = solver_fluid_qsys_analyzer(sn, options);
    AN = sn_get_arvr_from_tput(sn, TN, self.getAvgTputHandles());
    self.setAvgResults(QN,UN,RN,TN,AN,[],CN,XN,runtime,resolved,1);
    self.setTranAvgResults(Qt,Ut,{},Tt,{},{},runtime);
    return
end

M = sn.nstations;
K = sn.nclasses;

%%
lastSol= [];
Q = zeros(M,K); R = zeros(M,K); T = zeros(M,K);
U = zeros(M,K); C = zeros(1,K); X = zeros(1,K);
Qt=[];
s0 = sn.space;
s0prior = sn.stateprior;
s0_sz = cellfun(@(x) size(x,1), s0)';
s0_id = pprod(s0_sz-1);
cur_state = sn.state;
while s0_id>=0 % for all possible initial states
    s0prior_val = 1;
    for ind=1:sn.nnodes
        if sn.isstateful(ind)
            isf = sn.nodeToStateful(ind);
            s0prior_val = s0prior_val * s0prior{isf}(1+s0_id(isf)); % update prior
            %sn.state{isf} = s0{isf}(1+s0_id(isf),:); % assign initial state to network
            self.model.nodes{ind}.setState(s0{isf}(1+s0_id(isf),:));
        end
    end
    sn = self.model.getStruct;
    % SSA draws a sample path and the fluid ODEs read mu*phi as a flow, so their
    % surrogate must be a genuine phase-type: a matrix exponential has neither.
    options.config.phfit = 'ph';
    sn = sn_nonmarkov_toph(sn, options);  % Re-apply conversion after fresh struct
    if s0prior_val > 0
        %useJLine = false;
        %if useJLine
        %    [Qfull, Ufull, Rfull, Tfull, Cfull, Xfull, t, Qfull_t, Ufull_t, Tfull_t, lastSol] = JLINE.runFluidAnalyzer(self.model, options);
        %else
        % A non-hyperbolic fluid fixed point (balanced bottlenecks, a saturated
        % multiclass station, an overloaded open station) leaves the linear
        % noise approximation with no stationary covariance. It cannot be seen
        % before the mean is solved, so FLUID_MINNORMAL_APPLICABLE cannot
        % decline it and the moment closure raises at the Lyapunov step.
        %
        % THE LADDER HAS TWO RUNGS, AND THE FIRST ONE KEEPS THE CLOSURE. Most of
        % these failures are not a property of the model at all: SOLVER_FLUID_-
        % MOMENTS must start its alternation at sigma2 = 0, where min(n,c) has no
        % derivative, so a saturated or balanced model's first-order fixed point
        % lands on the kink, sits on a continuum of equilibria, and the Jacobian
        % there is neutral. SOLVER_FLUID_DAE seeds the variance POSITIVE and
        % never adopts sigma2 = 0 as an iterate, so the smoothed E[min(X,c)]
        % breaks the degeneracy and the fixed point is isolated and hyperbolic --
        % it answers the same closure, with a covariance, where the alternation
        % cannot. Dropping straight to first order instead is not merely a lost
        % second moment: on a balanced two-station PS cycle at N=10 it returns
        % [9 1] against the exact [5 5], because a first-order method has no
        % reason to prefer one point of the continuum over another.
        %
        % The second rung is the first-order method, taken when 'dae' declines
        % the model in advance (FLUID_DAE_APPLICABLE) or fails on the same
        % identifier, which is the genuinely non-hyperbolic case: an unstable
        % open station has no stationary distribution to approximate under any
        % closure. Fall back whether 'minnormal' was RESOLVED from 'default' or
        % REQUESTED outright -- the closure has no stationary covariance either
        % way, so refusing an explicit request would only deny the caller the
        % mean still available.
        try
            [Qfull, Ufull, Rfull, Tfull, Cfull, Xfull, t, Qfull_t, Ufull_t, Tfull_t, lastSol, iter, aoiResults] = solver_fluid_analyzer(sn, options);
        catch ME
            if ~strcmp(ME.identifier,'LINE:FluidNonHyperbolic')
                rethrow(ME);
            end
            solved = false;
            % Skip the rung when 'dae' is what just failed, so the ladder cannot
            % run the same solve twice.
            if ~any(strcmp(options.method,{'dae','fluid.dae'}))
                [daeOk, daeReason] = fluid_dae_applicable(sn, options);
                if daeOk
                    daeOptions = options;
                    daeOptions.method = 'dae';
                    try
                        [Qfull, Ufull, Rfull, Tfull, Cfull, Xfull, t, Qfull_t, Ufull_t, Tfull_t, lastSol, iter, aoiResults] = solver_fluid_analyzer(sn, daeOptions);
                        options = daeOptions;
                        solved = true;
                        line_debug(options, 'FLD minnormal declined at the Lyapunov step (%s); falling back to dae', ...
                            ME.message);
                    catch MEdae
                        if ~strcmp(MEdae.identifier,'LINE:FluidNonHyperbolic')
                            rethrow(MEdae);
                        end
                        line_debug(options, 'FLD dae also declined at the Lyapunov step (%s)', MEdae.message);
                    end
                else
                    line_debug(options, 'FLD dae not applicable as a fallback (%s)', daeReason);
                end
            end
            if ~solved
                if sn_has_dps(sn)
                    options.method = 'closing';
                else
                    options.method = 'matrix';
                end
                line_debug(options, 'FLD moment closure declined at the Lyapunov step (%s); falling back to %s', ...
                    ME.message, options.method);
                [Qfull, Ufull, Rfull, Tfull, Cfull, Xfull, t, Qfull_t, Ufull_t, Tfull_t, lastSol, iter, aoiResults] = solver_fluid_analyzer(sn, options);
            end
        end
        %end

        [t,uniqueIdx] = unique(t);
        if isempty(lastSol) % if solution fails
            Q = NaN*ones(M,K); R = NaN*ones(M,K);
            T = NaN*ones(M,K); U = NaN*ones(M,K);
            C = NaN*ones(1,K); X = NaN*ones(1,K);
            Qt = cell(M,K); Ut = cell(M,K); Tt = cell(M,K);
            for ist=1:M
                for r=1:K
                    Qt{ist,r} = [NaN,NaN];
                    Ut{ist,r} = [NaN,NaN];
                    Tt{ist,r} = [NaN,NaN];
                end
            end
        else
            if isempty(self.result) && max(size(Qt))==0 %~exist('Qt','var')
                Q = Qfull*s0prior_val;
                R = Rfull*s0prior_val;
                T = Tfull*s0prior_val;
                U = Ufull*s0prior_val;
                C = Cfull*s0prior_val;
                X = Xfull*s0prior_val;
                Qt = cell(M,K);
                Ut = cell(M,K);
                Tt = cell(M,K);
                for ist=1:M
                    for r=1:K
                        if ~isempty(Qfull_t{ist,r}) && length(Qfull_t{ist,r}) >= max(uniqueIdx)
                            Qfull_t{ist,r} = Qfull_t{ist,r}(uniqueIdx);
                            Ufull_t{ist,r} = Ufull_t{ist,r}(uniqueIdx);
                            Tfull_t{ist,r} = Tfull_t{ist,r}(uniqueIdx);
                            Qt{ist,r} = [Qfull_t{ist,r} * s0prior_val,t];
                            Ut{ist,r} = [Ufull_t{ist,r} * s0prior_val,t];
                            Tt{ist,r} = [Tfull_t{ist,r} * s0prior_val,t];
                        else
                            Qt{ist,r} = [NaN,NaN];
                            Ut{ist,r} = [NaN,NaN];
                            Tt{ist,r} = [NaN,NaN];
                        end
                    end
                end
            else
                Q = Q + Qfull*s0prior_val;
                R = R + Rfull*s0prior_val;
                T = T + Tfull*s0prior_val;
                U = U + Ufull*s0prior_val;
                C = C + Cfull*s0prior_val;
                X = X + Xfull*s0prior_val;
                for ist=1:M
                    for r=1:K
                        [t,uniqueIdx] = unique(t);
                        Qfull_t{ist,r} = Qfull_t{ist,r}(uniqueIdx);
                        Ufull_t{ist,r} = Ufull_t{ist,r}(uniqueIdx);
                        %                                  Tfull_t{i,r} = Tfull_t{i,r}(uniqueIdx);

                        tunion = union(Qt{ist,r}(:,2), t);
                        dataOld = interp1(Qt{ist,r}(:,2),Qt{ist,r}(:,1),tunion);
                        dataNew = interp1(t,Qfull_t{ist,r},tunion);
                        Qt{ist,r} = [dataOld + s0prior_val * dataNew, tunion];

                        dataOld = interp1(Ut{ist,r}(:,2),Ut{ist,r}(:,1),tunion);
                        dataNew = interp1(t,Ufull_t{ist,r},tunion);
                        Ut{ist,r} = [dataOld + s0prior_val * dataNew, tunion];

                        %                                 dataOld = interp1(Tt{i,r}(:,2),Tt{i,r}(:,1),tunion);
                        %                                 dataNew = interp1(t,Tfull_t{i,r},tunion);
                        %                                 Tt{i,r} = [dataOld + s0prior_val * dataNew, tunion];
                    end
                end
            end
        end
    end
    s0_id=pprod(s0_id,s0_sz-1); % update initial state
end
% Now we restore the original state 
for ind=1:sn.nnodes
    if sn.isstateful(ind)
        isf = sn.nodeToStateful(ind);
        self.model.nodes{ind}.setState(cur_state{isf});
    end
end
runtime = toc(T0);
self.result.solverSpecific = lastSol;
% Store AoI results if available
if exist('aoiResults', 'var') && ~isempty(aoiResults)
    self.result.solverSpecific.aoiResults = aoiResults;
end
QN = Q; UN=U; RN=R; TN=T; CN=C; XN=X;

% For cache models, update hit/miss probs in the model. Both fluid routes to
% the decomposition analyzer produce them: 'rmf' with the first-order network
% step, 'minnormal' with the moment closure in its place.
if any(strcmp(options.method, {'rmf','fluid.rmf','minnormal','fluid.minnormal'})) && any(sn.nodetype == NodeType.Cache)
    caches = find(sn.nodetype == NodeType.Cache);
    % Retrieve hitprob/missprob from the cacheqn result stored in lastSol
    if isfield(lastSol, 'cacheHitProb')
        for cIdx = 1:length(caches)
            ind = caches(cIdx);
            self.model.nodes{ind}.setResultHitProb(lastSol.cacheHitProb(cIdx,:));
            self.model.nodes{ind}.setResultMissProb(lastSol.cacheMissProb(cIdx,:));
        end
        self.model.refreshStruct(true);
        sn = self.model.getStruct(true);
    end
end

% Compute average arrival rate at steady-state
AN = sn_get_arvr_from_tput(sn, TN, self.getAvgTputHandles());
if strcmp(orig_method,'default') && ~strcmp(options.method,'default')
    self.setAvgResults(QN,UN,RN,TN,AN,[],CN,XN,runtime,['default/',options.method],iter);
elseif ~isempty(reportedMethod)
    self.setAvgResults(QN,UN,RN,TN,AN,[],CN,XN,runtime,reportedMethod,iter);
else
    self.setAvgResults(QN,UN,RN,TN,AN,[],CN,XN,runtime,options.method,iter);
end
Rt={}; Xt={}; Ct={};
self.setTranAvgResults(Qt,Ut,Rt,Tt,Ct,Xt,runtime);
if lineTimeoutExceeded(options)
    self.result.Avg.timedOut = true;
    line_warning(mfilename,'Solver stopped after the wall-clock time budget (options.timeout=%gs) was exceeded; returning the interim solution.\n', options.timeout);
end
end
