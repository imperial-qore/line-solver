function runtime = runAnalyzer(self, options)
% RUNTIME = RUN()
% Run the solver

T0=tic;

if nargin<2
    options = self.getOptions;
end

% Chain mode: the generator is user-supplied, so there is no state space to
% generate and no performance metric to derive, only the stationary vector.
if self.isChainSolver()
    [pi, infGen, stateSpace, runtime] = solver_ctmc_chain(self.chainModel, options);
    self.result = struct();
    self.result.('solver') = getName(self);
    self.result.infGen = infGen;
    self.result.space = stateSpace;
    self.result.Prob.joint = pi;
    self.result.runtime = runtime;
    return
end
% Wall-clock time-budget launch marker (see options.timeout / lineTimeoutExceeded)
options.timeout_tic = T0;
% Session-level deadline for budget checkpoints in deep utilities that have no
% options argument (e.g. multichoose during state-space generation). Cleared
% automatically when this analyzer returns or errors.
if isfield(options,'timeout') && ~isempty(options.timeout) && isfinite(options.timeout) && options.timeout > 0
    setappdata(0, 'LINEtimeoutDeadline', struct('tic', T0, 'budget', options.timeout));
    timeoutDeadlineCleanup = onCleanup(@() setappdata(0, 'LINEtimeoutDeadline', [])); %#ok<NASGU>
end

[pyHandled, options, pyRuntime] = self.runAnalyzerPreamble(options, 'CTMC');
if pyHandled
    runtime = pyRuntime;
    return
end

% QRF (Quadratic/Linear Reduction) LP-based bounds moved to SolverBA. The text
% lives in unsupportedMethodReason, which runAnalyzerChecks asks BEFORE it
% reports an unlisted method, so a caller carrying an old options.method is
% told where they went rather than only that the method is unsupported. Kept
% here for the enableChecks=false path, which skips that gate entirely.
moved = self.unsupportedMethodReason(options.method);
if ~isempty(moved)
    line_error(mfilename, moved);
end

if ~isinf(options.timespan(1)) && (options.timespan(1) == options.timespan(2))
    line_warning(mfilename,'%s: timespan is a single point, spacing by options.tol (%e).\n',mfilename, options.tol);
    options.timespan(2) = options.timespan(1) + options.tol;
end


verboseGuard = self.runAnalyzerChecks(options); %#ok<NASGU> restores the caller verbosity on return
% Finite Capacity Region: enforced in solver_ctmc.m by filtering the state space
% to states within the aggregate per-region job/memory/linear caps (blocking-
% before-entry). Per-station setCapacity is also honored.
Solver.resetRandomGeneratorSeed(options.seed);


% see _kb/06-solver-catalog.md (CTMC section, feature gate and LN redirect) for rationale

% Inform user about reducible routing handling
if self.enableChecks
    [isErg, ergInfo] = self.model.isRoutingErgodic();
    if ~isErg && ~isempty(ergInfo.absorbingStations)
        absNames = strjoin(ergInfo.absorbingStations, ', ');
        line_printf([...
            'Note: Model has reducible routing with absorbing stations: %s\n' ...
            '  Results represent limiting/absorption probabilities.\n' ...
            '  Use model.getReducibilityInfo() for detailed analysis.\n'], ...
            absNames);
    end
end

sn = getStruct(self);
line_debug(options, 'CTMC: using lang=matlab');

% Chain aggregation, opt-in through options.config.chain_aggregation. The
% state space of a multiclass model grows with the per-class populations, so
% collapsing every chain onto a single class is the standard way to make an
% otherwise intractable model solvable; ModelAdapter.aggregateChains builds the
% collapsed model and sn_deaggregate_chain_results maps its metrics back, both
% of which existed with no solver consumer until this branch. The aggregation is
% EXACT on a product-form model and an approximation otherwise, because one
% aggregate service law replaces the per-class ones weighted by alpha.
% Flow-equivalent server aggregation, opt-in through
% options.config.fes_stations. ModelAdapter.aggregateFES collapses the named
% station subset into one load-dependent station and the collapsed stations'
% own metrics are recovered by conditioning on its population; see
% ctmcFesAggregation. Exact when the subnetwork is product-form.
if isfield(options.config,'fes_stations') && ~isempty(options.config.fes_stations)
    [QN,UN,RN,TN,CN,XN] = ctmcFesAggregation(self, sn, options);
    runtime = toc(T0);
    T = getAvgTputHandles(self);
    AN = sn_get_arvr_from_tput(sn, TN, T);
    self.setAvgResults(QN,UN,RN,TN,AN,[],CN,XN,runtime,[options.method '/fes']);
    return
end

if isfield(options.config,'chain_aggregation') && options.config.chain_aggregation ...
        && sn.nchains < sn.nclasses
    % Driven by @NetworkSolver/transformSolve, so the aggregate is solved by an
    % instance of THIS solver rather than a hard-wired SolverCTMC. Clearing the
    % flag states that the aggregate must not be re-aggregated, rather than
    % relying on its nchains == nclasses guard to decline.
    options.config.chain_aggregation = false;
    options.config.transform = 'chains';
    tr = transformSolve(self, options);
    QN = tr.QN; UN = tr.UN; RN = tr.RN; TN = tr.TN; CN = tr.CN; XN = tr.XN;
    runtime = toc(T0);
    T = getAvgTputHandles(self);
    AN = sn_get_arvr_from_tput(sn, TN, T);
    self.setAvgResults(QN,UN,RN,TN,AN,[],CN,XN,runtime,[options.method '/chainaggr']);
    return
end

% The 'mdd' method never builds the explicit generator, so it returns before
% the state-space path below and leaves result.infGen/space empty by design.
if strcmpi(options.method, 'mdd')
    [QN,UN,RN,TN,CN,XN,mddinfo] = solver_ctmc_mdd_analyzer(sn, options, self.model);
    runtime = toc(T0);
    self.result.mdd = mddinfo;
    T = getAvgTputHandles(self);
    [TN,~,RN] = sn_pn_avg_rates(sn, QN, TN, [], RN);
    AN = sn_get_arvr_from_tput(sn, TN, T);
    self.setAvgResults(QN,UN,RN,TN,AN,[],CN,XN,runtime,options.method);
    if lineTimeoutExceeded(options)
        self.result.Avg.timedOut = true;
    end
    return
end

% Perfect sampling replaces enumeration: intercepted before the state space is
% built, so the memory gate below never applies to it.
if startsWith(lower(options.method), 'cftp')
    line_debug(options, 'CTMC: using perfect sampling (%s), %d samples', options.method, options.samples);
    [QN,UN,RN,TN,CN,XN,Xs,Ts,pAggr,SSq] = solver_ctmc_cftp(sn, options);
    runtime = toc(T0);
    self.result.space = SSq;
    self.result.spaceAggr = SSq;
    self.result.Prob.sampled = pAggr;
    self.result.cftp.samples = Xs;
    self.result.cftp.horizon = Ts;
    T = getAvgTputHandles(self);
    AN = sn_get_arvr_from_tput(sn, TN, T);
    self.setAvgResults(QN,UN,RN,TN,AN,[],CN,XN,runtime,options.method);
    return
end

% Native fork-join support: solve the tag-augmented copy exactly and fold
% the auxiliary sibling classes back into the original classes at the end
isFJ = any(sn.nodetype == NodeType.Fork) || any(sn.nodetype == NodeType.Join);
if isFJ
    % same predicate supportsModelMethod asks, kept for enableChecks=false
    [tranOk, tranWhy] = SolverCTMC.transientSupports(sn, options);
    if ~tranOk
        line_error(mfilename, tranWhy);
    end
    [sn, fjctx] = solver_tr_fjtag_analyzer(self, 'expand', sn, options);
    options.config.state_space_gen = 'reachable';
end

% Convert non-Markovian distributions to PH
sn = sn_nonmarkov_toph(sn, options);
line_debug(options, 'CTMC: converted non-Markovian distributions to PH (nstations=%d, nclasses=%d)', sn.nstations, sn.nclasses);

M = sn.nstations;
K = sn.nclasses;
NK = sn.njobs;
sizeEstimator = 0;
for k=1:K
    sizeEstimator = sizeEstimator + gammaln(1+NK(k)+M-1) - gammaln(1+M-1) - gammaln(1+NK(k)); % worst-case estimate of the state space
end

if any(isinf(sn.njobs))
    if isinf(options.cutoff)
        line_warning(mfilename,sprintf('The model has open chains, it is recommended to specify a finite cutoff value, e.g., SolverCTMC(model,''cutoff'',1).\n'));
        self.options.cutoff= ceil(6000^(1/(M*K)));
        options.cutoff= ceil(6000^(1/(M*K)));
        line_debug(options, 'Open/mixed model: auto-setting cutoff=%d for %d stations, %d classes', options.cutoff, M, K);
        line_warning(mfilename,sprintf('Setting cutoff=%d.\n',self.options.cutoff));
    end
    % Mandatory truncation warning for open/mixed models
    line_printf('CTMC solver using state space cutoff = %d for open/mixed model.\n', options.cutoff);
    line_warning(mfilename,'State space truncation may cause inaccurate results. Consider varying cutoff to assess sensitivity.\n');
end

% see _kb/06-solver-catalog.md (CTMC section, memory pre-gate) for rationale
logNstates = ctmc_state_space_logsize(sn, options);
forceFlag = isfield(options,'force') && ~isempty(options.force) && options.force;
if isfield(options,'memorySafetyFraction') && ~isempty(options.memorySafetyFraction)
    safetyFraction = options.memorySafetyFraction;
else
    safetyFraction = 0.6;
end
line_debug(options, 'State space size estimate: exp(%f)', logNstates);
gateVerbose = isfield(options,'verbose') && options.verbose == VerboseLevel.DEBUG;
[gateOk, gateMsg] = ctmc_memory_gate(logNstates, forceFlag, gateVerbose, safetyFraction);
if ~gateOk
    line_error(mfilename, sprintf('%s Stopping SolverCTMC.\n', gateMsg));
    return
end

% we compute all metrics anyway because CTMC has essentially
% the same cost
if isinf(options.timespan(1))
    line_debug(options, 'Using standard CTMC method for steady-state analysis');
    s0 = sn.state;
    s0prior = sn.stateprior;
    for ind=1:sn.nnodes
        if sn.isstateful(ind)
            isf = sn.nodeToStateful(ind);
            sn.state{isf} = s0{isf}(maxpos(s0prior{1}),:); % pick one particular initial state
        end
    end
    [QN,UN,RN,TN,CN,XN,Q,SS,SSq,Dfilt,~,~,sn,DfiltAux,StartN,PreemptN] = solver_ctmc_analyzer(sn, options);
    if ~isFJ
        % update initial state if this has been corrected by the state space
        % generator (skipped on fork-join models: the analyzed struct is the
        % tag-augmented copy, whose states do not fit the original model)
        for isf=1:sn.nstateful
            ind = sn.statefulToNode(isf);
            self.model.nodes{ind}.setState(sn.state{isf});
            switch class(self.model.nodes{sn.statefulToNode(isf)})
                case 'Cache'
                    self.model.nodes{sn.statefulToNode(isf)}.setResultHitProb(sn.nodeparam{ind}.actualhitprob);
                    self.model.nodes{sn.statefulToNode(isf)}.setResultMissProb(sn.nodeparam{ind}.actualmissprob);
                    if isfield(sn.nodeparam{ind}, 'actualresidt')
                        self.model.nodes{sn.statefulToNode(isf)}.setResultResidT(sn.nodeparam{ind}.actualresidt);
                    end
                    if isfield(sn.nodeparam{ind}, 'actualdelayedhitprob')
                        self.model.nodes{sn.statefulToNode(isf)}.setResultDelayedHitProb(sn.nodeparam{ind}.actualdelayedhitprob);
                    end
                    if isfield(sn.nodeparam{ind}, 'actualitemprob')
                        self.model.nodes{sn.statefulToNode(isf)}.setResultItemProb(sn.nodeparam{ind}.actualitemprob);
                    end
                    if isfield(sn.nodeparam{ind}, 'delayedhitqlen')
                        self.model.nodes{sn.statefulToNode(isf)}.setResultDelayedHitQLen(...
                            sn.nodeparam{ind}.delayedhitqlen, sn.nodeparam{ind}.delayedhitqlenfull);
                    end
                    self.model.refreshChains();
            end
        end
    end
    %sn.space = SS;
    self.result.infGen = Q;
    self.result.space = SS;
    self.result.spaceAggr = SSq;
    self.result.nodeSpace = sn.space;
    self.result.eventFilt = Dfilt;
    % Derived START/PREEMPT filtration and its rates. Kept in their own
    % fields: eventFilt is paired with sn.sync one-to-one and is summed as
    % D1 in @SolverCTMC/sample.m, so a derived filtration that rides on the
    % same arcs must not be appended to it.
    self.result.auxFilt = DfiltAux;
    self.result.startRate = StartN;
    self.result.preemptRate = PreemptN;
    runtime = toc(T0);
    sn.space = {};
    T = getAvgTputHandles(self);
    if isFJ
        [QN,UN,RN,TN,AN,CN,XN] = solver_tr_fjtag_analyzer(self, 'lift', fjctx, QN,UN,RN,TN,CN,XN, T);
    else
        [TN,~,RN] = sn_pn_avg_rates(sn, QN, TN, [], RN);
        AN=sn_get_arvr_from_tput(sn, TN, T);
    end
    self.setAvgResults(QN,UN,RN,TN,AN,[],CN,XN,runtime,options.method);
    if lineTimeoutExceeded(options)
        self.result.Avg.timedOut = true;
        line_warning(mfilename,'Solver exceeded the wall-clock time budget (options.timeout=%gs).\n', options.timeout);
    end
else
    line_debug(options, 'Transient analysis: timespan=[%f,%f]', options.timespan(1), options.timespan(2));
    lastSol= [];
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
                sn.state{isf} = s0{isf}(1+s0_id(isf),:); % assign initial state to network
            end
        end
        if s0prior_val > 0
            [t,pit,QNt,UNt,~,TNt,~,~,Q,SS,SSq,Dfilt,runtime_t] = solver_ctmc_transient_analyzer(sn, options);
            self.result.space = SS;
            self.result.spaceAggr = SSq;
            self.result.infGen = Q;
            self.result.eventFilt = Dfilt;
            %sn.space = SS;
            setTranProb(self,t,pit,SS,runtime_t);
            if isempty(self.result) || ~isfield(self.result,'Tran') || ~isfield(self.result.Tran,'Avg') || ~isfield(self.result.Tran.Avg,'Q')
                self.result.Tran.Avg.Q = cell(M,K);
                self.result.Tran.Avg.U = cell(M,K);
                self.result.Tran.Avg.T = cell(M,K);
                for ist=1:M
                    for r=1:K
                        self.result.Tran.Avg.Q{ist,r} = [QNt{ist,r} * s0prior_val,t];
                        self.result.Tran.Avg.U{ist,r} = [UNt{ist,r} * s0prior_val,t];
                        self.result.Tran.Avg.T{ist,r} = [TNt{ist,r} * s0prior_val,t];
                    end
                end
            else
                for ist=1:M
                    for r=1:K
                        tunion = union(self.result.Tran.Avg.Q{ist,r}(:,2), t);
                        dataOld = interp1(self.result.Tran.Avg.Q{ist,r}(:,2),self.result.Tran.Avg.Q{ist,r}(:,1),tunion);
                        dataNew = interp1(t,QNt{ist,r},tunion);
                        self.result.Tran.Avg.Q{ist,r} = [dataOld+s0prior_val*dataNew,tunion];
                        dataOld = interp1(self.result.Tran.Avg.U{ist,r}(:,2),self.result.Tran.Avg.U{ist,r}(:,1),tunion);
                        dataNew = interp1(t,UNt{ist,r},tunion);
                        self.result.Tran.Avg.U{ist,r} = [dataOld+s0prior_val*dataNew,tunion];

                        dataOld = interp1(self.result.Tran.Avg.T{ist,r}(:,2),self.result.Tran.Avg.T{ist,r}(:,1),tunion);
                        dataNew = interp1(t,TNt{ist,r},tunion);
                        self.result.Tran.Avg.T{ist,r} = [dataOld+s0prior_val*dataNew,tunion];
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
    sn.space = {};
    self.result.('solver') = getName(self);
    self.result.runtime = runtime;
    self.result.solverSpecific = lastSol;
end
end