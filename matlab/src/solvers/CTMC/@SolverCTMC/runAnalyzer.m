function runtime = runAnalyzer(self, options)
% RUNTIME = RUN()
% Run the solver

T0=tic;

if nargin<2
    options = self.getOptions;
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

% QRF (Quadratic/Linear Reduction) LP-based bounds moved to SolverBA.
if startsWith(options.method, 'qrf')
    line_error(mfilename, ['The QRF reduction bounds (method ''%s'') moved to SolverBA. ' ...
        'Use SolverBA(model,''method'',''%s'') (or the ''qr''/''lr'' aliases).'], ...
        options.method, options.method);
end

if ~isinf(options.timespan(1)) && (options.timespan(1) == options.timespan(2))
    line_warning(mfilename,'%s: timespan is a single point, spacing by options.tol (%e).\n',mfilename, options.tol);
    options.timespan(2) = options.timespan(1) + options.tol;
end


self.runAnalyzerChecks(options);
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

% Native fork-join support: solve the tag-augmented copy exactly and fold
% the auxiliary sibling classes back into the original classes at the end
isFJ = any(sn.nodetype == NodeType.Fork) || any(sn.nodetype == NodeType.Join);
if isFJ
    if ~isinf(options.timespan(1))
        line_error(mfilename,'Transient analysis of fork-join models is not supported by SolverCTMC.\n');
    end
    sn_orig = sn;
    Korig = sn.nclasses;
    [~, fjsn, fjclassmap] = ModelAdapter.fjtag(self.model);
    sn = fjsn;
    options.config.state_space_gen = 'reachable';
    line_debug(options, 'CTMC: fork-join tag augmentation, %d classes (%d auxiliary), %d fork firings', sn.nclasses, sn.nclasses-Korig, length(sn.fjsync));
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
logNstates = 0;
nkEff = zeros(1,K);
for k = 1:K
    if isinf(NK(k))
        nk = options.cutoff;
        if numel(nk) > 1; nk = max(nk(:)); end
    else
        nk = NK(k);
    end
    nkEff(k) = nk;
    logNstates = logNstates + gammaln(1+nk+M-1) - gammaln(1+M-1) - gammaln(1+nk);
end
% see _kb/06-solver-catalog.md (CTMC section, memory pre-gate) for rationale
if isfield(sn,'phasessz') && ~isempty(sn.phasessz)
    shareSched = [SchedStrategy.INF, SchedStrategy.PS, SchedStrategy.DPS, ...
        SchedStrategy.GPS, SchedStrategy.PSPRIO, SchedStrategy.DPSPRIO, ...
        SchedStrategy.GPSPRIO, SchedStrategy.LPS];
    for i = 1:min(M, size(sn.phasessz,1))
        for k = 1:min(K, size(sn.phasessz,2))
            p = sn.phasessz(i,k);
            if ~isfinite(p) || p <= 1
                continue
            end
            if sn.sched(i) == SchedStrategy.EXT
                m = 1;
            elseif any(sn.sched(i) == shareSched)
                m = nkEff(k);
            else
                m = min(nkEff(k), sn.nservers(i));
            end
            if ~isfinite(m)
                m = nkEff(k);
            end
            logNstates = logNstates + gammaln(1+m+p-1) - gammaln(1+p-1) - gammaln(1+m);
        end
    end
end
% Routing factor: each (node,class) doing RROBIN/WRROBIN adds a pointer over
% that node's outgoing links.
if isfield(sn,'routing') && ~isempty(sn.routing) && isfield(sn,'connmatrix') && ~isempty(sn.connmatrix)
    for ind = 1:min(size(sn.routing,1), size(sn.connmatrix,1))
        nout = nnz(sn.connmatrix(ind,:));
        if nout <= 1
            continue
        end
        nrr = sum(sn.routing(ind,:) == RoutingStrategy.RROBIN | ...
                  sn.routing(ind,:) == RoutingStrategy.WRROBIN);
        if nrr > 0
            logNstates = logNstates + nrr * log(nout);
        end
    end
end
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
    [QN,UN,RN,TN,CN,XN,Q,SS,SSq,Dfilt,~,~,sn] = solver_ctmc_analyzer(sn, options);
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
    runtime = toc(T0);
    sn.space = {};
    T = getAvgTputHandles(self);
    if isFJ
        [QN,UN,RN,TN,CN,XN] = sn_fj_foldback(QN,UN,RN,TN,CN,XN,fjclassmap,Korig);
        % A Place counts tokens, not firings: rescale before the arrival rates
        % are derived, so that everything downstream sees one convention.
        [TN,~,RN] = sn_pn_avg_rates(sn_orig, QN, TN, [], RN);
        AN=sn_get_arvr_from_tput(sn_orig, TN, T);
        % Join stations report the per-sibling waiting time (JMT convention):
        % QLen over the sibling arrival rate rather than the join firing rate
        for ist=1:sn_orig.nstations
            if sn_orig.nodetype(sn_orig.stationToNode(ist)) == NodeType.Join
                for r=1:Korig
                    if AN(ist,r) > 0
                        RN(ist,r) = QN(ist,r)/AN(ist,r);
                    end
                end
            end
        end
        self.result.fjclassmap = fjclassmap;
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