
function [runtime, tranSysState, tranSync] = runAnalyzer(self, options)
% [RUNTIME, TRANSYSSTATE] = RUN()

T0=tic;
if nargin<2 %~exist('options','var')
    options = self.getOptions;
end
% Wall-clock time-budget launch marker (see options.timeout / lineTimeoutExceeded)
options.timeout_tic = T0;

% options.events (DES event budget) overrides options.samples when set;
% samples remains accepted as a deprecated alias for the event budget.
if isfield(options,'events') && ~isempty(options.events) && isfinite(options.events) && options.events > 0
    options.samples = round(options.events);
end

[pyHandled, options, pyRuntime] = self.runAnalyzerPreamble(options, 'SSA');
if pyHandled
    runtime = pyRuntime;
    return
end

verboseGuard = self.runAnalyzerChecks(options); %#ok<NASGU> restores the caller verbosity on return
Solver.resetRandomGeneratorSeed(options.seed);


% Check if confidence intervals are requested
[confintEnabled, confintLevel] = Solver.parseConfInt(options.confint);

%options.lang = 'java';

sn = getStruct(self);
% see _kb/06-solver-catalog.md for rationale (SSA FCR)

% Native fork-join support: simulate the tag-augmented copy and fold the
% auxiliary sibling classes back into the original classes at the end
isFJ = any(sn.nodetype == NodeType.Fork) || any(sn.nodetype == NodeType.Join);
if isFJ
    if strcmp(options.lang,'java')
        line_warning(mfilename,'Fork-join models are not supported by the JLINE SSA backend, switching to lang=matlab.\n');
        options.lang = 'matlab';
    end
    if strcmp(options.method,'parallel')
        line_warning(mfilename,'The parallel method does not support fork-join models, switching to the serial method.\n');
    end
    options.method = 'serial';
    [sn, fjctx] = solver_tr_fjtag_analyzer(self, 'expand', sn, options);
end

switch options.lang
    case 'java'
        line_debug(options, 'SSA: using lang=java, delegating to JLINE');
        % see _kb/06-solver-catalog.md for rationale (SSA lang=java transient trajectory)
        if nargout > 1
            options.method = 'serial';
        end
        switch options.method
            case {'default','serial','parallel'}
                switch options.method
                    case 'default'
                        options.verbose = VerboseLevel.SILENT;
                        actualmethod = 'parallel';
                end
                jmodel = LINE2JLINE(self.model);
                M = jmodel.getNumberOfStations;
                R = jmodel.getNumberOfClasses;
                tic;
                jsolver = JLINE.SolverSSA(jmodel, options);
                % getAvgTable(true) is the UNFILTERED grid, and the reshape below needs it:
                % the no-argument getter DROPS every (station,class) cell whose six metrics
                % are all zero, so on a model with a disabled pair it returns fewer than M*R
                % entries and reshape(...,R,M) errors out. MATLAB applies its own filter when
                % the table is PRINTED, so the bridge must carry the whole grid, zeros included.
                [QN,UN,RN,WN,AN,TN] = JLINE.arrayListToResults(jsolver.getAvgTable(true));

                % see _kb/06-solver-catalog.md for rationale (SSA lang=java transient trajectory)
                if nargout > 1
                    jres = jsolver.result;
                    if isempty(jres.tranSysState)
                        line_error(mfilename,'SSA transient trajectory unavailable from JLINE (serial analyzer required).');
                    end
                    tranSysState = cell(1, sn.nstateful + 1);
                    tranSysState{1} = JLINE.from_jline_matrix(jres.tranSysState.get(java.lang.Integer(0)));
                    for isf = 1:sn.nstateful
                        tranSysState{1+isf} = JLINE.from_jline_matrix(jres.tranSysState.get(java.lang.Integer(isf)));
                    end
                    tranSync = JLINE.from_jline_matrix(jres.tranSync);
                    tranSync = tranSync(:)';
                    % see _kb/06-solver-catalog.md for rationale (SSA lang=java transient trajectory)
                    runtime = jsolver.result.runtime;
                    return;
                end
                CN = JLINE.from_jline_matrix(jsolver.getAvgSysRespT());
                XN = JLINE.from_jline_matrix(jsolver.getAvgSysTput());
                runtime = jsolver.result.runtime;
                QN = reshape(QN',R,M)';
                UN = reshape(UN',R,M)';
                RN = reshape(RN',R,M)';
                TN = reshape(TN',R,M)';
                WN = reshape(WN',R,M)';
                AN = reshape(AN',R,M)';
                % Extract cache hit/miss probabilities from Java model
                for ind = 1:sn.nnodes
                    if sn.nodetype(ind) == NodeType.Cache
                        hitRatioVec = JLINE.from_jline_matrix(jmodel.getNodeByIndex(ind-1).getHitRatio());
                        missRatioVec = JLINE.from_jline_matrix(jmodel.getNodeByIndex(ind-1).getMissRatio());
                        hitClass = self.model.nodes{ind}.getHitClass;
                        nk = length(hitClass);
                        hitprob = zeros(1, nk);
                        missprob = zeros(1, nk);
                        for k = 1:nk
                            if hitClass(k) > 0 && k <= length(hitRatioVec)
                                hitprob(k) = hitRatioVec(k);
                                missprob(k) = missRatioVec(k);
                            end
                        end
                        self.model.nodes{ind}.setResultHitProb(hitprob);
                        self.model.nodes{ind}.setResultMissProb(missprob);
                    end
                end
                if any(sn.nodetype == NodeType.Cache)
                    self.model.refreshStruct(true);
                end
                self.setAvgResults(QN,UN,RN,TN,AN,WN,CN,XN,runtime);

                % Extract confidence intervals from Java solver results
                if confintEnabled
                    result = jsolver.result;
                    % Helper to check for non-null Java objects
                    isValidMatrix = @(x) ~isempty(x) && isa(x, 'jline.util.matrix.Matrix');
                    if isValidMatrix(result.QNCI)
                        QNCI = JLINE.from_jline_matrix(result.QNCI);
                        QNCI = reshape(QNCI', R, M)';
                    else
                        QNCI = [];
                    end
                    if isValidMatrix(result.UNCI)
                        UNCI = JLINE.from_jline_matrix(result.UNCI);
                        UNCI = reshape(UNCI', R, M)';
                    else
                        UNCI = [];
                    end
                    if isValidMatrix(result.RNCI)
                        RNCI = JLINE.from_jline_matrix(result.RNCI);
                        RNCI = reshape(RNCI', R, M)';
                    else
                        RNCI = [];
                    end
                    if isValidMatrix(result.TNCI)
                        TNCI = JLINE.from_jline_matrix(result.TNCI);
                        TNCI = reshape(TNCI', R, M)';
                    else
                        TNCI = [];
                    end
                    if isValidMatrix(result.ANCI)
                        ANCI = JLINE.from_jline_matrix(result.ANCI);
                        ANCI = reshape(ANCI', R, M)';
                    else
                        ANCI = [];
                    end
                    if isValidMatrix(result.WNCI)
                        WNCI = JLINE.from_jline_matrix(result.WNCI);
                        WNCI = reshape(WNCI', R, M)';
                    else
                        WNCI = [];
                    end
                    % Store CI results
                    self.setAvgResultsCI(QNCI, UNCI, RNCI, TNCI, ANCI, WNCI, [], []);
                end
            otherwise
                line_error(mfilename, ['the ',options.method',' method is not available.']);
        end
    case 'matlab'
        line_debug(options, 'SSA: using lang=matlab');
        [QN,UN,RN,TN,CN,XN,~,actualmethod,tranSysState, tranSync, sn, QNCI, UNCI, RNCI, TNCI, ANCI, WNCI, StartN, PreemptN, tranStartTag, tranPreemptTag] = solver_ssa_analyzer(sn, options);

        for isf=1:sn.nstateful
            ind = sn.statefulToNode(isf);
            switch sn.nodetype(sn.statefulToNode(isf))
                case NodeType.Cache
                    self.model.nodes{sn.statefulToNode(isf)}.setResultHitProb(sn.nodeparam{ind}.actualhitprob);
                    self.model.nodes{sn.statefulToNode(isf)}.setResultMissProb(sn.nodeparam{ind}.actualmissprob);
                    if isfield(sn.nodeparam{ind}, 'actualdelayedhitprob')
                        self.model.nodes{sn.statefulToNode(isf)}.setResultDelayedHitProb(sn.nodeparam{ind}.actualdelayedhitprob);
                    end
                    if isfield(sn.nodeparam{ind}, 'actualresidt')
                        self.model.nodes{sn.statefulToNode(isf)}.setResultResidT(sn.nodeparam{ind}.actualresidt);
                    end
                    self.model.refreshChains();
            end
        end
        line_debug(options, 'SSA analysis complete: extracting results (nstations=%d, nclasses=%d)', sn.nstations, sn.nclasses);
        runtime = toc(T0);
        T = getAvgTputHandles(self);
        if isFJ
            [QN,UN,RN,TN,AN,CN,XN] = solver_tr_fjtag_analyzer(self, 'lift', fjctx, QN,UN,RN,TN,CN,XN, T);
        else
            [TN,~,RN] = sn_pn_avg_rates(sn, QN, TN, [], RN);
            AN = sn_get_arvr_from_tput(sn, TN, T);
        end
        if strcmp(options.method,'default') && exist('actualmethod','var')
            self.setAvgResults(QN,UN,RN,TN,AN,[],CN,XN,runtime,['default/' actualmethod]);
        else
            self.setAvgResults(QN,UN,RN,TN,AN,[],CN,XN,runtime,options.method);
        end
        self.result.space = sn.space;
        % Derived START/PREEMPT rates and the per-step tags of the sampled
        % path. Kept in their own fields: they are annotations on existing
        % transitions, not metrics, so they add no getAvgTable column.
        self.result.startRate = StartN;
        self.result.preemptRate = PreemptN;
        self.result.startTag = tranStartTag;
        self.result.preemptTag = tranPreemptTag;

        % Store CI data if computed
        if confintEnabled && ~isempty(QNCI)
            self.setAvgResultsCI(QNCI, UNCI, RNCI, TNCI, ANCI, WNCI, [], []);
            % HOW LONG THE RUN SHOULD HAVE BEEN, when the caller asked. The
            % batch-means half-width h at confidence 1-alpha over a run of n
            % samples pins the ASYMPTOTIC variance, sigma^2 = (h/z)^2 n, which is
            % the quantity a run length is planned from -- not the stationary
            % variance, which on M/M/1 differs from it by a factor blowing up
            % like (1-rho)^-2. SIM_RUNLENGTH then turns it into the sample count
            % that reaches the requested RELATIVE precision.
            self.result.runLengthPlan = ssa_plan_run_length(options, QN, QNCI, ...
                confintLevel, numel(tranSysState{1}));
        end
        if lineTimeoutExceeded(options)
            self.result.Avg.timedOut = true;
            line_warning(mfilename,'Solver stopped after the wall-clock time budget (options.timeout=%gs) was exceeded; returning the interim solution.\n', options.timeout);
        end
end
end

function plan = ssa_plan_run_length(options, QN, QNCI, confintLevel, nSamples)
% PLAN = SSA_PLAN_RUN_LENGTH(OPTIONS, QN, QNCI, CONFINTLEVEL, NSAMPLES)
%
% The run length the caller would need for the precision they asked for.
% Returns [] unless options.config.runLengthPlan is set; it is either a scalar
% target RELATIVE precision or a struct with fields relprecision and confidence.
%
% See also SIM_RUNLENGTH_PLAN.
plan = [];
if ~isfield(options,'config') || ~isstruct(options.config) ...
        || ~isfield(options.config,'runLengthPlan') || isempty(options.config.runLengthPlan)
    return
end
spec = options.config.runLengthPlan;
relprecision = 0.05;
confidence = confintLevel;
if isstruct(spec)
    if isfield(spec,'relprecision') && ~isempty(spec.relprecision)
        relprecision = spec.relprecision;
    end
    if isfield(spec,'confidence') && ~isempty(spec.confidence)
        confidence = spec.confidence;
    end
elseif isnumeric(spec) && isscalar(spec) && spec > 0
    relprecision = spec;
end
plan = sim_runlength_plan(QN, QNCI, nSamples, 'relprecision', relprecision, ...
    'confidence', confidence);
end
