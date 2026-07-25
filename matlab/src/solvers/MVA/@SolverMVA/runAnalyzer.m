function runtime = runAnalyzer(self, options)
% RUNTIME = RUNANALYZER(OPTIONS)
% Run the solver

if nargin<2
    options = self.getOptions;
end

% Wall-clock time-budget launch marker (see options.timeout / lineTimeoutExceeded)
options.timeout_tic = tic;

[pyHandled, options, pyRuntime] = self.runAnalyzerPreamble(options, 'MVA');
if pyHandled
    runtime = pyRuntime;
    return
end

self.runAnalyzerChecks(options);

sn = self.getStruct();
% Finite Capacity Region: MVA does not enforce the aggregate per-region job
% limit and would silently return the unconstrained product-form answer.
if isfield(sn,'nregions') && sn.nregions > 0
    line_error(mfilename,'This model uses a Finite Capacity Region (addRegion), which is not supported by SolverMVA. Use SolverJMT, or setCapacity for a single-station limit.');
end
if isfield(sn,'immfeed') && ~isempty(sn.immfeed) && any(sn.immfeed(:))
    line_warning(mfilename,'SolverMVA does not handle immediate feedback (immfeed); the solver will treat self-loops as class-switching with re-queueing.\n');
end

Solver.resetRandomGeneratorSeed(options.seed);


iter = 0;
%options.lang='java';

switch options.lang
    case 'java'
        line_debug(options, 'MVA: using lang=java, delegating to JLINE');
        sn = getStruct(self); % doesn't need initial state
        jmodel = LINE2JLINE(self.model);
        %M = jmodel.getNumberOfStatefulNodes;
        M = jmodel.getNumberOfStations;
        R = jmodel.getNumberOfClasses;
        jsolver = JLINE.SolverMVA(jmodel, options);
        % carry the JAR-side fork-join (MMT) iterate across the rebuilt JLINE solver;
        % see _kb/06-solver-catalog.md (MVA section) on the fork-join fixed point
        if options.config.fj_warmstart && ~isempty(self.fjForkLambda)
            jsolver.setForkWarmStart(JLINE.from_line_matrix(self.fjForkLambda));
        end
        [QN,UN,RN,WN,AN,TN] = JLINE.arrayListToResults(jsolver.getAvgTable);
        if self.model.hasFork
            self.fjForkLambda = JLINE.from_jline_matrix(jsolver.getForkWarmStart());
        end
        runtime = jsolver.result.runtime;
        CN = [];
        XN = [];
        QN = reshape(QN',R,M)';
        UN = reshape(UN',R,M)';
        RN = reshape(RN',R,M)';
        TN = reshape(TN',R,M)';
        WN = reshape(WN',R,M)';
        AN = reshape(AN',R,M)';
        lG = NaN;
        lastiter = NaN;
        for ind = 1:sn.nnodes
            if sn.nodetype(ind) == NodeType.Cache
                jnode = jmodel.getNodeByIndex(ind-1);
                self.model.nodes{ind}.setResultHitProb(JLINE.from_jline_matrix(jnode.getHitRatio()));
                self.model.nodes{ind}.setResultMissProb(JLINE.from_jline_matrix(jnode.getMissRatio()));
                % Retrieval-cache extras (delayed-hit ratio, per-list hit ratio
                % and expected latency) so getAvgCacheTable matches the native path.
                self.model.nodes{ind}.setResultDelayedHitProb(JLINE.from_jline_matrix(jnode.getDelayedHitRatio()));
                self.model.nodes{ind}.setResultHitProbList(JLINE.from_jline_matrix(jnode.getHitRatioByList()));
                self.model.nodes{ind}.setResultItemProb(JLINE.from_jline_matrix(jnode.getItemProb()));
                self.model.nodes{ind}.setResultResidT(JLINE.from_jline_matrix(jnode.getResidT()));
            end
        end
        %self.model.refreshChains();
        self.model.refreshStruct(true);
        self.setAvgResults(QN,UN,RN,TN,AN,WN,CN,XN,runtime,options.method,lastiter);
        self.result.Prob.logNormConstAggr = lG;
        return
    case 'matlab'
        line_debug(options, 'MVA: using lang=matlab');
        % solver-agnostic fork-join fixed point (fjFixedPoint) driven with
        % mvaDispatch as callback; see _kb/06-solver-catalog.md (MVA section)
        fjres = self.fjFixedPoint(options, @(sn_, opt_) self.mvaDispatch(sn_, opt_));
        QN = fjres.QN; UN = fjres.UN; RN = fjres.RN; TN = fjres.TN;
        CN = fjres.CN; XN = fjres.XN; lG = fjres.lG;
        runtime = fjres.runtime; iter = fjres.iter; method = fjres.method;
        actualmethod = fjres.actualmethod;

        sn = self.model.getStruct();

        % Compute average residence time at steady-state
        AN = sn_get_arvr_from_tput(sn, TN, self.getAvgTputHandles());
        WN = sn_get_residt_from_respt(sn, RN, self.getAvgResidTHandles());
        if strcmp(method,'default') && ~isempty(actualmethod)
            self.setAvgResults(QN,UN,RN,TN,AN,WN,CN,XN,runtime,['default/' actualmethod],iter);
        else
            self.setAvgResults(QN,UN,RN,TN,AN,WN,CN,XN,runtime,method,iter);
        end
        self.result.Prob.logNormConstAggr = lG;
        if lineTimeoutExceeded(options)
            self.result.Avg.timedOut = true;
            line_warning(mfilename,'Solver stopped after the wall-clock time budget (options.timeout=%gs) was exceeded; returning the interim solution.\n', options.timeout);
        end
end
end


