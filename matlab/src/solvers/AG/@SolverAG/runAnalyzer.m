function runtime = runAnalyzer(self, options)
% RUNTIME = RUNANALYZER(OPTIONS)
% Run the agent-based (RCAT) solver

T0=tic;
if nargin<2
    options = self.getOptions;
end

[pyHandled, options, pyRuntime] = self.runAnalyzerPreamble(options, 'AG');
if pyHandled
    runtime = pyRuntime;
    return
end

% runAnalyzerChecks subsumes the coarse supports() test; see _kb/06-solver-catalog.md for rationale
verboseGuard = self.runAnalyzerChecks(options); %#ok<NASGU> restores the caller verbosity on return

sn = self.getStruct();
% Finite Capacity Region: the RCAT decomposition has no aggregate per-region
% job limit and would silently return the unconstrained answer.
if isfield(sn,'nregions') && sn.nregions > 0
    line_error(mfilename,'This model uses a Finite Capacity Region (addRegion), which is not supported by SolverAG. Use SolverCTMC, SolverJMT, SolverSSA or SolverLDES, or setCapacity for a single-station limit.');
end

Solver.resetRandomGeneratorSeed(options.seed);

switch options.lang
    case 'java'
        line_debug(options, 'AG: using lang=java, delegating to JLINE');
        jmodel = LINE2JLINE(self.model);
        M = jmodel.getNumberOfStations;
        R = jmodel.getNumberOfClasses;
        jsolver = JLINE.SolverAG(jmodel, options);
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
        self.setAvgResults(QN,UN,RN,TN,AN,WN,CN,XN,runtime,options.method,NaN);
        self.result.Prob.logNormConstAggr = NaN;
        return
    case 'matlab'
        line_debug(options, 'AG: calling solver_ag_analyzer (method=%s)', options.method);
        [QN,UN,RN,TN,CN,XN,~,actualmethod,iter,percResults] = solver_ag_analyzer(sn, options);
        T = getAvgTputHandles(self);
        AN=sn_get_arvr_from_tput(sn, TN, T);

        runtime=toc(T0);
        if strcmp(options.method,'default')
            self.setAvgResults(QN,UN,RN,TN,AN,[],CN,XN,runtime,['default/' actualmethod],iter);
        else
            self.setAvgResults(QN,UN,RN,TN,AN,[],CN,XN,runtime,options.method,iter);
        end

        if ~isempty(percResults)
            self.result.Percentile = percResults;
        end
end
end
