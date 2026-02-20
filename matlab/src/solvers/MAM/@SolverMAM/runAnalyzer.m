function runtime = runAnalyzer(self, options)
% RUNTIME = RUN()
% Run the solver

T0=tic;
if nargin<2
    options = self.getOptions;
end

if self.enableChecks && ~self.supports(self.model)
    line_error(mfilename,'This model contains features not supported by the solver.\n');
end

self.runAnalyzerChecks(options);
Solver.resetRandomGeneratorSeed(options.seed);

% Show library attribution if verbose and not yet shown
if options.verbose ~= VerboseLevel.SILENT && ~GlobalConstants.isLibraryAttributionShown()
    sn = self.getStruct();
    libs = SolverMAM.getLibrariesUsed(sn, options);
    if ~isempty(libs)
        line_printf('The solver will leverage %s.\n', strjoin(libs, ', '));
        GlobalConstants.setLibraryAttributionShown(true);
    end
end

%options.lang = 'java';

switch options.lang
    case 'java'
        jmodel = LINE2JLINE(self.model);
        %M = jmodel.getNumberOfStatefulNodes;
        M = jmodel.getNumberOfStations;
        R = jmodel.getNumberOfClasses;
        jsolver = JLINE.SolverMAM(jmodel, options);
        [QN,UN,RN,WN,AN,TN] = JLINE.arrayListToResults(jsolver.getAvgTable);
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
        return
    case 'matlab'
        sn = getStruct(self);
        if true%~snHasMultipleClosedClasses(sn)
            % Call solver_mam_analyzer - percResults is optional 10th output
            [QN,UN,RN,TN,CN,XN,~,actualmethod,iter,percResults] = solver_mam_analyzer(sn, options);
            T = getAvgTputHandles(self);
            AN=sn_get_arvr_from_tput(sn, TN, T);

            runtime=toc(T0);
            if strcmp(options.method,'default')
                self.setAvgResults(QN,UN,RN,TN,AN,[],CN,XN,runtime,['default/' actualmethod],iter);
            else
                self.setAvgResults(QN,UN,RN,TN,AN,[],CN,XN,runtime,options.method,iter);
            end

            % Store percentile results if FJ_codes was used
            if ~isempty(percResults)
                self.result.Percentile = percResults;
            end
        else
            line_warning(mfilename,'SolverMAM supports at most a single closed class.\n');
            runtime=toc(T0);
            self.setAvgResults([],[],[],[],[],[],[],[],runtime,options.method,0);
        end
end
end