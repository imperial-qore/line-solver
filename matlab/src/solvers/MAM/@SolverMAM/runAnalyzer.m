function runtime = runAnalyzer(self, options)
% RUNTIME = RUN()
% Run the solver

T0=tic;
if nargin<2
    options = self.getOptions;
end

[pyHandled, options, pyRuntime] = self.runAnalyzerPreamble(options, 'MAM');
if pyHandled
    runtime = pyRuntime;
    return
end

% runAnalyzerChecks subsumes the coarse supports() test; see _kb/06-solver-catalog.md for rationale
verboseGuard = self.runAnalyzerChecks(options); %#ok<NASGU> restores the caller verbosity on return

sn = self.getStruct();
% Finite Capacity Region: MAM does not enforce the aggregate per-region job
% limit and would silently return the unconstrained answer.
if isfield(sn,'nregions') && sn.nregions > 0
    line_error(mfilename,'This model uses a Finite Capacity Region (addRegion), which is not supported by SolverMAM. Use SolverCTMC, SolverJMT, SolverSSA or SolverLDES, or setCapacity for a single-station limit.');
end

Solver.resetRandomGeneratorSeed(options.seed);


%options.lang = 'java';

switch options.lang
    case 'java'
        line_debug(options, 'MAM: using lang=java, delegating to JLINE');
        jmodel = LINE2JLINE(self.model);
        %M = jmodel.getNumberOfStatefulNodes;
        M = jmodel.getNumberOfStations;
        R = jmodel.getNumberOfClasses;
        jsolver = JLINE.SolverMAM(jmodel, options);
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
        return
    case 'matlab'
        line_debug(options, 'MAM: using lang=matlab');
        sn = getStruct(self);

        % Check if transient analysis is requested
        isTran = length(options.timespan) >= 2 && ~isinf(options.timespan(2));

        if isTran && strcmp(options.method, 'ldqbd')
            if mam_transient_qbd_applicable(sn)
                % Correlated MAP arrival/service or non-Poisson arrival: use the
                % Laplace-domain transient QBD solver on the true MAP blocks.
                line_debug(options, 'MAM: transient analysis via Laplace transient QBD');
                [Qt, Ut, Tt] = solver_mam_transient_qbd(sn, options);
            else
                line_debug(options, 'MAM: transient analysis via standard QBD');
                sn = sn_nonmarkov_toph(sn, options);
                [Qt, Ut, Tt] = solver_mam_ldqbd_transient(sn, options);
            end
            Rt = cell(size(Qt));
            Xt = cell(1, sn.nclasses);
            Ct = cell(1, sn.nclasses);
            runtime = toc(T0);
            self.setTranAvgResults(Qt, Ut, Rt, Tt, Ct, Xt, runtime);
        elseif true%~snHasMultipleClosedClasses(sn)
            line_debug(options, 'MAM: calling solver_mam_analyzer (method=%s)', options.method);
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
