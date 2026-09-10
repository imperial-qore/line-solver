function runtime = runAnalyzer(self, options)
% RUNTIME = RUNANALYZER(OPTIONS)
% Run the solver

if nargin<2
    options = self.options;
end

% Reset random seed
Solver.resetRandomGeneratorSeed(options.seed);

line_debug('AUTO solver starting: method=%s, lang=%s, model=%s', options.method, options.lang, class(self.model));

% The bridges receive AUTO's selection method name, not the delegate's method: the
% constructor resets options.method, so without this every mode would arrive as
% 'default' on the JAR and native-python sides.
if ~strcmp(options.lang,'matlab') && ~strcmp(self.selectionMode,'default')
    options.method = self.selectionMode;
end

% Check if using Java backend
switch options.lang
    case 'python'
        line_debug(options, 'AUTO: using lang=python, delegating to native line_solver');
        [QN,UN,RN,TN,AN,WN,runtime] = PYLINE.getAvg('SolverAuto', self.model, options);
        self.setAvgResults(QN,UN,RN,TN,AN,WN,[],[],runtime,options.method,NaN);
        return
    case 'cpp'
        % The CHOSEN candidate is what runs under lang='cpp'. line-cli's own
        % `-s auto` re-runs its chooser on the C++ side and can pick a different
        % engine than the one this heuristic reported, which would leave the
        % selection this solver announces unattributable to the numbers.
        line_debug(options, 'AUTO: using lang=cpp, delegating the chosen solver to the C++ line-cli');
        chosenSolver = self.chooseSolver('getAvg');
        CPPLINE.solverToken(chosenSolver.getName()); % refuse a wrapper by name
        if ~chosenSolver.supports(self.model)
            line_error(mfilename, sprintf(['AUTO chose %s for this model, which does not ' ...
                'support it; lang=''cpp'' delegates the chosen solver rather than trying ' ...
                'the candidates in turn, so there is no fallback here. Use lang=''matlab''.'], ...
                chosenSolver.getName()));
        end
        cppOptions = chosenSolver.getOptions();
        cppOptions.lang = 'cpp';
        if isfield(options,'arith')
            cppOptions.arith = options.arith;
        end
        runtime = chosenSolver.runAnalyzer(cppOptions);
        self.result = chosenSolver.result;
        self.result.('solver') = self.getName();
        self.result.SelectedSolver = chosenSolver.getName();
        return
    case 'java'
        % Use Java SolverAUTO
        jmodel = LINE2JLINE(self.model);
        M = jmodel.getNumberOfStations;
        R = jmodel.getNumberOfClasses;
        jsolver = JLINE.SolverAuto(jmodel, options);
        % getAvgTable(true) is the UNFILTERED grid, and the reshape below needs it:
        % the no-argument getter DROPS every (station,class) cell whose six metrics
        % are all zero, so on a model with a disabled pair it returns fewer than M*R
        % entries and reshape(...,R,M) errors out. MATLAB applies its own filter when
        % the table is PRINTED, so the bridge must carry the whole grid, zeros included.
        [QN, UN, RN, WN, AN, TN] = JLINE.arrayListToResults(jsolver.getAvgTable(true));
        runtime = jsolver.result.runtime;
        CN = [];
        XN = [];
        QN = reshape(QN', R, M)';
        UN = reshape(UN', R, M)';
        RN = reshape(RN', R, M)';
        TN = reshape(TN', R, M)';
        WN = reshape(WN', R, M)';
        AN = reshape(AN', R, M)';
        self.setAvgResults(QN, UN, RN, TN, AN, WN, CN, XN, runtime, options.method, NaN);
        self.result.SelectedSolver = char(jsolver.result.selectedSolver);
        return
    case 'matlab'
        % Fall through to MATLAB implementation below
end

% Delegate to the appropriate solver
switch class(self.model)
    case 'Network'
        % Select the best solver
        if length(self.candidates) > 1
            chosenSolver = self.chooseSolver('getAvg');
            if chosenSolver.supports(self.model)
                proposedSolvers = {chosenSolver, self.candidates{:}};
            else
                proposedSolvers = self.candidates;
            end
        else
            proposedSolvers = {self.solvers{:}};
        end
        
        % Try each proposed solver
        line_debug('AUTO trying %d candidate solvers', length(proposedSolvers));
        for s = 1:length(proposedSolvers)
            try
                solver = proposedSolvers{s};
                line_debug('AUTO attempting solver: %s', solver.getName());
                runtime = solver.runAnalyzer(options);
                
                % Copy results from the successful solver
                if ~isempty(solver.result)
                    self.result = solver.result;
                    % Update the solver name to reflect AUTO solver
                    self.result.('solver') = self.getName();
                    self.result.SelectedSolver = solver.getName();
                end
                
                if self.options.verbose
                    line_printf('AUTO solver: analysis completed successfully by %s.\n', solver.getName());
                end
                return;
            catch ME
                if self.options.verbose
                    line_printf('AUTO solver: %s failed with error: %s\n', proposedSolvers{s}.getName(), ME.message);
                end
            end
        end
        
        % If we get here, all solvers failed
        line_error(mfilename, 'All candidate solvers failed to analyze the model.');
        
    case 'LayeredNetwork'
        % Similar logic for LayeredNetwork
        line_debug('AUTO handling LayeredNetwork model');
        chosenSolver = self.chooseSolver('getAvg');
        if chosenSolver.supports(self.model)
            proposedSolvers = {chosenSolver, self.candidates{:}};
        else
            proposedSolvers = self.candidates;
        end

        line_debug('AUTO trying %d candidate solvers for LayeredNetwork', length(proposedSolvers));
        for s = 1:length(proposedSolvers)
            try
                solver = proposedSolvers{s};
                line_debug('AUTO attempting solver: %s', solver.getName());
                runtime = solver.runAnalyzer(options);
                
                % Copy results from the successful solver
                if ~isempty(solver.result)
                    self.result = solver.result;
                    % Update the solver name to reflect AUTO solver
                    self.result.('solver') = self.getName();
                    self.result.SelectedSolver = solver.getName();
                end
                
                if self.options.verbose
                    line_printf('AUTO solver: analysis completed successfully by %s.\n', solver.getName());
                end
                return;
            catch ME
                if self.options.verbose
                    line_printf('AUTO solver: %s failed with error: %s\n', proposedSolvers{s}.getName(), ME.message);
                end
            end
        end
        
        % If we get here, all solvers failed
        line_error(mfilename, 'All candidate solvers failed to analyze the model.');
end

end