classdef EnsembleSolver < Solver
    % Abstract class for solvers applicable to Ensemble models
    %
    % Copyright (c) 2012-2024, Imperial College London
    % All rights reserved.

    properties
        ensemble;
        solvers;
        results;
    end

    methods (Hidden)
        function self = EnsembleSolver(ensmodel, name, options)
            % SELF = ENSEMBLESOLVER(MODEL, NAME, OPTIONS)

            self@Solver(ensmodel, name);
            if nargin>=3 %exist('options','var')
                self.setOptions(options);
            else
                self.setOptions(Solver.defaultOptions);
            end
            %self.ensemble = ensmodel.getEnsemble;
            self.solvers = {};
        end
    end
    
    methods %(Abstract) % implemented with errors for Octave compatibility
        
        function [QN,UN,RT,TT] = getEnsembleAvg(self)
            % [QN,UN,RT,TT] = GETENSEMBLEAVG()

            line_error(mfilename,'An abstract method was called. The function needs to be overridden by a subclass.');

        end

        function init(self) % operations before starting to iterate
            % INIT() % OPERATIONS BEFORE STARTING TO ITERATE

            line_error(mfilename,'An abstract method was called. The function needs to be overridden by a subclass.');

        end

        function pre(self, it) % operations before an iteration
            % PRE(IT) % OPERATIONS BEFORE AN ITERATION

            line_error(mfilename,'An abstract method was called. The function needs to be overridden by a subclass.');

        end

        function [results, runtime] = analyze(self, it, e) % operations within an iteration
            % [RESULTS, RUNTIME] = ANALYZE(IT, E) % OPERATIONS WITHIN AN ITERATION

            line_error(mfilename,'An abstract method was called. The function needs to be overridden by a subclass.');

        end

        function post(self, it) % operations after an iteration
            % POST(IT) % OPERATIONS AFTER AN ITERATION

            line_error(mfilename,'An abstract method was called. The function needs to be overridden by a subclass.');

        end

        function finish(self) % operations after interations are completed
            % FINISH() % OPERATIONS AFTER INTERATIONS ARE COMPLETED
            line_error(mfilename,'An abstract method was called. The function needs to be overridden by a subclass.');

        end

        function bool = converged(self, it) % convergence test at iteration it
            % BOOL = CONVERGED(IT) % CONVERGENCE TEST AT ITERATION IT

            line_error(mfilename,'An abstract method was called. The function needs to be overridden by a subclass.');

        end
    end

    methods % default implementations
        function submodels = list(self, it)
            % SUBMODELS = LIST(IT)

            % submodels to be considered at iteration it
            submodels = 1:self.getNumberOfModels;
        end

        function it = getIteration(self)
            % IT = GETITERATION()

            it = size(results,1);
        end
    end

    methods
        function solver = getSolver(self, e) % solver for ensemble model e
            % SOLVER = GETSOLVER(E) % SOLVER FOR ENSEMBLE MODEL E

            solver = self.solvers{e};
        end

        % setSolver(solvers) : solver cell array is stored as such
        % setSolver(solver) : solver is assigned to all stages
        % setSolver(solver, e) : solver is assigned to stage e
        function solver = setSolver(self, solver, e)
            % SOLVER = SETSOLVER(SOLVER, E)
            %solver.options.verbose = self.options.verbose;
            if iscell(solver)
                self.solvers = solver;
            else
                if nargin<3 %~exist('e','var')
                    for e=1:self.getNumberOfModels
                        self.solvers{e} = solver;
                    end
                else
                    self.solvers{e} = solver;
                end
            end
        end

        function E = getNumberOfModels(self)
            % E = GETNUMBEROFMODELS()

            E = length(self.ensemble);
        end

        function [runtime, sruntime, results] = iterate(self, options)
            % [RUNTIME, SRUNTIME, RESULTS] = ITERATE()
            T0 = tic;
            it = 0;
            options = self.options;
            E = getNumberOfModels(self);
            results = cell(1,E);
            sruntime = zeros(1,E); % solver runtimes
            init(self);
            % nearly identical, but parfor based
            while ~self.converged(it) && it < options.iter_max
                it = it + 1;
                self.pre(it);
                sruntime(it,1:E) = 0;
                T1=tic;
                switch options.method
                    case {'para'}
                        parfor e = self.list(it)
                            [results{it,e}, solverTime] = self.analyze(it,e);
                            sruntime(it,e) = sruntime(it,e) + solverTime;
                        end
                    otherwise
                        for e = self.list(it)
                            [results{it,e}, solverTime] = self.analyze(it,e);
                            sruntime(it,e) = sruntime(it,e) + solverTime;
                        end
                end
                self.results = results;
                if options.verbose
                    Tsolve(it)=toc(T1);
                    Ttot=toc(T0);
                    if it==1
                        line_printf('Iter %2d. ',it);
                    else
                        if it>2
                            line_printf('\nIter %2d. ',it);
                        else
                            line_printf('Iter %2d. ',it);
                        end
                    end
                end
                T2=tic;
                self.post(it);
                Tsynch(it)=toc(T2);
                if options.verbose
                    line_printf('Analyze time: %.3fs. Update time: %.3fs. Runtime: %.3fs. ',Tsolve(it),Tsynch(it),Ttot);
                end
            end
            finish(self);
            runtime = toc(T0);
            if options.verbose
                line_printf('\nSummary: Analyze avg time: %.3fs. Update avg time: %.3fs. Total runtime: %.3fs. ',mean(Tsolve),mean(Tsynch),runtime);
            end
        end

        function AvgTables = getEnsembleAvgTables(self)
            E = getNumberOfModels(self);
            AvgTables = cell(1,E);
            for e=1:E
                AvgTables{1,e} = self.solvers{e}.getAvgTable();
            end
        end
    end

end
