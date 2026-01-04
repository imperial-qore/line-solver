classdef EnsembleSolver < Solver
    % EnsembleSolver Abstract base class for ensemble model solvers
    %
    % EnsembleSolver extends the base Solver class to handle ensemble models that
    % consist of multiple interconnected or related submodels. It provides the
    % framework for iterative solution algorithms that decompose complex systems
    % into manageable subproblems and coordinate their solution.
    %
    % @brief Abstract base class for solvers handling ensemble/multi-model systems
    %
    % Key characteristics:
    % - Manages collections of interconnected models
    % - Iterative decomposition and solution algorithms
    % - Convergence checking across multiple submodels
    % - Coordinated analysis of model ensembles
    % - Support for hierarchical and environmental models
    %
    % Ensemble solver workflow:
    % 1. Initialize ensemble and submodel solvers
    % 2. Iterate until convergence:
    %    a. Pre-iteration operations
    %    b. Analyze each submodel
    %    c. Post-iteration operations
    %    d. Check convergence
    % 3. Finalize and compute ensemble metrics
    %
    % Subclasses must implement:
    % - init(): Initialize before iteration
    % - pre()/post(): Pre/post iteration operations
    % - analyze(): Analyze individual submodels
    % - converged(): Check iteration convergence
    % - getEnsembleAvg(): Compute ensemble metrics
    %
    % Copyright (c) 2012-2026, Imperial College London
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
            self.solvers = {};
        end
    end

    methods (Abstract)
        init(self); % operations before starting to iterate
        pre(self, it); % operations before an iteration
        [results, runtime] = analyze(self, it, e); % operations within an iteration
        post(self, it); % operations after an iteration
        finish(self); % operations after ending to iterate
        bool = converged(self, it); % iteration convergence test
        [QN,UN,RN,TN,AN,WN] = getEnsembleAvg(self); % get average metrics for each ensemble model
    end

    methods % default implementations
        function submodels = list(self, it)
            % SUBMODELS = LIST(IT)

            % submodels to be considered at iteration it
            submodels = 1:self.getNumberOfModels;
        end

        function it = getIteration(self)
            % IT = GETITERATION()

            it = size(self.results,1);
        end
    end

    methods
        function solver = getSolver(self, e) % solver for ensemble model e
            % SOLVER = GETSOLVER(E)
            %
            % Return solver for ensemble model E

            solver = self.solvers{e};
        end

        % setSolver(solvers)   : solver cell array is stored as such
        % setSolver(solver)    : solver is assigned to all stages
        % setSolver(solver, e) : solver is assigned to stage e
        function solver = setSolver(self, solver, e)
            % SOLVER = SETSOLVER(SOLVER, E)
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
            %
            % Return number of ensemble models
            E = length(self.ensemble);
        end

        function AvgTables = getEnsembleAvgTables(self)
            E = getNumberOfModels(self);
            AvgTables = cell(1,E);
            for e=1:E
                AvgTables{1,e} = self.solvers{e}.getAvgTable();
            end
        end
        
        % Kotlin-style aliases for get* methods
        function solver = solver(self, e)
            % SOLVER Kotlin-style alias for getSolver
            solver = self.getSolver(e);
        end
        
        function it = iteration(self)
            % ITERATION Kotlin-style alias for getIteration
            it = self.getIteration();
        end
        
        function e = numberOfModels(self)
            % NUMBEROFMODELS Kotlin-style alias for getNumberOfModels
            e = self.getNumberOfModels();
        end
        
        function avg_tables = ensembleAvgTables(self)
            % ENSEMBLEAVGTABLES Kotlin-style alias for getEnsembleAvgTables
            avg_tables = self.getEnsembleAvgTables();
        end

        % Table -> T aliases
        function avg_tables = getEnsembleAvgTs(self)
            % GETENSEMBLEAVGTS Short alias for getEnsembleAvgTables
            avg_tables = self.getEnsembleAvgTables();
        end

        function avg_tables = ensembleAvgTs(self)
            % ENSEMBLEAVGTS Short alias for ensembleAvgTables
            avg_tables = self.ensembleAvgTables();
        end

        function varargout = ensembleAvg(self)
            % ENSEMBLEAVG Kotlin-style alias for getEnsembleAvg
            [varargout{1:nargout}] = self.getEnsembleAvg();
        end
    end

    methods
        [runtime, sruntime, results] = iterate(self, options); % core iteration
    end

end
