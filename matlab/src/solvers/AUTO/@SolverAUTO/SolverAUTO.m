classdef SolverAUTO < handle
    % SolverAuto Automatic solver selection based on model characteristics
    %
    % SolverAuto is an intelligent solver that automatically selects the most
    % appropriate solution method based on the characteristics of the queueing
    % network model. It analyzes model properties such as network type, class
    % types, scheduling policies, and size to determine the optimal solver.
    %
    % @brief Intelligent automatic solver selection for queueing network models
    %
    % The solver can choose from multiple candidate solvers including:
    % - MVA (Mean Value Analysis) for product-form networks
    % - NC (Normalizing Constant) for closed networks  
    % - MAM (Matrix Analytic Methods) for non-product-form features
    % - Fluid approximation for large-scale models
    % - LDES simulation for complex models
    % - SSA (Stochastic State-space Analysis) for detailed analysis
    % - CTMC (Continuous Time Markov Chain) for small state spaces
    % - LQNS for layered queueing networks
    %
    % The selection process considers model complexity, solver accuracy,
    % and computational efficiency to provide optimal performance.

    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties (Hidden, Access = public)
        enableChecks;
        % Dispatch mode: 'default'/'heur' heuristic, 'exact',
        % 'sim', 'fast', 'accurate', 'bound'. Held apart from options.method
        % because that field is handed to the delegated solver, which does not
        % know AUTO's selection modes.
        selectionMode = 'default';
    end

    properties (Hidden)
        % Population at or below which an exact solver is preferred over an
        % approximation, when one is available.
        EXACT_POPULATION_MAX = 5;
        % Network solvers
        CANDIDATE_MVA = 1;
        CANDIDATE_NC = 2;
        CANDIDATE_MAM = 3;
        CANDIDATE_FLUID = 4;
        % JMT is not an AUTO candidate: LDES subsumes its feature set, so
        % automatic selection never dispatches to the external simulator.
        % SolverJMT stays reachable through the explicit 'jmt' method name.
        CANDIDATE_SSA = 5;
        CANDIDATE_CTMC = 6;
        CANDIDATE_LDES = 7;
        % LayeredNetwork solvers
        CANDIDATE_LQNS = 1;
        CANDIDATE_LN_NC = 2;
        CANDIDATE_LN_MVA = 3;
        CANDIDATE_LN_MAM = 4;
        CANDIDATE_LN_FLUID = 5;
        % Environment solvers
        CANDIDATE_ENV_MVA = 1;
        CANDIDATE_ENV_NC = 2;
        CANDIDATE_ENV_FLUID = 3;
    end

    properties (Hidden)
        candidates; % feasible solvers
        solvers;
        options;
    end

    properties
        model;
        name;
        % Result of the last analysis, copied from the delegate that ran it.
        % Every other solver carries one and runAnalyzer assigns to it, so
        % without the property every candidate appeared to fail after having
        % run successfully.
        result;
    end

    methods
        % Constructor
        function self = SolverAUTO(model, varargin)
            % SOLVERAUTO Create an automatic solver instance
            %
            % @brief Creates a SolverAUTO instance that automatically selects solvers
            % @param model Network or LayeredNetwork model to be solved
            % @param varargin Optional parameters for solver configuration
            % @return self SolverAuto instance configured for the given model
            self.options = Solver.parseOptions(varargin, Solver.defaultOptions);
            self.model = model;
            self.name = 'SolverAuto';
            if self.options.verbose
                %line_printf('Running LINE version %s',model.getVersion);
            end
            % A method name is either a selection intent or a request for a
            % specific algorithm. Resolving it here keeps the qualified form
            % 'family.submethod' intact, and turns an unknown method name into an
            % error rather than into an empty candidate list.
            [tokenKind, tokenFamily, tokenSubmethod] = self.resolveMethodToken(model, self.options.method);
            if strcmp(tokenKind,'family')
                % The token stays in selectionMode for the lang bridges, which
                % must send the family rather than the delegate's own method.
                self.selectionMode = tokenFamily;
                self.options.method = tokenSubmethod;
                self.solvers{1,1} = self.buildFamilySolver(tokenFamily, model, self.options);
                return
            end
            if SolverAUTO.hasPriorParameters(model)
                % Uncertain parameters expand the model into one instance per
                % alternative; the intent applies to each of them.
                self.selectionMode = tokenFamily;
                inner = self.options;
                inner.verbose = 0;
                self.solvers{1,1} = SolverUQ(model, @(m) LINE(m, inner), self.options);
                return
            end
            if strcmp(tokenFamily,'bound')
                % Bound analysis keeps the method 'auto' because runAnalyzer
                % hands AUTO's options to the delegate, so resetting it to
                % 'default' would silently downgrade the bound to gb.upper.
                self.selectionMode = 'bound';
                self.options.method = 'auto';
                self.solvers{1,1} = SolverBA(model,self.options);
                return
            end
            self.options.method = tokenFamily;
            switch self.options.method
                case {'default','heur','sim','exact','fast','accurate'}
                    % Selection modes: all of them keep the full candidate
                    % list, so chooseSolver can rank per requested metric.
                    self.selectionMode = self.options.method;
                    self.options.method = 'default';
                    %solvers sorted from fastest to slowest
                    self.solvers = {};
                    switch class(model)
                        case 'Network'
                            self.solvers{1,self.CANDIDATE_MAM} = SolverMAM(model);
                            self.solvers{1,self.CANDIDATE_MVA} = SolverMVA(model);
                            self.solvers{1,self.CANDIDATE_NC} = SolverNC(model);
                            self.solvers{1,self.CANDIDATE_FLUID} = SolverFluid(model);
                            self.solvers{1,self.CANDIDATE_SSA} = SolverSSA(model);
                            self.solvers{1,self.CANDIDATE_CTMC} = SolverCTMC(model);
                            self.solvers{1,self.CANDIDATE_LDES} = SolverLDES(model);
                            boolSolver = false(length(self.solvers),1);
                            % A candidate keeps its OWN defaults and receives
                            % only what the caller actually overrode on AUTO.
                            % Handing over AUTO's whole options struct dropped
                            % every solver-specific field the generic defaults
                            % do not carry (CTMC config.state_space_gen and
                            % config.hide_immediate, SSA config.nreplicas, LDES
                            % lang='java') and overwrote the ones they do (CTMC
                            % cutoff 10 -> Inf). CTMC then failed on the missing
                            % field, and AUTO reported that as a warning and
                            % fell through to an approximation.
                            userOptions = SolverAUTO.optionsDelta(self.options, Solver.defaultOptions);
                            userOptions = SolverAUTO.withNamedOptions(userOptions, self.options, varargin);
                            for s=1:length(self.solvers)
                                boolSolver(s) = self.solvers{s}.supports(self.model);
                                self.solvers{s}.setOptions(Solver.mergeOptions(userOptions, self.solvers{s}.getOptions()));
                            end
                            self.candidates = {self.solvers{find(boolSolver)}}; %#ok<FNDSB>
                        case 'LayeredNetwork'
                            self.solvers{1,self.CANDIDATE_LQNS} = SolverLQNS(model,self.options);
                            self.solvers{1,self.CANDIDATE_LN_NC} = SolverLN(model,@(m) SolverNC(m,'verbose', false),self.options);
                            self.solvers{1,self.CANDIDATE_LN_MVA} = SolverLN(model,@(m) SolverMVA(m,'verbose', false),self.options);
                            self.solvers{1,self.CANDIDATE_LN_MAM} = SolverLN(model,@(m) SolverMAM(m,'verbose', false),self.options);
                            self.solvers{1,self.CANDIDATE_LN_FLUID} = SolverLN(model,@(m) SolverFluid(m,'verbose', false),self.options);
                            self.candidates = self.solvers;
                        case 'Environment'
                            self.solvers{1,self.CANDIDATE_ENV_MVA} = SolverENV(model,@(m) SolverMVA(m,'verbose', false),self.options);
                            self.solvers{1,self.CANDIDATE_ENV_NC} = SolverENV(model,@(m) SolverNC(m,'verbose', false),self.options);
                            self.solvers{1,self.CANDIDATE_ENV_FLUID} = SolverENV(model,@(m) SolverFluid(m,'verbose', false),self.options);
                            self.candidates = self.solvers;
                    end
            end
            %turn off warnings temporarily
            wstatus = warning('query');
            warning off;
            warning(wstatus);
        end

        function sn = getStruct(self)
            % QN = GETSTRUCT()

            % Get data structure summarizing the model
            sn = self.model.getStruct(true);
        end

        function out = getName(self)
            % OUT = GETNAME()
            % Get solver name
            out = self.name;
        end

        function model = getModel(self)
            % MODEL = GETMODEL()
            % Get the model being solved
            model = self.model;
        end

        function results = getResults(self)
            % RESULTS = GETRESULTS()
            % Return results data structure from chosen solver
            results = self.delegate('getResults', 1);
        end

        function bool = hasResults(self)
            % BOOL = HASRESULTS()
            % Check if the solver, or any delegate it dispatched to, has
            % computed results.
            bool = ~isempty(self.result);
            if bool
                return
            end
            for s = 1:length(self.solvers)
                if ~isempty(self.solvers{s})
                    try
                        if self.solvers{s}.hasResults()
                            bool = true;
                            return
                        end
                    catch
                        % A delegate without results contributes nothing.
                    end
                end
            end
        end

        function options = getOptions(self)
            % OPTIONS = GETOPTIONS()
            % Return options data structure
            options = self.options;
        end

        function self = setOptions(self, options)
            % SELF = SETOPTIONS(OPTIONS)
            % Set a new options data structure for all candidate solvers
            self.options = options;
            for s = 1:length(self.solvers)
                if ~isempty(self.solvers{s})
                    self.solvers{s}.setOptions(options);
                end
            end
        end

        function self = setAvgResults(self,Q,U,R,T,A,W,C,X,runtime,method,iter)
            % SELF = SETAVGRESULTS(SELF,Q,U,R,T,A,W,C,X,RUNTIME,METHOD,ITER)
            % Store average metrics at steady-state. The lang bridges return
            % arrays rather than a delegate, so they need this entry point.
            if nargin < 11
                method = self.options.method;
            end
            if nargin < 12
                iter = NaN;
            end
            self.result.('solver') = self.getName();
            self.result.Avg.('method') = method;
            self.result.Avg.('iter') = iter;
            if isnan(Q), Q=[]; end
            if isnan(R), R=[]; end
            if isnan(T), T=[]; end
            if isnan(U), U=[]; end
            if isnan(X), X=[]; end
            if isnan(C), C=[]; end
            if isnan(A), A=[]; end
            if isnan(W), W=[]; end
            if isa(self.model,'Network')
                [Q,U] = NetworkSolver.zeroSourceMetrics(self.model.getStruct(),Q,U);
            end
            self.result.Avg.Q = real(Q);
            self.result.Avg.R = real(R);
            self.result.Avg.X = real(X);
            self.result.Avg.U = real(U);
            self.result.Avg.T = real(T);
            self.result.Avg.C = real(C);
            self.result.Avg.A = real(A);
            self.result.Avg.W = real(W);
            self.result.Avg.runtime = runtime;
            if ~isfield(self.result.Avg,'timedOut')
                self.result.Avg.timedOut = false;
            end
        end

        function [name, solver] = getSelectedMethod(self)
            % [NAME, SOLVER] = GETSELECTEDMETHOD()
            % Name of the method family that answered the last request, and the
            % delegate that ran it. Before any run, the family that would be
            % chosen for the mean metrics is reported instead.
            solver = self.getSelectedSolver();
            if isempty(solver)
                name = '';
            else
                name = SolverAUTO.familyOfSolver(solver);
            end
        end

        function solver = getSelectedSolver(self)
            % SOLVER = GETSELECTEDSOLVER()
            % The delegate that holds results, or the one that the selection
            % rule would use for the mean performance metrics.
            solver = [];
            for s = 1:length(self.solvers)
                if ~isempty(self.solvers{s}) && self.solvers{s}.hasResults()
                    solver = self.solvers{s};
                    return
                end
            end
            if numel(self.solvers) == 1
                solver = self.solvers{1};
                return
            end
            try
                solver = self.chooseSolver('getAvg');
            catch
                % No feasible delegate: the caller learns this from the
                % analysis itself, not from an introspection call.
                solver = [];
            end
        end

        function bool = isStochastic(self)
            % BOOL = ISSTOCHASTIC()
            % If a delegate solver has already produced results, classify
            % by the solver that actually ran (which knows the method it
            % resolved at runtime). Before any run, the delegate choice is
            % unknown, so classify conservatively: true if any feasible
            % candidate is stochastic.
            bool = false;
            ran = false;
            for s = 1:length(self.solvers)
                if ~isempty(self.solvers{s}) && self.solvers{s}.hasResults()
                    ran = true;
                    if self.solvers{s}.isStochastic()
                        bool = true;
                        return;
                    end
                end
            end
            if ~ran
                if ~isempty(self.candidates)
                    list = self.candidates;
                else
                    list = self.solvers;
                end
                for s = 1:length(list)
                    if ~isempty(list{s}) && list{s}.isStochastic()
                        bool = true;
                        return;
                    end
                end
            end
        end

        function allMethods = listAllMethods(self)
            % ALLMETHODS = LISTALLMETHODS()
            % Every token the CONSTRUCTOR accepts, INDEPENDENT of the model:
            % the selection intents, the method families, and every family
            % method in its qualified form. The unqualified form of a family
            % method is accepted too, and is left out here only to keep the
            % list unambiguous.
            %
            % THIS is the list a method-NAME check must gate on.
            % LISTVALIDMETHODS narrows it to the model in hand, and gating a
            % name check on that would replace a rejection the delegate would
            % have EXPLAINED with a flat "the method is unsupported by this
            % solver" -- the same distinction SolverBA draws between its own
            % listAllMethods and listValidMethods.
            allMethods = [SolverAUTO.selectionIntents(), {'auto'}]; % 'ai' not yet available
            families = SolverAUTO.familyNames();
            probeOptions = self.options;
            probeOptions.verbose = 0;
            probeOptions.method = 'default';
            for f = 1:length(families)
                try
                    probe = self.buildFamilySolver(families{f}, self.model, probeOptions);
                    if ismethod(probe, 'listAllMethods')
                        declared = probe.listAllMethods();
                    else
                        declared = probe.listValidMethods();
                    end
                catch
                    % A family that cannot even be INSTANTIATED here (SolverLQNS
                    % without the lqns binary, SolverLN on a flat Network)
                    % contributes nothing.
                    continue
                end
                allMethods{end+1} = families{f}; %#ok<AGROW>
                for m = 1:length(declared)
                    allMethods{end+1} = [families{f},'.',declared{m}]; %#ok<AGROW>
                end
            end
            allMethods = unique(allMethods);
        end

        function allMethods = listValidMethods(self)
            % ALLMETHODS = LISTVALIDMETHODS()
            % The method names of LISTALLMETHODS that THIS MODEL can actually run.
            %
            % IT IS THE RUNNABLE ROWS OF FINDSOLVER, projected onto their
            % Method column. The narrowing used to be written out a second
            % time here, and a second copy of one gate is how two answers to
            % one question start to differ; FINDSOLVER owns it now, and this
            % adds only the method names that name no single method: the selection
            % intents, which name a RANKING rather than an algorithm, and each
            % family's bare method name.
            %
            % The gate FINDSOLVER applies is the one CHOOSESOLVERRANKED
            % applies before delegating, asked of every candidate instead of
            % the first feasible one: a family whose FEATURE SET refuses the
            % model contributes nothing, and a method whose SUPPORTSMODELMETHOD
            % refuses the model is not offered. Both are needed, because the
            % rules a flat feature set cannot express (product form for
            % 'exact', a binding finite buffer, NC 'mem' applicability) live
            % only in the method-level gate.
            %
            % Without them this returned the token universe regardless of the
            % model: 190 entries for the two-station BAS-blocking model of
            % cqn_bas_blocking, naming all 37 SolverNC methods although
            % SolverNC refuses that model method by method, and every
            % SolverBA method although SolverBA's feature set refuses it
            % outright. A caller enumerating the list was being invited to ask
            % for an analysis no candidate would perform.
            %
            % A family whose every method is refused loses its bare method name too:
            % 'nc' alone delegates to SolverNC, which is exactly the rejection
            % the per-method gate just returned. That falls out of the
            % projection, since a family with no runnable row contributes no
            % row to take its name from.
            allMethods = [SolverAUTO.selectionIntents(), {'auto'}]; % 'ai' not yet available
            T = self.findSolver();
            for k = 1:height(T)
                allMethods{end+1} = T.Method{k}; %#ok<AGROW>
                allMethods{end+1} = T.Solver{k}; %#ok<AGROW>
            end
            allMethods = unique(allMethods);
        end

    end

    methods
        % findSolver: which solvers and methods can analyze this model
        T = findSolver(self, metric, showAll);

        % delegate execution of method to chosen solver
        varargout = delegate(self, method, nretout, varargin);

        % chooseSolver: choses a solver from static properties of the model
        solver = chooseSolver(self, method);
        % Heuristic choice of solver
        solver = chooseSolverHeur(self, method);
        % Heuristic choice of solver for Avg* methods
        solver = chooseAvgSolverHeur(self);

    end

    methods
        function reset(self)
            for s=1:length(self.solvers)
                self.solvers{s}.reset();
            end
        end
        function setChecks(self,bool)
            for s=1:length(self.solvers)
                self.solvers{s}.setChecks(bool);
            end
        end
    end

    methods
        function varargout = getAvgChainTable(self, varargin)
            % [AVGCHAINTABLE] = GETAVGCHAINTABLE(self)
            % The result recorder captures the returned table together with the solver
            % that produced it -- see LineResultRecorder. Recording an ensemble here
            % rather than in the member solver it delegates to is what keeps an
            % AUTO/LN/ENV/UQ answer from being filed under the member's name.
            [scope, scopeGuard] = LineResultRecorder.enter(); %#ok<ASGLU>
            [varargout{1:max(nargout,1)}] = self.getAvgChainTable_impl(varargin{:});
            LineResultRecorder.capture(scope, self, 'chain', varargout{1});
        end

        function AvgChainTable = getAvgChainTable_impl(self)
            % GETAVGCHAINTABLE_IMPL Implementation of GETAVGCHAINTABLE; see the wrapper above.
            try
                AvgChainTable = self.delegate('getAvgChainTable', 1);
            catch
                line_error(mfilename, 'Fatal error in getAvgChainTable.')
            end
        end

        function varargout = getAvgQLenTable(self, varargin)
            % [AVGQLENTABLE, QT] = GETAVGQLENTABLE(self, Q, keepDisabled)
            % The result recorder captures the returned table together with the solver
            % that produced it -- see LineResultRecorder. Recording an ensemble here
            % rather than in the member solver it delegates to is what keeps an
            % AUTO/LN/ENV/UQ answer from being filed under the member's name.
            [scope, scopeGuard] = LineResultRecorder.enter(); %#ok<ASGLU>
            [varargout{1:max(nargout,1)}] = self.getAvgQLenTable_impl(varargin{:});
            LineResultRecorder.capture(scope, self, 'qlen', varargout{1});
        end

        function [AvgQLenTable, QT] = getAvgQLenTable_impl(self, Q, keepDisabled)
            % GETAVGQLENTABLE_IMPL Implementation of GETAVGQLENTABLE; see the wrapper above.
            if nargin < 2
                Q = [];
            end
            if nargin < 3
                keepDisabled = false;
            end
            [AvgQLenTable, QT] = self.delegate('getAvgQLenTable', 2, Q, keepDisabled);
        end

        function varargout = getAvgTputTable(self, varargin)
            % [AVGTPUTTABLE, TT] = GETAVGTPUTTABLE(self, T, keepDisabled)
            % The result recorder captures the returned table together with the solver
            % that produced it -- see LineResultRecorder. Recording an ensemble here
            % rather than in the member solver it delegates to is what keeps an
            % AUTO/LN/ENV/UQ answer from being filed under the member's name.
            [scope, scopeGuard] = LineResultRecorder.enter(); %#ok<ASGLU>
            [varargout{1:max(nargout,1)}] = self.getAvgTputTable_impl(varargin{:});
            LineResultRecorder.capture(scope, self, 'tput', varargout{1});
        end

        function [AvgTputTable, TT] = getAvgTputTable_impl(self, T, keepDisabled)
            % GETAVGTPUTTABLE_IMPL Implementation of GETAVGTPUTTABLE; see the wrapper above.
            if nargin < 2
                T = [];
            end
            if nargin < 3
                keepDisabled = false;
            end
            [AvgTputTable, TT] = self.delegate('getAvgTputTable', 2, T, keepDisabled);
        end

        function varargout = getAvgRespTTable(self, varargin)
            % [AVGRESPTTABLE, RT] = GETAVGRESPTTABLE(self, R, keepDisabled)
            % The result recorder captures the returned table together with the solver
            % that produced it -- see LineResultRecorder. Recording an ensemble here
            % rather than in the member solver it delegates to is what keeps an
            % AUTO/LN/ENV/UQ answer from being filed under the member's name.
            [scope, scopeGuard] = LineResultRecorder.enter(); %#ok<ASGLU>
            [varargout{1:max(nargout,1)}] = self.getAvgRespTTable_impl(varargin{:});
            LineResultRecorder.capture(scope, self, 'respt', varargout{1});
        end

        function [AvgRespTTable, RT] = getAvgRespTTable_impl(self, R, keepDisabled)
            % GETAVGRESPTTABLE_IMPL Implementation of GETAVGRESPTTABLE; see the wrapper above.
            if nargin < 2
                R = [];
            end
            if nargin < 3
                keepDisabled = false;
            end
            [AvgRespTTable, RT] = self.delegate('getAvgRespTTable', 2, R, keepDisabled);
        end

        function varargout = getAvgUtilTable(self, varargin)
            % [AVGUTILTABLE, UT] = GETAVGUTILTABLE(self, U, keepDisabled)
            % The result recorder captures the returned table together with the solver
            % that produced it -- see LineResultRecorder. Recording an ensemble here
            % rather than in the member solver it delegates to is what keeps an
            % AUTO/LN/ENV/UQ answer from being filed under the member's name.
            [scope, scopeGuard] = LineResultRecorder.enter(); %#ok<ASGLU>
            [varargout{1:max(nargout,1)}] = self.getAvgUtilTable_impl(varargin{:});
            LineResultRecorder.capture(scope, self, 'util', varargout{1});
        end

        function [AvgUtilTable, UT] = getAvgUtilTable_impl(self, U, keepDisabled)
            % GETAVGUTILTABLE_IMPL Implementation of GETAVGUTILTABLE; see the wrapper above.
            if nargin < 2
                U = [];
            end
            if nargin < 3
                keepDisabled = false;
            end
            [AvgUtilTable, UT] = self.delegate('getAvgUtilTable', 2, U, keepDisabled);
        end

        function varargout = getAvgSysTable(self, varargin)
            % [AVGSYSTABLE] = GETAVGSYSTABLE(self)
            % The result recorder captures the returned table together with the solver
            % that produced it -- see LineResultRecorder. Recording an ensemble here
            % rather than in the member solver it delegates to is what keeps an
            % AUTO/LN/ENV/UQ answer from being filed under the member's name.
            [scope, scopeGuard] = LineResultRecorder.enter(); %#ok<ASGLU>
            [varargout{1:max(nargout,1)}] = self.getAvgSysTable_impl(varargin{:});
            LineResultRecorder.capture(scope, self, 'sys', varargout{1});
        end

        function AvgSysTable = getAvgSysTable_impl(self)
            % GETAVGSYSTABLE_IMPL Implementation of GETAVGSYSTABLE; see the wrapper above.
            AvgSysTable = self.delegate('getAvgSysTable', 1);
        end

        function varargout = getAvgNodeTable(self, varargin)
            % [AVGNODETABLE] = GETAVGNODETABLE(self)
            % The result recorder captures the returned table together with the solver
            % that produced it -- see LineResultRecorder. Recording an ensemble here
            % rather than in the member solver it delegates to is what keeps an
            % AUTO/LN/ENV/UQ answer from being filed under the member's name.
            [scope, scopeGuard] = LineResultRecorder.enter(); %#ok<ASGLU>
            [varargout{1:max(nargout,1)}] = self.getAvgNodeTable_impl(varargin{:});
            LineResultRecorder.capture(scope, self, 'node', varargout{1});
        end

        function AvgNodeTable = getAvgNodeTable_impl(self)
            % GETAVGNODETABLE_IMPL Implementation of GETAVGNODETABLE; see the wrapper above.
            AvgNodeTable = self.delegate('getAvgNodeTable', 1);
        end

        function varargout = getAvgTable(self, varargin)
            % [AVGTABLE] = GETAVGTABLE(self)
            % The result recorder captures the returned table together with the solver
            % that produced it -- see LineResultRecorder. Recording an ensemble here
            % rather than in the member solver it delegates to is what keeps an
            % AUTO/LN/ENV/UQ answer from being filed under the member's name.
            [scope, scopeGuard] = LineResultRecorder.enter(); %#ok<ASGLU>
            [varargout{1:max(nargout,1)}] = self.getAvgTable_impl(varargin{:});
            LineResultRecorder.capture(scope, self, 'avg', varargout{1});
        end

        function AvgTable = getAvgTable_impl(self)
            % GETAVGTABLE_IMPL Implementation of GETAVGTABLE; see the wrapper above.
            AvgTable = self.delegate('getAvgTable', 1);
        end

        function [QN,UN,RN,TN,AN,WN] = getAvg(self,Q,U,R,T)
            %[QN,UN,RN,TN] = GETAVG(SELF,Q,U,R,T)

            if nargin>1
                [QN,UN,RN,TN,AN,WN] = self.delegate('getAvg', 6, Q,U,R,T);
            else
                [QN,UN,RN,TN,AN,WN] = self.delegate('getAvg', 6);
            end
        end

        function [QNc,UNc,RNc,TNc] = getAvgChain(self,Q,U,R,T)
            %[QNC,UNC,RNC,TNC] = GETAVGCHAIN(SELF,Q,U,R,T)

            if nargin>1
                [QNc,UNc,RNc,TNc] = self.delegate('getAvgChain', 4, Q,U,R,T);
            else
                [QNc,UNc,RNc,TNc] = self.delegate('getAvgChain', 4);
            end
        end

        function [CNc,XNc] = getAvgSys(self,R,T)
            %[CNC,XNC] = GETAVGSYS(SELF,R,T)

            if nargin>1
                [CNc,XNc] = self.delegate('getAvgSys', 2, R,T);
            else
                [CNc,XNc] = self.delegate('getAvgSys', 2);
            end
        end

        function [QN,UN,RN,TN,AN,WN] = getAvgNode(self,Q,U,R,T,A)
            if nargin>1
                [QN,UN,RN,TN,AN,WN] = self.delegate('getAvgNode', 6, Q,U,R,T,A);
            else
                [QN,UN,RN,TN,AN,WN] = self.delegate('getAvgNode', 6);
            end
        end

        function [AN] = getAvgArvRChain(self,A)
            if nargin>1
                AN = self.delegate('getAvgArvRChain', 1, A);
            else
                AN = self.delegate('getAvgArvRChain', 1);
            end
        end

        function [QN] = getAvgQLenChain(self,Q)
            if nargin>1
                QN = self.delegate('getAvgQLenChain', 1, Q);
            else
                QN = self.delegate('getAvgQLenChain', 1);
            end
        end

        function [UN] = getAvgUtilChain(self,U)
            if nargin>1
                UN = self.delegate('getAvgUtilChain', 1, U);
            else
                UN = self.delegate('getAvgUtilChain', 1);
            end
        end

        function [RN] = getAvgRespTChain(self,R)
            if nargin>1
                RN = self.delegate('getAvgRespTChain', 1, R);
            else
                RN = self.delegate('getAvgRespTChain', 1);
            end
        end

        function [TN] = getAvgTputChain(self,T)
            if nargin>1
                TN = self.delegate('getAvgTputChain', 1, T);
            else
                TN = self.delegate('getAvgTputChain', 1);
            end
        end

        function [RN] = getAvgSysRespT(self,R)
            if nargin>1
                RN = self.delegate('getAvgSysRespT', 1, R);
            else
                RN = self.delegate('getAvgSysRespT', 1);
            end
        end

        function [TN] = getAvgSysTput(self,T)
            if nargin>1
                TN = self.delegate('getAvgSysTput', 1, T);
            else
                TN = self.delegate('getAvgSysTput', 1);
            end
        end

        function [QNt,UNt,TNt] = getTranAvg(self,Qt,Ut,Tt)
            % [QNT,UNT,TNT] = GETTRANAVG(SELF,QT,UT,TT)

            if nargin>1
                [QNt,UNt,TNt] = self.delegate('getTranAvg', 3, Qt,Ut,Tt);
            else
                [QNt,UNt,TNt] = self.delegate('getTranAvg', 3);
            end
        end

        function RD = getTranCdfPassT(self, R)
            % RD = GETTRANCDFPASST(R)

            if nargin>1
                RD = self.delegate('getTranCdfPassT', 1, R);
            else
                RD = self.delegate('getTranCdfPassT', 1);
            end
        end

        function RD = getTranCdfRespT(self, R)
            % RD = GETTRANCDFRESPT(R)

            if nargin>1
                RD = self.delegate('getTranCdfRespT', 1, R);
            else
                RD = self.delegate('getTranCdfRespT', 1);
            end
        end

        function [Pi_t, SSnode] = getTranProb(self, node)
            % [PI, SS] = GETTRANPROB(NODE)
            [Pi_t, SSnode] = self.delegate('getTranProb', 2, node);
        end

        function [Pi_t, SSnode_a] = getTranProbAggr(self, node)
            % [PI, SS] = GETTRANPROBAGGR(NODE)
            [Pi_t, SSnode_a] = self.delegate('getTranProbAggr', 2, node);
        end

        function [Pi_t, SSsys] = getTranProbSys(self)
            % [PI, SS] = GETTRANPROBSYS()
            [Pi_t, SSsys] = self.delegate('getTranProbSys', 2);
        end

        function [Pi_t, SSsysa] = getTranProbSysAggr(self)
            % [PI, SS] = GETTRANPROBSYSAGGR()
            [Pi_t, SSsysa] = self.delegate('getTranProbSysAggr', 2);
        end

        function sampleNodeState = sample(self, node, numEvents)
            sampleNodeState = self.delegate('sample', 1, node, numEvents);
        end

        function stationStateAggr = sampleAggr(self, node, numEvents)
            stationStateAggr = self.delegate('sampleAggr', 1, node, numEvents);
        end

        function tranSysState = sampleSys(self, numEvents)
            tranSysState = self.delegate('sampleSys', 1, numEvents);
        end

        function sysStateAggr = sampleSysAggr(self, numEvents)
            sysStateAggr = self.delegate('sampleSysAggr', 1, numEvents);
        end

        function RD = getCdfRespT(self, R)
            if nargin>1
                RD = self.delegate('getCdfRespT', 1, R);
            else
                RD = self.delegate('getCdfRespT', 1);
            end
        end

        function varargout = getProb(self, varargin)
            % The result recorder captures the scalar together with the solver that
            % produced it -- see LineResultRecorder. Recording the ensemble here
            % rather than in the member it delegates to keeps an AUTO answer from
            % being filed under the member's name.
            [scope, scopeGuard] = LineResultRecorder.enter(); %#ok<ASGLU>
            [varargout{1:max(nargout,1)}] = self.getProb_impl(varargin{:});
            LineResultRecorder.captureScalar(scope, self, 'prob', varargout{1});
        end

        function Pnir = getProb_impl(self, node, state)
            % GETPROB_IMPL Implementation of GETPROB; see the wrapper above.
            Pnir = self.delegate('getProb', 1, node, state);
        end

        function varargout = getProbAggr(self, varargin)
            % The result recorder captures the scalar together with the solver that
            % produced it -- see LineResultRecorder. Recording the ensemble here
            % rather than in the member it delegates to keeps an AUTO answer from
            % being filed under the member's name.
            [scope, scopeGuard] = LineResultRecorder.enter(); %#ok<ASGLU>
            [varargout{1:max(nargout,1)}] = self.getProbAggr_impl(varargin{:});
            LineResultRecorder.captureScalar(scope, self, 'probAggr', varargout{1});
        end

        function Pnir = getProbAggr_impl(self, node, state_a)
            % GETPROBAGGR_IMPL Implementation of GETPROBAGGR; see the wrapper above.
            Pnir = self.delegate('getProbAggr', 1, node, state_a);
        end

        function varargout = getProbSys(self, varargin)
            % The result recorder captures the scalar together with the solver that
            % produced it -- see LineResultRecorder. Recording the ensemble here
            % rather than in the member it delegates to keeps an AUTO answer from
            % being filed under the member's name.
            [scope, scopeGuard] = LineResultRecorder.enter(); %#ok<ASGLU>
            [varargout{1:max(nargout,1)}] = self.getProbSys_impl(varargin{:});
            LineResultRecorder.captureScalar(scope, self, 'probSys', varargout{1});
        end

        function Pn = getProbSys_impl(self)
            % GETPROBSYS_IMPL Implementation of GETPROBSYS; see the wrapper above.
            Pn = self.delegate('getProbSys',1);
        end

        function varargout = getProbSysAggr(self, varargin)
            % The result recorder captures the scalar together with the solver that
            % produced it -- see LineResultRecorder. Recording the ensemble here
            % rather than in the member it delegates to keeps an AUTO answer from
            % being filed under the member's name.
            [scope, scopeGuard] = LineResultRecorder.enter(); %#ok<ASGLU>
            [varargout{1:max(nargout,1)}] = self.getProbSysAggr_impl(varargin{:});
            LineResultRecorder.captureScalar(scope, self, 'probSysAggr', varargout{1});
        end

        function Pn = getProbSysAggr_impl(self)
            % GETPROBSYSAGGR_IMPL Implementation of GETPROBSYSAGGR; see the wrapper above.
            Pn = self.delegate('getProbSysAggr',1);
        end

        function [logNormConst] = getProbNormConstAggr(self)
            logNormConst = self.delegate('getProbNormConstAggr',1);
        end

        % Basic metric methods
        function QN = getAvgQLen(self)
            % QN = GETAVGQLEN()
            % Compute average queue-lengths at steady-state
            QN = self.delegate('getAvgQLen', 1);
        end

        function UN = getAvgUtil(self)
            % UN = GETAVGUTIL()
            % Compute average utilizations at steady-state
            UN = self.delegate('getAvgUtil', 1);
        end

        function RN = getAvgRespT(self)
            % RN = GETAVGRESPT()
            % Compute average response times at steady-state
            RN = self.delegate('getAvgRespT', 1);
        end

        function WN = getAvgResidT(self)
            % WN = GETAVGRESIDT()
            % Compute average residence times at steady-state
            WN = self.delegate('getAvgResidT', 1);
        end

        function WT = getAvgWaitT(self)
            % WT = GETAVGWAITT()
            % Compute average waiting time in queue excluding service
            WT = self.delegate('getAvgWaitT', 1);
        end

        function TN = getAvgTput(self)
            % TN = GETAVGTPUT()
            % Compute average throughputs at steady-state
            TN = self.delegate('getAvgTput', 1);
        end

        function AN = getAvgArvR(self)
            % AN = GETAVGARVR()
            % Compute average arrival rate at steady-state
            AN = self.delegate('getAvgArvR', 1);
        end

        % Additional chain methods
        function [WN] = getAvgResidTChain(self, W)
            if nargin > 1
                WN = self.delegate('getAvgResidTChain', 1, W);
            else
                WN = self.delegate('getAvgResidTChain', 1);
            end
        end

        function [QN] = getAvgNodeQLenChain(self, Q)
            if nargin > 1
                QN = self.delegate('getAvgNodeQLenChain', 1, Q);
            else
                QN = self.delegate('getAvgNodeQLenChain', 1);
            end
        end

        function [UN] = getAvgNodeUtilChain(self, U)
            if nargin > 1
                UN = self.delegate('getAvgNodeUtilChain', 1, U);
            else
                UN = self.delegate('getAvgNodeUtilChain', 1);
            end
        end

        function [RN] = getAvgNodeRespTChain(self, R)
            if nargin > 1
                RN = self.delegate('getAvgNodeRespTChain', 1, R);
            else
                RN = self.delegate('getAvgNodeRespTChain', 1);
            end
        end

        function [WN] = getAvgNodeResidTChain(self, W)
            if nargin > 1
                WN = self.delegate('getAvgNodeResidTChain', 1, W);
            else
                WN = self.delegate('getAvgNodeResidTChain', 1);
            end
        end

        function [TN] = getAvgNodeTputChain(self, T)
            if nargin > 1
                TN = self.delegate('getAvgNodeTputChain', 1, T);
            else
                TN = self.delegate('getAvgNodeTputChain', 1);
            end
        end

        function [AN] = getAvgNodeArvRChain(self, A)
            if nargin > 1
                AN = self.delegate('getAvgNodeArvRChain', 1, A);
            else
                AN = self.delegate('getAvgNodeArvRChain', 1);
            end
        end

        % Additional table method
        function varargout = getAvgNodeChainTable(self, varargin)
            % The result recorder captures the returned table together with the solver
            % that produced it -- see LineResultRecorder. Recording an ensemble here
            % rather than in the member solver it delegates to is what keeps an
            % AUTO/LN/ENV/UQ answer from being filed under the member's name.
            [scope, scopeGuard] = LineResultRecorder.enter(); %#ok<ASGLU>
            [varargout{1:max(nargout,1)}] = self.getAvgNodeChainTable_impl(varargin{:});
            LineResultRecorder.capture(scope, self, 'nodechain', varargout{1});
        end

        function [AvgNodeChainTable, QTc, UTc, RTc, WTc, ATc, TTc] = getAvgNodeChainTable_impl(self, Q, U, R, T)
            % GETAVGNODECHAINTABLE_IMPL Implementation of GETAVGNODECHAINTABLE; see the wrapper above.
            if nargin > 1
                [AvgNodeChainTable, QTc, UTc, RTc, WTc, ATc, TTc] = self.delegate('getAvgNodeChainTable', 7, Q, U, R, T);
            else
                [AvgNodeChainTable, QTc, UTc, RTc, WTc, ATc, TTc] = self.delegate('getAvgNodeChainTable', 7);
            end
        end

        % Probability method
        function Pmarg = getProbMarg(self, node, jobclass, state_m)
            % PMARG = GETPROBMARG(NODE, JOBCLASS, STATE_M)
            % Return marginalized state probability for station and class
            Pmarg = self.delegate('getProbMarg', 1, node, jobclass, state_m);
        end

        % Distribution method
        function RD = getCdfPassT(self, R)
            % RD = GETCDFPASST(R)
            % Return cumulative distribution of passage times at steady-state
            if nargin > 1
                RD = self.delegate('getCdfPassT', 1, R);
            else
                RD = self.delegate('getCdfPassT', 1);
            end
        end

        % Percentile method
        function [PercRT, PercTable] = getPerctRespT(self, percentiles, jobclass)
            % [PERCRT, PERCTABLE] = GETPERCTRESPT(SELF, PERCENTILES, JOBCLASS)
            % Extract response time percentiles from CDF or solver-specific results
            if nargin < 3
                [PercRT, PercTable] = self.delegate('getPerctRespT', 2, percentiles);
            else
                [PercRT, PercTable] = self.delegate('getPerctRespT', 2, percentiles, jobclass);
            end
        end

        % Handle methods
        function [Q, U, R, T, A, W] = getAvgHandles(self)
            % [Q,U,R,T,A,W] = GETAVGHANDLES()
            [Q, U, R, T, A, W] = self.delegate('getAvgHandles', 6);
        end

        function [Qt, Ut, Tt] = getTranHandles(self)
            % [QT,UT,TT] = GETTRANHANDLES()
            [Qt, Ut, Tt] = self.delegate('getTranHandles', 3);
        end

        function Q = getAvgQLenHandles(self)
            % Q = GETAVGQLENHANDLES()
            Q = self.delegate('getAvgQLenHandles', 1);
        end

        function U = getAvgUtilHandles(self)
            % U = GETAVGUTILHANDLES()
            U = self.delegate('getAvgUtilHandles', 1);
        end

        function R = getAvgRespTHandles(self)
            % R = GETAVGRESPTHANDLES()
            R = self.delegate('getAvgRespTHandles', 1);
        end

        function T = getAvgTputHandles(self)
            % T = GETAVGTPUTHANDLES()
            T = self.delegate('getAvgTputHandles', 1);
        end

        function A = getAvgArvRHandles(self)
            % A = GETAVGARVRHANDLES()
            A = self.delegate('getAvgArvRHandles', 1);
        end

        function W = getAvgResidTHandles(self)
            % W = GETAVGRESIDTHANDLES()
            W = self.delegate('getAvgResidTHandles', 1);
        end

        % Aliases for table methods
        function avg_table = avgTable(self)
            % AVGTABLE Alias for getAvgTable
            avg_table = self.getAvgTable();
        end

        function avg_sys_table = avgSysTable(self)
            % AVGSYSTABLE Alias for getAvgSysTable
            avg_sys_table = self.getAvgSysTable();
        end

        function avg_node_table = avgNodeTable(self)
            % AVGNODETABLE Alias for getAvgNodeTable
            avg_node_table = self.getAvgNodeTable();
        end

        function avg_chain_table = avgChainTable(self)
            % AVGCHAINTABLE Alias for getAvgChainTable
            avg_chain_table = self.getAvgChainTable();
        end

        function avg_node_chain_table = avgNodeChainTable(self)
            % AVGNODECHAINTABLE Alias for getAvgNodeChainTable
            avg_node_chain_table = self.getAvgNodeChainTable();
        end

        % Table -> T short aliases
        function avg_table = avgT(self)
            % AVGT Short alias for avgTable
            avg_table = self.avgTable();
        end

        function avg_sys_table = avgSysT(self)
            % AVGSYST Short alias for avgSysTable
            avg_sys_table = self.avgSysTable();
        end

        function avg_node_table = avgNodeT(self)
            % AVGNODET Short alias for avgNodeTable
            avg_node_table = self.avgNodeTable();
        end

        function avg_chain_table = avgChainT(self)
            % AVGCHAINT Short alias for avgChainTable
            avg_chain_table = self.avgChainTable();
        end

        function avg_node_chain_table = avgNodeChainT(self)
            % AVGNODECHAINT Short alias for avgNodeChainTable
            avg_node_chain_table = self.avgNodeChainTable();
        end

        % Aliases for composite methods
        function varargout = avgChain(self, varargin)
            % AVGCHAIN Alias for getAvgChain
            [varargout{1:nargout}] = self.getAvgChain(varargin{:});
        end

        function varargout = avgSys(self, varargin)
            % AVGSYS Alias for getAvgSys
            [varargout{1:nargout}] = self.getAvgSys(varargin{:});
        end

        function varargout = avgNode(self, varargin)
            % AVGNODE Alias for getAvgNode
            [varargout{1:nargout}] = self.getAvgNode(varargin{:});
        end

        function varargout = avg(self, varargin)
            % AVG Alias for getAvg
            [varargout{1:nargout}] = self.getAvg(varargin{:});
        end

        function sys_resp_time = avgSysRespT(self, varargin)
            % AVGSYSRESPT Alias for getAvgSysRespT
            sys_resp_time = self.getAvgSysRespT(varargin{:});
        end

        function sys_tput = avgSysTput(self, varargin)
            % AVGSYSTPUT Alias for getAvgSysTput
            sys_tput = self.getAvgSysTput(varargin{:});
        end

        % Aliases for chain methods
        function arvr_chain = avgArvRChain(self, varargin)
            % AVGARVCHAIN Alias for getAvgArvRChain
            arvr_chain = self.getAvgArvRChain(varargin{:});
        end

        function qlen_chain = avgQLenChain(self, varargin)
            % AVGQLENCHAIN Alias for getAvgQLenChain
            qlen_chain = self.getAvgQLenChain(varargin{:});
        end

        function util_chain = avgUtilChain(self, varargin)
            % AVGUTILCHAIN Alias for getAvgUtilChain
            util_chain = self.getAvgUtilChain(varargin{:});
        end

        function resp_t_chain = avgRespTChain(self, varargin)
            % AVGRESPTCHAIN Alias for getAvgRespTChain
            resp_t_chain = self.getAvgRespTChain(varargin{:});
        end

        function resid_t_chain = avgResidTChain(self, varargin)
            % AVGRESIDTCHAIN Alias for getAvgResidTChain
            resid_t_chain = self.getAvgResidTChain(varargin{:});
        end

        function tput_chain = avgTputChain(self, varargin)
            % AVGTPUTCHAIN Alias for getAvgTputChain
            tput_chain = self.getAvgTputChain(varargin{:});
        end

        % Aliases for node chain methods
        function node_arvr_chain = avgNodeArvRChain(self, varargin)
            % AVGNODERVRCHAIN Alias for getAvgNodeArvRChain
            node_arvr_chain = self.getAvgNodeArvRChain(varargin{:});
        end

        function node_qlen_chain = avgNodeQLenChain(self, varargin)
            % AVGNODEQLENCHAIN Alias for getAvgNodeQLenChain
            node_qlen_chain = self.getAvgNodeQLenChain(varargin{:});
        end

        function node_util_chain = avgNodeUtilChain(self, varargin)
            % AVGNODEUTILCHAIN Alias for getAvgNodeUtilChain
            node_util_chain = self.getAvgNodeUtilChain(varargin{:});
        end

        function node_resp_t_chain = avgNodeRespTChain(self, varargin)
            % AVGNODERESPTCHAIN Alias for getAvgNodeRespTChain
            node_resp_t_chain = self.getAvgNodeRespTChain(varargin{:});
        end

        function node_resid_t_chain = avgNodeResidTChain(self, varargin)
            % AVGNODERESIDTCHAIN Alias for getAvgNodeResidTChain
            node_resid_t_chain = self.getAvgNodeResidTChain(varargin{:});
        end

        function node_tput_chain = avgNodeTputChain(self, varargin)
            % AVGNODETPUTCHAIN Alias for getAvgNodeTputChain
            node_tput_chain = self.getAvgNodeTputChain(varargin{:});
        end

        % Aliases for transient methods
        function varargout = tranAvg(self, varargin)
            % TRANAVG Alias for getTranAvg
            [varargout{1:nargout}] = self.getTranAvg(varargin{:});
        end

        function rd = tranCdfRespT(self, varargin)
            % TRANCDFRESPT Alias for getTranCdfRespT
            rd = self.getTranCdfRespT(varargin{:});
        end

        function rd = tranCdfPassT(self, varargin)
            % TRANCDFPASST Alias for getTranCdfPassT
            rd = self.getTranCdfPassT(varargin{:});
        end

        function varargout = tranProb(self, varargin)
            % TRANPROB Alias for getTranProb
            [varargout{1:nargout}] = self.getTranProb(varargin{:});
        end

        function varargout = tranProbAggr(self, varargin)
            % TRANPROBAGGR Alias for getTranProbAggr
            [varargout{1:nargout}] = self.getTranProbAggr(varargin{:});
        end

        function varargout = tranProbSys(self)
            % TRANPROBSYS Alias for getTranProbSys
            [varargout{1:nargout}] = self.getTranProbSys();
        end

        function varargout = tranProbSysAggr(self)
            % TRANPROBSYSAGGR Alias for getTranProbSysAggr
            [varargout{1:nargout}] = self.getTranProbSysAggr();
        end

        % Aliases for CDF methods
        function rd = cdfRespT(self, varargin)
            % CDFRESPT Alias for getCdfRespT
            rd = self.getCdfRespT(varargin{:});
        end

        function rd = cdfPassT(self, varargin)
            % CDFPASST Alias for getCdfPassT
            rd = self.getCdfPassT(varargin{:});
        end

        function varargout = perctRespT(self, varargin)
            % PERCTRESPT Alias for getPerctRespT
            [varargout{1:nargout}] = self.getPerctRespT(varargin{:});
        end

        % Aliases for probability methods
        function pstate = prob(self, varargin)
            % PROB Alias for getProb
            pstate = self.getProb(varargin{:});
        end

        function psysstate = probSys(self)
            % PROBSYS Alias for getProbSys
            psysstate = self.getProbSys();
        end

        function pnir = probAggr(self, varargin)
            % PROBAGGR Alias for getProbAggr
            pnir = self.getProbAggr(varargin{:});
        end

        function pnjoint = probSysAggr(self)
            % PROBSYSAGGR Alias for getProbSysAggr
            pnjoint = self.getProbSysAggr();
        end

        function pmarg = probMarg(self, varargin)
            % PROBMARG Alias for getProbMarg
            pmarg = self.getProbMarg(varargin{:});
        end

        function lnormconst = probNormConstAggr(self)
            % PROBNORMCONSTAGGR Alias for getProbNormConstAggr
            lnormconst = self.getProbNormConstAggr();
        end

        % Aliases for handle methods
        function varargout = avgHandles(self)
            % AVGHANDLES Alias for getAvgHandles
            [varargout{1:nargout}] = self.getAvgHandles();
        end

        function varargout = tranHandles(self)
            % TRANHANDLES Alias for getTranHandles
            [varargout{1:nargout}] = self.getTranHandles();
        end

        function q = avgQLenHandles(self)
            % AVGQLENHANDLES Alias for getAvgQLenHandles
            q = self.getAvgQLenHandles();
        end

        function u = avgUtilHandles(self)
            % AVGUTILHANDLES Alias for getAvgUtilHandles
            u = self.getAvgUtilHandles();
        end

        function r = avgRespTHandles(self)
            % AVGRESPTHANDLES Alias for getAvgRespTHandles
            r = self.getAvgRespTHandles();
        end

        function t = avgTputHandles(self)
            % AVGTPUTHANDLES Alias for getAvgTputHandles
            t = self.getAvgTputHandles();
        end

        function a = avgArvRHandles(self)
            % AVGARVRHANDLES Alias for getAvgArvRHandles
            a = self.getAvgArvRHandles();
        end

        function w = avgResidTHandles(self)
            % AVGRESIDTHANDLES Alias for getAvgResidTHandles
            w = self.getAvgResidTHandles();
        end

        % Aliases for basic metric methods
        function qn = avgQLen(self)
            % AVGQLEN Alias for getAvgQLen
            qn = self.getAvgQLen();
        end

        function un = avgUtil(self)
            % AVGUTIL Alias for getAvgUtil
            un = self.getAvgUtil();
        end

        function rn = avgRespT(self)
            % AVGRESPT Alias for getAvgRespT
            rn = self.getAvgRespT();
        end

        function wn = avgResidT(self)
            % AVGRESIDT Alias for getAvgResidT
            wn = self.getAvgResidT();
        end

        function wt = avgWaitT(self)
            % AVGWAITT Alias for getAvgWaitT
            wt = self.getAvgWaitT();
        end

        function tn = avgTput(self)
            % AVGTPUT Alias for getAvgTput
            tn = self.getAvgTput();
        end

        function an = avgArvR(self)
            % AVGARVR Alias for getAvgArvR
            an = self.getAvgArvR();
        end

        % Aliases for table methods with parameters
        function varargout = avgQLenTable(self, varargin)
            % AVGQLENTABLE Alias for getAvgQLenTable
            [varargout{1:nargout}] = self.getAvgQLenTable(varargin{:});
        end

        function varargout = avgUtilTable(self, varargin)
            % AVGUTILTABLE Alias for getAvgUtilTable
            [varargout{1:nargout}] = self.getAvgUtilTable(varargin{:});
        end

        function varargout = avgRespTTable(self, varargin)
            % AVGRESPTTABLE Alias for getAvgRespTTable
            [varargout{1:nargout}] = self.getAvgRespTTable(varargin{:});
        end

        function varargout = avgTputTable(self, varargin)
            % AVGTPUTTABLE Alias for getAvgTputTable
            [varargout{1:nargout}] = self.getAvgTputTable(varargin{:});
        end

        %% LayeredNetwork / EnsembleSolver methods
        % These methods are available when solving LayeredNetwork models

        function [QN, UN, RN, TN, AN, WN] = getEnsembleAvg(self)
            % [QN, UN, RN, TN, AN, WN] = GETENSEMBLEAVG()
            % Get average performance metrics for LayeredNetwork ensemble models
            [QN, UN, RN, TN, AN, WN] = self.delegate('getEnsembleAvg', 6);
        end

        function solver = getSolver(self, e)
            % SOLVER = GETSOLVER(E)
            % Get solver for ensemble model e (LayeredNetwork only)
            solver = self.delegate('getSolver', 1, e);
        end

        function solver = setSolver(self, solver, e)
            % SOLVER = SETSOLVER(SOLVER, E)
            % Set solver for ensemble model e (LayeredNetwork only)
            if nargin < 3
                solver = self.delegate('setSolver', 1, solver);
            else
                solver = self.delegate('setSolver', 1, solver, e);
            end
        end

        function E = getNumberOfModels(self)
            % E = GETNUMBEROFMODELS()
            % Get number of ensemble models (LayeredNetwork only)
            E = self.delegate('getNumberOfModels', 1);
        end

        function it = getIteration(self)
            % IT = GETITERATION()
            % Get current iteration number (LayeredNetwork only)
            it = self.delegate('getIteration', 1);
        end

        function AvgTables = getEnsembleAvgTables(self)
            % AVGTABLES = GETENSEMBLEAVGTABLES()
            % Get average tables for all ensemble models (LayeredNetwork only)
            AvgTables = self.delegate('getEnsembleAvgTables', 1);
        end

        function state = get_state(self)
            % STATE = GET_STATE()
            % Export current solver state for continuation (SolverLN only)
            state = self.delegate('get_state', 1);
        end

        function set_state(self, state)
            % SET_STATE(STATE)
            % Import solution state for continuation (SolverLN only)
            self.delegate('set_state', 0, state);
        end

        function update_solver(self, solverFactory)
            % UPDATE_SOLVER(SOLVERFACTORY)
            % Change the solver for all layers (SolverLN only)
            self.delegate('update_solver', 0, solverFactory);
        end

        % Aliases for LayeredNetwork/EnsembleSolver methods
        function varargout = ensembleAvg(self)
            % ENSEMBLEAVG Alias for getEnsembleAvg
            [varargout{1:nargout}] = self.getEnsembleAvg();
        end

        function solver = solver(self, e)
            % SOLVER Alias for getSolver
            solver = self.getSolver(e);
        end

        function it = iteration(self)
            % ITERATION Alias for getIteration
            it = self.getIteration();
        end

        function e = numberOfModels(self)
            % NUMBEROFMODELS Alias for getNumberOfModels
            e = self.getNumberOfModels();
        end

        function avg_tables = ensembleAvgTables(self)
            % ENSEMBLEAVGTABLES Alias for getEnsembleAvgTables
            avg_tables = self.getEnsembleAvgTables();
        end

        function avg_tables = getEnsembleAvgTs(self)
            % GETENSEMBLEAVGTS Short alias for getEnsembleAvgTables
            avg_tables = self.getEnsembleAvgTables();
        end

        function avg_tables = ensembleAvgTs(self)
            % ENSEMBLEAVGTS Short alias for ensembleAvgTables
            avg_tables = self.ensembleAvgTables();
        end

        % Specialized result tables. SolverAUTO is not a NetworkSolver
        % subclass, so every accessor a delegate exposes must be repeated here
        % or it is simply unreachable through AUTO.
        function ctmc = ctmcSolver(self)
            % CTMC = CTMCSOLVER()
            % The CTMC candidate, built on demand. State-space and generator
            % accessors are CTMC-only concepts, so they resolve here rather
            % than through the ranked selection.
            ctmc = [];
            if numel(self.solvers) >= self.CANDIDATE_CTMC
                ctmc = self.solvers{self.CANDIDATE_CTMC};
            end
            if isempty(ctmc) || ~isa(ctmc,'SolverCTMC')
                ctmc = SolverCTMC(self.model);
            end
        end

        function varargout = getGenerator(self, varargin)
            % [INFGEN, EVENTFILT, SYNCHINFO] = GETGENERATOR(OPTIONS)
            ctmcs = self.ctmcSolver();
            [varargout{1:max(nargout,1)}] = ctmcs.getGenerator(varargin{:});
        end

        function varargout = getSymbolicGenerator(self, varargin)
            % [INFGEN, EVENTFILT, SYNCHINFO, STATESPACE, NODESTATESPACE] = GETSYMBOLICGENERATOR(INVERTSYMBOL)
            ctmcs = self.ctmcSolver();
            [varargout{1:max(nargout,1)}] = ctmcs.getSymbolicGenerator(varargin{:});
        end

        function varargout = getStateSpace(self, varargin)
            % [STATESPACE, NODESTATESPACE] = GETSTATESPACE(OPTIONS)
            ctmcs = self.ctmcSolver();
            [varargout{1:max(nargout,1)}] = ctmcs.getStateSpace(varargin{:});
        end

        function varargout = getStateSpaceAggr(self, varargin)
            % [STATESPACEAGGR, NODESTATESPACE] = GETSTATESPACEAGGR(OPTIONS)
            ctmcs = self.ctmcSolver();
            [varargout{1:max(nargout,1)}] = ctmcs.getStateSpaceAggr(varargin{:});
        end

        function varargout = getInfGen(self, varargin)
            % [INFGEN, EVENTFILT, EV] = GETINFGEN(OPTIONS)
            ctmcs = self.ctmcSolver();
            [varargout{1:max(nargout,1)}] = ctmcs.getInfGen(varargin{:});
        end

        function varargout = getTransMat(self, varargin)
            % P = GETTRANSMAT()
            ctmcs = self.ctmcSolver();
            [varargout{1:max(nargout,1)}] = ctmcs.getTransMat(varargin{:});
        end

        function varargout = getMarkedCTMC(self, varargin)
            % MCTMC = GETMARKEDCTMC(OPTIONS)
            ctmcs = self.ctmcSolver();
            [varargout{1:max(nargout,1)}] = ctmcs.getMarkedCTMC(varargin{:});
        end

        function varargout = getSymbolicSolution(self, varargin)
            % [PI, SS, SYMBOLS] = GETSYMBOLICSOLUTION(...)
            ctmcs = self.ctmcSolver();
            [varargout{1:max(nargout,1)}] = ctmcs.getSymbolicSolution(varargin{:});
        end

        function varargout = getSensitivity(self, varargin)
            % [SENS, ...] = GETSENSITIVITY(...)
            ctmcs = self.ctmcSolver();
            [varargout{1:max(nargout,1)}] = ctmcs.getSensitivity(varargin{:});
        end

        function varargout = getSensitivityRanking(self, varargin)
            % RANKING = GETSENSITIVITYRANKING(...)
            ctmcs = self.ctmcSolver();
            [varargout{1:max(nargout,1)}] = ctmcs.getSensitivityRanking(varargin{:});
        end

        function varargout = getCdfSysRespT(self, varargin)
            % RD = GETCDFSYSRESPT(R)
            ctmcs = self.ctmcSolver();
            [varargout{1:max(nargout,1)}] = ctmcs.getCdfSysRespT(varargin{:});
        end

        % Fluid-only accessors. Like the state-space ones above, these are
        % properties of the ODE representation rather than of the model, so
        % they resolve on the fluid family instead of the ranked selection.
        function fld = fldSolver(self)
            % FLD = FLDSOLVER()
            fld = [];
            if numel(self.solvers) >= self.CANDIDATE_FLUID
                fld = self.solvers{self.CANDIDATE_FLUID};
            end
            if isempty(fld) || ~isa(fld,'SolverFluid')
                fld = SolverFluid(self.model);
            end
        end

        function varargout = getCdfPT(self, varargin)
            % RD = GETCDFPT(R)
            [varargout{1:max(nargout,1)}] = self.fldSolver().getCdfPT(varargin{:});
        end

        function varargout = getAvgAoI(self, varargin)
            % [AOI, PAOI, AOITABLE] = GETAVGAOI()
            [varargout{1:max(nargout,1)}] = self.fldSolver().getAvgAoI(varargin{:});
        end

        function varargout = getCdfAoI(self, varargin)
            % [AOI_CDF, PAOI_CDF] = GETCDFAOI(T)
            [varargout{1:max(nargout,1)}] = self.fldSolver().getCdfAoI(varargin{:});
        end

        function varargout = getMoments(self, varargin)
            % MOMENTS = GETMOMENTS()
            [varargout{1:max(nargout,1)}] = self.fldSolver().getMoments(varargin{:});
        end

        function varargout = getTranAvgVar(self, varargin)
            % [T, QVART, SIGMAT] = GETTRANAVGVAR()
            [varargout{1:max(nargout,1)}] = self.fldSolver().getTranAvgVar(varargin{:});
        end

        function varargout = getJacobian(self, varargin)
            % [J, RHS, VARS, EQUILIBRIA] = GETJACOBIAN(OPTIONS)
            [varargout{1:max(nargout,1)}] = self.fldSolver().getJacobian(varargin{:});
        end

        function varargout = exportODEs(self, varargin)
            % [TEX, SYS] = EXPORTODES(FILENAME, NOTATION)
            [varargout{1:max(nargout,1)}] = self.fldSolver().exportODEs(varargin{:});
        end

        function varargout = getSymbolicDrift(self, varargin)
            % [RHS, VARS, SYS] = GETSYMBOLICDRIFT(OPTIONS)
            [varargout{1:max(nargout,1)}] = self.fldSolver().getSymbolicDrift(varargin{:});
        end

        function nc = ncSolver(self)
            % NC = NCSOLVER()
            nc = [];
            if numel(self.solvers) >= self.CANDIDATE_NC
                nc = self.solvers{self.CANDIDATE_NC};
            end
            if isempty(nc) || ~isa(nc,'SolverNC')
                nc = SolverNC(self.model);
            end
        end

        function varargout = getNormalizingConstant(self, varargin)
            % [NORMCONST, LNORMCONST] = GETNORMALIZINGCONSTANT()
            [varargout{1:max(nargout,1)}] = self.ncSolver().getNormalizingConstant(varargin{:});
        end

        function ba = baSolver(self)
            % BA = BASOLVER()
            % Bounds are never in the ranked candidate list, since they return
            % an interval rather than a point estimate, so they are built here.
            ba = [];
            if numel(self.solvers) == 1 && isa(self.solvers{1},'SolverBA')
                ba = self.solvers{1};
            end
            if isempty(ba)
                opt = self.options;
                if ~SolverAUTO.isBoundMethod(opt.method)
                    opt.method = 'auto';
                end
                ba = SolverBA(self.model, opt);
            end
        end

        function bounds = getBounds(self)
            % BOUNDS = GETBOUNDS()
            bounds = self.baSolver().getBounds();
        end

        function varargout = getBoundsTable(self, varargin)
            % BOUNDSTABLE = GETBOUNDSTABLE(KEEPDISABLED)
            [varargout{1:max(nargout,1)}] = self.baSolver().getBoundsTable(varargin{:});
        end

        function uq = uqSolver(self)
            % UQ = UQSOLVER()
            uq = [];
            if numel(self.solvers) == 1 && isa(self.solvers{1},'UQ')
                uq = self.solvers{1};
            end
            if isempty(uq)
                inner = self.options;
                inner.verbose = 0;
                uq = SolverUQ(self.model, @(m) LINE(m, inner), self.options);
            end
        end

        function varargout = getPosteriorTable(self, varargin)
            % PTABLE = GETPOSTERIORTABLE()
            [varargout{1:max(nargout,1)}] = self.uqSolver().getPosteriorTable(varargin{:});
        end

        function varargout = getPosteriorDist(self, varargin)
            % EMPDIST = GETPOSTERIORDIST(METRIC, STATION, CLASS)
            [varargout{1:max(nargout,1)}] = self.uqSolver().getPosteriorDist(varargin{:});
        end

        % Metrics that several families implement: the ranked selection picks
        % the one that supports the model, exactly as for the mean metrics.
        function varargout = getAvgReward(self, varargin)
            % [REWARDS, TABLE] = GETAVGREWARD(...)
            [varargout{1:max(nargout,1)}] = self.delegate('getAvgReward', max(nargout,1), varargin{:});
        end

        function varargout = getTranReward(self, varargin)
            % [T, RT] = GETTRANREWARD(...)
            [varargout{1:max(nargout,1)}] = self.delegate('getTranReward', max(nargout,1), varargin{:});
        end

        function varargout = getSjrnT(self, varargin)
            % RD = GETSJRNT(R)
            [varargout{1:max(nargout,1)}] = self.delegate('getSjrnT', max(nargout,1), varargin{:});
        end

        function varargout = sjrnT(self, varargin)
            % SJRNT Alias for getSjrnT
            [varargout{1:max(nargout,1)}] = self.delegate('sjrnT', max(nargout,1), varargin{:});
        end

        function bool = supportsTransientAnalysis(self)
            % BOOL = SUPPORTSTRANSIENTANALYSIS()
            % True when any feasible candidate can run a transient analysis,
            % since that is the one the fallback would reach for it.
            bool = false;
            if numel(self.solvers) == 1
                try
                    bool = self.solvers{1}.supportsTransientAnalysis();
                catch
                    bool = false;
                end
                return
            end
            for s = 1:length(self.solvers)
                if isempty(self.solvers{s})
                    continue
                end
                try
                    if self.solvers{s}.supportsTransientAnalysis()
                        bool = true;
                        return
                    end
                catch
                    % A family without the method cannot answer transiently.
                end
            end
        end

        % Layered and random-environment accessors, resolved on the delegate
        % that decomposes the model.
        function varargout = getStageAvg(self, varargin)
            % [QN,UN,RN,TN] = GETSTAGEAVG(...)
            [varargout{1:max(nargout,1)}] = self.delegate('getStageAvg', max(nargout,1), varargin{:});
        end

        function varargout = getTranAvgCoupled(self, varargin)
            % [QT,UT,TT] = GETTRANAVGCOUPLED(QT,UT,TT)
            [varargout{1:max(nargout,1)}] = self.delegate('getTranAvgCoupled', max(nargout,1), varargin{:});
        end

        function varargout = getTranAvgDecoupled(self, varargin)
            % [QT,UT,TT] = GETTRANAVGDECOUPLED(QT,UT,TT)
            [varargout{1:max(nargout,1)}] = self.delegate('getTranAvgDecoupled', max(nargout,1), varargin{:});
        end

        function varargout = computeSojournCdf(self, varargin)
            % CDF = COMPUTESOJOURNCDF(E, T)
            [varargout{1:max(nargout,1)}] = self.delegate('computeSojournCdf', max(nargout,1), varargin{:});
        end

        function varargout = getSamplePathTable(self, varargin)
            % [T, SEGMENTRESULTS] = GETSAMPLEPATHTABLE(SAMPLEPATH)
            [varargout{1:max(nargout,1)}] = self.delegate('getSamplePathTable', max(nargout,1), varargin{:});
        end

        function varargout = getAvgCacheTable(self, varargin)
            % CACHEAVGTABLE = GETAVGCACHETABLE()
            % The result recorder captures the returned table together with the solver
            % that produced it -- see LineResultRecorder. Recording an ensemble here
            % rather than in the member solver it delegates to is what keeps an
            % AUTO/LN/ENV/UQ answer from being filed under the member's name.
            [scope, scopeGuard] = LineResultRecorder.enter(); %#ok<ASGLU>
            [varargout{1:max(nargout,1)}] = self.getAvgCacheTable_impl(varargin{:});
            LineResultRecorder.capture(scope, self, 'cache', varargout{1});
        end

        function varargout = getAvgCacheTable_impl(self, varargin)
            % GETAVGCACHETABLE_IMPL Implementation of GETAVGCACHETABLE; see the wrapper above.
            [varargout{1:max(nargout,1)}] = self.delegate('getAvgCacheTable', max(nargout,1), varargin{:});
        end

        function varargout = getAvgItemTable(self, varargin)
            % ITEMAVGTABLE = GETAVGITEMTABLE()
            % The result recorder captures the returned table together with the solver
            % that produced it -- see LineResultRecorder. Recording an ensemble here
            % rather than in the member solver it delegates to is what keeps an
            % AUTO/LN/ENV/UQ answer from being filed under the member's name.
            [scope, scopeGuard] = LineResultRecorder.enter(); %#ok<ASGLU>
            [varargout{1:max(nargout,1)}] = self.getAvgItemTable_impl(varargin{:});
            LineResultRecorder.capture(scope, self, 'item', varargout{1});
        end

        function varargout = getAvgItemTable_impl(self, varargin)
            % GETAVGITEMTABLE_IMPL Implementation of GETAVGITEMTABLE; see the wrapper above.
            [varargout{1:max(nargout,1)}] = self.delegate('getAvgItemTable', max(nargout,1), varargin{:});
        end

        function varargout = getAvgLossTable(self, varargin)
            % LOSSTABLE = GETAVGLOSSTABLE()
            % The result recorder captures the returned table together with the solver
            % that produced it -- see LineResultRecorder. Recording an ensemble here
            % rather than in the member solver it delegates to is what keeps an
            % AUTO/LN/ENV/UQ answer from being filed under the member's name.
            [scope, scopeGuard] = LineResultRecorder.enter(); %#ok<ASGLU>
            [varargout{1:max(nargout,1)}] = self.getAvgLossTable_impl(varargin{:});
            LineResultRecorder.capture(scope, self, 'loss', varargout{1});
        end

        function varargout = getAvgLossTable_impl(self, varargin)
            % GETAVGLOSSTABLE_IMPL Implementation of GETAVGLOSSTABLE; see the wrapper above.
            [varargout{1:max(nargout,1)}] = self.delegate('getAvgLossTable', max(nargout,1), varargin{:});
        end

        function varargout = getAvgRegionLossTable(self, varargin)
            % LOSSTABLE = GETAVGREGIONLOSSTABLE()
            % The result recorder captures the returned table together with the solver
            % that produced it -- see LineResultRecorder. Recording an ensemble here
            % rather than in the member solver it delegates to is what keeps an
            % AUTO/LN/ENV/UQ answer from being filed under the member's name.
            [scope, scopeGuard] = LineResultRecorder.enter(); %#ok<ASGLU>
            [varargout{1:max(nargout,1)}] = self.getAvgRegionLossTable_impl(varargin{:});
            LineResultRecorder.capture(scope, self, 'regionloss', varargout{1});
        end

        function varargout = getAvgRegionLossTable_impl(self, varargin)
            % GETAVGREGIONLOSSTABLE_IMPL Implementation of GETAVGREGIONLOSSTABLE; see the wrapper above.
            [varargout{1:max(nargout,1)}] = self.delegate('getAvgRegionLossTable', max(nargout,1), varargin{:});
        end

        function varargout = getAvgOrbitTable(self, varargin)
            % ORBITTABLE = GETAVGORBITTABLE()
            % The result recorder captures the returned table together with the solver
            % that produced it -- see LineResultRecorder. Recording an ensemble here
            % rather than in the member solver it delegates to is what keeps an
            % AUTO/LN/ENV/UQ answer from being filed under the member's name.
            [scope, scopeGuard] = LineResultRecorder.enter(); %#ok<ASGLU>
            [varargout{1:max(nargout,1)}] = self.getAvgOrbitTable_impl(varargin{:});
            LineResultRecorder.capture(scope, self, 'orbit', varargout{1});
        end

        function varargout = getAvgOrbitTable_impl(self, varargin)
            % GETAVGORBITTABLE_IMPL Implementation of GETAVGORBITTABLE; see the wrapper above.
            [varargout{1:max(nargout,1)}] = self.delegate('getAvgOrbitTable', max(nargout,1), varargin{:});
        end

        function ON = getAvgOrbit(self)
            % ON = GETAVGORBIT()
            ON = self.delegate('getAvgOrbit', 1);
        end

        function varargout = getAvgTableLayered(self, varargin)
            % [AVGTABLE,QT,UT,RT,WT,AT,TT] = GETAVGTABLELAYERED()
            % The result recorder captures the returned table together with the solver
            % that produced it -- see LineResultRecorder. Recording an ensemble here
            % rather than in the member solver it delegates to is what keeps an
            % AUTO/LN/ENV/UQ answer from being filed under the member's name.
            [scope, scopeGuard] = LineResultRecorder.enter(); %#ok<ASGLU>
            [varargout{1:max(nargout,1)}] = self.getAvgTableLayered_impl(varargin{:});
            LineResultRecorder.capture(scope, self, 'avg', varargout{1});
        end

        function varargout = getAvgTableLayered_impl(self, varargin)
            % GETAVGTABLELAYERED_IMPL Implementation of GETAVGTABLELAYERED; see the wrapper above.
            [varargout{1:max(nargout,1)}] = self.delegate('getAvgTableLayered', max(nargout,1), varargin{:});
        end

        function [QN,UN,RN,WN,AN,TN] = getAvgNodeChain(self,Q,U,R,T)
            % [QN,UN,RN,WN,AN,TN] = GETAVGNODECHAIN(Q,U,R,T)
            if nargin == 1
                [QN,UN,RN,WN,AN,TN] = self.delegate('getAvgNodeChain', 6);
            else
                [QN,UN,RN,WN,AN,TN] = self.delegate('getAvgNodeChain', 6, Q,U,R,T);
            end
        end

        function varargout = getMomentTable(self, varargin)
            % [MOMENTTABLE, MOM] = GETMOMENTTABLE(ORDER)
            [varargout{1:max(nargout,1)}] = self.delegate('getMomentTable', max(nargout,1), varargin{:});
        end

        function varargout = getMomentChainTable(self, varargin)
            % [MOMENTCHAINTABLE, MOM] = GETMOMENTCHAINTABLE(ORDER)
            [varargout{1:max(nargout,1)}] = self.delegate('getMomentChainTable', max(nargout,1), varargin{:});
        end

        function varargout = getMomentStationTable(self, varargin)
            % [MOMENTSTATIONTABLE, MOM] = GETMOMENTSTATIONTABLE(ORDER)
            [varargout{1:max(nargout,1)}] = self.delegate('getMomentStationTable', max(nargout,1), varargin{:});
        end

        function varargout = getSensitivityTable(self, varargin)
            % [SENSTABLE, SENS] = GETSENSITIVITYTABLE(...)
            [varargout{1:max(nargout,1)}] = self.delegate('getSensitivityTable', max(nargout,1), varargin{:});
        end

        % Solver introspection and third-party attribution, resolved on the
        % delegate actually chosen for the model.
        function featSupported = getMethodFeatureSet(self, method)
            % FEATSUPPORTED = GETMETHODFEATURESET(METHOD)
            featSupported = self.delegate('getMethodFeatureSet', 1, method);
        end

        function [bool, reason] = supportsModelMethod(self, method)
            % [BOOL, REASON] = SUPPORTSMODELMETHOD(METHOD)
            [bool, reason] = self.delegate('supportsModelMethod', 2, method);
        end

        function tf = supportsExactSensitivity(self)
            % TF = SUPPORTSEXACTSENSITIVITY()
            tf = self.delegate('supportsExactSensitivity', 1);
        end

        function bool = isStochasticMethod(self, method)
            % BOOL = ISSTOCHASTICMETHOD(METHOD)
            bool = self.delegate('isStochasticMethod', 1, method);
        end

        function libs = libraries(self, options)
            % LIBS = LIBRARIES(OPTIONS)
            if nargin < 2
                libs = self.delegate('libraries', 1);
            else
                libs = self.delegate('libraries', 1, options);
            end
        end

        function libs = getLibrariesUsed(self, options)
            % LIBS = GETLIBRARIESUSED(OPTIONS)
            % Instance-level form: the solver-level twin is static, so it is
            % not reachable through an AUTO object without this wrapper.
            if nargin < 2
                libs = self.delegate('libraries', 1);
            else
                libs = self.delegate('libraries', 1, options);
            end
        end

        function showLibraryAttribution(self, options)
            % SHOWLIBRARYATTRIBUTION(OPTIONS)
            if nargin < 2
                self.delegate('showLibraryAttribution', 0);
            else
                self.delegate('showLibraryAttribution', 0, options);
            end
        end

        function self = initFromSolver(self, initSolver)
            % SELF = INITFROMSOLVER(INITSOLVER)
            % Warm start every candidate: which one runs is decided later.
            for s = 1:length(self.solvers)
                if ~isempty(self.solvers{s}) && ismethod(self.solvers{s}, 'initFromSolver')
                    self.solvers{s}.initFromSolver(initSolver);
                end
            end
        end

        % Aliases for the tables above, mirroring @NetworkSolver.
        function varargout = getAvgCacheT(self, varargin)
            % GETAVGCACHET Alias for getAvgCacheTable
            [varargout{1:nargout}] = self.getAvgCacheTable(varargin{:});
        end

        function varargout = cacheAvgT(self, varargin)
            % CACHEAVGT Alias for getAvgCacheTable
            [varargout{1:nargout}] = self.getAvgCacheTable(varargin{:});
        end

        function varargout = aCaT(self, varargin)
            % ACAT Alias for getAvgCacheTable
            [varargout{1:nargout}] = self.getAvgCacheTable(varargin{:});
        end

        function varargout = getAvgItemT(self, varargin)
            % GETAVGITEMT Alias for getAvgItemTable
            [varargout{1:nargout}] = self.getAvgItemTable(varargin{:});
        end

        function varargout = itemAvgT(self, varargin)
            % ITEMAVGT Alias for getAvgItemTable
            [varargout{1:nargout}] = self.getAvgItemTable(varargin{:});
        end

        function varargout = aIT(self, varargin)
            % AIT Alias for getAvgItemTable
            [varargout{1:nargout}] = self.getAvgItemTable(varargin{:});
        end

        function varargout = getAvgLossT(self, varargin)
            % GETAVGLOSST Alias for getAvgLossTable
            [varargout{1:nargout}] = self.getAvgLossTable(varargin{:});
        end

        function varargout = lossAvgT(self, varargin)
            % LOSSAVGT Alias for getAvgLossTable
            [varargout{1:nargout}] = self.getAvgLossTable(varargin{:});
        end

        function varargout = aLT(self, varargin)
            % ALT Alias for getAvgLossTable
            [varargout{1:nargout}] = self.getAvgLossTable(varargin{:});
        end

        function varargout = getAvgRegionLossT(self, varargin)
            % GETAVGREGIONLOSST Alias for getAvgRegionLossTable
            [varargout{1:nargout}] = self.getAvgRegionLossTable(varargin{:});
        end

        function varargout = regionLossAvgT(self, varargin)
            % REGIONLOSSAVGT Alias for getAvgRegionLossTable
            [varargout{1:nargout}] = self.getAvgRegionLossTable(varargin{:});
        end

        function varargout = aRLT(self, varargin)
            % ARLT Alias for getAvgRegionLossTable
            [varargout{1:nargout}] = self.getAvgRegionLossTable(varargin{:});
        end

        function varargout = getAvgOrbitT(self, varargin)
            % GETAVGORBITT Alias for getAvgOrbitTable
            [varargout{1:nargout}] = self.getAvgOrbitTable(varargin{:});
        end

        function varargout = orbitAvgT(self, varargin)
            % ORBITAVGT Alias for getAvgOrbitTable
            [varargout{1:nargout}] = self.getAvgOrbitTable(varargin{:});
        end

        function varargout = aOT(self, varargin)
            % AOT Alias for getAvgOrbitTable
            [varargout{1:nargout}] = self.getAvgOrbitTable(varargin{:});
        end

        function varargout = getMomentT(self, varargin)
            % GETMOMENTT Alias for getMomentTable
            [varargout{1:nargout}] = self.getMomentTable(varargin{:});
        end

        function varargout = momentT(self, varargin)
            % MOMENTT Alias for getMomentTable
            [varargout{1:nargout}] = self.getMomentTable(varargin{:});
        end

        function varargout = mT(self, varargin)
            % MT Alias for getMomentTable
            [varargout{1:nargout}] = self.getMomentTable(varargin{:});
        end

        function varargout = getMomentChainT(self, varargin)
            % GETMOMENTCHAINT Alias for getMomentChainTable
            [varargout{1:nargout}] = self.getMomentChainTable(varargin{:});
        end

        function varargout = momentChainT(self, varargin)
            % MOMENTCHAINT Alias for getMomentChainTable
            [varargout{1:nargout}] = self.getMomentChainTable(varargin{:});
        end

        function varargout = mCT(self, varargin)
            % MCT Alias for getMomentChainTable
            [varargout{1:nargout}] = self.getMomentChainTable(varargin{:});
        end

        function varargout = getMomentStationT(self, varargin)
            % GETMOMENTSTATIONT Alias for getMomentStationTable
            [varargout{1:nargout}] = self.getMomentStationTable(varargin{:});
        end

        function varargout = momentStationT(self, varargin)
            % MOMENTSTATIONT Alias for getMomentStationTable
            [varargout{1:nargout}] = self.getMomentStationTable(varargin{:});
        end

        function varargout = mST(self, varargin)
            % MST Alias for getMomentStationTable
            [varargout{1:nargout}] = self.getMomentStationTable(varargin{:});
        end

        function varargout = getSensitivityT(self, varargin)
            % GETSENSITIVITYT Alias for getSensitivityTable
            [varargout{1:nargout}] = self.getSensitivityTable(varargin{:});
        end

        function varargout = sensitivityT(self, varargin)
            % SENSITIVITYT Alias for getSensitivityTable
            [varargout{1:nargout}] = self.getSensitivityTable(varargin{:});
        end

        function varargout = sT(self, varargin)
            % ST Alias for getSensitivityTable
            [varargout{1:nargout}] = self.getSensitivityTable(varargin{:});
        end

        % @NetworkSolver short aliases of the average tables. These delegate
        % by their own name: the delegate's alias accepts the full
        % (Q,U,R,T,A,W,keepDisabled) signature that AUTO's getAvgTable does not.
        function varargout = getAvgT(self, varargin)
            % GETAVGT Alias for getAvgTable
            [varargout{1:max(nargout,1)}] = self.delegate('getAvgT', max(nargout,1), varargin{:});
        end

        function varargout = aT(self, varargin)
            % AT Alias for getAvgTable
            [varargout{1:max(nargout,1)}] = self.delegate('aT', max(nargout,1), varargin{:});
        end

        function varargout = getChainAvgT(self, varargin)
            % GETCHAINAVGT Alias for getAvgChainTable
            [varargout{1:max(nargout,1)}] = self.delegate('getChainAvgT', max(nargout,1), varargin{:});
        end

        function varargout = chainAvgT(self, varargin)
            % CHAINAVGT Alias for getAvgChainTable
            [varargout{1:max(nargout,1)}] = self.delegate('chainAvgT', max(nargout,1), varargin{:});
        end

        function varargout = aCT(self, varargin)
            % ACT Alias for getAvgChainTable
            [varargout{1:max(nargout,1)}] = self.delegate('aCT', max(nargout,1), varargin{:});
        end

        function varargout = getNodeAvgT(self, varargin)
            % GETNODEAVGT Alias for getAvgNodeTable
            [varargout{1:max(nargout,1)}] = self.delegate('getNodeAvgT', max(nargout,1), varargin{:});
        end

        function varargout = nodeAvgT(self, varargin)
            % NODEAVGT Alias for getAvgNodeTable
            [varargout{1:max(nargout,1)}] = self.delegate('nodeAvgT', max(nargout,1), varargin{:});
        end

        function varargout = aNT(self, varargin)
            % ANT Alias for getAvgNodeTable
            [varargout{1:max(nargout,1)}] = self.delegate('aNT', max(nargout,1), varargin{:});
        end

        function varargout = getNodeChainAvgT(self, varargin)
            % GETNODECHAINAVGT Alias for getAvgNodeChainTable
            [varargout{1:max(nargout,1)}] = self.delegate('getNodeChainAvgT', max(nargout,1), varargin{:});
        end

        function varargout = nodeChainAvgT(self, varargin)
            % NODECHAINAVGT Alias for getAvgNodeChainTable
            [varargout{1:max(nargout,1)}] = self.delegate('nodeChainAvgT', max(nargout,1), varargin{:});
        end

        function varargout = aNCT(self, varargin)
            % ANCT Alias for getAvgNodeChainTable
            [varargout{1:max(nargout,1)}] = self.delegate('aNCT', max(nargout,1), varargin{:});
        end

        function varargout = getSysAvgT(self, varargin)
            % GETSYSAVGT Alias for getAvgSysTable
            [varargout{1:max(nargout,1)}] = self.delegate('getSysAvgT', max(nargout,1), varargin{:});
        end

        function varargout = sysAvgT(self, varargin)
            % SYSAVGT Alias for getAvgSysTable
            [varargout{1:max(nargout,1)}] = self.delegate('sysAvgT', max(nargout,1), varargin{:});
        end

        function varargout = aST(self, varargin)
            % AST Alias for getAvgSysTable
            [varargout{1:max(nargout,1)}] = self.delegate('aST', max(nargout,1), varargin{:});
        end

    end

    methods (Static)
        function options = defaultOptions()
            % OPTIONS = DEFAULTOPTIONS()
            options = Solver.defaultOptions;
        end

        function delta = optionsDelta(options, baseOptions)
            % DELTA = OPTIONSDELTA(OPTIONS, BASEOPTIONS)
            % The fields of OPTIONS that differ from BASEOPTIONS: what the
            % caller actually asked for, as opposed to what the generic
            % defaults supply. Struct-valued fields (config, odesolvers) are
            % compared one level down, so an untouched config contributes
            % nothing and a candidate's own config survives the overlay.
            delta = struct();
            fn = fieldnames(options);
            for i = 1:numel(fn)
                f = fn{i};
                if ~isfield(baseOptions, f)
                    delta.(f) = options.(f);
                elseif isstruct(options.(f)) && isstruct(baseOptions.(f)) && ...
                        isscalar(options.(f)) && isscalar(baseOptions.(f))
                    sub = SolverAUTO.optionsDelta(options.(f), baseOptions.(f));
                    if ~isempty(fieldnames(sub))
                        delta.(f) = sub;
                    end
                elseif ~isequaln(options.(f), baseOptions.(f))
                    delta.(f) = options.(f);
                end
            end
        end

        function delta = withNamedOptions(delta, options, args)
            % DELTA = WITHNAMEDOPTIONS(DELTA, OPTIONS, ARGS)
            % Add to DELTA every option the caller NAMED in ARGS, whatever its
            % value. OPTIONSDELTA compares values, so an explicit setting that
            % happens to equal the generic default is invisible to it --
            % 'samples',1e4 IS the generic default, while SolverLDES declares
            % 2e5 and SolverNC 1e5 -- and the candidate would keep its own
            % value against an explicit request. An options STRUCT is not
            % inspected here: every field of it would count as named, which is
            % the wholesale override this exists to avoid.
            for i = 1:numel(args)
                a = args{i};
                if ~(ischar(a) || (isstring(a) && isscalar(a)))
                    continue
                end
                name = char(a);
                if startsWith(name, 'config.')
                    sub = name(numel('config.')+1:end);
                    if isfield(options, 'config') && isfield(options.config, sub)
                        if ~isfield(delta, 'config')
                            delta.config = struct();
                        end
                        delta.config.(sub) = options.config.(sub);
                    end
                elseif isfield(options, name)
                    delta.(name) = options.(name);
                end
            end
        end

        function intents = selectionIntents()
            % INTENTS = SELECTIONINTENTS()
            % Tokens that state what the solution is for, rather than naming an
            % algorithm.
            intents = {'default','heur','sim','exact','fast','accurate','bound'};
        end

        function names = familyNames()
            % NAMES = FAMILYNAMES()
            % Method families, in the order in which an unqualified method name
            % is looked up.
            %
            % 'ag' SITS AFTER 'mam', whose RCAT names it took over, and is a
            % family for the METHOD NAME and REPORT tables only: LINE(model,'ag.inap')
            % and model.help() reach SolverAG through it. It is deliberately
            % NOT a candidate of the automatic ranking (the constructor builds
            % that list by hand) and NOT in the feature-set union of
            % getFeatureSet, so a G-network is still refused by 'default' and
            % has to be asked for by name; see _kb/06-solver-catalog.md,
            % "SolverAG owns the RCAT methods".
            names = {'mva','nc','ctmc','fluid','mam','ag','ba','ssa','ldes','jmt','qns','ln','env','lqns','uq'};
        end

        function groups = metricGroups()
            % GROUPS = METRICGROUPS()
            %
            % The measure groups FINDSOLVER reports on, in report order. A
            % group is a family of accessors that stand or fall together: a
            % solver that returns getCdfRespT returns getCdfPassT and
            % getPerctRespT as well, because all three read the same passage
            % time, so listing the three separately would say nothing extra.
            groups = {'avg','tran','cdf','prob','tranprob','sample', ...
                'cache','loss','orbit','moment','sens'};
        end

        function group = metricGroupOf(name)
            % GROUP = METRICGROUPOF(NAME)
            %
            % The measure group an accessor belongs to, '' when NAME names
            % none. A group name maps to itself, so findSolver('cdf') and
            % findSolver('getCdfRespT') ask the same question.
            %
            % THIS IS NOT CHOOSESOLVERHEUR'S TABLE, although both are keyed by
            % accessor name. That one maps an accessor to a RANKING, i.e. which
            % candidate should be preferred; this one maps it to a
            % CAPABILITY question, i.e. which candidates can answer it at all.
            % The two differ wherever a family can serve a measure but is never
            % the one AUTO would pick for it.
            group = '';
            if isempty(name) || ~(ischar(name) || isstring(name))
                return
            end
            name = char(name);
            if any(strcmp(name, SolverAUTO.metricGroups()))
                group = name;
                return
            end
            if any(strcmpi(name, {'any','all',''}))
                return
            end
            switch name
                case {'getTranAvg','getTranAvgVar','tranAvg'}
                    group = 'tran';
                case {'getCdfRespT','getCdfPassT','getPerctRespT', ...
                        'getTranCdfPassT','getTranCdfRespT','getCdfSysRespT'}
                    group = 'cdf';
                case {'getTranProb','getTranProbSys','getTranProbAggr','getTranProbSysAggr'}
                    group = 'tranprob';
                case {'getProb','getProbAggr','getProbSys','getProbSysAggr', ...
                        'getProbMarg','getProbNormConstAggr'}
                    group = 'prob';
                case {'sample','sampleSys','sampleAggr','sampleSysAggr'}
                    group = 'sample';
                case {'getAvgCacheTable','getAvgCacheT','getAvgItemTable','getAvgItemT', ...
                        'cacheAvgT','itemAvgT','aCaT','aIT'}
                    group = 'cache';
                case {'getAvgLossTable','getAvgLossT','getAvgRegionLossTable', ...
                        'getAvgRegionLossT','lossAvgT','regionLossAvgT','aLT','aRLT'}
                    group = 'loss';
                case {'getAvgOrbitTable','getAvgOrbitT','getAvgOrbit','orbitAvgT','aOT'}
                    group = 'orbit';
                case {'getMomentTable','getMomentChainTable','getMomentStationTable', ...
                        'getMomentT','getMomentChainT','getMomentStationT', ...
                        'momentT','momentChainT','momentStationT','mT','mCT','mST'}
                    group = 'moment';
                case {'getSensitivityTable','getSensitivityT','sensitivityT','sT', ...
                        'getSensitivity','getSensitivityRanking'}
                    group = 'sens';
                otherwise
                    % Everything else in the accessor surface is a mean
                    % measure: getAvg, its chain, node and system forms, their
                    % handles and their short aliases.
                    if strncmp(name,'getAvg',6) || strncmp(name,'avg',3) || ...
                            any(strcmp(name,{'getAvgSysRespT','getAvgSysTput','aT','aNT','aCT','aST','aNCT'}))
                        group = 'avg';
                    end
            end
        end

        function groups = familyMetrics(family)
            % GROUPS = FAMILYMETRICS(FAMILY)
            %
            % The measure groups a method family can answer. Every family
            % answers 'avg', which is what a solver is for; the rest is the
            % capability declaration this class owns.
            %
            % SOURCES, so that a claim here can be checked rather than trusted:
            % 'tran' is SUPPORTSTRANSIENTANALYSIS, which FLD, CTMC, LDES and
            % JMT override to true and no one else does. 'cdf', 'prob',
            % 'tranprob' and 'sample' are the families that carry an
            % implementation of the corresponding accessor (@SolverX/getCdfRespT.m
            % and its siblings) rather than inheriting the base refusal. The
            % remaining five groups are computed by NETWORKSOLVER from a
            % solver's own results, so no per-solver file marks them: their
            % lists are CHOOSESOLVERHEUR's rankings for the same accessors,
            % which is where AUTO already records who can serve them.
            %
            % A family that gains or loses a measure must be edited here in the
            % same change, the way a solver that gains a feature is edited into
            % its feature set: an omission here does not fail, it silently
            % hides the family from a caller asking for that measure.
            switch family
                case 'mva'
                    groups = {'avg','prob','cache','orbit','moment','sens'};
                case 'nc'
                    groups = {'avg','cdf','prob','cache','moment','sens'};
                case 'ctmc'
                    groups = {'avg','tran','cdf','prob','tranprob','sample', ...
                        'cache','loss','orbit','moment'};
                case 'fluid'
                    groups = {'avg','tran','cdf','prob','cache','sens'};
                case 'mam'
                    groups = {'avg','cdf'};
                case 'ag'
                    % The RCAT fixed point reports means only; the passage-time
                    % law it answers is the base exponential fit, not its own.
                    groups = {'avg'};
                case 'ba'
                    % A bound brackets the mean measures and nothing else.
                    groups = {'avg'};
                case 'ssa'
                    groups = {'avg','cdf','prob','sample','loss'};
                case 'ldes'
                    groups = {'avg','tran','cdf','prob','sample','cache','loss','orbit'};
                case 'jmt'
                    groups = {'avg','tran','cdf','prob','tranprob','sample'};
                case 'qns'
                    groups = {'avg'};
                case 'ln'
                    groups = {'avg','tran','cdf','sens'};
                case 'env'
                    groups = {'avg','tran'};
                case 'lqns'
                    groups = {'avg'};
                case 'uq'
                    groups = {'avg'};
                otherwise
                    groups = {'avg'};
            end
        end

        function cls = methodClass(family, method, isStochastic, isProductForm, isQbdShape, hasCache)
            % CLS = METHODCLASS(FAMILY, METHOD, ISSTOCHASTIC, ISPRODUCTFORM, ISQBDSHAPE, HASCACHE)
            %
            % What KIND of answer a method returns: 'exact', 'approx', 'bound'
            % or 'simulation'.
            %
            % 'simulation' is not decided here: ISSTOCHASTIC is the solver's own
            % ISSTOCHASTICMETHOD, which already tokenizes qualified and
            % runtime-resolved names and is the only place that knowledge lives.
            %
            % 'exact' IS CLAIMED ONLY WHERE IT IS TRUE OF THIS MODEL, never of
            % the algorithm in the abstract. Exactness of a normalizing constant
            % or of mean value analysis is a property of the product-form
            % model it is computed on, and of the QBD shape for the matrix
            % analytic methods, so both conditions are passed in and a method
            % that needs one reports 'approx' without it. The bias is
            % deliberate: an under-claimed 'approx' costs a user a better
            % method they could have had, an over-claimed 'exact' costs them a
            % wrong number they trusted.
            %
            % A CACHE IS THE THIRD CONDITION, and it was the over-claim the bias
            % above exists to prevent. SN_HAS_PRODUCT_FORM answers about the
            % QUEUEING network and knows nothing of a cache: the hit/miss split
            % is a class switch whose probabilities are not routing data but the
            % output of a cache model, so a network holding one reads as product
            % form and 'mva.exact' was labelled exact on it. Measured on the
            % tut06 shape with an LRU cache: exact MVA returns QLen 0.2516 at
            % the hit station where the CTMC returns 0.3022 and simulation
            % 0.3023, a 17% error under a label that says there is none. The
            % analytic families are conditioned on it; SolverCTMC is NOT,
            % because its state space carries the cache contents and it is exact
            % there, which is what the two numbers above show.
            if nargin < 6
                hasCache = false;
            end
            if strcmp(family,'ba')
                % Bounds are what SolverBA is for; every one of its methods
                % returns a bracket rather than an estimate.
                cls = 'bound';
                return
            end
            if isStochastic
                cls = 'simulation';
                return
            end
            cls = 'approx';
            switch family
                case 'ctmc'
                    % The generator is solved as written, so every state-space
                    % route is exact. 'cftp.approx' says in its own name that
                    % it is not, and 'mdd' is exact on a product-form model and
                    % an approximation otherwise (solver_ctmc_mdd_analyzer).
                    if strcmp(method,'cftp.approx')
                        cls = 'approx';
                    elseif strcmp(method,'mdd')
                        cls = SolverAUTO.exactIf(isProductForm);
                    else
                        cls = 'exact';
                    end
                case 'nc'
                    % The normalizing-constant routes that evaluate G exactly
                    % rather than expanding or estimating it. The asymptotic
                    % expansions (pana, le, kt, bk, gm, ...) and the
                    % non-product-form 'morrison' are approximations by
                    % construction and are left out.
                    if any(strcmp(method,{'exact','divdiff','ca','comom','comomld','rec','ms','cub','rgf'}))
                        cls = SolverAUTO.exactIf(isProductForm && ~hasCache);
                    end
                case 'mva'
                    % Exact MVA; every 'amva.*' arm is an approximation, and so
                    % are the open-network QNA transforms.
                    if any(strcmp(method,{'exact','mva'}))
                        cls = SolverAUTO.exactIf(isProductForm && ~hasCache);
                    end
                case 'jmt'
                    % JMVA's exact algorithms. 'jsim' and 'replication' are
                    % simulation and never reach here.
                    if any(strcmp(method,{'jmva.mva','jmva.recal','jmva.comom','jmva.treeconv'}))
                        cls = SolverAUTO.exactIf(isProductForm);
                    end
                case 'mam'
                    % The QBD is solved exactly on the shape it is stated for,
                    % one queueing station fed by a Source. Everything named
                    % 'dec.*' is a decomposition of a larger network into such
                    % queues and is therefore an approximation of it.
                    if any(strcmp(method,{'default','mna','ldqbd','bgchain','retrial'}))
                        cls = SolverAUTO.exactIf(isQbdShape);
                    end
                case 'ag'
                    % Every RCAT arm estimates the reversed rate of each
                    % synchronising action and iterates to a fixed point, an
                    % approximation by construction; SolverAG's 'exact' is a
                    % vestigial alias that warns and runs 'inap', so nothing
                    % here is claimed exact.
            end
        end

        function cls = exactIf(bool)
            % CLS = EXACTIF(BOOL)
            % 'exact' when the model meets the condition the method's exactness
            % rests on, 'approx' when it does not.
            if bool
                cls = 'exact';
            else
                cls = 'approx';
            end
        end

        function bool = familyAcceptsModelClass(family, model)
            % BOOL = FAMILYACCEPTSMODELCLASS(FAMILY, MODEL)
            % Is FAMILY defined over models of MODEL's kind at all?
            %
            % The composite families are not universal: SolverLN and SolverLQNS
            % analyze a LayeredNetwork and SolverENV an Environment, which is
            % the partition the CONSTRUCTOR already uses when it builds the
            % candidate list. Their feature sets describe what they accept
            % inside a LAYER, so asking one whether it supports a flat Network
            % gets a yes to a question it was never asked: SolverLQNS declares
            % Queue, Exp and SchedStrategy_FCFS, which every closed exponential
            % network uses, and would otherwise be offered for one.
            %
            % THE PARTITION HAS TO BE READ IN BOTH DIRECTIONS, and reading it in
            % one was a bug. Naming the model class each composite family takes
            % says nothing about what the FLAT-NETWORK families take, so every
            % one of them answered yes for an Environment and the report offered
            % seven SolverQNS methods on the tut12 random-environment model, half
            % the table. SolverQNS does not refuse such a model either: its
            % constructor leaves self.model empty rather than raising, so the
            % gate is asked about nothing at all and says yes, and the refusal
            % arrives only at getAvgTable as a bare 'nstations' error.
            %
            % A LayeredNetwork needs no such rule and deliberately does not get
            % one: SolverLDES analyzes an LQN natively and belongs in that
            % report, and the families that do not take one fail to CONSTRUCT
            % over it, which findSolver already treats as "contributes nothing".
            % An Environment is the class where construction does not filter.
            switch family
                case {'ln','lqns'}
                    bool = isa(model,'LayeredNetwork');
                case 'env'
                    bool = isa(model,'Environment');
                otherwise
                    bool = ~isa(model,'Environment');
            end
        end

        function prefixes = methodAliasPrefixes(family)
            % PREFIXES = METHODALIASPREFIXES(FAMILY)
            %
            % The prefixes under which FAMILY advertises a SECOND SPELLING of a
            % method it already declares plainly.
            %
            % THIS IS A DECLARATION, not a derivation, and belongs beside
            % FAMILYMETRICS and METHODCLASS for the same reason: the knowledge
            % lives in the solver's own dispatch (SolverMVA strips a leading
            % 'amva.' before selecting an algorithm) and no accessor exposes it,
            % so a family that gains or loses an alias spelling must be edited
            % into all four copies in the SAME change. An omission does not
            % fail; it puts the same algorithm in the report twice.
            switch family
                case 'mva'
                    prefixes = {'amva.'};
                otherwise
                    prefixes = {};
            end
        end

        function bool = isMethodAlias(family, name, declared)
            % BOOL = ISMETHODALIAS(FAMILY, NAME, DECLARED)
            %
            % Is NAME a second spelling of another method this family declares?
            %
            % The remainder has to be declared too, which is what keeps the rule
            % from eating a genuine method that merely starts with the prefix:
            % it is an alias only when the thing it aliases is there beside it.
            bool = false;
            prefixes = SolverAUTO.methodAliasPrefixes(family);
            for p = 1:length(prefixes)
                prefix = prefixes{p};
                if length(name) > length(prefix) && strncmp(name, prefix, length(prefix)) ...
                        && any(strcmp(name(length(prefix)+1:end), declared))
                    bool = true;
                    return
                end
            end
        end

        function fam = familyAlias(name)
            % FAM = FAMILYALIAS(NAME)
            % Canonical family of a method name, or '' when the method name names none.
            if isempty(name) || ~(ischar(name) || isstring(name))
                fam = '';
                return
            end
            switch lower(char(name))
                case 'mam'
                    fam = 'mam';
                case 'ag'
                    fam = 'ag';
                case 'mva'
                    fam = 'mva';
                case 'nc'
                    fam = 'nc';
                case {'fluid','fld'}
                    fam = 'fluid';
                case 'jmt'
                    fam = 'jmt';
                case 'ssa'
                    fam = 'ssa';
                case 'ctmc'
                    fam = 'ctmc';
                case {'ldes','des'}
                    fam = 'ldes';
                case 'ba'
                    fam = 'ba';
                case 'env'
                    fam = 'env';
                case 'ln'
                    fam = 'ln';
                case {'lqns','lqsim'}
                    fam = 'lqns';
                case 'qns'
                    fam = 'qns';
                case 'uq'
                    fam = 'uq';
                otherwise
                    fam = '';
            end
        end

        function fam = familyOfSolver(solver)
            % FAM = FAMILYOFSOLVER(SOLVER)
            fam = '';
            if isempty(solver)
                return
            end
            switch class(solver)
                case {'SolverMVA','MVA'}
                    fam = 'mva';
                case {'SolverNC','NC'}
                    fam = 'nc';
                case {'SolverCTMC','CTMC'}
                    fam = 'ctmc';
                case {'SolverFluid','SolverFLD','FLD','Fluid'}
                    fam = 'fluid';
                case {'SolverMAM','MAM'}
                    fam = 'mam';
                case {'SolverAG','AG'}
                    fam = 'ag';
                case {'SolverSSA','SSA'}
                    fam = 'ssa';
                case {'SolverLDES','LDES'}
                    fam = 'ldes';
                case {'SolverJMT','JMT'}
                    fam = 'jmt';
                case {'SolverBA','BA'}
                    fam = 'ba';
                case {'SolverLN','LN'}
                    fam = 'ln';
                case {'SolverENV','ENV'}
                    fam = 'env';
                case {'SolverLQNS','LQNS'}
                    fam = 'lqns';
                case {'SolverQNS','QNS'}
                    fam = 'qns';
                case {'SolverUQ','UQ'}
                    fam = 'uq';
                otherwise
                    fam = lower(class(solver));
            end
        end

        function bool = isBoundMethod(method)
            % BOOL = ISBOUNDMETHOD(METHOD)
            % True when the method name names a bound family rather than a point
            % estimator, so that a pinned bound is not reset to the composite.
            bool = false;
            if isempty(method) || ~(ischar(method) || isstring(method))
                return
            end
            method = lower(char(method));
            roots = {'aba','bjb','gb','pb','sb','mwba','pbh','pbk','bjbk','cbh', ...
                'ssd','cub','mbjb','sib','scb','ldbcmp','lr','qr','qrf','auto'};
            head = strtok(method,'.');
            bool = any(strcmp(head, roots));
        end

        function featSupported = getFeatureSet()
            % FEATSUPPORTED = GETFEATURESET()
            % Union of the feature sets of the method families: what the single
            % solver accepts is what at least one of its algorithms accepts.
            featSupported = SolverFeatureSet;
            classes = {'SolverMVA','SolverNC','SolverCTMC','SolverFluid','SolverMAM', ...
                'SolverSSA','SolverLDES','SolverJMT'};
            fields = fieldnames(featSupported.list);
            for c = 1:length(classes)
                try
                    fs = feval([classes{c},'.getFeatureSet']);
                catch
                    continue
                end
                for f = 1:length(fields)
                    if isfield(fs.list, fields{f}) && fs.list.(fields{f})
                        featSupported.setTrue(fields{f});
                    end
                end
            end
        end

        function [bool, featSupported, featUsed] = supports(model)
            % [BOOL, FEATSUPPORTED, FEATUSED] = SUPPORTS(MODEL)
            % True when at least one family can analyze the model. Composite
            % models are delegated to the family that decomposes them.
            featSupported = [];
            featUsed = [];
            switch class(model)
                case 'LayeredNetwork'
                    bool = false;
                    try
                        bool = SolverLN.supports(model);
                    catch
                        bool = true;
                    end
                    return
                case 'Environment'
                    try
                        [bool, featSupported] = SolverENV.supports(model);
                    catch
                        bool = true;
                    end
                    return
            end
            featUsed = model.getUsedLangFeatures();
            featSupported = SolverAUTO.getFeatureSet();
            bool = SolverFeatureSet.supports(featSupported, featUsed);
        end

        function tf = hasPriorParameters(model)
            % TF = HASPRIORPARAMETERS(MODEL)
            % True when the model carries uncertain parameters, which expand
            % into one model instance per alternative.
            %
            % Asked of UQ, which owns the question: UQ.supports gates on the
            % same answer, so the constructor's routing here and the gate there
            % cannot drift apart.
            tf = UQ.modelHasPrior(model);
        end

        function [ok, msg, logNstates] = isStateSpaceTractable(model, options)
            % [OK, MSG, LOGNSTATES] = ISSTATESPACETRACTABLE(MODEL, OPTIONS)
            % Same gate the exact analysis runs, exposed so that a caller can
            % rank state enumeration out before paying for it.
            if nargin < 2
                [ok, msg, logNstates] = SolverCTMC.isStateSpaceTractable(model);
            else
                [ok, msg, logNstates] = SolverCTMC.isStateSpaceTractable(model, options);
            end
        end

        function printInfGen(Q, SS)
            % PRINTINFGEN(Q, SS)
            SolverCTMC.printInfGen(Q, SS);
        end
    end
end
