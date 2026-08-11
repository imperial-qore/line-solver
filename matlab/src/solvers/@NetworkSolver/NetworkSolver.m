classdef NetworkSolver < Solver
    % NetworkSolver Abstract base class for queueing network solvers
    %
    % NetworkSolver provides the common interface and functionality for all
    % solvers that can analyze Network models. It handles performance metric
    % computation, language switching between MATLAB and Java implementations,
    % and provides standardized solver initialization and execution patterns.
    %
    % @brief Abstract base class for all queueing network analysis solvers
    %
    % Key characteristics:
    % - Abstract base for all network-specific solvers
    % - Manages performance metric handles and computation
    % - Supports both MATLAB native and Java/JLINE implementations
    % - Provides standardized solver options and configuration
    % - Handles model language switching and delegation
    %
    % NetworkSolver serves as the foundation for:
    % - Analytical solvers (MVA, NC, CTMC, etc.)
    % - Simulation solvers (JMT, SSA)
    % - Approximation methods (Fluid, MAM, NN)
    % - Automatic solver selection (AUTO)
    %
    % Example usage pattern:
    % @code
    % solver = SolverMVA(model, 'MyMVASolver');
    % solver.getAvg();  % Get average performance metrics
    % @endcode
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.


    properties (Access = protected)
        handles; % performance metric handles
    end

    properties (Hidden)
        % Cached MMT fork-join transformation, reused across SolverLN outer
        % iterations (tagged by structVersion, NOT cleared by reset()).
        % see _kb/05-solvers-overview.md for rationale
        mmtCache = [];
    end

    properties (Hidden)
        % Percentile extraction method of the last getPerctRespT call
        % ('cdf', 'forktail', ...), so that citations() can report the paper
        % behind a tail estimate the user actually asked for.
        lastPerctMethod = '';
    end

    properties (Access = public)
        % Auxiliary-class arrival rates of the fork-join (MMT) fixed point,
        % retained across runAnalyzer calls so that an outer iteration (e.g.
        % SolverLN) restarts the MMT loop from the previous converged point
        % instead of from GlobalConstants.FineTol. Like options.init_sol it
        % survives reset(); use resetForkWarmStart to discard it.
        fjForkLambda = [];
    end

    methods

        function self = NetworkSolver(model, name, options)
            % NETWORKSOLVER Create a NetworkSolver instance
            %
            % @brief Creates a NetworkSolver with the specified model and options
            % @param model Network model to be analyzed
            % @param name String identifier for the solver instance
            % @param options Optional SolverOptions structure for configuration
            % @return self NetworkSolver instance ready for analysis
            %
            % The constructor initializes the solver with the provided model,
            % validates inputs, sets options, and prepares performance metric handles.
            self@Solver(model, name);
            if isempty(model)
                line_error(mfilename,'The model supplied in input is empty');
            end
            if nargin>=3 %exist('options','var'),
                self.setOptions(options);
            end
            self.result = [];

            if isempty(model.obj) % not a Java object
                initHandles(self)
            end
        end

        function resetForkWarmStart(self)
            % RESETFORKWARMSTART()
            % Discard the retained MMT fixed point. Mirrors options.init_sol:
            % the iterate survives reset() and is invalidated explicitly by the
            % caller when the chain basis changes.
            self.fjForkLambda = [];
        end

        function setLang(self)
            % NOTE: self.model is passed by reference so it affects the
            % user model content
            switch self.options.lang %
                case 'matlab' % matlab solver
                    if self.model.isJavaNative() % java model
                        matlab_model = JLINE.jline_to_line(self.model.obj);
                        matlab_model.obj = self.model.obj;
                        self.model = matlab_model;
                        self.initHandles;
                    else  % matlab model
                        % no-op
                    end
                case 'java' % solver lang
                    joptions = self.options;
                    if ~isempty(self.model.obj) % java model
                        % no-op
                    else  % matlab model
                        self.model.obj = JLINE.line_to_jline(self.model);
                    end
                    switch self.name
                        case 'SolverCTMC'
                            self.obj = JLINE.SolverCTMC(self.model.obj, joptions);
                        case 'SolverLDES'
                            self.obj = JLINE.SolverLDES(self.model.obj, joptions);
                        case {'SolverFluid','SolverFLD'}
                            self.obj = JLINE.SolverFluid(self.model.obj, joptions);
                        case 'SolverJMT'
                            self.obj = JLINE.SolverJMT(self.model.obj, joptions);
                        case 'SolverMAM'
                            self.obj = JLINE.SolverMAM(self.model.obj, joptions);
                        case 'SolverMVA'
                            self.obj = JLINE.SolverMVA(self.model.obj, joptions);
                        case 'SolverNC'
                            self.obj = JLINE.SolverNC(self.model.obj, joptions);
                        case 'SolverSSA'
                            self.obj = JLINE.SolverSSA(self.model.obj, joptions);
                        case 'SolverQNS'
                            self.obj = JLINE.SolverQNS(self.model.obj, joptions);
                    end
            end
        end

        function initHandles(self)
            [Q,U,R,T,A,W] = self.model.getAvgHandles();

            % Get tardiness handles if available
            if ismethod(self.model, 'getAvgTardHandles')
                Tard = self.model.getAvgTardHandles();
            else
                Tard = [];
            end
            if ismethod(self.model, 'getAvgSysTardHandles')
                SysTard = self.model.getAvgSysTardHandles();
            else
                SysTard = [];
            end

            self.setAvgHandles(Q,U,R,T,A,W,Tard,SysTard);

            [Qt,Ut,Tt] = self.model.getTranHandles;
            self.setTranHandles(Qt,Ut,Tt);

            % refreshStruct is called lazily by model.getStruct() when first needed
        end

        function self = runAnalyzerChecks(self, options)
            % Single, method-aware feature gate shared by every solver.
            % It resolves the concrete method that will actually run (for
            % options.method='default' this may map to a specific method via
            % feature-driven selection), then gates the model against that
            % method's per-method feature set rather than the solver-level
            % union. Solvers whose methods are all equivalent inherit the base
            % behavior transparently (getMethodFeatureSet defaults to the
            % coarse solver set), so this replaces the previous bespoke
            % per-solver overrides without changing their behavior.

            % Propagate solver verbose level to global so that model-level
            % messages (e.g., priority info in refreshStruct) respect it
            GlobalConstants.setVerbose(options.verbose);
            if ~self.enableChecks
                return;
            end
            if ~any(cellfun(@(s) strcmp(s,options.method),self.listValidMethods))
                line_error(mfilename,sprintf('The ''%s'' method is unsupported by this solver.\n',options.method));
            end
            method = self.resolveMethod(options);
            [bool, reason] = self.supportsModelMethod(method);
            if ~bool
                if strcmp(method, options.method)
                    line_error(mfilename, sprintf('This model contains features not supported by the solver. %s', reason));
                else
                    line_error(mfilename, sprintf('This model contains features not supported by the solver''s ''%s'' method. %s', method, reason));
                end
            end
        end

        function [bool, reason] = supportsModelMethod(self, method)
            % Fine, method-aware gate: does the model fit the concrete METHOD?
            % Base behavior derives the answer from getMethodFeatureSet(method):
            % when that is empty the solver does not diverge per method and the
            % solver's own supports(model) is used (preserving structural checks
            % such as LQNS layers / LDES LayeredNetwork special-casing).
            % Solvers with non-feature-set structural per-method rules (e.g. NC
            % 'mem' via solver_nc_mem_supports) override this.
            featSupported = self.getMethodFeatureSet(method);
            if isempty(featSupported)
                bool = self.supports(self.model);
                reason = '';
            else
                featUsed = self.model.getUsedLangFeatures();
                [bool, reason] = SolverFeatureSet.supports(featSupported, featUsed);
            end
        end

        function method = resolveMethod(self, options)
            % Resolve the concrete method that will run. Base behavior is a
            % no-op (returns options.method unchanged). Solvers that perform
            % feature-driven selection for options.method='default' override
            % this (typically via selectMethod).
            method = options.method;
        end

        function featSupported = getMethodFeatureSet(self, method) %#ok<INUSD>
            % Per-method feature set. Returning empty ([]) signals "this solver
            % does not diverge per method": the gate then uses the solver's own
            % supports(model), preserving any structural checks it carries.
            % Divergent solvers (e.g. MVA, MAM, NC) override this to return the
            % base envelope with per-method add/remove deltas applied.
            featSupported = [];
        end

        function method = selectMethod(self, preferenceList)
            % Feature-driven method selection: returns the first method in
            % PREFERENCELIST whose per-method feature set covers the model.
            % Falls back to the last entry when none fully covers (the base
            % gate then reports the precise unsupported features).
            method = preferenceList{end};
            for k = 1:numel(preferenceList)
                cand = preferenceList{k};
                [ok, ~] = self.supportsModelMethod(cand);
                if ok
                    method = cand;
                    return;
                end
            end
        end

        function self = setTranHandles(self,Qt,Ut,Tt)
            self.handles.Qt = Qt;
            self.handles.Ut = Ut;
            self.handles.Tt = Tt;
        end

        function self = setAvgHandles(self,Q,U,R,T,A,W,Tard,SysTard)
            self.handles.Q = Q;
            self.handles.U = U;
            self.handles.R = R;
            self.handles.T = T;
            self.handles.A = A;
            self.handles.W = W;
            if nargin >= 8
                self.handles.Tard = Tard;
            end
            if nargin >= 9
                self.handles.SysTard = SysTard;
            end
        end

        function [Qt,Ut,Tt] = getTranHandles(self)
            Qt = self.handles.Qt;
            Ut = self.handles.Ut;
            Tt = self.handles.Tt;
        end

        function [Q,U,R,T,A,W] = getAvgHandles(self)
            if isempty(self.handles)
                self.handles.Q = [];
                self.handles.U = [];
                self.handles.R = [];
                self.handles.T = [];
                self.handles.A = [];
                self.handles.W = [];
                self.handles.Tard = [];
                self.handles.SysTard = [];
                initHandles(self);
            end
            Q = self.handles.Q;
            U = self.handles.U;
            R = self.handles.R;
            T = self.handles.T;
            A = self.handles.A;
            W = self.handles.W;
        end

        function Q = getAvgQLenHandles(self)
            if isempty(self.handles) || ~isstruct(self.handles)
                self.getAvgHandles();
            end
            Q = self.handles.Q;
        end

        function U = getAvgUtilHandles(self)
            if isempty(self.handles) || ~isstruct(self.handles)
                self.getAvgHandles();
            end
            U = self.handles.U;
        end

        function R = getAvgRespTHandles(self)
            if isempty(self.handles) || ~isstruct(self.handles)
                self.getAvgHandles();
            end
            R = self.handles.R;
        end

        function T = getAvgTputHandles(self)
            if isempty(self.handles) || ~isstruct(self.handles)
                self.getAvgHandles();
            end
            T = self.handles.T;
        end

        function A = getAvgArvRHandles(self)
            if isempty(self.handles) || ~isstruct(self.handles)
                self.getAvgHandles();
            end
            A = self.handles.A;
        end

        function W = getAvgResidTHandles(self)
            if isempty(self.handles) || ~isstruct(self.handles)
                self.getAvgHandles();
            end
            W = self.handles.W;
        end
    end


    methods (Access = 'protected')
        function bool = hasAvgResults(self)
            % BOOL = HASAVGRESULTS()

            % Returns true if the solver has computed steady-state average metrics.
            bool = false;
            if self.hasResults
                if isfield(self.result,'Avg')
                    bool = true;
                end
            end
        end

        function bool = hasTranResults(self)
            % BOOL = HASTRANRESULTS()

            % Return true if the solver has computed transient average metrics.
            bool = false;
            if self.hasResults
                if isfield(self.result,'Tran')
                    if isfield(self.result.Tran,'Avg')
                        bool = isfield(self.result.Tran.Avg,'Q');
                    end
                end
            end
        end

        function bool = hasDistribResults(self)
            % BOOL = HASDISTRIBRESULTS()

            % Return true if the solver has computed steady-state distribution metrics.
            bool = false;
            if self.hasResults
                bool = isfield(self.result.Distribution,'C');
            end
        end
    end

    methods (Sealed)

        function self = setModel(self, model)
            % SELF = SETMODEL(MODEL)

            % Assign the model to be solved.
            self.model = model;
        end

        function QN = getAvgQLen(self)
            % QN = GETAVGQLEN()

            % Compute average queue-lengths at steady-state
            Q = getAvgQLenHandles(self);
            [QN,~,~,~] = self.getAvg(Q,[],[],[],[],[]);
        end

        function UN = getAvgUtil(self)
            % UN = GETAVGUTIL()

            % Compute average utilizations at steady-state
            U = getAvgUtilHandles(self);
            [~,UN,~,~] = self.getAvg([],U,[],[],[],[]);
        end

        function RN = getAvgRespT(self)
            % RN = GETAVGRESPT()

            % Compute average response times at steady-state
            R = getAvgRespTHandles(self);
            [~,~,RN,~] = self.getAvg([],[],R,[],[],[]);
        end

        function WN = getAvgResidT(self)
            % WN = GETAVGRESIDT()

            % Compute average residence times at steady-state
            R = getAvgRespTHandles(self);
            W = getAvgResidTHandles(self);
            [~,~,~,~,~,WN] = self.getAvg([],[],R,[],[],W);
        end

        function WT = getAvgWaitT(self)
            % RN = GETAVGWAITT()
            % Compute average waiting time in queue excluding service
            R = getAvgRespTHandles(self);
            [~,~,RN,~] = self.getAvg([],[],R,[],[],[]);
            if isempty(RN)
                WT = [];
                return
            end
            sn = self.model.getStruct;
            WT = RN - 1./ sn.rates(:);
            WT(sn.nodetype==NodeType.Source) = 0;
        end

        function TN = getAvgTput(self)
            % TN = GETAVGTPUT()

            % Compute average throughputs at steady-state
            T = getAvgTputHandles(self);
            [~,~,~,TN] = self.getAvg([],[],[],T,[],[]);
        end

        function AN = getAvgArvR(self)
            % AN = GETAVGARVR()
            sn = self.model.getStruct();

            % Compute average arrival rate at steady-state
            TH = getAvgTputHandles(self);
            [~,~,~,TN] = self.getAvg([],[],[],TH,[],[]);
            AN = sn_get_arvr_from_tput(sn, TN, TH);
        end

        % also accepts a cell array with the handlers in it
        [QN,UN,RN,TN,AN,WN]       = getAvg(self,Q,U,R,T,A,W);
        [QN,UN,RN,TN,AN,WN]       = getAvgNode(self,Q,U,R,T,A,W);

        [AvgTable,QT,UT,RT,WT,AT,TT] = getAvgTable(self,Q,U,R,T,A,keepDisabled);

        [AvgTable,QT] = getAvgQLenTable(self,Q,keepDisabled);
        [AvgTable,UT] = getAvgUtilTable(self,U,keepDisabled);
        [AvgTable,RT] = getAvgRespTTable(self,R,keepDisabled);
        [AvgTable,TT] = getAvgTputTable(self,T,keepDisabled);

        [NodeAvgTable,QTn,UTn,RTn,WTn,ATn,TTn] = getAvgNodeTable(self,Q,U,R,T,A,keepDisabled);
        [CacheAvgTable] = getAvgCacheTable(self);
        [ItemAvgTable] = getAvgItemTable(self);
        [AvgChain,QTc,UTc,RTc,WTc,ATc,TTc] = getAvgChainTable(self,Q,U,R,T);
        [AvgChain,QTc,UTc,RTc,WTc,ATc,TTc] = getAvgNodeChainTable(self,Q,U,R,T);

        [QNc,UNc,RNc,WNc,ANc,TNc]   = getAvgChain(self,Q,U,R,T);
        [QNc]                       = getAvgQLenChain(self,Q);
        [UNc]                       = getAvgUtilChain(self,U);
        [RNc]                       = getAvgRespTChain(self,R);
        [WNc]                       = getAvgResidTChain(self,W);
        [TNc]                       = getAvgTputChain(self,T);
        [ANc]                       = getAvgArvRChain(self,A);
        [QNc]                       = getAvgNodeQLenChain(self,Q);
        [UNc]                       = getAvgNodeUtilChain(self,U);
        [RNc]                       = getAvgNodeRespTChain(self,R);
        [WNc]                       = getAvgNodeResidTChain(self,W);
        [TNc]                       = getAvgNodeTputChain(self,T);
        [ANc]                       = getAvgNodeArvRChain(sef,A);

        [CNc,XNc]           = getAvgSys(self,R,T);
        [CT,XT]             = getAvgSysTable(self,R,T);
        [RN]                = getAvgSysRespT(self,R);
        [TN]                = getAvgSysTput(self,T);


        function self = setAvgResults(self,Q,U,R,T,A,W,C,X,runtime,method,iter)
            % SELF = SETAVGRESULTS(SELF,Q,U,R,T,A,W,C,X,RUNTIME,METHOD,ITER)
            % Store average metrics at steady-state
            self.result.('solver') = getName(self);
            if nargin<11 %~exist('method','var')
                method = getOptions(self).method;
            end
            if nargin<12 %~exist('iter','var')
                iter = NaN;
            end
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
                % Wall-clock time-budget flag; solvers that stop early on
                % options.timeout overwrite this with true after calling setAvgResults.
                self.result.Avg.timedOut = false;
            end
            if getOptions(self).verbose
                try
                    solvername = erase(self.result.solver,'Solver');
                catch
                    solvername = self.result.solver(7:end);
                end
                if isnan(iter) || iter==1 || strcmp(solvername,'LDES') || strcmp(solvername,'SSA')
                    line_printf('%s analysis [method: %s, lang: %s, env: %s] completed in %fs.',solvername,self.result.Avg.method,self.options.lang,version("-release"),runtime);
                else
                    line_printf('%s analysis [method: %s, lang: %s, env: %s] completed in %fs. Iterations: %d.',solvername,self.result.Avg.method,self.options.lang,version("-release"),runtime,iter);
                end
                line_printf('\n');
            end
        end

        function self = setAvgResultsCI(self, QCI, UCI, RCI, TCI, ACI, WCI, CCI, XCI)
            % SELF = SETAVGRESULTSCI(SELF, QCI, UCI, RCI, TCI, ACI, WCI, CCI, XCI)
            % Store confidence interval bounds for average metrics
            % Each CI parameter is an [M x K x 2] array with lower/upper bounds
            % or an [M x K] array with half-widths (mean ± halfwidth)
            if nargin >= 2 && ~isempty(QCI)
                self.result.Avg.QCI = QCI;
            end
            if nargin >= 3 && ~isempty(UCI)
                self.result.Avg.UCI = UCI;
            end
            if nargin >= 4 && ~isempty(RCI)
                self.result.Avg.RCI = RCI;
            end
            if nargin >= 5 && ~isempty(TCI)
                self.result.Avg.TCI = TCI;
            end
            if nargin >= 6 && ~isempty(ACI)
                self.result.Avg.ACI = ACI;
            end
            if nargin >= 7 && ~isempty(WCI)
                self.result.Avg.WCI = WCI;
            end
            if nargin >= 8 && ~isempty(CCI)
                self.result.Avg.CCI = CCI;
            end
            if nargin >= 9 && ~isempty(XCI)
                self.result.Avg.XCI = XCI;
            end
        end

        function self = setDistribResults(self,Cd,runtime)
            % SELF = SETDISTRIBRESULTS(SELF,CD,RUNTIME)

            % Store distribution metrics at steady-state
            self.result.('solver') = getName(self);
            self.result.Distribution.('method') = getOptions(self).method;
            self.result.Distribution.C = Cd;
            self.result.Distribution.runtime = runtime;
        end

        function self = setTranProb(self,t,pi_t,SS,runtimet)
            % SELF = SETTRANPROB(SELF,T,PI_T,SS,RUNTIMET)

            % Store transient average metrics
            self.result.('solver') = getName(self);
            self.result.Tran.Prob.('method') = getOptions(self).method;
            self.result.Tran.Prob.t = t;
            self.result.Tran.Prob.pi_t = pi_t;
            self.result.Tran.Prob.SS = SS;
            self.result.Tran.Prob.runtime = runtimet;
        end

        function self = setTranAvgResults(self,Qt,Ut,Rt,Tt,Ct,Xt,runtimet)
            % SELF = SETTRANAVGRESULTS(SELF,QT,UT,RT,TT,CT,XT,RUNTIMET)

            % Store transient average metrics
            self.result.('solver') = getName(self);
            self.result.Tran.Avg.('method') = getOptions(self).method;
            % Clear individual cells that are scalar NaN (not entire array)
            for i=1:size(Qt,1), for r=1:size(Qt,2), if isscalar(Qt{i,r}) && any(isnan(Qt{i,r})), Qt{i,r}=[]; end, end, end
            for i=1:size(Rt,1), for r=1:size(Rt,2), if isscalar(Rt{i,r}) && any(isnan(Rt{i,r})), Rt{i,r}=[]; end, end, end
            for i=1:size(Ut,1), for r=1:size(Ut,2), if isscalar(Ut{i,r}) && any(isnan(Ut{i,r})), Ut{i,r}=[]; end, end, end
            for i=1:size(Tt,1), for r=1:size(Tt,2), if isscalar(Tt{i,r}) && any(isnan(Tt{i,r})), Tt{i,r}=[]; end, end, end
            for i=1:size(Xt,1), for r=1:size(Xt,2), if isscalar(Xt{i,r}) && any(isnan(Xt{i,r})), Xt{i,r}=[]; end, end, end
            for i=1:size(Ct,1), for r=1:size(Ct,2), if isscalar(Ct{i,r}) && any(isnan(Ct{i,r})), Ct{i,r}=[]; end, end, end
            self.result.Tran.Avg.Q = Qt;
            self.result.Tran.Avg.R = Rt;
            self.result.Tran.Avg.U = Ut;
            self.result.Tran.Avg.T = Tt;
            self.result.Tran.Avg.X = Xt;
            self.result.Tran.Avg.C = Ct;
            self.result.Tran.Avg.runtime = runtimet;
        end

    end

    methods
        [QNt,UNt,TNt] = getTranAvg(self,Qt,Ut,Tt);

        function [lNormConst] = getProbNormConstAggr(self)
            % [LNORMCONST] = GETPROBNORMCONST()

            % Return normalizing constant of state probabilities
            line_error(mfilename,sprintf('getProbNormConstAggr is not supported by %s',class(self)));
        end

        function Pstate = getProb(self, node, state)
            % PSTATE = GETPROBSTATE(NODE, STATE)

            % Return marginal state probability for station ist state
            line_error(mfilename,sprintf('getProb is not supported by %s',class(self)));
        end

        function Psysstate = getProbSys(self)
            % PSYSSTATE = GETPROBSYSSTATE()

            % Return joint state probability
            line_error(mfilename,sprintf('getProbSys is not supported by %s',class(self)));
        end

        function Pnir = getProbAggr(self, node, state_a)
            % PNIR = GETPROBSTATEAGGR(NODE, STATE_A)

            % Return marginal state probability for station ist state
            line_error(mfilename,sprintf('getProbAggr is not supported by %s',class(self)));
        end

        function Pnjoint = getProbSysAggr(self)
            % PNJOINT = GETPROBSYSSTATEAGGR()

            % Return joint state probability
            line_error(mfilename,sprintf('getProbSysAggr is not supported by %s',class(self)));
        end

        function Pmarg = getProbMarg(self, node, jobclass, state_m)
            % PMARG = GETPROBMARG(NODE, JOBCLASS, STATE_M)

            % Return marginalized state probability for station and class
            % This computes probabilities marginalized on a given class, in contrast to
            % getProbAggr which aggregates over all classes.
            line_error(mfilename,sprintf('getProbMarg is not supported by %s',class(self)));
        end

        function tstate = sample(self, node, numEvents)
            % TSTATE = SAMPLE(NODE, numEvents)

            % Return marginal state probability for station ist state
            line_error(mfilename,sprintf('sample is not supported by %s',class(self)));
        end

        function tstate = sampleAggr(self, node, numEvents)
            % TSTATE = SAMPLEAGGR(NODE, numEvents)

            % Return marginal state probability for station ist state
            line_error(mfilename,sprintf('sampleAggr is not supported by %s',class(self)));
        end

        function tstate = sampleSys(self, numEvents)
            % TSTATE = SAMPLESYS(numEvents)

            % Return joint state probability
            line_error(mfilename,sprintf('sampleSys is not supported by %s',class(self)));
        end

        function tstate = sampleSysAggr(self, numEvents)
            % TSTATE = SAMPLESYSAGGR(numEvents)

            % Return joint state probability
            line_error(mfilename,sprintf('sampleSysAggr is not supported by %s',class(self)));
        end

        function RD = getCdfRespT(self, R)
            % RD = GETCDFRESPT(R)

            % Return cumulative distribution of response times at steady-state
            % This uses a trivial approximation that assumes exponential
            % distributions everywhere with mean as RN(i,r)

            sn = self.getStruct;
            RD = cell(sn.nstations,sn.nclasses);
            if GlobalConstants.DummyMode
                return
            end
            T0 = tic;
            if nargin<2 %~exist('R','var')
                R = self.getAvgRespTHandles;
                % to do: check if some R are disabled
            end
            if ~self.hasAvgResults
                self.getAvg; % get steady-state solution
            end
            for i=1:sn.nstations
                if sn.nodetype(sn.stationToNode(i)) ~= NodeType.Source
                    for c=1:sn.nclasses
                        if isfinite(self.result.Avg.R(i,c)) && self.result.Avg.R(i,c)>0
                            lambda = 1/self.result.Avg.R(i,c);
                            n = 100; % number of points
                            quantiles = linspace(0.001, 0.999, n);
                            RD{i,c} = [quantiles;-log(1 - quantiles) / lambda]';
                        else
                            RD{i,c} = [1,0];
                        end
                    end
                end
            end
            runtime = toc(T0);
            self.setDistribResults(RD, runtime);
        end

        function RD = getTranCdfRespT(self, R)
            % RD = GETTRANCDFRESPT(R)

            % Return cumulative distribution of response times during transient
            line_error(mfilename,sprintf('getTranCdfRespT is not supported by %s',class(self)));
        end

        function RD = getCdfPassT(self, R)
            % RD = GETCDFPASST(R)

            % Return cumulative distribution of passage times at steady-state
            line_error(mfilename,sprintf('getCdfPassT is not supported by %s',class(self)));
        end

        function RD = getTranCdfPassT(self, R)
            % RD = GETTRANCDFPASST(R)

            % Return cumulative distribution of passage times during transient
            line_error(mfilename,sprintf('getTranCdfPassT is not supported by %s',class(self)));
        end
        
        % Kotlin-style aliases for getAvg* methods
        function avg_table = avgTable(self)
            % AVGTABLE Kotlin-style alias for getAvgTable
            avg_table = self.getAvgTable();
        end
        
        function avg_sys_table = avgSysTable(self)
            % AVGSYSTABLE Kotlin-style alias for getAvgSysTable
            avg_sys_table = self.getAvgSysTable();
        end
        
        function avg_node_table = avgNodeTable(self)
            % AVGNODETABLE Kotlin-style alias for getAvgNodeTable
            avg_node_table = self.getAvgNodeTable();
        end
        
        function avg_chain_table = avgChainTable(self)
            % AVGCHAINTABLE Kotlin-style alias for getAvgChainTable
            avg_chain_table = self.getAvgChainTable();
        end
        
        function avg_node_chain_table = avgNodeChainTable(self)
            % AVGNODECHAINTABLE Kotlin-style alias for getAvgNodeChainTable
            avg_node_chain_table = self.getAvgNodeChainTable();
        end

        % Table -> T aliases
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

        function varargout = avgChain(self, varargin)
            % AVGCHAIN Kotlin-style alias for getAvgChain
            [varargout{1:nargout}] = self.getAvgChain(varargin{:});
        end
        
        function varargout = avgSys(self, varargin)
            % AVGSYS Kotlin-style alias for getAvgSys
            [varargout{1:nargout}] = self.getAvgSys(varargin{:});
        end
        
        function varargout = avgNode(self, varargin)
            % AVGNODE Kotlin-style alias for getAvgNode
            [varargout{1:nargout}] = self.getAvgNode(varargin{:});
        end
        
        function sys_resp_time = avgSysRespT(self, varargin)
            % AVGSYSRESPT Kotlin-style alias for getAvgSysRespT
            sys_resp_time = self.getAvgSysRespT(varargin{:});
        end
        
        function sys_tput = avgSysTput(self, varargin)
            % AVGSYSTPUT Kotlin-style alias for getAvgSysTput
            sys_tput = self.getAvgSysTput(varargin{:});
        end
        
        function arvr_chain = avgArvRChain(self, varargin)
            % AVGARVCHAIN Kotlin-style alias for getAvgArvRChain
            arvr_chain = self.getAvgArvRChain(varargin{:});
        end
        
        function qlen_chain = avgQLenChain(self, varargin)
            % AVGQLENCHAIN Kotlin-style alias for getAvgQLenChain
            qlen_chain = self.getAvgQLenChain(varargin{:});
        end
        
        function util_chain = avgUtilChain(self, varargin)
            % AVGUTILCHAIN Kotlin-style alias for getAvgUtilChain
            util_chain = self.getAvgUtilChain(varargin{:});
        end
        
        function resp_t_chain = avgRespTChain(self, varargin)
            % AVGRESPTCHAIN Kotlin-style alias for getAvgRespTChain
            resp_t_chain = self.getAvgRespTChain(varargin{:});
        end
        
        function resid_t_chain = avgResidTChain(self, varargin)
            % AVGRESIDTCHAIN Kotlin-style alias for getAvgResidTChain
            resid_t_chain = self.getAvgResidTChain(varargin{:});
        end
        
        function tput_chain = avgTputChain(self, varargin)
            % AVGTPUTCHAIN Kotlin-style alias for getAvgTputChain
            tput_chain = self.getAvgTputChain(varargin{:});
        end
        
        function node_arvr_chain = avgNodeArvRChain(self, varargin)
            % AVGNODERVRCHAIN Kotlin-style alias for getAvgNodeArvRChain
            node_arvr_chain = self.getAvgNodeArvRChain(varargin{:});
        end
        
        function node_qlen_chain = avgNodeQLenChain(self, varargin)
            % AVGNODEQLENCHAIN Kotlin-style alias for getAvgNodeQLenChain
            node_qlen_chain = self.getAvgNodeQLenChain(varargin{:});
        end
        
        function node_util_chain = avgNodeUtilChain(self, varargin)
            % AVGNODEUTILCHAIN Kotlin-style alias for getAvgNodeUtilChain
            node_util_chain = self.getAvgNodeUtilChain(varargin{:});
        end
        
        function node_resp_t_chain = avgNodeRespTChain(self, varargin)
            % AVGNODERESPTCHAIN Kotlin-style alias for getAvgNodeRespTChain
            node_resp_t_chain = self.getAvgNodeRespTChain(varargin{:});
        end
        
        function node_resid_t_chain = avgNodeResidTChain(self, varargin)
            % AVGNODERESIDTCHAIN Kotlin-style alias for getAvgNodeResidTChain
            node_resid_t_chain = self.getAvgNodeResidTChain(varargin{:});
        end
        
        function node_tput_chain = avgNodeTputChain(self, varargin)
            % AVGNODETPUTCHAIN Kotlin-style alias for getAvgNodeTputChain
            node_tput_chain = self.getAvgNodeTputChain(varargin{:});
        end
        
        % Kotlin-style aliases for getTran* methods
        function varargout = tranAvg(self, varargin)
            % TRANAVG Kotlin-style alias for getTranAvg
            [varargout{1:nargout}] = self.getTranAvg(varargin{:});
        end
        
        function rd = tranCdfRespT(self, varargin)
            % TRANCDFRESPT Kotlin-style alias for getTranCdfRespT
            rd = self.getTranCdfRespT(varargin{:});
        end
        
        function rd = tranCdfPassT(self, varargin)
            % TRANCDFPASST Kotlin-style alias for getTranCdfPassT
            rd = self.getTranCdfPassT(varargin{:});
        end
        
        % Kotlin-style aliases for getCdf* methods
        function rd = cdfRespT(self, varargin)
            % CDFRESPT Kotlin-style alias for getCdfRespT
            rd = self.getCdfRespT(varargin{:});
        end
        
        function rd = cdfPassT(self, varargin)
            % CDFPASST Kotlin-style alias for getCdfPassT
            rd = self.getCdfPassT(varargin{:});
        end
        
        % Kotlin-style aliases for getProb* methods
        function pstate = prob(self, varargin)
            % PROB Kotlin-style alias for getProb
            pstate = self.getProb(varargin{:});
        end
        
        function psysstate = probSys(self)
            % PROBSYS Kotlin-style alias for getProbSys
            psysstate = self.getProbSys();
        end
        
        function tf = supportsExactSensitivity(self) %#ok<MANU>
            % TF = SUPPORTSEXACTSENSITIVITY()
            % True when the solver evaluates a product-form recursion that
            % getSensitivityTable can differentiate analytically. False here,
            % so that a solver reaching this base implementation obtains its
            % sensitivities by finite differences on its own predictions.
            % Overridden by SolverMVA and SolverNC.
            tf = false;
        end

        function pnir = probAggr(self, varargin)
            % PROBAGGR Kotlin-style alias for getProbAggr
            pnir = self.getProbAggr(varargin{:});
        end
        
        function pnjoint = probSysAggr(self)
            % PROBSYSAGGR Kotlin-style alias for getProbSysAggr
            pnjoint = self.getProbSysAggr();
        end
        
        function pmarg = probMarg(self, varargin)
            % PROBMARG Kotlin-style alias for getProbMarg
            pmarg = self.getProbMarg(varargin{:});
        end
        
        function lnormconst = probNormConstAggr(self)
            % PROBNORMCONSTAGGR Kotlin-style alias for getProbNormConstAggr
            lnormconst = self.getProbNormConstAggr();
        end

        % Table -> T shorthand aliases for the auxiliary result tables
        % (moment/sensitivity/cache/item/orbit), mirroring aT for getAvgTable.
        function varargout = momentT(self, varargin)
            % MOMENTT Alias for getMomentTable
            [varargout{1:nargout}] = self.getMomentTable(varargin{:});
        end
        function varargout = mT(self, varargin)
            % MT Alias for getMomentTable
            [varargout{1:nargout}] = self.getMomentTable(varargin{:});
        end
        function varargout = getMomentT(self, varargin)
            % GETMOMENTT Alias for getMomentTable
            [varargout{1:nargout}] = self.getMomentTable(varargin{:});
        end
        function varargout = momentChainT(self, varargin)
            % MOMENTCHAINT Alias for getMomentChainTable
            [varargout{1:nargout}] = self.getMomentChainTable(varargin{:});
        end
        function varargout = mCT(self, varargin)
            % MCT Alias for getMomentChainTable
            [varargout{1:nargout}] = self.getMomentChainTable(varargin{:});
        end
        function varargout = getMomentChainT(self, varargin)
            % GETMOMENTCHAINT Alias for getMomentChainTable
            [varargout{1:nargout}] = self.getMomentChainTable(varargin{:});
        end
        function varargout = momentStationT(self, varargin)
            % MOMENTSTATIONT Alias for getMomentStationTable
            [varargout{1:nargout}] = self.getMomentStationTable(varargin{:});
        end
        function varargout = mST(self, varargin)
            % MST Alias for getMomentStationTable
            [varargout{1:nargout}] = self.getMomentStationTable(varargin{:});
        end
        function varargout = getMomentStationT(self, varargin)
            % GETMOMENTSTATIONT Alias for getMomentStationTable
            [varargout{1:nargout}] = self.getMomentStationTable(varargin{:});
        end
        function varargout = sensitivityT(self, varargin)
            % SENSITIVITYT Alias for getSensitivityTable
            [varargout{1:nargout}] = self.getSensitivityTable(varargin{:});
        end
        function varargout = sT(self, varargin)
            % ST Alias for getSensitivityTable
            [varargout{1:nargout}] = self.getSensitivityTable(varargin{:});
        end
        function varargout = getSensitivityT(self, varargin)
            % GETSENSITIVITYT Alias for getSensitivityTable
            [varargout{1:nargout}] = self.getSensitivityTable(varargin{:});
        end
        function varargout = cacheAvgT(self, varargin)
            % CACHEAVGT Alias for getAvgCacheTable
            [varargout{1:nargout}] = self.getAvgCacheTable(varargin{:});
        end
        function varargout = aCaT(self, varargin)
            % ACAT Alias for getAvgCacheTable
            [varargout{1:nargout}] = self.getAvgCacheTable(varargin{:});
        end
        function varargout = getAvgCacheT(self, varargin)
            % GETAVGCACHET Alias for getAvgCacheTable
            [varargout{1:nargout}] = self.getAvgCacheTable(varargin{:});
        end
        function varargout = itemAvgT(self, varargin)
            % ITEMAVGT Alias for getAvgItemTable
            [varargout{1:nargout}] = self.getAvgItemTable(varargin{:});
        end
        function varargout = aIT(self, varargin)
            % AIT Alias for getAvgItemTable
            [varargout{1:nargout}] = self.getAvgItemTable(varargin{:});
        end
        function varargout = getAvgItemT(self, varargin)
            % GETAVGITEMT Alias for getAvgItemTable
            [varargout{1:nargout}] = self.getAvgItemTable(varargin{:});
        end
        function varargout = orbitAvgT(self, varargin)
            % ORBITAVGT Alias for getAvgOrbitTable
            [varargout{1:nargout}] = self.getAvgOrbitTable(varargin{:});
        end
        function varargout = aOT(self, varargin)
            % AOT Alias for getAvgOrbitTable
            [varargout{1:nargout}] = self.getAvgOrbitTable(varargin{:});
        end
        function varargout = getAvgOrbitT(self, varargin)
            % GETAVGORBITT Alias for getAvgOrbitTable
            [varargout{1:nargout}] = self.getAvgOrbitTable(varargin{:});
        end
        function varargout = lossAvgT(self, varargin)
            % LOSSAVGT Alias for getAvgLossTable
            [varargout{1:nargout}] = self.getAvgLossTable(varargin{:});
        end
        function varargout = aLT(self, varargin)
            % ALT Alias for getAvgLossTable
            [varargout{1:nargout}] = self.getAvgLossTable(varargin{:});
        end
        function varargout = getAvgLossT(self, varargin)
            % GETAVGLOSST Alias for getAvgLossTable
            [varargout{1:nargout}] = self.getAvgLossTable(varargin{:});
        end
        function varargout = regionLossAvgT(self, varargin)
            % REGIONLOSSAVGT Alias for getAvgRegionLossTable
            [varargout{1:nargout}] = self.getAvgRegionLossTable(varargin{:});
        end
        function varargout = aRLT(self, varargin)
            % ARLT Alias for getAvgRegionLossTable
            [varargout{1:nargout}] = self.getAvgRegionLossTable(varargin{:});
        end
        function varargout = getAvgRegionLossT(self, varargin)
            % GETAVGREGIONLOSST Alias for getAvgRegionLossTable
            [varargout{1:nargout}] = self.getAvgRegionLossTable(varargin{:});
        end

        % Kotlin-style aliases for get*Handles methods
        function varargout = avgHandles(self)
            % AVGHANDLES Kotlin-style alias for getAvgHandles
            [varargout{1:nargout}] = self.getAvgHandles();
        end
        
        function varargout = tranHandles(self)
            % TRANHANDLES Kotlin-style alias for getTranHandles
            [varargout{1:nargout}] = self.getTranHandles();
        end
        
        function q = avgQLenHandles(self)
            % AVGQLENHANDLES Kotlin-style alias for getAvgQLenHandles
            q = self.getAvgQLenHandles();
        end
        
        function u = avgUtilHandles(self)
            % AVGUTILHANDLES Kotlin-style alias for getAvgUtilHandles
            u = self.getAvgUtilHandles();
        end
        
        function r = avgRespTHandles(self)
            % AVGRESPTHANDLES Kotlin-style alias for getAvgRespTHandles
            r = self.getAvgRespTHandles();
        end
        
        function t = avgTputHandles(self)
            % AVGTPUTHANDLES Kotlin-style alias for getAvgTputHandles
            t = self.getAvgTputHandles();
        end
        
        function a = avgArvRHandles(self)
            % AVGARVRHANDLES Kotlin-style alias for getAvgArvRHandles
            a = self.getAvgArvRHandles();
        end
        
        function w = avgResidTHandles(self)
            % AVGRESIDTHANDLES Kotlin-style alias for getAvgResidTHandles
            w = self.getAvgResidTHandles();
        end
        
        % Kotlin-style aliases for basic get* methods
        function qn = avgQLen(self)
            % AVGQLEN Kotlin-style alias for getAvgQLen
            qn = self.getAvgQLen();
        end
        
        function un = avgUtil(self)
            % AVGUTIL Kotlin-style alias for getAvgUtil
            un = self.getAvgUtil();
        end
        
        function rn = avgRespT(self)
            % AVGRESPT Kotlin-style alias for getAvgRespT
            rn = self.getAvgRespT();
        end
        
        function wn = avgResidT(self)
            % AVGRESIDT Kotlin-style alias for getAvgResidT
            wn = self.getAvgResidT();
        end
        
        function wt = avgWaitT(self)
            % AVGWAITT Kotlin-style alias for getAvgWaitT
            wt = self.getAvgWaitT();
        end
        
        function tn = avgTput(self)
            % AVGTPUT Kotlin-style alias for getAvgTput
            tn = self.getAvgTput();
        end
        
        function an = avgArvR(self)
            % AVGARVR Kotlin-style alias for getAvgArvR
            an = self.getAvgArvR();
        end

    end

    methods (Static)
        % Integer placement decided by an auxiliary solver steady state
        placement = warmStartPlacement(initSolver, sn)

        function [bool, reason] = checkBindingCapacity(model, solverName)
            % [BOOL, REASON] = CHECKBINDINGCAPACITY(MODEL, SOLVERNAME)
            % Shared structural gate for finite station capacity
            % (setCapacity) and finite per-class buffers (classCap), used by
            % the product-form solvers (MVA, NC). A product-form solver has no
            % representation of a finite buffer, so without this gate it
            % silently returns the UNCONSTRAINED answer (e.g. QLen=4 instead
            % of the M/M/1/2 value 0.8525). There is no registry feature name
            % for plain capacity, hence the structural test; this mirrors the
            % native Python check in solvers/solver_mva/solver_mva.py.
            %
            % The test reads the node-level cap/classCap set by the user, NOT
            % sn.cap/sn.classcap: refreshCapacity derives a FINITE sn.classcap
            % (= the chain population) for every closed model, so an sn-level
            % test would reject every closed model.
            %
            % Only a capacity that can actually BIND is rejected. A closed
            % model whose station capacity is at least the total population
            % can never block a job, so the declaration is a no-op and the
            % product-form answer stays exact (a common idiom: setCap(N) on an
            % order-independent station of an N-job closed model). njobs is Inf
            % for an open class, so any finite capacity reachable by an open
            % class binds.
            %
            % Cache models are exempt: Cache.m sets classCap=1 on the
            % retrieval queues it builds, and MVA/NC solve those through their
            % dedicated cache/retrieval analyzers rather than as a buffer
            % constraint.
            bool = true;
            reason = '';
            if ~isa(model, 'Network')
                return
            end
            nodes = model.getNodes();
            for i = 1:numel(nodes)
                if isa(nodes{i}, 'Cache')
                    return
                end
            end
            njobs = model.getStruct().njobs(:)';
            totalJobs = sum(njobs); % Inf as soon as one class is open
            for i = 1:numel(nodes)
                node = nodes{i};
                if ~isa(node, 'Station') || isa(node, 'Source') || isa(node, 'Sink')
                    continue
                end
                if ~isempty(node.cap) && ~isinf(node.cap) && node.cap >= 0 && node.cap < totalJobs
                    bool = false;
                    reason = sprintf('Finite station capacity (setCapacity=%g) at station ''%s'' is not supported by %s. Use SolverCTMC, SolverJMT or SolverLDES.', node.cap, node.getName(), solverName);
                    return
                end
                ccap = node.classCap;
                for r = 1:min(numel(ccap), numel(njobs))
                    if ~isinf(ccap(r)) && ccap(r) > 0 && ccap(r) < njobs(r)
                        bool = false;
                        reason = sprintf('Finite per-class capacity (classCap=%g for class %d) at station ''%s'' is not supported by %s. Use SolverCTMC, SolverJMT or SolverLDES.', ccap(r), r, node.getName(), solverName);
                        return
                    end
                end
            end
        end

        function solvers = getAllSolvers(model, options)
            % SOLVERS = GETALLSOLVERS(MODEL, OPTIONS)

            % Return a cell array with all Network solvers
            if nargin<2 %~exist('options','var')
                options = Solver.defaultOptions;
            end
            solvers = {};
            solvers{end+1} = SolverCTMC(model, options);
            solvers{end+1} = SolverFluid(model, options);
            solvers{end+1} = SolverJMT(model, options);
            solvers{end+1} = SolverMAM(model, options);
            solvers{end+1} = SolverMVA(model, options);
            solvers{end+1} = SolverNC(model, options);
            solvers{end+1} = SolverSSA(model, options);
        end
    end

end
