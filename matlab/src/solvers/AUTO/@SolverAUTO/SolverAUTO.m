classdef SolverAUTO
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
    % - JMT simulation for complex models
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
    end

    properties (Hidden)
        % Network solvers
        CANDIDATE_MVA = 1;
        CANDIDATE_NC = 2;
        CANDIDATE_MAM = 3;
        CANDIDATE_FLUID = 4;
        CANDIDATE_JMT = 5;
        CANDIDATE_SSA = 6;
        CANDIDATE_CTMC = 7;
        CANDIDATE_DES = 8;
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
            switch self.options.method
                case 'sim'
                    % Best simulator - JMT for most cases
                    self.options.method = 'default';
                    switch class(model)
                        case 'Network'
                            self.solvers{1,1} = SolverJMT(model,self.options);
                        case 'LayeredNetwork'
                            self.solvers{1,1} = SolverLQNS(model,self.options);
                        case 'Environment'
                            self.solvers{1,1} = SolverENV(model,@(m) SolverJMT(m,'verbose',false),self.options);
                    end
                case 'exact'
                    % Best exact method - NC for closed, CTMC for small state space
                    self.options.method = 'default';
                    switch class(model)
                        case 'Network'
                            sn = model.getStruct(true);
                            if all(isfinite(sn.njobs)) % closed network
                                self.solvers{1,1} = SolverNC(model,self.options);
                            else
                                self.solvers{1,1} = SolverCTMC(model,self.options);
                            end
                        case 'LayeredNetwork'
                            self.solvers{1,1} = SolverLN(model,@(m) SolverNC(m,'verbose',false),self.options);
                        case 'Environment'
                            self.solvers{1,1} = SolverENV(model,@(m) SolverNC(m,'verbose',false),self.options);
                    end
                case 'fast'
                    % Fast approximate method - MVA (fastest analytical)
                    self.options.method = 'default';
                    switch class(model)
                        case 'Network'
                            self.solvers{1,1} = SolverMVA(model,self.options);
                        case 'LayeredNetwork'
                            self.solvers{1,1} = SolverLN(model,@(m) SolverMVA(m,'verbose',false),self.options);
                        case 'Environment'
                            self.solvers{1,1} = SolverENV(model,@(m) SolverMVA(m,'verbose',false),self.options);
                    end
                case 'accurate'
                    % Accurate approximate method - Fluid or MAM
                    self.options.method = 'default';
                    switch class(model)
                        case 'Network'
                            self.solvers{1,1} = SolverFluid(model,self.options);
                        case 'LayeredNetwork'
                            self.solvers{1,1} = SolverLN(model,@(m) SolverFluid(m,'verbose',false),self.options);
                        case 'Environment'
                            self.solvers{1,1} = SolverENV(model,@(m) SolverFluid(m,'verbose',false),self.options);
                    end
                case 'mam'
                    self.options.method = 'default';
                    self.solvers{1,1} = SolverMAM(model,self.options);
                case 'mva'
                    self.options.method = 'default';
                    self.solvers{1,1} = SolverMVA(model,self.options);
                case 'nc'
                    self.options.method = 'default';
                    self.solvers{1,1} = SolverNC(model,self.options);
                case 'fluid'
                    self.options.method = 'default';
                    self.solvers{1,1} = SolverFluid(model,self.options);
                case 'jmt'
                    self.options.method = 'default';
                    self.solvers{1,1} = SolverJMT(model,self.options);
                case 'ssa'
                    self.options.method = 'default';
                    self.solvers{1,1} = SolverSSA(model,self.options);
                case 'ctmc'
                    self.options.method = 'default';
                    self.solvers{1,1} = SolverCTMC(model,self.options);
                case 'des'
                    self.options.method = 'default';
                    self.solvers{1,1} = SolverDES(model,self.options);
                case 'env'
                    self.options.method = 'default';
                    self.solvers{1,1} = SolverENV(model,@(m) SolverMVA(m,'verbose',false),self.options);
                case 'ln'
                    self.options.method = 'default';
                    self.solvers{1,1} = SolverLN(model,self.options);
                case 'lqns'
                    self.options.method = 'default';
                    self.solvers{1,1} = SolverLQNS(model,self.options);
                case {'default','heur'} % 'ai' method not yet available
                    %solvers sorted from fastest to slowest
                    self.solvers = {};
                    switch class(model)
                        case 'Network'
                            self.solvers{1,self.CANDIDATE_MAM} = SolverMAM(model);
                            self.solvers{1,self.CANDIDATE_MVA} = SolverMVA(model);
                            self.solvers{1,self.CANDIDATE_NC} = SolverNC(model);
                            self.solvers{1,self.CANDIDATE_FLUID} = SolverFluid(model);
                            self.solvers{1,self.CANDIDATE_JMT} = SolverJMT(model);
                            self.solvers{1,self.CANDIDATE_SSA} = SolverSSA(model);
                            self.solvers{1,self.CANDIDATE_CTMC} = SolverCTMC(model);
                            self.solvers{1,self.CANDIDATE_DES} = SolverDES(model);
                            boolSolver = false(length(self.solvers),1);
                            for s=1:length(self.solvers)
                                boolSolver(s) = self.solvers{s}.supports(self.model);
                                self.solvers{s}.setOptions(self.options);
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
            % Check if the solver has computed results
            bool = self.delegate('hasResults', 1);
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

        function bool = supports(self, model)
            % BOOL = SUPPORTS(MODEL)
            % Check if any candidate solver supports the given model
            bool = false;
            for s = 1:length(self.solvers)
                if ~isempty(self.solvers{s}) && self.solvers{s}.supports(model)
                    bool = true;
                    return;
                end
            end
        end

        function allMethods = listValidMethods(self)
            % ALLMETHODS = LISTVALIDMETHODS()
            % List valid methods from all candidate solvers
            allMethods = {'default', 'auto', 'heur', 'sim', 'exact', 'fast', 'accurate'}; % 'ai' not yet available
            for s = 1:length(self.solvers)
                if ~isempty(self.solvers{s})
                    try
                        solverMethods = self.solvers{s}.listValidMethods();
                        allMethods = unique([allMethods, solverMethods]);
                    catch
                        % Solver doesn't implement listValidMethods
                    end
                end
            end
        end

    end

    methods
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
        function AvgChainTable = getAvgChainTable(self)
            % [AVGCHAINTABLE] = GETAVGCHAINTABLE(self)
            try
                AvgChainTable = self.delegate('getAvgChainTable', 1);
            catch
                line_error(mfilename, 'Fatal error in getAvgChainTable.')
            end
        end

        function [AvgQLenTable, QT] = getAvgQLenTable(self, Q, keepDisabled)
            % [AVGQLENTABLE, QT] = GETAVGQLENTABLE(self, Q, keepDisabled)
            if nargin < 2
                Q = [];
            end
            if nargin < 3
                keepDisabled = false;
            end
            [AvgQLenTable, QT] = self.delegate('getAvgQLenTable', 2, Q, keepDisabled);
        end

        function [AvgTputTable, TT] = getAvgTputTable(self, T, keepDisabled)
            % [AVGTPUTTABLE, TT] = GETAVGTPUTTABLE(self, T, keepDisabled)
            if nargin < 2
                T = [];
            end
            if nargin < 3
                keepDisabled = false;
            end
            [AvgTputTable, TT] = self.delegate('getAvgTputTable', 2, T, keepDisabled);
        end

        function [AvgRespTTable, RT] = getAvgRespTTable(self, R, keepDisabled)
            % [AVGRESPTTABLE, RT] = GETAVGRESPTTABLE(self, R, keepDisabled)
            if nargin < 2
                R = [];
            end
            if nargin < 3
                keepDisabled = false;
            end
            [AvgRespTTable, RT] = self.delegate('getAvgRespTTable', 2, R, keepDisabled);
        end

        function [AvgUtilTable, UT] = getAvgUtilTable(self, U, keepDisabled)
            % [AVGUTILTABLE, UT] = GETAVGUTILTABLE(self, U, keepDisabled)
            if nargin < 2
                U = [];
            end
            if nargin < 3
                keepDisabled = false;
            end
            [AvgUtilTable, UT] = self.delegate('getAvgUtilTable', 2, U, keepDisabled);
        end

        function AvgSysTable = getAvgSysTable(self)
            % [AVGSYSTABLE] = GETAVGSYSTABLE(self)
            AvgSysTable = self.delegate('getAvgSysTable', 1);
        end

        function AvgNodeTable = getAvgNodeTable(self)
            % [AVGNODETABLE] = GETAVGNODETABLE(self)
            AvgNodeTable = self.delegate('getAvgNodeTable', 1);
        end

        function AvgTable = getAvgTable(self)
            % [AVGTABLE] = GETAVGTABLE(self)
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

        function sampleNodeState = sample(self, node, numSamples)
            sampleNodeState = self.delegate('sample', 1, node, numSamples);
        end

        function stationStateAggr = sampleAggr(self, node, numSamples)
            stationStateAggr = self.delegate('sampleAggr', 1, node, numSamples);
        end

        function tranSysState = sampleSys(self, numSamples)
            tranSysState = self.delegate('sampleSys', 1, numSamples);
        end

        function sysStateAggr = sampleSysAggr(self, numSamples)
            sysStateAggr = self.delegate('sampleSysAggr', 1, numSamples);
        end

        function RD = getCdfRespT(self, R)
            if nargin>1
                RD = self.delegate('getCdfRespT', 1, R);
            else
                RD = self.delegate('getCdfRespT', 1);
            end
        end

        function Pnir = getProb(self, node, state)
            Pnir = self.delegate('getProb', 1, node, state);
        end

        function Pnir = getProbAggr(self, node, state_a)
            Pnir = self.delegate('getProbAggr', 1, node, state_a);
        end

        function Pn = getProbSys(self)
            Pn = self.delegate('getProbSys',1);
        end

        function Pn = getProbSysAggr(self)
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
        function [AvgNodeChainTable, QTc, UTc, RTc, WTc, ATc, TTc] = getAvgNodeChainTable(self, Q, U, R, T)
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

        % Kotlin-style aliases for table methods
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

        % Kotlin-style aliases for composite methods
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

        function varargout = avg(self, varargin)
            % AVG Kotlin-style alias for getAvg
            [varargout{1:nargout}] = self.getAvg(varargin{:});
        end

        function sys_resp_time = avgSysRespT(self, varargin)
            % AVGSYSRESPT Kotlin-style alias for getAvgSysRespT
            sys_resp_time = self.getAvgSysRespT(varargin{:});
        end

        function sys_tput = avgSysTput(self, varargin)
            % AVGSYSTPUT Kotlin-style alias for getAvgSysTput
            sys_tput = self.getAvgSysTput(varargin{:});
        end

        % Kotlin-style aliases for chain methods
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

        % Kotlin-style aliases for node chain methods
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

        % Kotlin-style aliases for transient methods
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

        function varargout = tranProb(self, varargin)
            % TRANPROB Kotlin-style alias for getTranProb
            [varargout{1:nargout}] = self.getTranProb(varargin{:});
        end

        function varargout = tranProbAggr(self, varargin)
            % TRANPROBAGGR Kotlin-style alias for getTranProbAggr
            [varargout{1:nargout}] = self.getTranProbAggr(varargin{:});
        end

        function varargout = tranProbSys(self)
            % TRANPROBSYS Kotlin-style alias for getTranProbSys
            [varargout{1:nargout}] = self.getTranProbSys();
        end

        function varargout = tranProbSysAggr(self)
            % TRANPROBSYSAGGR Kotlin-style alias for getTranProbSysAggr
            [varargout{1:nargout}] = self.getTranProbSysAggr();
        end

        % Kotlin-style aliases for CDF methods
        function rd = cdfRespT(self, varargin)
            % CDFRESPT Kotlin-style alias for getCdfRespT
            rd = self.getCdfRespT(varargin{:});
        end

        function rd = cdfPassT(self, varargin)
            % CDFPASST Kotlin-style alias for getCdfPassT
            rd = self.getCdfPassT(varargin{:});
        end

        function varargout = perctRespT(self, varargin)
            % PERCTRESPT Kotlin-style alias for getPerctRespT
            [varargout{1:nargout}] = self.getPerctRespT(varargin{:});
        end

        % Kotlin-style aliases for probability methods
        function pstate = prob(self, varargin)
            % PROB Kotlin-style alias for getProb
            pstate = self.getProb(varargin{:});
        end

        function psysstate = probSys(self)
            % PROBSYS Kotlin-style alias for getProbSys
            psysstate = self.getProbSys();
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

        % Kotlin-style aliases for handle methods
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

        % Kotlin-style aliases for basic metric methods
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

        % Kotlin-style aliases for table methods with parameters
        function varargout = avgQLenTable(self, varargin)
            % AVGQLENTABLE Kotlin-style alias for getAvgQLenTable
            [varargout{1:nargout}] = self.getAvgQLenTable(varargin{:});
        end

        function varargout = avgUtilTable(self, varargin)
            % AVGUTILTABLE Kotlin-style alias for getAvgUtilTable
            [varargout{1:nargout}] = self.getAvgUtilTable(varargin{:});
        end

        function varargout = avgRespTTable(self, varargin)
            % AVGRESPTTABLE Kotlin-style alias for getAvgRespTTable
            [varargout{1:nargout}] = self.getAvgRespTTable(varargin{:});
        end

        function varargout = avgTputTable(self, varargin)
            % AVGTPUTTABLE Kotlin-style alias for getAvgTputTable
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

        % Kotlin-style aliases for LayeredNetwork/EnsembleSolver methods
        function varargout = ensembleAvg(self)
            % ENSEMBLEAVG Kotlin-style alias for getEnsembleAvg
            [varargout{1:nargout}] = self.getEnsembleAvg();
        end

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

        function avg_tables = getEnsembleAvgTs(self)
            % GETENSEMBLEAVGTS Short alias for getEnsembleAvgTables
            avg_tables = self.getEnsembleAvgTables();
        end

        function avg_tables = ensembleAvgTs(self)
            % ENSEMBLEAVGTS Short alias for ensembleAvgTables
            avg_tables = self.ensembleAvgTables();
        end

    end
end
