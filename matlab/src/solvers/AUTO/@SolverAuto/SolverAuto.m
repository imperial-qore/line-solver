classdef SolverAuto
    %A solver that selects the solution method automatically based on the
    %model characteristics.

    %Copyright (c) 2012-2023, Imperial College London
    %All rights reserved.

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
        % LayeredNetwork solvers
        CANDIDATE_LQNS = 1;
        CANDIDATE_LN_NC = 2;
        CANDIDATE_LN_MVA = 3;
        CANDIDATE_LN_MAM = 4;
        CANDIDATE_LN_FLUID = 5;
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
        %Constructor
        function self = SolverAuto(model, varargin)
            %SELF = LINE(MODEL, VARARGIN)
            self.options = Solver.parseOptions(varargin, Solver.defaultOptions);
            self.model = model;
            self.name = 'SolverAuto';
            if self.options.verbose
                %line_printf('Running LINE version %s',model.getVersion);
            end
            switch self.options.method
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
                case 'ln'
                    self.options.method = 'default';
                    self.solvers{1,1} = SolverLN(model,self.options);
                case 'lqns'
                    self.options.method = 'default';
                    self.solvers{1,1} = SolverLQNS(model,self.options);
                case {'default','ai','nn','heur'}
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
                            %boolSolver = [];
                            %for s=1:length(self.solvers)
                            %    boolSolver(s) = self.solvers{s}.supports(self.model);
                            %    self.solvers{s}.setOptions(self.options);
                            %end
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

    end

    methods
        % delegate execution of method to chosen solver
        varargout = delegate(self, method, nretout, varargin);

        % chooseSolver: choses a solver from static properties of the model
        solver = chooseSolver(self, method);
        % AI-based choice solver
        solver = chooseSolverAI(self, method);
        % Heuristic choice of solver
        solver = chooseSolverHeur(self, method);
        % AI-based choice solver for Avg* methods
        solver = chooseAvgSolverAI(self);
        solver = chooseAvgSolverNN(self); % AI v2
        % Heuristic choice of solver for Avg* methods
        solver = chooseAvgSolverHeur(self);

    end

    methods
        function reset(self)
            for s=1:length(self.solvers)
                self.solvers{s}.reset();
            end
        end
        function setDoChecks(self,bool)
            for s=1:length(self.solvers)
                self.solvers{s}.setDoChecks(bool);
            end
        end
    end

    methods
        function AvgChainTable = getAvgChainTable(self)
            % [AVGCHAINTABLE] = GETAVGCHAINTABLE(self)
            try
            AvgChainTable = self.delegate('getAvgChainTable', 1);
            catch
                keyboard
            end
        end

        function AvgQlenTable = getAvgQLenTable(self)
            AvgQlenTable = self.delegate('getAvgQLenTable', 1);
        end

        function AvgTputTable = getAvgTputTable(self)
            AvgTputTable = self.delegate('getAvgTputTable', 1);
        end

        function AvgRespTTable = getAvgRespTTable(self)
            AvgRespTTable = self.delegate('getAvgRespTTable', 1);
        end

        function AvgUtilTable = getAvgUtilTable(self)
            AvgUtilTable = self.delegate('getAvgUtilTable', 1);
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

    end
end
