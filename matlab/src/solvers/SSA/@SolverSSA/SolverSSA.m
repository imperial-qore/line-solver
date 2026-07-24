classdef SolverSSA < NetworkSolver
    % Stochastic Simulation Analysis solver
    %
    % Implements discrete-event stochastic simulation for queueing network analysis.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    methods
        function self = SolverSSA(model,varargin)
            % SOLVERSSA Create an SSA solver instance
            %
            % @brief Creates a Stochastic Simulation Analysis solver
            % @param model Network model to be analyzed via simulation
            % @param varargin Optional parameters (method, samples, seed, etc.)
            % @return self SolverSSA instance configured for simulation

            % An auxiliary solver passed as first optional argument requests a
            % warm start: its steady-state solution decides the initial state
            % of the simulated trajectory (see NetworkSolver.initFromSolver).
            initSolver = [];
            if ~isempty(varargin) && isa(varargin{1}, 'NetworkSolver')
                initSolver = varargin{1};
                varargin(1) = [];
            end
            self@NetworkSolver(model, mfilename);
            self.setOptions(Solver.parseOptions(varargin, self.defaultOptions));
            self.setLang();
            if ~isempty(initSolver)
                self.initFromSolver(initSolver);
            end
        end
        
        function sn = getStruct(self)
            % QN = GETSTRUCT()
            
            % Get data structure summarizing the model
            sn = self.model.getStruct(true);
        end
        
        [runtime, tranSysState, tranSync] = run(self, options);        
        Prob = getProb(self, node, state);
        ProbAggr = getProbAggr(self, node, state);
        ProbSys = getProbSys(self);
        ProbSysAggr = getProbSysAggr(self);
        tranNodeState = sample(self, node, numEvents, markActivePassive);
        tranNodeStateAggr = sampleAggr(self, node, numEvents, markActivePassive);
        tranSysStateAggr = sampleSysAggr(self, numEvents, markActivePassive);
        tranSysState = sampleSys(self, numEvents, markActivePassive);

        function [allMethods] = listValidMethods(self)
            % allMethods = LISTVALIDMETHODS()
            % List valid methods for this solver
            
            %sn = self.model.getStruct();
            allMethods = {'default','ssa','serial', 'para','parallel','nrm'};
        end

        function bool = isStochasticMethod(self, method) %#ok<INUSD>
            % BOOL = ISSTOCHASTICMETHOD(METHOD)
            % All SSA methods are stochastic simulation.
            bool = true;
        end
    end
    
    methods (Static)

        function featSupported = getFeatureSet()
            % FEATSUPPORTED = GETFEATURESET()

            featSupported = SolverFeatureSet;
            featSupported.setTrue({'Sink','Source','Router',...
                'ClassSwitch','Delay','DelayStation','Queue',...
                'Cache','CacheClassSwitcher',...
                'MAP','MMPP2','MMAP', 'APH', 'PH', 'Replayer',...
                'Coxian','Erlang','Exp','HyperExp',...
                'Det','Gamma','Lognormal','Pareto','Uniform','Weibull',...
                'StatelessClassSwitcher','InfiniteServer',...
                'SharedServer','Buffer','Dispatcher',...
                ... % Finite capacity regions: the NRM sample path honours them.
                'Region', ...
                'Server','JobSink','RandomSource','ServiceTunnel',...
                'SchedStrategy_INF','SchedStrategy_PS',...
                'SchedStrategy_DPS','SchedStrategy_FCFS',...
                'SchedStrategy_GPS','SchedStrategy_LPS','SchedStrategy_SIRO',...
                'SchedStrategy_HOL','SchedStrategy_LCFS',...
                'SchedStrategy_SEPT','SchedStrategy_LEPT',...
                'SchedStrategy_LCFSPR',...
                'SchedStrategy_PSPRIO','SchedStrategy_DPSPRIO','SchedStrategy_GPSPRIO',...
                'SchedStrategy_LCFSPRPRIO','SchedStrategy_FCFSPRPRIO',...
                'SchedStrategy_PAS',...
                'SchedStrategy_OI',...
                'SchedStrategy_POLLING',...
                'RoutingStrategy_RROBIN',...
                'RoutingStrategy_WRROBIN',...
                'RoutingStrategy_JSQ',...
                'RoutingStrategy_SQ',...
                'RoutingStrategy_RL',...
                'RoutingStrategy_PROB','RoutingStrategy_RAND',...
                'ReplacementStrategy_RR', 'ReplacementStrategy_FIFO','ReplacementStrategy_SFIFO','ReplacementStrategy_LRU',...
                'ReplacementStrategy_HLRU','ReplacementStrategy_CLIMB','ReplacementStrategy_QLRU',...
                'SchedStrategy_EXT','ClosedClass','SelfLoopingClass','OpenClass',...
                'OpenSignal','ClosedSignal',...
                'SignalType_NEGATIVE','SignalType_CATASTROPHE',...
                'SignalBatchRemoval','SignalRemovalPolicy',...
                'Fork','Join','Forker','Joiner',...
                'Place', 'Transition', 'Linkage', 'Enabling', 'Inhibiting', 'Timing', 'Firing', 'Storage',...
                'Balking','Reneging','Retrial',...
                'LoadDependence','ClassDependence','JointDependence'});
        end
        
        function [bool, featSupported] = supports(model)
            % [BOOL, FEATSUPPORTED] = SUPPORTS(MODEL)
            
            featUsed = model.getUsedLangFeatures();
            featSupported = SolverSSA.getFeatureSet();
            bool = SolverFeatureSet.supports(featSupported, featUsed);
        end        
               
        function options = defaultOptions()
            % OPTIONS = DEFAULTOPTIONS()

            options = SolverOptions('SSA');
        end

        function libs = getLibrariesUsed(sn, options)
            % GETLIBRARIESUSED Get list of external libraries used by SSA solver
            % SSA uses internal simulation, no external libraries needed
            libs = {};
        end

    end
end
