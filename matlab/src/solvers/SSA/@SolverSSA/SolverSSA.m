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
            
            self@NetworkSolver(model, mfilename);
            self.setOptions(Solver.parseOptions(varargin, self.defaultOptions));
            self.setLang();
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
        tranNodeState = sample(self, node, numSamples, markActivePassive);
        tranNodeStateAggr = sampleAggr(self, node, numSamples, markActivePassive);
        tranSysStateAggr = sampleSysAggr(self, numSamples, markActivePassive);
        tranSysState = sampleSys(self, numSamples, markActivePassive);

        function [allMethods] = listValidMethods(self)
            % allMethods = LISTVALIDMETHODS()
            % List valid methods for this solver
            
            %sn = self.model.getStruct();
            allMethods = {'default','ssa','serial', 'para','parallel','nrm'};
        end
    end
    
    methods (Static)

        function featSupported = getFeatureSet()
            % FEATSUPPORTED = GETFEATURESET()

            featSupported = SolverFeatureSet;
            featSupported.setTrue({'Sink','Source','Router',...
                'ClassSwitch','Delay','DelayStation','Queue',...
                'Cache','CacheClassSwitcher',...
                'MAP','MMPP2', 'APH', 'PH',...
                'Coxian','Erlang','Exp','HyperExp',...
                'Det','Gamma','Lognormal','Pareto','Uniform','Weibull',...
                'StatelessClassSwitcher','InfiniteServer',...
                'SharedServer','Buffer','Dispatcher',...
                'Server','JobSink','RandomSource','ServiceTunnel',...
                'SchedStrategy_INF','SchedStrategy_PS',...
                'SchedStrategy_DPS','SchedStrategy_FCFS',...
                'SchedStrategy_GPS','SchedStrategy_LPS','SchedStrategy_SIRO',...
                'SchedStrategy_HOL','SchedStrategy_LCFS',...
                'SchedStrategy_SEPT','SchedStrategy_LEPT',...
                'SchedStrategy_LCFSPR',...
                'SchedStrategy_PSPRIO','SchedStrategy_DPSPRIO','SchedStrategy_GPSPRIO',...
                'RoutingStrategy_RROBIN',...
                'RoutingStrategy_PROB','RoutingStrategy_RAND',...
                'ReplacementStrategy_RR', 'ReplacementStrategy_FIFO','ReplacementStrategy_SFIFO','ReplacementStrategy_LRU',...
                'SchedStrategy_EXT','ClosedClass','SelfLoopingClass','OpenClass',...
                'Place', 'Transition', 'Linkage', 'Enabling', 'Timing', 'Firing', 'Storage'});
            %                'Fork','Join','Forker','Joiner',...
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
