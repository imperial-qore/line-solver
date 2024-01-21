classdef SolverSSA < NetworkSolver
    % A solver based on discrete-event stochastic simulation analysis.
    %
    % Copyright (c) 2012-2024, Imperial College London
    % All rights reserved.
    methods
        function self = SolverSSA(model,varargin)
            % SELF = SOLVERSSA(MODEL,VARARGIN)
            
            self@NetworkSolver(model, mfilename);
            self.setOptions(Solver.parseOptions(varargin, self.defaultOptions));
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
    end
    
    methods (Static)
        function [allMethods] = listValidMethods()
            % allMethods = LISTVALIDMETHODS()
            % List valid methods for this solver
            allMethods = {'default','ssa','serial.hash','serial',...
            'para','parallel','para.hash','parallel.hash','hashed',...
            'taussa','tauleap'};

        end

        function featSupported = getFeatureSet()
            % FEATSUPPORTED = GETFEATURESET()
            
            featSupported = SolverFeatureSet;                
            featSupported.setTrue({'Sink','Source','Router',...
                'ClassSwitch','DelayStation','Queue',...
                'Cache','CacheClassSwitcher',...
                'MAP','MMPP2', 'APH', 'PH',...
                'Coxian','Erlang','Exp','HyperExp',...
                'StatelessClassSwitcher','InfiniteServer',...
                'SharedServer','Buffer','Dispatcher',...
                'Server','JobSink','RandomSource','ServiceTunnel',...
                'SchedStrategy_INF','SchedStrategy_PS',...
                'SchedStrategy_DPS','SchedStrategy_FCFS',...
                'SchedStrategy_GPS','SchedStrategy_SIRO',...
                'SchedStrategy_HOL','SchedStrategy_LCFS',...
                'SchedStrategy_SEPT','SchedStrategy_LEPT',...
                'SchedStrategy_LCFSPR',...
                'RoutingStrategy_RROBIN',...
                'RoutingStrategy_PROB','RoutingStrategy_RAND',...
                'SchedStrategy_EXT','ClosedClass','OpenClass'});
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
            
            options = lineDefaults('SSA');
        end
        
    end
end
