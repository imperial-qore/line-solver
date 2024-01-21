classdef SolverFluid < NetworkSolver
    % A solver based on fluid and mean-field approximation methods.
    %
    % Copyright (c) 2012-2024, Imperial College London
    % All rights reserved.
    
    methods
        function self = SolverFluid(model,varargin)
            % SELF = SOLVERFLUID(MODEL,VARARGIN)
            
            self@NetworkSolver(model, mfilename);
            self.setOptions(Solver.parseOptions(varargin, self.defaultOptions));
        end
    end
    
    methods
        RD = getTranCdfPassT(self, R);
        [Pnir,logPnir] = getProbAggr(self, ist);
        RD = getCdfRespT(self, R);
        
        function sn = getStruct(self)
            % QN = GETSTRUCT()
            
            % Get data structure summarizing the model
            sn = self.model.getStruct(true);
        end
        
        [QNt,UNt,TNt] = getTranAvg(self,Qt,Ut,Tt);
        
        % solve method is supplied by Solver superclass
        runtime = runAnalyzer(self, options);
    end
    
    methods (Static)
        function [allMethods] = listValidMethods()
            % allMethods = LISTVALIDMETHODS()
            % List valid methods for this solver
            allMethods = {'default','softmin','statedep','closing','matrix'};
        end

        function featSupported = getFeatureSet()
            % FEATSUPPORTED = GETFEATURESET()
            
            featSupported = SolverFeatureSet;
            featSupported.setTrue({
                'ClassSwitch','DelayStation','Queue',...
                'Cox2','Erlang','Exp','HyperExp',...
                'APH', 'Det',...
                'StatelessClassSwitcher','InfiniteServer','SharedServer','Buffer','Dispatcher',...
                'Server','ServiceTunnel',...
                'SchedStrategy_INF','SchedStrategy_PS',...
                'SchedStrategy_DPS','SchedStrategy_FCFS',...
                'SchedStrategy_SIRO',...
                'RoutingStrategy_PROB','RoutingStrategy_RAND',...
                'ClosedClass','Replayer', ...
                'RandomSource','Sink','Source','OpenClass','JobSink'});
            %SolverFluid has very weak performance on open models
        end
        
        function [bool, featSupported] = supports(model)
            % [BOOL, FEATSUPPORTED] = SUPPORTS(MODEL)
            
            featUsed = model.getUsedLangFeatures();
            featSupported = SolverFluid.getFeatureSet();
            bool = SolverFeatureSet.supports(featSupported, featUsed);
        end        
        
        function options = defaultOptions()
            % OPTIONS = DEFAULTOPTIONS()
            options = lineDefaults('Fluid');
        end
    end
end
