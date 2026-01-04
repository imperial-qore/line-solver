classdef SolverFLD < NetworkSolver
    % FLD - Fluid/Mean-Field Approximation solver for large-scale network analysis
    %
    % Implements fluid approximation by replacing discrete populations with continuous levels.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    methods
        function self = SolverFLD(model,varargin)
            % SOLVERFLUID Create a Fluid solver instance
            %
            % @brief Creates a Fluid solver for continuous approximation analysis
            % @param model Network model to be analyzed via fluid approximation
            % @param varargin Optional parameters (method, tolerances, etc.)
            % @return self SolverFLD instance configured for fluid analysis
            
            self@NetworkSolver(model, mfilename);
            self.setOptions(Solver.parseOptions(varargin, self.defaultOptions));
            self.setLang();
        end
    end
    
    methods
        RD = getTranCdfPassT(self, R);
        [Pnir,logPnir] = getProbAggr(self, ist);
        RD = getCdfRespT(self, R);
        [AoI, PAoI, aoiTable] = getAvgAoI(self);
        [AoI_cdf, PAoI_cdf] = getCdfAoI(self, t_values);

        function sn = getStruct(self)
            % QN = GETSTRUCT()
            
            % Get data structure summarizing the model
            sn = self.model.getStruct(true);
        end
        
        [QNt,UNt,TNt] = getTranAvg(self,Qt,Ut,Tt);
        
        % solve method is supplied by Solver superclass
        runtime = runAnalyzer(self, options);

        function [allMethods] = listValidMethods(self)
            % allMethods = LISTVALIDMETHODS()
            % List valid methods for this solver
            sn = self.model.getStruct();
            allMethods = {'default','softmin','pnorm','statedep','closing','matrix','diffusion','mfq'};
        end
    end
    
    methods (Static)

        function featSupported = getFeatureSet()
            % FEATSUPPORTED = GETFEATURESET()
            
            featSupported = SolverFeatureSet;
            featSupported.setTrue({
                'ClassSwitch','Delay','DelayStation','Queue',...
                'Cox2','Erlang','Exp','HyperExp',...
                'APH', 'Det',...
                'StatelessClassSwitcher','InfiniteServer','SharedServer','Buffer','Dispatcher',...
                'Server','ServiceTunnel',...
                'SchedStrategy_INF','SchedStrategy_PS',...
                'SchedStrategy_DPS','SchedStrategy_FCFS',...
                'SchedStrategy_SIRO','SchedStrategy_LCFS','SchedStrategy_LCFSPR',...
                'RoutingStrategy_PROB','RoutingStrategy_RAND',...
                'ClosedClass','SelfLoopingClass','Replayer', ...
                'RandomSource','Sink','Source','OpenClass','JobSink'});
            %SolverFLD has very weak performance on open models
        end
        
        function [bool, featSupported] = supports(model)
            % [BOOL, FEATSUPPORTED] = SUPPORTS(MODEL)
            
            featUsed = model.getUsedLangFeatures();
            featSupported = SolverFLD.getFeatureSet();
            bool = SolverFeatureSet.supports(featSupported, featUsed);
        end        
        
        function options = defaultOptions()
            % OPTIONS = DEFAULTOPTIONS()
            options = SolverOptions('Fluid');
        end

        function libs = getLibrariesUsed(sn, options)
            % GETLIBRARIESUSED Get list of external libraries used by FLD solver
            % FLD uses internal fluid approximation algorithms, no external libraries needed
            libs = {};
        end
    end
end
