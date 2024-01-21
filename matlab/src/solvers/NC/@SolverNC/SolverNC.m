classdef SolverNC < NetworkSolver
    % A solver based on normalizing constant methods.
    %
    % Copyright (c) 2012-2024, Imperial College London
    % All rights reserved.
    
    methods
        function self = SolverNC(model,varargin)
            % SELF = SOLVERNC(MODEL,VARARGIN)
            
            self@NetworkSolver(model, mfilename);
            self.setOptions(Solver.parseOptions(varargin, self.defaultOptions));
        end
        
        runtime = runAnalyzer(self, options)
        Pnir = getProb(self, node, state)
        Pnir = getProbAggr(self, node, state_a)
        Pn   = getProbSys(self)
        Pn   = getProbSysAggr(self)
        RD = getCdfRespT(self, R);
        
        function [normConst,lNormConst] = getNormalizingConstant(self)
            normConst = exp(getProbNormConstAggr(self));
            lNormConst = getProbNormConstAggr(self);
        end
        
        [lNormConst] = getProbNormConstAggr(self)
        
        function sn = getStruct(self)
            % QN = GETSTRUCT()
            
            % Get data structure summarizing the model
            sn = getStruct(self.model, false); %no need for initial state
        end
        
    end
    
    methods (Static)
        function [allMethods] = listValidMethods()
            % allMethods = LISTVALIDMETHODS()
            % List valid methods for this solver
            allMethods = {'default','exact','imci','ls',...
                'le','mmint2','panacea','pana',...
                'kt','deterministic','sampling',...
                'propfair','comom','mom','cub','brute',...
                'rd', 'nrp','nr.probit','nrl','nr.logit','gm'};
        end
        
        function featSupported = getFeatureSet()
            % FEATSUPPORTED = GETFEATURESET()
            
            featSupported = SolverFeatureSet;
            featSupported.setTrue({'Sink','Source',...
                'ClassSwitch','DelayStation','Queue',...
                'APH','Coxian','Erlang','Det','Exp','HyperExp',...
                'StatelessClassSwitcher','InfiniteServer',...
                'SharedServer','Buffer','Dispatcher',...
                'Server','JobSink','RandomSource','ServiceTunnel',...
                'SchedStrategy_INF','SchedStrategy_PS','SchedStrategy_SIRO',...
                'RoutingStrategy_PROB','RoutingStrategy_RAND',...
                'SchedStrategy_FCFS','ClosedClass','ClosedClass',...
                'Cache','CacheClassSwitcher','OpenClass'});
            %'OpenClass',...
        end
        
        function [bool, featSupported] = supports(model)
            % [BOOL, FEATSUPPORTED] = SUPPORTS(MODEL)
            
            featUsed = model.getUsedLangFeatures();
            featSupported = SolverNC.getFeatureSet();
            bool = SolverFeatureSet.supports(featSupported, featUsed);
        end        
        
        function options = defaultOptions()
            % OPTIONS = DEFAULTOPTIONS()
            options = lineDefaults('NC');
        end
    end
end
