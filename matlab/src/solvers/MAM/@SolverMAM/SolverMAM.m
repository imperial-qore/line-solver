classdef SolverMAM < NetworkSolver
    % A solver based on matrix-analytic methods.
    %
    % Copyright (c) 2012-2024, Imperial College London
    % All rights reserved.
    
    methods
        function self = SolverMAM(model,varargin)
            % SELF = SOLVERMAM(MODEL,VARARGIN)
            
            self@NetworkSolver(model, mfilename);
            self.setOptions(Solver.parseOptions(varargin, self.defaultOptions));
        end
                
        function sn = getStruct(self)
            % QN = GETSTRUCT()
            
            % Get data structure summarizing the model
            sn = self.model.getStruct(true);
        end
        
        runtime = runAnalyzer(self, options);
        RD = getCdfRespT(self, R);
    end
    
    methods (Static)

        function [allMethods] = listValidMethods()
            % allMethods = LISTVALIDMETHODS()
            % List valid methods for this solver
            allMethods = {'default','dec.source','dec.mmap','dec.poisson','qna'};
        end

        function featSupported = getFeatureSet()
            % FEATSUPPORTED = GETFEATURESET()
            
            featSupported = SolverFeatureSet;
            featSupported.setTrue({'Sink','Source',...
                'DelayStation','Queue',...
                'APH','Coxian','Erlang','Exp','HyperExp','MMPP2','MAP',...
                'StatelessClassSwitcher','InfiniteServer',...
                'ClassSwitch', ...
                'SharedServer','Buffer','Dispatcher',...
                'Server','JobSink','RandomSource','ServiceTunnel',...
                'SchedStrategy_INF','SchedStrategy_PS','SchedStrategy_HOL',...
                'SchedStrategy_FCFS',...
                'RoutingStrategy_PROB','RoutingStrategy_RAND',...
                'ClosedClass',...
                'OpenClass'});       
                %'Fork','Join','Forker','Joiner',...
        end
        
        function [bool, featSupported] = supports(model)
            % [BOOL, FEATSUPPORTED] = SUPPORTS(MODEL)
            
            featUsed = model.getUsedLangFeatures();
            featSupported = SolverMAM.getFeatureSet();
            bool = SolverFeatureSet.supports(featSupported, featUsed);
        end
        
        function options = defaultOptions()
            % OPTIONS = DEFAULTOPTIONS()
            options = lineDefaults('MAM');
        end
    end
end
