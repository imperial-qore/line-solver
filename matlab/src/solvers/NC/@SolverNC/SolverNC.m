classdef SolverNC < NetworkSolver
    % SolverNC Normalizing Constant solver for product-form networks
    %
    % SolverNC implements normalizing constant algorithms for analyzing closed
    % product-form queueing networks. It computes the normalizing constant and
    % associated performance measures efficiently without explicitly enumerating
    % all network states, making it suitable for medium to large closed networks.
    %
    % @brief Normalizing constant solver for efficient closed network analysis
    %
    % Key characteristics:
    % - Normalizing constant computation for product-form networks
    % - Avoids explicit state enumeration
    % - Efficient algorithms for closed networks
    % - Multiple computational methods (exact, approximation)
    % - State probability computation via normalization
    %
    % NC solver methods:
    % - Exact normalizing constant computation
    % - IMCI (Improved Modular Computer Implementation)
    % - Linearizer methods (LS, LE)
    % - Interpolation methods (MMINT2, GLEINT)
    % - Approximation methods (CA, Panacea)
    %
    % SolverNC is ideal for:
    % - Closed product-form networks
    % - Medium to large population networks
    % - Systems requiring efficient exact solutions
    % - Networks with complex routing patterns
    % - Performance analysis requiring state probabilities
    %
    % Example:
    % @code
    % solver = SolverNC(model, 'method', 'exact');
    % solver.getProbAggr();          % State probabilities
    % solver.getNormalizingConstant(); % Normalizing constant
    % @endcode
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    methods
        function self = SolverNC(model,varargin)
            % SOLVERNC Create a Normalizing Constant solver instance
            %
            % @brief Creates an NC solver for product-form network analysis
            % @param model Network model to be analyzed via normalizing constant methods
            % @param varargin Optional parameters (method, tolerance, etc.)
            % @return self SolverNC instance configured for NC analysis
            
            self@NetworkSolver(model, mfilename);
            self.setOptions(Solver.parseOptions(varargin, self.defaultOptions));
            self.setLang();
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
            sn = self.model.getStruct(false); %no need for initial state
        end
        
        function [allMethods] = listValidMethods(self)
            % allMethods = LISTVALIDMETHODS()
            % List valid methods for this solver
            sn = self.model.getStruct();
            allMethods = {'default','exact','imci','ls',...
                'le','mmint2','gleint','panacea','ca',...
                'kt','sampling',...
                'propfair','comom','cub',...
                'rd', 'nrp','nrl','gm','mem'};
        end
    end
    
    methods (Static)
        
        function featSupported = getFeatureSet()
            % FEATSUPPORTED = GETFEATURESET()
            
            featSupported = SolverFeatureSet;
            featSupported.setTrue({'Sink','Source',...
                'ClassSwitch','Delay','DelayStation','Queue',...
                'APH','Coxian','Erlang','Det','Exp','HyperExp',...
                'StatelessClassSwitcher','InfiniteServer',...
                'SharedServer','Buffer','Dispatcher',...
                'Server','JobSink','RandomSource','ServiceTunnel',...
                'SchedStrategy_INF','SchedStrategy_PS','SchedStrategy_SIRO',...
                'SchedStrategy_LCFS','SchedStrategy_LCFSPR',...
                'RoutingStrategy_PROB','RoutingStrategy_RAND',...
                'SchedStrategy_FCFS','ClosedClass','SelfLoopingClass',...
                'Cache','CacheClassSwitcher','OpenClass',            ...
                'ReplacementStrategy_RR', 'ReplacementStrategy_FIFO'});
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
            options = SolverOptions('NC');
        end

        function libs = getLibrariesUsed(sn, options)
            % GETLIBRARIESUSED Get list of external libraries used by NC solver
            % NC uses internal normalizing constant algorithms, no external libraries needed
            libs = {};
        end
    end
end
