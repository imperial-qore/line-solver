classdef SolverMAM < NetworkSolver
    % Matrix-Analytic and RCAT methods solver
    %
    % Implements matrix-analytic methods and RCAT for structured Markov chain analysis.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    methods
        function self = SolverMAM(model,varargin)
            % SOLVERMAM Create a Matrix-Analytic Methods solver instance
            %
            % @brief Creates a MAM solver for structured Markov chain analysis
            % @param model Network model to be analyzed via matrix-analytic methods
            % @param varargin Optional parameters (method, tolerance, etc.)
            % @return self SolverMAM instance configured for MAM analysis
            
            self@NetworkSolver(model, mfilename);
            self.setOptions(Solver.parseOptions(varargin, self.defaultOptions));
            self.setLang();
        end
                
        function sn = getStruct(self)
            % QN = GETSTRUCT()
            
            % Get data structure summarizing the model
            sn = self.model.getStruct(true);
        end
        
        runtime = runAnalyzer(self, options);
        RD = getCdfRespT(self, R);

            function [allMethods] = listValidMethods(self)
            % allMethods = LISTVALIDMETHODS()
            % List valid methods for this solver
            sn = self.model.getStruct();
            % Note: Method order must match expected values in test files.
            % New methods should be added at the end to preserve index alignment.
            % 'exact' method removed - autocat moved to line-legacy.git
            allMethods = {'default','dec.source','dec.mmap','dec.poisson','mna','inap','ldqbd','inapplus'};
        end
end
    
    methods (Static)


        function featSupported = getFeatureSet()
            % FEATSUPPORTED = GETFEATURESET()

            featSupported = SolverFeatureSet;
            % MAM features
            featSupported.setTrue({'Sink','Source',...
                'Fork','Join','Forker','Joiner',... % Fork-Join support (via FJ_codes)
                'Delay','DelayStation','Queue',...
                'APH','Coxian','Erlang','Exp','HyperExp','MMPP2','MAP',...
                'Det','Gamma','Lognormal','Pareto','Uniform','Weibull',...
                'StatelessClassSwitcher','InfiniteServer',...
                'ClassSwitch', ...
                'SharedServer','Buffer','Dispatcher',...
                'Server','JobSink','RandomSource','ServiceTunnel',...
                'SchedStrategy_INF','SchedStrategy_PS','SchedStrategy_HOL',...
                'SchedStrategy_FCFS',...
                'RoutingStrategy_PROB','RoutingStrategy_RAND',...
                'ClosedClass','SelfLoopingClass',...
                'OpenClass'});
            % Add RCAT (AG) features
            featSupported.setTrue({'Sink', 'Source', ...
                'Fork','Join','Forker','Joiner',... % Fork-Join support (via FJ_codes)
                'Delay', 'DelayStation', 'Queue', ...
                'APH', 'Coxian', 'Erlang', 'Exp', 'HyperExp', ...
                'Det','Gamma','Lognormal','Pareto','Uniform','Weibull',...
                'StatelessClassSwitcher', 'InfiniteServer', ...
                'SharedServer', 'Buffer', 'Dispatcher', ...
                'Server', 'JobSink', 'RandomSource', 'ServiceTunnel', ...
                'SchedStrategy_INF', 'SchedStrategy_PS', ...
                'SchedStrategy_FCFS', ...
                'RoutingStrategy_PROB', 'RoutingStrategy_RAND', ...
                'ClosedClass', ...
                'OpenClass'});
            % Add BMAP/PH/N/N retrial queue features
            featSupported.setTrue({'Retrial', 'BMAP', 'PH'});
        end
        
        function [bool, featSupported] = supports(model)
            % [BOOL, FEATSUPPORTED] = SUPPORTS(MODEL)
            
            featUsed = model.getUsedLangFeatures();
            featSupported = SolverMAM.getFeatureSet();
            bool = SolverFeatureSet.supports(featSupported, featUsed);
        end
        
        function options = defaultOptions()
            % OPTIONS = DEFAULTOPTIONS()
            options = SolverOptions('MAM');
        end

        function libs = getLibrariesUsed(sn, options)
            % GETLIBRARIESUSED Get list of external libraries used by MAM solver
            % Detect libraries used by MAM solver based on topology and method
            libs = {};

            % MAMSolver used for matrix-analytic methods (M/G/1, GI/M/1 types)
            % This includes default and decomposition methods
            if ismember(options.method, {'default', 'dec.source', 'dec.mmap', 'dec.poisson'})
                libs{end+1} = 'MAMSolver';
            end

            % Q-MAM used for specific RCAT-based methods
            if ismember(options.method, {'mna', 'inap', 'inapplus'})
                libs{end+1} = 'Q-MAM';
            end

            % SMCSolver used for QBD (Quasi-Birth-Death) analysis
            % Currently available but not actively used in default paths
            % Uncomment when QBD methods are activated:
            % if ismember(options.method, {'qbd'})
            %     libs{end+1} = 'SMCSolver';
            % end

            % BUTools usage detection (via KPCToolbox MAP functions)
            % BUTools is used when analyzing MAP/PH distributions
            if ~isempty(sn) && isfield(sn, 'proc') && ~isempty(sn.proc) && any(~cellfun(@isempty, sn.proc(:)))
                libs{end+1} = 'BUTools';
            end

            % Remove duplicates and maintain order
            libs = unique(libs, 'stable');
        end
    end
end
