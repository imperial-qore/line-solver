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

            % An auxiliary solver passed as first optional argument requests a
            % warm start: its steady-state solution decides the initial state
            % of the ODE integration (see NetworkSolver.initFromSolver).
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
    end
    
    methods
        RD = getTranCdfPassT(self, R);
        [Pnir,logPnir] = getProbAggr(self, ist);
        RD = getCdfRespT(self, R);
        RD = getCdfPassT(self, R);
        RD = getCdfPT(self, R);
        [AoI, PAoI, aoiTable] = getAvgAoI(self);
        [AoI_cdf, PAoI_cdf] = getCdfAoI(self, t_values);

        function sn = getStruct(self)
            % QN = GETSTRUCT()
            
            % Get data structure summarizing the model
            sn = self.model.getStruct(true);
        end
        
        [QNt,UNt,TNt] = getTranAvg(self,Qt,Ut,Tt);

        % export the mean-field ODE system to LaTeX (scalar or matrix notation)
        [tex, sys] = exportODEs(self, filename, notation);

        % solve method is supplied by Solver superclass
        runtime = runAnalyzer(self, options);

        function [allMethods] = listValidMethods(self)
            % allMethods = LISTVALIDMETHODS()
            % List valid methods for this solver
            allMethods = {'default','softmin','pnorm','statedep','closing','matrix','diffusion','mfq','rmf','tbi'};
        end
        function featSupported = getMethodFeatureSet(self, method) %#ok<INUSD>
            % All FLD methods share the solver-level feature envelope.
            %
            % Defining this is what lets NetworkSolver.supportsModelMethod name
            % the offending features: with no method feature set it falls back
            % to the coarse supports(model), which returns an empty reason, so
            % the gate could only report "features not supported" without
            % saying which ones.
            %
            % A non-Network model (e.g. a LayeredNetwork) has no
            % getUsedLangFeatures, so it keeps the coarse path and any
            % structural checks or redirects that operate on such models.
            if ~isa(self.model, 'Network')
                featSupported = [];
                return;
            end
            featSupported = SolverFLD.getFeatureSet();
        end


    end
    
    methods (Static)

        function featSupported = getFeatureSet()
            % FEATSUPPORTED = GETFEATURESET()
            
            featSupported = SolverFeatureSet;
            featSupported.setTrue({
                'ClassSwitch','Delay','DelayStation','Queue',...
                'Cache','CacheClassSwitcher',...
                ... % 'CacheRetrieval' is deliberately NOT declared: no fluid
                ... % code anywhere implements delayed-hit retrieval, and on
                ... % examples/basic/cacheModel/retrieval_simple the ODE
                ... % returned zero QLen, Util and Tput on every row while jobs
                ... % arrived at rate 1, i.e. flow was not conserved. Refusing
                ... % the model is the honest answer; the JAR SolverFluid does
                ... % the same.
                'Cox2','Coxian','Erlang','Exp','HyperExp',...
                'APH', 'Det','MAP','MMPP2','NHPP',...
                ... % Non-Markovian renewal distributions: converted to acyclic PH
                ... % by sn_nonmarkov_toph in runAnalyzer, so the fluid ODE solves them.
                'Gamma','Lognormal','Pareto','Uniform','Weibull',...
                'StatelessClassSwitcher','InfiniteServer','SharedServer','Buffer','Dispatcher',...
                'Server','ServiceTunnel',...
                'SchedStrategy_INF','SchedStrategy_PS',...
                'SchedStrategy_DPS','SchedStrategy_FCFS','SchedStrategy_HOL',...
                'SchedStrategy_SIRO','SchedStrategy_LCFS','SchedStrategy_LCFSPR',...
                'ReplacementStrategy_RR','ReplacementStrategy_FIFO',... % RAND(m)/FIFO(m) refined mean field
                'ReplacementStrategy_SFIFO',... % strict FIFO(m) position-resolved mean field; LRU/HLRU/CLIMB/QLRU rejected at runtime
                ...
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
            libs = {};
            if nargin >= 2 && isfield(struct(options), 'method')
                if any(strcmp(options.method, {'rmf', 'fluid.rmf'}))
                    libs{end+1} = 'rmf_tool';
                end
            end
        end
    end
end
