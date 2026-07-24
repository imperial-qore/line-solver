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
        
        function tf = supportsExactSensitivity(self) %#ok<MANU>
            % TF = SUPPORTSEXACTSENSITIVITY()
            % The normalizing-constant solver is exact on the same
            % product-form class that pfqn_sens differentiates, so
            % getSensitivityTable uses the analytic branch.
            tf = true;
        end

        function [allMethods] = listValidMethods(self)
            % allMethods = LISTVALIDMETHODS()
            % List valid methods for this solver
            sn = self.model.getStruct();
            allMethods = {'default','exact','erlangfp','mci','imci','ls',...
                'le','mmint2','gleint','panacea','ca',...
                'clw','kt','sampling','is',...
                'propfair','comom','cub',...
                'rd', 'nrp','nrl','gm','mem'};
        end

        function method = resolveMethod(self, options)
            % Feature-driven resolution of options.method='default'. An open
            % network with non-Markovian (any non-unit SCV) variability within
            % the MEM feature set is solved by the Maximum Entropy Method by
            % default, since the normalizing-constant path would silently
            % exponentialize it; plain Markovian models keep the exact
            % product-form path. Mirrors the dispatch in runAnalyzer.
            method = options.method;
            if strcmp(options.method, 'default')
                sn = self.model.getStruct();
                if solver_nc_mem_supports(sn)
                    scvv = sn.scv(isfinite(sn.scv));
                    if ~isempty(scvv) && any(abs(scvv - 1) > GlobalConstants.FineTol)
                        method = 'mem';
                    end
                end
            end
        end

        function featSupported = getMethodFeatureSet(self, method) %#ok<INUSD>
            % All NC methods share the solver-level feature envelope.
            %
            % Defining this is what lets NetworkSolver.supportsModelMethod name
            % the offending features: with no method feature set it falls back
            % to the coarse supports(model) and returns an empty reason, so the
            % gate could only report "features not supported" without saying
            % which ones.
            %
            % A non-Network model (e.g. a LayeredNetwork) has no
            % getUsedLangFeatures, so it keeps the coarse path and the
            % structural checks/redirects that operate on such models.
            if ~isa(self.model, 'Network')
                featSupported = [];
                return;
            end
            featSupported = SolverNC.getFeatureSet();
        end

        function [bool, reason] = supportsModelMethod(self, method)
            % MEM (Kouvatsos maximum entropy) has structural applicability
            % rules beyond a flat feature set (open-only, no class switching,
            % non-priority scheduling); delegate to solver_nc_mem_supports,
            % which returns a precise reason. All other NC methods inherit the
            % coarse product-form feature gate.
            memblocking = false;
            if strcmp(method, 'mem')
                sn = self.model.getStruct();
                [bool, reason, memblocking] = solver_nc_mem_supports(sn);
                if bool
                    % a mem-admissible model still has to clear the feature gate;
                    % see _kb/06-solver-catalog.md (mem.blocking note)
                    [bool, reason] = supportsModelMethod@NetworkSolver(self, method);
                end
            else
                [bool, reason] = supportsModelMethod@NetworkSolver(self, method);
            end
            % structural finite-capacity gate (no registry feature name); mem is
            % the one exception (it represents the buffer as a GE/GE/c/0;N queue).
            % see _kb/06-solver-catalog.md (finite capacity gate + mem.blocking)
            if bool && isa(self.model, 'Network') && ~memblocking
                [bool, reason] = NetworkSolver.checkBindingCapacity(self.model, 'SolverNC');
            end
        end

        function bool = isStochasticMethod(self, method) %#ok<INUSL>
            % BOOL = ISSTOCHASTICMETHOD(METHOD)
            % NC is deterministic except for the Monte Carlo integration
            % methods (mci/imci), logistic sampling (ls), the importance
            % sampling method (is), and the sampling method, whose estimates
            % depend on the random seed. Method names are tokenized so that
            % runtime-resolved names such as 'default/imci' and prefixed names
            % such as 'nc.ls' classify correctly.
            tokens = regexp(lower(method), '[./]', 'split');
            bool = any(ismember(tokens, {'mci','imci','ls','sampling','is'}));
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
                ... % Finite capacity regions: NC solves the OPEN single-Delay
                ... % loss-network case exactly (Erlang fixed point,
                ... % solver_nc_lossn_analyzer). It cannot do an FCR on queueing
                ... % stations -- a boolean feature cannot express that split, so
                ... % runAnalyzer keeps the imperative check for the queueing case,
                ... % the same pattern as SolverMVA.supportsFiniteCapacity.
                'Region', ...
                'Server','JobSink','RandomSource','ServiceTunnel',...
                'SchedStrategy_INF','SchedStrategy_PS','SchedStrategy_SIRO',...
                'SchedStrategy_LCFS','SchedStrategy_LCFSPR',...
                'RoutingStrategy_PROB','RoutingStrategy_RAND',...
                'SchedStrategy_FCFS',...
                'Fork','Join','Forker','Joiner',... % fork-join via the MMT transformation (fjFixedPoint)
                'ClosedClass','SelfLoopingClass',...
                'Cache','CacheClassSwitcher','OpenClass',            ...
                'ReplacementStrategy_RR', 'ReplacementStrategy_FIFO',...
                'ReplacementStrategy_HLRU',...
                'LoadDependence','ClassDependence','JointDependence'});
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
