classdef SolverMVA < NetworkSolver
    % Mean Value Analysis solver for queueing networks
    %
    % Implements MVA algorithms for analyzing closed and open queueing networks.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods
        function self = SolverMVA(model,varargin)
            % SOLVERMVA Create an MVA solver instance
            %
            % @brief Creates a Mean Value Analysis solver for the given model
            % @param model Network model to be analyzed
            % @param varargin Optional solver options (method, tolerance, etc.)
            % @return self SolverMVA instance configured with specified options

            self@NetworkSolver(model, mfilename);
            self.setOptions(Solver.parseOptions(varargin, SolverMVA.defaultOptions));
            self.setLang();
        end

        function sn = getStruct(self)
            % GETSTRUCT Get model data structure for analysis
            %
            % @brief Returns the internal data structure representing the model
            % @return sn Structured data representing the queueing network
            sn = self.model.getStruct(false);
        end

        function tf = supportsExactSensitivity(self) %#ok<MANU>
            % TF = SUPPORTSEXACTSENSITIVITY()
            % MVA differentiates its own recursion: getSensitivityTable uses
            % the analytic branch (pfqn_sens) wherever the model is in scope.
            tf = true;
        end

        [runtime, analyzer] = runAnalyzer(self, options);
        [lNormConst] = getProbNormConstAggr(self);
        [Pnir,logPnir] = getProbAggr(self, ist);
        [Pnir,logPn] = getProbSysAggr(self);

        function [allMethods] = listValidMethods(self)
            % LISTVALIDMETHODS Get all valid MVA solution methods
            %
            % @brief Returns cell array of valid MVA methods for the current model
            % @return allMethods Cell array of method names available for this model

            sn = self.model.getStruct;
            % base set of methods
            allMethods = {'default',...
                'mva','exact','amva', ...
                'sum','esum', ...
                'qdlin','amva.qdlin', ...
                'bs','amva.bs', ...
                'sqni', ...
                'qd','amva.qd', ...
                'qli','amva.qli', ...
                'fli','amva.fli', ...
                'ab','amva.ab', ...
                'schmidt','amva.schmidt', ...
                'schmidt-ext','amva.schmidt-ext'};

            % QNA/RQNA advertised for fully open models only; see
            % _kb/06-solver-catalog.md (MVA section) for the default-dispatch rules
            if sn_is_open_model(sn)
                allMethods = [allMethods(1:4), {'qna','rqna'}, allMethods(5:end)];
            end

            % SQD (Smith Queue Decomposition) is only valid for closed
            % single-chain Blocking-After-Service networks; kept at this
            % position to preserve BAS regression-baseline ordering.
            if sn_is_bas_model(sn)
                allMethods{end+1} = 'sqd'; %#ok<AGROW>
            end

            % MVAC (exact mean value analysis by chain, pfqn_mvac): closed
            % single-server product-form networks only; rejects open/mixed and
            % multiserver at solve time.
            if ~any(isinf(sn.njobs))
                allMethods{end+1} = 'mvac'; %#ok<AGROW>
            end

            % Marie withheld for open models and for class-dependent routing;
            % see _kb/06-solver-catalog.md (MVA section) for the rationale
            if ~sn_is_open_model(sn) && ~sn_has_classdep_routing(sn)
                allMethods = [allMethods, {'marie','amva.marie'}]; %#ok<AGROW>
            end

            allMethods = {allMethods{:}, ...
                'lin','egflin','gflin','amva.lin'}; %#ok<CCAT>

            % Bound methods (aba/bjb/gb/pb/sb/mwba) are NOT listed here: they
            % moved to SolverBA, and runAnalyzer raises line_error for the whole
            % family. Listing them made this a false claim, since every listed
            % name errored when run. See SolverBA.listValidMethods.

            if sn_is_open_model(sn) && sn.nstations == 2 && sn.nclasses == 1
                % methods to add for queueing systems
                qsys = {'mm1','mmk','mg1','mgi1','gm1','gig1','gim1','gig1.kingman', ...
                    'gigk','gigk.kingman_approx', ...
                    'gig1.gelenbe','gig1.heyman','gig1.kimura','gig1.allen', ...
                    'gig1.kobayashi','gig1.klb','gig1.marchal'};
                % append, keeping original order and avoiding duplicates
                allMethods = {allMethods{:}, qsys{:}}; %#ok<CCAT>
            end
        end

        function method = resolveMethod(self, options)
            % Feature-driven resolution of options.method='default'. A bursty
            % single-class open network has a non-renewal (MAP/MMPP) arrival
            % process whose autocorrelation a two-moment method cannot capture,
            % so the default dispatch selects RQNA (robust queueing network
            % analyzer, indices of dispersion). All other cases keep 'default',
            % which the analyzer expands with its own heuristics. Expressed as
            % RQNA's MAP-family coverage rather than a bespoke gate.
            method = options.method;
            if strcmp(options.method, 'default')
                sn = self.model.getStruct();
                if (sn.nclasses == 1) && all(isinf(sn.njobs)) && sn_has_bursty_arrival(sn)
                    method = 'rqna';
                end
            end
        end

        function featSupported = getMethodFeatureSet(self, method) %#ok<INUSL>
            % Per-method feature deltas applied to the base MVA envelope.
            % RQNA natively consumes non-renewal MAP/MMPP/RAP arrival and
            % service processes (open only); QNA is a two-moment open-network
            % method. The queueing-system and bounds methods are already
            % structurally restricted by listValidMethods, so they inherit the
            % base envelope unchanged.
            featSupported = SolverMVA.getFeatureSet();
            switch method
                case 'rqna'
                    featSupported.setTrue({'MAP','MMPP2','MMAP','RAP'});
                    featSupported.setFalse({'ClosedClass','SelfLoopingClass'});
                case 'qna'
                    featSupported.setFalse({'ClosedClass','SelfLoopingClass'});
            end
        end

        function [bool, reason] = supportsModelMethod(self, method)
            % Finite station/class capacity has no registry feature name, so
            % the coarse per-method feature gate cannot see it. Apply the
            % structural capacity check on top of it, otherwise MVA silently
            % returns the unconstrained product-form answer for models built
            % with setCapacity / a finite classCap (BUG-39).
            [bool, reason] = supportsModelMethod@NetworkSolver(self, method);
            if bool && isa(self.model, 'Network')
                [bool, reason] = SolverMVA.supportsFiniteCapacity(self.model);
            end
        end

    end

    methods(Static)
        function [bool, reason] = supportsFiniteCapacity(model)
            % [BOOL, REASON] = SUPPORTSFINITECAPACITY(MODEL)
            % MVA-specific finite-capacity gate: Blocking-After-Service models
            % are exempt because MVA offers the Smith queue-decomposition
            % method 'sqd', and solver_mva_analyzer routes a BAS model to
            % solver_sqd under the default method too, so the finite buffers
            % ARE honoured on every MVA path. Everything else defers to the
            % shared product-form gate.
            bool = true;
            reason = '';
            if ~isa(model, 'Network')
                return
            end
            if sn_is_bas_model(model.getStruct())
                return
            end
            [bool, reason] = NetworkSolver.checkBindingCapacity(model, 'SolverMVA');
        end

        function featSupported = getFeatureSet()
            % FEATSUPPORTED = GETFEATURESET()

            featSupported = SolverFeatureSet;
            featSupported.setTrue({'Sink','Source',...
                'ClassSwitch','Delay','DelayStation','Queue',...
                'APH','Coxian','Erlang','Exp','HyperExp','BMAP',...
                'Pareto','Weibull','Lognormal','Uniform','Det', ...
                'StatelessClassSwitcher','InfiniteServer','SharedServer','Buffer','Dispatcher',...
                'CacheClassSwitcher','Cache', ...
                'Server','JobSink','RandomSource','ServiceTunnel',...
                'SchedStrategy_INF','SchedStrategy_PS',...
                'SchedStrategy_DPS','SchedStrategy_FCFS','SchedStrategy_SIRO','SchedStrategy_HOL',...
                'SchedStrategy_LCFS','SchedStrategy_LCFSPR','SchedStrategy_POLLING',...
                'Fork','Forker','Join','Joiner',...
                'RoutingStrategy_PROB','RoutingStrategy_RAND',...
                'ReplacementStrategy_RR', 'ReplacementStrategy_FIFO', 'ReplacementStrategy_LRU',...
                'ReplacementStrategy_HLRU',...
                'MMAP',...  % marked MAP sources (cache LRU via cache_ttl_lrum_map)
                'ClosedClass','SelfLoopingClass','OpenClass','Replayer',...
                'LoadDependence','ClassDependence','JointDependence'});
        end

        function [bool, featSupported] = supports(model)
            % [BOOL, FEATSUPPORTED] = SUPPORTS(MODEL)

            featUsed = model.getUsedLangFeatures();
            featSupported = SolverMVA.getFeatureSet();
            bool = SolverFeatureSet.supports(featSupported, featUsed);
            if bool
                % Registry inclusion cannot see finite capacity (no feature
                % name); apply the structural gate too, so that the static
                % supports() agrees with the runAnalyzer feature gate.
                [bool, reason] = SolverMVA.supportsFiniteCapacity(model);
                if ~bool
                    line_warning(mfilename, '%s\n', reason);
                end
            end
        end

        function options = defaultOptions
            % OPTIONS = DEFAULTOPTIONS()

            options = SolverOptions('MVA');
        end

        function libs = getLibrariesUsed(sn, options)
            % GETLIBRARIESUSED Get list of external libraries used by MVA solver
            % MVA uses internal algorithms, no external library attribution needed
            libs = {};
        end

    end
end
