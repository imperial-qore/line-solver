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
                'mva','exact','amva','qna', ...
                'qdlin','amva.qdlin', ...
                'bs','amva.bs', ...
                'sqni', ...
                'qd','amva.qd', ...
                'qli','amva.qli', ...
                'fli','amva.fli', ...
                'ab','amva.ab', ...
                'schmidt','amva.schmidt', ...
                'schmidt-ext','amva.schmidt-ext', ...
                'lin','egflin','gflin','amva.lin'};

            if ~sn_is_open_model(sn) && sn.nclasses == 1
                bounds = {'aba.upper','aba.lower','bjb.upper','bjb.lower', ...
                    'gb.upper','gb.lower','pb.upper','pb.lower','sb.upper','sb.lower'};
                allMethods = {allMethods{:}, bounds{:}}; %#ok<CCAT>
            end

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

    end

    methods(Static)
        function featSupported = getFeatureSet()
            % FEATSUPPORTED = GETFEATURESET()

            featSupported = SolverFeatureSet;
            featSupported.setTrue({'Sink','Source',...
                'ClassSwitch','Delay','DelayStation','Queue',...
                'APH','Coxian','Erlang','Exp','HyperExp',...
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
                'ClosedClass','SelfLoopingClass','OpenClass','Replayer'});
        end

        function [bool, featSupported] = supports(model)
            % [BOOL, FEATSUPPORTED] = SUPPORTS(MODEL)

            featUsed = model.getUsedLangFeatures();
            featSupported = SolverMVA.getFeatureSet();
            bool = SolverFeatureSet.supports(featSupported, featUsed);
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
