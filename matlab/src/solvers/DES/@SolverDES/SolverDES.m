classdef SolverDES < NetworkSolver
    % SolverDES Discrete Event Simulation solver using SSJ library
    %
    % SolverDES implements a discrete-event simulation solver that uses the SSJ
    % (Stochastic Simulation in Java) library to analyze queueing networks.
    % It supports open and closed networks with various service distributions,
    % scheduling strategies, and advanced node types.
    %
    % For LayeredNetwork (LQN) models, SolverDES also supports DES simulation.
    %
    % @brief Discrete-event simulation solver using SSJ library
    %
    % Example:
    % @code
    % solver = SolverDES(model, 'samples', 1000000, 'seed', 23000);
    % solver.getAvg();  % Run DES simulation
    % @endcode
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods
        function self = SolverDES(model, varargin)
            % SOLVERDES Create a DES solver instance
            %
            % @brief Creates a Discrete Event Simulation solver
            % @param model Network model to be analyzed via DES simulation
            % @param varargin Optional parameters (samples, seed, method, etc.)
            % @return self SolverDES instance configured for simulation

            self@NetworkSolver(model, mfilename);
            self.setOptions(Solver.parseOptions(varargin, self.defaultOptions));
            % SolverDES only supports Java backend - force lang to 'java'
            self.options.lang = 'java';
            self.setLang();
        end

        function sn = getStruct(self)
            % QN = GETSTRUCT()

            % Get data structure summarizing the model
            sn = self.model.getStruct(true);
        end

        function [allMethods] = listValidMethods(self)
            % allMethods = LISTVALIDMETHODS()
            % List valid methods for this solver

            allMethods = {'default'};
        end
    end

    methods (Static)

        function featSupported = getFeatureSet()
            % FEATSUPPORTED = GETFEATURESET()

            featSupported = SolverFeatureSet;
            featSupported.setTrue({'Sink', 'Source', ...
                'Queue', 'Delay', ...
                'Fork', 'Join', 'Forker', 'Joiner', ... % Fork-Join node support
                'Linkage', 'Enabling', 'Timing', 'Firing', 'Storage', ... % Petri net section support
                'Logger', 'LogTunnel', ... % Logger node support
                'Buffer', ...  % Finite buffer capacity support
                'Region', ...  % Finite capacity region support
                'Exp', 'Erlang', 'HyperExp', 'PH', 'APH', 'Coxian', 'Cox2', 'MAP', 'MMPP2', 'Replayer', ...
                'Det', 'Uniform', 'Gamma', 'Pareto', 'Weibull', 'Lognormal', ... % Additional continuous distributions
                'Server', 'RandomSource', ...
                'SchedStrategy_FCFS', 'SchedStrategy_INF', ...
                'SchedStrategy_HOL', ... % Priority scheduling (FCFS with priorities)
                'SchedStrategy_PS', ... % Processor Sharing
                'SchedStrategy_DPS', ... % Discriminatory Processor Sharing
                'SchedStrategy_GPS', ... % Generalized Processor Sharing
                'SchedStrategy_LCFS', ... % Last Come First Served (non-preemptive)
                'SchedStrategy_LCFSPR', ... % LCFS Preemptive Resume
                'SchedStrategy_LCFSPI', ... % LCFS Preemptive Independent
                'SchedStrategy_SIRO', ...
                'SchedStrategy_SJF', 'SchedStrategy_LJF', ...
                'SchedStrategy_LEPT', ...
                'SchedStrategy_SEPT', ...
                'SchedStrategy_SRPT', ... % Shortest Remaining Processing Time (preemptive)
                'SchedStrategy_SRPTPRIO', ... % SRPT with priorities
                'SchedStrategy_PSJF', ... % Preemptive Shortest Job First
                'SchedStrategy_FB', ... % Feedback / Least Attained Service
                'SchedStrategy_LRPT', ... % Longest Remaining Processing Time
                'SchedStrategy_EXT', ...
                'SchedStrategy_POLLING', ... % Polling scheduling (GATED, EXHAUSTIVE, KLIMITED)
                'Router', ... % Router node support
                'ClassSwitch', 'StatelessClassSwitcher', ... % Class switching node support
                'Cache', 'CacheClassSwitcher', ... % Cache node support with replacement policies (LRU, FIFO, Strict FIFO, RR)
                'RoutingStrategy_PROB', 'RoutingStrategy_RAND', ...
                'RoutingStrategy_RROBIN', 'RoutingStrategy_WRROBIN', ...
                'RoutingStrategy_KCHOICES', ... % Power of K Choices routing
                'OpenClass', ...
                'ClosedClass', ...
                'ReplacementStrategy_RR', 'ReplacementStrategy_FIFO', 'ReplacementStrategy_SFIFO', 'ReplacementStrategy_LRU'});
        end

        function [bool, featSupported] = supports(model)
            % [BOOL, FEATSUPPORTED] = SUPPORTS(MODEL)

            if isa(model, 'LayeredNetwork')
                % LayeredNetwork not supported - use SolverLDES from line-apps
                bool = false;
                featSupported = SolverFeatureSet;
                line_warning(mfilename, 'SolverDES does not support LayeredNetwork models. Use SolverLDES from line-apps.');
                return;
            end

            % Regular Network support
            featUsed = model.getUsedLangFeatures();
            featSupported = SolverDES.getFeatureSet();
            bool = SolverFeatureSet.supports(featSupported, featUsed);
        end

        function options = defaultOptions()
            % OPTIONS = DEFAULTOPTIONS()

            options = SolverOptions('DES');
        end

    end
end
