classdef LineOptSolverOptions < handle
    % LineOptSolverOptions  Configuration for opt.LineOptSolver, mirroring
    % native-Python LineOptSolver.defaultOptions(). Chained setters return self.

    properties
        strategy = 'best1bin';
        popsize = 15;
        mutationLow = 0.5;
        mutationHigh = 1.0;
        recombination = 0.7;
        tol = 0.01;
        maxIterations = 100;
        timeLimit = 300.0;
        seed = [];
        verbose = false;
        penaltyWeight = 1e6;
        scenarioAggregation = 'worst';
        optimizer = 'evolution';
        fdStep = 1e-6;
        gradientRestarts = 4;
        % LayeredNetwork (LQN) gradient source (used only when the model is a
        % LayeredNetwork and the gradient path is taken):
        %   'fd'              finite-difference the whole LayeredNetwork per
        %                     parameter -- correct total derivative, robust
        %                     default.
        %   'partial_sens'    assemble the direction from SolverLN's per-layer
        %                     WITHIN-LAYER partial derivatives -- cheap, biased.
        %   'partial_plus_fd' partial_sens direction, corrected by a full
        %                     LayeredNetwork finite difference every fdRefresh
        %                     gradient evaluations.
        lqnGradient = 'fd';
        fdRefresh = 5;          % partial_plus_fd full-FD correction period
        % Explicit layer freezing (LQN only): a cell array of layer names (host
        % or task layers) whose variables are held at the model's current value
        % instead of being optimized. [] freezes nothing.
        frozenLayers = [];
    end

    methods
        function obj = setStrategy(obj, v), obj.strategy = v; end
        function obj = setPopsize(obj, v), obj.popsize = v; end
        function obj = setMutation(obj, low, high), obj.mutationLow = low; obj.mutationHigh = high; end
        function obj = setRecombination(obj, v), obj.recombination = v; end
        function obj = setTolerance(obj, v), obj.tol = v; end
        function obj = setMaxIterations(obj, v), obj.maxIterations = v; end
        function obj = setTimeLimit(obj, v), obj.timeLimit = v; end
        function obj = setSeed(obj, v), obj.seed = v; end
        function obj = setVerbose(obj, v), obj.verbose = v; end
        function obj = setPenaltyWeight(obj, v), obj.penaltyWeight = v; end
        function obj = setScenarioAggregation(obj, v), obj.scenarioAggregation = v; end
        function obj = setOptimizer(obj, v), obj.optimizer = v; end
        function obj = setFdStep(obj, v), obj.fdStep = v; end
        function obj = setGradientRestarts(obj, v), obj.gradientRestarts = v; end
        function obj = setLqnGradient(obj, v)
            if ~any(strcmp(v, {'fd','partial_sens','partial_plus_fd'}))
                line_error(mfilename, ['lqnGradient must be ''fd'', ' ...
                    '''partial_sens'', or ''partial_plus_fd''']);
            end
            obj.lqnGradient = v;
        end
        function obj = setFdRefresh(obj, v), obj.fdRefresh = round(v); end
        function obj = setFrozenLayers(obj, v)
            if isempty(v), obj.frozenLayers = []; else, obj.frozenLayers = v; end
        end
    end
end
