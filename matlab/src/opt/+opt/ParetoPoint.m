classdef ParetoPoint < handle
    % ParetoPoint  One point of a cost-performance tradeoff curve produced by
    % opt.ParetoSweep. Mirrors native-Python ParetoPoint.
    properties
        epsilon = 0.0;
        objectiveValue = inf;
        feasible = false;
        result = [];
    end
    methods
        function obj = ParetoPoint(epsilon, objectiveValue, feasible, result)
            if nargin >= 1, obj.epsilon = epsilon; end
            if nargin >= 2, obj.objectiveValue = objectiveValue; end
            if nargin >= 3, obj.feasible = feasible; end
            if nargin >= 4, obj.result = result; end
        end
    end
end
