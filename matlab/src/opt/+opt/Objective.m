classdef Objective < handle
    % Objective  Abstract base for line-opt objectives. Mirrors native-Python
    % line_solver.opt.objectives.Objective: defines the scalar to minimize,
    % with attached constraints folded in as penalties by evaluateWithPenalty.

    properties
        constraints = {};   % cell array of opt.Constraint
    end

    methods
        function c = getConstraints(obj)
            c = obj.constraints;
        end
        function v = evaluateWithPenalty(obj, result, variableValues, penaltyWeight)
            if nargin < 4, penaltyWeight = 1e6; end
            v = obj.evaluate(result, variableValues);
            for i = 1:numel(obj.constraints)
                v = v + obj.constraints{i}.evaluate(result, variableValues) * penaltyWeight;
            end
        end
    end

    methods (Abstract)
        v = evaluate(obj, result, variableValues)
        tf = isMinimization(obj)
    end

    methods (Static)
        function d = numericValue(value)
            if isempty(value)
                d = 0.0;
            elseif isnumeric(value)
                d = sum(value(:));
            else
                d = 0.0;
            end
        end
        function tf = isScalarNumeric(value)
            tf = isnumeric(value) && isscalar(value);
        end
    end
end
