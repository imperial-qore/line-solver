classdef Constraint < handle
    % Constraint  Abstract base for line-opt constraints. Each computes a
    % non-negative violation (0 if satisfied). Mirrors native-Python
    % line_solver.opt.objectives.Constraint.

    properties
        name = '';
    end

    methods
        function obj = Constraint(name)
            if nargin >= 1 && ~isempty(name)
                obj.name = name;
            end
        end
        function n = getName(obj)
            if isempty(obj.name)
                n = obj.generateName();
            else
                n = obj.name;
            end
        end
        function tf = isSatisfied(obj, result, variableValues, tol)
            if nargin < 4, tol = 1e-6; end
            tf = obj.evaluate(result, variableValues) <= tol;
        end
    end

    methods (Abstract)
        n = generateName(obj)
        v = evaluate(obj, result, variableValues)
    end

    methods (Static)
        function v = upperBoundViolation(actual, bound)
            if ~isfinite(actual)
                v = inf;
            else
                v = max(0.0, actual - bound);
            end
        end
        function v = lowerBoundViolation(actual, bound)
            if ~isfinite(actual)
                v = inf;
            else
                v = max(0.0, bound - actual);
            end
        end
        function nm = nameOf(x)
            if ischar(x) || isstring(x)
                nm = char(x);
            elseif isempty(x)
                nm = '';
            else
                nm = x.getName();
            end
        end
    end
end
