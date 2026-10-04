classdef OptimizationResult < handle
    % OptimizationResult  Result from a single optimization run. Mirrors
    % native-Python line_solver.opt.results.OptimizationResult.

    properties
        objectiveValue = inf;
        variableValues            % Map name -> value (scalar or vector)
        constraintViolations      % Map name -> violation
        feasible = false;
        iterations = 0;
        solveTime = 0.0;
        modelEvaluations = 0;
        convergenceHistory = [];
        terminatedBy = '';
    end

    methods
        function obj = OptimizationResult()
            obj.variableValues = configureDictionary('string', 'cell');
            obj.constraintViolations = configureDictionary('string', 'double');
        end

        function tf = isFeasible(obj)
            tf = obj.feasible;
        end
        function v = getObjectiveValue(obj)
            v = obj.objectiveValue;
        end
        function v = getVariableValue(obj, name)
            if isKey(obj.variableValues, name), v = obj.variableValues{name}; else, v = []; end
        end
        function v = getConstraintViolation(obj, name)
            if isKey(obj.constraintViolations, name), v = obj.constraintViolations(name); else, v = 0.0; end
        end
        function v = getTotalViolation(obj)
            if numEntries(obj.constraintViolations) == 0, v = 0.0;
            else, v = sum(values(obj.constraintViolations)); end
        end
    end
end
