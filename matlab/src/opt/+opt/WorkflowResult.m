classdef WorkflowResult < handle
    % WorkflowResult  Result from a decomposed workflow optimization. Mirrors
    % native-Python line_solver.opt.results.WorkflowResult.

    properties
        finalObjective = inf;
        subproblemResults         % Map name -> opt.SubProblemResult
        cyclesCompleted = 0;
        converged = false;
        totalSolveTime = 0.0;
        objectiveHistory = [];
        finalVariableValues       % Map name -> value
        % LQN layer-wise decomposition diagnostics (solveLayered): the final
        % frozen layer set and the total LINE solve count.
        frozenLayers = {};
        modelEvaluations = 0;
    end

    methods
        function obj = WorkflowResult()
            obj.subproblemResults = configureDictionary('string', 'cell');
            obj.finalVariableValues = configureDictionary('string', 'cell');
        end

        function tf = isConverged(obj)
            tf = obj.converged;
        end
        function r = getSubProblemResult(obj, name)
            if isKey(obj.subproblemResults, name), r = obj.subproblemResults{name}; else, r = []; end
        end
        function v = getFinalVariableValue(obj, name)
            if isKey(obj.finalVariableValues, name), v = obj.finalVariableValues{name}; else, v = []; end
        end
    end
end
