classdef SubProblemResult < handle
    % SubProblemResult  Result from solving one decomposition subproblem.
    % Mirrors native-Python line_solver.opt.results.SubProblemResult.

    properties
        name = '';
        result                    % opt.OptimizationResult
        variablesFixed            % Map name -> value
    end

    methods
        function obj = SubProblemResult(name, result)
            obj.variablesFixed = configureDictionary('string', 'cell');
            if nargin >= 1, obj.name = name; end
            if nargin >= 2, obj.result = result; else, obj.result = opt.OptimizationResult(); end
        end
    end
end
