classdef SubProblem < handle
    % SubProblem  A subset of variables to optimize while others are fixed.
    % Mirrors native-Python line_solver.opt.decomposition.SubProblem.
    properties
        name
        variableType
        variables       % cell of opt.DecisionVariable
        fixedValues
    end
    methods
        function obj = SubProblem(name, variableType, variables)
            obj.name = name;
            obj.variableType = variableType;
            obj.variables = variables;
            obj.fixedValues = containers.Map('KeyType','char','ValueType','any');
        end
        function names = getVariableNames(obj)
            names = cell(1, numel(obj.variables));
            for i = 1:numel(obj.variables), names{i} = obj.variables{i}.getName(); end
        end
    end
end
