classdef BudgetConstraint < opt.Constraint
    % BudgetConstraint  Budget constraint: total cost <= budget.
    properties
        budget
        costCoefficients   % dictionary name -> cost
    end
    methods
        function obj = BudgetConstraint(budget, costCoefficients, name)
            if nargin < 3, name = ''; end
            obj@opt.Constraint(name);
            obj.budget = budget;
            if nargin >= 2 && ~isempty(costCoefficients)
                obj.costCoefficients = costCoefficients;
            else
                obj.costCoefficients = configureDictionary('string', 'double');
            end
        end
        function n = generateName(obj)
            n = sprintf('Budget_le_%g', obj.budget);
        end
        function c = computeCost(obj, variableValues)
            c = 0.0;
            ks = keys(obj.costCoefficients);
            for i = 1:numel(ks)
                if isKey(variableValues, ks(i))
                    c = c + obj.costCoefficients(ks(i)) * opt.Objective.numericValue(variableValues{ks(i)});
                end
            end
        end
        function v = evaluate(obj, ~, variableValues)
            v = max(0.0, obj.computeCost(variableValues) - obj.budget);
        end
    end
end
