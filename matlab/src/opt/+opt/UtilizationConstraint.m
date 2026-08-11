classdef UtilizationConstraint < opt.Constraint
    % UtilizationConstraint  Utilization constraint: U <= maxValue.
    properties
        station
        maxValue
    end
    methods
        function obj = UtilizationConstraint(station, maxValue, name)
            if nargin < 3, name = ''; end
            obj@opt.Constraint(name);
            obj.station = opt.Constraint.nameOf(station);
            obj.maxValue = maxValue;
        end
        function n = generateName(obj)
            n = sprintf('Util_%s_le_%g', obj.station, obj.maxValue);
        end
        function v = evaluate(obj, result, ~)
            v = opt.Constraint.upperBoundViolation(result.getUtilization(obj.station), obj.maxValue);
        end
    end
end
