classdef ResponseTimeConstraint < opt.Constraint
    % ResponseTimeConstraint  Per-station response time constraint: RT <= maxValue.
    properties
        station
        jobclass = '';
        maxValue
    end
    methods
        function obj = ResponseTimeConstraint(station, jobclass, maxValue, name)
            if nargin < 4, name = ''; end
            obj@opt.Constraint(name);
            obj.station = opt.Constraint.nameOf(station);
            if nargin >= 2 && ~isempty(jobclass)
                obj.jobclass = opt.Constraint.nameOf(jobclass);
            end
            obj.maxValue = maxValue;
        end
        function n = generateName(obj)
            if ~isempty(obj.jobclass)
                n = sprintf('RT_%s_%s_le_%g', obj.station, obj.jobclass, obj.maxValue);
            else
                n = sprintf('RT_%s_le_%g', obj.station, obj.maxValue);
            end
        end
        function v = evaluate(obj, result, ~)
            if ~isempty(obj.jobclass)
                actual = result.getResponseTime(obj.station, obj.jobclass);
            else
                actual = result.getResponseTime(obj.station);
            end
            v = opt.Constraint.upperBoundViolation(actual, obj.maxValue);
        end
    end
end
