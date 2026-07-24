classdef ThroughputConstraint < opt.Constraint
    % ThroughputConstraint  Throughput constraint: Tput >= minValue.
    properties
        station
        jobclass = '';
        minValue
    end
    methods
        function obj = ThroughputConstraint(station, jobclass, minValue, name)
            if nargin < 4, name = ''; end
            obj@opt.Constraint(name);
            obj.station = opt.Constraint.nameOf(station);
            if nargin >= 2 && ~isempty(jobclass)
                obj.jobclass = opt.Constraint.nameOf(jobclass);
            end
            obj.minValue = minValue;
        end
        function n = generateName(obj)
            if ~isempty(obj.jobclass)
                n = sprintf('Tput_%s_%s_ge_%g', obj.station, obj.jobclass, obj.minValue);
            else
                n = sprintf('Tput_%s_ge_%g', obj.station, obj.minValue);
            end
        end
        function v = evaluate(obj, result, ~)
            if ~isempty(obj.jobclass)
                actual = result.getThroughput(obj.station, obj.jobclass);
            else
                actual = result.getThroughput(obj.station);
            end
            v = opt.Constraint.lowerBoundViolation(actual, obj.minValue);
        end
    end
end
