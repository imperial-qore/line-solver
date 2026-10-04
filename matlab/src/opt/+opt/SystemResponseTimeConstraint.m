classdef SystemResponseTimeConstraint < opt.Constraint
    % SystemResponseTimeConstraint  End-to-end response time: SysRespT <= maxValue.
    properties
        jobclass = '';
        maxValue
    end
    methods
        function obj = SystemResponseTimeConstraint(jobclass, maxValue, name)
            if nargin < 3, name = ''; end
            obj@opt.Constraint(name);
            if nargin >= 1 && ~isempty(jobclass)
                obj.jobclass = opt.Constraint.nameOf(jobclass);
            end
            obj.maxValue = maxValue;
        end
        function n = generateName(obj)
            if ~isempty(obj.jobclass)
                n = sprintf('SysRT_%s_le_%g', obj.jobclass, obj.maxValue);
            else
                n = sprintf('SysRT_le_%g', obj.maxValue);
            end
        end
        function v = evaluate(obj, result, ~)
            if ~isempty(obj.jobclass)
                actual = result.getSystemResponseTime(obj.jobclass);
            else
                actual = result.getSystemResponseTime();
            end
            v = opt.Constraint.upperBoundViolation(actual, obj.maxValue);
        end
    end
end
