classdef MinimizeSystemResponseTime < opt.Objective
    % MinimizeSystemResponseTime  Minimize end-to-end (system) response time.
    % Mirrors native-Python MinimizeSystemResponseTime.

    properties
        jobclass = '';
    end

    methods
        function obj = MinimizeSystemResponseTime(jobclass, subjectTo)
            if nargin >= 1 && ~isempty(jobclass)
                obj.jobclass = opt.Constraint.nameOf(jobclass);
            end
            if nargin >= 2 && ~isempty(subjectTo)
                obj.constraints = subjectTo;
            end
        end
        function tf = isMinimization(~), tf = true; end
        function v = evaluate(obj, result, ~)
            if ~isempty(obj.jobclass)
                v = result.getSystemResponseTime(obj.jobclass);
            else
                v = result.getSystemResponseTime();
            end
        end
    end
end
