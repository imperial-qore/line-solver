classdef JobPopulation < opt.DecisionVariable
    % JobPopulation  Optimize the fixed circulating population of a closed
    % class. Encodes an integer count in [minJobs, maxJobs]. Mirrors
    % native-Python JobPopulation.

    properties
        jobclass
        minJobs
        maxJobs
    end

    methods
        function obj = JobPopulation(jobclass, bounds, name)
            if nargin < 3 || isempty(name)
                name = [jobclass.getName() '_population'];
            end
            obj@opt.DecisionVariable(name);
            obj.jobclass = jobclass;
            obj.minJobs = round(bounds(1));
            obj.maxJobs = round(bounds(2));
        end

        function c = getJobClass(obj), c = obj.jobclass; end
        function b = getBounds(obj), b = opt.DecisionVariable.unitBounds(1); end

        function v = decode(obj, x)
            continuous = obj.minJobs + x(1) * (obj.maxJobs - obj.minJobs);
            v = min(max(round(continuous), obj.minJobs), obj.maxJobs);
        end

        function apply(obj, model, value)
            cls = opt.DecisionVariable.resolveClass(model, obj.jobclass);
            if isempty(cls), return; end
            cls.setPopulation(value);
        end

        function t = getVariableType(obj), t = 'job_population'; end
    end
end
