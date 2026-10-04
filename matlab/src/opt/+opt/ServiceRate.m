classdef ServiceRate < opt.DecisionVariable
    % ServiceRate  Optimize the exponential processing rate of a station for a
    % job class. Continuous (differentiable): exposes paramKey/decodeJacobian
    % for the analytic-gradient path. Mirrors native-Python ServiceRate.

    properties
        station
        jobclass
        minRate
        maxRate
    end

    methods
        function obj = ServiceRate(station, jobclass, bounds, name)
            if nargin < 4 || isempty(name)
                name = [station.getName() '_' jobclass.getName() '_rate'];
            end
            obj@opt.DecisionVariable(name);
            obj.station = station;
            obj.jobclass = jobclass;
            obj.minRate = bounds(1);
            obj.maxRate = bounds(2);
        end

        function s = getStation(obj), s = obj.station; end
        function c = getJobClass(obj), c = obj.jobclass; end
        function b = getBounds(obj), b = opt.DecisionVariable.unitBounds(1); end

        function v = decode(obj, x)
            v = obj.minRate + x(1) * (obj.maxRate - obj.minRate);
        end

        function apply(obj, model, value)
            cls = opt.DecisionVariable.resolveClass(model, obj.jobclass);
            if isempty(cls), cls = obj.jobclass; end
            nodes = model.getNodes();
            for i = 1:numel(nodes)
                if strcmp(nodes{i}.getName(), obj.station.getName())
                    nodes{i}.setService(cls, Exp(value));
                    return;
                end
            end
        end

        function t = getVariableType(obj), t = 'service_rate'; end

        function k = paramKey(obj)
            k = opt.SensitivityData.paramKey(obj.station.getName(), obj.jobclass.getName());
        end
        function j = decodeJacobian(obj, ~)
            j = obj.maxRate - obj.minRate;
        end
    end
end
