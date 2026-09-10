classdef StationReplicas < opt.DecisionVariable
    % StationReplicas  Optimize the number of identical station copies. N
    % replicas are represented as one multiserver station with N times the base
    % server count, keeping topology and names fixed. Mirrors native-Python
    % StationReplicas.

    properties
        station
        minReplicas
        maxReplicas
    end

    methods
        function obj = StationReplicas(station, bounds, name)
            if nargin < 3 || isempty(name)
                name = [station.getName() '_replicas'];
            end
            obj@opt.DecisionVariable(name);
            obj.station = station;
            obj.minReplicas = round(bounds(1));
            obj.maxReplicas = round(bounds(2));
        end

        function s = getStation(obj), s = obj.station; end
        function b = getBounds(obj), b = opt.DecisionVariable.unitBounds(1); end

        function v = decode(obj, x)
            continuous = obj.minReplicas + x(1) * (obj.maxReplicas - obj.minReplicas);
            v = min(max(round(continuous), obj.minReplicas), obj.maxReplicas);
        end

        function apply(obj, model, value)
            nodes = model.getNodes();
            for i = 1:numel(nodes)
                if strcmp(nodes{i}.getName(), obj.station.getName())
                    base = nodes{i}.getNumberOfServers();
                    if ~isfinite(base) || base < 1
                        base = 1;
                    end
                    nodes{i}.setNumberOfServers(base * value);
                    return;
                end
            end
        end

        function t = getVariableType(obj), t = 'station_replicas'; end
    end
end
