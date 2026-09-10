classdef ServerAllocation < opt.DecisionVariable
    % ServerAllocation  Optimize the number of servers at a station. Encodes an
    % integer server count in [minServers, maxServers]. Mirrors native-Python
    % ServerAllocation.

    properties
        station
        minServers
        maxServers
    end

    methods
        function obj = ServerAllocation(station, bounds, name)
            if nargin < 3 || isempty(name)
                name = [station.getName() '_servers'];
            end
            obj@opt.DecisionVariable(name);
            obj.station = station;
            obj.minServers = round(bounds(1));
            obj.maxServers = round(bounds(2));
        end

        function s = getStation(obj), s = obj.station; end
        function b = getBounds(obj), b = opt.DecisionVariable.unitBounds(1); end

        function v = decode(obj, x)
            continuous = obj.minServers + x(1) * (obj.maxServers - obj.minServers);
            v = min(max(round(continuous), obj.minServers), obj.maxServers);
        end

        function apply(obj, model, value)
            nodes = model.getNodes();
            for i = 1:numel(nodes)
                if strcmp(nodes{i}.getName(), obj.station.getName())
                    nodes{i}.setNumberOfServers(value);
                    return;
                end
            end
        end

        function t = getVariableType(obj), t = 'server_allocation'; end
    end
end
