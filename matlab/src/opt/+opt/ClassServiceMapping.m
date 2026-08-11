classdef ClassServiceMapping < opt.DecisionVariable
    % ClassServiceMapping  Optimize the class-to-station mapping by rerouting a
    % job class through a selected candidate station and bypassing the others,
    % preserving default routing of every other class. Mirrors native-Python
    % ClassServiceMapping.

    properties
        jobclass
        stations     % cell array of Station
    end

    methods
        function obj = ClassServiceMapping(jobclass, stations, name)
            if nargin < 3 || isempty(name)
                name = [jobclass.getName() '_mapping'];
            end
            obj@opt.DecisionVariable(name);
            obj.jobclass = jobclass;
            obj.stations = stations;
        end

        function c = getJobClass(obj), c = obj.jobclass; end
        function s = getStations(obj), s = obj.stations; end
        function b = getBounds(obj), b = opt.DecisionVariable.unitBounds(1); end

        function v = decode(obj, x)
            n = numel(obj.stations);
            idx = floor(x(1) * n);
            v = min(idx, n - 1);   % 0-based index
        end

        function apply(obj, model, value)
            if isempty(obj.stations), return; end
            v = min(max(value, 0), numel(obj.stations) - 1);

            nodes = model.getNodes();
            n = numel(nodes);
            conn = opt.DecisionVariable.connectionMatrix(model);

            selectedName = obj.stations{v + 1}.getName();
            blocked = [];
            for k = 1:numel(obj.stations)
                if ~strcmp(obj.stations{k}.getName(), selectedName)
                    idx = opt.DecisionVariable.indexOfNode(nodes, obj.stations{k}.getName());
                    if idx > 0, blocked(end+1) = idx; end %#ok<AGROW>
                end
            end

            mapped = opt.DecisionVariable.resolveClass(model, obj.jobclass);
            if isempty(mapped), return; end

            rt = model.initRoutingMatrix();
            classes = model.getClasses();
            for ci = 1:numel(classes)
                c = classes{ci};
                adj = conn;
                if strcmp(c.getName(), mapped.getName())
                    adj(:, blocked) = 0.0;
                end
                for i = 1:n
                    outDegree = sum(adj(i, :));
                    if outDegree <= 0, continue; end
                    for j = 1:n
                        if adj(i, j) > 0
                            rt.set(c, c, nodes{i}, nodes{j}, adj(i, j) / outDegree);
                        end
                    end
                end
            end
            model.link(rt);
        end

        function t = getVariableType(obj), t = 'class_mapping'; end
    end
end
