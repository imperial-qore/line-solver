classdef MinimizeCost < opt.Objective
    % MinimizeCost  Minimize infrastructure cost subject to SLA constraints.
    % Cost = server costs ('name_servers') + rate costs (variable names
    % containing the station name and 'rate') + replica costs
    % ('name_replicas'). Mirrors native-Python MinimizeCost.
    %
    % serverCost/rateCost/replicaCost are containers.Map keyed by station name
    % (char) -> cost, or []. subjectTo is a cell array of opt.Constraint.

    properties
        serverCost
        rateCost
        replicaCost
    end

    methods
        function obj = MinimizeCost(serverCost, rateCost, replicaCost, subjectTo)
            obj.serverCost = opt.MinimizeCost.suffixMap(serverCost, '_servers');
            obj.rateCost = opt.MinimizeCost.suffixMap(rateCost, '');
            obj.replicaCost = opt.MinimizeCost.suffixMap(replicaCost, '_replicas');
            if nargin >= 4 && ~isempty(subjectTo)
                obj.constraints = subjectTo;
            end
        end

        function tf = isMinimization(~), tf = true; end

        function total = evaluate(obj, ~, variableValues)
            total = 0.0;
            % server costs
            ks = keys(obj.serverCost);
            for i = 1:numel(ks)
                if isKey(variableValues, ks{i})
                    val = variableValues(ks{i});
                    if opt.Objective.isScalarNumeric(val)
                        total = total + obj.serverCost(ks{i}) * val;
                    end
                end
            end
            % rate costs (pattern match)
            rks = keys(obj.rateCost);
            vks = keys(variableValues);
            for i = 1:numel(rks)
                pattern = rks{i};
                for j = 1:numel(vks)
                    vname = vks{j};
                    if ~isempty(strfind(vname, pattern)) && ~isempty(strfind(lower(vname), 'rate')) %#ok<STREMP>
                        val = variableValues(vname);
                        if opt.Objective.isScalarNumeric(val)
                            total = total + obj.rateCost(pattern) * val;
                        end
                    end
                end
            end
            % replica costs
            pks = keys(obj.replicaCost);
            for i = 1:numel(pks)
                if isKey(variableValues, pks{i})
                    val = variableValues(pks{i});
                    if opt.Objective.isScalarNumeric(val)
                        total = total + obj.replicaCost(pks{i}) * val;
                    end
                end
            end
        end
    end

    methods (Static)
        function out = suffixMap(inMap, suffix)
            out = containers.Map('KeyType', 'char', 'ValueType', 'double');
            if nargin < 1 || isempty(inMap)
                return;
            end
            ks = keys(inMap);
            for i = 1:numel(ks)
                out([ks{i} suffix]) = inMap(ks{i});
            end
        end
    end
end
