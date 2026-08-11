classdef MaximizePerformance < opt.Objective
    % MaximizePerformance  Maximize weighted throughput + 1/RespT + 1/QLen,
    % subject to an optional budget constraint. Returns negated performance
    % (DE minimizes). Mirrors native-Python MaximizePerformance.

    properties
        throughputWeight = 1.0;
        responseTimeWeight = 1.0;
        queueLengthWeight = 0.0;
        stations = [];   % cell of station names, or [] = all
    end

    methods
        function obj = MaximizePerformance(throughputWeight, responseTimeWeight, ...
                queueLengthWeight, stations, budget, budgetTerms)
            if nargin >= 1 && ~isempty(throughputWeight), obj.throughputWeight = throughputWeight; end
            if nargin >= 2 && ~isempty(responseTimeWeight), obj.responseTimeWeight = responseTimeWeight; end
            if nargin >= 3 && ~isempty(queueLengthWeight), obj.queueLengthWeight = queueLengthWeight; end
            if nargin >= 4 && ~isempty(stations)
                names = cell(1, numel(stations));
                for i = 1:numel(stations)
                    names{i} = opt.Constraint.nameOf(stations{i});
                end
                obj.stations = names;
            end
            if nargin >= 5 && ~isempty(budget)
                if nargin < 6, budgetTerms = []; end
                obj.constraints = {opt.BudgetConstraint(budget, budgetTerms)};
            end
        end

        function tf = isMinimization(~), tf = true; end

        function v = evaluate(obj, result, ~)
            performance = 0.0;
            if isempty(obj.stations)
                % unique station names from throughput keys 'station||class'
                ks = keys(result.throughputs);
                sts = {};
                for i = 1:numel(ks)
                    parts = strsplit(ks{i}, '||');
                    sts{end+1} = parts{1}; %#ok<AGROW>
                end
                sts = unique(sts);
            else
                sts = obj.stations;
            end
            for i = 1:numel(sts)
                station = sts{i};
                if obj.throughputWeight > 0
                    performance = performance + obj.throughputWeight * result.getThroughput(station);
                end
                if obj.responseTimeWeight > 0
                    rt = result.getResponseTime(station);
                    if rt > 0 && isfinite(rt)
                        performance = performance + obj.responseTimeWeight * (1.0 / rt);
                    end
                end
                if obj.queueLengthWeight > 0
                    qlen = result.getQueueLength(station);
                    if qlen > 0
                        performance = performance + obj.queueLengthWeight * (1.0 / qlen);
                    end
                end
            end
            v = -performance;
        end
    end
end
