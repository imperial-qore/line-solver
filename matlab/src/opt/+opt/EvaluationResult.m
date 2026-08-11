classdef EvaluationResult < handle
    % EvaluationResult  Metrics from evaluating a LINE model via SolverAuto.
    % Mirrors native-Python line_solver.opt.results.EvaluationResult. Per-
    % (station, class) metrics use a containers.Map keyed by 'station||class';
    % utilizations by station; system metrics by chain/class name.

    properties
        feasible = true;
        responseTimes             % Map 'station||class' -> value
        throughputs
        queueLengths
        utilizations              % Map 'station' -> value
        systemResponseTimes       % Map chain/class -> value
        systemThroughputs
        solveTime = 0.0;
        solverUsed = '';
        sensitivities = [];       % opt.SensitivityData or []
    end

    methods
        function obj = EvaluationResult()
            obj.responseTimes = containers.Map('KeyType', 'char', 'ValueType', 'double');
            obj.throughputs = containers.Map('KeyType', 'char', 'ValueType', 'double');
            obj.queueLengths = containers.Map('KeyType', 'char', 'ValueType', 'double');
            obj.utilizations = containers.Map('KeyType', 'char', 'ValueType', 'double');
            obj.systemResponseTimes = containers.Map('KeyType', 'char', 'ValueType', 'double');
            obj.systemThroughputs = containers.Map('KeyType', 'char', 'ValueType', 'double');
        end

        function k = key(~, station, jobclass)
            k = [station '||' jobclass];
        end

        function setResponseTime(obj, station, jobclass, v)
            obj.responseTimes(obj.key(station, jobclass)) = v;
        end
        function setThroughput(obj, station, jobclass, v)
            obj.throughputs(obj.key(station, jobclass)) = v;
        end
        function setQueueLength(obj, station, jobclass, v)
            obj.queueLengths(obj.key(station, jobclass)) = v;
        end

        function v = getResponseTime(obj, station, jobclass)
            if nargin >= 3 && ~isempty(jobclass)
                k = obj.key(station, jobclass);
                if isKey(obj.responseTimes, k), v = obj.responseTimes(k); else, v = inf; end
            else
                v = obj.aggregate(obj.responseTimes, station, true, inf);
            end
        end
        function v = getThroughput(obj, station, jobclass)
            if nargin >= 3 && ~isempty(jobclass)
                k = obj.key(station, jobclass);
                if isKey(obj.throughputs, k), v = obj.throughputs(k); else, v = 0.0; end
            else
                v = obj.aggregate(obj.throughputs, station, false, 0.0);
            end
        end
        function v = getQueueLength(obj, station, jobclass)
            if nargin >= 3 && ~isempty(jobclass)
                k = obj.key(station, jobclass);
                if isKey(obj.queueLengths, k), v = obj.queueLengths(k); else, v = 0.0; end
            else
                v = obj.aggregate(obj.queueLengths, station, false, 0.0);
            end
        end
        function v = getUtilization(obj, station)
            if isKey(obj.utilizations, station), v = obj.utilizations(station); else, v = 0.0; end
        end

        function v = aggregate(~, m, station, doMean, defaultVal)
            ks = keys(m);
            total = 0.0; count = 0;
            prefix = [station '||'];
            for i = 1:numel(ks)
                if strncmp(ks{i}, prefix, numel(prefix))
                    total = total + m(ks{i});
                    count = count + 1;
                end
            end
            if count == 0
                v = defaultVal;
            elseif doMean
                v = total / count;
            else
                v = total;
            end
        end

        function v = getSystemResponseTime(obj, jobclass)
            if nargin >= 2 && ~isempty(jobclass)
                if isKey(obj.systemResponseTimes, jobclass), v = obj.systemResponseTimes(jobclass); else, v = inf; end
                return;
            end
            ks = keys(obj.systemResponseTimes);
            if isempty(ks), v = inf; return; end
            totalTput = 0.0; weighted = 0.0;
            for i = 1:numel(ks)
                rt = obj.systemResponseTimes(ks{i});
                if isKey(obj.systemThroughputs, ks{i}), tp = obj.systemThroughputs(ks{i}); else, tp = 0.0; end
                weighted = weighted + rt * tp; totalTput = totalTput + tp;
            end
            if totalTput > 0
                v = weighted / totalTput;
            else
                vals = cell2mat(values(obj.systemResponseTimes));
                v = sum(vals) / numel(vals);
            end
        end

        function v = getSystemThroughput(obj, jobclass)
            if nargin >= 2 && ~isempty(jobclass)
                if isKey(obj.systemThroughputs, jobclass), v = obj.systemThroughputs(jobclass); else, v = 0.0; end
                return;
            end
            if obj.systemThroughputs.Count == 0, v = 0.0; else, v = sum(cell2mat(values(obj.systemThroughputs))); end
        end
    end
end
