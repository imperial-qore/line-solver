classdef LineEvaluator < handle
    % LineEvaluator  Interface between the optimizer and SolverAuto. Applies
    % decision variable values to a per-evaluation model copy, solves, and
    % extracts per-(station,class) and system metrics. Mirrors native-Python
    % line_solver.opt.evaluator.LineEvaluator.
    %
    % A copied model carries a cached NetworkStruct that scalar setters do not
    % invalidate, so evaluateValues forces refreshStruct after applying
    % variables.

    properties
        baseModel
        variables       % cell of opt.DecisionVariable
        fixedVariables  % cell of {var, value}
        evaluationCount = 0;
        totalDimension
        varOffsets
        isLayered = false;
    end

    methods
        function obj = LineEvaluator(model, variables, fixedVariables)
            obj.baseModel = model;
            obj.variables = variables;
            obj.isLayered = opt.Layered.isLayered(model);
            if nargin >= 3 && ~isempty(fixedVariables)
                obj.fixedVariables = fixedVariables;
            else
                obj.fixedVariables = {};
            end
            dim = 0;
            obj.varOffsets = zeros(1, numel(variables));
            for i = 1:numel(variables)
                obj.varOffsets(i) = dim;
                dim = dim + variables{i}.getDimension();
            end
            obj.totalDimension = dim;
        end

        function n = getEvaluationCount(obj), n = obj.evaluationCount; end

        function b = getBounds(obj)
            b = [];
            for i = 1:numel(obj.variables)
                b = [b; obj.variables{i}.getBounds()]; %#ok<AGROW>
            end
        end

        function values = decodeVariables(obj, x)
            values = containers.Map('KeyType', 'char', 'ValueType', 'any');
            for i = 1:numel(obj.variables)
                var = obj.variables{i};
                offset = obj.varOffsets(i);
                dim = var.getDimension();
                slice = x(offset+1 : offset+dim);
                values(var.getName()) = var.decode(slice);
            end
        end

        function applyVariables(obj, model, values)
            for i = 1:numel(obj.variables)
                var = obj.variables{i};
                if isKey(values, var.getName())
                    var.apply(model, values(var.getName()));
                end
            end
        end

        function m = copyModel(obj)
            m = obj.baseModel.copy();
        end

        function result = evaluateValues(obj, values)
            obj.evaluationCount = obj.evaluationCount + 1;
            result = opt.EvaluationResult();
            try
                model = obj.copyModel();
                for k = 1:numel(obj.fixedVariables)
                    fv = obj.fixedVariables{k};
                    fv{1}.apply(model, fv{2});
                end
                obj.applyVariables(model, values);

                if obj.isLayered
                    % LQN path: solve with SolverLN, read the per-node LQN
                    % table. Force a fresh struct so the mutated LQN elements
                    % are read (the base is never solved, so a cached struct
                    % would only appear via a deep-copied lsn).
                    model.lsn = [];
                    [~, avgTable] = opt.Layered.solveAvg(model);
                    result.feasible = true;
                    result.solverUsed = 'SolverLN';
                    obj.extractLayeredMetrics(avgTable, result);
                    obj.extractLayeredSystemMetrics(model, result);
                    % LQN sensitivities are expensive (they re-solve every
                    % layer), so the gradient path computes them lazily.
                else
                    model.refreshStruct();

                    solver = SolverAuto(model);
                    [QN, UN, RN, ~, ~, TN] = solver.getAvg();
                    result.feasible = true;
                    try
                        result.solverUsed = solver.getSelectedSolverName();
                    catch
                        result.solverUsed = '';
                    end
                    obj.extractMetrics(model, QN, UN, RN, TN, result);
                    obj.extractSystemMetrics(solver, model, result);
                    result.sensitivities = opt.sens.computeModelSensitivities(model);
                end
            catch ME %#ok<NASGU>
                result.feasible = false;
            end
        end

        function extractMetrics(~, model, QN, UN, RN, TN, result)
            sn = model.getStruct();
            R = sn.nclasses;
            for ist = 1:sn.nstations
                nodeIdx = sn.stationToNode(ist);
                st = sn.nodenames{nodeIdx};
                for r = 1:R
                    cl = sn.classnames{r};
                    if ~isempty(RN), result.setResponseTime(st, cl, RN(ist, r)); end
                    if ~isempty(TN), result.setThroughput(st, cl, TN(ist, r)); end
                    if ~isempty(QN), result.setQueueLength(st, cl, QN(ist, r)); end
                    if ~isempty(UN)
                        if isKey(result.utilizations, st)
                            result.utilizations(st) = result.utilizations(st) + UN(ist, r);
                        else
                            result.utilizations(st) = UN(ist, r);
                        end
                    end
                end
            end
        end

        function extractSystemMetrics(~, solver, model, result)
            try
                [SysRespT, SysTput] = solver.getAvgSys();
            catch
                return;
            end
            SysRespT = SysRespT(:); SysTput = SysTput(:);
            sn = model.getStruct();
            chains = sn.chains; classnames = sn.classnames;
            njobs = sn.njobs(:);
            nchains = size(chains, 1); nclasses = size(chains, 2);
            for ci = 1:nchains
                members = {};
                isOpen = false;
                for k = 1:nclasses
                    if chains(ci, k) > 0
                        members{end+1} = classnames{k}; %#ok<AGROW>
                        if isinf(njobs(k)), isOpen = true; end
                    end
                end
                if isempty(members), continue; end
                tput = 0.0; if ci <= numel(SysTput), tput = SysTput(ci); end
                respt = NaN; if ci <= numel(SysRespT), respt = SysRespT(ci); end
                if isOpen && tput > 0 && result.queueLengths.Count > 0
                    jobsInSystem = 0.0;
                    ks = keys(result.queueLengths);
                    for kk = 1:numel(ks)
                        parts = strsplit(ks{kk}, '||');
                        if any(strcmp(parts{2}, members))
                            jobsInSystem = jobsInSystem + result.queueLengths(ks{kk});
                        end
                    end
                    respt = jobsInSystem / tput;
                end
                for mi = 1:numel(members)
                    if ~isnan(respt), result.systemResponseTimes(members{mi}) = respt; end
                    result.systemThroughputs(members{mi}) = tput;
                end
            end
        end

        function extractLayeredMetrics(~, avgTable, result)
            % Extract per-LQN-node metrics from SolverLN's average table (one
            % row per Processor/Task/Entry/Activity; columns Node, NodeType,
            % QLen, Util, RespT, ResidT, ArvR, Tput). Throughput/queue-length/
            % response-time keyed by (node,node); utilization by node, matching
            % the flat EvaluationResult convention. Non-finite cells skipped.
            if isempty(avgTable) || height(avgTable) == 0
                return;
            end
            nodes = cellstr(avgTable.Node);
            for r = 1:numel(nodes)
                node = nodes{r};
                v = avgTable.RespT(r);
                if isfinite(v), result.setResponseTime(node, node, v); end
                v = avgTable.Tput(r);
                if isfinite(v), result.setThroughput(node, node, v); end
                v = avgTable.QLen(r);
                if isfinite(v), result.setQueueLength(node, node, v); end
                v = avgTable.Util(r);
                if isfinite(v), result.utilizations(node) = v; end
            end
        end

        function extractLayeredSystemMetrics(~, model, result)
            % Derive end-to-end (system) metrics from the reference task(s): a
            % closed LQN's system throughput is the reference task's
            % throughput; its end-to-end response time is the sum of response
            % times over the reference task's entries. Keyed by the reference
            % task's name so MinimizeSystemResponseTime /
            % SystemResponseTimeConstraint resolve without a chain concept.
            tasks = model.getTasks();
            for i = 1:numel(tasks)
                task = tasks{i};
                if ~opt.Layered.isRefTask(task)
                    continue;
                end
                tname = task.getName();
                tputKey = [tname '||' tname];
                if isKey(result.throughputs, tputKey)
                    tput = result.throughputs(tputKey);
                    if tput > 0
                        result.systemThroughputs(tname) = tput;
                    end
                end
                entries = task.entries;
                totalRt = 0.0; haveRt = false;
                for j = 1:numel(entries)
                    ename = entries(j).getName();
                    ekey = [ename '||' ename];
                    if isKey(result.responseTimes, ekey)
                        totalRt = totalRt + result.responseTimes(ekey);
                        haveRt = true;
                    end
                end
                if haveRt
                    result.systemResponseTimes(tname) = totalRt;
                end
            end
        end

        function sens = evaluateLayeredSensitivities(obj, values)
            % Compute LQN per-layer service-rate partial sensitivities on
            % demand: rebuild the configured model copy, solve it with
            % SolverLN, and reshape SolverLN.getSensitivityTable into a
            % containers.Map 'Station||JobClass' -> struct(Tput,RespT,QLen,
            % Util). Called only by the partial-sensitivity gradient path.
            % Returns [] on any failure (the caller then finite-differences).
            sens = [];
            if ~obj.isLayered
                return;
            end
            try
                model = obj.copyModel();
                for k = 1:numel(obj.fixedVariables)
                    fv = obj.fixedVariables{k};
                    fv{1}.apply(model, fv{2});
                end
                obj.applyVariables(model, values);
                model.lsn = [];
                solver = opt.Layered.makeSolver(model);
                sens = opt.Layered.computeSensitivities(solver);
            catch
                sens = [];
            end
        end

        function result = evaluateValuesWithCache(obj, values, cache)
            key = opt.LineEvaluator.valuesKey(values);
            if isKey(cache, key)
                result = cache(key);
            else
                result = obj.evaluateValues(values);
                cache(key) = result; %#ok<NASGU>
            end
        end
    end

    methods (Static)
        function key = valuesKey(values)
            ks = sort(keys(values));
            parts = cell(1, numel(ks));
            for i = 1:numel(ks)
                v = values(ks{i});
                if isnumeric(v)
                    parts{i} = [ks{i} '=' mat2str(round(v(:).' * 1e9) / 1e9)];
                else
                    parts{i} = [ks{i} '=' num2str(v)];
                end
            end
            key = strjoin(parts, ';');
        end
    end
end
