classdef DecompositionWorkflow < handle
    % DecompositionWorkflow  Decomposes a joint problem into per-variable-type
    % subproblems solved via Gauss-Seidel cycling with fixed-value propagation.
    % An internal topological sort orders subproblems when dependencies are set.
    % Mirrors native-Python DecompositionWorkflow.

    properties (Constant)
        % Flat-network variable types first, then LayeredNetwork (LQN) types;
        % only types present in a given problem produce subproblems.
        DEFAULT_ORDER = {'server_allocation','station_replicas','service_rate', ...
            'job_population','class_priority','routing','class_mapping', ...
            'processor_multiplicity','task_multiplicity','task_replication', ...
            'host_demand','think_time'};
    end

    properties
        problem
        subproblems = {};
        dependencyGraph      % Map toNode -> cell of fromNodes
        solverOptions
    end

    methods
        function obj = DecompositionWorkflow(problem)
            obj.problem = problem;
            obj.dependencyGraph = containers.Map('KeyType','char','ValueType','any');
            obj.solverOptions = opt.LineOptSolverOptions();
        end

        function p = getProblem(obj), p = obj.problem; end
        function s = getSubProblems(obj), s = obj.subproblems; end
        function obj = setSolverOptions(obj, options), obj.solverOptions = options; end

        function obj = autoDecompose(obj)
            byType = containers.Map('KeyType','char','ValueType','any');
            vars = obj.problem.getVariables();
            for i = 1:numel(vars)
                t = vars{i}.getVariableType();
                if isKey(byType, t), lst = byType(t); else, lst = {}; end
                lst{end+1} = vars{i}; %#ok<AGROW>
                byType(t) = lst;
            end
            obj.subproblems = {};
            used = containers.Map('KeyType','char','ValueType','logical');
            for k = 1:numel(obj.DEFAULT_ORDER)
                t = obj.DEFAULT_ORDER{k};
                if isKey(byType, t)
                    obj.subproblems{end+1} = opt.SubProblem(t, t, byType(t)); %#ok<AGROW>
                    used(t) = true;
                end
            end
            bk = keys(byType);
            for k = 1:numel(bk)
                if ~isKey(used, bk{k})
                    obj.subproblems{end+1} = opt.SubProblem(bk{k}, bk{k}, byType(bk{k})); %#ok<AGROW>
                end
            end
        end

        function obj = setDependency(obj, fromProblem, toProblem)
            if isKey(obj.dependencyGraph, toProblem), s = obj.dependencyGraph(toProblem); else, s = {}; end
            s{end+1} = fromProblem;
            obj.dependencyGraph(toProblem) = s;
        end

        function obj = addSubProblem(obj, name, variables, after)
            if isempty(variables), vt = 'custom'; else, vt = variables{1}.getVariableType(); end
            obj.subproblems{end+1} = opt.SubProblem(name, vt, variables);
            if nargin >= 4 && ~isempty(after)
                for i = 1:numel(after), obj.setDependency(after{i}, name); end
            end
        end

        function ordered = getExecutionOrder(obj)
            if obj.dependencyGraph.Count == 0
                ordered = obj.subproblems; return;
            end
            names = cellfun(@(sp) sp.name, obj.subproblems, 'UniformOutput', false);
            indeg = containers.Map(names, num2cell(zeros(1, numel(names))));
            adj = containers.Map('KeyType','char','ValueType','any');
            for i = 1:numel(names), adj(names{i}) = {}; end
            tks = keys(obj.dependencyGraph);
            for i = 1:numel(tks)
                to = tks{i}; froms = obj.dependencyGraph(to);
                for j = 1:numel(froms)
                    from = froms{j};
                    if isKey(adj, from) && isKey(indeg, to)
                        a = adj(from); a{end+1} = to; adj(from) = a; %#ok<AGROW>
                        indeg(to) = indeg(to) + 1;
                    end
                end
            end
            q = {};
            for i = 1:numel(names), if indeg(names{i}) == 0, q{end+1} = names{i}; end; end %#ok<AGROW>
            orderNames = {}; head = 1;
            while head <= numel(q)
                n = q{head}; head = head + 1;
                orderNames{end+1} = n; %#ok<AGROW>
                succ = adj(n);
                for j = 1:numel(succ)
                    m = succ{j}; indeg(m) = indeg(m) - 1;
                    if indeg(m) == 0, q{end+1} = m; end %#ok<AGROW>
                end
            end
            if numel(orderNames) ~= numel(obj.subproblems)
                ordered = obj.subproblems; return;   % cycle: fall back
            end
            byName = containers.Map(names, obj.subproblems);
            ordered = cell(1, numel(orderNames));
            for i = 1:numel(orderNames), ordered{i} = byName(orderNames{i}); end
        end

        function result = solveSequential(obj, maxCycles, tolerance)
            if nargin < 2, maxCycles = 10; end
            if nargin < 3, tolerance = 0.01; end
            t0 = tic;
            result = opt.WorkflowResult();
            if isempty(obj.subproblems)
                result.converged = true; return;
            end
            ordered = obj.getExecutionOrder();
            fixedValues = containers.Map('KeyType','char','ValueType','any');
            prevObjective = inf;
            for cycle = 1:maxCycles
                for si = 1:numel(ordered)
                    sp = ordered{si};
                    partial = obj.createPartialProblem(sp, fixedValues);
                    spResult = opt.LineOptSolver(partial, obj.solverOptions).solve();
                    spr = opt.SubProblemResult(sp.name, spResult);
                    fk = keys(fixedValues);
                    for i = 1:numel(fk), spr.variablesFixed(fk{i}) = fixedValues(fk{i}); end
                    result.subproblemResults(sp.name) = spr;
                    vk = keys(spResult.variableValues);
                    for i = 1:numel(vk), fixedValues(vk{i}) = spResult.variableValues(vk{i}); end
                end
                currentObjective = obj.evaluateFullObjective(fixedValues);
                result.objectiveHistory(end+1) = currentObjective;
                if abs(currentObjective - prevObjective) < tolerance
                    result.converged = true; break;
                end
                prevObjective = currentObjective;
                result.cyclesCompleted = cycle;
            end
            if isempty(result.objectiveHistory)
                result.finalObjective = inf;
            else
                result.finalObjective = result.objectiveHistory(end);
            end
            fk = keys(fixedValues);
            for i = 1:numel(fk), result.finalVariableValues(fk{i}) = fixedValues(fk{i}); end
            result.totalSolveTime = toc(t0);
        end

        function result = solveHierarchical(obj)
            result = obj.solveSequential(1, 0.0);
        end

        function partial = createPartialProblem(obj, subproblem, fixedValues)
            partial = opt.OptimizationProblem(obj.problem.getModel());
            for i = 1:numel(subproblem.variables), partial.addVariable(subproblem.variables{i}); end
            subNames = subproblem.getVariableNames();
            fixedPairs = {};
            vars = obj.problem.getVariables();
            for i = 1:numel(vars)
                nm = vars{i}.getName();
                if ~any(strcmp(nm, subNames)) && isKey(fixedValues, nm)
                    fixedPairs{end+1} = {vars{i}, fixedValues(nm)}; %#ok<AGROW>
                end
            end
            partial.setFixedVariables(fixedPairs);
            partial.setObjective(obj.problem.getObjective());
            cons = obj.problem.getConstraints();
            for i = 1:numel(cons), partial.addConstraint(cons{i}); end
            scen = obj.problem.getScenarios();
            for i = 1:numel(scen), partial.addScenario(scen{i}{1}, scen{i}{2}); end
        end

        function value = evaluateFullObjective(obj, variableValues)
            penaltyWeight = obj.solverOptions.penaltyWeight;
            evaluator = opt.LineEvaluator(obj.problem.getModel(), obj.problem.getVariables(), {});
            res = evaluator.evaluateValues(variableValues);
            if ~res.feasible, value = inf; return; end
            objective = obj.problem.getObjective();
            value = objective.evaluateWithPenalty(res, variableValues, penaltyWeight);
            cons = obj.problem.getConstraints();
            for i = 1:numel(cons)
                value = value + cons{i}.evaluate(res, variableValues) * penaltyWeight;
            end
        end

        % ---- LayeredNetwork (LQN) layer-wise decomposition ---------------

        function result = solveLayered(obj, maxCycles, tolerance, autoFreeze, ...
                freezeTol, frozenLayers)
            % Solve an LQN by layer, optionally freezing converged layers.
            % Groups decision variables by the LQN layer they perturb (host or
            % task layer) and cycles Gauss-Seidel over the layer groups, fixing
            % every other layer's variables at their current values while one
            % layer is optimized. The LQN analogue of solveSequential, but the
            % subproblems are LAYERS rather than variable types.
            %
            % Freezing has two composable sources: frozenLayers (an explicit
            % seed set held fixed throughout) and autoFreeze (adaptive: after
            % each cycle a layer whose representative node metrics moved less
            % than freezeTol relative is frozen and skipped; unfrozen again if
            % any still-active layer later moves by more than freezeTol).
            % Convergence is on the full penalized objective delta, or when
            % every layer is frozen. Falls back to solveSequential for a flat
            % network. The WorkflowResult carries frozenLayers (final frozen
            % set) and modelEvaluations (total LINE solves).
            if nargin < 2 || isempty(maxCycles), maxCycles = 10; end
            if nargin < 3 || isempty(tolerance), tolerance = 0.01; end
            if nargin < 4 || isempty(autoFreeze), autoFreeze = true; end
            if nargin < 5 || isempty(freezeTol), freezeTol = 1e-3; end
            if nargin < 6 || isempty(frozenLayers), frozenLayers = {}; end

            t0 = tic;
            model = obj.problem.getModel();
            if ~opt.Layered.isLayered(model)
                result = obj.solveSequential(maxCycles, tolerance);
                return;
            end

            % Group variables by their primary (first) layer.
            groups = containers.Map('KeyType', 'char', 'ValueType', 'any');
            groupOrder = {};
            vars = obj.problem.getVariables();
            for i = 1:numel(vars)
                layers = vars{i}.getLayer(model);
                if isempty(layers), lkey = '_nolayer'; else, lkey = layers{1}; end
                if isKey(groups, lkey)
                    lst = groups(lkey);
                else
                    lst = {}; groupOrder{end+1} = lkey; %#ok<AGROW>
                end
                lst{end+1} = vars{i}; %#ok<AGROW>
                groups(lkey) = lst;
            end

            objective = obj.problem.getObjective();
            penaltyWeight = obj.solverOptions.penaltyWeight;
            evaluator = opt.LineEvaluator(model, vars, {});

            frozen = frozenLayers;
            fixedValues = containers.Map('KeyType', 'char', 'ValueType', 'any');
            prevSig = containers.Map('KeyType', 'char', 'ValueType', 'any');
            havePrevSig = false;
            prevObjective = inf;
            modelEvaluations = 0;

            result = opt.WorkflowResult();
            for cycle = 1:maxCycles
                for gi = 1:numel(groupOrder)
                    layer = groupOrder{gi};
                    if any(strcmp(layer, frozen)), continue; end
                    layerVars = groups(layer);
                    partial = opt.OptimizationProblem(model);
                    for i = 1:numel(layerVars), partial.addVariable(layerVars{i}); end
                    subNames = cellfun(@(v) v.getName(), layerVars, 'UniformOutput', false);
                    fixedPairs = {};
                    for i = 1:numel(vars)
                        nm = vars{i}.getName();
                        if ~any(strcmp(nm, subNames)) && isKey(fixedValues, nm)
                            fixedPairs{end+1} = {vars{i}, fixedValues(nm)}; %#ok<AGROW>
                        end
                    end
                    partial.setFixedVariables(fixedPairs);
                    partial.setObjective(objective);
                    cons = obj.problem.getConstraints();
                    for i = 1:numel(cons), partial.addConstraint(cons{i}); end

                    spResult = opt.LineOptSolver(partial, obj.solverOptions).solve();
                    modelEvaluations = modelEvaluations + spResult.modelEvaluations;
                    spr = opt.SubProblemResult(layer, spResult);
                    fk = keys(fixedValues);
                    for i = 1:numel(fk), spr.variablesFixed(fk{i}) = fixedValues(fk{i}); end
                    result.subproblemResults(layer) = spr;
                    vk = keys(spResult.variableValues);
                    for i = 1:numel(vk), fixedValues(vk{i}) = spResult.variableValues(vk{i}); end
                end

                % Full evaluation: objective + per-layer signatures for freezing.
                evalResult = evaluator.evaluateValues(fixedValues);
                modelEvaluations = modelEvaluations + 1;
                if ~evalResult.feasible
                    currentObjective = inf;
                    sig = containers.Map('KeyType', 'char', 'ValueType', 'any');
                else
                    currentObjective = objective.evaluateWithPenalty( ...
                        evalResult, fixedValues, penaltyWeight);
                    cons = obj.problem.getConstraints();
                    for i = 1:numel(cons)
                        currentObjective = currentObjective + ...
                            cons{i}.evaluate(evalResult, fixedValues) * penaltyWeight;
                    end
                    sig = opt.DecompositionWorkflow.layerSignatures(evalResult, groupOrder);
                end

                if autoFreeze && havePrevSig
                    moved = containers.Map('KeyType', 'char', 'ValueType', 'double');
                    for gi = 1:numel(groupOrder)
                        L = groupOrder{gi};
                        a = []; b = [];
                        if isKey(prevSig, L), a = prevSig(L); end
                        if isKey(sig, L), b = sig(L); end
                        moved(L) = opt.DecompositionWorkflow.sigDelta(a, b);
                    end
                    activeMoved = false;
                    for gi = 1:numel(groupOrder)
                        L = groupOrder{gi};
                        if ~any(strcmp(L, frozen)) && moved(L) > freezeTol
                            activeMoved = true; break;
                        end
                    end
                    for gi = 1:numel(groupOrder)
                        L = groupOrder{gi};
                        if any(strcmp(L, frozen))
                            if activeMoved, frozen(strcmp(frozen, L)) = []; end
                        elseif moved(L) < freezeTol
                            frozen{end+1} = L; %#ok<AGROW>
                        end
                    end
                end

                prevSig = sig; havePrevSig = true;
                result.objectiveHistory(end+1) = currentObjective;
                if abs(currentObjective - prevObjective) < tolerance
                    result.converged = true; break;
                end
                prevObjective = currentObjective;
                result.cyclesCompleted = cycle;
                % All layers frozen: nothing left to optimize.
                if numel(frozen) >= numel(groupOrder)
                    result.converged = true; break;
                end
            end

            if isempty(result.objectiveHistory)
                result.finalObjective = inf;
            else
                result.finalObjective = result.objectiveHistory(end);
            end
            fk = keys(fixedValues);
            for i = 1:numel(fk), result.finalVariableValues(fk{i}) = fixedValues(fk{i}); end
            result.totalSolveTime = toc(t0);
            result.frozenLayers = sort(frozen);
            result.modelEvaluations = modelEvaluations;
        end
    end

    methods (Static)
        function sig = layerSignatures(evalResult, layers)
            % Representative [Util, QLen, Tput, RespT] per layer, keyed by its
            % node. A layer named after a processor or task has a same-named
            % node in the LQN average table; its metrics are the layer's
            % convergence signature. Layers without a matching node (e.g. the
            % '_nolayer' bucket) get an empty signature so they never auto-freeze.
            sig = containers.Map('KeyType', 'char', 'ValueType', 'any');
            for i = 1:numel(layers)
                L = layers{i};
                util = evalResult.getUtilization(L);
                qlen = evalResult.getQueueLength(L);
                tput = evalResult.getThroughput(L);
                respt = evalResult.getResponseTime(L);
                probe = [util, qlen, tput];
                if any(probe ~= 0 & ~isinf(probe))
                    sig(L) = [util, qlen, tput, respt];
                else
                    sig(L) = [];
                end
            end
        end

        function d = sigDelta(a, b)
            % Max relative change between two layer signatures (inf if unknown).
            if isempty(a) || isempty(b), d = inf; return; end
            eps0 = 1e-12; d = 0.0;
            for i = 1:min(numel(a), numel(b))
                ai = a(i); bi = b(i);
                if ~(isfinite(ai) && isfinite(bi)), continue; end
                d = max(d, abs(bi - ai) / (abs(ai) + eps0));
            end
        end
    end
end
