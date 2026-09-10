classdef LineOptSolver < handle
    % LineOptSolver  Main line-opt solver. Minimizes a penalized scalar
    % objective (constraints as penalties) over the decision variables,
    % aggregating across scenarios, using the self-contained
    % opt.de.DifferentialEvolution engine (numpy-exact RNG) or an analytic/FD
    % projected-gradient path. Mirrors native-Python LineOptSolver.

    properties
        problem
        opt
        fixedValueMap
        freeVariables     % cell of opt.DecisionVariable (post-freeze)
        evaluators        % cell of opt.LineEvaluator
        scenarioWeights
        caches            % cell of dictionary
        iterations = 0;
        bestValue = inf;
        bestX = [];
        convergenceHistory = [];
        startTime
        deadline = inf;
        lqnSensCache      % dict valuesKey -> sensitivity dict (or [])
        gradCalls = 0;
    end

    methods
        function obj = LineOptSolver(problem, options)
            obj.problem = problem;
            obj.opt = options;
            obj.fixedValueMap = configureDictionary('string', 'cell');
            fixed = problem.getFixedVariables();
            for k = 1:numel(fixed)
                obj.fixedValueMap{fixed{k}{1}.getName()} = fixed{k}{2};
            end

            % Explicit layer freezing: variables whose layer is in frozenLayers
            % are held at the model's current parameter value (moved from free
            % to fixed) instead of being optimized. See freezeLayers.
            freeVars = problem.getVariables();
            [freeVars, fixed] = obj.freezeLayers(freeVars, fixed);
            obj.freeVariables = freeVars;

            obj.evaluators = {opt.LineEvaluator(problem.getModel(), freeVars, fixed)};
            scen = problem.getScenarios();
            w = zeros(1, 1 + numel(scen));
            w(1) = 1.0;
            for i = 1:numel(scen)
                obj.evaluators{end+1} = opt.LineEvaluator(scen{i}{1}, freeVars, fixed);
                w(i+1) = scen{i}{2};
            end
            obj.scenarioWeights = w;
            obj.lqnSensCache = configureDictionary('string', 'cell');
        end

        function [freeVars, fixed] = freezeLayers(obj, freeVars, fixed)
            % Partition variables by the frozenLayers option (LQN only). A
            % variable whose layer set (DecisionVariable.getLayer) intersects
            % the frozen set is moved from free to fixed, held at its current
            % model parameter value (DecisionVariable.currentValue). If the
            % current value cannot be read the variable is dropped from the
            % optimization, leaving the model's built-in value untouched.
            % No-op for a flat model or when no layers are frozen.
            frozen = obj.opt.frozenLayers;
            if isempty(frozen), return; end
            model = obj.problem.getModel();
            if ~opt.Layered.isLayered(model), return; end

            already = {};
            for k = 1:numel(fixed), already{end+1} = fixed{k}{1}.getName(); end %#ok<AGROW>
            kept = {};
            for i = 1:numel(freeVars)
                var = freeVars{i};
                layers = var.getLayer(model);
                if ~isempty(layers) && any(ismember(layers, frozen))
                    if any(strcmp(var.getName(), already))
                        continue;
                    end
                    value = var.currentValue(model);
                    if ~isempty(value)
                        fixed{end+1} = {var, value}; %#ok<AGROW>
                        obj.fixedValueMap{var.getName()} = value;
                    end
                    % else: drop; the model keeps its built-in value
                else
                    kept{end+1} = var; %#ok<AGROW>
                end
            end
            freeVars = kept;
        end

        function result = solve(obj)
            obj.startTime = tic;
            obj.iterations = 0;
            obj.bestValue = inf;
            obj.bestX = [];
            obj.convergenceHistory = [];
            obj.caches = cell(1, numel(obj.evaluators));
            for i = 1:numel(obj.evaluators)
                obj.caches{i} = configureDictionary('string', 'cell');
            end
            obj.lqnSensCache = configureDictionary('string', 'cell');
            obj.gradCalls = 0;
            obj.deadline = obj.opt.timeLimit;

            bounds = obj.evaluators{1}.getBounds();
            if isempty(bounds)
                result = obj.buildEmptyResult();
                return;
            end
            if obj.shouldUseGradient()
                result = obj.solveGradient(bounds);
            else
                result = obj.solveEvolution(bounds);
            end
        end

        function result = solveEvolution(obj, bounds)
            low = bounds(:, 1).';
            high = bounds(:, 2).';
            if isempty(obj.opt.seed)
                seed = mod(int64(feature('timing','cpucount')), 2^31);
            else
                seed = obj.opt.seed;
            end
            de = opt.de.DifferentialEvolution(@(x) obj.objectiveFunction(x), low, high, ...
                obj.opt.strategy, obj.opt.popsize, obj.opt.maxIterations, ...
                obj.opt.mutationLow, obj.opt.mutationHigh, obj.opt.recombination, obj.opt.tol, seed);
            de.callback = @(bx, nit) obj.deCallback(nit);

            timedOut = false;
            try
                r = de.solve();
                resX = obj.engineBestVector(r);
                resFun = r.fun;
            catch ME
                if strcmp(ME.identifier, 'LineOpt:TimeLimit')
                    timedOut = true;
                    if isempty(obj.bestX)
                        result = obj.buildEmptyResult(); return;
                    end
                    resX = obj.bestX; resFun = obj.bestValue;
                else
                    rethrow(ME);
                end
            end
            solveTime = toc(obj.startTime);
            result = obj.buildResult(resX, resFun, solveTime);
            if timedOut, result.terminatedBy = 'time_limit'; end
        end

        function stop = deCallback(obj, nit)
            obj.iterations = nit;
            obj.convergenceHistory(end+1) = obj.bestValue;
            stop = toc(obj.startTime) >= obj.deadline;
        end

        function x = engineBestVector(obj, r)
            if ~isempty(obj.bestX) && obj.bestValue <= r.fun
                x = obj.bestX;
            else
                x = r.x;
            end
        end

        function total = objectiveFunction(obj, x)
            if toc(obj.startTime) >= obj.deadline
                error('LineOpt:TimeLimit', 'time limit reached');
            end
            values = obj.evaluators{1}.decodeVariables(x);
            allValues = obj.mergeValues(values);
            objective = obj.problem.getObjective();
            pw = obj.opt.penaltyWeight;
            cons = obj.problem.getConstraints();
            scenarioValues = zeros(1, numel(obj.evaluators));
            for i = 1:numel(obj.evaluators)
                [res, obj.caches{i}] = obj.evaluators{i}.evaluateValuesWithCache(values, obj.caches{i});
                if ~res.feasible
                    total = inf; return;
                end
                value = objective.evaluateWithPenalty(res, allValues, pw);
                for c = 1:numel(cons)
                    value = value + cons{c}.evaluate(res, allValues) * pw;
                end
                scenarioValues(i) = value;
            end
            total = obj.aggregateScenarios(scenarioValues);
            if total < obj.bestValue
                obj.bestValue = total;
                obj.bestX = x;
            end
        end

        function av = mergeValues(obj, values)
            % merge: start from fixed, overlay values
            av = configureDictionary('string', 'cell');
            fk = keys(obj.fixedValueMap);
            for i = 1:numel(fk), av{fk(i)} = obj.fixedValueMap{fk(i)}; end
            vk = keys(values);
            for i = 1:numel(vk), av{vk(i)} = values{vk(i)}; end
        end

        function total = aggregateScenarios(obj, values)
            if numel(values) == 1
                total = values(1); return;
            end
            if strcmp(obj.opt.scenarioAggregation, 'mean')
                total = sum(obj.scenarioWeights .* values) / sum(obj.scenarioWeights);
            else
                total = max(values);
            end
        end

        function tf = shouldUseGradient(obj)
            if strcmp(obj.opt.optimizer, 'gradient'), tf = true; return; end
            if strcmp(obj.opt.optimizer, 'evolution'), tf = false; return; end
            tf = obj.allContinuous();
        end

        function tf = allContinuous(obj)
            % Continuous (differentiable) types, incl. the continuous LQN knobs
            % host demand and think time.
            vars = obj.freeVariables;
            if isempty(vars), tf = false; return; end
            tf = true;
            for i = 1:numel(vars)
                vt = vars{i}.getVariableType();
                if ~any(strcmp(vt, {'service_rate', 'routing', ...
                        'host_demand', 'think_time'}))
                    tf = false; return;
                end
            end
        end

        function result = solveGradient(obj, bounds)
            dim = size(bounds, 1);
            nStart = max(1, obj.opt.gradientRestarts);
            if isempty(obj.opt.seed), sd = 0; else, sd = obj.opt.seed; end
            rs = RandStream('mt19937ar', 'Seed', sd);
            starts = {0.5 * ones(1, dim)};
            for s = 2:nStart, starts{end+1} = rand(rs, 1, dim); end %#ok<AGROW>
            best = []; bestFun = inf; timedOut = false;
            for si = 1:numel(starts)
                if toc(obj.startTime) >= obj.deadline, timedOut = true; break; end
                try
                    xx = obj.projectedGradientDescent(starts{si});
                    fx = obj.objectiveFunction(xx);
                    if isfinite(fx) && fx < bestFun, bestFun = fx; best = xx; end
                catch ME
                    if strcmp(ME.identifier, 'LineOpt:TimeLimit'), timedOut = true; break; else, rethrow(ME); end
                end
            end
            if isempty(best) || (~isempty(obj.bestX) && obj.bestValue < bestFun)
                if ~isempty(obj.bestX), best = obj.bestX; bestFun = obj.bestValue;
                elseif isempty(best), result = obj.buildEmptyResult(); return; end
            end
            solveTime = toc(obj.startTime);
            result = obj.buildResult(best, bestFun, solveTime);
            if timedOut, result.terminatedBy = 'time_limit'; end
        end

        function x = projectedGradientDescent(obj, x0)
            x = min(max(x0, 0), 1);
            f = obj.objectiveFunction(x);
            for iter = 1:obj.opt.maxIterations
                g = obj.objectiveGradient(x);
                if norm(g) < 1e-9, break; end
                step = 1.0; improved = false;
                for ls = 1:30
                    xn = min(max(x - step * g, 0), 1);
                    fn = obj.objectiveFunction(xn);
                    if isfinite(fn) && fn < f - 1e-12
                        x = xn; f = fn; improved = true; break;
                    end
                    step = step * 0.5;
                end
                if ~improved, break; end
            end
        end

        function g = objectiveGradient(obj, x)
            % Gradient of the penalized objective in encoded space. For a
            % LayeredNetwork the source is selected by lqnGradient: 'fd'
            % (whole-model finite difference, robust default), 'partial_sens'
            % (SolverLN per-layer partial derivatives, cheap/biased), or
            % 'partial_plus_fd' (partial with a periodic full-FD correction).
            % A flat network always finite-differences. All paths fall back to
            % finite differences, which always work.
            if obj.evaluators{1}.isLayered
                mode = obj.opt.lqnGradient;
                if any(strcmp(mode, {'partial_sens', 'partial_plus_fd'}))
                    obj.gradCalls = obj.gradCalls + 1;
                    refresh = max(1, obj.opt.fdRefresh);
                    % partial_plus_fd periodically replaces the biased partial
                    % direction with the correct whole-model finite difference.
                    if ~(strcmp(mode, 'partial_plus_fd') && mod(obj.gradCalls, refresh) == 0)
                        g = obj.lqnAnalyticGradient(x);
                        if ~isempty(g), return; end
                    end
                end
                g = obj.finiteDifferenceGradient(x);
                return;
            end
            g = obj.finiteDifferenceGradient(x);
        end

        function g = finiteDifferenceGradient(obj, x)
            % Central finite-difference gradient of the penalized scalar
            % objective. Works for any model or solver (for a LayeredNetwork
            % each perturbed evaluation re-solves the whole ensemble, giving
            % the correct total derivative). One-sided differences near an
            % infeasible/unstable boundary where a two-sided value is non-finite.
            % A layered evaluation re-solves an iterative fixed point, so the
            % step must clear its noise floor (see fdStepLayered).
            if obj.evaluators{1}.isLayered
                h = obj.opt.fdStepLayered;
            else
                h = obj.opt.fdStep;
            end
            dim = numel(x);
            g = zeros(1, dim);
            f0 = [];
            for i = 1:dim
                xp = x; xm = x;
                xp(i) = min(1.0, x(i) + h);
                xm(i) = max(0.0, x(i) - h);
                fp = obj.objectiveFunction(xp);
                fm = obj.objectiveFunction(xm);
                if isfinite(fp) && isfinite(fm) && xp(i) > xm(i)
                    g(i) = (fp - fm) / (xp(i) - xm(i)); continue;
                end
                if isempty(f0), f0 = obj.objectiveFunction(x); end
                if isfinite(fp) && isfinite(f0) && xp(i) > x(i)
                    g(i) = (fp - f0) / (xp(i) - x(i));
                elseif isfinite(fm) && isfinite(f0) && x(i) > xm(i)
                    g(i) = (f0 - fm) / (x(i) - xm(i));
                else
                    g(i) = 0.0;
                end
            end
        end

        function g = lqnAnalyticGradient(obj, x)
            % Partial-sensitivity gradient for a LayeredNetwork ([] to fall back
            % to the whole-model finite difference). Assembles d(objective)/dx
            % from SolverLN's per-layer WITHIN-LAYER service-rate derivatives,
            % WITHOUT re-solving per parameter: (i) read d(layer metric)/d(rate)
            % at the variable's host-layer row, (ii) map each layer metric to
            % the LQN node metric it approximates and take d(penalized scalar)/
            % d(that node metric) by cheap metric-space finite differences,
            % (iii) chain through d(rate)/d(demand) = -1/D^2 and the linear
            % decode. Returns [] when multi-scenario, any variable is not a
            % host-demand variable, or the sensitivity table is unavailable.
            % BIASED (omits cross-layer coupling); fd/partial_plus_fd correct it.
            g = [];
            vars = obj.freeVariables;
            if numel(obj.evaluators) ~= 1 || isempty(vars), return; end
            for i = 1:numel(vars)
                if ~strcmp(vars{i}.getVariableType(), 'host_demand') ...
                        || ~ismethod(vars{i}, 'sensKey')
                    return;
                end
            end

            values = obj.evaluators{1}.decodeVariables(x);
            allValues = obj.mergeValues(values);
            [res, obj.caches{1}] = obj.evaluators{1}.evaluateValuesWithCache(values, obj.caches{1});
            if ~res.feasible, return; end

            skey = opt.LineEvaluator.valuesKey(values);
            if isKey(obj.lqnSensCache, skey)
                sens = obj.lqnSensCache{skey};
            else
                sens = obj.evaluators{1}.evaluateLayeredSensitivities(values);
                obj.lqnSensCache{skey} = sens;
            end
            if isempty(sens), return; end

            objective = obj.problem.getObjective();
            cons = obj.problem.getConstraints();
            pw = obj.opt.penaltyWeight;
            h = obj.opt.fdStep;
            model = obj.problem.getModel();
            kinds = {'Tput', 'RespT', 'QLen', 'Util'};

            grad = zeros(1, numel(x));
            for i = 1:numel(vars)
                var = vars{i};
                kv = var.sensKey(model);
                if isempty(kv) || ~isKey(sens, kv)
                    % No sensitivity row: leave 0 (other components still make
                    % progress; the FD modes cover it fully).
                    continue;
                end
                row = sens{kv};
                targets = var.sensMetricTargets(model);
                dS_drate = 0.0;
                for kk = 1:numel(kinds)
                    kind = kinds{kk};
                    dmetric_drate = row.(kind);
                    if dmetric_drate == 0, continue; end
                    mkey = targets.(kind);
                    dS_dmetric = obj.scalarMetricDerivative(res, allValues, ...
                        kind, mkey, objective, cons, pw, h);
                    dS_drate = dS_drate + dS_dmetric * dmetric_drate;
                end
                demand = allValues{var.getName()};
                if isnumeric(demand)
                    rateJac = var.rateJacobian(demand);
                else
                    rateJac = 0.0;
                end
                grad(i) = dS_drate * rateJac * var.decodeJacobian(x(i));
            end
            g = grad;
        end

        function d = scalarMetricDerivative(obj, res, allValues, kind, mkey, ...
                objective, cons, pw, h)
            % d(penalized scalar)/d(node metric[kind][mkey]) by central FD in
            % metric space (pure arithmetic, no solving).
            d = 0.0;
            if isempty(mkey), return; end
            fld = obj.metricFieldFor(kind);
            if isempty(fld) || ~isKey(res.(fld), mkey), return; end
            base = res.(fld)(mkey);
            res.(fld)(mkey) = base + h;
            fp = obj.scalarObjective(res, allValues, objective, cons, pw);
            res.(fld)(mkey) = base - h;
            fm = obj.scalarObjective(res, allValues, objective, cons, pw);
            res.(fld)(mkey) = base;
            d = (fp - fm) / (2.0 * h);
        end

        function v = scalarObjective(~, res, allValues, objective, cons, pw)
            v = objective.evaluateWithPenalty(res, allValues, pw);
            for c = 1:numel(cons)
                v = v + cons{c}.evaluate(res, allValues) * pw;
            end
        end

        function fld = metricFieldFor(~, kind)
            % property name, not the dictionary itself: a dictionary is a value
            % type, so the perturbation must be written back through res
            switch kind
                case 'RespT', fld = 'responseTimes';
                case 'QLen',  fld = 'queueLengths';
                case 'Tput',  fld = 'throughputs';
                case 'Util',  fld = 'utilizations';
                otherwise,    fld = '';
            end
        end

        function result = buildResult(obj, x, objectiveValue, solveTime)
            result = opt.OptimizationResult();
            result.objectiveValue = objectiveValue;
            vals = obj.evaluators{1}.decodeVariables(x);
            vk = keys(vals);
            for i = 1:numel(vk), result.variableValues{vk(i)} = vals{vk(i)}; end
            result.iterations = obj.iterations;
            result.solveTime = solveTime;
            evals = 0;
            for i = 1:numel(obj.evaluators), evals = evals + obj.evaluators{i}.getEvaluationCount(); end
            result.modelEvaluations = evals;
            result.convergenceHistory = obj.convergenceHistory;

            allValues = obj.mergeValues(result.variableValues);
            objective = obj.problem.getObjective();
            allCons = [objective.getConstraints(), obj.problem.getConstraints()];
            result.feasible = true;
            for i = 1:numel(obj.evaluators)
                [er, obj.caches{i}] = obj.evaluators{i}.evaluateValuesWithCache(result.variableValues, obj.caches{i});
                if ~er.feasible, result.feasible = false; continue; end
                for c = 1:numel(allCons)
                    viol = allCons{c}.evaluate(er, allValues);
                    if viol > 0
                        nm = allCons{c}.getName();
                        prev = 0.0; if isKey(result.constraintViolations, nm), prev = result.constraintViolations(nm); end
                        result.constraintViolations(nm) = max(viol, prev);
                        result.feasible = false;
                    end
                end
            end
            elapsed = toc(obj.startTime);
            if elapsed >= obj.opt.timeLimit
                result.terminatedBy = 'time_limit';
            elseif obj.iterations >= obj.opt.maxIterations
                result.terminatedBy = 'iterations';
            else
                result.terminatedBy = 'convergence';
            end
        end

        function result = buildEmptyResult(obj) %#ok<MANU>
            result = opt.OptimizationResult();
            result.objectiveValue = 0.0;
            result.feasible = true;
            result.terminatedBy = 'empty';
        end
    end
end
