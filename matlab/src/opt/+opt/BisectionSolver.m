classdef BisectionSolver < handle
    % BisectionSolver  Exact O(log n) solver for a single integer decision
    % variable with monotone feasibility. direction='min_feasible' finds the
    % smallest feasible value (server sizing); 'max_feasible' the largest
    % (population sizing). Mirrors native-Python
    % line_solver.opt.sizing.BisectionSolver.

    properties
        problem
        direction = 'min_feasible';
        variable
        lo
        hi
        fixedValueDict
        evaluators
        baseResultAt = [];
        violationsAt
        probeCache
    end

    methods
        function obj = BisectionSolver(problem, direction)
            obj.problem = problem;
            if nargin >= 2 && ~isempty(direction), obj.direction = direction; end
            vars = problem.getVariables();
            if numel(vars) ~= 1 || vars{1}.getDimension() ~= 1
                line_error(mfilename, 'BisectionSolver requires exactly one dimension-1 decision variable');
            end
            obj.variable = vars{1};
            loV = obj.variable.decode(0.0);
            hiV = obj.variable.decode(1.0);
            if mod(loV,1) ~= 0 || mod(hiV,1) ~= 0
                line_error(mfilename, 'BisectionSolver requires an integer-valued variable');
            end
            obj.lo = round(loV); obj.hi = round(hiV);

            obj.fixedValueDict = configureDictionary('string','cell');
            fixed = problem.getFixedVariables();
            for k = 1:numel(fixed), obj.fixedValueDict{fixed{k}{1}.getName()} = fixed{k}{2}; end
            obj.evaluators = {opt.LineEvaluator(problem.getModel(), vars, fixed)};
            scen = problem.getScenarios();
            for i = 1:numel(scen)
                obj.evaluators{end+1} = opt.LineEvaluator(scen{i}{1}, vars, fixed); %#ok<AGROW>
            end
            obj.violationsAt = configureDictionary('string','double');
            obj.probeCache = configureDictionary('double','cell');
        end

        function cons = allConstraints(obj)
            objective = obj.problem.getObjective();
            cons = {};
            if ~isempty(objective), cons = objective.getConstraints(); end
            cons = [cons, obj.problem.getConstraints()];
        end

        function feasible = probe(obj, value)
            if isKey(obj.probeCache, value)
                c = obj.probeCache{value};
                obj.baseResultAt = c{2}; obj.violationsAt = c{3};
                feasible = c{1}; return;
            end
            values = configureDictionary('string','cell');
            values{obj.variable.getName()} = value;
            allValues = configureDictionary('string','cell');
            fk = keys(obj.fixedValueDict);
            for i = 1:numel(fk), allValues{fk(i)} = obj.fixedValueDict{fk(i)}; end
            allValues{obj.variable.getName()} = value;
            cons = obj.allConstraints();

            feasible = true;
            violations = configureDictionary('string','double');
            baseResult = [];
            for e = 1:numel(obj.evaluators)
                res = obj.evaluators{e}.evaluateValues(values);
                if isempty(baseResult), baseResult = res; end
                if ~res.feasible, feasible = false; continue; end
                for c = 1:numel(cons)
                    v = cons{c}.evaluate(res, allValues);
                    if v > 0
                        feasible = false;
                        nm = cons{c}.getName();
                        prev = 0.0; if isKey(violations, nm), prev = violations(nm); end
                        violations(nm) = max(v, prev);
                    end
                end
            end
            obj.baseResultAt = baseResult;
            obj.violationsAt = violations;
            obj.probeCache{value} = {feasible, baseResult, violations};
        end

        function result = solve(obj)
            t0 = tic;
            loV = obj.lo; hiV = obj.hi; iterations = 0;
            if strcmp(obj.direction, 'min_feasible')
                while loV < hiV
                    mid = floor((loV + hiV) / 2);
                    feasible = obj.probe(mid); iterations = iterations + 1;
                    if feasible, hiV = mid; else, loV = mid + 1; end
                end
            elseif strcmp(obj.direction, 'max_feasible')
                while loV < hiV
                    mid = floor((loV + hiV + 1) / 2);
                    feasible = obj.probe(mid); iterations = iterations + 1;
                    if feasible, loV = mid; else, hiV = mid - 1; end
                end
            else
                line_error(mfilename, ['Unknown direction: ' obj.direction]);
            end

            chosen = loV;
            feasible = obj.probe(chosen);

            result = opt.OptimizationResult();
            result.variableValues{obj.variable.getName()} = chosen;
            result.feasible = feasible;
            vk = keys(obj.violationsAt);
            for i = 1:numel(vk), result.constraintViolations(vk(i)) = obj.violationsAt(vk(i)); end
            result.iterations = iterations;
            evals = 0;
            for e = 1:numel(obj.evaluators), evals = evals + obj.evaluators{e}.getEvaluationCount(); end
            result.modelEvaluations = evals;
            result.terminatedBy = 'bisection';

            objective = obj.problem.getObjective();
            if ~isempty(objective) && ~isempty(obj.baseResultAt)
                allValues = configureDictionary('string','cell');
                fk = keys(obj.fixedValueDict);
                for i = 1:numel(fk), allValues{fk(i)} = obj.fixedValueDict{fk(i)}; end
                allValues{obj.variable.getName()} = chosen;
                result.objectiveValue = objective.evaluate(obj.baseResultAt, allValues);
            end
            result.solveTime = toc(t0);
        end
    end
end
