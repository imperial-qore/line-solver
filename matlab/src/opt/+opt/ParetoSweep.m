classdef ParetoSweep < handle
    % ParetoSweep  Epsilon-constraint sweep for bi-objective tradeoff analysis.
    % Solves the problem once per epsilon, each time adding
    % constraintFactory(epsilon), and filters to the non-dominated cost
    % frontier. Mirrors native-Python ParetoSweep.
    %
    % constraintFactory is a function handle mapping an epsilon to an
    % opt.Constraint. solver is 'de' (default) or 'bisection'.

    properties
        problem
        constraintFactory
        epsilons
        solver = 'de';
        points = {};
    end

    methods
        function obj = ParetoSweep(problem, constraintFactory, epsilons, solver)
            obj.problem = problem;
            obj.constraintFactory = constraintFactory;
            obj.epsilons = epsilons;
            if nargin >= 4 && ~isempty(solver), obj.solver = solver; end
        end

        function p = getPoints(obj), p = obj.points; end

        function clone = cloneProblem(obj, epsilon)
            clone = opt.OptimizationProblem(obj.problem.getModel());
            vars = obj.problem.getVariables();
            for i = 1:numel(vars), clone.addVariable(vars{i}); end
            clone.setObjective(obj.problem.getObjective());
            cons = obj.problem.getConstraints();
            for i = 1:numel(cons), clone.addConstraint(cons{i}); end
            clone.setFixedVariables(obj.problem.getFixedVariables());
            scen = obj.problem.getScenarios();
            for i = 1:numel(scen), clone.addScenario(scen{i}{1}, scen{i}{2}); end
            clone.addConstraint(obj.constraintFactory(epsilon));
        end

        function pts = solve(obj, options)
            if nargin < 2, options = opt.LineOptSolverOptions(); end
            obj.points = {};
            for i = 1:numel(obj.epsilons)
                epsilon = obj.epsilons(i);
                clone = obj.cloneProblem(epsilon);
                if strcmp(obj.solver, 'bisection')
                    result = opt.BisectionSolver(clone).solve();
                else
                    result = clone.solve(options);
                end
                obj.points{end+1} = opt.ParetoPoint(epsilon, result.objectiveValue, ...
                    result.feasible, result); %#ok<AGROW>
            end
            pts = obj.points;
        end

        function frontier = getFrontier(obj)
            feasible = {};
            for i = 1:numel(obj.points)
                if obj.points{i}.feasible, feasible{end+1} = obj.points{i}; end %#ok<AGROW>
            end
            frontier = {};
            for i = 1:numel(feasible)
                p = feasible{i};
                dominated = false;
                for j = 1:numel(feasible)
                    q = feasible{j};
                    if q.objectiveValue <= p.objectiveValue && q.epsilon <= p.epsilon && ...
                            (q.objectiveValue < p.objectiveValue || q.epsilon < p.epsilon)
                        dominated = true; break;
                    end
                end
                if ~dominated, frontier{end+1} = p; end %#ok<AGROW>
            end
            % sort by epsilon
            if ~isempty(frontier)
                eps = cellfun(@(x) x.epsilon, frontier);
                [~, ord] = sort(eps);
                frontier = frontier(ord);
            end
        end
    end
end
