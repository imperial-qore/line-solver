classdef OptimizationProblem < handle
    % OptimizationProblem  Declarative specification of a queueing-network
    % optimization problem. Mirrors native-Python
    % line_solver.opt.problem.OptimizationProblem.

    properties
        model
        variables = {};        % cell of opt.DecisionVariable
        objective = [];
        constraints = {};      % cell of opt.Constraint
        fixedVariables = {};   % cell of {var, value}
        scenarios = {};        % cell of {model, weight}
        isLayeredModel = false;
    end

    properties (Constant)
        % Decision-variable types operating on a LayeredNetwork vs a flat one.
        LQN_VAR_TYPES = {'host_demand','think_time','task_multiplicity', ...
            'task_replication','processor_multiplicity'};
        FLAT_VAR_TYPES = {'server_allocation','station_replicas', ...
            'service_rate','job_population','class_priority','routing', ...
            'class_mapping'};
    end

    methods
        function obj = OptimizationProblem(model)
            obj.model = model;
            obj.isLayeredModel = opt.Layered.isLayered(model);
        end

        function tf = isLayered(obj), tf = obj.isLayeredModel; end
        function m = getModel(obj), m = obj.model; end
        function v = getVariables(obj), v = obj.variables; end
        function o = getObjective(obj), o = obj.objective; end
        function c = getConstraints(obj), c = obj.constraints; end

        function obj = addVariable(obj, variable)
            obj.variables{end+1} = variable;
        end
        function obj = setObjective(obj, objective)
            obj.objective = objective;
        end
        function obj = addConstraint(obj, constraint)
            obj.constraints{end+1} = constraint;
        end
        function obj = setFixedVariables(obj, pairs)
            obj.fixedVariables = pairs;
        end
        function p = getFixedVariables(obj), p = obj.fixedVariables; end
        function obj = addScenario(obj, scenarioModel, weight)
            if nargin < 3, weight = 1.0; end
            obj.scenarios{end+1} = {scenarioModel, weight};
        end
        function s = getScenarios(obj), s = obj.scenarios; end

        function errors = validate(obj)
            errors = {};
            if isempty(obj.model), errors{end+1} = 'Model is not set'; end
            if isempty(obj.variables), errors{end+1} = 'No decision variables defined'; end
            if isempty(obj.objective), errors{end+1} = 'Objective function is not set'; end

            % Model/variable-kind consistency: LQN models take LQN variables and
            % flat models take flat variables; mixing them silently produces
            % no-ops (a variable whose element is never found in the copy).
            for i = 1:numel(obj.variables)
                vt = obj.variables{i}.getVariableType();
                nm = obj.variables{i}.getName();
                if obj.isLayeredModel && any(strcmp(vt, obj.FLAT_VAR_TYPES))
                    errors{end+1} = sprintf(['Variable ''%s'' (%s) is a ' ...
                        'flat-network variable but the model is a ' ...
                        'LayeredNetwork'], nm, vt); %#ok<AGROW>
                elseif ~obj.isLayeredModel && any(strcmp(vt, obj.LQN_VAR_TYPES))
                    errors{end+1} = sprintf(['Variable ''%s'' (%s) is a ' ...
                        'LayeredNetwork variable but the model is a flat ' ...
                        'Network'], nm, vt); %#ok<AGROW>
                end
            end
        end

        function tf = isValid(obj)
            tf = isempty(obj.validate());
        end

        function result = solve(obj, options)
            errors = obj.validate();
            if ~isempty(errors)
                line_error(mfilename, ['Invalid problem: ' strjoin(errors, ', ')]);
            end
            if nargin < 2 || isempty(options)
                options = opt.LineOptSolverOptions();
            end
            solver = opt.LineOptSolver(obj, options);
            result = solver.solve();
        end

        function workflow = decompose(obj)
            workflow = opt.DecompositionWorkflow(obj);
        end
    end
end
