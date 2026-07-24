classdef TaskReplication < opt.DecisionVariable
    % TaskReplication  Optimize the replication (fan-out replicas) of an LQN
    % Task (integer). Mirrors native-Python TaskReplication.

    properties
        task
        minValue
        maxValue
    end

    methods
        function obj = TaskReplication(task, bounds, name)
            taskName = opt.Layered.elemName(task);
            if nargin < 3 || isempty(name)
                name = [taskName '_replication'];
            end
            obj@opt.DecisionVariable(name);
            obj.task = taskName;
            obj.minValue = round(bounds(1));
            obj.maxValue = round(bounds(2));
        end

        function b = getBounds(obj), b = opt.DecisionVariable.unitBounds(1); end

        function v = decode(obj, x)
            continuous = obj.minValue + x(1) * (obj.maxValue - obj.minValue);
            v = min(max(round(continuous), obj.minValue), obj.maxValue);
        end

        function apply(obj, model, value)
            task = opt.Layered.resolveTask(model, obj.task);
            if ~isempty(task)
                task.setReplication(round(value));
            end
        end

        function t = getVariableType(obj), t = 'task_replication'; end

        function layers = getLayer(obj, model)
            layers = {obj.task};
            proc = opt.Layered.taskProcessorName(model, obj.task);
            if ~isempty(proc), layers{end+1} = proc; end
        end

        function v = currentValue(obj, model)
            v = [];
            task = opt.Layered.resolveTask(model, obj.task);
            if isempty(task), return; end
            r = task.getReplication();
            if isnumeric(r) && isfinite(r), v = round(r); end
        end
    end
end
