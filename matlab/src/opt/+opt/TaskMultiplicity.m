classdef TaskMultiplicity < opt.DecisionVariable
    % TaskMultiplicity  Optimize the multiplicity (thread/instance count) of an
    % LQN Task (integer). Mirrors native-Python TaskMultiplicity.

    properties
        task
        minValue
        maxValue
    end

    methods
        function obj = TaskMultiplicity(task, bounds, name)
            taskName = opt.Layered.elemName(task);
            if nargin < 3 || isempty(name)
                name = [taskName '_multiplicity'];
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
                % Tasks store multiplicity as a public property (no setter);
                % setting it directly mirrors how the LQN builder assigns it.
                task.multiplicity = round(value);
            end
        end

        function t = getVariableType(obj), t = 'task_multiplicity'; end

        function layers = getLayer(obj, model)
            layers = {obj.task};
            proc = opt.Layered.taskProcessorName(model, obj.task);
            if ~isempty(proc), layers{end+1} = proc; end
        end

        function v = currentValue(obj, model)
            v = [];
            task = opt.Layered.resolveTask(model, obj.task);
            if isempty(task), return; end
            m = task.multiplicity;
            if isnumeric(m) && isfinite(m), v = round(m); end
        end
    end
end
