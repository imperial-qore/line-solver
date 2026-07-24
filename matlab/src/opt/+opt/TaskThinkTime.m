classdef TaskThinkTime < opt.DecisionVariable
    % TaskThinkTime  Optimize the think time of an LQN Task (continuous).
    % Mirrors native-Python TaskThinkTime.

    properties
        task
        minValue
        maxValue
    end

    methods
        function obj = TaskThinkTime(task, bounds, name)
            taskName = opt.Layered.elemName(task);
            if nargin < 3 || isempty(name)
                name = [taskName '_thinktime'];
            end
            obj@opt.DecisionVariable(name);
            obj.task = taskName;
            obj.minValue = bounds(1);
            obj.maxValue = bounds(2);
        end

        function b = getBounds(obj), b = opt.DecisionVariable.unitBounds(1); end

        function v = decode(obj, x)
            v = obj.minValue + x(1) * (obj.maxValue - obj.minValue);
        end

        function apply(obj, model, value)
            task = opt.Layered.resolveTask(model, obj.task);
            if ~isempty(task)
                task.setThinkTime(double(value));
            end
        end

        function t = getVariableType(obj), t = 'think_time'; end

        function layers = getLayer(obj, model)
            layers = {obj.task};
            proc = opt.Layered.taskProcessorName(model, obj.task);
            if ~isempty(proc), layers{end+1} = proc; end
        end

        function v = currentValue(obj, model)
            v = [];
            task = opt.Layered.resolveTask(model, obj.task);
            if ~isempty(task)
                v = opt.Layered.distMean(task.thinkTime);
            end
        end
    end
end
