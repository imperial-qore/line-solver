classdef ActivityThinkTime < opt.DecisionVariable
    % ActivityThinkTime  Optimize the activity-level think time of an LQN
    % Activity (continuous). Mirrors native-Python ActivityThinkTime.

    properties
        activity
        minValue
        maxValue
    end

    methods
        function obj = ActivityThinkTime(activity, bounds, name)
            actName = opt.Layered.elemName(activity);
            if nargin < 3 || isempty(name)
                name = [actName '_thinktime'];
            end
            obj@opt.DecisionVariable(name);
            obj.activity = actName;
            obj.minValue = bounds(1);
            obj.maxValue = bounds(2);
        end

        function b = getBounds(obj), b = opt.DecisionVariable.unitBounds(1); end

        function v = decode(obj, x)
            v = obj.minValue + x(1) * (obj.maxValue - obj.minValue);
        end

        function apply(obj, model, value)
            act = opt.Layered.resolveActivity(model, obj.activity);
            if ~isempty(act)
                act.setThinkTime(double(value));
            end
        end

        function t = getVariableType(obj), t = 'think_time'; end

        function layers = getLayer(obj, model)
            layers = {};
            task = opt.Layered.taskOfActivity(model, obj.activity);
            if ~isempty(task), layers{end+1} = task.getName(); end
            proc = opt.Layered.activityProcessorName(model, obj.activity);
            if ~isempty(proc), layers{end+1} = proc; end
            if isempty(layers), layers = {}; end
        end

        function v = currentValue(obj, model)
            v = [];
            act = opt.Layered.resolveActivity(model, obj.activity);
            if ~isempty(act)
                v = opt.Layered.distMean(act.thinkTime);
            end
        end
    end
end
